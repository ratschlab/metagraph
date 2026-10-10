#include <algorithm>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <filesystem>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <tuple>
#include <vector>

#include <netinet/in.h>
#include <netinet/tcp.h>
#if defined(__APPLE__)
#include <netinet/tcp_fsm.h>
#endif
#include <sys/socket.h>
#include <zlib.h>
#include <json/json.h>
#include <server_http.hpp>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "server_utils.hpp"
#include "traverse.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;

/**
 * Compress a STL string using zlib with given compression level and return
 * the binary data.
 * Source: https://panthema.net/2007/0328-ZLibString.html
 */
std::string compress_string(const std::string &str, int compressionlevel, bool gzip,
                            const std::function<void()> &check, size_t max_piece) {
    z_stream zs; // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    // zlib's deflate stream by default; a gzip container (window bits + 16) when the
    // client asked for gzip, which is what most HTTP clients send by default
    const int window_bits = 15 + (gzip ? 16 : 0);
    if (deflateInit2(&zs, compressionlevel, Z_DEFLATED, window_bits, 8,
                     Z_DEFAULT_STRATEGY) != Z_OK)
        throw std::runtime_error("deflateInit failed while compressing.");

    // zlib counts its input in 32 bits (uInt avail_in): the text is handed over in pieces of
    // at most that many bytes, the last one with Z_FINISH. Assigned whole, a text of 4 GiB or
    // more would be cut to its size modulo 2^32, and the server would answer 200 with a
    // well-formed stream of that prefix. A text of one piece — every text below 4 GiB — is
    // compressed exactly as by a single call
    const size_t piece_limit = std::numeric_limits<uInt>::max();
    const size_t piece_max = max_piece ? std::min(max_piece, piece_limit) : piece_limit;
    const char *next = str.data();
    size_t left = str.size();

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // retrieve the compressed bytes blockwise
    do {
        if (check) {
            // a stop between blocks (an attempt at its bound, a client gone): the stream is
            // released before the exception leaves
            try {
                check();
            } catch (...) {
                deflateEnd(&zs);
                throw;
            }
        }
        if (!zs.avail_in && left) {
            // the next piece, once zlib consumed the last one (no input is added after
            // Z_FINISH: that flush is used only once nothing is left)
            const size_t piece = std::min(left, piece_max);
            zs.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(next));
            zs.avail_in = static_cast<uInt>(piece);
            next += piece;
            left -= piece;
        }
        zs.next_out = reinterpret_cast<Bytef *>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, left ? Z_NO_FLUSH : Z_FINISH);

        // append the block to the output string (what this call wrote: total_out is a uLong,
        // 32 bits on some platforms)
        outstring.append(outbuffer, sizeof(outbuffer) - zs.avail_out);
    } while (ret == Z_OK);

    deflateEnd(&zs);

    if (ret != Z_STREAM_END) { // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
        throw std::runtime_error(oss.str());
    }

    return outstring;
}

namespace {

// Request::connection is private and the pinned Simple-Web-Server has no accessor for it. An
// explicit instantiation may name a private member ([temp.spec.general]/6): the one
// conforming way to reach the connection's socket without patching the submodule. A
// submodule bump that renames or retypes the member fails to compile here, loudly.
auto request_connection(const HttpServer::Request &request);
template <auto Member>
struct ConnectionOf {
    friend auto request_connection(const HttpServer::Request &request) {
        return (request.*Member).lock();
    }
};
template struct ConnectionOf<&HttpServer::Request::connection>;

// The checks of one delivered text: |check| called whenever |out| reached the next multiple of
// the interval, and the longest time between two of them measured (from the start to the
// first, and from the last to finish(): what nothing interrupts, a token's preparation
// included)
class DeliveryChecks {
  public:
    DeliveryChecks(const std::string &out, const std::function<void()> &check, size_t interval,
                   double *max_gap_ms)
          : out_(out), check_(check), interval_(interval), next_(out.size() + interval),
            max_gap_ms_(max_gap_ms) {}
    // the bytes |out| may grow by before the next check is due (at least 1)
    size_t room() const { return next_ > out_.size() ? next_ - out_.size() : 1; }
    void tick() {
        if (out_.size() < next_)
            return;
        next_ = out_.size() + interval_;
        gap();
        check_();
        since_.reset();
    }
    void finish() { gap(); }

  private:
    void gap() {
        if (max_gap_ms_)
            *max_gap_ms_ = std::max(*max_gap_ms_, since_.elapsed() * 1000);
    }
    const std::string &out_;
    const std::function<void()> &check_;
    size_t interval_;
    size_t next_;
    double *max_gap_ms_;
    Timer since_;
};

// A streambuf appending to a string under DeliveryChecks: the JSON text of a large response is
// written under the caller's check, byte for byte what an std::ostringstream receives (the same
// writer writes into either). A piece the writer hands over is copied in pieces up to the next
// check, so that a check comes every |interval| bytes however large the piece (otherwise one
// 16 MiB string value would be appended whole, with one check after it)
class CheckedStringBuf : public std::streambuf {
  public:
    CheckedStringBuf(std::string *out, DeliveryChecks *checks) : out_(out), checks_(checks) {}

  protected:
    int_type overflow(int_type c) override {
        if (traits_type::eq_int_type(c, traits_type::eof()))
            return traits_type::not_eof(c);
        out_->push_back(traits_type::to_char_type(c));
        checks_->tick();
        return c;
    }
    std::streamsize xsputn(const char_type *s, std::streamsize n) override {
        for (size_t left = static_cast<size_t>(n); left; ) {
            const size_t piece = std::min(left, checks_->room());
            out_->append(s, piece);
            s += piece;
            left -= piece;
            checks_->tick();
        }
        return n;
    }

  private:
    std::string *out_;
    DeliveryChecks *checks_;
};

std::string trimmed(const std::string &s) {
    const size_t from = s.find_first_not_of(" \t");
    if (from == std::string::npos)
        return "";
    return s.substr(from, s.find_last_not_of(" \t") - from + 1);
}

std::string lowercase(std::string s) {
    for (char &c : s) {
        c = std::tolower(static_cast<unsigned char>(c));
    }
    return s;
}

// An RFC 9110 qvalue: "0" or "1", optionally with a point and up to three decimals,
// at most 1
std::optional<double> parse_qvalue(const std::string &s) {
    if (s.empty() || s.size() > 5 || (s[0] != '0' && s[0] != '1'))
        return std::nullopt;
    if (s.size() > 1 && s[1] != '.')
        return std::nullopt;
    for (size_t i = 2; i < s.size(); ++i) {
        if (!std::isdigit(static_cast<unsigned char>(s[i])))
            return std::nullopt;
    }
    const double q = std::stod(s);
    if (q > 1)
        return std::nullopt;
    return q;
}

#ifdef ASIO_STANDALONE
using TransportBuffer = asio::streambuf;
#else
using TransportBuffer = boost::asio::streambuf;
#endif

// Sizes the transport's buffer of |response| for |bytes| at once. Response::write copies the
// response into the Response's asio::streambuf, which grows 128 bytes at a time through
// std::vector::resize, that is by doubling: its last step holds an old and a new buffer of up
// to the body's size beside the body (a 68 MB body costs 129 MB in that copy). The pinned
// Simple-Web-Server has no call to size it; its Response is the std::ostream over that
// streambuf, so it is reached through rdbuf(). Called after the response's last check, so
// nothing reaches the transport before every check passed (a 503 at the bound or a client
// gone still writes nothing). A streambuf of another type (a submodule bump) is left to grow
// by itself
void reserve_transport(std::ostream &response, size_t bytes) {
    if (auto *buffer = dynamic_cast<TransportBuffer *>(response.rdbuf()))
        buffer->prepare(bytes);
}

// What Response::write puts before the body: the status line, the header fields and the
// Content-Length it adds, with room to spare (a shortfall would only grow the buffer once)
size_t response_head_bytes(const SimpleWeb::CaseInsensitiveMultimap &header) {
    size_t bytes = 256;
    for (const auto &[name, value] : header) {
        bytes += name.size() + value.size() + 4;
    }
    return bytes;
}

} // namespace

// The content encoding to send, "" for none (RFC 9110 §12.5.3). Accept-Encoding is a
// comma-separated list of case-insensitive codings, each with an optional weight ";q="
// (1 when absent, 0 = not acceptable), `*` standing for every coding not listed by name;
// several header lines form one list, and an element with a malformed weight is ignored.
// A coding is acceptable at its own weight, else at the weight of `*`, else not at all.
// gzip wins over deflate at equal weight, the higher weight otherwise; the response is
// uncompressed when neither is acceptable or the client weights identity above both
// (uncompressed is also the fallback after "identity;q=0": refusing to answer is not what
// such a client asked for). A substring test would send gzip for "gzip;q=0, deflate;q=1"
// and for "identity, gzip;q=0".
std::string requested_encoding(const std::shared_ptr<HttpServer::Request> &request) {
    const auto [from, to] = request->header.equal_range("Accept-Encoding");
    if (from == to)
        return "";
    std::map<std::string, double> weight;   // the first weight given for each coding
    for (auto it = from; it != to; ++it) {
        std::istringstream list(it->second);
        std::string element;
        while (std::getline(list, element, ',')) {
            std::istringstream parts(element);
            std::string coding, param;
            std::getline(parts, coding, ';');
            coding = lowercase(trimmed(coding));
            if (coding.empty())
                continue;
            std::optional<double> q = 1.0;
            while (q && std::getline(parts, param, ';')) {
                param = trimmed(param);
                if (param.size() >= 2 && std::tolower(static_cast<unsigned char>(param[0])) == 'q'
                        && param[1] == '=') {
                    q = parse_qvalue(trimmed(param.substr(2)));
                }
            }
            if (q)
                weight.emplace(coding == "x-gzip" ? "gzip" : coding, *q);
        }
    }
    auto weight_of = [&](const std::string &coding) {
        if (auto it = weight.find(coding); it != weight.end())
            return it->second;
        if (auto it = weight.find("*"); it != weight.end())
            return it->second;
        return 0.0;
    };
    const double gzip = weight_of("gzip");
    const double deflate = weight_of("deflate");
    const double best = std::max(gzip, deflate);
    if (best <= 0)
        return "";
    if (auto it = weight.find("identity"); it != weight.end() && it->second > best)
        return "";
    return gzip >= deflate ? "gzip" : "deflate";
}

bool is_compression_requested(const std::shared_ptr<HttpServer::Request> &request) {
    return !requested_encoding(request).empty();
}

bool client_gone(const HttpServer::Request &request) {
    auto connection = request_connection(request);
    if (!connection)
        return true;
    return peer_closed(connection->socket->lowest_layer().native_handle());
}

// Whether the TCP connection |fd| is past its established state, with data still waiting
// before the end: a peek returns that data, never the end behind it. The kernel's connection
// state says it. A handler's connection leaves ESTABLISHED only by a FIN or a reset of the
// peer (CLOSE_WAIT, CLOSED) or by the server's own shutdown — the HTTP server's content
// timeout shuts the connection (FIN_WAIT1/2, CLOSING, TIME_WAIT, LAST_ACK) — and in every one
// of these no response can be delivered. Reading only the peer's states is not enough: on
// Linux, whose shutdown keeps the waiting bytes (macOS discards them, so the peek reads the
// end), a walk that outlived the content timeout with bytes waiting would compute on to its
// end, against the SPEC's "stopped by this too". SYN_RECV stays connected: an accepted
// connection is there only under TCP Fast Open, which this server's listener does not enable,
// and a live client must never read as gone. Not a TCP socket, or a platform without the
// query: false (connected)
static bool tcp_peer_finished(int fd) {
#if defined(__linux__)
    struct tcp_info info;
    socklen_t length = sizeof(info);
    if (::getsockopt(fd, IPPROTO_TCP, TCP_INFO, &info, &length) != 0)
        return false;
    return info.tcpi_state != TCP_ESTABLISHED && info.tcpi_state != TCP_SYN_RECV;
#elif defined(__APPLE__)
    struct tcp_connection_info info;
    socklen_t length = sizeof(info);
    if (::getsockopt(fd, IPPROTO_TCP, TCP_CONNECTION_INFO, &info, &length) != 0)
        return false;
    return info.tcpi_state != TCPS_ESTABLISHED && info.tcpi_state != TCPS_SYN_RECEIVED;
#else
    (void)fd;
    return false;
#endif
}

bool peer_closed(int fd) {
    if (fd < 0)
        return true;
    char byte;
    ssize_t n;
    do {
        n = ::recv(fd, &byte, 1, MSG_PEEK | MSG_DONTWAIT);
    } while (n < 0 && errno == EINTR);
    if (n == 0)
        return true;     // closed: an orderly close or a half-close
    if (n > 0) {
        // Data waiting — a pipelined request, or bytes past the request such as a trailing
        // CRLF (RFC 9112 §2.2) — hides a close behind it from the peek: the connection's
        // state tells (otherwise a client that closed after a trailing CRLF would be
        // walked to the end and answered into a dead socket)
        return tcp_peer_finished(fd);
    }
    // EBADF and ENOTSOCK: the descriptor holds no connection (closed, or not a socket), so no
    // client can be answered through it. The server's descriptor is always its connection's
    // socket while the handler holds the request, so a request never meets these
    return errno == ECONNRESET || errno == ENOTCONN || errno == EPIPE || errno == ETIMEDOUT
        || errno == EBADF || errno == ENOTSOCK;
}

std::string json_text(const Json::Value &value, bool compact, const std::function<void()> &check,
                      double *max_gap_ms) {
    Json::StreamWriterBuilder builder;
    if (compact)
        builder["indentation"] = "";   // the traversal routes: half the bytes of the indented form
    if (!check)
        return Json::writeString(builder, value);
    std::string out;
    DeliveryChecks checks(out, check, kDeliveryCheckBytes, max_gap_ms);
    CheckedStringBuf buf(&out, &checks);
    std::ostream stream(&buf);
    // an exception of the check reaches the caller rather than setting badbit
    stream.exceptions(std::ios::badbit);
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    writer->write(value, &stream);
    checks.finish();
    return out;
}

std::string assemble_traverse_response(const Json::Value &envelope,
                                       std::vector<std::string> results,
                                       const std::function<void()> &check,
                                       double *max_gap_ms) {
    Json::Value before(Json::objectValue);
    Json::Value after(Json::objectValue);
    for (const std::string &name : envelope.getMemberNames()) {
        if (name == "results")
            throw std::logic_error("assemble_traverse_response: the envelope holds results");
        (name < "results" ? before : after)[name] = envelope[name];
    }
    std::string head = json_text(before, true, check, max_gap_ms);     // {...}
    // the envelope's members after "results" are written first (the usage, the timing: a few
    // KiB), so that the whole response's size is known before a byte of the results is copied
    const std::string tail = json_text(after, true, check, max_gap_ms);  // {...}
    head.pop_back();
    const bool head_members = head.size() > 1;
    const bool tail_members = tail.size() > 2;
    static const std::string kResults = "\"results\":[";
    // the exact size: `{` and the members before, a comma, "results":[, the texts and their
    // commas, `]`, then a comma and the members after with their `}`, or the `}` alone. Reserved
    // at once, the text never grows by doubling, whose last step held an old and a new buffer
    // of up to its size beside the texts
    size_t size = head.size() + head_members + kResults.size() + 1
                + (tail_members ? tail.size() : 1);
    for (size_t i = 0; i < results.size(); ++i) {
        size += results[i].size() + (i > 0);
    }
    std::string out;
    out.reserve(size);
    out += head;
    std::string().swap(head);
    if (head_members)
        out += ',';
    out += kResults;
    // each text is freed once copied: what was copied and what is left to copy are the
    // response once; texts kept alive beside it until the handler returns would stay through
    // the compression and the transport's copy (2.5x the text held with gzip, 3.5x without)
    if (!check) {
        for (size_t i = 0; i < results.size(); ++i) {
            if (i)
                out += ',';
            out += results[i];
            std::string().swap(results[i]);
        }
    } else {
        // each text copied in pieces up to the next check
        DeliveryChecks checks(out, check, kDeliveryCheckBytes, max_gap_ms);
        for (size_t i = 0; i < results.size(); ++i) {
            if (i)
                out += ',';
            const std::string &text = results[i];
            for (size_t at = 0; at < text.size(); ) {
                const size_t piece = std::min(text.size() - at, checks.room());
                out.append(text, at, piece);
                at += piece;
                checks.tick();
            }
            std::string().swap(results[i]);
        }
        checks.finish();
    }
    out += ']';
    if (tail_members) {
        out += ',';
        out.append(tail, 1, std::string::npos);
    } else {
        out += '}';
    }
    // the size computed is the size written: a mismatch would only cost a reallocation, never
    // a byte, so it is asserted where tests run unoptimised
    assert(out.size() == size);
    return out;
}

GraphListEntry parse_graph_list_line(const std::string &text, size_t line) {
    std::vector<std::string> columns;
    size_t from = 0;
    while (true) {
        const size_t comma = text.find(',', from);
        columns.push_back(text.substr(from, comma == std::string::npos ? std::string::npos
                                                                       : comma - from));
        if (comma == std::string::npos)
            break;
        from = comma + 1;
    }
    auto bad = [&](const std::string &what) {
        return std::invalid_argument("line " + std::to_string(line) + " of the graph list ('"
                                     + text + "'): " + what);
    };
    if (columns.size() < 3)
        throw bad("expected name,graph_path,annotation_path[,manifest_path[,index_ns]]");
    if (columns.size() > 5) {
        throw bad(std::to_string(columns.size()) + " columns, at most five are read "
                  "(name,graph_path,annotation_path,manifest_path,index_ns): a column the "
                  "server would not read is refused rather than ignored");
    }
    GraphListEntry e;
    e.line = line;
    e.name = columns[0];
    e.graph_path = columns[1];
    e.annotation_path = columns[2];
    if (columns.size() > 3)
        e.manifest_path = columns[3];
    if (columns.size() > 4) {
        e.index_ns = columns[4];
        // the name is an MGT token and goes into file names (as --index-name)
        bool ok = !e.index_ns.empty();
        for (unsigned char c : e.index_ns) {
            ok &= std::isalnum(c) || c == '.' || c == '_' || c == '-';
        }
        if (!e.index_ns.empty() && !ok)
            throw bad("index_ns '" + e.index_ns + "' does not match [A-Za-z0-9._-]+");
    }
    return e;
}

std::map<std::pair<std::string, std::string>, std::pair<std::string, std::string>>
graph_list_identities(const std::vector<GraphListEntry> &entries,
                      const std::function<std::string(const GraphListEntry &)> &fingerprint,
                      const std::function<std::vector<std::string>(const GraphListEntry &)> &inventory) {
    using Pair = std::pair<std::string, std::string>;
    struct Stated {
        std::string value;
        size_t line = 0;
    };
    std::map<Pair, std::pair<Stated, Stated>> stated;     // (ns, fp) and the lines stating them
    std::map<std::pair<Pair, std::string>, std::string> digests;   // (pair, manifest) -> fp
    auto agree = [](Stated *have, const std::string &value, size_t line, const char *what,
                    const Pair &pair) {
        if (value.empty())
            return;
        if (have->value.empty()) {
            *have = { value, line };
            return;
        }
        if (have->value != value) {
            throw std::invalid_argument(
                    std::string("lines ") + std::to_string(have->line) + " and "
                    + std::to_string(line) + " of the graph list state different " + what
                    + " for the same index (" + pair.first + ", " + pair.second + "): '"
                    + have->value + "' and '" + value + "'; one index has one identity");
        }
    };
    // One index is one set of loaded files, however its paths are spelled: the lines are
    // grouped by the real paths of the pair's COMPLETE loader inventory (index_load_inventory:
    // the main files and every identity sidecar, derived from the spelling as the loaders derive
    // them; the graph's derived data, its mask and Bloom filter, is not part of an identity),
    // so that two spellings of one index agree, and two pairs whose main files are symlinks to
    // the same files but whose sidecars differ are two indexes, each validated against its
    // manifest (grouped by the main files' real paths, the second pair's own .seqs would never
    // be checked and both would state one index_fp)
    auto real = [](const std::string &path) {
        std::error_code ec;
        const std::filesystem::path p = std::filesystem::weakly_canonical(path, ec);
        return ec ? path : p.string();
    };
    // the bundle of a listed pair: its files' real paths in the inventory's order, and the
    // real main files (for messages)
    using Bundle = std::vector<std::string>;
    std::map<Pair, Bundle> bundle_of;                      // the pair as listed -> its bundle
    std::map<Bundle, Pair> pair_of;                        // a bundle -> the first pair listing it
    for (const GraphListEntry &e : entries) {
        const Pair listed { e.graph_path, e.annotation_path };
        auto [b, fresh_pair] = bundle_of.try_emplace(listed);
        if (fresh_pair) {
            const std::vector<std::string> files
                = inventory ? inventory(e) : index_bundle_files(e.graph_path, e.annotation_path);
            for (const std::string &f : files) {
                b->second.push_back(real(f));
            }
            pair_of.try_emplace(b->second, listed);
        }
        const Bundle &bundle = b->second;
        // keyed by the bundle's first listing: the messages name a pair as listed
        const Pair pair = pair_of.at(bundle);
        auto &[ns, fp] = stated[pair];
        agree(&ns, e.index_ns, e.line, "index_ns", listed);
        if (e.manifest_path.empty())
            continue;
        // every distinct bundle is validated against each manifest a line gives for it
        auto [it, fresh] = digests.try_emplace({ pair, e.manifest_path });
        if (fresh)
            it->second = fingerprint(e);
        agree(&fp, it->second, e.line, "index_fp (manifest digests)", listed);
    }
    // Two different indexes never state one fingerprint: the server checks sizes, not
    // contents, so one manifest whose files have the base names and sizes of another pair's
    // (two annotations of one graph with swapped memberships, written by one tool run) would
    // let a client take one index for the other. Copies of one index under
    // two paths are refused too (telling them apart would need hashing): list one path
    std::map<std::string, std::pair<Pair, size_t>> by_fp;
    for (const auto &[pair, s] : stated) {
        if (s.second.value.empty())
            continue;
        auto [it, fresh] = by_fp.try_emplace(s.second.value, pair, s.second.line);
        if (!fresh) {
            throw std::invalid_argument(
                    "lines " + std::to_string(it->second.second) + " and "
                    + std::to_string(s.second.line) + " of the graph list state one index_fp "
                    + s.second.value + " for two different indexes (" + it->second.first.first
                    + ", " + it->second.first.second + ") and (" + pair.first + ", "
                    + pair.second + "): a manifest states the identity of one graph with one "
                    "annotation; give each pair its own manifest (scripts/traversal/"
                    "index_manifest.py --server-csv writes one per pair), or list copies of one "
                    "index by one path");
        }
    }
    std::map<Pair, std::pair<std::string, std::string>> out;
    for (const auto &[listed, bundle] : bundle_of) {
        const auto &s = stated.at(pair_of.at(bundle));
        out[listed] = { s.first.value, s.second.value };
    }
    return out;
}

Json::Value parse_json_string(const std::string &msg) {
    Json::Value json;

    Json::CharReaderBuilder rbuilder;
    std::unique_ptr<Json::CharReader> reader { rbuilder.newCharReader() };
    std::string errors;

    if (!reader->parse(msg.data(), msg.data() + msg.size(), &json, &errors))
        throw std::invalid_argument("Bad json received: " + errors);

    return json;
}

std::string json_str_with_error_msg(const std::string &msg) {
    Json::Value root;
    root["error"] = msg;
    return Json::writeString(Json::StreamWriterBuilder(), root);
}

RequestAnswer answer_request(const std::string &content, const std::string &encoding,
                             size_t request_id,
                             const std::function<Json::Value(const std::string &)> &process,
                             bool compact,
                             const ResponseControl *control) {
    RequestAnswer answer;
    std::string &ret = answer.body;
    static const std::function<void()> kNoCheck;
    const std::function<void()> &check = control ? control->check : kNoCheck;

    try {
        // Return JSON string
        answer.status = 200;
        if (control && control->write) {
            ret = control->write(process(content), check);
        } else {
            ret = json_text(process(content), compact, check);
        }
        if (!encoding.empty()) {
            Timer compressing;
            const size_t text_bytes = ret.size();
            ret = compress_string(ret, control ? control->compression_level : Z_BEST_COMPRESSION,
                                  encoding == "gzip", check);
            if (control && control->on_compressed)
                control->on_compressed(text_bytes, compressing.elapsed());
            answer.header.emplace_back("Content-Encoding", encoding);
            answer.header.emplace_back("Content-Length", std::to_string(ret.size()));
        }
        // once more before the response is handed to the transport: a body shorter than the
        // check's interval is otherwise never checked after it was built
        if (check)
            check();
    } catch (const ClientGone &e) {
        // nobody to answer: nothing is written, and the connection is closed rather than
        // kept for a next request
        logger->info("[Server] Request {}: {}; no response written", request_id, e.what());
        return RequestAnswer();
    } catch (const HttpError &e) {
        logger->warn("[Server] Error on request {} ({}): {}", request_id, e.status(), e.what());
        answer.status = e.status();
        answer.header.clear();
        // the body as the route wrote it, uncompressed like every error
        ret = json_text(e.body(), compact);
    } catch (const CurrentlyInitializingError& e) {
        logger->info("[Server] Got a request during initialization. Asked to come back later");
        answer.status = 503;
        answer.header.emplace_back("Retry-After", "60"); // ask to come back in 60 seconds
        ret = json_str_with_error_msg("Server is currently initializing, please come back later.");
    } catch (const std::exception& e) {
        logger->warn("[Server] Error on request {}: {}", request_id, e.what());
        answer.status = 400;
        ret = json_str_with_error_msg(e.what());
    } catch (...) {
        logger->warn("[Server] Error on request {}", request_id);
        answer.status = 500;
        ret = json_str_with_error_msg("Internal server error");
    }
    // a route that answers nobody who left asks once more, whatever was built: its errors
    // too (the success path's check asked already; an error never passed one). Apart from
    // the deadline, so that a 503 at the deadline still reaches a client that is there
    if (control && control->gone && control->gone()) {
        logger->info("[Server] Request {}: the client is gone or the server stops; the {} "
                     "answer is not written", request_id, answer.status);
        return RequestAnswer();
    }
    return answer;
}

void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     size_t request_id,
                     const std::function<Json::Value(const std::string &)> &process,
                     bool compact,
                     const ResponseControl *control) {
    logger->info("[Server] {} request {} from {}", request->path, request_id,
                 request->remote_endpoint().address().to_string());
    Timer timer;
    const RequestAnswer answer = answer_request(request->content.string(),
                                                requested_encoding(request), request_id,
                                                process, compact, control);
    if (!answer.status) {
        // nothing is written, and the connection is closed rather than kept for a next
        // request
        response->close_connection_after_response = true;
        if (control && control->on_written)
            control->on_written(0, std::nullopt);
        return;
    }
    const auto status = static_cast<SimpleWeb::StatusCode>(answer.status);
    // the Content-Type first, then the fields in the order added (the multimap built as it
    // always was, so that its fields are written in the same order)
    SimpleWeb::CaseInsensitiveMultimap header({ { "Content-Type", "application/json" } });
    for (const auto &field : answer.header) {
        header.insert(field);
    }
    const std::string &ret = answer.body;
    double processing_time = timer.elapsed();
    // after the last check: the transport's copy is made into a buffer of its final size
    reserve_transport(*response, ret.size() + response_head_bytes(header));
    response->write(status, ret, header);
    if (control && control->on_written)
        control->on_written(static_cast<int>(status), ret.size());
    logger->info("[Server] Request {} processing time: {:.3f} sec, response size: {:.1f} KB, "
                 "finished in {:.3f} sec",
                 request_id, processing_time, (double)ret.size() / 1000, timer.elapsed());
}

// ---------------------------------------------------------------- multi-graph servers

bool LoadReservations::reserve(size_t bytes, const std::function<bool()> &gone,
                               uint64_t poll_ms) {
    assert(bytes <= capacity_);
    std::unique_lock<std::mutex> lock(mutex_);
    if (!gone) {
        freed_.wait(lock, [&]() { return left_ >= bytes; });
    } else {
        while (left_ < bytes) {
            // asked without the lock: a socket's peek, a flag
            lock.unlock();
            const bool abandoned = gone();
            lock.lock();
            if (abandoned)
                return false;
            if (left_ >= bytes)
                break;
            freed_.wait_for(lock, std::chrono::milliseconds(poll_ms));
        }
    }
    left_ -= bytes;
    return true;
}

void LoadReservations::release(size_t bytes) {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        left_ += bytes;
        assert(left_ <= capacity_);
    }
    freed_.notify_all();
}

size_t LoadReservations::left() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return left_;
}

InRamPlan in_ram_plan(bool in_ram, bool server_on_mmap, size_t bytes, size_t capacity) {
    if (!in_ram)
        return InRamPlan::RESIDENT;
    // a server that loaded its indexes into RAM has them there already
    if (!server_on_mmap)
        return InRamPlan::RESIDENT_IN_RAM;
    if (bytes > capacity)
        return InRamPlan::RESIDENT_TOO_LARGE;
    return InRamPlan::LOAD;
}

const char* to_string(InRamPlan plan) {
    switch (plan) {
        case InRamPlan::RESIDENT: return "resident";
        case InRamPlan::RESIDENT_IN_RAM: return "resident_in_ram";
        case InRamPlan::RESIDENT_TOO_LARGE: return "resident_too_large";
        case InRamPlan::LOAD: return "load";
    }
    return "";
}

std::optional<bool> in_ram_field(const Json::Value &request) {
    if (!request.isObject() || !request.isMember("in_ram"))
        return std::nullopt;
    if (!request["in_ram"].isBool())
        throw std::invalid_argument("request.in_ram: expected a boolean");
    return request["in_ram"].asBool();
}

GraphSelection traverse_graph_selection(const Json::Value &request, bool multi_graph) {
    GraphSelection selection;
    if (!request.isObject())
        return selection;
    if (!multi_graph) {
        if (request.isMember("graph") || request.isMember("graph_path")) {
            throw std::invalid_argument("Bad request: this server hosts a single graph; "
                                        "remove the 'graph' / 'graph_path' field");
        }
        if (request.isMember("graphs")) {
            throw std::invalid_argument("Bad request: this server hosts a single graph; "
                                        "remove the 'graphs' field");
        }
        return selection;
    }
    if (request.isMember("graphs")) {
        // /search's field, one name: a traversal reads one graph
        const Json::Value &graphs = request["graphs"];
        if (request.isMember("graph")) {
            throw std::invalid_argument("Bad request: give the graph as 'graph' or as "
                                        "'graphs', not both");
        }
        if (!graphs.isArray() || graphs.size() != 1 || !graphs[0].isString()) {
            throw std::invalid_argument("Bad request: 'graphs' names the one graph a traversal "
                                        "reads: expected [name]");
        }
        selection.name = graphs[0].asString();
    } else {
        if (!request.isMember("graph") || !request["graph"].isString())
            throw std::invalid_argument("Bad request: 'graph' (index name) is required in "
                                        "multi-graph mode");
        selection.name = request["graph"].asString();
    }
    if (request.isMember("graph_path") && request["graph_path"].isString())
        selection.graph_path = request["graph_path"].asString();
    return selection;
}

std::vector<std::string> pattern_graph_names(const Json::Value &request,
                                             const std::vector<std::string> &known,
                                             size_t max_without_graphs) {
    std::vector<std::string> names;
    if (request.isObject() && request.isMember("graphs")) {
        const Json::Value &graphs = request["graphs"];
        const std::string expected = "request.graphs: expected a non-empty array of graph names "
                                     "(GET /capabilities lists them)";
        if (!graphs.isArray() || graphs.empty())
            throw std::invalid_argument(expected);
        for (const Json::Value &name : graphs) {
            if (!name.isString())
                throw std::invalid_argument(expected);
            names.push_back(name.asString());
        }
        for (const std::string &name : names) {
            if (std::find(known.begin(), known.end(), name) == known.end()) {
                throw std::invalid_argument("request.graphs: unknown graph '" + name
                                            + "' (GET /capabilities lists the graphs)");
            }
        }
    } else {
        if (known.size() > max_without_graphs) {
            throw std::invalid_argument("request.graphs: required on this server, which hosts "
                                        + std::to_string(known.size()) + " graph names (more "
                                        "than " + std::to_string(max_without_graphs)
                                        + "; GET /capabilities lists them)");
        }
        names = known;
    }
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    return names;
}

ColumnOverlap column_overlap(const std::vector<std::vector<std::string>> &columns,
                             const std::function<size_t(const std::string&)> &hash_of) {
    // the names by their hash, every collision resolved on the names themselves: 16 bytes per
    // column beside the names the annotations hold anyway
    struct Entry {
        size_t hash;
        uint32_t pair;
        uint32_t column;
    };
    ColumnOverlap result;
    std::vector<Entry> entries;
    for (const auto &names : columns) {
        result.columns += names.size();
    }
    entries.reserve(result.columns);
    const std::function<size_t(const std::string&)> hash
            = hash_of ? hash_of : std::function<size_t(const std::string&)>(std::hash<std::string>());
    for (size_t p = 0; p < columns.size(); ++p) {
        for (size_t c = 0; c < columns[p].size(); ++c) {
            entries.push_back({ hash(columns[p][c]), static_cast<uint32_t>(p),
                                static_cast<uint32_t>(c) });
        }
    }
    std::sort(entries.begin(), entries.end(), [](const Entry &a, const Entry &b) {
        return std::tie(a.hash, a.pair, a.column) < std::tie(b.hash, b.pair, b.column);
    });
    for (size_t i = 0; i < entries.size(); ) {
        size_t j = i + 1;
        while (j < entries.size() && entries[j].hash == entries[i].hash) {
            ++j;
        }
        if (entries[j - 1].pair != entries[i].pair) {
            // names of several pairs share this hash: which of them are one name
            std::map<std::string, std::vector<uint32_t>> pairs_of;
            for (size_t e = i; e < j; ++e) {
                pairs_of[columns[entries[e].pair][entries[e].column]].push_back(entries[e].pair);
            }
            for (auto &[name, pairs] : pairs_of) {
                std::sort(pairs.begin(), pairs.end());
                if (pairs.front() == pairs.back())
                    continue;
                ++result.shared;
                if (result.example.empty() || name < result.example)
                    result.example = name;
            }
        }
        i = j;
    }
    result.disjoint = !result.shared;
    return result;
}

} // namespace cli
} // namespace mtg
