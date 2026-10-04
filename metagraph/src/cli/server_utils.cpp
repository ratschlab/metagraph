#include <cctype>
#include <cerrno>
#include <filesystem>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
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


namespace mtg {
namespace cli {

using mtg::common::logger;

/**
 * Compress a STL string using zlib with given compression level and return
 * the binary data.
 * Source: https://panthema.net/2007/0328-ZLibString.html
 */
std::string compress_string(const std::string &str,
                            int compressionlevel = Z_BEST_COMPRESSION,
                            bool gzip = false,
                            const std::function<void()> &check = nullptr) {
    z_stream zs; // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    // zlib's deflate stream by default; a gzip container (window bits + 16) when the
    // client asked for gzip, which is what most HTTP clients send by default
    const int window_bits = 15 + (gzip ? 16 : 0);
    if (deflateInit2(&zs, compressionlevel, Z_DEFLATED, window_bits, 8,
                     Z_DEFAULT_STRATEGY) != Z_OK)
        throw std::runtime_error("deflateInit failed while compressing.");

    zs.next_in = (Bytef *)(str.data());
    zs.avail_in = str.size(); // set the z_stream's input

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
        zs.next_out = reinterpret_cast<Bytef *>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, Z_FINISH);

        if (outstring.size() < zs.total_out) {
            // append the block to the output string
            outstring.append(outbuffer, zs.total_out - outstring.size());
        }
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

// A streambuf appending to a string that calls |check| every |interval| bytes: the JSON text
// of a large response is written under the caller's check, byte for byte what an
// std::ostringstream receives (the same writer writes into either)
class CheckedStringBuf : public std::streambuf {
  public:
    CheckedStringBuf(std::string *out, const std::function<void()> &check, size_t interval)
          : out_(out), check_(check), interval_(interval), next_(interval) {}

  protected:
    int_type overflow(int_type c) override {
        if (traits_type::eq_int_type(c, traits_type::eof()))
            return traits_type::not_eof(c);
        out_->push_back(traits_type::to_char_type(c));
        tick();
        return c;
    }
    std::streamsize xsputn(const char_type *s, std::streamsize n) override {
        out_->append(s, static_cast<size_t>(n));
        tick();
        return n;
    }

  private:
    void tick() {
        if (out_->size() < next_)
            return;
        next_ = out_->size() + interval_;
        check_();
    }

    std::string *out_;
    const std::function<void()> &check_;
    size_t interval_;
    size_t next_;
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

} // namespace

// The content encoding to send, "" for none (RFC 9110 §12.5.3). Accept-Encoding is a
// comma-separated list of case-insensitive codings, each with an optional weight ";q="
// (1 when absent, 0 = not acceptable), `*` standing for every coding not listed by name;
// several header lines form one list, and an element with a malformed weight is ignored.
// A coding is acceptable at its own weight, else at the weight of `*`, else not at all.
// gzip wins over deflate at equal weight, the higher weight otherwise; the response is
// uncompressed when neither is acceptable or the client weights identity above both
// (uncompressed is also the fallback after "identity;q=0": refusing to answer is not what
// such a client asked for). A substring test sent gzip for "gzip;q=0, deflate;q=1" and
// for "identity, gzip;q=0" (review round 4, finding 4).
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

// Whether the TCP connection |fd| received the peer's FIN, or was reset, with data still
// waiting before the end: a peek returns that data, never the end behind it. The kernel's
// connection state says it (CLOSE_WAIT: the FIN was received; CLOSED: reset). Not a TCP socket,
// or a platform without the query: false (connected, as before)
static bool tcp_peer_finished(int fd) {
#if defined(__linux__)
    struct tcp_info info;
    socklen_t length = sizeof(info);
    if (::getsockopt(fd, IPPROTO_TCP, TCP_INFO, &info, &length) != 0)
        return false;
    return info.tcpi_state == TCP_CLOSE_WAIT || info.tcpi_state == TCP_CLOSE;
#elif defined(__APPLE__)
    struct tcp_connection_info info;
    socklen_t length = sizeof(info);
    if (::getsockopt(fd, IPPROTO_TCP, TCP_CONNECTION_INFO, &info, &length) != 0)
        return false;
    return info.tcpi_state == TCPS_CLOSE_WAIT || info.tcpi_state == TCPS_CLOSED;
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
        // state tells (review of the stage-4 backend, F3: a client that closed after a
        // trailing CRLF was walked to the end and answered with 68 MB into a dead socket)
        return tcp_peer_finished(fd);
    }
    return errno == ECONNRESET || errno == ENOTCONN || errno == EPIPE || errno == ETIMEDOUT;
}

std::string json_text(const Json::Value &value, bool compact, const std::function<void()> &check) {
    Json::StreamWriterBuilder builder;
    if (compact)
        builder["indentation"] = "";   // the traversal routes: half the bytes of the indented form
    if (!check)
        return Json::writeString(builder, value);
    std::string out;
    CheckedStringBuf buf(&out, check, size_t(1) << 16);
    std::ostream stream(&buf);
    // an exception of the check reaches the caller rather than setting badbit
    stream.exceptions(std::ios::badbit);
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    writer->write(value, &stream);
    return out;
}

std::string assemble_traverse_response(const Json::Value &envelope,
                                       const std::vector<std::string> &results,
                                       const std::function<void()> &check) {
    Json::Value before(Json::objectValue);
    Json::Value after(Json::objectValue);
    for (const std::string &name : envelope.getMemberNames()) {
        if (name == "results")
            throw std::logic_error("assemble_traverse_response: the envelope holds results");
        (name < "results" ? before : after)[name] = envelope[name];
    }
    std::string out = json_text(before, true, check);      // {...}
    out.pop_back();
    if (out.size() > 1)
        out += ',';
    out += "\"results\":[";
    size_t next_check = out.size() + (size_t(1) << 16);
    for (size_t i = 0; i < results.size(); ++i) {
        if (i)
            out += ',';
        out += results[i];
        if (check && out.size() >= next_check) {
            next_check = out.size() + (size_t(1) << 16);
            check();
        }
    }
    out += ']';
    const std::string tail = json_text(after, true, check);  // {...}
    if (tail.size() > 2) {
        out += ',';
        out.append(tail, 1, std::string::npos);
    } else {
        out += '}';
    }
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
                      const std::function<std::string(const GraphListEntry &)> &fingerprint) {
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
    // One index is one pair of files, however its paths are spelled: the lines are grouped by
    // the pair's real paths (the files themselves), so that two spellings of one pair agree too
    auto real = [](const std::string &path) {
        std::error_code ec;
        const std::filesystem::path p = std::filesystem::weakly_canonical(path, ec);
        return ec ? path : p.string();
    };
    std::map<Pair, Pair> real_of;                          // the pair as listed -> real paths
    for (const GraphListEntry &e : entries) {
        const Pair listed { e.graph_path, e.annotation_path };
        const Pair pair = real_of.emplace(listed, Pair { real(e.graph_path),
                                                         real(e.annotation_path) })
                                 .first->second;
        auto &[ns, fp] = stated[pair];
        agree(&ns, e.index_ns, e.line, "index_ns", listed);
        if (e.manifest_path.empty())
            continue;
        auto [it, fresh] = digests.try_emplace({ pair, e.manifest_path });
        if (fresh)
            it->second = fingerprint(e);
        agree(&fp, it->second, e.line, "index_fp (manifest digests)", listed);
    }
    // Two different indexes never state one fingerprint: the server checks sizes, not
    // contents, so one manifest whose files have the base names and sizes of another pair's
    // (two annotations of one graph with swapped memberships, written by one tool run) would
    // let a client take one index for the other (review of pass 5). Copies of one index under
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
    for (const auto &[listed, pair] : real_of) {
        const auto &s = stated.at(pair);
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

void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     size_t request_id,
                     const std::function<Json::Value(const std::string &)> &process,
                     bool compact,
                     const ResponseControl *control) {
    logger->info("[Server] {} request {} from {}", request->path, request_id,
                 request->remote_endpoint().address().to_string());
    Timer timer;
    // Retrieve string:
    std::string content = request->content.string();
    SimpleWeb::CaseInsensitiveMultimap header({ { "Content-Type", "application/json" } });
    SimpleWeb::StatusCode status;
    std::string ret;
    static const std::function<void()> kNoCheck;
    const std::function<void()> &check = control ? control->check : kNoCheck;

    try {
        // Return JSON string
        status = SimpleWeb::StatusCode::success_ok;
        if (control && control->write) {
            ret = control->write(process(content), check);
        } else {
            ret = json_text(process(content), compact, check);
        }
        const std::string encoding = requested_encoding(request);
        if (!encoding.empty()) {
            Timer compressing;
            const size_t text_bytes = ret.size();
            ret = compress_string(ret, control ? control->compression_level : Z_BEST_COMPRESSION,
                                  encoding == "gzip", check);
            if (control && control->on_compressed)
                control->on_compressed(text_bytes, compressing.elapsed());
            header.insert(std::make_pair("Content-Encoding", encoding));
            header.insert(std::make_pair("Content-Length", std::to_string(ret.size())));
        }
        // once more before the response is handed to the transport: a body shorter than the
        // check's interval is otherwise never checked after it was built
        if (check)
            check();
    } catch (const ClientGone &e) {
        // nobody to answer: nothing is written, and the connection is closed rather than
        // kept for a next request
        logger->info("[Server] Request {}: {}; no response written", request_id, e.what());
        response->close_connection_after_response = true;
        if (control && control->on_written)
            control->on_written(0, std::nullopt);
        return;
    } catch (const HttpError &e) {
        logger->warn("[Server] Error on request {} ({}): {}", request_id, e.status(), e.what());
        status = static_cast<SimpleWeb::StatusCode>(e.status());
        header = SimpleWeb::CaseInsensitiveMultimap({ { "Content-Type", "application/json" } });
        // the body as the route wrote it, uncompressed like every error
        ret = json_text(e.body(), compact);
    } catch (const CurrentlyInitializingError& e) {
        logger->info("[Server] Got a request during initialization. Asked to come back later");
        status = SimpleWeb::StatusCode::server_error_service_unavailable;
        header.insert(std::make_pair("Retry-After", "60")); // ask to come back in 60 seconds
        ret = json_str_with_error_msg("Server is currently initializing, please come back later.");
    } catch (const std::exception& e) {
        logger->warn("[Server] Error on request {}: {}", request_id, e.what());
        status = SimpleWeb::StatusCode::client_error_bad_request;
        ret = json_str_with_error_msg(e.what());
    } catch (...) {
        logger->warn("[Server] Error on request {}", request_id);
        status = SimpleWeb::StatusCode::server_error_internal_server_error;
        ret = json_str_with_error_msg("Internal server error");
    }
    double processing_time = timer.elapsed();
    response->write(status, ret, header);
    if (control && control->on_written)
        control->on_written(static_cast<int>(status), ret.size());
    logger->info("[Server] Request {} processing time: {:.3f} sec, response size: {:.1f} KB, "
                 "finished in {:.3f} sec",
                 request_id, processing_time, (double)ret.size() / 1000, timer.elapsed());
}

} // namespace cli
} // namespace mtg
