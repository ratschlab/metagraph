#include <algorithm>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <filesystem>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <netinet/in.h>
#include <netinet/tcp.h>
#if defined(__APPLE__)
#include <netinet/tcp_fsm.h>
#endif
#include <sys/socket.h>
#include <sys/stat.h>
#include <zlib.h>
#include <spdlog/fmt/fmt.h>
#include <json/json.h>
#include <server_http.hpp>

#include "annotation/coord_to_header.hpp"
#include "cli/json_helpers.hpp"
#include "common/logger.hpp"
#include "common/threads/threading.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/string_utils.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "server_envelope.hpp"
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

namespace {

#ifdef ASIO_STANDALONE
namespace net = ::asio;
#else
namespace net = boost::asio;
#endif

// What the drain of a refused body reads at once (refuse_body_too_large): its whole buffer,
// whatever the client sends
constexpr size_t kDiscardBufferBytes = 16 << 10;

// The drain of a refused body: its buffer, and how many more bytes it discards before it closes
struct Drain {
    std::vector<char> buffer;
    uint64_t left;
};

// Reads and discards from |connection| until the client closes it, the connection's timer shuts
// it (both end the pending read with an error) or |drain|'s bytes are spent; then closes it (the
// last reference, held by the pending read, frees it)
template <class Connection>
void discard_until_closed(std::shared_ptr<Connection> connection, std::shared_ptr<Drain> drain) {
    // the socket and the buffer's view taken before the handler is made: it takes both pointers
    // over, and the order in which a call's arguments are evaluated is unspecified (GCC makes
    // the handler first, which left a dereference of the moved-from buffer)
    auto *socket = connection->socket.get();
    const auto view = net::buffer(drain->buffer);
    socket->async_read_some(view,
            [connection = std::move(connection), drain = std::move(drain)](
                    const SimpleWeb::error_code &ec, size_t bytes) mutable {
        auto lock = connection->handler_runner->continue_lock();
        if (!lock)
            return;
        if (ec || bytes >= drain->left) {
            connection->close();
            return;
        }
        drain->left -= bytes;
        discard_until_closed(std::move(connection), std::move(drain));
    });
}

} // namespace

bool refuse_body_too_large(const std::shared_ptr<HttpServer::Request> &request,
                           uint64_t max_mb) {
    auto connection = request_connection(*request);
    if (!connection)
        return false;
    std::optional<uint64_t> length;
    if (auto it = request->header.find("Content-Length"); it != request->header.end()) {
        // the library read it with std::stoull, which takes "-1" as 2^64 - 1 and "3145728junk"
        // as 3145728: a Content-Length that is no byte count (content_length_value) is
        // dropped as the library drops one it cannot read at all ("abc", 2^64)
        length = content_length_value(it->second);
        if (!length) {
            connection->close();
            return false;
        }
    }
    const std::string body = error_first_json_text(body_too_large_json(length, max_mb));
    auto text = std::make_shared<std::string>(
            "HTTP/1.1 413 Payload Too Large\r\nContent-Type: application/json\r\n"
            "Content-Length: " + std::to_string(body.size()) + "\r\nConnection: close\r\n\r\n"
            + body);
    // the whole of what follows (the answer and the drain) within the connection's own timer,
    // which shuts the connection when it ends (and so ends a pending read or write)
    connection->set_timeout(kBodyTooLargeDrainS);
    // what the drain reads at most: the larger of kBodyTooLargeDrainBytes and the cap
    auto drain = std::make_shared<Drain>(Drain {
        std::vector<char>(kDiscardBufferBytes),
        std::max<uint64_t>(kBodyTooLargeDrainBytes, max_mb << 20) });
    // (the socket and the text's view taken before the handler takes the connection over)
    auto *socket = connection->socket.get();
    const auto view = net::buffer(*text);
    net::async_write(*socket, view,
            [connection = std::move(connection), text, drain = std::move(drain)](
                    const SimpleWeb::error_code &ec, size_t) mutable {
        auto lock = connection->handler_runner->continue_lock();
        if (!lock)
            return;
        if (ec) {
            connection->close();
            return;
        }
        // the answer is complete (FIN after it); what the client still sends is read and
        // dropped, so that its connection is not reset while it sends
        SimpleWeb::error_code ignored;
        connection->socket->lowest_layer().shutdown(net::ip::tcp::socket::shutdown_send,
                                                    ignored);
        discard_until_closed(std::move(connection), std::move(drain));
    });
    return true;
}

// ---------------------------------------------------------------- the connection cap (H3)

namespace {

using IoContext = decltype(HttpServer::io_service)::element_type;

// the shortest time between two warnings of the accept loop pausing at the cap
constexpr auto kCapLogInterval = std::chrono::seconds(10);

} // namespace

struct CappedHttpServer::CapState {
    std::mutex mutex;
    size_t max = 0;
    // the connections open, the one awaiting an accept included
    size_t open = 0;
    // the accept loop stopped at the cap: the next close resumes it
    bool paused = false;
    std::chrono::steady_clock::time_point paused_at;
    size_t pauses = 0;
    // the warnings are throttled: the pauses since the last one, and whether the current
    // pause was warned of (its resume is then logged too)
    std::chrono::steady_clock::time_point warned_at;
    bool warned = false;
    size_t pauses_unwarned = 0;
    bool episode_warned = false;
    // what a close needs to resume the loop: the server's io_context (weak: a close while it
    // is destroyed posts nothing), its handler runner (stopped with the server: a resume
    // posted before then does nothing) and the server itself
    std::weak_ptr<IoContext> io;
    std::shared_ptr<SimpleWeb::ScopeRunner> runner;
    CappedHttpServer *server = nullptr;

    // a connection closed: resume the accept loop if it stopped at the cap
    void closed() {
        std::shared_ptr<IoContext> context;
        std::shared_ptr<SimpleWeb::ScopeRunner> run;
        CappedHttpServer *to_resume = nullptr;
        bool log_resume = false;
        double waited_ms = 0;
        {
            std::lock_guard<std::mutex> lock(mutex);
            --open;
            if (!paused || open >= max)
                return;
            paused = false;
            context = io.lock();
            run = runner;
            to_resume = server;
            log_resume = episode_warned;
            waited_ms = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - paused_at).count();
        }
        if (!context || !run || !to_resume)
            return;
        SimpleWeb::post(*context, [run, to_resume, log_resume, waited_ms]() {
            auto lock = run->continue_lock();
            if (!lock)
                return;
            if (log_resume) {
                logger->info("[Server] Accepting connections again after {:.0f} ms at the "
                             "connection cap", waited_ms);
            }
            to_resume->accept();
        });
    }
};

CappedHttpServer::CappedHttpServer() : cap_(std::make_shared<CapState>()) {}

void CappedHttpServer::set_max_connections(size_t max_connections) {
    std::lock_guard<std::mutex> lock(cap_->mutex);
    cap_->max = max_connections;
}

size_t CappedHttpServer::max_connections() const {
    std::lock_guard<std::mutex> lock(cap_->mutex);
    return cap_->max;
}

size_t CappedHttpServer::open_connections() const {
    std::lock_guard<std::mutex> lock(cap_->mutex);
    return cap_->open;
}

size_t CappedHttpServer::pauses() const {
    std::lock_guard<std::mutex> lock(cap_->mutex);
    return cap_->pauses;
}

void CappedHttpServer::accept() {
    // stopped (the acceptor closed): no accept, as the library's own loop ends at its
    // operation_aborted
    if (!acceptor || !acceptor->is_open())
        return;
    {
        std::lock_guard<std::mutex> lock(cap_->mutex);
        cap_->io = io_service;
        cap_->runner = handler_runner;
        cap_->server = this;
        if (cap_->max && cap_->open >= cap_->max) {
            // at the cap: nothing awaits an accept until a connection closes (CapState::
            // closed); the kernel's listen queue holds the connections meanwhile
            const auto now = std::chrono::steady_clock::now();
            cap_->paused = true;
            cap_->paused_at = now;
            ++cap_->pauses;
            cap_->episode_warned = !cap_->warned || now - cap_->warned_at >= kCapLogInterval;
            if (cap_->episode_warned) {
                logger->warn("[Server] {} connections open (--max-connections {}): new "
                             "connections wait in the listen queue until one closes{}",
                             cap_->open, cap_->max,
                             cap_->pauses_unwarned
                                 ? fmt::format(" ({} more pauses at the cap since the last "
                                               "such line)", cap_->pauses_unwarned)
                                 : std::string());
                cap_->warned = true;
                cap_->warned_at = now;
                cap_->pauses_unwarned = 0;
            } else {
                ++cap_->pauses_unwarned;
            }
            return;
        }
        ++cap_->open;
    }
    // as the library's create_connection, with a deleter that also counts the close
    auto connections = this->connections;
    std::shared_ptr<CapState> cap = cap_;
    std::shared_ptr<Connection> connection(
            new Connection(handler_runner, *io_service),
            [connections, cap](Connection *c) {
                {
                    SimpleWeb::LockGuard lock(connections->mutex);
                    auto it = connections->set.find(c);
                    if (it != connections->set.end())
                        connections->set.erase(it);
                }
                delete c;
                cap->closed();
            });
    {
        SimpleWeb::LockGuard lock(connections->mutex);
        connections->set.emplace(connection.get());
    }
    // as the library's Server<HTTP>::accept
    acceptor->async_accept(*connection->socket,
                           [this, connection](const SimpleWeb::error_code &ec) {
        auto lock = connection->handler_runner->continue_lock();
        if (!lock)
            return;
        // immediately accept the next one (unless the server stops), or stop at the cap
        if (ec != SimpleWeb::error::operation_aborted)
            this->accept();
        auto session = std::make_shared<Session>(config.max_request_streambuf_size, connection);
        if (!ec) {
            net::ip::tcp::no_delay option(true);
            SimpleWeb::error_code ignored;
            session->connection->socket->set_option(option, ignored);
            this->read(session);
        } else if (this->on_error) {
            this->on_error(session->request, ec);
        }
    });
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

namespace {

// The compact text of |object| cut at its member |key|, which it does not hold: |*head| gets
// `{`, the members before |key| (jsoncpp writes an object's members in byte order of their
// names) and their comma when there are any, then `"key":`; |*tail| the members after |key|,
// each after a comma, and the closing `}`. Together around |key|'s value they are byte for
// byte json_text of |object| holding it. |check| and |max_gap_ms| as json_text's
void split_object_text(const Json::Value &object, const std::string &key,
                       const std::function<void()> &check, double *max_gap_ms,
                       std::string *head, std::string *tail) {
    Json::Value before(Json::objectValue);
    Json::Value after(Json::objectValue);
    for (const std::string &name : object.getMemberNames()) {
        if (name == key) {
            throw std::logic_error("assemble_traverse_response: the envelope holds " + key
                                   + " already");
        }
        (name < key ? before : after)[name] = object[name];
    }
    *head = json_text(before, true, check, max_gap_ms);         // {...}
    const std::string after_text = json_text(after, true, check, max_gap_ms);  // {...}
    head->pop_back();
    if (head->size() > 1)
        *head += ',';
    *head += '"' + key + "\":";
    // `{}` without members after |key|, `{...}` with them: the first `{` becomes a comma
    *tail = after_text.size() > 2 ? ',' + after_text.substr(1) : std::string("}");
}

/**
 * |prefix|, the texts of |results| separated by commas, then |suffix|: the text of a /traverse
 * response whose results (each seed's, written once built) are an array's elements between
 * the two. The exact size is reserved at once, so the text never grows by doubling, whose last
 * step held an old and a new buffer of up to its size beside the texts; and each text is freed
 * once copied: what was copied and what is left to copy are the response once, where texts
 * kept alive beside it until the handler returned would stay through the compression and the
 * transport's copy (2.5x the text held with gzip, 3.5x without). |check| and |max_gap_ms| as
 * json_text's, the texts copied in pieces up to the next check
 */
std::string splice_results(std::string prefix, std::vector<std::string> results,
                           const std::string &suffix, const std::function<void()> &check,
                           double *max_gap_ms) {
    size_t size = prefix.size() + suffix.size();
    for (size_t i = 0; i < results.size(); ++i) {
        size += results[i].size() + (i > 0);
    }
    std::string out;
    out.reserve(size);
    out += prefix;
    std::string().swap(prefix);
    if (!check) {
        for (size_t i = 0; i < results.size(); ++i) {
            if (i)
                out += ',';
            out += results[i];
            std::string().swap(results[i]);
        }
    } else {
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
    out += suffix;
    // the size computed is the size written: a mismatch would only cost a reallocation, never
    // a byte, so it is asserted where tests run unoptimised
    assert(out.size() == size);
    return out;
}

} // namespace

std::string assemble_traverse_envelope(const Json::Value &envelope,
                                       std::vector<std::string> results,
                                       const std::function<void()> &check,
                                       double *max_gap_ms) {
    const Json::Value &answers = envelope["answers"];
    if (!answers.isArray() || answers.size() != 1 || !answers[0].isObject()) {
        throw std::logic_error("assemble_traverse_envelope: the envelope has no single entry "
                               "to hold the results");
    }
    // the envelope around its one entry, then the entry around its results: the members after
    // either (the usage, the timing, a few KiB) are written before a byte of the results is
    // copied, so that the whole response's size is known first
    Json::Value outer(Json::objectValue);
    for (const std::string &name : envelope.getMemberNames()) {
        if (name != "answers")
            outer[name] = envelope[name];
    }
    std::string outer_head, outer_tail, entry_head, entry_tail;
    split_object_text(outer, "answers", check, max_gap_ms, &outer_head, &outer_tail);
    split_object_text(answers[0], "results", check, max_gap_ms, &entry_head, &entry_tail);
    // {...,"answers":[{...,"results":[ the texts ]...}]...}
    return splice_results(outer_head + '[' + entry_head + '[', std::move(results),
                          ']' + entry_tail + ']' + outer_tail, check, max_gap_ms);
}

std::string assemble_traverse_response(const Json::Value &envelope,
                                       std::vector<std::string> results,
                                       const std::function<void()> &check,
                                       double *max_gap_ms) {
    // the envelope's members after "results" are written first (the usage, the timing: a few
    // KiB), so that the whole response's size is known before a byte of the results is copied
    std::string head, tail;
    split_object_text(envelope, "results", check, max_gap_ms, &head, &tail);
    // {...,"results":[ the texts ]...}
    return splice_results(head + '[', std::move(results), ']' + tail, check, max_gap_ms);
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
    // made absolute first: weakly_canonical of a relative path whose first element does not
    // exist is left relative by libstdc++ ("g1"), while "./g1" is made absolute, so the two
    // spellings of one missing file would be two indexes on Linux (libc++ makes both absolute)
    auto real = [](const std::string &path) {
        std::error_code ec;
        const std::filesystem::path abs = std::filesystem::absolute(path, ec);
        if (ec)
            return path;
        const std::filesystem::path p = std::filesystem::weakly_canonical(abs, ec);
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

RequestAnswer busy_answer(const ServerBusy &busy, bool compact) {
    RequestAnswer answer;
    answer.status = 503;
    answer.header.emplace_back("Retry-After", std::to_string(kBusyRetryAfterS));
    Json::Value body;
    body["code"] = "busy";
    body["error"] = busy.what();
    body["server_instance"] = busy.server_instance();
    answer.body = json_text(body, compact);
    return answer;
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
            // in pieces of kCompressPieceBytes, so that a stop (an attempt at its bound, a
            // client gone) is seen at least every 4 MiB of text: the check per 32 KiB block of
            // output alone runs once per tens of MiB on a text that compresses strongly. The
            // compressed bytes are those of one call over the whole text
            // (GraphletServer.CompressionPiecesOfTheServerKeepTheBytes)
            ret = compress_string(ret, control ? control->compression_level : Z_BEST_COMPRESSION,
                                  encoding == "gzip", check, kCompressPieceBytes);
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
        ret = e.text().empty() ? json_text(e.body(), compact) : e.text();
    } catch (const CurrentlyInitializingError& e) {
        logger->info("[Server] Got a request during initialization. Asked to come back later");
        answer.status = 503;
        answer.header.emplace_back("Retry-After", "60"); // ask to come back in 60 seconds
        ret = json_str_with_error_msg("Server is currently initializing, please come back later.");
    } catch (const ServerBusy &e) {
        // a route at its cap (E19): the request was not parsed, nothing of it registered (the
        // handlers refuse it before process_request; this is the case of a slot taken after the
        // index finished loading)
        logger->info("[Server] Request {}: busy; 503", request_id);
        answer = busy_answer(e, compact);
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
                     const ResponseControl *control,
                     std::function<void(const SimpleWeb::error_code &)> on_sent) {
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
    if (on_sent)
        response->send(std::move(on_sent));
    logger->info("[Server] Request {} processing time: {:.3f} sec, response size: {:.1f} KB, "
                 "finished in {:.3f} sec",
                 request_id, processing_time, (double)ret.size() / 1000, timer.elapsed());
}

void write_answer(std::shared_ptr<HttpServer::Response> &response, const RequestAnswer &answer) {
    // as process_request writes an answer: the Content-Type first, then the answer's fields
    SimpleWeb::CaseInsensitiveMultimap header({ { "Content-Type", "application/json" } });
    for (const auto &field : answer.header) {
        header.insert(field);
    }
    response->write(static_cast<SimpleWeb::StatusCode>(answer.status), answer.body, header);
}

void process_serialized_request(std::shared_ptr<HttpServer::Response> &response,
                                const std::shared_ptr<HttpServer::Request> &request,
                                size_t request_id,
                                const std::function<std::string(const std::string &)> &process,
                                std::function<void(const SimpleWeb::error_code &)> on_sent) {
    std::string body;
    ResponseControl control;
    control.write = [&](const Json::Value &, const std::function<void()> &) {
        return std::move(body);
    };
    process_request(response, request, request_id, [&](const std::string &content) {
        body = process(content);
        return Json::Value();
    }, false, &control, std::move(on_sent));
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
        // `pairs` (release 10) is no field of a single-graph request: the parser refuses it as
        // an unknown field, as before release 10
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
        }
        // each name checked once, in a set of the known names: the check runs on the HTTP
        // thread before any pair is queued, and a list may repeat a name many times over
        // (the body is unlimited by default). The refusal names the first unknown name in
        // the list's order
        const std::unordered_set<std::string> known_set(known.begin(), known.end());
        std::unordered_set<std::string> seen;
        for (const Json::Value &name : graphs) {
            std::string value = name.asString();
            if (!seen.insert(value).second)
                continue;
            if (!known_set.count(value)) {
                throw std::invalid_argument("request.graphs: unknown graph '"
                                            + ::utils::clip_for_message(value)
                                            + "' (GET /capabilities lists the graphs)");
            }
            names.push_back(std::move(value));
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

Json::Value max_graphs_without_selection_json(bool multi_graph, uint64_t threshold) {
    return multi_graph ? Json::Value(Json::UInt64(threshold)) : Json::Value();
}

Json::Value max_request_body_mb_json(uint64_t mb) {
    return mb ? Json::Value(Json::UInt64(mb)) : Json::Value();
}

std::optional<uint64_t> content_length_value(const std::string &text) {
    std::optional<uint64_t> count;
    size_t from = 0;
    while (true) {
        size_t to = text.find(',', from);
        if (to == std::string::npos)
            to = text.size();
        size_t a = from, b = to;
        while (a < b && (text[a] == ' ' || text[a] == '\t'))
            ++a;
        while (b > a && (text[b - 1] == ' ' || text[b - 1] == '\t'))
            --b;
        // leading zeros are the same count; at most 20 digits remain below 2^64
        while (a + 1 < b && text[a] == '0')
            ++a;
        if (a == b || b - a > 20 || !std::all_of(text.begin() + a, text.begin() + b,
                                                 [](unsigned char c) { return std::isdigit(c); }))
            return std::nullopt;
        uint64_t value;
        try {
            value = std::stoull(text.substr(a, b - a));
        } catch (const std::exception &) {
            return std::nullopt;     // 2^64 and more
        }
        if (count && *count != value)
            return std::nullopt;     // a list of differing counts
        count = value;
        if (to == text.size())
            return count;
        from = to + 1;
    }
}

Json::Value body_too_large_json(std::optional<uint64_t> content_length, uint64_t max_mb) {
    Json::Value body(Json::objectValue);
    body["code"] = "body_too_large";
    const std::string cap = std::to_string(max_mb) + " MiB (--max-request-body-mb)";
    body["error"] = content_length
        ? "the request body (" + std::to_string(*content_length) + " bytes) is larger than "
          "this server reads, " + cap + ": it was not read and nothing of it ran"
        : "the request body (chunked) grew past what this server reads, " + cap
          + ": the rest was not read and nothing of it ran";
    body["max_request_body_mb"] = Json::UInt64(max_mb);
    return body;
}

Json::Value too_many_json_values_json(uint64_t max_values) {
    Json::Value body(Json::objectValue);
    body["code"] = "too_many_json_values";
    body["error"] = TooManyJsonValues(max_values).what();
    body["max_request_json_values"] = Json::UInt64(max_values);
    return body;
}

std::string error_first_json_text(const Json::Value &body) {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    std::string text = "{\"error\":" + Json::writeString(builder, body["error"]);
    for (const std::string &name : body.getMemberNames()) {
        if (name != "error") {
            text += "," + Json::writeString(builder, Json::Value(name)) + ":"
                    + Json::writeString(builder, body[name]);
        }
    }
    return text + "}";
}

Json::Value request_intake_json(uint64_t max_body_mb, uint64_t max_json_values,
                                uint64_t max_connections) {
    auto refusal = [](int status, const char *code, const char *cap) {
        Json::Value r(Json::objectValue);
        r["status"] = status;
        r["code"] = code;
        Json::Value members(Json::arrayValue);
        for (const char *member : { "code", "error", cap }) {
            members.append(member);
        }
        r["members"] = std::move(members);
        return r;
    };
    Json::Value intake(Json::objectValue);
    intake["max_request_body_mb"] = max_request_body_mb_json(max_body_mb);
    intake["body_too_large"] = refusal(413, "body_too_large", "max_request_body_mb");
    intake["max_request_json_values"]
            = max_json_values ? Json::Value(Json::UInt64(max_json_values)) : Json::Value();
    Json::Value routes(Json::arrayValue);
    for (const char *route : kJsonValuesRoutes) {
        routes.append(route);
    }
    intake["json_values_routes"] = std::move(routes);
    intake["too_many_json_values"]
            = refusal(400, "too_many_json_values", "max_request_json_values");
    // H3: the connections held open at once (a further one waits in the listen queue)
    intake["max_connections"]
            = max_connections ? Json::Value(Json::UInt64(max_connections)) : Json::Value();
    return intake;
}

// ---------------------------------------------------------------- GET /column_labels (H3)

ColumnLabelList::ColumnLabelList(std::vector<const std::vector<std::string> *> parts)
      : parts_(std::move(parts)) {
    starts_.reserve(parts_.size());
    for (const auto *part : parts_) {
        starts_.push_back(total_);
        total_ += part->size();
    }
}

namespace {

// What jsoncpp quotes as it is (isAnyCharRequiredQuoting): printable ASCII but `"` and `\`
bool plain_json_string(const std::string &label) {
    for (char c : label) {
        const auto u = static_cast<unsigned char>(c);
        if (u < 0x20 || u >= 0x80 || c == '"' || c == '\\')
            return false;
    }
    return true;
}

// The code point jsoncpp's writer reads at |*s| (its utf8ToCodepoint, byte for byte: a lead byte
// is taken with as many following bytes as it announces, whatever they are, and anything it
// cannot read is U+FFFD), |*s| left on the code point's last byte
unsigned int jsoncpp_codepoint(const char *&s, const char *e) {
    constexpr unsigned int kReplacement = 0xFFFD;
    const unsigned int b0 = static_cast<unsigned char>(*s);
    auto cont = [&](int i) { return static_cast<unsigned int>(s[i]) & 0x3F; };
    if (b0 < 0x80)
        return b0;
    if (b0 < 0xE0) {
        if (e - s < 2)
            return kReplacement;
        const unsigned int c = ((b0 & 0x1F) << 6) | cont(1);
        s += 1;
        return c < 0x80 ? kReplacement : c;
    }
    if (b0 < 0xF0) {
        if (e - s < 3)
            return kReplacement;
        const unsigned int c = ((b0 & 0x0F) << 12) | (cont(1) << 6) | cont(2);
        s += 2;
        if (c >= 0xD800 && c <= 0xDFFF)
            return kReplacement;
        return c < 0x800 ? kReplacement : c;
    }
    if (b0 < 0xF8) {
        if (e - s < 4)
            return kReplacement;
        const unsigned int c = ((b0 & 0x07) << 18) | (cont(1) << 12) | (cont(2) << 6) | cont(3);
        s += 3;
        return c < 0x10000 ? kReplacement : c;
    }
    return kReplacement;
}

// The JSON string jsoncpp writes for |label| (valueToQuotedStringN without emitUTF8, as
// json_text and Json::writeString write it), passed piece by piece to |emit|(data, size):
// quoted, `"` `\` \b \f \n \r \t as two characters, every other code point outside
// [0x20, 0x7F] as \uxxxx (a surrogate pair beyond the BMP). Nothing is built: |emit| counts
// the length or appends the text
template <class Emit>
void jsoncpp_json_string(const std::string &label, Emit &&emit) {
    emit("\"", 1);
    if (plain_json_string(label)) {
        emit(label.data(), label.size());
        emit("\"", 1);
        return;
    }
    static const char kHex[] = "0123456789abcdef";
    auto u16 = [&](unsigned int x) {
        const char t[6] = { '\\', 'u', kHex[(x >> 12) & 0xF], kHex[(x >> 8) & 0xF],
                            kHex[(x >> 4) & 0xF], kHex[x & 0xF] };
        emit(t, 6);
    };
    const char *const end = label.data() + label.size();
    for (const char *c = label.data(); c != end; ++c) {
        switch (*c) {
            case '"': emit("\\\"", 2); break;
            case '\\': emit("\\\\", 2); break;
            case '\b': emit("\\b", 2); break;
            case '\f': emit("\\f", 2); break;
            case '\n': emit("\\n", 2); break;
            case '\r': emit("\\r", 2); break;
            case '\t': emit("\\t", 2); break;
            default: {
                unsigned int cp = jsoncpp_codepoint(c, end);
                if (cp >= 0x20 && cp <= 0x7F) {
                    const char ch = static_cast<char>(cp);
                    emit(&ch, 1);
                } else if (cp < 0x10000) {
                    u16(cp);
                } else {
                    cp -= 0x10000;
                    u16((cp >> 10) + 0xD800);
                    u16((cp & 0x3FF) + 0xDC00);
                }
            }
        }
    }
    emit("\"", 1);
}

} // namespace

size_t json_string_length(const std::string &label) {
    if (plain_json_string(label))
        return label.size() + 2;
    size_t n = 0;
    jsoncpp_json_string(label, [&n](const char *, size_t k) { n += k; });
    return n;
}

size_t append_json_string(std::string *out, const std::string &label) {
    const size_t before = out->size();
    jsoncpp_json_string(label, [out](const char *data, size_t k) { out->append(data, k); });
    return out->size() - before;
}

std::optional<std::string> column_labels_text(const ColumnLabelList &list, uint64_t max_bytes) {
    // GPT review 6 #13: the labels' JSON strings are measured first, nothing built, so that a
    // refusal allocates nothing and an answer is built at its final size
    uint64_t bytes = 0;
    bool over = false;
    list.visit_from(0, [&](const std::string &label) {
        bytes += json_string_length(label);
        over = max_bytes && bytes > max_bytes;
        return !over;
    });
    if (over)
        return std::nullopt;
    // json_text's styled array: "[", then each element on its own line after a tab, the
    // elements separated by ",", then "\n]"; an empty array "[]"
    std::string out;
    out.reserve(bytes + 3 * list.total() + 2);
    out.push_back('[');
    bool first = true;
    list.visit_from(0, [&](const std::string &label) {
        out.append(first ? "\n\t" : ",\n\t");
        first = false;
        append_json_string(&out, label);
        return true;
    });
    out.append(list.total() ? "\n]" : "]");
    return out;
}

namespace {

// a parameter of GET /column_labels: decimal digits only, in [|min|, 2^53 - 1]
uint64_t column_labels_integer(const std::string &key, const std::string &value, uint64_t min) {
    constexpr uint64_t kMax = (uint64_t(1) << 53) - 1;
    uint64_t v = 0;
    bool ok = !value.empty() && value.size() <= 16;
    for (char c : value) {
        if (!ok || c < '0' || c > '9') {
            ok = false;
            break;
        }
        v = v * 10 + static_cast<uint64_t>(c - '0');
    }
    if (!ok || v < min || v > kMax) {
        // the value quoted as every error text quotes a client's string (clipped past 1 KiB)
        throw std::invalid_argument("Bad request: GET /column_labels: '" + key
                                    + "' must be an integer in [" + std::to_string(min) + ", "
                                    + std::to_string(kMax) + "], got '"
                                    + ::utils::clip_for_message(value) + "'");
    }
    return v;
}

} // namespace

ColumnLabelsPage column_labels_parameters(
        const std::vector<std::pair<std::string, std::string>> &query) {
    ColumnLabelsPage page;
    bool offset_given = false;
    bool limit_given = false;
    for (const auto &[key, value] : query) {
        bool *given = key == "offset" ? &offset_given : key == "limit" ? &limit_given : nullptr;
        if (!given) {
            throw std::invalid_argument("Bad request: unknown parameter '"
                                        + ::utils::clip_for_message(key)
                                        + "' of GET /column_labels (offset, limit)");
        }
        if (*given)
            throw std::invalid_argument("Bad request: the parameter '" + key + "' is given twice");
        *given = true;
        if (key == "offset") {
            page.offset = column_labels_integer(key, value, 0);
        } else {
            page.limit = column_labels_integer(key, value, 1);
        }
    }
    page.paged = offset_given || limit_given;
    return page;
}

std::string column_labels_page_text(const ColumnLabelList &list, const ColumnLabelsPage &page,
                                    uint64_t max_bytes) {
    // compact, the members in byte order, as json_text(compact) writes an object
    std::string out = "{\"labels\":[";
    uint64_t bytes = 0;
    uint64_t count = 0;
    list.visit_from(page.offset, [&](const std::string &label) {
        if (page.limit && count >= *page.limit)
            return false;
        // measured before it is written (GPT review 6 #13: a label that does not fit is never
        // built); at least one label a page, so that paging moves on past a label longer than
        // the cap
        const size_t n = json_string_length(label);
        if (count && max_bytes && bytes + n > max_bytes)
            return false;
        if (count)
            out.push_back(',');
        append_json_string(&out, label);
        bytes += n;
        ++count;
        return true;
    });
    const uint64_t next = page.offset + count;
    out.append("],\"next_offset\":");
    out.append(count && next < list.total() ? std::to_string(next) : "null");
    out.append(",\"offset\":" + std::to_string(page.offset));
    out.append(",\"total\":" + std::to_string(list.total()) + "}");
    return out;
}

Json::Value column_labels_too_large_json(uint64_t total, uint64_t max_bytes) {
    Json::Value body(Json::objectValue);
    body["code"] = "column_labels_too_large";
    body["error"] = "the column labels (" + std::to_string(total) + " labels) take more than "
                    + std::to_string(max_bytes) + " bytes of JSON strings, more than one answer "
                    "holds (--max-column-labels-bytes): page through them with GET "
                    "/column_labels?offset=0&limit=<n>, following next_offset";
    body["max_column_labels_bytes"] = Json::UInt64(max_bytes);
    body["total"] = Json::UInt64(total);
    return body;
}

Json::Value column_labels_capabilities_json(uint64_t max_bytes, std::optional<uint64_t> total) {
    Json::Value c(Json::objectValue);
    c["max_column_labels_bytes"] = max_bytes ? Json::Value(Json::UInt64(max_bytes))
                                             : Json::Value();
    Json::Value paging(Json::arrayValue);
    paging.append("limit");
    paging.append("offset");
    c["paging"] = std::move(paging);
    Json::Value refusal(Json::objectValue);
    refusal["status"] = 400;
    refusal["code"] = "column_labels_too_large";
    Json::Value members(Json::arrayValue);
    for (const char *member : { "code", "error", "max_column_labels_bytes", "total" }) {
        members.append(member);
    }
    refusal["members"] = std::move(members);
    c["too_large"] = std::move(refusal);
    c["total"] = total ? Json::Value(Json::UInt64(*total)) : Json::Value();
    return c;
}

// The entries of /pattern's envelope are the shared envelope's (server_envelope.hpp), which
// tags a pair and its refusal the same way on every route
Json::Value pattern_pair_answered(const std::string &name,
                                  const std::string &graph_path,
                                  const std::string &annotation_path,
                                  const std::string &index_fp,
                                  Json::Value answer) {
    return envelope_entry_answered({ name, { graph_path, annotation_path } }, index_fp,
                                   std::move(answer));
}

Json::Value pattern_pair_refused(const std::string &name,
                                 const std::string &graph_path,
                                 const std::string &annotation_path,
                                 const std::string &index_fp,
                                 int http_status,
                                 const Json::Value &body) {
    return envelope_entry_refused({ name, { graph_path, annotation_path } }, index_fp,
                                  http_status, body);
}

Json::Value pattern_envelope(const std::vector<std::string> &names,
                             std::vector<Json::Value> entries,
                             double elapsed_ms) {
    Json::Value out;
    out["pattern_contract_version"] = graph::pattern::kPatternContractVersion;
    Json::Value graphs(Json::arrayValue);
    for (const std::string &name : names) {
        graphs.append(name);
    }
    out["graphs"] = std::move(graphs);
    Json::UInt answered = 0;
    Json::UInt refused = 0;
    Json::Value answers(Json::arrayValue);
    for (Json::Value &entry : entries) {
        (entry["outcome"].asString() == "refused" ? refused : answered)++;
        answers.append(std::move(entry));
    }
    out["answered"] = answered;
    out["refused"] = refused;
    out["answers"] = std::move(answers);
    out["timing"]["elapsed_ms"] = elapsed_ms;
    return out;
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

Json::Value headers_unique_json(const std::vector<HeaderPair> &pairs, size_t num_threads,
                                size_t max_examples, size_t max_listed) {
    std::vector<const annot::CoordToHeader*> indexes;
    bool every_pair = true;
    for (const HeaderPair &pair : pairs) {
        indexes.push_back(pair.coord_to_header);
        every_pair &= pair.coord_to_header != nullptr;
    }
    // the byte-smallest repeated headers, so that the examples do not depend on the hash
    std::map<std::string, Json::Value> examples;
    const uint64_t collisions = annot::CoordToHeader::find_repeated_headers(indexes,
            [&](std::string_view header, uint64_t count,
                const std::vector<annot::CoordToHeader::HeaderOccurrence> &listed) {
                if (examples.size() == max_examples
                        && (!max_examples || header >= examples.rbegin()->first)) {
                    return;
                }
                Json::Value example;
                example["header"] = std::string(header);
                example["occurrences"] = Json::Value(static_cast<Json::UInt64>(count));
                Json::Value occurrences(Json::arrayValue);
                for (const auto &o : listed) {
                    const HeaderPair &pair = pairs[o.index];
                    Json::Value p;
                    p["graph"] = pair.graph;
                    p["graph_path"] = pair.graph_path;
                    p["annotation_path"] = pair.annotation_path;
                    p["column"] = pair.columns && o.column < pair.columns->size()
                            ? Json::Value((*pair.columns)[o.column])
                            : Json::Value();
                    occurrences.append(std::move(p));
                }
                example["pairs"] = std::move(occurrences);
                examples[std::string(header)] = std::move(example);
                if (examples.size() > max_examples)
                    examples.erase(std::prev(examples.end()));
            }, max_listed, num_threads);
    Json::Value result;
    result["value"] = every_pair ? Json::Value(collisions == 0) : Json::Value();
    result["collisions"] = Json::Value(static_cast<Json::UInt64>(collisions));
    result["examples"] = Json::Value(Json::arrayValue);
    for (auto &[header, example] : examples) {
        result["examples"].append(std::move(example));
    }
    return result;
}

Json::Value labels_disjoint_json(const std::vector<HeaderPair> &pairs, size_t num_threads,
                                 size_t max_examples, size_t max_listed) {
    bool every_pair = true;
    std::vector<size_t> with_seqs;
    for (size_t p = 0; p < pairs.size(); ++p) {
        if (pairs[p].coord_to_header) {
            with_seqs.push_back(p);
        } else {
            every_pair = false;
        }
    }
    // every column name with its pair, by name, then pair (views into the annotations' labels:
    // 24 bytes per column); none is needed when no pair has headers
    std::vector<std::pair<std::string_view, size_t>> names;
    if (!with_seqs.empty()) {
        for (size_t p = 0; p < pairs.size(); ++p) {
            if (!pairs[p].columns)
                continue;
            for (const std::string &name : *pairs[p].columns) {
                names.emplace_back(name, p);
            }
        }
        std::sort(names.begin(), names.end());
    }
    // the distinct names (their first entry), and whether each is a header of a pair
    std::vector<size_t> distinct;
    for (size_t i = 0; i < names.size(); ++i) {
        if (!i || names[i].first != names[i - 1].first)
            distinct.push_back(i);
    }
    std::vector<uint8_t> is_header(distinct.size(), false);
    #pragma omp parallel for num_threads(num_threads ? num_threads : get_num_threads()) schedule(dynamic, 1024)
    for (size_t d = 0; d < distinct.size(); ++d) {
        const std::string_view name = names[distinct[d]].first;
        for (size_t p : with_seqs) {
            if (pairs[p].coord_to_header->find_header(name)) {
                is_header[d] = true;
                break;
            }
        }
    }
    auto tags = [&](size_t p) {
        Json::Value t;
        t["graph"] = pairs[p].graph;
        t["graph_path"] = pairs[p].graph_path;
        t["annotation_path"] = pairs[p].annotation_path;
        return t;
    };
    uint64_t collisions = 0;
    Json::Value examples(Json::arrayValue);
    for (size_t d = 0; d < distinct.size(); ++d) {
        if (!is_header[d])
            continue;
        ++collisions;
        // the names come in byte order: the first ones are the byte-smallest
        if (examples.size() >= max_examples)
            continue;
        const std::string_view name = names[distinct[d]].first;
        Json::Value example;
        example["header"] = std::string(name);
        Json::Value records(Json::arrayValue);
        for (size_t p : with_seqs) {
            if (records.size() >= max_listed)
                break;
            if (auto found = pairs[p].coord_to_header->find_header(name)) {
                Json::Value r = tags(p);
                const auto column = found->first;
                r["column"] = pairs[p].columns && column < pairs[p].columns->size()
                        ? Json::Value((*pairs[p].columns)[column])
                        : Json::Value();
                records.append(std::move(r));
            }
        }
        Json::Value columns(Json::arrayValue);
        const size_t end = d + 1 < distinct.size() ? distinct[d + 1] : names.size();
        for (size_t i = distinct[d]; i < end && columns.size() < max_listed; ++i) {
            // (a pair listed twice for one name would repeat: a name labels one column of a pair)
            if (i == distinct[d] || names[i].second != names[i - 1].second)
                columns.append(tags(names[i].second));
        }
        example["records"] = std::move(records);
        example["columns"] = std::move(columns);
        examples.append(std::move(example));
    }
    Json::Value result;
    result["value"] = every_pair ? Json::Value(collisions == 0) : Json::Value();
    result["collisions"] = Json::Value(static_cast<Json::UInt64>(collisions));
    result["examples"] = std::move(examples);
    return result;
}

std::optional<FileStamp> file_stamp(const std::string &path) {
    struct stat st;
    if (::stat(path.c_str(), &st))
        return std::nullopt;
#if defined(__APPLE__)
    const struct timespec mtime = st.st_mtimespec;
    const struct timespec ctime = st.st_ctimespec;
#else
    const struct timespec mtime = st.st_mtim;
    const struct timespec ctime = st.st_ctim;
#endif
    FileStamp stamp;
    stamp.device = static_cast<uint64_t>(st.st_dev);
    stamp.inode = static_cast<uint64_t>(st.st_ino);
    stamp.size = static_cast<uint64_t>(st.st_size);
    stamp.mtime_ns = static_cast<int64_t>(mtime.tv_sec) * 1000000000 + mtime.tv_nsec;
    stamp.ctime_ns = static_cast<int64_t>(ctime.tv_sec) * 1000000000 + ctime.tv_nsec;
    return stamp;
}

std::vector<PairFileStamp> pair_file_stamps(const std::string &graph,
                                            const std::string &annotation, bool coord_mapping) {
    std::vector<PairFileStamp> stamps;
    for (const IndexFile &f : index_load_inventory(graph, annotation, coord_mapping)) {
        stamps.push_back({ f.path, f.role, file_stamp(f.path) });
    }
    for (const IndexDerivedFile &f : index_derived_files(graph)) {
        stamps.push_back({ f.path, f.role, file_stamp(f.path) });
    }
    return stamps;
}

std::string pair_files_changed(const std::vector<PairFileStamp> &stamps,
                               const std::string &graph, const std::string &annotation) {
    for (const PairFileStamp &f : stamps) {
        const std::optional<FileStamp> now = file_stamp(f.path);
        if (now == f.stamp)
            continue;
        const char *how = !now ? "it is gone or cannot be read"
                : !f.stamp ? "it appeared"
                : now->device != f.stamp->device || now->inode != f.stamp->inode
                        ? "it is another file"
                : now->size != f.stamp->size ? "its size"
                : now->mtime_ns != f.stamp->mtime_ns ? "its modification time"
                : "its status change time";
        return fmt::format("The file '{}' ({}) of the pair ('{}', '{}') changed on disk since the "
                           "server started ({}): a request does not load it into RAM, since the "
                           "answer would carry the start-up identity and certificates; restart "
                           "onto the new files, or send the request without in_ram", f.path,
                           f.role, graph, annotation, how);
    }
    return "";
}

} // namespace cli
} // namespace mtg
