#include "server_utils.hpp"
#include "common/logger.hpp"
#include "common/unix_tools.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <sstream>
#include <netinet/in.h>
#include <netinet/tcp.h>
#if defined(__APPLE__)
#include <netinet/tcp_fsm.h>
#endif
#include <sys/socket.h>

namespace mtg::cli {
using common::logger;
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
    response->send(std::move(on_sent));
    logger->info("[Server] Request {} processing time: {:.3f} sec, response size: {:.1f} KB, "
                 "finished in {:.3f} sec",
                 request_id, processing_time, (double)ret.size() / 1000, timer.elapsed());
}


} // namespace mtg::cli
