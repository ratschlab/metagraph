#include <cctype>
#include <cerrno>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <streambuf>

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
        ret = json_text(process(content), compact, check);
        const std::string encoding = requested_encoding(request);
        if (!encoding.empty()) {
            ret = compress_string(ret, Z_BEST_COMPRESSION, encoding == "gzip", check);
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
