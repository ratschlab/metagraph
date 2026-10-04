#ifndef __METAGRAPH_SERVER_UTILS_HPP__
#define __METAGRAPH_SERVER_UTILS_HPP__

#include <functional>
#include <optional>
#include <stdexcept>

#include <json/json.h>
#include <server_http.hpp>

#include "server_checks.hpp"


namespace mtg {
namespace cli {

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

// Answers a request with |status| and the JSON |body| (an error that states more than its
// message: the usage of a ledger-managed /traverse, the state of a conflicting attempt)
class HttpError : public std::runtime_error {
  public:
    HttpError(int status, Json::Value body)
          : std::runtime_error(body.get("error", "").asString()), status_(status),
            body_(std::move(body)) {}
    int status() const { return status_; }
    const Json::Value& body() const { return body_; }

  private:
    int status_;
    Json::Value body_;
};

// The client of the request is gone: nothing is written, and the connection is closed
class ClientGone : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

/**
 * Whether the client of |request| is gone: peer_closed() on its socket. Safe from the
 * handler's thread: while a handler runs, the HTTP server has no read pending on the socket
 * (the next read starts after the response was sent), the only other thread that can touch it
 * is the content timeout's, whose shutdown is safe beside a peek, and the connection, with its
 * descriptor, lives as long as the handler holds the Request.
 */
bool client_gone(const HttpServer::Request &request);

// What process_request does besides writing the result (null: nothing, as before)
struct ResponseControl {
    // called while the body is written (every 64 KiB of JSON text) and compressed (every
    // block): throws to stop — ClientGone (nothing is written) or HttpError (that answer)
    std::function<void()> check;
    // after the handler is done: the status written and the body's size in bytes, or 0 and
    // nullopt when nothing was written (the client is gone)
    std::function<void(int status, std::optional<size_t> bytes)> on_written;
    // the zlib level of a compressed body (1-9; 9, the best compression, for every route
    // that does not choose: the traversal routes choose a faster one, see
    // Config::traverse_compression_level)
    int compression_level = 9;
    // after a body of |text_bytes| was compressed in |seconds| (a delivery rate measured)
    std::function<void(size_t text_bytes, double seconds)> on_compressed;
    // writes the result as text instead of json_text(result, compact, check) — a route that
    // wrote parts of it already (the /traverse results, each written once built) assembles
    // them here; must give the same bytes json_text would
    std::function<std::string(const Json::Value &result, const std::function<void()> &check)> write;
};

// Runs |process| on the request body and writes its JSON result; compressed with gzip or
// deflate when the client's Accept-Encoding makes one acceptable (RFC 9110 weights; gzip
// on a tie), uncompressed otherwise. |compact|
// writes the JSON without indentation (the traversal routes, whose bodies are large).
// |process| may throw HttpError (that status and body) or ClientGone (no response); with
// |control|, the text and the compression are written under its check.
void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     size_t request_id,
                     const std::function<Json::Value(const std::string &)> &process,
                     bool compact = false,
                     const ResponseControl *control = nullptr);

class CurrentlyInitializingError : public std::runtime_error {
  public:
    CurrentlyInitializingError()
        : std::runtime_error("Server is currently initializing") {}
};

Json::Value parse_json_string(const std::string &msg);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_UTILS_HPP__
