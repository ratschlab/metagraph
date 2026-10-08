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

// HttpError, ClientGone, CurrentlyInitializingError, ResponseControl and answer_request, the
// parts that need no HTTP library, are in server_checks.hpp

/**
 * Whether the client of |request| is gone: peer_closed() on its socket. Safe from the
 * handler's thread: while a handler runs, the HTTP server has no read pending on the socket
 * (the next read starts after the response was sent), the only other thread that can touch it
 * is the content timeout's, whose shutdown is safe beside a peek, and the connection, with its
 * descriptor, lives as long as the handler holds the Request.
 */
bool client_gone(const HttpServer::Request &request);

// Runs |process| on the request body and writes its JSON result; compressed with gzip or
// deflate when the client's Accept-Encoding makes one acceptable (RFC 9110 weights; gzip
// on a tie), uncompressed otherwise. |compact|
// writes the JSON without indentation (the traversal routes, whose bodies are large).
// |process| may throw HttpError (that status and body) or ClientGone (no response); with
// |control|, the text and the compression are written under its check, and nothing at all
// is written when its |gone| answers true (answer_request, server_checks.hpp).
void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     size_t request_id,
                     const std::function<Json::Value(const std::string &)> &process,
                     bool compact = false,
                     const ResponseControl *control = nullptr);

Json::Value parse_json_string(const std::string &msg);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_UTILS_HPP__
