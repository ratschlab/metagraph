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

// How long the server reads and discards what the client of a 413 body_too_large still sends,
// at most, before it closes the connection (refuse_body_too_large)
constexpr long kBodyTooLargeDrainS = 30;
// ... and how many bytes, at least: the drain ends after the larger of this and the cap
constexpr uint64_t kBodyTooLargeDrainBytes = 64ull << 20;

/**
 * Answers |request|, whose body is longer than the server reads (--max-request-body-mb
 * |max_mb| MiB: the HTTP library's message_size, reported to on_error), with the 413 of
 * body_too_large_json, on its connection, then closes the connection; true when it answered.
 * The library refuses such a request before reading its body (a Content-Length above the
 * bound) or at the chunk that passes it (a chunked body); it prepares a 413 of its own but does
 * not send it, so the connection would be closed with the body unread, a reset for the client.
 * The answer is written here instead, on the connection reached as client_gone reaches it,
 * with Connection: close; then the server half-closes and reads and discards (never keeps)
 * what the client still sends, until the client closes, kBodyTooLargeDrainS pass (the
 * connection's own timer) or the larger of kBodyTooLargeDrainBytes and the cap is read, so that
 * a client that sends its whole body before it reads finds the 413 rather than a reset. A
 * Content-Length that content_length_value reads as no byte count ("-1", "3145728junk", which
 * the library's std::stoull reads as 2^64 - 1 and 3145728, or a list of differing counts) is
 * closed without an answer (false), as the library drops one it cannot read at all. Every step
 * runs under the library's handler runner, so a shutdown ends it. Called from on_error only,
 * with no operation of the library pending on the connection.
 */
bool refuse_body_too_large(const std::shared_ptr<HttpServer::Request> &request,
                           uint64_t max_mb);

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
                     const ResponseControl *control = nullptr,
                     std::function<void(const SimpleWeb::error_code &)> on_sent = nullptr);

// Writes |answer| (a status other than 0) as process_request writes one, without its log lines:
// the busy refusal of a handler that refuses its request before reading it (E19)
void write_answer(std::shared_ptr<HttpServer::Response> &response, const RequestAnswer &answer);

// Send an already serialized JSON body through the same compression and error handling.
void process_serialized_request(std::shared_ptr<HttpServer::Response> &response,
                                const std::shared_ptr<HttpServer::Request> &request,
                                size_t request_id,
                                const std::function<std::string(const std::string &)> &process,
                                std::function<void(const SimpleWeb::error_code &)> on_sent);


Json::Value parse_json_string(const std::string &msg);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_UTILS_HPP__
