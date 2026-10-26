#ifndef __METAGRAPH_SERVER_UTILS_HPP__
#define __METAGRAPH_SERVER_UTILS_HPP__

#include <server_http.hpp>
#include "server_checks.hpp"


namespace mtg {
namespace cli {

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

/**
 * Run a handler returning serialized JSON, then write the response.
 * on_sent receives the socket write's error code; success does not confirm
 * receipt by the client application. /search uses failures to prioritize
 * retaining results for retries.
 */
void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     size_t request_id,
                     const std::function<std::string(const std::string &)> &process,
                     std::function<void(const SimpleWeb::error_code &)> on_sent = nullptr);

// JSON routes with bounded serialization and cancellation.
void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     size_t request_id,
                     const std::function<Json::Value(const std::string &)> &process,
                     bool compact,
                     const ResponseControl *control = nullptr,
                     std::function<void(const SimpleWeb::error_code &)> on_sent = nullptr);

std::string requested_encoding(const std::shared_ptr<HttpServer::Request> &request);
bool client_gone(const HttpServer::Request &request);

Json::Value parse_json_string(const std::string &msg);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_UTILS_HPP__
