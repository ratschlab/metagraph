#ifndef __METAGRAPH_SERVER_UTILS_HPP__
#define __METAGRAPH_SERVER_UTILS_HPP__

#include <server_http.hpp>


namespace mtg {
namespace cli {

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

/**
 * Run the JSON-in / JSON-out request handler `process` and write the result back
 * to the client.
 *
 * If `on_sent` is non-null it's invoked with the asio error_code returned by
 * the async network write — empty error_code means the bytes left this host
 * without a socket-level error; a non-empty code means the connection failed
 * before/during delivery. Use this to drive a result-cache state machine
 * (e.g. mark FAILED → keep result for retries).
 */
void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     size_t request_id,
                     const std::function<Json::Value(const std::string &)> &process,
                     std::function<void(const SimpleWeb::error_code &)> on_sent = nullptr);

class CurrentlyInitializingError : public std::runtime_error {
  public:
    CurrentlyInitializingError()
        : std::runtime_error("Server is currently initializing") {}
};

Json::Value parse_json_string(const std::string &msg);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_UTILS_HPP__
