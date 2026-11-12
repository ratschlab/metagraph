#ifndef __METAGRAPH_SERVER_UTILS_HPP__
#define __METAGRAPH_SERVER_UTILS_HPP__

#include <initializer_list>
#include <string_view>

#include <server_http.hpp>


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

class CurrentlyInitializingError : public std::runtime_error {
  public:
    CurrentlyInitializingError()
        : std::runtime_error("Server is currently initializing") {}
};

Json::Value parse_json_string(const std::string &msg);

// Refuse a request that is not an object or has a field outside |known|, so that a
// misspelled field is never ignored
void check_request_fields(const Json::Value &request,
                          std::initializer_list<std::string_view> known);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_UTILS_HPP__
