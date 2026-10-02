#ifndef __METAGRAPH_SERVER_UTILS_HPP__
#define __METAGRAPH_SERVER_UTILS_HPP__

#include <server_http.hpp>


namespace mtg {
namespace cli {

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

// Runs |process| on the request body and writes its JSON result; compressed (gzip
// preferred, else deflate) when the client's Accept-Encoding allows it. |compact|
// writes the JSON without indentation (the traversal routes, whose bodies are large).
void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     size_t request_id,
                     const std::function<Json::Value(const std::string &)> &process,
                     bool compact = false);

class CurrentlyInitializingError : public std::runtime_error {
  public:
    CurrentlyInitializingError()
        : std::runtime_error("Server is currently initializing") {}
};

Json::Value parse_json_string(const std::string &msg);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_UTILS_HPP__
