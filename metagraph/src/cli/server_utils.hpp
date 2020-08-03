#ifndef __METAGRAPH_SERVER_UTILS_HPP__
#define __METAGRAPH_SERVER_UTILS_HPP__

#include <server_http.hpp>


namespace mtg {
namespace cli {

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

void process_request(std::shared_ptr<HttpServer::Response> &response,
                     const std::shared_ptr<HttpServer::Request> &request,
                     const bool execute_separate_thread,
                     const float timeout,
                     const std::function<std::string(const std::string &)> &process);

Json::Value parse_json_string(const std::string &msg);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_UTILS_HPP__
