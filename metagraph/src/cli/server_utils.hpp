#ifndef METAGRAPH_SERVER_UTILS_HPP
#define METAGRAPH_SERVER_UTILS_HPP

#include "server_http.hpp"

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

void process_request(std::shared_ptr<HttpServer::Response> &response,
                     std::shared_ptr<HttpServer::Request> &request,
                     const std::function<std::string(const std::string &)> &process);

Json::Value parse_json_string(const std::string &msg);

#endif // METAGRAPH_SERVER_UTILS_HPP