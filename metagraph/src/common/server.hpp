#ifndef __SERVER_HPP__
#define __SERVER_HPP__

#include <string>
#include <memory>

#include <asio.hpp>


class Server {
  public:
    Server(asio::io_context &io_context, short port,
           std::function<std::string(const std::string &input)> generate_out
               = [](const auto &input) { return "Echo: " + input; });

  private:
    void do_accept();

    asio::ip::tcp::acceptor acceptor_;
    std::function<std::string(const std::string &input)> generate_response_;
};

#endif // __SERVER_HPP__
