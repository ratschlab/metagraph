#include "server.hpp"

#include <ctime>
#include <iostream>
#include <string>
#include <utility>

#include <asio.hpp>

using asio::ip::tcp;

// Based on example http://think-async.com/Asio/asio-1.12.2/src/examples/cpp11/echo/async_tcp_echo_server.cpp


class Connection : public std::enable_shared_from_this<Connection> {
    static constexpr size_t MAX_MESSAGE_SIZE = 1'000'000;

  public:
    Connection(tcp::socket socket,
               std::function<std::string(const std::string &input)> generate_response)
          : socket_(std::move(socket)),
            generate_response_(generate_response) {
        std::cout << "New connection with "
                  << socket_.remote_endpoint().address()
                  << ":"
                  << socket_.remote_endpoint().port()
                  << std::endl;
    }

    ~Connection() {
        std::cout << "Closing connection with "
                  << socket_.remote_endpoint().address()
                  << ":"
                  << socket_.remote_endpoint().port()
                  << std::endl;
    }

    void start() { do_read(); }

  private:
    void do_read() {
        auto self(shared_from_this());
        socket_.async_read_some(asio::buffer(data_, MAX_MESSAGE_SIZE),
            [this, self](std::error_code ec, size_t length) {
                if (!ec)
                    do_write(length);
            }
        );
    }

    void do_write(size_t length) {
        auto self(shared_from_this());

        response_ = generate_response_(std::string(data_, length));

        asio::async_write(socket_, asio::buffer(response_),
            [this, self](std::error_code ec, size_t /*length*/) {
                if (!ec)
                    do_read();
            }
        );
    }

    tcp::socket socket_;
    char data_[MAX_MESSAGE_SIZE];
    std::function<std::string(const std::string &input)> generate_response_;
    std::string response_;
};


Server::Server(asio::io_context &io_context, short port,
               std::function<std::string(const std::string &input)> generate_out)
      : acceptor_(io_context, tcp::endpoint(tcp::v4(), port)),
        generate_response_(generate_out) {
    do_accept();
}

void Server::do_accept() {
    acceptor_.async_accept([this](std::error_code ec, tcp::socket socket) {
        if (!ec)
            std::make_shared<Connection>(std::move(socket), generate_response_)->start();

        do_accept();
    });
}
