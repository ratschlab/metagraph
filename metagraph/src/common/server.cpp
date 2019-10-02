#include "server.hpp"

#include <ctime>
#include <iostream>
#include <string>
#include <utility>

#include <asio.hpp>
#include <json/json.h>

using asio::ip::tcp;
using asio::detail::socket_ops::host_to_network_long;
using asio::detail::socket_ops::network_to_host_long;

const size_t MAX_MESSAGE_SIZE = 1'000'000'000;


class Connection : public std::enable_shared_from_this<Connection> {
  public:
    Connection(tcp::socket socket,
               std::function<std::string(const std::string &input)> generate_response)
          : socket_(std::move(socket)),
            generate_response_(generate_response) {
        try {
            address_ = socket_.remote_endpoint().address();
            port_ = socket_.remote_endpoint().port();
        } catch (...) {
            std::cerr << "Error in remote endpoint" << std::endl;
        }
        std::cout << "New connection " << address_ << ":" << port_ << std::endl;
    }

    ~Connection() {
        std::cout << "Closing connection " << address_ << ":" << port_ << std::endl;
    }

    void start() {
        // start the routine asynchronously
        auto self(shared_from_this());
        socket_.get_io_service().post([this, self]() { do_read(); });
    }

  private:
    void do_read() {
        uint64_t message_size;
        try {
            message_size = receive_number();
        } catch (const asio::system_error &e) {
            if (e.code() != asio::error::eof)
                print_reading_error();
            return;
        }

        if (message_size > MAX_MESSAGE_SIZE) {
            std::cerr << "Error: Trying to reveive message of size "
                      << message_size << " bytes from client "
                      << address_ << ":" << port_ << ". Maximum message size: "
                      << MAX_MESSAGE_SIZE << std::endl;
            return;
        }

        data_.resize(message_size);

        auto self(shared_from_this());
        asio::async_read(socket_, asio::buffer(&data_[0], message_size),
            [this, self](std::error_code ec, size_t /*length*/) {
                if (ec) {
                    print_reading_error();
                    return;
                }

                try {
                    data_ = generate_response_(data_);
                    do_write();
                } catch (const std::exception &ex) {
                    std::cerr << "Error: Exception thrown when trying to "
                              << "generate a response for client "
                              << address_ << ":" << port_ << std::endl;

                    Json::Value message;
                    message["error"] = ex.what();

                    Json::StreamWriterBuilder builder;
                    data_ = Json::writeString(builder, message);

                    do_write();
                    return;
                }
            }
        );
    }

    void do_write() {
        try {
            send_number(data_.size());
        } catch (const asio::system_error &e) {
            print_writing_error();
            return;
        }

        auto self(shared_from_this());
        asio::async_write(socket_, asio::buffer(data_),
            [this, self](std::error_code ec, size_t /*length*/) {
                if (!ec) {
                    do_read();
                } else {
                    print_writing_error();
                }
            }
        );
    }

    uint64_t receive_number() {
        char small_buffer_[8];
        asio::read(socket_, asio::buffer(small_buffer_));

        static_assert(sizeof(uint32_t) * 2 == sizeof(small_buffer_));
        const uint32_t high = *reinterpret_cast<const uint32_t *>(small_buffer_);
        const uint32_t low = *(reinterpret_cast<const uint32_t *>(small_buffer_) + 1);

        return (static_cast<uint64_t>(host_to_network_long(high)) << 32) | host_to_network_long(low);
    }

    void send_number(uint64_t value) {
        const uint32_t high = host_to_network_long(static_cast<uint32_t>(value >> 32));
        const uint32_t low = host_to_network_long(static_cast<uint32_t>(value & 0xFFFFFFFFLL));

        char small_buffer_[8];
        *reinterpret_cast<uint32_t *>(small_buffer_) = high;
        *(reinterpret_cast<uint32_t *>(small_buffer_) + 1) = low;

        asio::write(socket_, asio::buffer(small_buffer_));
    }

    void print_reading_error() const {
        std::cerr << "Error: can't receive data from "
                  << address_ << ":" << port_ << std::endl;
    }

    void print_writing_error() const {
        std::cerr << "Error: can't send data to "
                  << address_ << ":" << port_ << std::endl;
    }

    tcp::socket socket_;
    asio::ip::address address_;
    unsigned short port_;

    std::function<std::string(const std::string &input)> generate_response_;
    std::string data_;
};


Server::Server(asio::io_context &io_context, short port,
               std::function<std::string(const std::string &input)> get_output)
      : acceptor_(io_context, tcp::endpoint(tcp::v4(), port)),
        get_output_(get_output) {
    do_accept();
}

void Server::do_accept() {
    acceptor_.async_accept([this](std::error_code ec, tcp::socket socket) {
        if (!ec)
            std::make_shared<Connection>(std::move(socket), get_output_)->start();

        do_accept();
    });
}
