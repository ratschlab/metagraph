#ifndef __FILE_UTILS_HPP__
#define __FILE_UTILS_HPP__

#include <cassert>
#include <fstream>
#include <future>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>


namespace utils {

void set_swap_path(std::filesystem::path dir_path);
std::filesystem::path get_swap_path();

std::filesystem::path create_temp_dir(std::filesystem::path path,
                                      const std::string &name = "");
void remove_temp_dir(std::filesystem::path dir_name);

bool check_if_writable(const std::string &filename);


class TempFile {
  public:
    // The state flow:
    //    init -> APPEND -> READ -> deinit
    enum State { APPEND, READ };

    TempFile(const std::string &tmp_dir = "");
    ~TempFile();

    std::ofstream& ofstream();
    std::ifstream& ifstream();

    const std::string& name() const;

  private:
    std::string tmp_file_name_;
    State state_;
    std::unique_ptr<std::ofstream> tmp_ostream_;
    std::unique_ptr<std::ifstream> tmp_istream_;
};


template <typename T>
class BufferedAsyncWriter {
    static constexpr uint32_t CAPACITY = 100'000;

  public:
    BufferedAsyncWriter(const std::string &name, std::ofstream *fos)
          : name_(name), fos_(fos) {
        assert(fos);

        buf_.reserve(CAPACITY);
        buf_dump_.reserve(CAPACITY);
    }

    ~BufferedAsyncWriter() { flush(); }

    void push(const T &v) {
        if (buf_.size() == CAPACITY) {
            wait_for_write();
            buf_.swap(buf_dump_);
            write_future_ = std::async(std::launch::async, flush_to_stream,
                                       buf_dump_, fos_, name_);
            buf_.resize(0);
        }
        buf_.push_back(v);
    }

    void flush() {
        // wait for the other to finish
        wait_for_write();

        // dump the current buffer
        flush_to_stream(buf_, fos_, name_);
        buf_.resize(0);

        fos_->flush();
    }

  private:
    void wait_for_write() {
        if (write_future_.valid())
            write_future_.wait();
    }

    static void flush_to_stream(const std::vector<T> &buf,
                                std::ofstream *fos,
                                const std::string &name) {
        if (!fos->write(reinterpret_cast<const char *>(buf.data()), sizeof(T) * buf.size())) {
            std::cerr << "Error: Writing to '" << name << "' failed" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    std::vector<T> buf_;
    std::vector<T> buf_dump_;

    std::shared_future<void> write_future_;

    std::string name_;
    std::ofstream *fos_;
};


class NamedIfstream : public std::ifstream {
  public:
    NamedIfstream(const std::string &fname, std::ios_base::openmode mode)
        : std::ifstream(fname, mode), fname_(fname) {}

    const std::string& get_name() const { return fname_; }

  private:
    std::string fname_;
};

} // namespace utils

#endif // __FILE_UTILS_HPP__
