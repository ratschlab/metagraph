#ifndef __FILE_UTILS_HPP__
#define __FILE_UTILS_HPP__

#include <fstream>
#include <future>
#include <iostream>
#include <memory>
#include <string>


namespace utils {

bool check_if_writable(const std::string &filename);

/**
 *  Returns the input file type, given a filename
 *  One of: [KMC|VCF|FASTQ|FASTA].
 *  If filetype is unknown, return empty string "".
 */
std::string get_filetype(const std::string &fname);

class TempFile {
  public:
    // The state flow:
    //    init -> APPEND -> READ -> deinit
    enum State { APPEND, READ };

    TempFile(const std::string &tmp_dir = "");
    ~TempFile();

    std::ofstream &ofstream();
    std::ifstream &ifstream();

  private:
    std::string tmp_file_name_;
    State state_;
    std::unique_ptr<std::ofstream> tmp_ostream_;
    std::unique_ptr<std::ifstream> tmp_istream_;
};

template <typename T>
class BufferedAsyncWriter {
    static constexpr uint32_t capacity = 100000;

  public:
    BufferedAsyncWriter(const std::string &name, std::fstream *f)
        : buf1_(capacity), buf2_(capacity), buf_(&buf1_), name_(name), f_(f) {}
    static void
    flush_to_file(const std::string name, std::fstream *f, const std::vector<T> *buf) {
        if (!f->write((char *)&((*buf)[0]), sizeof((*buf)[0]) * buf->size())) {
            std::cerr << "Error: Writing to '{}' failed" << name;
            std::exit(EXIT_FAILURE);
        }
    }

    void push(const T &v) {
        if (buf_->size() == capacity) {
            wait_for_write();
            write_future = std::async(flush_to_file, name_, f_, buf_);
            buf_ = (buf_ == &buf1_ ? &buf2_ : &buf1_);
            buf_->resize(0);
        }
        buf_->push_back(v);
    }

    void flush() {
        flush_to_file(name_, f_, buf_);
        wait_for_write();
    }

  private:
    void wait_for_write() {
        if (write_future.valid()) {
            write_future.wait();
        }
    }

    std::vector<T> buf1_, buf2_, *buf_;
    std::future<void> write_future;
    std::string name_;
    std::fstream *f_;
};

} // namespace utils

#endif // __FILE_UTILS_HPP__
