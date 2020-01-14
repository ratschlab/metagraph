#include "file_utils.hpp"

#include <cassert>
#include <filesystem>
#include <algorithm>

#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#else
#include <stdio.h>
#endif

namespace utils {

bool check_if_writable(const std::string &filename) {
    std::ifstream ifstream(filename, std::ios::binary);
    bool existed = ifstream.good();
    ifstream.close();

    std::ofstream ofstream(filename, std::ios::binary
                                        | std::ofstream::ios_base::app);
    bool can_write = ofstream.good();
    ofstream.close();

    if (!can_write)
        return false;

    if (!existed)
        std::remove(filename.c_str());

    return true;
}


TempFile::TempFile(const std::string &tmp_dir)
      : tmp_file_name_((tmp_dir.size()
                          ? tmp_dir
                          : std::filesystem::temp_directory_path().string())
                                                + std::string("/tmp.XXXXXX")) {
    // create a file
    int fd = mkstemp(tmp_file_name_.data());
    if (fd == -1)
        throw std::runtime_error("Error: temp file "
                                    + tmp_file_name_ + " creation failed");
    // close the file descriptor
    close(fd);

    tmp_ostream_.reset(new std::ofstream(tmp_file_name_,
                                         std::ios::binary | std::ios::app));
    if (!tmp_ostream_->good()) {
        unlink(tmp_file_name_.c_str());
        throw std::runtime_error("Error: temp file "
                                    + tmp_file_name_ + " open failed");
    }
    state_ = APPEND;
}

TempFile::~TempFile() {
    unlink(tmp_file_name_.c_str());
}

std::ofstream& TempFile::ofstream() {
    assert(state_ == APPEND && "Can't write after reading");
    return *tmp_ostream_;
}

std::ifstream& TempFile::ifstream() {
    if (!tmp_istream_.get()) {
        state_ = READ;
        tmp_ostream_.reset();
        tmp_istream_.reset(new std::ifstream(tmp_file_name_, std::ios::binary));
    }
    return *tmp_istream_;
}

} // namespace utils
