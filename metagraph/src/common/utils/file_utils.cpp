#include "file_utils.hpp"

#include <cassert>
#include <csignal>
#include <cstdlib>
#include <filesystem>
#include <algorithm>

#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#else
#include <stdio.h>
#endif

#include <sdsl/memory_management.hpp>

#include "common/logger.hpp"

namespace utils {

using mtg::common::logger;

std::filesystem::path SWAP_PATH;
std::vector<std::string> TMP_DIRS;
std::mutex TMP_DIRS_MUTEX;

static bool WITH_MMAP = false;

bool with_mmap(bool set_bit) {
    if (set_bit) {
        // TODO: is there a good way to check this, like when opening outstreams?
        // It could be logger->info() but that prints to stdout. Print info() to stderr?
        logger->warn("Memory mapping enabled. Make sure all output files are"
                     " different from the input to avoid errors.");
        WITH_MMAP = true;
    }
    return WITH_MMAP;
}

std::unique_ptr<std::ifstream>
open_ifstream(const std::string &filename, bool mmap_stream) {
    std::unique_ptr<std::ifstream> in;
    if (mmap_stream) {
        in.reset(new sdsl::mmap_ifstream(filename, std::ios_base::binary));
    } else {
        in.reset(new std::ifstream(filename, std::ios_base::binary));
    }
    return in;
}

void set_swap_path(std::filesystem::path dir_path) {
    SWAP_PATH = dir_path;
}

std::filesystem::path get_swap_path() {
    return SWAP_PATH;
}

void cleanup_tmp_dir_on_signal(int sig) {
    logger->trace("Got signal {}. Exiting...", sig);
    // call std::exit to invoke handlers registered in std::atexit
    std::exit(sig);
}

void cleanup_temp_dir_nolock(const std::filesystem::path &tmp_dir) {
    logger->trace("Cleaning up temporary directory {}", tmp_dir);
    try {
        std::filesystem::remove_all(tmp_dir);
    } catch (...) {
        logger->error("Failed to clean up temporary directory {}", tmp_dir);
    }
}

void cleanup_tmp_dir_on_exit() {
    std::for_each(TMP_DIRS.begin(), TMP_DIRS.end(), cleanup_temp_dir_nolock);
}

std::filesystem::path create_temp_dir(std::filesystem::path path,
                                      const std::string &name) {
    if (path.empty())
        path = "./";

    std::string tmp_dir_str(path/("temp_" + name + "_XXXXXX"));
    if (!mkdtemp(tmp_dir_str.data())) {
        logger->error("Failed to create a temporary directory in {}", path);
        exit(1);
    }

    std::lock_guard<std::mutex> lock(TMP_DIRS_MUTEX);

    if (TMP_DIRS.empty()) {
        if (std::signal(SIGINT, cleanup_tmp_dir_on_signal) == SIG_ERR)
            logger->error("Couldn't reset the signal handler for SIGINT");
        if (std::signal(SIGTERM, cleanup_tmp_dir_on_signal) == SIG_ERR)
            logger->error("Couldn't reset the signal handler for SIGTERM");
        if (std::atexit(cleanup_tmp_dir_on_exit))
            logger->error("Couldn't reset the atexit handler");
    }

    logger->trace("Registered temporary directory {}", tmp_dir_str);

    TMP_DIRS.push_back(tmp_dir_str);

    return tmp_dir_str;
}

void remove_temp_dir(std::filesystem::path dir_name) {
    {
        std::lock_guard<std::mutex> lock(TMP_DIRS_MUTEX);
        auto it = std::find(TMP_DIRS.begin(), TMP_DIRS.end(), dir_name);
        assert(it != TMP_DIRS.end());
        TMP_DIRS.erase(it);
    }
    cleanup_temp_dir_nolock(dir_name);
}


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
                          ? std::filesystem::path(tmp_dir)
                          : std::filesystem::temp_directory_path()
                        )/"tmp.XXXXXX") {
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

const std::string& TempFile::name() const {
    return tmp_file_name_;
}

} // namespace utils
