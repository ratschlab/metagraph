#include "file_utils.hpp"

#include <cassert>
#include <csignal>
#include <cstdlib>
#include <cerrno>
#include <filesystem>
#include <algorithm>
#include <system_error>
#include <string>

#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#else
#include <stdio.h>
#endif

#include <sdsl/memory_management.hpp>

#include "common/logger.hpp"

namespace utils {

using mtg::common::logger;
namespace fs = std::filesystem;

// Function-local statics instead of namespace-scope globals to sidestep a
// static-initialization-order hazard: test_dump_dir() in test_helpers.hpp
// calls create_temp_dir() during a test TU's dynamic init, which push_backs
// into tmp_dirs. If that push runs before this TU's own dynamic init, the
// vector's default-member-initializers (_M_start = nullptr, ...) later
// re-run via the default constructor and silently drop the pushed entry.
// The mutex (constexpr default ctor) and pid_t were not strictly at risk;
// they're wrapped here for consistency.
static std::filesystem::path& swap_path() {
    static std::filesystem::path instance;
    return instance;
}
static std::vector<std::string>& tmp_dirs() {
    static std::vector<std::string> instance;
    return instance;
}
static std::mutex& tmp_dirs_mutex() {
    static std::mutex instance;
    return instance;
}
// The pid of the process that first registered a temp dir. Used to skip
// cleanup in forked children — e.g. gtest's threadsafe death-test style forks
// the test process, and the child's atexit would otherwise wipe the parent's
// temp dirs out from under it.
static pid_t& tmp_dirs_owner_pid() {
    static pid_t instance = 0;
    return instance;
}

static bool WITH_MMAP = false;
static bool WITH_MADVISE = false;

void set_mmap(bool set_bit) {
    WITH_MMAP = set_bit;
}

bool with_mmap() {
    return WITH_MMAP;
}

void set_madvise(bool set_bit) {
    WITH_MADVISE = set_bit;
}

bool with_madvise() {
    return WITH_MADVISE;
}

std::unique_ptr<std::ifstream> open_ifstream(const std::string &filename, bool mmap_stream) {
    std::unique_ptr<std::ifstream> in;
    if (mmap_stream) {
        in.reset(new sdsl::mmap_ifstream(filename, std::ios_base::binary));
    } else {
        in.reset(new named_ifstream(filename, std::ios_base::binary));
    }
    return in;
}

const std::string& get_filename(std::istream &f) {
    if (auto *mmap_in = dynamic_cast<sdsl::mmap_ifstream *>(&f))
        return mmap_in->get_filename();
    if (auto *named_in = dynamic_cast<named_ifstream *>(&f))
        return named_in->get_filename();
    throw std::runtime_error(
            "get_filename: stream is neither sdsl::mmap_ifstream nor utils::named_ifstream");
}

void madvise_random_range(std::istream &f, std::streamoff start, std::streamoff length) {
    if (!with_madvise())
        return;
    auto *mmap_in = dynamic_cast<sdsl::mmap_ifstream *>(&f);
    if (!mmap_in)
        return;
    if (length < 0)
        length = mmap_in->get_mmap_context()->file_size_bytes() - start;
    const std::streamoff pagesize = sysconf(_SC_PAGESIZE);
    std::streamoff aligned_start = start - start % pagesize;
    if (madvise(mmap_in->get_mmap_context()->data() + aligned_start,
                length + (start - aligned_start),
                MADV_RANDOM)) {
        logger->warn("madvise(MADV_RANDOM) failed for range [{}, {})", start, start + length);
    }
}

void load_mmap_random(const std::string &filename, std::streamoff offset,
                      const std::function<void(std::istream &)> &fn) {
    sdsl::mmap_ifstream in(filename);
    in.seekg(offset, std::ios_base::beg);
    fn(in);
    madvise_random_range(in, offset);
}

std::string file_read_failure_detail(const fs::path &path) {
    std::error_code ec;
    fs::file_status st = fs::status(path, ec);
    // fs::status sets |ec| for errors it can observe via stat(2): nonexistent
    // path, or a parent directory that we can't traverse (EACCES bubbles up
    // as std::errc::permission_denied). It does NOT report the case where
    // the file itself exists and is stat-able but isn't readable — stat
    // checks directory-traversal perms, not file-read perms. The access()
    // call below fills that gap.
    if (ec) {
        if (ec == std::errc::permission_denied)
            return std::string(ec.message())
                    + " (check read permissions on this path and its parent directories)";
        return ec.message();
    }
    if (!fs::exists(st))
        return "no such file or directory";
#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    if (::access(path.c_str(), R_OK) != 0) {
        const char *err = std::strerror(errno);
        return std::string(err ? err : "access failed")
                + " (missing read permission on this file)";
    }
#endif
    return "could not read contents (possibly corrupted or incompatible)";
}

std::ofstream open_new_ofstream(const std::string &filename) {
    // If file already exists, remove it.
    // This ensures that if that old file was open for reading (e.g., memory mapped),
    // the readers can keep reading it and the new out stream points to a new inode.
    std::filesystem::remove(filename);
    return std::ofstream(filename, std::ios_base::binary);
}

void set_swap_path(std::filesystem::path dir_path) {
    swap_path() = dir_path;
}

std::filesystem::path get_swap_path() {
    return swap_path();
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
    if (getpid() != tmp_dirs_owner_pid())
        return;
    auto &dirs = tmp_dirs();
    std::for_each(dirs.begin(), dirs.end(), cleanup_temp_dir_nolock);
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

    std::lock_guard<std::mutex> lock(tmp_dirs_mutex());

    auto &dirs = tmp_dirs();
    if (dirs.empty()) {
        tmp_dirs_owner_pid() = getpid();
        if (std::signal(SIGINT, cleanup_tmp_dir_on_signal) == SIG_ERR)
            logger->error("Couldn't reset the signal handler for SIGINT");
        if (std::signal(SIGTERM, cleanup_tmp_dir_on_signal) == SIG_ERR)
            logger->error("Couldn't reset the signal handler for SIGTERM");
        if (std::atexit(cleanup_tmp_dir_on_exit))
            logger->error("Couldn't reset the atexit handler");
    }

    logger->trace("Registered temporary directory {}", tmp_dir_str);

    dirs.push_back(tmp_dir_str);

    return tmp_dir_str;
}

void remove_temp_dir(std::filesystem::path dir_name) {
    {
        std::lock_guard<std::mutex> lock(tmp_dirs_mutex());
        auto &dirs = tmp_dirs();
        auto it = std::find(dirs.begin(), dirs.end(), dir_name);
        assert(it != dirs.end());
        dirs.erase(it);
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


void rename_or_move_file(const std::string &old_fname, const std::string &fname) {
    std::error_code ec;
    fs::rename(old_fname, fname, ec);
    // check if rename was successful
    if (!ec)
        return;
    // if cross-device, fallback to copy + remove
    if (ec == std::errc::cross_device_link) {
        fs::copy_file(old_fname, fname, fs::copy_options::overwrite_existing);
        fs::remove(old_fname);
        return;
    }
    // other errors
    throw fs::filesystem_error("rename failed", old_fname, fname, ec);
}


void append_file(const std::string &source_fname, const std::string &target_fname) {
    std::string concat_command = "cat " + source_fname + " >> " + target_fname;
    if (std::system(concat_command.c_str())) {
        logger->error("Error while cat-ing files: {}", concat_command);
        exit(1);
    }
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
