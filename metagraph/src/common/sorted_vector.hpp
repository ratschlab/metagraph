#ifndef __SORTED_VECTOR__
#define __SORTED_VECTOR__

#include <vector>
#include <filesystem>

#include <ips4o.hpp>

#include "common/utils/file_utils.hpp"
#include "common/elias_fano/elias_fano.hpp"
#include "common/elias_fano/elias_fano_merger.hpp"


namespace mtg {
namespace common {

template <typename T>
class SortedVector {
  public:
    SortedVector(size_t num_threads,
                 uint64_t buffer_size,
                 const std::string &swap_dir,
                 size_t max_chunks_open = 1000)
          : num_threads_(num_threads),
            max_chunks_open_(max_chunks_open) {
        buffer_.reserve(buffer_size);
        if (swap_dir.size())
            tmp_dir_ = utils::create_temp_dir(swap_dir, "sorted_vector");
    }

    // Don't make the default move constructor and copy operator get removed
    // because of the destructor.
    SortedVector(SortedVector&& other) = default;
    SortedVector& operator=(SortedVector&& other) = default;

    ~SortedVector() {
        try {
            if (!tmp_dir_.empty())
                utils::remove_temp_dir(tmp_dir_);

        } catch (const std::exception &e) {
            std::cerr << "ERROR: Failed to destruct SortedVector: "
                      << e.what() << std::endl;
        } catch (...) {
            std::cerr << "ERROR: Failed to destruct SortedVector";
        }
    }

    template <class... Args>
    void emplace_back(Args&&... args) {
        buffer_.emplace_back(std::forward<Args>(args)...);
        total_size_++;
        if (buffer_.size() == buffer_.capacity())
            flush();
    }

    void push_back(T&& value) {
        buffer_.push_back(std::move(value));
        total_size_++;
        if (buffer_.size() == buffer_.capacity())
            flush();
    }

    void push_back(const T &value) {
        buffer_.push_back(value);
        total_size_++;
        if (buffer_.size() == buffer_.capacity())
            flush();
    }

    size_t size() const { return total_size_; }

    template <class Callback>
    void for_each(const Callback &callback) {
        if (tmp_dir_.empty()) {
            ips4o::parallel::sort(buffer_.begin(), buffer_.end(), std::less<>(), num_threads_);
            for (const auto &v : buffer_) {
                callback(v);
            }
            return;
        }

        if (merged_)
            throw std::runtime_error("ERROR: Trying to merge chunks more than once");

        flush();

        // release the buffer
        buffer_ = std::vector<T>();

        std::vector<std::string> filenames;
        // for stage 1, fwd bits are already counted, so we skip that chunk
        for (uint32_t chunk = 0; chunk < num_chunks_; ++chunk) {
            filenames.push_back(tmp_file(chunk));
        }

        const bool remove_chunks = true;
        elias_fano::merge_files<T>(filenames, callback, remove_chunks, max_chunks_open_);

        merged_ = true;
    }

    // release the buffer if `final = true`
    void flush(bool final = false) {
        if (tmp_dir_.empty() || (buffer_.empty() && num_chunks_))
            return;
        ips4o::parallel::sort(buffer_.begin(), buffer_.end(), std::less<>(), num_threads_);
        elias_fano::EliasFanoEncoderBuffered<T>::append_block(buffer_, tmp_file(num_chunks_++));
        if (final) {
            buffer_ = std::vector<T>();
        } else {
            buffer_.resize(0);
        }
    }

    std::vector<T>& get_buffer() { return buffer_; }

  private:
    std::string tmp_file(size_t chunk) const {
        return tmp_dir_/fmt::format("chunk_{}", chunk);
    };

    std::vector<T> buffer_;
    size_t num_threads_; // threads to use for sorting
    size_t max_chunks_open_; // maximum number of chunks open when merging
    std::filesystem::path tmp_dir_; // directory where to create a temp dir
    size_t num_chunks_ = 0;
    size_t total_size_ = 0;
    bool merged_ = false;
};

} // namespace common
} // namespace mtg

#endif // __SORTED_VECTOR__
