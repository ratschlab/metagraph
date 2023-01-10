#ifndef __DISK_UTILS__
#define __DISK_UTILS__

#include <cinttypes>
#include <vector>
#include <fstream>

namespace mtg {
namespace annot {
namespace matrix {

class DiskWriter {
  public:
    DiskWriter(std::ofstream &out, size_t buff_size) : out_(out) {
        buff_.reserve(buff_size);
    }

    void add(uint64_t v, uint8_t n_bits) {
        buff_write_.add(v, n_bits, [&](uint64_t x) { store(x); });
    }

    void flush() {
        buff_write_.flush([&](uint64_t x) { store(x); });
        flush_buffer();
    }

    ~DiskWriter() {
        flush();
    }

  private:
    class Uint64BuffWrite {
      public:
        template <class Store>
        void add(uint64_t v, uint8_t n_bits, Store store) {
            assert(n_bits && n_bits <= (uint8_t)64);
            if (pos_ == 64) {
                store(buffer_);
                buffer_ = v;
                pos_ = n_bits;
                return;
            }
            buffer_ |= v << pos_;
            auto free_bits = 64 - pos_;
            if (free_bits >= n_bits) {
                pos_ += n_bits;
            } else {
                store(buffer_);
                buffer_ = v >> free_bits;
                pos_ = n_bits - free_bits;
            }
        }

        template <class Store>
        void flush(Store store) {
            if (pos_) {
                store(buffer_);
                pos_ = 0;
                buffer_ = 0;
            }
        }

      private:
        uint64_t buffer_ = 0;
        uint8_t pos_ = 0;
    };

    std::vector<uint64_t> buff_;
    Uint64BuffWrite buff_write_;
    std::ofstream &out_;

    void flush_buffer() {
        for (uint64_t x : buff_) {
            serialize_number(out_, x);
        }
        buff_.resize(0);
    }

    void store(uint64_t x) {
        buff_.push_back(x);
        if (buff_.size() == buff_.capacity())
            flush_buffer();
    }
};


class DiskRandomReader {
  public:
    DiskRandomReader(std::ifstream &in, size_t buff_size) : in_(in) {
        offset_ = in_.tellg();

        in_.seekg(0, std::ios::end);

        file_size_ = in.tellg();

        in_.seekg(offset_, std::ios::beg);

        buffer_start_pos_ = offset_;
        buffer_.resize(buff_size);
        load(buffer_start_pos_);
    }

    void start_reading_at(uint64_t bit_pos) {
        i_ = bit_pos / 64 - (buffer_start_pos_ - offset_) / 8;
        reaload_if_needed();
        buff_read_.init(buffer_[i_++], bit_pos % 64);
    }

    uint64_t get(uint8_t n_bits) {
        uint64_t v;
        buff_read_.get(v, n_bits, [&](uint64_t& x) {
            reaload_if_needed();
            x = buffer_[i_++];
        });
        return v;
    }

  private:
    class Uint64BuffRead {
      public:
        template<class Load>
        void get(uint64_t &v, uint8_t n_bits, Load load) {
            assert(n_bits && n_bits <= (uint8_t)64);
            if (pos_ == 64) {
                load(buffer_);
                pos_ = 0;
            }
            auto have = 64 - pos_;
            if (have >= n_bits) {
                v = buffer_ >> pos_;
                if (n_bits != 64)
                    v &= (1ull << n_bits) - 1;
                pos_ += n_bits;
            } else {
                v = buffer_ >> pos_;
                load(buffer_);
                v |= buffer_ << have;
                if (n_bits != 64)
                    v &= (1ull << n_bits) - 1;
                pos_ = n_bits - have;
            }
        }

        void init(uint64_t buffer, uint8_t pos) {
            buffer_ = buffer;
            pos_ = pos;
        }

      private:
        uint64_t buffer_;
        uint8_t pos_ = 64;
    };

    std::vector<uint64_t> buffer_;
    size_t i_ = 0; // index in `buffer_`
    size_t buffer_start_pos_ = 0; // in bytes
    std::ifstream &in_;
    size_t offset_;
    Uint64BuffRead buff_read_;
    size_t file_size_;

    void load(size_t pos) {
        size_t to_load = std::min(buffer_.size(), (file_size_ - pos) / 8); // how many uint64_t words to load
        for (size_t i = 0; i < to_load; ++i)
            buffer_[i] = load_number(in_);
    }

    void reaload_if_needed() {
        auto pos = buffer_start_pos_ + i_ * 8; // convert to position in file
        if (pos < buffer_start_pos_ || pos >= buffer_start_pos_ + 8 * buffer_.size()) {
            in_.seekg(pos, std::ios::beg);
            buffer_start_pos_ = pos;
            load(buffer_start_pos_);
        }
    }
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __DISK_UTILS__
