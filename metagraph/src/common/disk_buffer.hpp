#ifndef __DISK_BUFFER__
#define __DISK_BUFFER__

#include <cinttypes>
#include <vector>
#include <fstream>

namespace mtg {
namespace common {

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

        buffer_start_ = 0;
        buffer_.resize(buff_size);
        load_buffer(offset_);
    }

    void start_reading_at(uint64_t bit_pos) {
        i_ = bit_pos / 64;
        buff_read_.init(get_next_word(), bit_pos % 64);
    }

    uint64_t get(uint8_t n_bits) {
        uint64_t v;
        buff_read_.get(v, n_bits, [&]() {
            return get_next_word();
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
                buffer_ = load();
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
                buffer_ = load();
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
    size_t i_ = 0; // index of the uint64_t word fed to `buff_read_`
    size_t buffer_start_ = 0; // in uint64_t words
    std::ifstream &in_;
    size_t offset_;
    Uint64BuffRead buff_read_;
    size_t file_size_;

    void load_buffer(size_t pos) {
        size_t to_load = std::min(buffer_.size(), (file_size_ - pos) / 8); // how many uint64_t words to load
        for (size_t i = 0; i < to_load; ++i) {
            buffer_[i] = load_number(in_);
        }
    }

    uint64_t get_next_word() {
        if (i_ < buffer_start_ || i_ >= buffer_start_ + buffer_.size()) {
            auto pos = offset_ + i_ * 8; // convert to position in file
            in_.seekg(pos, std::ios::beg);
            buffer_start_ = i_;
            load_buffer(pos);
        }
        return buffer_[i_++ - buffer_start_];
    }
};

} // namespace common
} // namespace mtg

#endif // __DISK_BUFFER__
