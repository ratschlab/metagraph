#ifndef __DISK_UTILS__
#define __DISK_UTILS__
#include <cinttypes>
#include <vector>
#include <fstream>

namespace mtg {
namespace annot {
namespace matrix {

class Uint64BuffWrite {    
public:
    template<class Store>
    void add(uint64_t v, uint8_t n_bits, Store store) {        
        if (pos_ == 64) {
            store(buffer_);
            buffer_ = v;
            pos_ = n_bits;
            return;
        }
        buffer_ += v << pos_;
        auto free_bits = 64 - pos_;
        if (free_bits >= n_bits) {            
            pos_ += n_bits;
        }
        else {
            store(buffer_);
            buffer_ = v >> free_bits;
            pos_ = n_bits - free_bits;
        }
    }
    template<class Store>
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

class Uint64BuffRead {
public:
    template<class Load>
    void get(uint64_t& v, uint8_t n_bits, Load load) {
        if (pos_ == 64) {
            load(buffer_);
            pos_ = 0;
        }
        auto have = 64 - pos_;
        if (have >= n_bits) {
            v = buffer_ >> pos_;
            if(n_bits != 64)
                v &= (1ull << n_bits) - 1;
            pos_ += n_bits;
        }
        else {
            v = buffer_ >> pos_;
            load(buffer_);
            v += buffer_ << have;
            if (n_bits != 64)
                v &= (1ull << n_bits) - 1;
            pos_ = n_bits - have;
        }
    }

    void set_state(uint64_t buffer, uint8_t pos) {
        buffer_ = buffer;
        pos_ = pos;
    }

private:
    uint64_t buffer_;
    uint8_t pos_ = 64;
};


class DiskWriter {
public:
    DiskWriter(std::ofstream& out, size_t buff_size) : out_(out) {
        buff_.reserve(buff_size);
    }
    void add(uint64_t v, uint8_t n_bits) {
        buff_write_.add(v, n_bits, [&](uint64_t x) {
            add(x);
        });
    }
    void flush() {
        buff_write_.flush([&](uint64_t x) {
            add(x);
        });
        if (buff_.size()) {            
            store_in_file();
        }
    }
    ~DiskWriter() {
        flush();
    }

private:
    std::vector<uint64_t> buff_;
    Uint64BuffWrite buff_write_;
    std::ofstream& out_;

    void store_in_file() {
        for (auto x : buff_)
            serialize_number(out_, x);

        buff_.clear();
    }

    void add(uint64_t x) {
        buff_.push_back(x);
        if (buff_.size() == buff_.capacity())
            store_in_file();
    }
};


class DiskRandomReader {
public:
    DiskRandomReader(std::ifstream& in, size_t buff_size) : in_(in){
        offset_ = in_.tellg();

        in_.seekg(0, std::ios::end);

        file_size_ = in.tellg();

        in_.seekg(offset_, std::ios::beg);
        buffer_start_pos_in_file_ = offset_;
        buffer_.resize(buff_size);
        load(buffer_start_pos_in_file_);
    }
    void start_reading_at(uint64_t bit_pos) {
        read_pos_ = bit_pos / 64;
        reaload_if_needed();
        buff_read_.set_state(buffer_[read_pos_ - (buffer_start_pos_in_file_ - offset_) / 8], bit_pos % 64);
        read_pos_++;
    }

    void get(uint64_t& v, uint8_t n_bits) {
        buff_read_.get(v, n_bits, [&](uint64_t& x) {
            reaload_if_needed();
            x = buffer_[read_pos_ - (buffer_start_pos_in_file_ - offset_) / 8];
            read_pos_++;
        });
    }
private:    
    std::vector<uint64_t> buffer_;
    size_t read_pos_{}; //global in uint64_t
    size_t buffer_start_pos_in_file_{}; //in bytes
    std::ifstream& in_;
    size_t offset_;
    Uint64BuffRead buff_read_;
    size_t file_size_;

    void load(size_t pos) {
        size_t uint64_t_elems_to_load = buffer_.size();
        if (file_size_ - pos < uint64_t_elems_to_load * 8)
            uint64_t_elems_to_load = (file_size_ - pos) / 8;
        
        for (size_t i = 0; i < uint64_t_elems_to_load; ++i)
            buffer_[i] = load_number(in_);
    }

    void reaload_if_needed() {
        auto read_pos_in_file = offset_ + read_pos_ * 8;
        if (read_pos_in_file < buffer_start_pos_in_file_ ||
            read_pos_in_file >= buffer_start_pos_in_file_ + 8 * buffer_.size()) {
            in_.seekg(read_pos_in_file, std::ios::beg);
            buffer_start_pos_in_file_ = read_pos_in_file;
            load(buffer_start_pos_in_file_);
        }
    }    
};

} // namespace matrix
} // namespace annot
} // namespace mtg

#endif // __DISK_UTILS__