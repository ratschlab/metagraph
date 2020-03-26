#include "elias_fano.hpp"

#include <algorithm>
#include <cassert>
#include <filesystem>

#include <sdsl/uint128_t.hpp>


namespace mg {
namespace common {

template <class T, class Enable = void>
struct Unaligned;

/**
 * Representation of an unaligned value of a POD type.
 */
template <class T>
struct Unaligned<T, typename std::enable_if<std::is_pod<T>::value>::type> {
    Unaligned() = default; // uninitialized
    /* implicit */ Unaligned(T v) : value(v) {}
    T value;
} __attribute__((__packed__));

/**
 * Read an unaligned value of type T and return it.
 */
template <class T>
inline T load_unaligned(const void *p) {
    static_assert(sizeof(Unaligned<T>) == sizeof(T), "Invalid unaligned size");
    static_assert(alignof(Unaligned<T>) == 1, "Invalid alignment");
    return static_cast<const Unaligned<T> *>(p)->value;
}

/**
 * Inform the compiler that the argument can be assumed true. It is
 * undefined behavior if the argument is not actually true, so use
 * with care.
 *
 * Implemented as a function instead of a macro because
 * __builtin_assume does not evaluate its argument at runtime, so it
 * cannot be used with expressions that have side-effects.
 */
inline __attribute__((__always_inline__)) void assume(bool cond) {
#if defined(__clang__) // Must go first because Clang also defines __GNUC__.
    __builtin_assume(cond);
#elif defined(__GNUC__)
    if (!cond) {
        __builtin_unreachable();
    }
#endif
}

/**
 * Write an unaligned value of type T.
 */
template <class T>
inline void store_unaligned(void *p, T value) {
    static_assert(sizeof(Unaligned<T>) == sizeof(T), "Invalid unaligned size");
    static_assert(alignof(Unaligned<T>) == 1, "Invalid alignment");
    // Prior to C++14, the spec says that a placement new like this
    // is required to check that p is not nullptr, and to do nothing
    // if p is a nullptr. By assuming it's not a nullptr, we get a
    // nice loud segfault in optimized builds if p is nullptr, rather
    // than just silently doing nothing.
    assume(p != nullptr);
    new (p) Unaligned<T>(value);
}

// ------------ EliasFanoEncoder ----------------
template <typename T>
EliasFanoEncoder<T>::EliasFanoEncoder(size_t size,
                                      T min_value,
                                      T max_value,
                                      const std::string &out_filename,
                                      bool is_append)
    : declared_size_(size), offset_(min_value) {
    std::ios_base::openmode open_flag = is_append ? std::ios::app : std::ios::out;
    sink_internal_ = std::ofstream(out_filename, std::ios::binary | open_flag);
    sink_ = &sink_internal_;
    if (!sink_->good()) {
        std::cerr << "Unable to write to " << out_filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    init(size, max_value);
}

template <typename T>
EliasFanoEncoder<T>::EliasFanoEncoder(const Vector<T> &data, std::ofstream &sink)
    : declared_size_(data.size()), sink_(&sink) {
    if (data.size() == 0) {
        return;
    }
    offset_ = data.front();
    init(data.size(), data.back());
    for (const auto &v : data) {
        add(v);
    }
}

template <typename T>
void EliasFanoEncoder<T>::add(T value) {
#ifndef NDEBUG
    assert(value >= last_value_);
#endif
    value -= offset_;
    const T upper_bits = value >> num_lower_bits_;

    // We are adding the size_-th element, so we have a 1 followed by upper_bits
    // zeros, plus the 1s for the previous size_ elements; this is not trivial to
    // understand, so spend some time thinking about why this is correct
    const T pos = upper_bits + size_;
    upper_[pos / 8] |= 1U << (pos % 8);

    // Append the #num_lower_bits_ bits of #value to #lower_
    if (num_lower_bits_ != 0) {
        const T lowerBits = value & lower_bits_mask_;
        size_t pos_bits = size_ * num_lower_bits_;
        if (pos_bits - cur_pos_lbits_ >= 64) { // first 64 bits are ready to be written
            cur_pos_lbits_ += 64;
            sink_->write(reinterpret_cast<char *>(lower_), sizeof(uint64_t));
            lower_[0] = lower_[1];
            lower_[1] = 0;
        }
        write_bits(reinterpret_cast<uint8_t *>(lower_), pos_bits % 64, lowerBits);
    }

#ifndef NDEBUG
    last_value_ = value;
#endif
    ++size_;
}

template <typename T>
size_t EliasFanoEncoder<T>::finish() {
    assert(size_ == declared_size_);
    if (size_ == 0) {
        if (sink_internal_.is_open()) {
            sink_internal_.close();
        }
        return 0;
    }
    // Append the remaining lower bits
    if (num_lower_bits_ != 0) {
        size_t cur_pos_bytes = (cur_pos_lbits_ + 7) / 8;
        assert(cur_pos_bytes <= num_lower_bytes_);
        assert(num_lower_bytes_ - cur_pos_bytes < 16);
        sink_->write(reinterpret_cast<char *>(lower_), num_lower_bytes_ - cur_pos_bytes);
    }
    sink_->write(upper_.data(), num_upper_bytes_);
    if (sink_internal_.is_open()) {
        sink_internal_.close();
    }
    return num_lower_bytes_ + num_upper_bytes_ + sizeof(size_) + sizeof(num_lower_bits_)
            + sizeof(T) + sizeof(num_upper_bytes_) + sizeof(num_lower_bytes_);
}

template <typename T>
void EliasFanoEncoder<T>::init(size_t size, T max_value) {
    if (size == 0) {
        return;
    }
    max_value -= offset_;
    // cap at 56 because #write_bits supports a max of 56 bits
    num_lower_bits_
            = std::min(get_num_lower_bits(max_value, size), static_cast<uint8_t>(56));
    lower_bits_mask_ = (T(1) << num_lower_bits_) - 1;
    // Number of 0-bits to be stored + 1-bits
    const uint64_t upper_size_bits
            = static_cast<uint64_t>(max_value >> num_lower_bits_) + size;
    num_upper_bytes_ = (upper_size_bits + 7) / 8;
    num_lower_bytes_ = (num_lower_bits_ * size + 7) / 8;

    // Current read/write logic assumes that the 7 bytes following the last byte of
    // lower and upper sequences are readable (the stored value doesn't matter and
    // won't be changed), so we reserve an additional 7 bytes for padding
    if (size > 0) {
        upper_.reserve(num_upper_bytes_ + 7);
        upper_.resize(num_upper_bytes_);
    }

    sink_->write(reinterpret_cast<char *>(&size), sizeof(size_t));
    sink_->write(reinterpret_cast<char *>(&offset_), sizeof(T));
    sink_->write(reinterpret_cast<char *>(&num_lower_bits_), 1);
    sink_->write(reinterpret_cast<char *>(&num_lower_bytes_), sizeof(size_t));
    sink_->write(reinterpret_cast<char *>(&num_upper_bytes_), sizeof(size_t));
}

template <typename T>
uint8_t EliasFanoEncoder<T>::get_num_lower_bits(T max_value, size_t size) {
    if (size == 0 || max_value < size) {
        return 0;
    }
    // Result that should be returned is "floor(log(upperBound / size))".
    // In order to avoid expensive division, we rely on
    // "floor(a) - floor(b) - 1 <= floor(a - b) <= floor(a) - floor(b)".
    // Assuming "candidate = floor(log(upperBound)) - floor(log(upperBound))",
    // then result is either "candidate - 1" or "candidate".
    size_t candidate
            = sdsl::bits::hi(static_cast<uint64_t>(max_value)) - sdsl::bits::hi(size);

    return (size > static_cast<uint64_t>(max_value >> candidate)) ? candidate - 1 : candidate;
}

template <typename T>
void EliasFanoEncoder<T>::write_bits(uint8_t *data, size_t pos, uint64_t value) {
    assert(sdsl::bits::hi(value) < 56);
    unsigned char *const ptr = data + (pos / 8);
    uint64_t ptrv = load_unaligned<uint64_t>(ptr);
    ptrv |= value << (pos % 8);
    store_unaligned<uint64_t>(ptr, ptrv);
}


// ------------------------------- EliasFanoDecoder ------------------------------------

template <typename T>
EliasFanoDecoder<T>::EliasFanoDecoder(const std::string &source_name) {
    if (std::filesystem::file_size(source_name) == 0) {
        position_ = size_ = 0;
        all_read_ = true;
        return;
    }
    source_internal_ = std::ifstream(source_name, std::ios::binary);
    source_ = &source_internal_;
    // profiling indicates that setting a larger buffer slightly increases performance
    source_->rdbuf()->pubsetbuf(buffer_, 1024 * 1024);
    if (!source_->good()) {
        std::cerr << "Unable to read from " << source_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
    source_->seekg(0, source_->end);
    file_end_ = source_->tellg();
    source_->seekg(0, source_->beg);
    init();
}


template <typename T>
EliasFanoDecoder<T>::EliasFanoDecoder(std::ifstream &source, std::streampos file_end)
    : source_(&source), file_end_(file_end) {}

template <typename T>
T EliasFanoDecoder<T>::next_upper() {
    // Skip to the first non-zero block.
    while (upper_block_ == 0) {
        upper_pos_ += sizeof(uint64_t);
        upper_block_ = load_unaligned<uint64_t>(upper_.data() + upper_pos_);
    }

    size_t trailing_zeros = __builtin_ctzll(upper_block_); // count trailing zeros
    upper_block_ = upper_block_ & (upper_block_ - 1); // reset the lowest 1 bit

    return static_cast<T>(8 * upper_pos_ + trailing_zeros - position_);
}

template <typename T>
T EliasFanoDecoder<T>::next_lower() {
    assert(position_ < size_);
    const size_t pos_bits = position_ * num_lower_bits_;
    if (pos_bits - cur_pos_bits_ >= 64) {
        // +16 because we read 16 bytes in #init()
        const size_t bytes_read = cur_pos_bits_ / 8 + 16;
        cur_pos_bits_ += 64;
        lower_[0] = lower_[1];
        lower_[1] = 0;
        if (num_lower_bytes_ > bytes_read) {
            source_->read(reinterpret_cast<char *>(&lower_[1]),
                          std::min(sizeof(uint64_t), num_lower_bytes_ - bytes_read));
        }
    }
    const size_t adjusted_pos = pos_bits - cur_pos_bits_;

    const uint8_t *ptr = reinterpret_cast<uint8_t *>(lower_) + (adjusted_pos / 8);
    const uint64_t ptrv = load_unaligned<uint64_t>(ptr);
    return T(clear_high_bits(ptrv >> (adjusted_pos % 8), num_lower_bits_));
}

template <typename T>
std::optional<T> EliasFanoDecoder<T>::next() {
    if (position_ == size_) {
        if (all_read_) {
            return {};
        }
        init(); // read the next chunk of compressed data
    }
    T result = next_lower() | (next_upper() << num_lower_bits_);
    position_++;
    if (end_of_chunk()) {
        // move to the beginning of next chunk
        source_->seekg(num_upper_bytes_, std::ios::cur);
    }
    return result + offset_;
}

template <typename T>
bool EliasFanoDecoder<T>::end_of_chunk() {
    return position_ == size_;
}

// TODO: make this public and avoid reconstruction
template <typename T>
void EliasFanoDecoder<T>::init() {
    position_ = 0;
    cur_pos_bits_ = 0;
    upper_block_ = 0;
    lower_[0] = lower_[1] = 0;
    // Initialized to a negative number to save on decrement instruction in
    // #next_upper.
    upper_pos_ = static_cast<size_t>(-sizeof(size_t));
    source_->read(reinterpret_cast<char *>(&size_), sizeof(size_t));
    source_->read(reinterpret_cast<char *>(&offset_), sizeof(T));
    source_->read(reinterpret_cast<char *>(&num_lower_bits_), 1);
    source_->read(reinterpret_cast<char *>(&num_lower_bytes_), sizeof(size_t));
    source_->read(reinterpret_cast<char *>(&num_upper_bytes_), sizeof(size_t));
    size_t low_bytes_read = std::min(2 * sizeof(uint64_t), num_lower_bytes_);
    source_->read(reinterpret_cast<char *>(lower_), low_bytes_read);

    std::streampos pos = source_->tellg();
    // to avoid jumping through the file, we read the relatively small upper_bytes
    // into memory;
    source_->seekg(num_lower_bytes_ - low_bytes_read, std::ios::cur);
    upper_.reserve(num_upper_bytes_ + 7);
    upper_.resize(num_upper_bytes_);
    source_->read(upper_.data(), num_upper_bytes_);
    assert(static_cast<uint32_t>(source_->gcount()) == num_upper_bytes_);
    if (source_->tellg() == file_end_) {
        all_read_ = true;
    }
    source_->seekg(pos, source_->beg);
}

template <typename T>
uint64_t EliasFanoDecoder<T>::clear_high_bits(uint64_t value, uint32_t index) {
    assert(index < 64);
    return value & ((uint64_t(1) << index) - 1);
}

// -------------------------- EliasFanoEncoder<std::pair> -------------------------------
template <typename T, typename C>
EliasFanoEncoder<std::pair<T, C>>::EliasFanoEncoder(size_t size,
                                                    const T &first_value,
                                                    const T &last_value,
                                                    const std::string &sink_name,
                                                    bool is_append)
    : ef_encoder(size, first_value, last_value, sink_name, is_append),
      sink_second_name_(sink_name + ".count") {
    std::ios_base::openmode open_flag = is_append ? std::ios::app : std::ios::out;
    sink_second_ = std::ofstream(sink_second_name_, std::ios::binary | open_flag);
}


template <typename T, typename C>
void EliasFanoEncoder<std::pair<T, C>>::add(const std::pair<T, C> &value) {
    ef_encoder.add(value.first);
    sink_second_.write(reinterpret_cast<const char *>(&value.second), sizeof(C));
}

template <typename T, typename C>
size_t EliasFanoEncoder<std::pair<T, C>>::finish() {
    size_t first_size = ef_encoder.finish();
    sink_second_.close();
    return first_size + std::filesystem::file_size(sink_second_name_);
}

// ------------------------- EliasFandDecoder<std::pair> --------------------------------
template <typename T, typename C>
EliasFanoDecoder<std::pair<T, C>>::EliasFanoDecoder(const std::string &source)
    : source_first_(source) {
    std::string source_second_name = source + ".count";
    source_second_ = std::ifstream(source_second_name, std::ios::binary);
    // profiling indicates that setting a larger buffer slightly increases performance
    source_second_.rdbuf()->pubsetbuf(buffer_, 1024 * 1024);
    if (!source_second_.good()) {
        std::cerr << "Unable to write to " << source_second_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

template <typename T, typename C>
std::optional<std::pair<T, C>> EliasFanoDecoder<std::pair<T, C>>::next() {
    std::optional<T> first = source_first_.next();
    C second;
    if (!first.has_value()) {
        assert(!source_second_.read(reinterpret_cast<char *>(&second), sizeof(C)));
        return {};
    }
    source_second_.read(reinterpret_cast<char *>(&second), sizeof(C));
    assert(source_second_);
    return std::make_pair(first.value(), second);
}

// ---------------------- EliasFandEncoder<sdsl::uint128_t> -----------------------------
EliasFanoEncoder<sdsl::uint128_t>::EliasFanoEncoder(size_t,
                                                    sdsl::uint128_t,
                                                    sdsl::uint128_t,
                                                    const std::string &sink_name,
                                                    bool is_append) {
    std::ios_base::openmode open_flag = is_append ? std::ios::app : std::ios::out;
    // open file for appending, as we may encode multiple compressed chunks in the same file
    sink_internal_ = std::ofstream(sink_name, std::ios::binary | open_flag);
    sink_ = &sink_internal_;
    if (!sink_->good()) {
        std::cerr << "Unable to write to " << sink_name << std::endl;
        std::exit(EXIT_FAILURE);
    }

}

EliasFanoEncoder<sdsl::uint128_t>::EliasFanoEncoder(const Vector<sdsl::uint128_t> &data,
                                                    std::ofstream &sink) : sink_(&sink) {
    for (const auto &v : data) {
        add(v);
    }
}


void EliasFanoEncoder<sdsl::uint128_t>::add(const sdsl::uint128_t &value) {
    uint64_t current_hi = static_cast<uint64_t>(value >> 64);
    if (buffer_.empty()) {
        last_hi_ = current_hi;
    }
    if (current_hi == last_hi_) {
        buffer_.push_back((uint64_t)value);
    } else {
        sink_->write(reinterpret_cast<const char *>(&last_hi_), sizeof(uint64_t));
        EliasFanoEncoder<uint64_t> encoder64_ = EliasFanoEncoder<uint64_t>(buffer_, *sink_);
        total_size_ += (sizeof(uint64_t) + encoder64_.finish());
        buffer_.resize(0);
        buffer_.push_back(value);
        last_hi_ = current_hi;
    }
    size_++;
}

size_t EliasFanoEncoder<sdsl::uint128_t>::finish() {
    if (!buffer_.empty()) {
        sink_->write(reinterpret_cast<const char *>(&last_hi_), sizeof(uint64_t));
        EliasFanoEncoder<uint64_t> encoder64_ = EliasFanoEncoder<uint64_t>(buffer_, *sink_);
        total_size_ += (sizeof(uint64_t) + encoder64_.finish());
    }
    sink_->close();
    return total_size_;
}

// --------------------- EliasFanoEncoder<sdsl::uint256_t> ------------------------------
EliasFanoEncoder<sdsl::uint256_t>::EliasFanoEncoder(size_t,
                                                    sdsl::uint256_t,
                                                    sdsl::uint256_t,
                                                    const std::string &sink_name,
                                                    bool is_append) {
    std::ios_base::openmode open_flag = is_append ? std::ios::app : std::ios::out;
    // open file for appending, as we may encode multiple compressed chunks in the same file
    sink_ = std::ofstream(sink_name, std::ios::binary | open_flag);
    if (!sink_.good()) {
        std::cerr << "Unable to write to " << sink_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

EliasFanoEncoder<sdsl::uint256_t>::EliasFanoEncoder(const Vector<sdsl::uint256_t> &data,
                                                    std::ofstream &sink) {
    sink.write(reinterpret_cast<const char *>(data.data()),
               data.size() * sizeof(sdsl::uint256_t));
    total_size_ += data.size() * sizeof(sdsl::uint256_t);
}

void EliasFanoEncoder<sdsl::uint256_t>::add(const sdsl::uint256_t &value) {
    sink_.write(reinterpret_cast<const char *>(&value), sizeof(sdsl::uint256_t));
    total_size_ += sizeof(sdsl::uint256_t);
    size_++;
}

size_t EliasFanoEncoder<sdsl::uint256_t>::finish() {
    sink_.close();
    return total_size_;
}

// ---------------------- EliasFanoDecoder<sdsl::uint128_t> -----------------------------
EliasFanoDecoder<sdsl::uint128_t>::EliasFanoDecoder(const std::string &source_name) {
    source_ = std::ifstream(source_name, std::ios::binary);
    source_.rdbuf()->pubsetbuf(buffer_, 1024 * 1024);
    if (!source_.good()) {
        std::cerr << "Unable to read from " << source_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
    source_.seekg(0, source_.end);
    std::streampos file_end = source_.tellg();
    source_.seekg(0, source_.beg);
    decoder64_ = EliasFanoDecoder<uint64_t>(source_, file_end);
}

std::optional<sdsl::uint128_t> EliasFanoDecoder<sdsl::uint128_t>::next() {
    if (new_chunk_) {
        if (!source_.read(reinterpret_cast<char *>(&last_hi_), sizeof(uint64_t))) {
            assert(!decoder64_.next().has_value());
            return {};
        }
    }
    std::optional<uint64_t> next_lo_ = decoder64_.next();
    assert(next_lo_.has_value());
    new_chunk_ = decoder64_.end_of_chunk();
    return ((sdsl::uint128_t)last_hi_ << 64) + next_lo_.value();
}

/** Creates a decoder that retrieves data from the given source */
EliasFanoDecoder<sdsl::uint256_t>::EliasFanoDecoder(const std::string &source_name) {
    source_ = std::ifstream(source_name, std::ios::binary);
    if (!source_.good()) {
        std::cerr << "Unable to read from " << source_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

/**
 * Returns the next compressed element or empty if all elements were read.
 */
std::optional<sdsl::uint256_t> EliasFanoDecoder<sdsl::uint256_t>::next() {
    sdsl::uint256_t result;
    if (source_.read(reinterpret_cast<char *>(&result), sizeof(sdsl::uint256_t))) {
        return result;
    }
    return {};
}


// ------------------------------ EliasFanoEncoderBuffered ----------------------------
template <typename T>
EliasFanoEncoderBuffered<T>::EliasFanoEncoderBuffered(const std::string &file_name,
                                                      size_t buffer_size)
    : file_name_(file_name) {
    std::filesystem::remove(file_name);
    buffer_.reserve(buffer_size);
}

template <typename T>
void EliasFanoEncoderBuffered<T>::add(const T &value) {
    buffer_.push_back(value);
    if (buffer_.size() == buffer_.capacity()) {
        encode_chunk();
    }
}

template <typename T>
size_t EliasFanoEncoderBuffered<T>::finish() {
    encode_chunk();
    return total_size_;
}

template <typename T>
void EliasFanoEncoderBuffered<T>::encode_chunk() {
    EliasFanoEncoder<T> encoder_
            = EliasFanoEncoder<T>(buffer_.size(), buffer_.empty() ? T(0) : buffer_.front(),
                                  buffer_.empty() ? T(0) : buffer_.back(), file_name_, true);
    for (const auto &v : buffer_) {
        encoder_.add(v);
    }
    total_size_ += encoder_.finish();
    buffer_.resize(0);
}

// ----------------------- EliasFanoEncoderBuffered<std::pair> --------------------------
template <typename T, typename C>
EliasFanoEncoderBuffered<std::pair<T, C>>::EliasFanoEncoderBuffered(const std::string &file_name,
                                                                    size_t buffer_size)
    : file_name_(file_name) {
    std::filesystem::remove(file_name);
    sink_ = std::ofstream(file_name, std::ios::binary | std::ios::out);
    sink_second_ = std::ofstream(file_name + ".count", std::ios::binary | std::ios::out);
    buffer_.reserve(buffer_size);
    buffer_second_.reserve(buffer_size);
}

template <typename T, typename C>
void EliasFanoEncoderBuffered<std::pair<T, C>>::add(const std::pair<T, C> &value) {
    buffer_.push_back(value.first);
    buffer_second_.push_back(value.second);
    if (buffer_.size() == buffer_.capacity()) {
        encode_chunk();
    }
}

template <typename T, typename C>
size_t EliasFanoEncoderBuffered<std::pair<T, C>>::finish() {
    encode_chunk();
    sink_.close();
    sink_second_.close();
    return total_size_;
}

template <typename T, typename C>
void EliasFanoEncoderBuffered<std::pair<T, C>>::encode_chunk() {
    EliasFanoEncoder<T> encoder_ = EliasFanoEncoder<T>(buffer_, sink_);

    sink_second_.write(reinterpret_cast<char *>(buffer_second_.data()),
                       buffer_second_.size() * sizeof(C));
    total_size_ += (encoder_.finish() + buffer_second_.size() * sizeof(C));
    buffer_.resize(0);
    buffer_second_.resize(0);
}

// instantiate used templates
template class EliasFanoEncoder<uint32_t>;
template class EliasFanoEncoder<uint64_t>;
template class EliasFanoEncoder<std::pair<uint64_t, uint8_t>>;
template class EliasFanoEncoder<std::pair<uint64_t, uint16_t>>;
template class EliasFanoEncoder<std::pair<uint64_t, uint32_t>>;
template class EliasFanoEncoder<std::pair<uint32_t, uint8_t>>;
template class EliasFanoEncoder<std::pair<uint32_t, uint16_t>>;
template class EliasFanoEncoder<std::pair<uint32_t, uint32_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint128_t, uint8_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint128_t, uint16_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint128_t, uint32_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint256_t, uint8_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint256_t, uint16_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint256_t, uint32_t>>;


template class EliasFanoDecoder<uint32_t>;
template class EliasFanoDecoder<uint64_t>;
template class EliasFanoDecoder<std::pair<uint32_t, uint8_t>>;
template class EliasFanoDecoder<std::pair<uint32_t, uint16_t>>;
template class EliasFanoDecoder<std::pair<uint32_t, uint32_t>>;
template class EliasFanoDecoder<std::pair<uint64_t, uint8_t>>;
template class EliasFanoDecoder<std::pair<uint64_t, uint16_t>>;
template class EliasFanoDecoder<std::pair<uint64_t, uint32_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint128_t, uint8_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint128_t, uint16_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint128_t, uint32_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint256_t, uint8_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint256_t, uint16_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint256_t, uint32_t>>;

template class EliasFanoEncoderBuffered<uint32_t>;
template class EliasFanoEncoderBuffered<uint64_t>;
template class EliasFanoEncoderBuffered<sdsl::uint128_t>;
template class EliasFanoEncoderBuffered<sdsl::uint256_t>;
template class EliasFanoEncoderBuffered<std::pair<uint32_t, uint8_t>>;
template class EliasFanoEncoderBuffered<std::pair<uint32_t, uint16_t>>;
template class EliasFanoEncoderBuffered<std::pair<uint32_t, uint32_t>>;
template class EliasFanoEncoderBuffered<std::pair<uint64_t, uint8_t>>;
template class EliasFanoEncoderBuffered<std::pair<uint64_t, uint16_t>>;
template class EliasFanoEncoderBuffered<std::pair<uint64_t, uint32_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint128_t, uint8_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint128_t, uint16_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint128_t, uint32_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint256_t, uint8_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint256_t, uint16_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint256_t, uint32_t>>;

} // namespace common
} // namespace mg
