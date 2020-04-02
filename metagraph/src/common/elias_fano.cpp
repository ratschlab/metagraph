#include "elias_fano.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <filesystem>

#include <sdsl/uint128_t.hpp>

#include <common/utils/template_utils.hpp>


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

/** Template specialization for an unaligned sdsl::uint256_t, which is not a POD */
template <>
struct Unaligned<sdsl::uint256_t> {
    Unaligned() = default; // uninitialized
    /* implicit */ Unaligned(sdsl::uint256_t v) : lo(v), hi(v >> 128) {}
    sdsl::uint128_t lo; // sdsl::uint256 is not a POD, so needs to be stored as 2 uint128_t's
    sdsl::uint128_t hi;
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
// TODO(ddanciu): define load_unaligned when MODE_TI is not present
template <>
inline sdsl::uint256_t load_unaligned(const void *p) {
    static_assert(sizeof(Unaligned<sdsl::uint256_t>) == sizeof(sdsl::uint256_t),
                  "Invalid unaligned size");
    static_assert(alignof(Unaligned<sdsl::uint256_t>) == 1, "Invalid alignment");
    return sdsl::uint256_t(static_cast<const Unaligned<sdsl::uint256_t> *>(p)->lo,
                           static_cast<const Unaligned<sdsl::uint256_t> *>(p)->hi);
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

/** Returns the rank of the highest non-null bit */
template <typename T>
inline uint32_t log2_floor(T x) {
    return sdsl::bits::hi(x);
}

template <>
inline uint32_t log2_floor(sdsl::uint128_t x) {
    if (x == 0U) {
        return 0;
    }

    uint64_t hi = static_cast<uint64_t>(x >> 64);
    if (hi == 0U) {
        uint64_t lo = static_cast<uint64_t>(x);
        return 63 - __builtin_clzll(lo);
    } else {
        return 127 - __builtin_clzll(hi);
    }
}

template <>
inline uint32_t log2_floor(sdsl::uint256_t x) {
    if (x == 0U) {
        return 0;
    }

    sdsl::uint128_t hi = static_cast<sdsl::uint128_t>(x >> 128);
    if (hi == 0U) {
        return log2_floor(static_cast<sdsl::uint128_t>(x));
    } else {
        return 128 + log2_floor(hi);
    }
}

/** Returns the trailing zeros in the binary representation of a number */
template <typename T>
inline uint32_t count_trailing_zeros(T v) {
    assert(v != 0U);
    return __builtin_ctzll(v);
}


template <>
inline uint32_t count_trailing_zeros(sdsl::uint128_t v) {
    assert(v != 0U);
    uint64_t lo = static_cast<uint64_t>(v);
    if (lo != 0) {
        return __builtin_ctzll(lo);
    } else {
        uint64_t hi = static_cast<uint64_t>(v >> 64);
        return 64 + __builtin_ctzll(hi);
    }
}

template <>
inline uint32_t count_trailing_zeros(sdsl::uint256_t v) {
    assert(v != 0U);
    sdsl::uint128_t lo = static_cast<sdsl::uint128_t>(v);
    if (lo != 0U) {
        return count_trailing_zeros(lo);
    } else {
        sdsl::uint128_t hi = static_cast<sdsl::uint128_t>(v >> 128);
        return 128 + count_trailing_zeros(hi);
    }
}

static void safe_close(std::ofstream &sink) {
    if (sink.is_open()) {
        sink.close();
    }
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
    sink_internal_upper_ = std::ofstream(out_filename + ".up", std::ios::binary | open_flag);
    sink_upper_ = &sink_internal_upper_;
    if (!sink_->good() || !sink_upper_->good()) {
        std::cerr << "Unable to write to " << out_filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    init(size, max_value);
}

template <typename T>
EliasFanoEncoder<T>::EliasFanoEncoder(const Vector<T> &data,
                                      std::ofstream *sink,
                                      std::ofstream *sink_upper)
    : declared_size_(data.size()), sink_(sink), sink_upper_(sink_upper) {
    if (data.size() == 0U) {
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
    // using ">> 3" for /8 and "&7" for %8 because sdsl::uint256_t doesn't define / and %
    upper_[(uint64_t)(pos >> 3)] |= 1U << (uint64_t)(pos & 7);

    // Append the #num_lower_bits_ bits of #value to #lower_
    if (num_lower_bits_ != 0) {
        const T lowerBits = value & lower_bits_mask_;
        size_t pos_bits = size_ * num_lower_bits_;
        if (pos_bits - cur_pos_lbits_ >= 8 * sizeof(T)) {
            // first sizeof(T)*8 bits are ready to be written
            cur_pos_lbits_ += 8 * sizeof(T);
            sink_->write(reinterpret_cast<char *>(lower_), sizeof(T));
            lower_[0] = lower_[1];
            lower_[1] = 0;
        }
        write_bits(reinterpret_cast<uint8_t *>(lower_), pos_bits % (8 * sizeof(T)),
                   lowerBits);
    }

#ifndef NDEBUG
    last_value_ = value;
#endif
    ++size_;
}

template <typename T>
size_t EliasFanoEncoder<T>::finish() {
    assert(size_ == declared_size_);
    if (size_ == 0U) {
        safe_close(sink_internal_);
        safe_close(sink_internal_upper_);
        return 0;
    }
    // Append the remaining lower bits
    if (num_lower_bits_ != 0) {
        size_t cur_pos_bytes = (cur_pos_lbits_ + 7) / 8;
        assert(cur_pos_bytes <= num_lower_bytes_);
        assert(num_lower_bytes_ - cur_pos_bytes < 2 * sizeof(T));
        sink_->write(reinterpret_cast<char *>(lower_), num_lower_bytes_ - cur_pos_bytes);
    }
    sink_upper_->write(upper_.data(), num_upper_bytes_);
    safe_close(sink_internal_);
    safe_close(sink_internal_upper_);
    return num_lower_bytes_ + num_upper_bytes_ + sizeof(size_) + sizeof(num_lower_bits_)
            + sizeof(T) + sizeof(num_upper_bytes_) + sizeof(num_lower_bytes_);
}

template <typename T>
void EliasFanoEncoder<T>::init(size_t size, T max_value) {
    if (size == 0U) {
        return;
    }
    max_value -= offset_;
    // #write_bits supports a max of sizeof(T)-1 bytes
    num_lower_bits_ = std::min(get_num_lower_bits(max_value, size),
                               static_cast<uint8_t>(8 * sizeof(T) - 8));
    lower_bits_mask_ = (T(1) << num_lower_bits_) - 1UL;
    // Number of 0-bits to be stored + 1-bits
    const uint64_t upper_size_bits
            = static_cast<uint64_t>(max_value >> num_lower_bits_) + size;
    num_upper_bytes_ = (upper_size_bits + 7) / 8;
    num_lower_bytes_ = (num_lower_bits_ * size + 7) / 8;

    // Current read/write logic assumes that the sizeof(T)-1 bytes following the last byte of
    // lower and upper sequences are readable (the stored value doesn't matter and
    // won't be changed), so we reserve an additional 7 bytes for padding
    if (size > 0) {
        upper_.reserve(num_upper_bytes_ + sizeof(T) - 1);
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
    size_t candidate = log2_floor(max_value) - log2_floor(size);

    return (size > static_cast<uint64_t>(max_value >> candidate)) ? candidate - 1 : candidate;
}

template <typename T>
void EliasFanoEncoder<T>::write_bits(uint8_t *data, size_t pos, T value) {
    unsigned char *const ptr = data + (pos / 8);
    if constexpr (sizeof(T) >= 16) {
        assert(log2_floor(value) < 8 * (sizeof(T) - 1));
        T ptrv = load_unaligned<T>(ptr);
        ptrv |= value << (pos % 8);
        store_unaligned<T>(ptr, ptrv);
    } else { // all types <=64 bits are stored in 64-bit chunks
        assert(log2_floor(value) < 56);
        uint64_t ptrv = load_unaligned<uint64_t>(ptr);
        ptrv |= value << (pos % 8);
        store_unaligned<uint64_t>(ptr, ptrv);
    }
}


// ------------------------------- EliasFanoDecoder ------------------------------------

template <typename T>
EliasFanoDecoder<T>::EliasFanoDecoder(const std::string &source_name) {
    source_internal_ = std::ifstream(source_name, std::ios::binary);
    source_ = &source_internal_;
    source_internal_upper_ = std::ifstream(source_name + ".up", std::ios::binary);
    source_upper_ = &source_internal_upper_;
    if (!source_->good() || !source_upper_->good()) {
        std::cerr << "Unable to read from " << source_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
    file_end_ = std::filesystem::file_size(source_name + ".up");
    init();
}


template <typename T>
EliasFanoDecoder<T>::EliasFanoDecoder(std::ifstream *source,
                                      std::ifstream *source_upper,
                                      std::streampos file_end_upper)
    : source_(source), source_upper_(source_upper), file_end_(file_end_upper) {}

template <typename T>
T EliasFanoDecoder<T>::next_upper() {
    // Skip to the first non-zero block.
    while (upper_block_ == 0U) {
        upper_pos_ += sizeof(T);
        upper_block_ = load_unaligned<T>(upper_.data() + upper_pos_);
    }

    size_t trailing_zeros = count_trailing_zeros(upper_block_);
    upper_block_ = upper_block_ & (upper_block_ - 1UL); // reset the lowest 1 bit

    return static_cast<T>(8 * upper_pos_ + trailing_zeros - position_);
}

template <typename T>
T EliasFanoDecoder<T>::next_lower() {
    assert(position_ < size_);
    const size_t pos_bits = position_ * num_lower_bits_;
    if (pos_bits - cur_pos_bits_ >= 8 * sizeof(T)) {
        cur_pos_bits_ += 8 * sizeof(T);
        lower_idx_++;
        if (lower_idx_ == READ_BUF_SIZE-1 && num_lower_bytes_ > 0) {
            const uint32_t to_read
                    = std::min(sizeof(lower_) - sizeof(T), num_lower_bytes_);
            lower_[0] = lower_[lower_idx_];
            source_->read(reinterpret_cast<char *>(&lower_[1]), to_read);
            num_lower_bytes_ -= to_read;
            lower_idx_ = 0;
        }
    }
    const size_t adjusted_pos = pos_bits - cur_pos_bits_;

    const uint8_t *ptr
            = reinterpret_cast<uint8_t *>(&lower_[lower_idx_]) + (adjusted_pos / 8);
    const T ptrv = load_unaligned<T>(ptr);
    return T(clear_high_bits(ptrv >> (adjusted_pos % 8), num_lower_bits_));
}

template <typename T>
std::optional<T> EliasFanoDecoder<T>::next() {
    if (position_ == size_) {
        if (!init()) { // read the next chunk of compressed data
            return {};
        }
    }
    T result = next_lower() | (next_upper() << num_lower_bits_);
    position_++;
    return result + offset_;
}

template <typename T>
bool EliasFanoDecoder<T>::end_of_chunk() {
    return position_ == size_;
}

// TODO: make this public and avoid reconstruction
template <typename T>
bool EliasFanoDecoder<T>::init() {
    position_ = 0;
    cur_pos_bits_ = 0;
    upper_block_ = 0;
    memset(lower_, 0, sizeof(lower_));
    // Initialized to a negative number to save on decrement instruction in
    // #next_upper.
    upper_pos_ = static_cast<uint64_t>(-sizeof(T));
    if (!source_->read(reinterpret_cast<char *>(&size_), sizeof(size_t))) {
        assert(!source_upper_->read(upper_.data(), 1));
        return false;
    }
    source_->read(reinterpret_cast<char *>(&offset_), sizeof(T));
    source_->read(reinterpret_cast<char *>(&num_lower_bits_), 1);
    source_->read(reinterpret_cast<char *>(&num_lower_bytes_), sizeof(size_t));
    source_->read(reinterpret_cast<char *>(&num_upper_bytes_), sizeof(size_t));
    size_t low_bytes_read = std::min(sizeof(lower_), num_lower_bytes_);
    source_->read(reinterpret_cast<char *>(lower_), low_bytes_read);
    num_lower_bytes_ -= low_bytes_read;

    upper_.reserve(num_upper_bytes_ + sizeof(T) - 1);
#ifndef NDEBUG
    // silence Valgrind uninitialized warnings. Because our reads are unaligned, we read
    // (and correctly ignore) from uninitilized memory
    memset(upper_.data(), 0, num_upper_bytes_ + sizeof(T) - 1);
#endif
    upper_.resize(num_upper_bytes_);
    source_upper_->read(upper_.data(), num_upper_bytes_);
    assert(static_cast<uint32_t>(source_upper_->gcount()) == num_upper_bytes_);

    return true;
}

template <typename T>
T EliasFanoDecoder<T>::clear_high_bits(T value, uint32_t index) {
    assert(index < 8 * sizeof(T));
    return value & ((T(1) << index) - 1UL);
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
    sink_upper_ = std::ofstream(file_name + ".up", std::ios::binary | std::ios::out);
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
    sink_upper_.close();
    sink_second_.close();
    return total_size_;
}

template <typename T, typename C>
void EliasFanoEncoderBuffered<std::pair<T, C>>::encode_chunk() {
    EliasFanoEncoder<T> encoder_ = EliasFanoEncoder<T>(buffer_, &sink_, &sink_upper_);

    sink_second_.write(reinterpret_cast<char *>(buffer_second_.data()),
                       buffer_second_.size() * sizeof(C));
    total_size_ += (encoder_.finish() + buffer_second_.size() * sizeof(C));
    buffer_.resize(0);
    buffer_second_.resize(0);
}

// instantiate used templates
template class EliasFanoEncoder<uint32_t>;
template class EliasFanoEncoder<uint64_t>;
template class EliasFanoEncoder<sdsl::uint128_t>;
template class EliasFanoEncoder<sdsl::uint256_t>;
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
template class EliasFanoDecoder<sdsl::uint128_t>;
template class EliasFanoDecoder<sdsl::uint256_t>;

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
