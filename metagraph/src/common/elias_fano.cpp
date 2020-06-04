#include "elias_fano.hpp"
#include "logger.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <filesystem>

#include <sdsl/uint128_t.hpp>

#include <common/utils/template_utils.hpp>


namespace mtg {
namespace common {

void concat(const std::vector<std::string> &files, const std::string &result) {
    if (files.empty())
        return;

    std::vector<std::string> suffixes = { "", ".up" };
    if (std::filesystem::exists(files[0] + ".count"))
        suffixes.push_back(".count");

    for (const std::string &suffix : suffixes) {
        std::string concat_command = "cat ";
        for (uint32_t i = 1; i < files.size(); ++i) {
            concat_command += files[i] + suffix + " ";
        }
        concat_command += " >> " + files[0] + suffix;

        if (std::system(concat_command.c_str()))
            throw std::runtime_error("Error while cat-ing files: " + concat_command);

        std::filesystem::rename(files[0] + suffix, result + suffix);
        for (const std::string &f : files) {
            std::filesystem::remove(f + suffix);
        }
    }
}

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
#ifndef MODE_TI
// need to specialize for sdsl::uint128_t which is not a POD if MODE_TI is not available
template <>
struct Unaligned<sdsl::uint128_t> {
    Unaligned() = default; // uninitialized
    /* implicit */ Unaligned(sdsl::uint128_t v) : lo(v), hi(v >> 64) {}
    uint64_t lo; // sdsl::uint128 is not a POD, so needs to be stored as 2 uint64_t's
    uint64_t hi;
} __attribute__((__packed__));

template <>
inline sdsl::uint128_t load_unaligned(const void *p) {
    static_assert(sizeof(Unaligned<sdsl::uint128_t>) == sizeof(sdsl::uint128_t),
                  "Invalid unaligned size");
    static_assert(alignof(Unaligned<sdsl::uint128_t>) == 1, "Invalid alignment");
    return sdsl::uint128_t(static_cast<const Unaligned<sdsl::uint128_t> *>(p)->lo,
                           static_cast<const Unaligned<sdsl::uint128_t> *>(p)->hi);
}
#endif

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
inline uint32_t log2_floor(const T &x) {
    return sdsl::bits::hi(x);
}

template <>
inline uint32_t log2_floor(const sdsl::uint128_t &x) {
#ifdef MODE_TI
    if (x == 0U) {
        return 0;
    }

    uint64_t hi = static_cast<uint64_t>(x >> 64);
    if (hi == 0U) {
        return sdsl::bits::hi(static_cast<uint64_t>(x));
    } else {
        return sdsl::bits::hi(hi) + 64;
    }
#else
    return x.hi();
#endif
}

template <>
inline uint32_t log2_floor(const sdsl::uint256_t &x) {
    return x.hi();
}

static void safe_close(std::ofstream &sink) {
    if (sink.is_open()) {
        sink.close();
    }
}

// ------------ EliasFanoEncoder --------------------------------------------------------
template <typename T>
EliasFanoEncoder<T>::EliasFanoEncoder(size_t size,
                                      T min_value,
                                      T max_value,
                                      const std::string &out_filename)
    : declared_size_(size), offset_(min_value) {
    sink_internal_ = std::ofstream(out_filename, std::ios::binary);
    sink_ = &sink_internal_;
    sink_internal_upper_ = std::ofstream(out_filename + ".up", std::ios::binary);
    sink_upper_ = &sink_internal_upper_;
    if (!sink_->good() || !sink_upper_->good()) {
        logger->error("Unable to write to {}", out_filename);
        std::exit(EXIT_FAILURE);
    }
    init(size, max_value);
}

template <typename T>
EliasFanoEncoder<T>::EliasFanoEncoder(const std::vector<T> &data,
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
EliasFanoEncoder<T>::~EliasFanoEncoder() {
    assert(!sink_internal_.is_open());
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
    assert(upper_bits + size_ < std::numeric_limits<uint64_t>::max());
    const uint64_t pos = upper_bits + size_;
    // using ">> 3" for /8 and "&7" for %8 because sdsl::uint256_t doesn't define / and %
    upper_[pos >> 3] |= 1U << (pos & 7);

    // Append the #num_lower_bits_ bits of #value to #lower_
    if (num_lower_bits_ != 0) {
        const T lower_bits = value & lower_bits_mask_;
        size_t pos_bits = size_ * num_lower_bits_;
        constexpr size_t max_buf_size = sizeof(lower_) - sizeof(T);
        if (pos_bits - cur_pos_lbits_ >= 8 * max_buf_size) {
            // first sizeof(T)*8 bits were written
            cur_pos_lbits_ += 8 * max_buf_size;
            sink_->write(lower_, max_buf_size);
            std::memcpy(lower_, lower_ + max_buf_size, sizeof(T));
            std::memset(lower_ + sizeof(T), 0, max_buf_size);
        }
        write_bits(lower_ + (pos_bits - cur_pos_lbits_) / 8, pos_bits % 8, lower_bits);
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
        size_t pos_bits = size_ * num_lower_bits_;
        size_t num_bytes = (pos_bits - cur_pos_lbits_ + 7) / 8;
        assert(cur_pos_lbits_ / 8 + num_bytes == num_lower_bytes_);
        sink_->write(lower_, num_bytes);
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
    std::memset(lower_, 0, sizeof(lower_));
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
        upper_.resize(num_upper_bytes_ + sizeof(T) - 1, 0);
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
void EliasFanoEncoder<T>::write_bits(char *data, size_t pos, T value) {
    char *const ptr = data + (pos / 8);
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
EliasFanoDecoder<T>::EliasFanoDecoder(const std::string &source_name, bool remove_source)
    : source_name_(source_name), remove_source_(remove_source) {
    source_.open(source_name, std::ios::binary);
    source_upper_.open(source_name + ".up", std::ios::binary);
    if (!source_.good() || !source_upper_.good()) {
        logger->error("Unable to read from {}", source_name);
        std::exit(EXIT_FAILURE);
    }
    init();
}

template <typename T>
size_t EliasFanoDecoder<T>::decompress_next_block() {
    buffer_pos_ = 0;
    buffer_end_ = 0;

    if (size_ == static_cast<size_t>(-1))
        return 0;

    size_t block_begin = position_;
    size_t block_end = std::min(size_, position_ + sizeof(buffer_) / sizeof(T));

    while (true) {
        if (position_ == block_end) {
            // rollback #position_ and read the upper bits
            buffer_end_ -= (block_end - block_begin);
            position_ = block_begin;
            while (position_ < block_end) {
                // Skip to the first non-zero block.
                while (upper_[upper_pos_] == 0U) {
                    upper_pos_++;
                }
                size_t trailing_zeros = sdsl::bits::lo(upper_[upper_pos_]);
                upper_[upper_pos_] &= (upper_[upper_pos_] - 1UL); // reset the lowest 1 bit
                T upper = 64 * upper_pos_ + trailing_zeros - position_;
                buffer_[buffer_end_] |= (upper << num_lower_bits_);
                buffer_[buffer_end_] += offset_;
                buffer_end_++;
                position_++;
            }
            if (buffer_end_ == sizeof(buffer_) / sizeof(T)) {
                break;
            }
            assert(position_ == size_);
            if (!init()) { // read the next chunk of compressed data
                break;
            }
            assert(position_ == 0U);
            block_begin = 0;
            block_end = std::min(size_, sizeof(buffer_) / sizeof(T) - buffer_end_);
        } else {
            assert(position_ < size_);
            assert(num_lower_bits_ < 8 * sizeof(T));
            const size_t pos_bits = position_ * num_lower_bits_;
            if (pos_bits - cur_pos_bits_ >= 8 * sizeof(T)) {
                cur_pos_bits_ += 8 * sizeof(T);
                lower_idx_++;
                if (lower_idx_ == READ_BUF_SIZE - 1 && num_lower_bytes_ > 0) {
                    const uint32_t to_read = std::min(sizeof(lower_) - sizeof(T), num_lower_bytes_);
                    lower_[0] = lower_[lower_idx_];
                    source_.read(reinterpret_cast<char *>(&lower_[1]), to_read);
                    num_lower_bytes_ -= to_read;
                    lower_idx_ = 0;
                }
            }
            const size_t adjusted_pos = pos_bits - cur_pos_bits_;
            const uint8_t *ptr
                    = reinterpret_cast<uint8_t *>(&lower_[lower_idx_]) + (adjusted_pos / 8);
            buffer_[buffer_end_++] = (load_unaligned<T>(ptr) >> (adjusted_pos % 8))
                                        & lower_bits_mask_;
            position_++;
        }
    }

    return buffer_end_;
}

// TODO: make this public and avoid reconstruction
template <typename T>
bool EliasFanoDecoder<T>::init() {
    position_ = 0;
    cur_pos_bits_ = 0;
    lower_idx_ = 0;
    memset(lower_, 0, sizeof(lower_));
    upper_pos_ = 0;
    if (!source_.read(reinterpret_cast<char *>(&size_), sizeof(size_t))) {
        if (remove_source_) {
            std::filesystem::remove(source_name_);
            std::filesystem::remove(source_name_ + ".up");
        }
        size_ = static_cast<size_t>(-1);
        return false;
    }
    source_.read(reinterpret_cast<char *>(&offset_), sizeof(T));
    source_.read(reinterpret_cast<char *>(&num_lower_bits_), 1);
    assert(num_lower_bits_ < 8 * sizeof(T));
    lower_bits_mask_ = (T(1) << num_lower_bits_) - 1UL;
    source_.read(reinterpret_cast<char *>(&num_lower_bytes_), sizeof(size_t));
    source_.read(reinterpret_cast<char *>(&num_upper_bytes_), sizeof(size_t));
    size_t low_bytes_read = std::min(sizeof(lower_), num_lower_bytes_);
    source_.read(reinterpret_cast<char *>(lower_), low_bytes_read);
    num_lower_bytes_ -= low_bytes_read;

    // Reserve a bit extra space for unaligned reads and set all to
    // zero to silence Valgrind uninitilized memory warnings
    upper_.resize((num_upper_bytes_ + 7) / 8, 0);
    source_upper_.read(reinterpret_cast<char *>(upper_.data()), num_upper_bytes_);
    assert(static_cast<uint32_t>(source_upper_.gcount()) == num_upper_bytes_);

    return true;
}

// -------------------------- EliasFanoEncoder<std::pair> -------------------------------
template <typename T, typename C>
EliasFanoEncoder<std::pair<T, C>>::EliasFanoEncoder(size_t size,
                                                    const T &first_value,
                                                    const T &last_value,
                                                    const std::string &sink_name)
    : ef_encoder(size, first_value, last_value, sink_name),
      sink_second_name_(sink_name + ".count"),
      sink_second_(sink_second_name_, std::ios::binary) {}

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
EliasFanoDecoder<std::pair<T, C>>::EliasFanoDecoder(const std::string &source,
                                                    bool remove_source)
    : source_first_(source, remove_source),
      source_second_name_(source + ".count"),
      remove_source_(remove_source) {
    source_second_ = std::ifstream(source_second_name_, std::ios::binary);
    if (!source_second_.good()) {
        logger->error("Unable to read from {}", source_second_name_);
        std::exit(EXIT_FAILURE);
    }
}

// ------------------------------ EliasFanoEncoderBuffered ----------------------------
template <typename T>
EliasFanoEncoderBuffered<T>::EliasFanoEncoderBuffered(const std::string &file_name,
                                                      size_t buffer_size)
    : file_name_(file_name) {
    sink_ = std::ofstream(file_name, std::ios::binary);
    sink_upper_ = std::ofstream(file_name + ".up", std::ios::binary);
    buffer_.reserve(buffer_size);
}

template <typename T>
EliasFanoEncoderBuffered<T>::~EliasFanoEncoderBuffered() {
    assert(!sink_.is_open());
}


template <typename T>
size_t EliasFanoEncoderBuffered<T>::finish() {
    encode_chunk();
    sink_.close();
    sink_upper_.close();
    return total_size_;
}

template <typename T>
void EliasFanoEncoderBuffered<T>::encode_chunk() {
    EliasFanoEncoder<T> encoder_(buffer_, &sink_, &sink_upper_);
    total_size_ += encoder_.finish();
    buffer_.resize(0);
}

// ----------------------- EliasFanoEncoderBuffered<std::pair> --------------------------
template <typename T, typename C>
EliasFanoEncoderBuffered<std::pair<T, C>>::EliasFanoEncoderBuffered(const std::string &file_name,
                                                                    size_t buffer_size)
    : sink_(file_name, std::ios::binary),
      sink_upper_(file_name + ".up", std::ios::binary),
      sink_second_(file_name + ".count", std::ios::binary),
      file_name_(file_name) {
    buffer_.reserve(buffer_size);
    buffer_second_.reserve(buffer_size);
}

template <typename T, typename C>
EliasFanoEncoderBuffered<std::pair<T, C>>::~EliasFanoEncoderBuffered() {
    assert(!sink_.is_open());
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
    EliasFanoEncoder<T> encoder_(buffer_, &sink_, &sink_upper_);

    sink_second_.write(reinterpret_cast<char *>(buffer_second_.data()),
                       buffer_second_.size() * sizeof(C));
    total_size_ += (encoder_.finish() + buffer_second_.size() * sizeof(C));
    buffer_.resize(0);
    buffer_second_.resize(0);
}

// instantiate used templates
template class EliasFanoEncoder<uint64_t>;
template class EliasFanoEncoder<sdsl::uint128_t>;
template class EliasFanoEncoder<sdsl::uint256_t>;
template class EliasFanoEncoder<std::pair<uint64_t, uint8_t>>;
template class EliasFanoEncoder<std::pair<uint64_t, uint16_t>>;
template class EliasFanoEncoder<std::pair<uint64_t, uint32_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint128_t, uint8_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint128_t, uint16_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint128_t, uint32_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint256_t, uint8_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint256_t, uint16_t>>;
template class EliasFanoEncoder<std::pair<sdsl::uint256_t, uint32_t>>;

template class EliasFanoDecoder<uint64_t>;
template class EliasFanoDecoder<sdsl::uint128_t>;
template class EliasFanoDecoder<sdsl::uint256_t>;
template class EliasFanoDecoder<std::pair<uint64_t, uint8_t>>;
template class EliasFanoDecoder<std::pair<uint64_t, uint16_t>>;
template class EliasFanoDecoder<std::pair<uint64_t, uint32_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint128_t, uint8_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint128_t, uint16_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint128_t, uint32_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint256_t, uint8_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint256_t, uint16_t>>;
template class EliasFanoDecoder<std::pair<sdsl::uint256_t, uint32_t>>;

template class EliasFanoEncoderBuffered<uint64_t>;
template class EliasFanoEncoderBuffered<sdsl::uint128_t>;
template class EliasFanoEncoderBuffered<sdsl::uint256_t>;
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
} // namespace mtg
