#include "elias_fano.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <filesystem>

#include <sdsl/uint128_t.hpp>

#include "common/logger.hpp"
#include "common/utils/template_utils.hpp"


namespace mtg {
namespace elias_fano {

using mtg::common::logger;

void concat(const std::vector<std::string> &files, const std::string &result) {
    if (files.empty())
        return;

    std::vector<std::string> suffixes = { "", ".up" };
    if (std::filesystem::exists(files[0] + ".count"))
        suffixes.push_back(".count");

    for (const std::string &suffix : suffixes) {
        if (files.size() > 1) {
            std::string concat_command = "cat ";
            for (uint32_t i = 1; i < files.size(); ++i) {
                concat_command += files[i] + suffix + " ";
            }
            concat_command += ">> " + files[0] + suffix;

            if (std::system(concat_command.c_str())) {
                logger->error("Error while cat-ing files: {}", concat_command);
                std::exit(EXIT_FAILURE);
            }
        }

        std::filesystem::rename(files[0] + suffix, result + suffix);
        for (size_t i = 1; i < files.size(); ++i) {
            std::string fname = files[i] + suffix;
            if (!std::filesystem::exists(fname))
                logger->error("Trying to remove non-existent file: {}", fname);

            std::filesystem::remove(fname);
        }
    }
}

void remove_chunks(const std::vector<std::string> &files) {
    for (const std::string &f : files) {
        std::vector<std::string> suffixes = { "", ".up" };
        if (std::filesystem::exists(f + ".count"))
            suffixes.push_back(".count");

        for (const std::string &suffix : suffixes) {
            std::string fname = f + suffix;
            if (!std::filesystem::exists(fname))
                logger->error("Trying to remove non-existent file: {}", fname);

            std::filesystem::remove(f + suffix);
        }
    }
}

uint64_t chunk_size(const std::string &file) {
    std::vector<std::string> suffixes = { "", ".up" };
    if (std::filesystem::exists(file + ".count"))
        suffixes.push_back(".count");

    uint64_t size = 0;
    for (const std::string &suffix : suffixes) {
        size += std::filesystem::file_size(file + suffix);
    }
    return size;
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

void close_and_check(std::ofstream &sink, const std::string &fname) {
    sink.close();
    if (sink.fail()) {
        logger->error("Unable to close file {}", fname);
        std::exit(EXIT_FAILURE);
    }
}

// ------------ EliasFanoEncoder --------------------------------------------------------
/**
 * Elias-Fano encoder that streams the encoded result into a file.
 * Loosely inspired  by
 * https://github.com/facebook/folly/blob/master/folly/experimental/EliasFanoCoding.h
 */
template <typename T>
class EliasFanoEncoder {
  public:
    static constexpr uint32_t WRITE_BUF_SIZE = 1024;

    /** Encodes the #data array */
    static size_t append_block(const std::vector<T> &data,
                               std::ofstream *sink,
                               std::ofstream *sink_upper);
  private:
    /** Constructs an encoder that encodes the #data array */
    EliasFanoEncoder(const std::vector<T> &data, std::ofstream *sink, std::ofstream *sink_upper);

    /** Encodes the next number */
    void add(T value);

    /** Dumps any pending data to the stream. Must be called exactly once when done #add-ing */
    size_t finish();

    /**
     * Returns the number of lower bits used in the Elias-Fano encoding of a sorted array
     * of size #size and maximum value max_value.
     */
    static uint8_t get_num_lower_bits(T max_value, size_t size);

    /** Writes #value (with len up to 56 bits) to #data starting at the #pos-th bit. */
    static void write_bits(char *data, size_t pos, T value);

    void init(size_t size, T max_value);

    /**
     * The lower bits of the encoded number, obtained by simply concatenating the
     * binary representation of the lower bits of each number.
     * To save memory, only the last 2*sizeof(T) bytes are kept in memory. As soon as a
     * chunk of 8 bytes is ready to be written, we flush it to #sink_ and shift the data
     * in #lower_ to the left by sizeof(T) bytes.
     */
    char lower_[WRITE_BUF_SIZE * sizeof(T)];

    /**
     * Upper bits of the encoded numbers. Upper bits are stored using unary delta
     * encoding, with a 1 followed by as many zeros as the value to encode. For example,
     * the  upper bits, (3 5 5 9) will be encoded as the deltas (3 2 0 4). The 3 is
     * encoded as 1000, the 2 as 100, the 0 as 1 and the 4 as 10000,  resulting in
     * 1000011001000 in base 2.
     */
    std::vector<char> upper_;

    /** Current number of elements added for encoding */
    size_t size_ = 0;

    /**
     * Number of elements the decoder was initialized with. When all elements are added
     * the #declared_size_ must equal size_.
     */
    size_t declared_size_ = 0;

    /**
     * Each encoded integer is split into a "lower" and an "upper" part. This is the
     * number of bits used for the "lower" part of the Elias-Fano encoding. It is
     * capped at 56, as this is the maximum value supported by #write_bits
     */
    uint8_t num_lower_bits_;

    /** Mask to extract the lower bits from a value T. Equal to 2^#num_lower_bits_-1. */
    T lower_bits_mask_;

    /** The size in bytes of lower_, without the 7 byte padding */
    size_t num_lower_bytes_;
    /** The size in bytes of upper_, without the 7 byte padding */
    size_t num_upper_bytes_;
#ifndef NDEBUG
    /**
     * The last value that was added to the encoder. Only used to assert that the
     * numbers are added in increasing order.
     */
    T last_value_ = T(0);
#endif

    /**
     * Sink to write the encoded values to (except the upper bytes). Points to
     * an externally provided sink.
     * */
    std::ofstream *sink_;
    /**
     * Sink to write the upper bytes to. Points to an externally provided sink.
     * Upper bytes are written to a different sink in order
     * to avoid costly seek operations within the file.
     * */
    std::ofstream *sink_upper_;

    /** Number of lower bits that were written to disk */
    size_t cur_pos_lbits_ = 0;

    /** Offset to add to each element when decoding (used for minimizing the range). */
    T offset_ = 0;
};

template <typename T>
EliasFanoEncoder<T>::EliasFanoEncoder(const std::vector<T> &data,
                                      std::ofstream *sink,
                                      std::ofstream *sink_upper)
    : declared_size_(data.size()), sink_(sink), sink_upper_(sink_upper) {
    assert(sink_->is_open() && sink_upper_->is_open());
    if (data.size() == 0U) {
        return;
    }
    offset_ = data.front();
    init(data.size(), data.back());
}

template <typename T>
size_t EliasFanoEncoder<T>::append_block(const std::vector<T> &data,
                                         std::ofstream *sink,
                                         std::ofstream *sink_upper) {
    EliasFanoEncoder<T> encoder(data, sink, sink_upper);
    for (const auto &v : data) {
        encoder.add(v);
    }
    return encoder.finish();
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
    if (size_ == 0U)
        return 0;

    // Append the remaining lower bits
    if (num_lower_bits_ != 0) {
        size_t pos_bits = size_ * num_lower_bits_;
        size_t num_bytes = (pos_bits - cur_pos_lbits_ + 7) / 8;
        assert(cur_pos_lbits_ / 8 + num_bytes == num_lower_bytes_);
        sink_->write(lower_, num_bytes);
    }
    sink_upper_->write(upper_.data(), num_upper_bytes_);
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

    sink_->write(reinterpret_cast<const char *>(&size), sizeof(size_t));
    sink_->write(reinterpret_cast<const char *>(&offset_), sizeof(T));
    sink_->write(reinterpret_cast<const char *>(&num_lower_bits_), 1);
    sink_->write(reinterpret_cast<const char *>(&num_lower_bytes_), sizeof(size_t));
    sink_->write(reinterpret_cast<const char *>(&num_upper_bytes_), sizeof(size_t));
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
    source_ = std::ifstream(source_name, std::ios::binary);
    if (!source_) {
        logger->error("Unable to open {}", source_name);
        std::exit(EXIT_FAILURE);
    }
    source_upper_ = std::ifstream(source_name + ".up", std::ios::binary);
    if (!source_upper_) {
        logger->error("Unable to open {}", source_name + ".up");
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
    source_.read(reinterpret_cast<char *>(&size_), sizeof(size_t));
    if (source_.bad()) {
        logger->error("Error while reading from {}", source_name_);
        std::exit(EXIT_FAILURE);
    } else if (source_.eof()) {
        source_.close();
        source_upper_.close();
        if (remove_source_) {
            if (!std::filesystem::exists(source_name_))
                logger->error("Trying to remove non-existent file: {}", source_name_);

            if (!std::filesystem::exists(source_name_ + ".up"))
                logger->error("Trying to remove non-existent file: {}.up", source_name_);

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

    if (!source_.good() || !source_upper_.good()) {
        logger->error("Failed reading from {}", source_name_);
        std::exit(EXIT_FAILURE);
    }

    return true;
}

// ------------------------- EliasFandDecoder<std::pair> --------------------------------
template <typename T, typename C>
EliasFanoDecoder<std::pair<T, C>>::EliasFanoDecoder(const std::string &source,
                                                    bool remove_source)
    : source_first_(source, remove_source),
      source_second_name_(source + ".count"),
      remove_source_(remove_source) {
    source_second_ = std::ifstream(source_second_name_, std::ios::binary);
    if (!source_second_) {
        logger->error("Unable to open {}", source_second_name_);
        std::exit(EXIT_FAILURE);
    }
}

// ------------------------------ EliasFanoEncoderBuffered ----------------------------
template <typename T>
EliasFanoEncoderBuffered<T>::EliasFanoEncoderBuffered(const std::string &file_name,
                                                      size_t buffer_size,
                                                      bool append)
    : file_name_(file_name) {
    auto mode = append ? (std::ios::binary|std::ios::app) : std::ios::binary;
    sink_ = std::ofstream(file_name, mode);
    if (!sink_) {
        logger->error("Unable to open {} for writing", file_name);
        std::exit(EXIT_FAILURE);
    }
    sink_upper_ = std::ofstream(file_name + ".up", mode);
    if (!sink_upper_) {
        logger->error("Unable to open {} for writing", file_name + ".up");
        std::exit(EXIT_FAILURE);
    }
    buffer_.reserve(buffer_size);
}

template <typename T>
EliasFanoEncoderBuffered<T>::~EliasFanoEncoderBuffered() {
    assert(!sink_.is_open());
    assert(!sink_upper_.is_open());
}

/** Append sorted array #data to EF-coded #out_fname */
template <typename T>
size_t EliasFanoEncoderBuffered<T>::append_block(const std::vector<T> &data,
                                                 const std::string &file_name) {
    std::ofstream sink(file_name, std::ios::binary | std::ios::app);
    if (!sink) {
        logger->error("Unable to open {} for writing", file_name);
        std::exit(EXIT_FAILURE);
    }
    std::ofstream sink_upper(file_name + ".up", std::ios::binary | std::ios::app);
    if (!sink_upper) {
        logger->error("Unable to open {} for writing", file_name + ".up");
        std::exit(EXIT_FAILURE);
    }
    size_t total_size = EliasFanoEncoder<T>::append_block(data, &sink, &sink_upper);
    close_and_check(sink, file_name);
    close_and_check(sink_upper, file_name + ".up");
    return total_size;
}

template <typename T>
size_t EliasFanoEncoderBuffered<T>::finish() {
    encode_chunk();
    close_and_check(sink_, file_name_);
    close_and_check(sink_upper_, file_name_ + ".up");
    return total_size_;
}

template <typename T>
void EliasFanoEncoderBuffered<T>::encode_chunk() {
    total_size_ += EliasFanoEncoder<T>::append_block(buffer_, &sink_, &sink_upper_);
    buffer_.resize(0);
}

// ----------------------- EliasFanoEncoderBuffered<std::pair> --------------------------
template <typename T, typename C>
EliasFanoEncoderBuffered<std::pair<T, C>>::EliasFanoEncoderBuffered(const std::string &file_name,
                                                                    size_t buffer_size,
                                                                    bool append)
    : encoder_first_(file_name, buffer_size, append) {
    sink_second_ = std::ofstream(file_name + ".count",
                                 append ? (std::ios::binary|std::ios::app) : std::ios::binary);
    if (!sink_second_) {
        logger->error("Unable to open {} for writing", file_name + ".count");
        std::exit(EXIT_FAILURE);
    }
    buffer_second_.reserve(buffer_size);
}

template <typename T, typename C>
EliasFanoEncoderBuffered<std::pair<T, C>>::~EliasFanoEncoderBuffered() {
    assert(!sink_second_.is_open());
}

/** Append sorted array #data to EF-coded #out_fname */
template <typename T, typename C>
size_t EliasFanoEncoderBuffered<std::pair<T, C>>
::append_block(const std::vector<std::pair<T, C>> &data,
               const std::string &file_name) {
    std::vector<T> ts;
    std::vector<C> cs;
    ts.reserve(data.size());
    cs.reserve(data.size());
    for (const auto &[v, c] : data) {
        ts.push_back(v);
        cs.push_back(c);
    }
    size_t total_size = EliasFanoEncoderBuffered<T>::append_block(ts, file_name);

    std::ofstream sink_second(file_name + ".count", std::ios::binary|std::ios::app);
    if (!sink_second) {
        logger->error("Unable to open {} for writing", file_name + ".count");
        std::exit(EXIT_FAILURE);
    }
    sink_second.write(reinterpret_cast<const char *>(cs.data()), cs.size() * sizeof(C));

    return total_size + cs.size() * sizeof(C);
}

template <typename T, typename C>
size_t EliasFanoEncoderBuffered<std::pair<T, C>>::finish() {
    encode_chunk();
    close_and_check(sink_second_, name() + ".count");
    return encoder_first_.finish() + size() * sizeof(C);
}

template <typename T, typename C>
void EliasFanoEncoderBuffered<std::pair<T, C>>::encode_chunk() {
    sink_second_.write(reinterpret_cast<const char *>(buffer_second_.data()),
                       buffer_second_.size() * sizeof(C));
    buffer_second_.resize(0);
}

// instantiate used templates
template class EliasFanoDecoder<uint64_t>;
template class EliasFanoDecoder<sdsl::uint128_t>;
template class EliasFanoDecoder<sdsl::uint256_t>;
template class EliasFanoDecoder<std::pair<uint64_t, uint8_t>>;
template class EliasFanoDecoder<std::pair<uint64_t, uint16_t>>;
template class EliasFanoDecoder<std::pair<uint64_t, uint32_t>>;
template class EliasFanoDecoder<std::pair<uint64_t, uint64_t>>;
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
template class EliasFanoEncoderBuffered<std::pair<uint64_t, uint64_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint128_t, uint8_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint128_t, uint16_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint128_t, uint32_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint256_t, uint8_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint256_t, uint16_t>>;
template class EliasFanoEncoderBuffered<std::pair<sdsl::uint256_t, uint32_t>>;

} // namespace elias_fano
} // namespace mtg
