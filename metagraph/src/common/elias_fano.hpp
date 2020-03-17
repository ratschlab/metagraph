#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <optional>

#include <sdsl/bits.hpp>
#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#include "common/vector.hpp"


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

/**
 * Elias-Fano encoder that streams the encoded result into a file.
 * Loosely inspired  by
 * https://github.com/facebook/folly/blob/master/folly/experimental/EliasFanoCoding.h
 */
template <typename T>
class EliasFanoEncoder {
  public:
    EliasFanoEncoder() {}

    /**
     * Constructs an Elias-Fano encoder of an array with the given size and given max
     * value. The encoded output is written to #sink.
     */
    EliasFanoEncoder(size_t size, T max_value, const std::string &sink_name)
        : declared_size_(size) {
        // open file for appending, as we may encode multiple compressed chunks in the same file
        init(size, sink_name, max_value);
    }

    /**
     * Adds a new value to be encoded.
     */
    void add(T value) {
        assert(value >= last_value_);

        const T upper_bits = value >> num_lower_bits_;

        // We are adding the size_-th element, so we have a 1 followed by upper_bits
        // zeros, plus the 1s for the previous size_ elements; this is not trivial to
        // understand, so spend some time thinking about why this is correct
        const T pos = upper_bits + size_;
        upper_[pos / 8] |= 1U << (pos % 8);

        // Append the #num_lower_bits_ bits of #value to #lower_
        if (num_lower_bits_ != 0) {
            const T lowerBits = value & ((T(1) << num_lower_bits_) - 1);
            size_t pos_bits = size_ * num_lower_bits_;
            if (pos_bits - cur_pos_lbits_ >= 64) { // first 64 bits are ready to be written
                cur_pos_lbits_ += 64;
                sink_.write(reinterpret_cast<char *>(lower_.data()), sizeof(uint64_t));
                lower_[0] = lower_[1];
                lower_[1] = 0;
            }
            write_bits(reinterpret_cast<uint8_t *>(lower_.data()), pos_bits % 64,
                       num_lower_bits_, lowerBits);
        }

        last_value_ = value;
        ++size_;
    }

    size_t finish() {
        assert(size_ == declared_size_);
        // Append the remaining lower bits
        if (num_lower_bits_ != 0) {
            size_t cur_pos_bytes = (cur_pos_lbits_ + 7) / 8;
            assert(cur_pos_bytes <= num_lower_bytes_);
            assert(num_lower_bytes_ - cur_pos_bytes < 16);
            sink_.write(reinterpret_cast<char *>(lower_.data()),
                        num_lower_bytes_ - cur_pos_bytes);
        }
        if (size_ > 0) {
            sink_.write(upper_.data(), num_upper_bytes_);
        }
        sink_.close();
        return num_lower_bytes_ + num_upper_bytes_ + sizeof(size_) + sizeof(num_lower_bits_)
                + sizeof(num_upper_bytes_) + sizeof(num_lower_bytes_);
    }

  private:
    void init(size_t size, const std::string &sink_name, T max_value) {
        sink_ = std::ofstream(sink_name, std::ios::binary | std::ios::app);
        if (!sink_.good()) {
            std::cerr << "Unable to write to " << sink_name << std::endl;
            std::exit(EXIT_FAILURE);
        }
        // cap at 56 because #write_bits supports a max of 56 bits
        num_lower_bits_
                = std::min(get_num_lower_bits(max_value, size), static_cast<uint8_t>(56));

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

        sink_.write(reinterpret_cast<char *>(&size), sizeof(size_t));
        sink_.write(reinterpret_cast<char *>(&num_lower_bits_), 1);
        sink_.write(reinterpret_cast<char *>(&num_lower_bytes_), sizeof(size_t));
        sink_.write(reinterpret_cast<char *>(&num_upper_bytes_), sizeof(size_t));
    }

    /**
     * Returns the number of lower bits used in the Elias-Fano encoding of a sorted array
     * of size #size and maximum value max_value.
     */
    static uint8_t get_num_lower_bits(T max_value, size_t size) {
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

        // NOTE: As size != 0, "candidate" is always < 64.
        return (size > static_cast<uint64_t>(max_value >> candidate)) ? candidate - 1
                                                                      : candidate;
    }

    /** Writes #value (with len up to 56 bits) to #data starting at the #pos-th bit. */
    static void write_bits(uint8_t *data,
                           size_t pos,
                           uint8_t __attribute__((unused)) len,
                           uint64_t value) {
        assert(uint32_t(len) < 56);
        assert(0 == (value & ~((uint64_t(1) << len) - 1)));
        unsigned char *const ptr = data + (pos / 8);
        uint64_t ptrv = load_unaligned<uint64_t>(ptr);
        ptrv |= value << (pos % 8);
        store_unaligned<uint64_t>(ptr, ptrv);
    }

  private:
    /**
     * The lower bits of the encoded number, obtained by simply concatenating the
     * binary representation of the lower bits of each number.
     * To save memory, only the last 16 bytes are kept in memory. As soon as a chunk of 8
     * bytes is ready to be written, we flush it to #sink_ and shift the data in
     * #lower_ to the left by 8 bytes.
     */
    std::array<uint64_t, 2> lower_ = { 0, 0 };

    /**
     * Upper bits of the encoded numbers. Upper bits are stored using unary delta
     * encoding, with a 1 followed by as many zeros as the value to encode. For example,
     * the  upper bits, (3 5 5 9) will be encoded as the deltas (3 2 0 4). The 3 is
     * encoded as 1000, the 2 as 100, the 0 as 1 and the 4 as 10000,  resulting in
     * 1000011001000 in base 2.
     */
    Vector<char> upper_;

    /**
     * Current number of elements added for encoding.
     */
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

    /**
     * The size in bytes of lower_, without the 7 byte padding.
     */
    size_t num_lower_bytes_;
    /**
     * The size in bytes of upper_, without the 7 byte padding.
     */
    size_t num_upper_bytes_;
    /**
     * The last value that was added to the encoder. Only used to assert that the
     * numbers are added in increasing order.
     */
    T last_value_ = T(0);

    /**
     * Sink to write the encoded values to.
     */
    std::ofstream sink_;

    /**
     * Number of lower bits that were written to disk.
     */
    size_t cur_pos_lbits_ = 0;
};

/**
 * Decodes a list of compressed sorted integers stored in a file using #EliasFanoEncoder.
 */
template <typename T>
class EliasFanoDecoder {
  public:
    /** Creates a decoder that retrieves data from the given source */
    EliasFanoDecoder(const std::string &source_name) {
        if (std::filesystem::file_size(source_name) == 0) {
            position_ = size_ = 0;
            all_read_ = true;
            return;
        }
        source_ = std::ifstream(source_name, std::ios::binary);
        // profiling indicates that setting a larger buffer slightly increases performance
        source_.rdbuf()->pubsetbuf(buffer_, 1024 * 1024);
        if (!source_.good()) {
            std::cerr << "Unable to read from " << source_name << std::endl;
            std::exit(EXIT_FAILURE);
        }
        source_.seekg(0, source_.end);
        file_end_ = source_.tellg();
        source_.seekg(0, source_.beg);
        init();
    }

    /**
     * Returns the upper part of the next compressed element.
     */
    T next_upper() {
        // Skip to the first non-zero block.
        while (upper_block_ == 0) {
            upper_pos_ += sizeof(uint64_t);
            upper_block_ = load_unaligned<uint64_t>(upper_.data() + upper_pos_);
        }

        size_t trailing_zeros = __builtin_ctzll(upper_block_); // count trailing zeros
        upper_block_ = upper_block_ & (upper_block_ - 1); // reset the lowest 1 bit

        return static_cast<T>(8 * upper_pos_ + trailing_zeros - position_);
    }

    /**
     * Returns the lower part of the next compressed element.
     */
    T next_lower() {
        assert(position_ < size_);
        const size_t pos_bits = position_ * num_lower_bits_;
        const size_t pos_bytes = (pos_bits + 7) / 8 + 8;
        if (pos_bits - cur_pos_bits_ >= 64) {
            cur_pos_bits_ += 64;
            lower_[0] = lower_[1];
            lower_[1] = 0;
            if (num_lower_bytes_ > pos_bytes) {
                source_.read(reinterpret_cast<char *>(&lower_[1]),
                             std::min(sizeof(uint64_t), num_lower_bytes_ - pos_bytes));
            }
        }
        const size_t adjusted_pos = pos_bits - cur_pos_bits_;

        const uint8_t *ptr = reinterpret_cast<uint8_t *>(lower_.data()) + (adjusted_pos / 8);
        const uint64_t ptrv = load_unaligned<uint64_t>(ptr);
        return T(clear_high_bits(ptrv >> (adjusted_pos % 8), num_lower_bits_));
    }

    /**
     * Returns the next compressed element or empty if all elements were read.
     */
    std::optional<T> next() {
        if (position_ == size_) {
            if (all_read_) {
                return {};
            }
            init(); // read the next chunk of compressed data
        }
        T result = next_lower() | (next_upper() << num_lower_bits_);
        position_++;
        return result;
    }

  private:
    void init() {
        position_ = 0;
        cur_pos_bits_ = 0;
        upper_block_ = 0;
        lower_ = { 0, 0 };
        // Initialized to a negative number to save on decrement instruction in
        // #next_upper.
        upper_pos_ = static_cast<size_t>(-sizeof(size_t));
        source_.seekg(num_upper_bytes_, std::ios::cur); // move to beg of chunk
        source_.read(reinterpret_cast<char *>(&size_), sizeof(size_t));
        source_.read(reinterpret_cast<char *>(&num_lower_bits_), 1);
        source_.read(reinterpret_cast<char *>(&num_lower_bytes_), sizeof(size_t));
        source_.read(reinterpret_cast<char *>(&num_upper_bytes_), sizeof(size_t));
        size_t low_bytes_read = std::min(2 * sizeof(uint64_t), num_lower_bytes_);
        source_.read(reinterpret_cast<char *>(lower_.data()), low_bytes_read);

        std::streampos pos = source_.tellg();
        // to avoid jumping through the file, we read the relatively small upper_bytes
        // into memory;
        source_.seekg(num_lower_bytes_ - low_bytes_read, std::ios::cur);
        upper_.reserve(num_upper_bytes_ + 7);
        upper_.resize(num_upper_bytes_);
        source_.read(upper_.data(), num_upper_bytes_);
        assert(source_.gcount() == num_upper_bytes_);
        if (source_.tellg() == file_end_) {
            all_read_ = true;
        }
        source_.seekg(pos, source_.beg);
    }

    /** Clear the high bits of #value starting at position #index. */
    static uint64_t clear_high_bits(uint64_t value, uint32_t index) {
        assert(index < 64);
        return value & ((uint64_t(1) << index) - 1);
    }

  private:
    /** Index of current element */
    size_t position_;
    /**
     * Current position in the #upper_ vector.
     */
    size_t upper_pos_ = static_cast<size_t>(-sizeof(size_t));

    /** The sequence of 8 upper bytes currently being processed */
    uint64_t upper_block_;

    /** Number of lower bits that were read from disk. */
    size_t cur_pos_bits_;

    /**
     * The lower bits of the encoded number, obtained by simply concatenating the
     * binary representation of the lower bits of each number. To save memory, only the
     * currently needed window of 16 bytes is read from the file.
     */
    std::array<uint64_t, 2> lower_;

    /**
     * Upper bits of the encoded numbers. Upper bits are stored using unary delta
     * encoding, with a 1 followed by as many zeros as the value to encode.
     */
    Vector<char> upper_;

    /** Total number of elements encoded. */
    size_t size_;

    /**
     * Each encoded integer is split into a "lower" and an "upper" part. This is the
     * number of bits used for the "lower" part of the Elias-Fano encoding.
     */
    uint8_t num_lower_bits_;

    /** The size in bytes of lower_. */
    size_t num_lower_bytes_;

    /** The size in bytes of upper_. */
    size_t num_upper_bytes_ = 0;

    /** Stream containing the compressed data. */
    std::ifstream source_;

    /* 1MB buffer for reading from #source_ */
    char buffer_[1024 * 1024];

    std::streampos file_end_;

    /** True if we read all bytes from source_ */
    bool all_read_ = false;
};

/**
 * Encoder specialization for an std::pair. The first member of the pair is assumed to be
 * in nondecreasing order and is compressed using Elias-Fano encoding, while the second
 * member of the pair is written to disk as is.
 */
template <typename T, typename C>
class EliasFanoEncoder<std::pair<T, C>> {
  public:
    EliasFanoEncoder() {}
    /**
     * Constructs an Elias-Fano encoder of an array with the given size and given last
     * value. The encoded output is dumped to #sink_name.
     */
    EliasFanoEncoder(size_t size, const std::pair<T, C> &last_value, const std::string &sink_name)
        : ef_encoder(size, last_value.first, sink_name),
          sink_second_name_(sink_name + ".count") {
        // open file for appending, as we may encode multiple compressed chunks in the same file
        sink_second_ = std::ofstream(sink_second_name_, std::ios::binary | std::ios::app);
    }

    /**
     * Adds a new value to be encoded.
     */
    void add(std::pair<T, C> value) {
        ef_encoder.add(value.first);
        sink_second_.write(reinterpret_cast<char *>(&value.second), sizeof(C));
    }

    size_t finish() {
        size_t first_size = ef_encoder.finish();
        sink_second_.close();
        return first_size + std::filesystem::file_size(sink_second_name_);
    }

  private:
    EliasFanoEncoder<T> ef_encoder;
    std::string sink_second_name_;
    std::ofstream sink_second_;
};

template <typename T, typename C>
class EliasFanoDecoder<std::pair<T, C>> {
  public:
    /** Creates a decoder that retrieves data from the given source */
    EliasFanoDecoder(const std::string &source) : source_first_(source) {
        std::string source_second_name = source + ".count";
        source_second_ = std::ifstream(source_second_name, std::ios::binary);
        // profiling indicates that setting a larger buffer slightly increases performance
        source_second_.rdbuf()->pubsetbuf(buffer_, 1024 * 1024);
        if (!source_second_.good()) {
            std::cerr << "Unable to write to " << source_second_name << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    /**
     * Returns the next compressed element or empty if all elements were read.
     */
    std::optional<std::pair<T, C>> next() {
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

  private:
    EliasFanoDecoder<T> source_first_;
    std::ifstream source_second_;
    /* 1MB buffer for reading from #source_second_ */
    char buffer_[1024 * 1024];
};

/**
 * Template specialization for 128 bit integers that simply writes the integers to file uncompressed.
 * TODO(dd): separate the 128 bit integer into 2 64 bit ones and compress in chunks,
 * taking care that in each cunk the lower 64 bits are in increasing order.
 */
template <>
class EliasFanoEncoder<sdsl::uint128_t> {
  public:
    EliasFanoEncoder() : EliasFanoEncoder<sdsl::uint128_t>(0, sdsl::uint128_t(0), "") {}

    /**
     * Constructs an Elias-Fano encoder of an array with the given size and given max
     * value. The encoded output is written to #sink.
     */
    EliasFanoEncoder(size_t size,
                     __attribute__((unused)) sdsl::uint128_t max_value,
                     const std::string &sink_name)
        : declared_size_(size) {
        // open file for appending, as we may encode multiple compressed chunks in the same file
        sink_ = std::ofstream(sink_name, std::ios::binary | std::ios::app);
        if (!sink_.good()) {
            std::cerr << "Unable to write to " << sink_name << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    void add(const sdsl::uint128_t &value) {
        sink_.write(reinterpret_cast<const char *>(&value), sizeof(sdsl::uint128_t));
        total_size_ += sizeof(sdsl::uint128_t);
        size_++;
    }

    size_t finish() {
        assert(size_ == declared_size_);
        sink_.close();
        return total_size_;
    }

  private:
    size_t declared_size_;

    /**
     * Sink to write the encoded values to.
     */
    std::ofstream sink_;

    size_t total_size_ = 0;

    size_t size_ = 0;
};

/**
 * Template specialization for 128 bit integers that simply writes the integers to file uncompressed.
 */
template <>
class EliasFanoEncoder<sdsl::uint256_t> {
  public:
    EliasFanoEncoder() : EliasFanoEncoder<sdsl::uint256_t>(0, sdsl::uint256_t(0), "") {}
    /**
     * Constructs an Elias-Fano encoder of an array with the given size and given max
     * value. The encoded output is written to #sink.
     */
    EliasFanoEncoder(size_t size,
                     __attribute__((unused)) sdsl::uint256_t max_value,
                     const std::string &sink_name)
        : declared_size_(size) {
        // open file for appending, as we may encode multiple compressed chunks in the same file
        sink_ = std::ofstream(sink_name, std::ios::binary | std::ios::app);
        if (!sink_.good()) {
            std::cerr << "Unable to write to " << sink_name << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    void add(const sdsl::uint256_t &value) {
        sink_.write(reinterpret_cast<const char *>(&value), sizeof(sdsl::uint256_t));
        total_size_ += sizeof(sdsl::uint256_t);
        size_++;
    }

    size_t finish() {
        assert(size_ == declared_size_);
        sink_.close();
        return total_size_;
    }

  private:
    size_t declared_size_;

    /**
     * Sink to write the encoded values to.
     */
    std::ofstream sink_;

    size_t total_size_ = 0;

    size_t size_ = 0;
};

/**
 * Decodes a list of compressed sorted integers stored in a file using #EliasFanoEncoder.
 */
template <>
class EliasFanoDecoder<sdsl::uint128_t> {
  public:
    /** Creates a decoder that retrieves data from the given source */
    EliasFanoDecoder(const std::string &source_name)  {
        source_ = std::ifstream(source_name, std::ios::binary);
        if (!source_.good()) {
            std::cerr << "Unable to read from " << source_name << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    /**
     * Returns the next compressed element or empty if all elements were read.
     */
    std::optional<sdsl::uint128_t> next() {
        sdsl::uint128_t result;
        if (source_.read(reinterpret_cast<char *>(&result), sizeof(sdsl::uint128_t))) {
            return result;
        }
        return {};
    }


  private:
    /** Stream containing the compressed data. */
    std::ifstream source_;
};

/**
 * Decodes a list of compressed sorted integers stored in a file using #EliasFanoEncoder.
 */
template <>
class EliasFanoDecoder<sdsl::uint256_t> {
  public:
    /** Creates a decoder that retrieves data from the given source */
    EliasFanoDecoder(const std::string &source_name)  {
        source_ = std::ifstream(source_name, std::ios::binary);
        if (!source_.good()) {
            std::cerr << "Unable to read from " << source_name << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    /**
     * Returns the next compressed element or empty if all elements were read.
     */
    std::optional<sdsl::uint256_t> next() {
        sdsl::uint256_t result;
        if (source_.read(reinterpret_cast<char *>(&result), sizeof(sdsl::uint256_t))) {
            return result;
        }
        return {};
    }


  private:
    /** Stream containing the compressed data. */
    std::ifstream source_;
};

/**
 * Specialization of #EliasFanoEncoder that can encode sequences of unknown size. It uses
 * a buffer to acummulate data and then dumps it in chunks to an EliasFanoEncoder.
 */
template <typename T>
class EliasFanoEncoderBuffered {
  public:
    EliasFanoEncoderBuffered(const std::string &file_name, size_t buffer_size)
        : file_name_(file_name) {
        std::ofstream ofs(file_name, std::ofstream::out | std::ofstream::trunc);
        ofs.close();
        buffer_.reserve(buffer_size);
    }

    void add(const T &value) {
        buffer_.push_back(value);
        if (buffer_.size() == buffer_.capacity()) {
            encode_chunk();
        }
    }

    size_t finish() {
        encode_chunk();
        return total_size_;
    }

  private:
    void encode_chunk() {
        if (buffer_.empty()) {
            return;
        }
        encoder_ = EliasFanoEncoder<T>(buffer_.size(), buffer_.back(), file_name_);
        for (const auto &v : buffer_) {
            encoder_.add(v);
        }
        total_size_ += encoder_.finish();
        buffer_.resize(0);
    }

  private:
    EliasFanoEncoder<T> encoder_;
    Vector<T> buffer_;
    std::string file_name_;
    size_t total_size_ = 0;
};

/**
 * Specialization of #EliasFanoEncoder that can encode sequences of pairs of unknown size.
 * It uses a buffer to acummulate data and then dumps it in chunks to an EliasFanoEncoder.
 */
template <typename T, typename C>
class EliasFanoEncoderBuffered<std::pair<T, C>> {
  public:
    EliasFanoEncoderBuffered(const std::string &file_name, size_t buffer_size)
        : file_name_(file_name) {
        std::ofstream ofs(file_name, std::ofstream::out | std::ofstream::trunc);
        ofs.close();
        buffer_.reserve(buffer_size);
    }

    void add(const std::pair<T, C> &value) {
        buffer_.push_back(value);
        if (buffer_.size() == buffer_.capacity()) {
            encode_chunk();
        }
    }

    size_t finish() {
        encode_chunk();
        return total_size_;
    }

  private:
    void encode_chunk() {
        if (buffer_.size() == 0) {
            return;
        }
        encoder_ = EliasFanoEncoder<std::pair<T, C>>(buffer_.size(), buffer_.back(),
                                                     file_name_);
        for (const auto &v : buffer_) {
            encoder_.add(v);
        }
        total_size_ += encoder_.finish();
        buffer_.resize(0);
    }

  private:
    EliasFanoEncoder<std::pair<T, C>> encoder_;
    Vector<std::pair<T, C>> buffer_;
    std::string file_name_;
    size_t total_size_ = 0;
};

} // namespace common
} // namespace mg
