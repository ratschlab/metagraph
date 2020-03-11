#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <optional>

#include <folly/lang/Bits.h>

#include "common/vector.hpp"

namespace mg {
namespace common {

/** Writes #value (with len up to 56 bits) to #data starting at the #pos-th bit. */
void write_bits(uint8_t *data, size_t pos, uint8_t len, uint64_t value);

/**
 * Returns the number of lower bits used in the Elias-Fano encoding of a sorted array of
 * size #size and maximum value max_value.
 */
uint8_t get_num_lower_bits(size_t max_value, size_t size);

/** Clear the high bits of #value starting at position #index. */
uint64_t clear_high_bits(uint64_t value, uint32_t index);

/**
 * Elias-Fano encoder that streams the encoded result into a file.
 * Loosely inspired  by
 * https://github.com/facebook/folly/blob/master/folly/experimental/EliasFanoCoding.h
 */
template <typename T>
class EliasFanoEncoder {
  public:
    /**
     * Constructs an Elias-Fano encoder of an array with the given size and given max
     * value. The encoded output is written to #sink.
     */
    EliasFanoEncoder(size_t size, T max_value, std::ofstream &sink)
        : declared_size_(size), sink_(sink) {
        init(size, max_value);
    }

    /**
     * Encodes the given vector using Elias-Fano encoding. The encoded output is
     * written to #sink. The vector elements must be non-decreasing.
     */
    EliasFanoEncoder(const Vector<T> &data, std::ofstream &sink)
        : EliasFanoEncoder(data.size(), data.back(), sink) {
        for (const auto &v : data) {
            add(v);
        }
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
        const size_t pos = upper_bits + size_;
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
    void init(size_t size, T max_value) {
        // cap at 56 because #write_bits supports a max of 56 bits
        num_lower_bits_
                = std::min(get_num_lower_bits(max_value, size), static_cast<uint8_t>(56));

        // Number of 0-bits to be stored + 1-bits
        const size_t upper_size_bits = (max_value >> num_lower_bits_) + size;
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
    size_t declared_size_;

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
    std::ofstream &sink_;

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
    EliasFanoDecoder(std::ifstream &source) : source_(source) {
        source.read(reinterpret_cast<char *>(&size_), sizeof(size_t));
        source.read(reinterpret_cast<char *>(&num_lower_bits_), 1);
        source.read(reinterpret_cast<char *>(&num_lower_bytes_), sizeof(size_t));
        source.read(reinterpret_cast<char *>(&num_upper_bytes_), sizeof(size_t));
        source.read(reinterpret_cast<char *>(lower_.data()),
                    std::min(2 * sizeof(uint64_t), num_lower_bytes_));
        std::streampos pos = source.tellg();
        // to avoid jumping through the file, we read the relatively small upper_bytes
        // into memory; TODO(dd): consider writing upperbytes to a separate file to avoid
        // loading it into memory
        source.seekg(-num_upper_bytes_, source.end);
        upper_.reserve(num_upper_bytes_ + 7);
        upper_.resize(num_upper_bytes_);
        source.read(upper_.data(), num_upper_bytes_);
        source.seekg(pos, source.beg);
    }

    /**
     * Returns the upper part of the next compressed element.
     */
    T next_upper() {
        // Skip to the first non-zero block.
        while (upper_block_ == 0) {
            upper_pos_ += sizeof(uint64_t);
            upper_block_ = folly::loadUnaligned<uint64_t>(upper_.data() + upper_pos_);
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
        const uint64_t ptrv = folly::loadUnaligned<uint64_t>(ptr);
        return clear_high_bits(ptrv >> (adjusted_pos % 8), num_lower_bits_);
    }

    /**
     * Returns the next compressed element or empty if all elements were read.
     */
    std::optional<T> next() {
        if (position_ == size_) {
            return {};
        }
        T result = next_lower() | (next_upper() << num_lower_bits_);
        position_++;
        return result;
    }

  private:
    /** Index of current element */
    size_t position_ = 0;
    /**
     * Current position in the #upper_ vector.
     * Initialized to a negative number to avoid one decrement instruction in #next_upper.
     */
    size_t upper_pos_ = static_cast<size_t>(-sizeof(size_t));

    /** The sequence of 8 upper bytes currently being processed */
    uint64_t upper_block_ = 0;

    /** Number of lower bits that were read from disk. */
    size_t cur_pos_bits_ = 0;

    /**
     * The lower bits of the encoded number, obtained by simply concatenating the
     * binary representation of the lower bits of each number. To save memory, only the
     * currently needed window of 16 bytes is read from the file.
     */
    std::array<uint64_t, 2> lower_ = { 0, 0 };

    /**
     * Upper bits of the encoded numbers. Upper bits are stored using unary delta
     * encoding, with a 1 followed by as many zeros as the value to encode.
     */
    Vector<char> upper_;

    /** Total number of elements encoded. */
    size_t size_ = 0;

    /**
     * Each encoded integer is split into a "lower" and an "upper" part. This is the
     * number of bits used for the "lower" part of the Elias-Fano encoding.
     */
    uint8_t num_lower_bits_;

    /** The size in bytes of lower_. */
    size_t num_lower_bytes_;

    /** The size in bytes of upper_. */
    size_t num_upper_bytes_;

    /** Stream containing the compressed data. */
    std::ifstream &source_;
};

} // namespace common
} // namespace mg
