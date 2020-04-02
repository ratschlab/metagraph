#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <optional>

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>

#include "common/vector.hpp"


namespace mg {
namespace common {

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
     * Constructs an Elias-Fano encoder of an array with the given #size and given
     * #max_value. The encoded output is written to #out_filename.
     */
    EliasFanoEncoder(size_t size,
                     T min_value,
                     T max_value,
                     const std::string &out_filename,
                     bool is_append = false);

    /** Constructs an encoder that encodes the #data array */
    EliasFanoEncoder(const Vector<T> &data, std::ofstream *sink, std::ofstream *sink_upper);

    /** Encodes the next number */
    void add(T value);

    /** Dumps any pending data to the stream. Must be called exactly once when done #add-ing */
    size_t finish();

  private:
    /**
     * Returns the number of lower bits used in the Elias-Fano encoding of a sorted array
     * of size #size and maximum value max_value.
     */
    static uint8_t get_num_lower_bits(T max_value, size_t size);

    /** Writes #value (with len up to 56 bits) to #data starting at the #pos-th bit. */
    static void write_bits(uint8_t *data, size_t pos, T value);

    void init(size_t size, T max_value);

  private:
    /**
     * The lower bits of the encoded number, obtained by simply concatenating the
     * binary representation of the lower bits of each number.
     * To save memory, only the last 2*sizeof(T) bytes are kept in memory. As soon as a
     * chunk of 8 bytes is ready to be written, we flush it to #sink_ and shift the data
     * in #lower_ to the left by sizeof(T) bytes.
     */
    T lower_[2] = { 0, 0 };

    /**
     * Upper bits of the encoded numbers. Upper bits are stored using unary delta
     * encoding, with a 1 followed by as many zeros as the value to encode. For example,
     * the  upper bits, (3 5 5 9) will be encoded as the deltas (3 2 0 4). The 3 is
     * encoded as 1000, the 2 as 100, the 0 as 1 and the 4 as 10000,  resulting in
     * 1000011001000 in base 2.
     */
    Vector<char> upper_;

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
     * Sink to write the encoded values to (except the upper bytes). Points to either
     * #sink_internal or to an externally provided sink
     * */
    std::ofstream *sink_;
    /**
     * Sink to write the upper bytes to. Points to either #sink_internal_upper_ or to
     * an externally provided sink. Upper bytes are written to a different sink in order
     * to avoid costly seek operations within the file
     * */
    std::ofstream *sink_upper_;

    /** Internal sink for EF encoding except the upper bytes */
    std::ofstream sink_internal_;

    /**
     * Internal sink for the upper_ bytes, which are saved in a separate file to avoid
     * costly seekg/tellg opreations.
     */
    std::ofstream sink_internal_upper_;

    /** Number of lower bits that were written to disk */
    size_t cur_pos_lbits_ = 0;

    /** Offset to add to each element when decoding (used for minimizing the range). */
    T offset_ = 0;
};

/**
 * Decodes a list of compressed sorted integers stored in a file using #EliasFanoEncoder.
 */
template <typename T>
class EliasFanoDecoder {
    static constexpr uint32_t READ_BUF_SIZE = 1024;
  public:
    EliasFanoDecoder() {}

    /** Creates a decoder that retrieves data from the given file */
    EliasFanoDecoder(const std::string &source_name);

    /**
     * Create a decoder that reads the data from the provided streams.
     * @param source the stream containing all encoded data, except the upper bytes
     * @param source_upper the stream containing the upper bytes
     * @param file_end_upper the position where the algorithm should stop reading. Note
     * that this is not necessarily the end of the file, see #EliasFanoEncoderBuffered.
     */
    EliasFanoDecoder(std::ifstream *source,
                     std::ifstream *source_upper,
                     std::streampos file_end_upper);

    /** Returns the next compressed element or empty if all elements were read */
    std::optional<T> next();

  private:
    bool init();

    /** Returns the upper part of the next compressed element */
    T next_upper();

    /** Returns the lower part of the next compressed element */
    T next_lower();

    /** Clear the high bits of #value starting at position #index. */
    static T clear_high_bits(T value, uint8_t index);

  private:
    /** Index of current element */
    size_t position_ = 0;

    /** Current position in the #upper_ vector */
    uint64_t upper_pos_;

    /** The sequence of 8 upper bytes currently being processed */
    T upper_block_;

    /** Number of lower bits that were read from disk. */
    size_t cur_pos_bits_;

    /**
     * The lower bits of the encoded number, obtained by simply concatenating the
     * binary representation of the lower bits of each number. To save memory, only the
     * currently needed window of 16 bytes is read from the file.
     */
    T lower_[READ_BUF_SIZE];

    /** Points to current element in #lower_ */
    uint32_t lower_idx_;

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

    /** The size in bytes of lower_ */
    size_t num_lower_bytes_;

    /** The size in bytes of upper_ */
    size_t num_upper_bytes_ = 0;

    /**
     * Stream containing the compressed data. Points to either #source_internal_ or to
     * a stream provided in the constructor
     */
    std::ifstream *source_;
    /**
     * Stream containing the upper bytes of the compressed data (saved separately to avoid
     * costly tellg/seekg operations in the file.
     */
    std::ifstream *source_upper_;

    std::ifstream source_internal_;
    std::ifstream source_internal_upper_;

    std::streampos file_end_;

    /** Value to add to each decoded element */
    T offset_;
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

    EliasFanoEncoder(size_t size,
                     const T &first_value,
                     const T &last_value,
                     const std::string &sink_name,
                     bool is_append = false);

    void add(const std::pair<T, C> &value);

    size_t finish();

  private:
    EliasFanoEncoder<T> ef_encoder;
    std::string sink_second_name_;
    std::ofstream sink_second_;
};

/** Decoder specialization for an std::pair */
template <typename T, typename C>
class EliasFanoDecoder<std::pair<T, C>> {
  public:
    EliasFanoDecoder(const std::string &source);


    std::optional<std::pair<T, C>> next();

  private:
    EliasFanoDecoder<T> source_first_;
    std::ifstream source_second_;
};

/**
 * Specialization of #EliasFanoEncoder that can encode sequences of unknown range. It uses
 * a buffer to accumulate data and then dumps it in chunks to an EliasFanoEncoder.
 */
template <typename T>
class EliasFanoEncoderBuffered {
  public:
    EliasFanoEncoderBuffered(const std::string &file_name, size_t buffer_size);

    void add(const T &value);

    size_t finish();
  private:
    void encode_chunk();
  private:
    Vector<T> buffer_;
    std::string file_name_;
    size_t total_size_ = 0;
};

/**
 * Specialization of #EliasFanoEncoder that can encode sequences of pairs of unknown size.
 * It uses a buffer to accumulate data and then dumps it in chunks to an EliasFanoEncoder.
 */
// TODO: Pack the C count, by passing max count to the constructor and then dump it to
// sdsl::int_vector_buffer with width = hi(max)+1. Add a CompressedInt wrapper that
// encodes T with EliasFano and C (if present) with sdsl packing
template <typename T, typename C>
class EliasFanoEncoderBuffered<std::pair<T, C>> {
  public:
    EliasFanoEncoderBuffered(const std::string &file_name, size_t buffer_size);

    void add(const std::pair<T, C> &value);

    size_t finish();

  private:
    void encode_chunk();

  private:
    Vector<T> buffer_;
    Vector<C> buffer_second_;
    std::ofstream sink_;
    std::ofstream sink_upper_;
    std::ofstream sink_second_;
    std::string file_name_;
    size_t total_size_ = 0;
};

} // namespace common
} // namespace mg
