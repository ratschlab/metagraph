#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <functional>
#include <optional>
#include <vector>

#include <sdsl/uint128_t.hpp>
#include <sdsl/uint256_t.hpp>


namespace mtg {
namespace common {

/**
 * Concatenates Elias-Fano files into result.
 * The files store data that is ordered and the values in a file are smaller than the
 * values in the next file.
 */
void concat(const std::vector<std::string> &files, const std::string &result);

void remove_chunks(const std::vector<std::string> &files);

/**
 * Elias-Fano encoder that streams the encoded result into a file.
 * Loosely inspired  by
 * https://github.com/facebook/folly/blob/master/folly/experimental/EliasFanoCoding.h
 */
template <typename T>
class EliasFanoEncoder {
  public:
    static constexpr uint32_t WRITE_BUF_SIZE = 1024;

    EliasFanoEncoder(EliasFanoEncoder&&) = delete;
    EliasFanoEncoder& operator=(EliasFanoEncoder&&) = delete;

    /**
     * Constructs an Elias-Fano encoder of an array with the given #size and given
     * #max_value. The encoded output is written to #out_filename.
     */
    EliasFanoEncoder(size_t size,
                     T min_value,
                     T max_value,
                     const std::string &out_filename,
                     bool append = false);

    /** Constructs an encoder that encodes the #data array */
    EliasFanoEncoder(const std::vector<T> &data, std::ofstream *sink, std::ofstream *sink_upper);

    ~EliasFanoEncoder();

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
    static void write_bits(char *data, size_t pos, T value);

    void init(size_t size, T max_value);

  private:
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
    static_assert( std::is_integral_v<T> || std::is_same_v<T, sdsl::uint256_t>);

  public:
    EliasFanoDecoder() {}

    /** Creates a decoder that retrieves data from the given file */
    EliasFanoDecoder(const std::string &source_name, bool remove_source = true);

    /** Returns the next compressed element or empty if all elements were read */
    inline std::optional<T> next() {
        if (buffer_pos_ == buffer_end_) {
            if (!decompress_next_block()) {
                std::optional<T> no_value;
                return no_value;
            }
        }
        return buffer_[buffer_pos_++];
    }

  private:
    bool init();

    /** Decompressed the next block into buffer and returns the number of elements in it */
    size_t decompress_next_block();

  private:
    /** Index of current element */
    size_t position_ = 0;

    /**
     * Buffer with uncompressed elements.
     */
    T buffer_[READ_BUF_SIZE];
    size_t buffer_pos_ = 0;
    size_t buffer_end_ = 0;

    /** Current position in the #upper_ vector */
    uint64_t upper_pos_;

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
    std::vector<uint64_t> upper_;

    /** Total number of elements encoded. */
    size_t size_ = 0;

    /**
     * Each encoded integer is split into a "lower" and an "upper" part. This is the
     * number of bits used for the "lower" part of the Elias-Fano encoding.
     */
    uint8_t num_lower_bits_;
    T lower_bits_mask_;

    /** The size in bytes of lower_ */
    size_t num_lower_bytes_;

    /** The size in bytes of upper_ */
    size_t num_upper_bytes_ = 0;

    /** Name of the source file. Only set if #source_internal_ is set */
    std::string source_name_;

    /**
     * Stream containing the compressed data.
     */
    std::ifstream source_;
    /**
     * Stream containing the upper bytes of the compressed data (saved separately to avoid
     * costly tellg/seekg operations in the file.
     */
    std::ifstream source_upper_;

    /** Value to add to each decoded element */
    T offset_;

    /** If true, the source file is removed after decompression */
    bool remove_source_;
};

/**
 * Encoder specialization for an std::pair. The first member of the pair is assumed to be
 * in nondecreasing order and is compressed using Elias-Fano encoding, while the second
 * member of the pair is written to disk as is.
 */
template <typename T, typename C>
class EliasFanoEncoder<std::pair<T, C>> {
  public:
    EliasFanoEncoder(size_t size,
                     const T &first_value,
                     const T &last_value,
                     const std::string &sink_name,
                     bool append = false);

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
    EliasFanoDecoder(const std::string &source, bool remove_source = true);

    inline std::optional<std::pair<T, C>> next() {
        std::optional<T> first = source_first_.next();
        C second;
        if (!first.has_value()) {
            assert(!source_second_.read(reinterpret_cast<char *>(&second), sizeof(C)));
            if (remove_source_) {
                std::filesystem::remove(source_second_name_);
            }
            return {};
        }
        source_second_.read(reinterpret_cast<char *>(&second), sizeof(C));
        assert(source_second_);
        return std::make_pair(first.value(), second);
    }

  private:
    EliasFanoDecoder<T> source_first_;
    std::string source_second_name_;
    std::ifstream source_second_;
    bool remove_source_;
};

/**
 * Specialization of #EliasFanoEncoder that can encode sequences of unknown range. It uses
 * a buffer to accumulate data and then dumps it in chunks to an EliasFanoEncoder.
 */
template <typename T>
class EliasFanoEncoderBuffered {
  public:
    EliasFanoEncoderBuffered(const std::string &file_name, size_t buffer_size);

    EliasFanoEncoderBuffered(EliasFanoEncoderBuffered&&) = default;
    EliasFanoEncoderBuffered& operator=(EliasFanoEncoderBuffered&&) = default;

    ~EliasFanoEncoderBuffered();

    inline void add(const T &value) {
        buffer_.push_back(value);
        if (buffer_.size() == buffer_.capacity()) {
            encode_chunk();
        }
        size_++;
    }

    const std::string& name() { return file_name_; }

    size_t size() const { return size_; }

    size_t finish();

  private:
    void encode_chunk();
  private:
    std::vector<T> buffer_;
    std::ofstream sink_;
    std::ofstream sink_upper_;
    std::string file_name_;
    size_t total_size_ = 0; // total encoded size in bytes
    size_t size_ = 0; // total number of encoded elements
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

    EliasFanoEncoderBuffered(EliasFanoEncoderBuffered&&) = default;
    EliasFanoEncoderBuffered& operator=(EliasFanoEncoderBuffered&&) = default;

    ~EliasFanoEncoderBuffered();

    inline void add(const std::pair<T, C> &value) {
        buffer_.push_back(value.first);
        buffer_second_.push_back(value.second);
        if (buffer_.size() == buffer_.capacity()) {
            encode_chunk();
        }
    }

    const std::string& name() { return file_name_; }

    size_t finish();

  private:
    void encode_chunk();
  private:
    std::vector<T> buffer_;
    std::vector<C> buffer_second_;
    std::ofstream sink_;
    std::ofstream sink_upper_;
    std::ofstream sink_second_;
    std::string file_name_;
    size_t total_size_ = 0;
};

} // namespace common
} // namespace mtg
