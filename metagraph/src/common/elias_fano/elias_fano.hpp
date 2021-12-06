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

#include "common/logger.hpp"


namespace mtg {
namespace elias_fano {

/**
 * Concatenates Elias-Fano files into result.
 * The files store data that is ordered and the values in a file are smaller than the
 * values in the next file.
 */
void concat(const std::vector<std::string> &files, const std::string &result);

void remove_chunks(const std::vector<std::string> &files);

// get size in bytes
uint64_t chunk_size(const std::string &file);


/**
 * Decodes a list of compressed sorted integers stored in a file using #EliasFanoEncoder.
 */
template <typename T>
class EliasFanoDecoder {
    static constexpr uint32_t READ_BUF_SIZE = 1024;
    static_assert(std::is_integral_v<T> || std::is_same_v<T, sdsl::uint256_t>);

  public:
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

    /** If reading fails, retry |max_num_retries_| times until crashing */
    size_t num_retries_ = 0;
    static constexpr size_t max_num_retries_ = 100;

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

/** Decoder specialization for an std::pair */
template <typename T, typename C>
class EliasFanoDecoder<std::pair<T, C>> {
  public:
    EliasFanoDecoder(const std::string &source, bool remove_source = true);

    inline std::optional<std::pair<T, C>> next() {
        std::optional<T> first = source_first_.next();
        C second;
        source_second_.read(reinterpret_cast<char *>(&second), sizeof(C));
        if (!first.has_value()) {
            if (!source_second_.eof()) {
                common::logger->error("EliasFanoDecoder error: file {} is not read to the end",
                                      source_second_name_);
                std::exit(EXIT_FAILURE);
            }
            source_second_.close();
            if (remove_source_)
                std::filesystem::remove(source_second_name_);

            return {};
        }
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
    EliasFanoEncoderBuffered(const std::string &file_name,
                             size_t buffer_size,
                             bool append = false);

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

    /** Append sorted array #data to EF-coded #out_fname */
    static size_t append_block(const std::vector<T> &data,
                               const std::string &file_name);

  private:
    void encode_chunk();

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
    EliasFanoEncoderBuffered(const std::string &file_name,
                             size_t buffer_size,
                             bool append = false);

    EliasFanoEncoderBuffered(EliasFanoEncoderBuffered&&) = default;
    EliasFanoEncoderBuffered& operator=(EliasFanoEncoderBuffered&&) = default;

    ~EliasFanoEncoderBuffered();

    inline void add(const std::pair<T, C> &value) {
        encoder_first_.add(value.first);
        buffer_second_.push_back(value.second);
        if (buffer_second_.size() == buffer_second_.capacity()) {
            encode_chunk();
        }
    }

    const std::string& name() { return encoder_first_.name(); }

    size_t size() const { return encoder_first_.size(); }

    size_t finish();

    /** Append sorted array #data to EF-coded #out_fname */
    static size_t append_block(const std::vector<std::pair<T, C>> &data,
                               const std::string &file_name);

  private:
    void encode_chunk();

    EliasFanoEncoderBuffered<T> encoder_first_;
    std::vector<C> buffer_second_;
    std::ofstream sink_second_;
};

} // namespace elias_fano
} // namespace mtg
