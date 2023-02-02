#ifndef __SEQUENCE_IO__
#define __SEQUENCE_IO__

#include <functional>
#include <vector>
#include <string>

#include <zlib.h>
#include <htslib/kseq.h>

#include "common/batch_accumulator.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/threads/threading.hpp"


namespace mtg {
namespace seq_io {

KSEQ_DECLARE(gzFile);


class FastaWriter {
  public:
    FastaWriter(const std::string &filebase,
                const std::string &header = "",
                bool enumerate_sequences = false,
                bool async = false,
                const char *mode = "w");

    ~FastaWriter();

    void write(const std::string &sequence);
    void write(std::string&& sequence);

    void flush();

  private:
    void write_to_disk(const std::string &sequence);

    gzFile gz_out_;
    const std::string header_;
    bool enumerate_sequences_;
    uint64_t count_ = 0;
    ThreadPool worker_;
    BatchAccumulator<std::string> batcher_;
};

template <typename T = uint32_t>
class ExtendedFastaWriter {
  public:
    typedef T feature_type;
    typedef std::pair<std::string, std::vector<feature_type>> value_type;

    /**
     * Write sequences to fasta file `<filebase>.fasta.gz` and dump
     * features for their respective k-mers (such as k-mer counts)
     * to `<filebase>.<feature_name>.gz`.
     * The features are dumped to a compressed array of `T`.
     * Each sequence dumped must contain at least one k-mer, that
     * is, must be at least |kmer_length| in length.
     */
    ExtendedFastaWriter(const std::string &filebase,
                        const std::string &feature_name,
                        uint32_t kmer_length,
                        const std::string &header = "",
                        bool enumerate_sequences = false,
                        bool async = false,
                        const char *mode = "w");

    ~ExtendedFastaWriter();

    /**
     * Features are compressed with gzip before dumping to the target file.
     */
    void write(const std::string &sequence,
               const std::vector<feature_type> &kmer_features);
    void write(std::string&& sequence,
               std::vector<feature_type>&& kmer_features);

    void flush();

  private:
    void write_to_disk(const value_type &value_pair);

    gzFile fasta_gz_out_;
    gzFile feature_gz_out_;
    uint32_t kmer_length_;
    const std::string header_;
    bool enumerate_sequences_;
    uint64_t count_ = 0;
    ThreadPool worker_;
    BatchAccumulator<value_type> batcher_;
};

bool write_fasta(gzFile gz_out, const kseq_t &kseq);

bool write_fasta(gzFile gz_out, const std::string &header,
                                const std::string &sequence);

bool write_fastq(gzFile gz_out, const kseq_t &kseq);

void read_fasta_file_critical(const std::string &filename,
                              std::function<void(kseq_t*)> callback,
                              bool with_reverse = false);

/**
 * Read fasta/fastq file with sequences and parse their corresponding k-mer's
 * reatures from <filebase>.<feature_name>.gz written by ExtendedFastaWriter.
 */
template <typename T>
void read_extended_fasta_file_critical(const std::string &filebase,
                                       const std::string &feature_name,
                                       std::function<void(size_t k, const kseq_t*, const T*)> callback,
                                       bool with_reverse = false);

void read_vcf_file_critical(const std::string &filename,
                            const std::string &ref_filename,
                            size_t k,
                            std::function<void(std::string_view)> callback,
                            bool with_reverse = false);

void read_vcf_file_with_annotations_critical(const std::string &filename,
                                             const std::string &ref_filename,
                                             size_t k,
                                             std::function<void(std::string&&,
                                                                const std::vector<std::string>&)> callback,
                                             bool with_reverse = false);

void read_fasta_from_string(const std::string &fasta_flat,
                            std::function<void(kseq_t*)> callback,
                            bool with_reverse = false);


class FastaParser {
  public:
    class iterator;

    FastaParser(const std::string &filename,
                bool with_reverse_complement = false)
      : filename_(filename),
        with_reverse_complement_(with_reverse_complement) {}

    inline iterator begin() const;
    inline iterator end() const;

    const std::string& get_filename() const { return filename_; }

  private:
    std::string filename_;
    bool with_reverse_complement_;
};

class FastaParser::iterator {
    friend FastaParser;

  public:
    using iterator_category = std::input_iterator_tag;
    using value_type = kseq_t;
    using difference_type = size_t;
    using pointer = kseq_t*;
    using reference = kseq_t&;

    iterator() {}

    iterator(const iterator &other) { *this = other; }
    iterator(iterator&& other) { *this = std::move(other); }

    iterator& operator=(const iterator &other);
    iterator& operator=(iterator&& other);

    ~iterator() { deinit_stream(); }

    kseq_t& operator*() { return *read_stream_; }
    const kseq_t& operator*() const { return *read_stream_; }

    kseq_t* operator->() { return read_stream_; }
    const kseq_t* operator->() const { return read_stream_; }

    iterator& operator++() {
        if (with_reverse_complement_ && !is_reverse_complement_) {
            reverse_complement(read_stream_->seq);
            is_reverse_complement_ = true;
            return *this;
        }

        if (kseq_read(read_stream_) < 0) {
            deinit_stream();
        }
        is_reverse_complement_ = false;
        return *this;
    }

    bool operator==(const iterator &other) const {
        return (read_stream_ && other.read_stream_
                    && is_reverse_complement_ == other.is_reverse_complement_
                    && read_stream_->f->seek_pos == other.read_stream_->f->seek_pos
                    && with_reverse_complement_ == other.with_reverse_complement_
                    && filename_ == other.filename_)
            || (!read_stream_ && !other.read_stream_);
    }

    bool operator!=(const iterator &other) const {
        return !(*this == other);
    }

  private:
    iterator(const std::string &filename, bool with_reverse_complement);
    void deinit_stream();

    std::string filename_;
    bool with_reverse_complement_;
    kseq_t *read_stream_ = NULL;
    bool is_reverse_complement_ = false;
};

FastaParser::iterator
FastaParser::begin() const { return iterator(filename_, with_reverse_complement_); }

FastaParser::iterator
FastaParser::end() const { return iterator(); }

} // namespace seq_io
} // namespace mtg

#endif // __SEQUENCE_IO__
