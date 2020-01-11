#ifndef __SEQUENCE_IO__
#define __SEQUENCE_IO__

#include <functional>
#include <vector>
#include <string>

#include <zlib.h>
#include <htslib/kseq.h>

KSEQ_INIT(gzFile, gzread);


class FastaWriter {
  public:
    FastaWriter(const std::string &filebase,
                const std::string &header = "",
                bool enumerate_sequences = false);

    ~FastaWriter();

    void write(const std::string &sequence);

  public:
    gzFile gz_out_;
    const std::string header_;
    bool enumerate_sequences_;
    uint64_t count_ = 0;
};

template <typename T = uint32_t>
class ExtendedFastaWriter {
  public:
    typedef T feature_type;

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
                        bool enumerate_sequences = false);

    ~ExtendedFastaWriter();

    /**
     * Features are compressed with gzip before dumping to the target file.
     */
    void write(const std::string &sequence,
               const std::vector<feature_type> &kmer_features);

  public:
    gzFile fasta_gz_out_;
    gzFile feature_gz_out_;
    uint32_t kmer_length_;
    const std::string header_;
    bool enumerate_sequences_;
    uint64_t count_ = 0;
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
                            std::function<void(std::string&&)> callback,
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

#endif // __SEQUENCE_IO__
