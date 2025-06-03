#ifndef __SEQUENCE_IO__
#define __SEQUENCE_IO__

#include <functional>
#include <vector>
#include <string>
#include <variant>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zstd.hpp>

#include <zlib.h>
#include <htslib/kseq.h>

#include "common/batch_accumulator.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/threads/threading.hpp"


namespace mtg {
namespace seq_io {

struct compFile {
    using bfi = boost::iostreams::filtering_istream;
    using bfo = boost::iostreams::filtering_ostream;

    ~compFile() { close(); }

    bool good() const {
        return std::visit([&](const auto &f) {
            if constexpr(std::is_same_v<std::decay_t<decltype(f)>, gzFile>) {
                return f != Z_NULL;
            } else {
                return f && f->good();
            }
        }, f_);
    }

    void match_offset(const compFile &other) {
        if (f_.index() != other.f_.index()) {
            throw std::runtime_error("Matching offsets for different stream types");
        }

        std::visit([&](auto &f) {
            std::visit([&](const auto &other_f) {
                if constexpr(std::is_same_v<std::decay_t<decltype(f)>, gzFile>) {
                    if constexpr(std::is_same_v<std::decay_t<decltype(other_f)>, gzFile>) {
                        gzseek(f, gztell(other_f), SEEK_SET);
                    } else {
                        assert(false && "Matching offsets for different stream types");
                    }
                } else if constexpr(std::is_same_v<std::decay_t<decltype(f)>, std::shared_ptr<bfi>>) {
                    if constexpr(std::is_same_v<std::decay_t<decltype(other_f)>, std::shared_ptr<bfi>>) {
                        f->seekg(other_f->tellg(), std::ios::beg);
                    } else {
                        assert(false && "Matching offsets for different stream types");
                    }
                }
            }, other.f_);
        }, f_);
    }

    int put(char c) {
        return std::visit([&](auto &f) {
            if constexpr(std::is_same_v<std::decay_t<decltype(f)>, gzFile>) {
                return gzputc(f, c);
            } else if constexpr(std::is_same_v<std::decay_t<decltype(f)>, std::shared_ptr<bfo>>) {
                f->put(c);
                return !f->eof() ? c : f->eof();
            } else {
                return -1;
            }
        }, f_);
    }

    int get() {
        return std::visit([&](auto &f) {
            if constexpr(std::is_same_v<std::decay_t<decltype(f)>, gzFile>) {
                return gzgetc(f);
            } else if constexpr(std::is_same_v<std::decay_t<decltype(f)>, std::shared_ptr<bfi>>) {
                return f->get();
            } else {
                return -1;
            }
        }, f_);
    }

    std::streamsize read(voidp buf, unsigned len) {
        return std::visit([&](auto &f) {
            if constexpr(std::is_same_v<std::decay_t<decltype(f)>, gzFile>) {
                return static_cast<std::streamsize>(gzread(f, buf, len));
            } else if constexpr(std::is_same_v<std::decay_t<decltype(f)>, std::shared_ptr<bfi>>) {
                f->read(reinterpret_cast<char*>(buf), len);
                return f->gcount();
            } else {
                return std::streamsize(0);
            }
        }, f_);
    }

    std::streamsize write(voidpc buf, unsigned len) {
        return std::visit([&](auto &f) {
            if constexpr(std::is_same_v<std::decay_t<decltype(f)>, gzFile>) {
                return static_cast<std::streamsize>(gzwrite(f, buf, len));
            } else if constexpr(std::is_same_v<std::decay_t<decltype(f)>, std::shared_ptr<bfo>>) {
                f->write(reinterpret_cast<const char*>(buf), len);
                return static_cast<std::streamsize>(len);
            } else {
                return std::streamsize(0);
            }
        }, f_);
    }

    int close() {
        return std::visit([&](auto &f) {
            if constexpr(std::is_same_v<std::decay_t<decltype(f)>, gzFile>) {
                return gzclose(f);
            } else {
                try {
                    f.reset();
                    if (cf_) {
                        int ret_val = fclose(cf_);
                        cf_ = NULL;
                        return ret_val;
                    }

                    return Z_OK;
                } catch (...) {
                    return Z_ERRNO;
                }
            }
        }, f_);
    }

    bool eof() const {
        return std::visit([&](const auto &f) {
            if constexpr(std::is_same_v<std::decay_t<decltype(f)>, gzFile>) {
                return static_cast<bool>(gzeof(f));
            } else {
                return f->eof();
            }
        }, f_);
    }

    std::string error() const {
        return std::visit([&](const auto &f) {
            int err_state;
            if constexpr(std::is_same_v<std::decay_t<decltype(f)>, gzFile>) {
                return gzerror(f, &err_state);
            }

            return "";
        }, f_);
    }

    static int read(compFile &f, voidp buf, unsigned len) {
        return f.read(buf, len);
    }

    static int write(compFile &f, voidpc buf, unsigned len) {
        return f.write(buf, len);
    }

    static compFile open_read(const char *path) {
        compFile f;
        const char *dotpos = strrchr(path, '.');
        if (dotpos == NULL)
            throw std::runtime_error("No file extension");

        std::string ext(dotpos);

        if (ext == ".zst") {
            auto fio = std::make_shared<bfi>();
            fio->push(boost::iostreams::zstd_decompressor());
            fio->push(boost::iostreams::file_descriptor(path));
            f.f_ = fio;
        } else {
            f.f_ = gzopen(path, "r");
        }

        return f;
    }

    static compFile open_write(const char *path) {
        compFile f;
        const char *dotpos = strrchr(path, '.');
        if (dotpos == NULL)
            throw std::runtime_error("No file extension");

        std::string ext(dotpos);

        if (ext == ".zst") {
            auto fio = std::make_shared<bfo>();
            fio->push(boost::iostreams::zstd_decompressor());
            fio->push(boost::iostreams::file_descriptor(path));
            f.f_ = fio;
        } else {
            f.f_ = gzopen(path, "w");
        }

        return f;
    }

    template <typename T>
    static compFile dopen_read(int fd) {
        compFile f;
        if constexpr(std::is_same_v<T, gzFile>) {
            f.f_ = gzdopen(fd, "r");
        } else {
            auto fio = std::make_shared<bfi>();
            f.cf_ = fdopen(fd, "r");
            fio->push(boost::iostreams::zstd_decompressor());
            fio->push(f.cf_);
            f.f_ = fio;
        }
        return f;
    }

    std::variant<gzFile, std::shared_ptr<bfi>, std::shared_ptr<bfo>> f_;
    FILE *cf_ = NULL;
};

KSEQ_DECLARE(compFile);


class FastaWriter {
  public:
    FastaWriter(const std::string &filebase,
                const std::string &header = "",
                bool enumerate_sequences = false,
                bool async = false,
                const char *mode = "w",
                bool zstd = false);

    ~FastaWriter();

    void write(const std::string &sequence);
    void write(std::string&& sequence);

    void flush();

    const std::string& get_fname() const { return fname_; }

  private:
    void write_to_disk(const std::string &sequence);

    std::string fname_;
    compFile out_;
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
                        const char *mode = "w",
                        bool zstd = false);

    ~ExtendedFastaWriter();

    /**
     * Features are compressed with gzip before dumping to the target file.
     */
    void write(const std::string &sequence,
               const std::vector<feature_type> &kmer_features);
    void write(std::string&& sequence,
               std::vector<feature_type>&& kmer_features);

    void flush();

    const std::string& get_fasta_fname() const { return fasta_fname_; }

  private:
    void write_to_disk(const value_type &value_pair);

    std::string fasta_fname_;
    compFile fasta_out_;
    compFile feature_out_;
    uint32_t kmer_length_;
    const std::string header_;
    bool enumerate_sequences_;
    uint64_t count_ = 0;
    ThreadPool worker_;
    BatchAccumulator<value_type> batcher_;
};

bool write_fasta(compFile &out, const kseq_t &kseq);

bool write_fasta(compFile &out, const std::string &header,
                                const std::string &sequence);

bool write_fastq(compFile &out, const kseq_t &kseq);

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
