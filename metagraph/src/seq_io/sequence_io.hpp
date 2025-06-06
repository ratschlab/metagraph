#ifndef __SEQUENCE_IO__
#define __SEQUENCE_IO__

#include <functional>
#include <vector>
#include <string>
#include <variant>
#include <fstream>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <htslib/kseq.h>

#include "common/batch_accumulator.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/threads/threading.hpp"


namespace mtg {
namespace seq_io {

class compFile {
  public:
    bool good() const {
        return std::visit([&](const auto &f) {
            return f && f->good();
        }, f_);
    }

    int put(char c) {
        return std::visit([&](auto &f) {
            if (!f)
                return -1;

            if constexpr(std::is_base_of_v<std::ostream, std::decay_t<decltype(*f)>>) {
                if (!*f)
                    return -1;

                f->put(c);
                return !*f ? -1 : c;
            }

            return -1;
        }, f_);
    }

    int get() {
        return std::visit([&](auto &f) {
            if (!f)
                return -1;

            if constexpr(std::is_base_of_v<std::istream, std::decay_t<decltype(*f)>>) {
                return !*f ? -1 : f->get();
            }

            return -1;
        }, f_);
    }

    int read(void *buf, unsigned len) {
        return std::visit([&](auto &f) -> int {
            if (!f)
                return -1;

            if constexpr(std::is_base_of_v<std::istream, std::decay_t<decltype(*f)>>) {
                if (f->eof())
                    return 0;

                if (f->bad() || f->fail())
                    return -1;

                f->read(reinterpret_cast<char*>(buf), len);
                int ret_val = f->gcount();
                if (f->good() || f->eof())
                    return ret_val;
            }

            return -1;
        }, f_);
    }

    int write(const void *buf, unsigned len) {
        return std::visit([&](auto &f) -> int {
            if (!f)
                return -1;

            if constexpr(std::is_base_of_v<std::ostream, std::decay_t<decltype(*f)>>) {
                if (*f) {
                    f->write(reinterpret_cast<const char*>(buf), len);
                    if (*f)
                        return static_cast<int>(len);
                }
            }

            return -1;
        }, f_);
    }

    bool eof() const {
        return std::visit([&](const auto &f) {
            return f && f->eof();
        }, f_);
    }

    std::string error() const { return ""; }

    static int read(compFile &f, void *buf, unsigned len) {
        return f.read(buf, len);
    }

    static int write(compFile &f, const void *buf, unsigned len) {
        return f.write(buf, len);
    }

    static compFile open_read(const char *path) {
        compFile f;
        const char *dotpos = strrchr(path, '.');
        if (dotpos == NULL)
            throw std::runtime_error("No file extension");

        std::string ext(dotpos);

        if (ext == ".zst") {
            auto fio = std::make_shared<boost::iostreams::filtering_istream>();
            fio->push(boost::iostreams::zstd_decompressor());
            fio->push(boost::iostreams::file_source(path, std::ios::binary));
            f.f_.emplace<std::shared_ptr<std::istream>>(fio);
        } else if (ext == ".gz" || ext == ".bgz") {
            auto fio = std::make_shared<boost::iostreams::filtering_istream>();
            fio->push(boost::iostreams::gzip_decompressor());
            fio->push(boost::iostreams::file_source(path, std::ios::binary));
            f.f_.emplace<std::shared_ptr<std::istream>>(fio);
        } else {
            f.f_.emplace<std::shared_ptr<std::istream>>(
                std::make_shared<std::ifstream>(path)
            );
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
            auto foo = std::make_shared<boost::iostreams::filtering_ostream>();
            foo->push(boost::iostreams::zstd_compressor());
            foo->push(boost::iostreams::file_sink(path, std::ios::binary));
            f.f_.emplace<std::shared_ptr<std::ostream>>(foo);
        } else if (ext == ".gz") {
            auto foo = std::make_shared<boost::iostreams::filtering_ostream>();
            foo->push(boost::iostreams::gzip_compressor());
            foo->push(boost::iostreams::file_sink(path, std::ios::binary));
            f.f_.emplace<std::shared_ptr<std::ostream>>(foo);
        } else {
            f.f_.emplace<std::shared_ptr<std::ostream>>(
                std::make_shared<std::ofstream>(path)
            );
        }

        return f;
    }

    static compFile open_read(std::istream &in) {
        compFile f;
        f.f_ = std::shared_ptr<std::istream>(std::shared_ptr<std::istream>{}, &in);
        return f;
    }

    std::shared_ptr<std::ios> data() {
        return std::visit([&](auto &f) -> std::shared_ptr<std::ios> {
            return f;
        }, f_);
    }

  private:
    std::variant<std::shared_ptr<std::istream>, std::shared_ptr<std::ostream>> f_;
};

KSEQ_DECLARE(compFile);


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

    struct deinit_stream { void operator()(kseq_t *read_stream); };

  public:
    using iterator_category = std::input_iterator_tag;
    using stream_type = std::unique_ptr<kseq_t, deinit_stream>;
    using value_type = stream_type::element_type;
    using difference_type = size_t;
    using pointer = stream_type::pointer;
    using reference = stream_type::element_type&;

    iterator() {}

    // some fasta formats are not random access, so we can't copy iterators
    iterator& operator=(const iterator&) = delete;
    iterator(const iterator&) = delete;

    iterator& operator=(iterator&&) = default;
    iterator(iterator&&) = default;

    kseq_t& operator*() { return *read_stream_; }
    const kseq_t& operator*() const { return *read_stream_; }

    kseq_t* operator->() { return read_stream_.get(); }
    const kseq_t* operator->() const { return read_stream_.get(); }

    iterator& operator++() {
        if (!read_stream_ || !read_stream_->f->f.data()) {
            *this = iterator();
            return *this;
        }

        if (with_reverse_complement_ && !is_reverse_complement_) {
            reverse_complement(read_stream_->seq);
            is_reverse_complement_ = true;
            return *this;
        }

        if (kseq_read(read_stream_.get()) < 0)
            *this = iterator();

        is_reverse_complement_ = false;
        return *this;
    }

    bool operator==(const iterator &other) const {
        return (read_stream_.get() && other.read_stream_.get()
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

    std::string filename_;
    bool with_reverse_complement_;
    stream_type read_stream_;
    bool is_reverse_complement_ = false;
};

FastaParser::iterator
FastaParser::begin() const { return iterator(filename_, with_reverse_complement_); }

FastaParser::iterator FastaParser::end() const { return iterator(); }

} // namespace seq_io
} // namespace mtg

#endif // __SEQUENCE_IO__
