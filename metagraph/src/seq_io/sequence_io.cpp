#include "sequence_io.hpp"

#include <iostream>
#include <fstream>
#include <thread>

#include "common/seq_tools/reverse_complement.hpp"
#include "common/utils/string_utils.hpp"
#include "vcf_parser.hpp"


namespace mtg {
namespace seq_io {

const char kDefaultFastQualityChar = 'I';

// Optimal values found from a grid search with the BM_WriteRandomSequences benchmark
const size_t kWorkerQueueSize = 1;
const size_t kBufferSize = 1'000'000;


FastaWriter::FastaWriter(const std::string &filebase,
                         const std::string &header,
                         bool enumerate_sequences,
                         bool async,
                         const char *mode)
      : header_(header),
        enumerate_sequences_(enumerate_sequences),
        worker_(async, kWorkerQueueSize) {
    auto filename = utils::remove_suffix(filebase, ".gz", ".fasta") + ".fasta.gz";

    gz_out_ = gzopen(filename.c_str(), mode);
    if (gz_out_ == Z_NULL) {
        std::cerr << "ERROR: Can't write to " << filename << std::endl;
        exit(1);
    }

    batcher_ = BatchAccumulator<std::string>(
        [&](std::vector<std::string>&& buffer) {
            worker_.enqueue([&](const auto &buffer) {
                                for (const std::string &sequence : buffer) {
                                    write_to_disk(sequence);
                                }
                            },
                            std::move(buffer));
        },
        kBufferSize / kWorkerQueueSize / sizeof(std::string),  // max size
        kBufferSize / kWorkerQueueSize  // max cumulative length
    );
}

FastaWriter::~FastaWriter() {
    flush();
    gzclose(gz_out_);
}

void FastaWriter::flush() {
    batcher_.process_all_buffered();
    worker_.join();
}

void FastaWriter::write(const std::string &sequence) {
    batcher_.push_and_pay(sequence.size(), sequence);
}

void FastaWriter::write(std::string&& sequence) {
    batcher_.push_and_pay(sequence.size(), std::move(sequence));
}

void FastaWriter::write_to_disk(const std::string &sequence) {
    if (!write_fasta(gz_out_,
                     enumerate_sequences_ ? header_ + std::to_string(++count_)
                                          : header_,
                     sequence)) {
        std::cerr << "ERROR: FastaWriter::write failed. Can't dump sequence to fasta" << std::endl;
        exit(1);
    }
}


template <typename T>
ExtendedFastaWriter<T>::ExtendedFastaWriter(const std::string &filebase,
                                            const std::string &feature_name,
                                            uint32_t kmer_length,
                                            const std::string &header,
                                            bool enumerate_sequences,
                                            bool async,
                                            const char *mode)
      : kmer_length_(kmer_length),
        header_(header),
        enumerate_sequences_(enumerate_sequences),
        worker_(async, kWorkerQueueSize) {
    assert(feature_name.size());

    auto filename = utils::remove_suffix(filebase, ".gz", ".fasta") + ".fasta.gz";

    fasta_gz_out_ = gzopen(filename.c_str(), mode);
    if (fasta_gz_out_ == Z_NULL) {
        std::cerr << "ERROR: Can't write to " << filename << std::endl;
        exit(1);
    }

    filename = utils::remove_suffix(filebase, ".gz", ".fasta") + "." + feature_name + ".gz";

    feature_gz_out_ = gzopen(filename.c_str(), "w");
    if (feature_gz_out_ == Z_NULL
            || gzwrite(feature_gz_out_, &kmer_length_, 4) != sizeof(kmer_length_)) {
        std::cerr << "ERROR: Can't write to " << filename << std::endl;
        exit(1);
    }

    batcher_ = BatchAccumulator<value_type>(
        [&](std::vector<value_type>&& buffer) {
            worker_.enqueue([&](const auto &buffer) {
                                for (const value_type &value_pair : buffer) {
                                    write_to_disk(value_pair);
                                }
                            },
                            std::move(buffer));
        },
        kBufferSize / kWorkerQueueSize / sizeof(value_type),  // max size
        kBufferSize / kWorkerQueueSize  // max cumulative length
    );
}

template <typename T>
ExtendedFastaWriter<T>::~ExtendedFastaWriter() {
    flush();
    gzclose(fasta_gz_out_);
    gzclose(feature_gz_out_);
}

template <typename T>
void ExtendedFastaWriter<T>::flush() {
    batcher_.process_all_buffered();
    worker_.join();
}

template <typename T>
void ExtendedFastaWriter<T>::write(const std::string &sequence,
                                   const std::vector<feature_type> &kmer_features) {
    assert(kmer_features.size() + kmer_length_ - 1 == sequence.size());

    batcher_.push_and_pay(sequence.size() + sizeof(feature_type) * kmer_features.size(),
                          std::make_pair(sequence, kmer_features));
}

template <typename T>
void ExtendedFastaWriter<T>::write(std::string&& sequence,
                                   std::vector<feature_type>&& kmer_features) {
    assert(kmer_features.size() + kmer_length_ - 1 == sequence.size());

    batcher_.push_and_pay(sequence.size() + sizeof(feature_type) * kmer_features.size(),
                          std::make_pair(std::move(sequence), std::move(kmer_features)));
}

template <typename T>
void ExtendedFastaWriter<T>::write_to_disk(const value_type &value_pair) {
    const auto &[sequence, kmer_features] = value_pair;
    if (!write_fasta(fasta_gz_out_,
                     enumerate_sequences_ ? header_ + std::to_string(++count_)
                                          : header_,
                     sequence)) {
        std::cerr << "ERROR: ExtendedFastaWriter::write failed. Can't dump sequence to fasta" << std::endl;
        exit(1);
    }

    if (gzwrite(feature_gz_out_, kmer_features.data(),
                                 kmer_features.size() * sizeof(feature_type))
            != static_cast<int>(kmer_features.size() * sizeof(feature_type))) {
        std::cerr << "ERROR: ExtendedFastaWriter::write failed. Can't dump k-mer features" << std::endl;
        exit(1);
    }
}

template class ExtendedFastaWriter<uint8_t>;
template class ExtendedFastaWriter<uint16_t>;
template class ExtendedFastaWriter<uint32_t>;
template class ExtendedFastaWriter<uint64_t>;


bool write_fasta(gzFile gz_out, const kseq_t &kseq) {
    return gzputc(gz_out, '>') == '>'
        && gzwrite(gz_out, kseq.name.s, kseq.name.l)
                        == static_cast<int>(kseq.name.l)
        && (!kseq.comment.l
            || (gzputc(gz_out, ' ') == ' '
                && gzwrite(gz_out, kseq.comment.s, kseq.comment.l)
                                    == static_cast<int>(kseq.comment.l)))
        && gzputc(gz_out, '\n') == '\n'
        && gzwrite(gz_out, kseq.seq.s, kseq.seq.l)
                        == static_cast<int>(kseq.seq.l)
        && gzputc(gz_out, '\n') == '\n';
}

bool write_fasta(gzFile gz_out, const std::string &header,
                                const std::string &sequence) {
    return gzputc(gz_out, '>') == '>'
        && gzwrite(gz_out, header.data(), header.length())
                        == static_cast<int>(header.length())
        && gzputc(gz_out, '\n') == '\n'
        && gzwrite(gz_out, sequence.data(), sequence.length())
                        == static_cast<int>(sequence.length())
        && gzputc(gz_out, '\n') == '\n';
}

bool write_fastq(gzFile gz_out, const kseq_t &kseq) {
    std::string qual(kseq.qual.s, kseq.qual.l);

    if (!kseq.qual.l && kseq.seq.l)
        qual.assign(kseq.seq.l, kDefaultFastQualityChar);

    return gzputc(gz_out, '@') == '@'
        && gzwrite(gz_out, kseq.name.s, kseq.name.l)
                        == static_cast<int>(kseq.name.l)
        && (kseq.comment.l == 0
            || (gzputc(gz_out, ' ') == ' '
                && gzwrite(gz_out, kseq.comment.s, kseq.comment.l)
                                    == static_cast<int>(kseq.comment.l)))
        && gzputc(gz_out, '\n') == '\n'
        && gzwrite(gz_out, kseq.seq.s, kseq.seq.l)
                        == static_cast<int>(kseq.seq.l)
        && gzputc(gz_out, '\n') == '\n'
        && gzputc(gz_out, '+') == '+'
        && gzwrite(gz_out, kseq.name.s, kseq.name.l)
                        == static_cast<int>(kseq.name.l)
        && (kseq.comment.l == 0
            || (gzputc(gz_out, ' ') == ' '
                && gzwrite(gz_out, kseq.comment.s, kseq.comment.l)
                                    == static_cast<int>(kseq.comment.l)))
        && gzputc(gz_out, '\n') == '\n'
        && gzwrite(gz_out, qual.data(), qual.size())
                        == static_cast<int>(qual.size())
        && gzputc(gz_out, '\n') == '\n';
}


FastaParser::iterator& FastaParser::iterator::operator=(const iterator &other) {
    if (!other.read_stream_) {
        // free memory, reset members and return
        deinit_stream();
        filename_.clear();
        with_reverse_complement_ = false;
        read_stream_ = NULL;
        is_reverse_complement_ = false;
        return *this;
    }

    if (!read_stream_) {
        *this = iterator(other.filename_, other.with_reverse_complement_);

    } else if (filename_ != other.filename_) {
        // The current iterator doesn't point to the target fasta file.
        // Thus, we need to open a new descriptor and seek to the target position.

        filename_ = other.filename_;

        // close the old file descriptor
        gzclose(read_stream_->f->f);

        // open a new file descriptor
        read_stream_->f->f = gzopen(filename_.c_str(), "r");
        if (read_stream_->f->f == Z_NULL) {
            std::cerr << "ERROR: Cannot read from file " << filename_ << std::endl;
            exit(1);
        }
    }

    gzseek(read_stream_->f->f, gztell(other.read_stream_->f->f), SEEK_SET);

#define KSTRING_COPY(kstr, other_kstr) \
            kstr.l = other_kstr.l; \
            if (kstr.m != other_kstr.m) { \
                kstr.m = other_kstr.m; \
                kstr.s = (char*)realloc(kstr.s, kstr.m); \
                if (!kstr.s) { \
                    std::cerr << "ERROR: realloc failed" << std::endl; \
                    exit(1); \
                } \
            } \
            memcpy(kstr.s, other_kstr.s, other_kstr.m); \

    // copy last cached record
    KSTRING_COPY(read_stream_->name, other.read_stream_->name);
    KSTRING_COPY(read_stream_->comment, other.read_stream_->comment);
    KSTRING_COPY(read_stream_->seq, other.read_stream_->seq);
    KSTRING_COPY(read_stream_->qual, other.read_stream_->qual);

    read_stream_->last_char = other.read_stream_->last_char;

    // copy stream state
    kstream_t *f = read_stream_->f;
    kstream_t *other_f = other.read_stream_->f;

    f->begin = other_f->begin;
    f->end = other_f->end;
    f->is_eof = other_f->is_eof;
    f->seek_pos = other_f->seek_pos;
    if (f->bufsize != other_f->bufsize) {
        f->bufsize = other_f->bufsize;
        f->buf = (unsigned char*)realloc(f->buf, f->bufsize);
        if (!f->buf) {
            std::cerr << "ERROR: realloc failed" << std::endl;
            exit(1);
        }
    }
    memcpy(f->buf, other_f->buf, other_f->bufsize);

    is_reverse_complement_ = other.is_reverse_complement_;

    return *this;
}

FastaParser::iterator& FastaParser::iterator::operator=(iterator&& other) {
    std::swap(filename_, other.filename_);
    with_reverse_complement_ = other.with_reverse_complement_;
    std::swap(read_stream_, other.read_stream_);
    is_reverse_complement_ = other.is_reverse_complement_;
    // the destructor in |other| will be responsible for freeing the memory now
    return *this;
}

FastaParser::iterator::iterator(const std::string &filename,
                                bool with_reverse_complement)
      : filename_(filename),
        with_reverse_complement_(with_reverse_complement) {
    gzFile input_p = gzopen(filename_.c_str(), "r");
    if (input_p == Z_NULL) {
        std::cerr << "ERROR: Cannot read from file " << filename_ << std::endl;
        exit(1);
    }

    read_stream_ = kseq_init(input_p);
    if (read_stream_ == NULL) {
        std::cerr << "ERROR: failed to initialize kseq file descriptor" << std::endl;
        exit(1);
    }

    if (kseq_read(read_stream_) < 0) {
        deinit_stream();
    }
}

void FastaParser::iterator::deinit_stream() {
    if (read_stream_) {
        gzFile input_p = read_stream_->f->f;
        kseq_destroy(read_stream_);
        gzclose(input_p);
    }
    read_stream_ = NULL;
}


template <class Callback>
void read_fasta_file_critical(gzFile input_p,
                              Callback callback,
                              bool with_reverse) {
    size_t seq_count = 0;

    if (input_p == Z_NULL) {
        std::cerr << "ERROR: Null file descriptor" << std::endl;
        exit(1);
    }
    //TODO: handle read_stream->qual
    kseq_t *read_stream = kseq_init(input_p);
    if (read_stream == NULL) {
        std::cerr << "ERROR: failed to initialize kseq file descriptor" << std::endl;
        exit(1);
    }

    while (kseq_read(read_stream) >= 0) {
        callback(read_stream);
        if (with_reverse) {
            reverse_complement(read_stream->seq);
            callback(read_stream);
        }

        seq_count++;
    }

    kseq_destroy(read_stream);
}

void read_fasta_file_critical(const std::string &filename,
                              std::function<void(kseq_t*)> callback,
                              bool with_reverse) {
    gzFile input_p = gzopen(filename.c_str(), "r");
    if (input_p == Z_NULL) {
        std::cerr << "ERROR: Cannot read file " << filename << std::endl;
        exit(1);
    }

    read_fasta_file_critical(input_p, callback, with_reverse);

    gzclose(input_p);
}

template <typename T>
void read_extended_fasta_file_critical(const std::string &filebase,
                                       const std::string &feature_name,
                                       std::function<void(size_t, const kseq_t*, const T*)> callback,
                                       bool with_reverse) {
    assert(feature_name.size());

    auto filename = utils::remove_suffix(filebase, ".gz", ".fasta") + ".fasta.gz";

    gzFile fasta_p = gzopen(filename.c_str(), "r");
    if (fasta_p == Z_NULL) {
        std::cerr << "ERROR: Cannot read from file " << filename << std::endl;
        exit(1);
    }

    filename = utils::remove_suffix(filebase, ".gz", ".fasta") + "." + feature_name + ".gz";
    uint32_t kmer_length;

    gzFile features_p = gzopen(filename.c_str(), "r");
    if (features_p == Z_NULL
            || gzread(features_p, &kmer_length, 4) != sizeof(kmer_length)) {
        std::cerr << "ERROR: Cannot read from file " << filename << std::endl;
        exit(1);
    }

    std::vector<T> counts;
    // read sequences from fasta file
    read_fasta_file_critical(fasta_p,
        [&](kseq_t *read_stream) {
            if (read_stream->seq.l < kmer_length) {
                std::cerr << "ERROR: Bad fasta file. Found sequence shorter"
                             " than k-mer length " << kmer_length << std::endl;
                exit(1);
            }
            counts.resize(read_stream->seq.l - kmer_length + 1);
            // read the features of k-mers in sequence |read_stream|
            if (gzread(features_p, counts.data(), sizeof(T) * counts.size())
                    != static_cast<int>(sizeof(T) * counts.size())) {
                std::cerr << "ERROR: Cannot read k-mer features" << std::endl;
                exit(1);
            }

            callback(kmer_length, read_stream, counts.data());

            if (with_reverse) {
                reverse_complement(read_stream->seq);
                std::reverse(counts.begin(), counts.end());
                callback(kmer_length, read_stream, counts.data());
            }
        },
        false
    );

    if (!gzeof(fasta_p) || gzgetc(features_p) != -1 || !gzeof(features_p)) {
        std::cerr << "ERROR: There are features left in extension unread" << std::endl;
        exit(1);
    }

    gzclose(fasta_p);
    gzclose(features_p);
}

template
void read_extended_fasta_file_critical(const std::string &filebase,
                                       const std::string &feature_name,
                                       std::function<void(size_t, const kseq_t*, const uint8_t*)> callback,
                                       bool with_reverse);
template
void read_extended_fasta_file_critical(const std::string &filebase,
                                       const std::string &feature_name,
                                       std::function<void(size_t, const kseq_t*, const uint16_t*)> callback,
                                       bool with_reverse);
template
void read_extended_fasta_file_critical(const std::string &filebase,
                                       const std::string &feature_name,
                                       std::function<void(size_t, const kseq_t*, const uint32_t*)> callback,
                                       bool with_reverse);
template
void read_extended_fasta_file_critical(const std::string &filebase,
                                       const std::string &feature_name,
                                       std::function<void(size_t, const kseq_t*, const uint64_t*)> callback,
                                       bool with_reverse);


void read_fasta_from_string(const std::string &fasta_flat,
                            std::function<void(kseq_t*)> callback,
                            bool with_reverse) {
    // create pipe
    int p[2];
    if (pipe(p) != 0) {
        std::cerr << "ERROR: opening for writing to pipe failed" << std::endl;
        exit(1);
    }

    // seems silly to spawn a thread for that, but for larger strings write would just block.
    // attempted with fmemopen as well, which, however, doesn't provide a file descriptor. But this seems necessary for gzdopen/gzopen
    std::thread writer([&]() {
        auto sent = write(p[1], fasta_flat.c_str(), fasta_flat.size());
        close(p[1]);

        if (sent != static_cast<int64_t>(fasta_flat.length())) {
            std::cerr << "ERROR: writing to pipe failed" << std::endl;

            close(p[0]);
            exit(1);
        }
    });

    // gzFile from pipe
    gzFile input_p = gzdopen(p[0], "r");
    if (input_p == Z_NULL) {
        std::cerr << "ERROR: opening for reading from pipe failed" << std::endl;
        close(p[0]);
        exit(1);
    }

    read_fasta_file_critical(input_p, callback, with_reverse);

    writer.join();

    gzclose(input_p);
    close(p[0]);
}


void read_vcf_file_critical(const std::string &filename,
                            const std::string &ref_filename,
                            size_t k,
                            std::function<void(std::string_view)> callback,
                            bool with_reverse) {
    VCFParser vcf(ref_filename, filename, k);

    if (with_reverse) {
        vcf.call_sequences(
            [&](auto&& sequence) {
                callback(sequence);
                reverse_complement(sequence.begin(), sequence.end());
                callback(sequence);
            }
        );
    } else {
        vcf.call_sequences(callback);
    }
}

void read_vcf_file_with_annotations_critical(
      const std::string &filename,
      const std::string &ref_filename,
      size_t k,
      std::function<void(std::string&&,
                         const std::vector<std::string>&)> callback,
      bool with_reverse
    ) {
    //TODO: make this a configurable option
    //default list of tokens to extract as annotations
    //TODO: extract these guys directly from vcf parsed
    const std::vector<std::string> annots = {
        "AC_AFR", "AC_EAS", "AC_AMR", "AC_ASJ",
        "AC_FIN", "AC_NFE", "AC_SAS", "AC_OTH"
    };

    VCFParser vcf(ref_filename, filename, k);

    if (with_reverse) {
        vcf.call_annotated_sequences(
            [&](auto&& sequence, const auto &annotation) {
                callback(std::string(sequence.begin(), sequence.end()), annotation);
                reverse_complement(sequence.begin(), sequence.end());
                callback(std::move(sequence), annotation);
            },
            annots
        );
    } else {
        vcf.call_annotated_sequences(callback, annots);
    }
}

} // namespace seq_io
} // namespace mtg
