#include "sequence_io.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <thread>

#include <unistd.h>

#include "common/seq_tools/reverse_complement.hpp"
#include "common/utils/string_utils.hpp"
#include "vcf_parser.hpp"


namespace mtg {
namespace seq_io {

// subset of KSTREAM_INIT
__KS_BASIC(/**/, compFile, 16384)
__KS_GETUNTIL(/**/, compFile::read)
__KS_INLINED(compFile::read)

// subset of KSEQ_INIT
__KSEQ_BASIC(/**/, compFile)
__KSEQ_READ(/**/)

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
    fname_ = utils::remove_suffix(filebase, ".zst");
    bool zstd = (fname_.size() < filebase.size());
    std::string comp_type = zstd ? ".zst" : ".gz";

    fname_ = utils::remove_suffix(filebase, comp_type, ".fasta") + ".fasta" + comp_type;

    out_ = compFile::open_write(fname_.c_str());
    if (!out_.good()) {
        std::cerr << "ERROR: Can't write to " << fname_ << std::endl;
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

FastaWriter::~FastaWriter() { flush(); }

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
    if (!write_fasta(out_,
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

    fasta_fname_ = utils::remove_suffix(filebase, ".zst");
    bool zstd = (fasta_fname_.size() < filebase.size());
    std::string comp_type = zstd ? ".zst" : ".gz";

    fasta_fname_ = utils::remove_suffix(filebase, comp_type, ".fasta") + ".fasta" + comp_type;

    fasta_out_ = compFile::open_write(fasta_fname_.c_str());
    if (!fasta_out_.good()) {
        std::cerr << "ERROR: Can't write to " << fasta_fname_ << std::endl;
        exit(1);
    }

    auto filename = utils::remove_suffix(filebase, comp_type, ".fasta") + "." + feature_name + comp_type;

    feature_out_ = compFile::open_write(filename.c_str());
    if (!feature_out_.good()
            || feature_out_.write(&kmer_length_, 4) != sizeof(kmer_length_)) {
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
ExtendedFastaWriter<T>::~ExtendedFastaWriter() { flush(); }

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
    if (!write_fasta(fasta_out_,
                     enumerate_sequences_ ? header_ + std::to_string(++count_)
                                          : header_,
                     sequence)) {
        std::cerr << "ERROR: ExtendedFastaWriter::write failed. Can't dump sequence to fasta" << std::endl;
        exit(1);
    }

    if (feature_out_.write(kmer_features.data(),
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


bool write_fasta(compFile &out, const kseq_t &kseq) {
    return out.put('>') == '>'
        && out.write(kseq.name.s, kseq.name.l) == static_cast<int>(kseq.name.l)
        && (!kseq.comment.l
            || (out.put(' ') == ' '
                && out.write(kseq.comment.s, kseq.comment.l)
                                    == static_cast<int>(kseq.comment.l)))
        && out.put('\n') == '\n'
        && out.write(kseq.seq.s, kseq.seq.l) == static_cast<int>(kseq.seq.l)
        && out.put('\n') == '\n';
}

bool write_fasta(compFile &out, const std::string &header,
                                const std::string &sequence) {
    return out.put('>') == '>'
        && out.write(header.data(), header.length())
                        == static_cast<int>(header.length())
        && out.put('\n') == '\n'
        && out.write(sequence.data(), sequence.length())
                        == static_cast<int>(sequence.length())
        && out.put('\n') == '\n';
}

bool write_fastq(compFile &out, const kseq_t &kseq) {
    std::string qual(kseq.qual.s, kseq.qual.l);

    if (!kseq.qual.l && kseq.seq.l)
        qual.assign(kseq.seq.l, kDefaultFastQualityChar);

    return out.put('@') == '@'
        && out.write(kseq.name.s, kseq.name.l) == static_cast<int>(kseq.name.l)
        && (kseq.comment.l == 0
            || (out.put(' ') == ' '
                && out.write(kseq.comment.s, kseq.comment.l)
                                    == static_cast<int>(kseq.comment.l)))
        && out.put('\n') == '\n'
        && out.write(kseq.seq.s, kseq.seq.l) == static_cast<int>(kseq.seq.l)
        && out.put('\n') == '\n'
        && out.put('+') == '+'
        && out.write(kseq.name.s, kseq.name.l) == static_cast<int>(kseq.name.l)
        && (kseq.comment.l == 0
            || (out.put(' ') == ' '
                && out.write(kseq.comment.s, kseq.comment.l)
                                    == static_cast<int>(kseq.comment.l)))
        && out.put('\n') == '\n'
        && out.write(qual.data(), qual.size()) == static_cast<int>(qual.size())
        && out.put('\n') == '\n';
}

FastaParser::iterator::iterator(const std::string &filename,
                                bool with_reverse_complement)
      : filename_(filename),
        with_reverse_complement_(with_reverse_complement) {
    compFile input_p = compFile::open_read(filename_.c_str());
    if (!input_p.good()) {
        std::cerr << "ERROR: Cannot read from file " << filename_ << std::endl;
        exit(1);
    }

    read_stream_.reset(kseq_init(input_p));
    if (!read_stream_) {
        std::cerr << "ERROR: failed to initialize kseq file descriptor" << std::endl;
        exit(1);
    }

    if (kseq_read(read_stream_.get()) < 0)
        *this = iterator();
}

void FastaParser::iterator::deinit_stream::operator()(kseq_t *read_stream) {
    if (read_stream) {
        // move the input file so it can be destroyed at the end of this scope
        // by ~compFile
        compFile input_p = std::move(read_stream->f->f);

        kseq_destroy(read_stream);
    }
}


template <class Callback>
void read_fasta_file_critical(compFile &input_p,
                              Callback callback,
                              bool with_reverse) {
    if (!input_p.good()) {
        std::cerr << "ERROR: Null file descriptor" << std::endl;
        exit(1);
    }
    //TODO: handle read_stream->qual
    FastaParser::iterator::stream_type read_stream;
    read_stream.reset(kseq_init(input_p));
    if (!read_stream) {
        std::cerr << "ERROR: failed to initialize kseq file descriptor" << std::endl;
        exit(1);
    }

    while (kseq_read(read_stream.get()) >= 0) {
        callback(read_stream.get());
        if (with_reverse) {
            reverse_complement(read_stream->seq);
            callback(read_stream.get());
        }
    }
}

void read_fasta_file_critical(const std::string &filename,
                              std::function<void(kseq_t*)> callback,
                              bool with_reverse) {
    compFile input_p = compFile::open_read(filename.c_str());
    if (!input_p.good()) {
        std::cerr << "ERROR: Cannot read file " << filename << std::endl;
        exit(1);
    }

    read_fasta_file_critical(input_p, callback, with_reverse);
}

template <typename T>
void read_extended_fasta_file_critical(const std::string &filebase,
                                       const std::string &feature_name,
                                       std::function<void(size_t, const kseq_t*, const T*)> callback,
                                       bool with_reverse) {
    assert(feature_name.size());

    const char *dirpos = strrchr(filebase.c_str(), '/');
    const char *dotpos = strrchr(filebase.c_str(), '.');

    std::string ext = (dotpos == NULL || (dirpos != NULL && dirpos > dotpos))
        ? ".gz"
        : dotpos;

    auto filename = filebase;

    if (ext != ".zst" && ext != ".gz") {
        ext = "";
    } else {
        filename = utils::remove_suffix(filename, ext);
    }

    auto feature_base = utils::remove_suffix(filebase, ".fasta");
    filename = feature_base + ".fasta" + ext;

    compFile fasta_p = compFile::open_read(filename.c_str());
    if (!fasta_p.good()) {
        std::cerr << "ERROR: Cannot read from file " << filename << std::endl;
        exit(1);
    }

    filename = feature_base + "." + feature_name + ext;
    uint32_t kmer_length;

    compFile features_p = compFile::open_read(filename.c_str());
    if (!features_p.good()
            || features_p.read(&kmer_length, 4) != sizeof(kmer_length)) {
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
            if (features_p.read(counts.data(), sizeof(T) * counts.size())
                    != static_cast<int>(sizeof(T) * counts.size())) {
                std::cerr << "ERROR: Cannot read k-mer features from file " << filename << std::endl;
                std::cerr << "ERROR: " << features_p.error() << std::endl;
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

    if (!fasta_p.eof() || !features_p.eof()) {
        std::cerr << "ERROR: There are features left in extension unread" << std::endl;
        exit(1);
    }
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
            exit(1);
        }
    });

    // gzFile from pipe
    // TODO: switch to zstd?
    auto input_p = compFile::dopen_read<boost::iostreams::gzip_decompressor>(p[0]);
    if (!input_p.good()) {
        std::cerr << "ERROR: opening for reading from pipe failed" << std::endl;
        exit(1);
    }

    read_fasta_file_critical(input_p, callback, with_reverse);

    writer.join();
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
