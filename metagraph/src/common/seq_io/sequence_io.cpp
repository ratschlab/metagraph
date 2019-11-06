#include "sequence_io.hpp"

#include "vcf_parser.hpp"

#include "../seq_tools/reverse_complement.hpp"
#include "../serialization.hpp"

#include <fstream>
#include <iostream>

const char kDefaultFastQualityChar = 'I';


FastaWriter::FastaWriter(const std::string &filename,
                         const std::string &header,
                         bool write_counts) : header_(header),
                                              write_counts_(write_counts) {
    gz_out_ = gzopen(filename.c_str(), "w");
    if (gz_out_ == Z_NULL) {
        std::cerr << "ERROR: Can't write to " << filename << std::endl;
        exit(1);
    }
}

FastaWriter::~FastaWriter() {
    gzclose(gz_out_);
}

void FastaWriter::write(const std::string &sequence) {
    if (!write_fasta(gz_out_,
                     write_counts_ ? header_ + std::to_string(++count_)
                                   : header_,
                     sequence)) {
        std::cerr << "ERROR: FastaWriter::write failed. Can't dump sequence to fasta" << std::endl;
        exit(1);
    }
}


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


void read_fasta_file_critical(gzFile input_p,
                              std::function<void(kseq_t*)> callback,
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

void read_fasta_from_string(const std::string &fasta_flat,
                            std::function<void(kseq_t*)> callback,
                            bool with_reverse) {
    // create pipe
    int p[2];
    if (pipe(p) != 0) {
        std::cerr << "ERROR: opening for writing to pipe failed" << std::endl;
        exit(1);
    }
    auto sent = write(p[1], fasta_flat.c_str(), fasta_flat.length());
    close(p[1]);

    if (sent != static_cast<int64_t>(fasta_flat.length())) {
        std::cerr << "ERROR: writing to pipe failed" << std::endl;
        close(p[0]);
        exit(1);
    }

    // gzFile from pipe
    gzFile input_p = gzdopen(p[0], "r");
    if (input_p == Z_NULL) {
        std::cerr << "ERROR: opening for reading from pipe failed" << std::endl;
        close(p[0]);
        exit(1);
    }

    read_fasta_file_critical(input_p, callback, with_reverse);

    gzclose(input_p);
    close(p[0]);
}


void read_vcf_file_critical(const std::string &filename,
                            const std::string &ref_filename,
                            size_t k,
                            std::function<void(std::string&&)> callback,
                            bool with_reverse) {
    VCFParser vcf(ref_filename, filename, k);

    if (with_reverse) {
        vcf.call_sequences(
            [&](auto&& sequence) {
                callback(std::string(sequence.begin(), sequence.end()));
                reverse_complement(sequence.begin(), sequence.end());
                callback(std::move(sequence));
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
