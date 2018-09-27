#include "sequence_io.hpp"

#include <iostream>
#include <fstream>

#include "vcf_parser.hpp"
#include "serialization.hpp"

const char kDefaultFastQualityChar = 'I';


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


void read_fasta_file_critical(const std::string &filename,
                              std::function<void(kseq_t*)> callback,
                              bool with_reverse,
                              Timer *timer,
                              const std::string &filter_filename) {
    std::vector<bool> filter;
    if (filter_filename.size()) {
        std::ifstream instream(filter_filename);
        try {
            filter = load_number_vector<bool>(instream);
        } catch (...) {
            std::cerr << "ERROR: Filter file " << filter_filename
                      << " is corrupted" << std::endl;
            exit(1);
        }
    }
    size_t seq_count = 0;

    gzFile input_p = gzopen(filename.c_str(), "r");
    if (input_p == Z_NULL) {
        std::cerr << "ERROR no such file " << filename << std::endl;
        exit(1);
    }
    //TODO: handle read_stream->qual
    kseq_t *read_stream = kseq_init(input_p);
    if (read_stream == NULL) {
        std::cerr << "ERROR while opening input file " << filename << std::endl;
        exit(1);
    }
    if (timer) {
        std::cout << "Start extracting sequences from file " << filename << std::endl;
    }
    while (kseq_read(read_stream) >= 0) {
        if (filter_filename.size() && filter.size() <= seq_count) {
            std::cerr << "ERROR: Filter file " << filter_filename
                      << " has fewer sequences" << std::endl;
            exit(1);
        }

        if (!filter_filename.size() || filter[seq_count]) {
            callback(read_stream);
            if (with_reverse) {
                reverse_complement(read_stream->seq);
                callback(read_stream);
            }
        }

        seq_count++;
    }
    if (filter_filename.size() && filter.size() != seq_count) {
        std::cerr << "ERROR: Filter file " << filter_filename
                  << " has more sequences" << std::endl;
        exit(1);
    }
    if (timer) {
        std::cout << "Finished extracting sequences from file " << filename
                  << " in " << timer->elapsed() << "sec"
                  << ", sequences extracted: " << seq_count << std::endl;
    }
    kseq_destroy(read_stream);
    gzclose(input_p);
}

void read_vcf_file_critical(const std::string &filename,
                            const std::string &ref_filename,
                            size_t k,
                            std::vector<std::string> *annotation,
                            std::function<void(std::string &,
                                               std::vector<std::string> *)> callback,
                            Timer *timer) {
    //TODO: make this a configurable option
    //default list of tokens to extract as annotations
    //TODO: extract these guys directly from vcf parsed
    const std::vector<std::string> annots = {
      "AC_AFR", "AC_EAS", "AC_AMR", "AC_ASJ",
      "AC_FIN", "AC_NFE", "AC_SAS", "AC_OTH"
    };

    vcf_parser vcf;
    if (!vcf.init(ref_filename, filename, k)) {
        std::cerr << "ERROR reading VCF " << filename << std::endl;
        exit(1);
    }
    if (timer) {
        std::cout << "Extracting sequences from file " << filename << std::endl;
    }
    size_t seq_count = 0;
    while (vcf.get_seq(annots, annotation)) {
        callback(vcf.seq, annotation);
        seq_count++;
    }
    if (timer) {
        std::cout << "Finished extracting sequences from file " << filename
                  << " in " << timer->elapsed() << "sec"
                  << ", sequences extracted: " << seq_count << std::endl;
    }
}
