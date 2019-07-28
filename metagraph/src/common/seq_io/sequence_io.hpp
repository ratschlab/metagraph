#ifndef __SEQUENCE_IO__
#define __SEQUENCE_IO__

#include <functional>
#include <vector>
#include <string>

#include <zlib.h>
#include <htslib/kseq.h>

#include "unix_tools.hpp"

KSEQ_INIT(gzFile, gzread);


bool write_fasta(gzFile gz_out, const kseq_t &kseq);

bool write_fasta(gzFile gz_out, const std::string &header,
                                const std::string &sequence);

bool write_fastq(gzFile gz_out, const kseq_t &kseq);

void read_fasta_file_critical(const std::string &filename,
                              std::function<void(kseq_t*)> callback,
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
