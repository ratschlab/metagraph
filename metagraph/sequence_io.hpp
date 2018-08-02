#ifndef __SEQUENCE_IO__
#define __SEQUENCE_IO__

#include <functional>
#include <vector>

#include <zlib.h>
#include <htslib/kseq.h>

#include "unix_tools.hpp"
#include "helpers.hpp"

KSEQ_INIT(gzFile, gzread);


// TODO: use BGZF from htslib
bool write_fasta(gzFile gz_out, const kseq_t &kseq);

bool write_fastq(gzFile gz_out, const kseq_t &kseq);


void read_fasta_file_critical(const std::string &filename,
                              std::function<void(kseq_t*)> callback,
                              bool with_reverse = false,
                              Timer *timer = NULL,
                              const std::string &filter_filename = "");

void read_vcf_file_critical(const std::string &filename,
                            const std::string &ref_filename,
                            size_t k,
                            std::vector<std::string> *annotation,
                            std::function<void(std::string&,
                                               std::vector<std::string>*)> callback,
                            Timer *timer = NULL);

#endif // __SEQUENCE_IO__
