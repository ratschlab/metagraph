#ifndef __SEQUENCE_READER_HPP__
#define __SEQUENCE_READER_HPP__

#include <string>
#include <vector>

#include <ips4o.hpp>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/template_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "seq_io/formats.hpp"
#include "seq_io/kmc_parser.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"


template <class Callback, class CountedKmer, class Loop>
void parse_sequences(const std::vector<std::string> &files,
                     const Config &config,
                     const Timer &timer,
                     Callback call_sequence,
                     CountedKmer call_kmer,
                     Loop call_sequences) {
    // iterate over input files
    for (const auto &file : files) {
        mg::common::logger->trace("Parsing '{}'", file);

        Timer data_reading_timer;

        if (file_format(file) == "VCF") {
            read_vcf_file_critical(file,
                                   config.refpath,
                                   config.k,
                                   [&](std::string&& sequence) {
                                       call_sequence(std::move(sequence));
                                   },
                                   config.forward_and_reverse);

        } else if (file_format(file) == "KMC") {
            bool warning_different_k = false;

            auto min_count = config.min_count;
            auto max_count = config.max_count;

            if (config.min_count_quantile > 0 || config.max_count_quantile < 1) {
                std::unordered_map<uint64_t, uint64_t> count_hist;
                kmc::read_kmers(
                    file,
                    [&](std::string&&, uint32_t count) {
                        count_hist[count] += (1 + config.forward_and_reverse);
                    },
                    !config.canonical && !config.forward_and_reverse
                );

                if (count_hist.size()) {
                    std::vector<std::pair<uint64_t, uint64_t>> count_hist_v(count_hist.begin(),
                                                                            count_hist.end());

                    ips4o::parallel::sort(count_hist_v.begin(), count_hist_v.end(),
                                          utils::LessFirst(), get_num_threads());

                    if (config.min_count_quantile > 0)
                        min_count = utils::get_quantile(count_hist_v, config.min_count_quantile);
                    if (config.max_count_quantile < 1)
                        max_count = utils::get_quantile(count_hist_v, config.max_count_quantile);

                    mg::common::logger->info("Used k-mer count thresholds:\n"
                                             "min (including): {}\n"
                                             "max (excluding): {}", min_count, max_count);
                }
            }

            kmc::read_kmers(
                file,
                [&](std::string&& sequence, uint32_t count) {
                    if (!warning_different_k && sequence.size() != config.k) {
                            mg::common::logger->warn("k-mers parsed from KMC database '{}' have "
                                                     "length {} but graph is constructed for k={}",
                                                     file, sequence.size(), config.k);
                            warning_different_k = true;
                        }
                        if (config.forward_and_reverse) {
                            std::string reverse = sequence;
                            reverse_complement(reverse.begin(), reverse.end());
                            call_kmer(std::move(sequence), count);
                            call_kmer(std::move(reverse), count);
                        } else {
                            call_kmer(std::move(sequence), count);
                        }
                    },
                    !config.canonical && !config.forward_and_reverse, min_count, max_count);

        } else if (file_format(file) == "FASTA"
                    || file_format(file) == "FASTQ") {
            if (files.size() >= get_num_threads()) {
                auto forward_and_reverse = config.forward_and_reverse;

                // capture all required values by copying to be able
                // to run task from other threads
                call_sequences([=](auto callback) {
                    read_fasta_file_critical(file, [=](kseq_t *read_stream) {
                        // add read to the graph constructor as a callback
                        callback(read_stream->seq.s);
                    }, forward_and_reverse);
                });
            } else {
                read_fasta_file_critical(file, [&](kseq_t *read_stream) {
                    // add read to the graph constructor as a callback
                    call_sequence(read_stream->seq.s);
                }, config.forward_and_reverse);
            }
        } else {
            mg::common::logger->error("File type unknown for '{}'", file);
            exit(1);
        }

        mg::common::logger->trace("Extracted all sequences from file '{}' in {} sec", file,
                                  timer.elapsed());
    }
}

#endif // __SEQUENCE_READER_HPP__
