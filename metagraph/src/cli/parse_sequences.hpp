#ifndef __SEQUENCE_READER_HPP__
#define __SEQUENCE_READER_HPP__

#include <string>
#include <vector>
#include <filesystem>

#include <ips4o.hpp>
#include <tsl/hopscotch_map.h>

#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/template_utils.hpp"
#include "common/utils/string_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "seq_io/formats.hpp"
#include "seq_io/kmc_parser.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"


namespace mtg {
namespace cli {

using namespace mtg::seq_io;

template <class Callback>
void parse_sequences(const std::string &file,
                     const Config &config,
                     Callback call_sequence);


template <class Callback>
void parse_sequences_and_log(const std::string &file,
                             const Config &config,
                             size_t k,
                             Callback call_sequence) {
    uint64_t num_bp = 0;
    uint64_t num_kmers = 0;
    Timer timer;

    parse_sequences(file, config,
        [&](std::string_view seq, auto count) {
            call_sequence(seq, count);
            num_bp += seq.size();
            num_kmers += seq.size() >= k ? seq.size() - k + 1 : 0;
        }
    );

    mtg::common::logger->trace("Extracted all sequences (total of {} bp, {} k-mers) from file {} in {} sec",
                               num_bp, num_kmers, file, timer.elapsed());
}


template <class Callback>
void parse_sequences(const std::string &file,
                     const Config &config,
                     Callback call_sequence) {
    mtg::common::logger->trace("Parsing {}", file);

    if (config.graph_mode == graph::DeBruijnGraph::PRIMARY && file_format(file) != "FASTA") {
        mtg::common::logger->error("Primary graphs can only be constructed from"
                                   " primary contigs");
        exit(1);
    }

    if (file_format(file) == "VCF") {
        read_vcf_file_critical(file,
                               config.refpath,
                               config.k,
                               [&](const auto &seq) { call_sequence(seq, 1); },
                               config.forward_and_reverse);

    } else if (file_format(file) == "KMC") {
        bool warning_different_k = false;

        auto min_count = config.min_count;
        auto max_count = config.max_count;

        // compute the quantiles to update the thresholds
        if (config.min_count_quantile > 0 || config.max_count_quantile < 1) {
            tsl::hopscotch_map<uint64_t, uint64_t> count_hist;
            read_kmers(file, [&](auto, uint32_t count) { count_hist[count]++; }, false);

            if (count_hist.size()) {
                std::vector<std::pair<uint64_t, uint64_t>> count_hist_v(count_hist.begin(),
                                                                        count_hist.end());

                ips4o::parallel::sort(count_hist_v.begin(), count_hist_v.end(),
                                      utils::LessFirst(), get_num_threads());

                if (config.min_count_quantile > 0)
                    min_count = utils::get_quantile(count_hist_v, config.min_count_quantile);
                if (config.max_count_quantile < 1)
                    max_count = utils::get_quantile(count_hist_v, config.max_count_quantile);

                mtg::common::logger->info("Used k-mer count thresholds:\n"
                                         "min (including): {}\n"
                                         "max (excluding): {}", min_count, max_count);
            }
        }

        read_kmers(
            file,
            [&](std::string_view sequence, uint32_t count) {
                if (!warning_different_k && sequence.size() != config.k) {
                    mtg::common::logger->warn("k-mers parsed from KMC database {} have "
                                             "length {} but graph is constructed for k={}",
                                             file, sequence.size(), config.k);
                    warning_different_k = true;
                }
                if (config.forward_and_reverse) {
                    std::string reverse(sequence);
                    reverse_complement(reverse.begin(), reverse.end());
                    call_sequence(sequence, count);
                    call_sequence(reverse, count);
                } else {
                    call_sequence(sequence, count);
                }
            },
            // For canonical graph the rev-compl sequences will be called automatically anyway.
            // Also, if forward_and_reverse = true, the rev-compl will be computed in the callback.
            config.graph_mode != graph::DeBruijnGraph::CANONICAL && !config.forward_and_reverse,
            min_count, max_count
        );

    } else if (file_format(file) == "FASTA"
                || file_format(file) == "FASTQ") {

        if (config.count_kmers
            && std::filesystem::exists(utils::remove_suffix(file, ".gz", ".fasta") + ".kmer_counts.gz")) {

            mtg::common::logger->trace("Parsing k-mer counts from {}",
                utils::remove_suffix(file, ".gz", ".fasta") + ".kmer_counts.gz"
            );
            read_extended_fasta_file_critical<uint32_t>(file, "kmer_counts",
                [&](size_t k, const kseq_t *read_stream, const uint32_t *kmer_counts) {
                    if (k != config.k) {
                        mtg::common::logger->error("File {} contains counts for k-mers of "
                                                  "length {} but graph is constructed with k={}",
                                                  file, k, config.k);
                        exit(1);
                    }
                    assert(read_stream->seq.l >= k && "sequences can't be shorter than k-mers");

                    const uint32_t *kmer_counts_end = kmer_counts + read_stream->seq.l - k + 1;
                    const uint32_t *same_counts_end;

                    const char *seq = read_stream->seq.s;

                    do {
                        same_counts_end = std::find_if(kmer_counts, kmer_counts_end,
                            [&](auto count) { return count != *kmer_counts; }
                        );
                        size_t segment_size = same_counts_end - kmer_counts;

                        call_sequence(std::string_view(seq, segment_size + k - 1), *kmer_counts);

                        kmer_counts += segment_size;
                        seq += segment_size;

                    } while (kmer_counts < kmer_counts_end);

                    assert(seq + k - 1 == read_stream->seq.s + read_stream->seq.l);
                },
                config.forward_and_reverse
            );

        } else {
            read_fasta_file_critical(file, [&](kseq_t *read_stream) {
                // add read to the graph constructor as a callback
                call_sequence(std::string_view(read_stream->seq.s, read_stream->seq.l), 1);
            }, config.forward_and_reverse);
        }
    } else {
        mtg::common::logger->error("File type unknown for '{}'", file);
        exit(1);
    }
}

} // namespace cli
} // namespace mtg

#endif // __SEQUENCE_READER_HPP__
