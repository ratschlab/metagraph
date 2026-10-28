// Benchmark driver for the aligner (graph/alignment): per-read latency, a
// per-read watchdog, crash isolation, and seed fan-out statistics.
//
// usage:
//   align_bench MODE TIMEOUT_SEC SKIP align -i GRAPH [-a ANNO] [align options] READS.fa
//
// MODE:
//   default | map     align (or, with map, only map via LabeledAligner::map_batch)
//                     each read with one aligner object reused across reads
//                     (as a binding would), print one TSV line per read:
//                     name, status, latency_us, score, cigar, orientation,
//                     query_begin, query_end, spelling, labels (';'-joined)
//   fanout            for every --align-min-seed-length window S of every read,
//                     count the real k-mers that end with S (suffix resolution,
//                     what the default aligner seeds from) and that start with S
//                     (forward resolution; counted as k-mers ending with rc(S), so
//                     the graph must contain both strands, i.e. CANONICAL mode);
//                     print the two histograms
//
// TIMEOUT_SEC: a watchdog prints "name\tTIMEOUT" and exits with code 3 when one
//              read takes longer; the caller resumes with SKIP set past that read.
// SKIP:        number of leading reads to skip.

#include <atomic>
#include <chrono>
#include <iostream>
#include <map>
#include <thread>
#include <unistd.h>

#include "cli/align.hpp"
#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "cli/load/load_graph.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/alignment/aligner_labeled.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "seq_io/sequence_io.hpp"

using namespace mtg;
using namespace mtg::graph;

namespace {

using Clock = std::chrono::steady_clock;

std::atomic<int64_t> read_start_us { -1 };
std::string current_read;
std::mutex current_read_mutex;

int64_t now_us() {
    return std::chrono::duration_cast<std::chrono::microseconds>(
            Clock::now().time_since_epoch()).count();
}

void watchdog(double timeout_sec) {
    const int64_t timeout_us = timeout_sec * 1e6;
    while (true) {
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        int64_t start = read_start_us.load();
        if (start >= 0 && now_us() - start > timeout_us) {
            std::lock_guard<std::mutex> lock(current_read_mutex);
            std::cout << current_read << "\tTIMEOUT\t" << (now_us() - start)
                      << "\t*\t*\t*\t*\t*\t*\t*" << std::endl;
            _exit(3);
        }
    }
}

void print_histogram(const std::string &name, const std::map<size_t, size_t> &hist) {
    for (const auto &[value, count] : hist) {
        std::cout << name << '\t' << value << '\t' << count << '\n';
    }
}

} // namespace

int main(int argc, char *argv[]) {
    if (argc < 6) {
        std::cerr << "usage: " << argv[0]
                  << " default|map|fanout TIMEOUT_SEC SKIP align [options] READS.fa\n";
        return 1;
    }
    const std::string mode = argv[1];
    const double timeout_sec = std::stod(argv[2]);
    const size_t skip = std::stoull(argv[3]);

    std::vector<char*> cli_argv { argv[0] };
    cli_argv.insert(cli_argv.end(), argv + 4, argv + argc);
    cli::Config config(cli_argv.size(), cli_argv.data());

    std::shared_ptr<DeBruijnGraph> graph = cli::load_critical_dbg(config.infbase);

    if (mode == "fanout") {
        const auto &dbg = dynamic_cast<const DBGSuccinct&>(*graph);
        if (dbg.get_mode() != DeBruijnGraph::CANONICAL) {
            std::cerr << "fanout needs a CANONICAL graph\n";
            return 1;
        }
        const size_t s = config.alignment_min_seed_length;
        std::map<size_t, size_t> suffix_hist, prefix_hist;
        size_t num_reads = 0;
        for (const auto &file : config.fnames) {
            seq_io::read_fasta_file_critical(file, [&](seq_io::kseq_t *read) {
                ++num_reads;
                std::string_view seq(read->seq.s, read->seq.l);
                std::string rc;
                for (size_t i = 0; i + s <= seq.size(); ++i) {
                    std::string_view window = seq.substr(i, s);
                    size_t num_suffix = 0;
                    dbg.call_nodes_with_suffix_matching_longest_prefix(
                        window, [&](auto, auto) { ++num_suffix; }, s);
                    rc.assign(window);
                    reverse_complement(rc.begin(), rc.end());
                    size_t num_prefix = 0;
                    dbg.call_nodes_with_suffix_matching_longest_prefix(
                        rc, [&](auto, auto) { ++num_prefix; }, s);
                    if (num_suffix || num_prefix) {
                        ++suffix_hist[num_suffix];
                        ++prefix_hist[num_prefix];
                    } else {
                        ++suffix_hist[0];
                        ++prefix_hist[0];
                    }
                }
            });
        }
        std::cout << "reads\t" << num_reads << '\n';
        print_histogram("suffix", suffix_hist);
        print_histogram("prefix", prefix_hist);
        return 0;
    }

    if (mode != "default" && mode != "map") {
        std::cerr << "unknown mode " << mode << '\n';
        return 1;
    }

    std::unique_ptr<AnnotatedDBG> anno;
    if (!config.infbase_annotators.empty())
        anno = cli::initialize_annotated_dbg(graph, config);

    // as in cli::align_to_graph
    auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(graph);
    if (dbg_succ)
        dbg_succ->reset_mask();
    std::shared_ptr<DeBruijnGraph> aln_graph = graph;
    if (aln_graph->get_mode() == DeBruijnGraph::PRIMARY)
        aln_graph = std::make_shared<CanonicalDBG>(aln_graph);
    if (dbg_succ)
        aln_graph->add_extension(std::make_shared<NodeFirstCache>(*dbg_succ));

    const auto old_config = cli::initialize_aligner_config(config, *aln_graph);
    std::unique_ptr<align::IDBGAligner> old_aligner;
    if (anno) {
        old_aligner = std::make_unique<align::LabeledAligner<>>(
                *aln_graph, old_config, anno->get_annotator());
    } else {
        old_aligner = std::make_unique<align::DBGAligner<>>(*aln_graph, old_config);
    }

    std::thread(watchdog, timeout_sec).detach();

    size_t read_id = 0;
    for (const auto &file : config.fnames) {
        seq_io::read_fasta_file_critical(file, [&](seq_io::kseq_t *read) {
            if (read_id++ < skip)
                return;

            const std::string seq(read->seq.s, read->seq.l);
            {
                std::lock_guard<std::mutex> lock(current_read_mutex);
                current_read = read->name.s;
            }
            std::string line;
            int64_t start = now_us();
            read_start_us = start;
            try {
                if (mode == "map") {
                    std::optional<align::LabeledAligner<>::LabelMapping> mapping;
                    dynamic_cast<const align::LabeledAligner<>&>(*old_aligner).map_batch(
                        { { std::string(read->name.s), seq } },
                        [&](const std::string &, auto&& result) { mapping = std::move(result); }
                    );
                    int64_t latency = now_us() - start;
                    if (mapping) {
                        std::vector<std::string> labels;
                        for (auto column : mapping->labels) {
                            labels.emplace_back(
                                anno->get_annotator().get_label_encoder().decode(column));
                        }
                        line = fmt::format("\tOK\t{}\t{}\t*\t{}\t{}\t{}\t*\t{}", latency,
                                           mapping->num_matches, mapping->orientation,
                                           mapping->query_begin, mapping->query_end,
                                           fmt::join(labels, ";"));
                    } else {
                        line = fmt::format("\tUNMAPPED\t{}\t*\t*\t*\t*\t*\t*\t*", latency);
                    }
                } else {
                    auto results = old_aligner->align(seq);
                    int64_t latency = now_us() - start;
                    const align::Alignment *best = nullptr;
                    for (const auto &aln : results) {
                        if (!best || aln.get_score() > best->get_score())
                            best = &aln;
                    }
                    if (best) {
                        std::vector<std::string> labels;
                        for (auto column : best->label_columns) {
                            labels.emplace_back(
                                anno->get_annotator().get_label_encoder().decode(column));
                        }
                        line = fmt::format("\tOK\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", latency,
                                           best->get_score(), best->get_cigar().to_string(),
                                           best->get_orientation(), best->get_clipping(),
                                           seq.size() - best->get_end_clipping(),
                                           best->get_sequence(), fmt::join(labels, ";"));
                    } else {
                        line = fmt::format("\tUNMAPPED\t{}\t*\t*\t*\t*\t*\t*\t*", latency);
                    }
                }
            } catch (const std::exception &e) {
                line = fmt::format("\tERROR\t{}\t*\t*\t*\t*\t*\t*\t{}", now_us() - start,
                                   e.what());
            }
            read_start_us = -1;
            std::lock_guard<std::mutex> lock(current_read_mutex);
            // flush per read, so that a crash pins down the read that caused it
            std::cout << read->name.s << line << std::endl;
        });
    }
    std::cout.flush();
    return 0;
}
