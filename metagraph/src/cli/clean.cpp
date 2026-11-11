#include "clean.hpp"

#include <ips4o.hpp>

#include "annotation/annotation_remap.hpp"
#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/graph_extensions/node_weights.hpp"
#include "graph/graph_cleaning.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using mtg::seq_io::FastaWriter;
using mtg::seq_io::ExtendedFastaWriter;


int clean_graph(Config *config) {
    assert(config);
    assert(config->fnames.size() == 1);

    if (!config->to_fasta) {
        logger->error("Clean graph can be serialized only in form "
                      "of its contigs or unitigs. Add flag --to-fasta");
        exit(1);
    }

    logger->trace("Loading graph...");

    return clean_graph(config, load_critical_dbg(config->fnames.at(0)));
}

int clean_graph(Config *config, std::shared_ptr<graph::DeBruijnGraph> graph,
                std::pair<sdsl::bit_vector, sdsl::bit_vector> *succ_pred) {
    assert(config);

    const auto &files = config->fnames;

    assert(files.size() == 1);
    assert(config->outfbase.size());
    assert(config->to_fasta);

    config->min_count = std::max(1u, config->min_count);

    Timer timer;

    if (config->row_map && graph->get_mode() == graph::DeBruijnGraph::CANONICAL) {
        logger->error("--row-map supports only graphs in basic or primary mode");
        exit(1);
    }
    // try loading k-mer counts
    auto node_weights = graph->load_extension<graph::NodeWeights>(files.at(0));

    if (node_weights) {
        if (auto *dbg_succ = dynamic_cast<graph::DBGSuccinct*>(graph.get()))
            dbg_succ->reset_mask();

        if (!node_weights->is_compatible(*graph)) {
            logger->error("Loaded k-mer counts not compatible with graph '{}'", files.at(0));
            exit(1);
        }
    }

    if (graph->get_mode() == graph::DeBruijnGraph::PRIMARY) {
        logger->warn("Cleaning primary graphs is not recommended."
                     " Consider building a canonical graph for cleaning instead.");
    }

    if (config->min_count > 1
            || config->max_count < std::numeric_limits<decltype(config->max_count)>::max()
            || config->min_unitig_median_kmer_abundance != 1
            || config->count_slice_quantiles[0] != 0
            || config->count_slice_quantiles[1] != 1) {
        if (!node_weights) {
            logger->error("Cannot load k-mer counts from file '{}'", files.at(0));
            exit(1);
        }

        if (config->min_unitig_median_kmer_abundance == 0) {
            // skip zero k-mer counts for dummy k-mers in DBGSuccinct
            const auto _graph = dynamic_cast<graph::DBGSuccinct*>(graph.get())
                    ? std::make_shared<graph::MaskedDeBruijnGraph>(graph,
                        [&](auto i) { return (*node_weights)[i] > 0; }, true)
                    : graph;

            uint64_t cutoff
                = estimate_min_kmer_abundance(*_graph, *node_weights,
                                              config->num_singleton_kmers,
                                              config->cleaning_threshold_percentile);

            if (cutoff != static_cast<uint64_t>(-1)) {
                config->min_unitig_median_kmer_abundance = cutoff;
            } else {
                if (config->fallback_abundance_cutoff == -1) {
                    logger->error("Cannot estimate expected minimum k-mer abundance"
                                  " and fallback is disabled (--fallback -1). Terminating.");
                    std::exit(129);
                }
                logger->warn("Cannot estimate expected minimum k-mer abundance."
                             " Using fallback value: {}", config->fallback_abundance_cutoff);
                config->min_unitig_median_kmer_abundance = config->fallback_abundance_cutoff;
            }
        }

        if (config->min_count > 1
                || config->max_count < std::numeric_limits<decltype(config->max_count)>::max()) {
            const auto &weights = *graph->get_extension<graph::NodeWeights>();

            graph = std::make_shared<graph::MaskedDeBruijnGraph>(graph,
                [&](auto i) { return weights[i] >= config->min_count
                                    && weights[i] <= config->max_count; },
                true,
                graph->get_mode()
            );
            graph->add_extension(node_weights);

            assert(node_weights->is_compatible(*graph));
        }
    }

    logger->trace("Graph loaded in {} sec", timer.elapsed());

    if (dynamic_cast<const graph::MaskedDeBruijnGraph *>(graph.get())) {
        logger->trace("Extracting sequences from subgraph...");
    } else {
        logger->trace("Extracting sequences from graph...");
    }

    timer.reset();

    auto call_clean_contigs = [&](auto callback, size_t num_threads) {
        if (config->min_unitig_median_kmer_abundance != 1) {
            assert(node_weights);
            if (!node_weights->is_compatible(*graph)) {
                logger->error("k-mer counts are not compatible with the subgraph");
                exit(1);
            }

            logger->info("Threshold for median k-mer abundance in unitigs: {}",
                         config->min_unitig_median_kmer_abundance);

            graph->call_unitigs([&](const std::string &unitig, const auto &path) {
                if (!is_unreliable_unitig(path,
                                          *node_weights,
                                          config->min_unitig_median_kmer_abundance))
                    callback(unitig, path);
            }, num_threads, config->min_tip_size, graph->get_mode() == graph::DeBruijnGraph::CANONICAL);

        } else if (config->unitigs || config->min_tip_size > 1 || config->smoothing_window > 1) {
            graph->call_unitigs(callback,
                                num_threads,
                                config->min_tip_size,
                                graph->get_mode() == graph::DeBruijnGraph::CANONICAL);

        } else {
            graph->call_sequences(callback, num_threads,
                                  graph->get_mode() == graph::DeBruijnGraph::CANONICAL);
        }
    };

    auto dump_contigs_to_fasta = [&](const std::string &outfbase, auto call_contigs) {
        std::mutex seq_mutex;
        uint64_t num_bp = 0;
        uint64_t num_kmers = 0;
        std::string fasta_fname;

        // An SSHash graph built from the written file numbers its k-mers in file
        // order, so the row of every k-mer in that graph is known here.
        std::unique_ptr<annot::RowMapWriter> row_map;
        std::filesystem::path row_map_tmp_dir;
        if (config->row_map) {
            row_map_tmp_dir = utils::create_temp_dir(
                    config->tmp_dir.empty() ? std::filesystem::path(outfbase).remove_filename()
                                            : config->tmp_dir,
                    "rowmap");
            // in primary mode, an SSHash graph has a second node for every reverse complement
            row_map = std::make_unique<annot::RowMapWriter>(graph->max_index(),
                                                            2 * graph->max_index(),
                                                            row_map_tmp_dir);
        }
        const graph::boss::BOSS *boss = nullptr;
        const bit_vector *boss_mask = nullptr;
        // The degrees in the source graph are those in the SSHash graph only if all
        // k-mers of a basic graph are written: no subgraph, no k-mer counts to filter
        // by, no pruned tips.
        auto *dbg_succ = dynamic_cast<const graph::DBGSuccinct*>(graph.get());
        if (!dbg_succ || node_weights || config->min_tip_size > 1
                || graph->get_mode() != graph::DeBruijnGraph::BASIC)
            succ_pred = nullptr;
        if (succ_pred) {
            boss = &dbg_succ->get_boss();
            boss_mask = dbg_succ->get_mask();
            // resized to the number of written k-mers after the extraction
            succ_pred->first = sdsl::bit_vector(graph->max_index() + 1, false);
            succ_pred->second = sdsl::bit_vector(graph->max_index() + 1, false);
        }
        // `first` is the number of k-mers written before this contig
        auto add_rows = [&](const auto &path, uint64_t first) {
            if (row_map) {
                uint64_t target = first;
                for (auto node : path) {
                    row_map->add(node - 1, ++target);
                }
            }
            if (succ_pred) {
                // K-mer j of this contig is SSHash node first + 1 + j. Edge path[j + 1]
                // leaves the BOSS node that edge path[j] enters, so the degrees of
                // that node are read here without fwd or bwd.
                for (size_t j = 0; j + 1 < path.size(); ++j) {
                    if (boss->is_single_outgoing(path[j + 1]))
                        set_bit(succ_pred->first.data(), first + 1 + j, true);

                    auto w = boss->get_W(path[j]);
                    bool single_incoming;
                    if (w < boss->alph_size) {
                        // path[j] is the first edge incoming to that node
                        single_incoming = boss->is_single_incoming(path[j], w);
                    } else if (!boss_mask) {
                        single_incoming = false;
                    } else {
                        // single if the only other incoming edge is a masked dummy
                        auto d = w - boss->alph_size;
                        auto x = boss->pred_W(path[j], d);
                        single_incoming = !(*boss_mask)[x]
                                    && boss->num_incoming_to_target(x, d) == 2;
                    }
                    if (single_incoming)
                        set_bit(succ_pred->second.data(), first + 2 + j, true);
                }
            }
        };
        // the end of the last written contig, for the edge between two contigs
        std::string prev_suffix;
        bool prev_single_outgoing = false;
        // call under the lock, before num_kmers is advanced
        auto link_contigs = [&](const std::string &contig, bool single_incoming,
                                bool single_outgoing) {
            if (!succ_pred)
                return;
            const size_t overlap = graph->get_k() - 1;
            if (num_kmers && prev_suffix == std::string_view(contig).substr(0, overlap)) {
                if (prev_single_outgoing)
                    set_bit(succ_pred->first.data(), num_kmers, true);
                if (single_incoming)
                    set_bit(succ_pred->second.data(), num_kmers + 1, true);
            }
            prev_suffix = contig.substr(contig.size() - overlap);
            prev_single_outgoing = single_outgoing;
        };

        if (node_weights) {
            if (!node_weights->is_compatible(*graph)) {
                logger->error("k-mer counts are not compatible with the subgraph");
                exit(1);
            }

            ExtendedFastaWriter<uint32_t> writer(outfbase,
                                                 "kmer_counts",
                                                 graph->get_k(),
                                                 config->header,
                                                 config->enumerate_out_sequences,
                                                 get_num_threads() > 1);
            call_contigs([&](const std::string &contig, const auto &path) {
                std::vector<uint32_t> kmer_counts;
                kmer_counts.reserve(path.size());
                for (auto node : path) {
                    kmer_counts.push_back((*node_weights)[node]);
                }
                // smooth k-mer counts in the unitig
                utils::smooth_vector(config->smoothing_window, &kmer_counts);

                uint64_t first;
                {
                    std::lock_guard<std::mutex> lock(seq_mutex);
                    writer.write(contig, kmer_counts);
                    num_bp += contig.size();
                    first = num_kmers;
                    num_kmers += contig.size() - graph->get_k() + 1;
                }
                add_rows(path, first);
            }, get_num_threads());

            fasta_fname = writer.get_fasta_fname();

        } else {
            FastaWriter writer(outfbase, config->header,
                               config->enumerate_out_sequences,
                               get_num_threads() > 1);

            call_contigs([&](const std::string &contig, const auto &path) {
                uint64_t first;
                const bool single_incoming = succ_pred && graph->has_single_incoming(path.front());
                const bool single_outgoing = succ_pred && graph->has_single_outgoing(path.back());
                {
                    std::lock_guard<std::mutex> lock(seq_mutex);
                    writer.write(contig);
                    link_contigs(contig, single_incoming, single_outgoing);
                    num_bp += contig.size();
                    first = num_kmers;
                    num_kmers += contig.size() - graph->get_k() + 1;
                }
                add_rows(path, first);
            }, get_num_threads());

            fasta_fname = writer.get_fname();
        }
        logger->trace("Extracted {} bp of sequences with {} k-mers and wrote to {}",
                      num_bp, num_kmers, fasta_fname);

        if (succ_pred) {
            succ_pred->first.resize(num_kmers + 1);
            succ_pred->second.resize(num_kmers + 1);
        }

        if (row_map) {
            const std::string row_map_fname = outfbase + ".rowmap";
            try {
                row_map->write(row_map_fname,
                               graph->get_mode() == graph::DeBruijnGraph::PRIMARY
                                    ? 2 * num_kmers : num_kmers,
                               get_num_threads());
            } catch (const std::exception &e) {
                logger->error("{}", e.what());
                row_map.reset();
                utils::remove_temp_dir(row_map_tmp_dir);
                std::error_code ec;
                std::filesystem::remove(row_map_fname, ec);
                exit(1);
            }
            row_map.reset();
            utils::remove_temp_dir(row_map_tmp_dir);
            logger->trace("Row map for an SSHash graph built from {} saved to {}",
                          fasta_fname, row_map_fname);
        }
    };

    assert(config->count_slice_quantiles.size() >= 2);

    if (config->count_slice_quantiles[0] == 0
            && config->count_slice_quantiles[1] == 1) {
        dump_contigs_to_fasta(config->outfbase, call_clean_contigs);

    } else {
        if (!node_weights) {
            logger->error("Need k-mer counts for binning k-mers by abundance");
            exit(1);
        }
        assert(node_weights->is_compatible(*graph));

        auto &weights = node_weights->get_data();

        assert(graph->max_index() + 1 == weights.size());

        // compute clean count histogram
        tsl::hopscotch_map<uint64_t, uint64_t> count_hist;

        if (config->min_unitig_median_kmer_abundance != 1 || config->min_tip_size > 1) {
            // cleaning required
            sdsl::bit_vector removed_nodes(weights.size(), 1);

            call_clean_contigs([&](const std::string&, const auto &path) {
                for (auto i : path) {
                    assert(weights[i]);
                    count_hist[weights[i]]++;
                    removed_nodes[i] = 0;
                }
            }, get_num_threads());

            call_ones(removed_nodes, [&weights](auto i) { weights[i] = 0; });

        } else if (auto dbg_succ = std::dynamic_pointer_cast<graph::DBGSuccinct>(graph)) {
            // use entire graph without dummy BOSS edges
            graph->call_nodes([&](auto i) {
                if (uint64_t count = weights[i])
                    count_hist[count]++;
            });
        } else {
            // use entire graph
            graph->call_nodes([&](auto i) {
                assert(weights[i]);
                count_hist[weights[i]]++;
            });
        }
        // must not have any zero weights
        assert(!count_hist.count(0));

        std::vector<std::pair<uint64_t, uint64_t>> count_hist_v(count_hist.begin(),
                                                                count_hist.end());
        count_hist.clear();

        ips4o::parallel::sort(count_hist_v.begin(), count_hist_v.end(),
                              utils::LessFirst(), get_num_threads());

        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 1; i < config->count_slice_quantiles.size(); ++i) {
            // extract sequences for k-mer counts bin |i|
            assert(config->count_slice_quantiles[i - 1] < config->count_slice_quantiles[i]);

            auto filebase = utils::remove_suffix(config->outfbase, ".gz", ".fasta")
                    + fmt::format(".{}.{}", config->count_slice_quantiles[i - 1],
                                            config->count_slice_quantiles[i]);

            if (!count_hist_v.size()) {
                dump_contigs_to_fasta(filebase, [](auto, auto) {});
                continue;
            }

            uint64_t min_count = config->count_slice_quantiles[i - 1] > 0
                ? utils::get_quantile(count_hist_v, config->count_slice_quantiles[i - 1])
                : 1;
            uint64_t max_count = config->count_slice_quantiles[i] < 1
                ? utils::get_quantile(count_hist_v, config->count_slice_quantiles[i])
                : std::numeric_limits<uint64_t>::max();

            logger->info("k-mer count thresholds:\n"
                         "min (including): {}\n"
                         "max (excluding): {}", min_count, max_count);

            assert(node_weights->is_compatible(*graph));

            graph::MaskedDeBruijnGraph graph_slice(graph,
                [&](auto i) { return weights[i] >= min_count && weights[i] < max_count; },
                true,
                graph->get_mode()
            );

            // the outer for loop is parallelized, so don't start another thread
            // pool here
            dump_contigs_to_fasta(filebase, [&](auto dump_sequence, auto /* num_threads */) {
                graph_slice.call_sequences(dump_sequence, 1,
                                           graph_slice.get_mode() == graph::DeBruijnGraph::CANONICAL);
            });
        }
    }

    logger->trace("Graph cleaning finished in {} sec", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg
