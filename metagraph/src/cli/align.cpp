#include "align.hpp"

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/alignment/aligner_methods.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;
using namespace mtg::graph::align;

using mtg::seq_io::kseq_t;
using mtg::common::logger;


DBGAlignerConfig initialize_aligner_config(const DeBruijnGraph &graph, const Config &config) {
    DBGAlignerConfig aligner_config;

    aligner_config.queue_size = config.alignment_queue_size;
    aligner_config.bandwidth = config.alignment_vertical_bandwidth;
    aligner_config.num_alternative_paths = config.alignment_num_alternative_paths;
    aligner_config.min_seed_length = config.alignment_min_seed_length;
    aligner_config.max_seed_length = config.alignment_max_seed_length;
    aligner_config.max_num_seeds_per_locus = config.alignment_max_num_seeds_per_locus;
    aligner_config.max_nodes_per_seq_char = config.alignment_max_nodes_per_seq_char;
    aligner_config.min_cell_score = config.alignment_min_cell_score;
    aligner_config.min_path_score = config.alignment_min_path_score;
    aligner_config.xdrop = config.alignment_xdrop;
    aligner_config.exact_kmer_match_fraction = config.discovery_fraction;
    aligner_config.gap_opening_penalty = -config.alignment_gap_opening_penalty;
    aligner_config.gap_extension_penalty = -config.alignment_gap_extension_penalty;
    aligner_config.forward_and_reverse_complement = config.align_both_strands;
    aligner_config.alignment_edit_distance = config.alignment_edit_distance;
    aligner_config.alignment_match_score = config.alignment_match_score;
    aligner_config.alignment_mm_transition_score = config.alignment_mm_transition_score;
    aligner_config.alignment_mm_transversion_score = config.alignment_mm_transversion_score;

    if (!aligner_config.min_seed_length)
        aligner_config.min_seed_length = graph.get_k();

    if (!aligner_config.max_seed_length)
        aligner_config.max_seed_length = graph.get_k();

    logger->trace("Alignment settings:");
    logger->trace("\t Alignments to report: {}", aligner_config.num_alternative_paths);
    logger->trace("\t Priority queue size: {}", aligner_config.queue_size);
    logger->trace("\t Min seed length: {}", aligner_config.min_seed_length);
    logger->trace("\t Max seed length: {}", aligner_config.max_seed_length);
    logger->trace("\t Max num seeds per locus: {}", aligner_config.max_num_seeds_per_locus);
    logger->trace("\t Max num nodes per sequence char: {}", aligner_config.max_nodes_per_seq_char);
    logger->trace("\t Gap opening penalty: {}", int64_t(aligner_config.gap_opening_penalty));
    logger->trace("\t Gap extension penalty: {}", int64_t(aligner_config.gap_extension_penalty));
    logger->trace("\t Min DP table cell score: {}", int64_t(aligner_config.min_cell_score));
    logger->trace("\t Min alignment score: {}", aligner_config.min_path_score);
    logger->trace("\t Bandwidth: {}", aligner_config.bandwidth);
    logger->trace("\t X drop-off: {}", aligner_config.xdrop);
    logger->trace("\t Exact k-mer match fraction: {}", aligner_config.exact_kmer_match_fraction);

    logger->trace("\t Scoring matrix: {}", config.alignment_edit_distance ? "unit costs" : "matrix");
    if (!config.alignment_edit_distance) {
        logger->trace("\t\t Match score: {}", int64_t(config.alignment_match_score));
        logger->trace("\t\t (DNA) Transition score: {}",
                      int64_t(config.alignment_mm_transition_score));
        logger->trace("\t\t (DNA) Transversion score: {}",
                      int64_t(config.alignment_mm_transversion_score));
    }

    aligner_config.set_scoring_matrix();

    return aligner_config;
}

std::unique_ptr<IDBGAligner> build_aligner(const DeBruijnGraph &graph, const Config &config) {
    assert(!config.canonical || graph.is_canonical_mode());

    DBGAlignerConfig aligner_config = initialize_aligner_config(graph, config);

    assert(aligner_config.min_seed_length <= aligner_config.max_seed_length);

    if (aligner_config.min_seed_length < graph.get_k()) {
        // seeds are ranges of nodes matching a suffix
        if (!dynamic_cast<const DBGSuccinct*>(&graph)) {
            logger->error("SuffixSeeder can be used only with succinct graph representation");
            exit(1);
        }

        // Use the seeder that seeds to node suffixes
        return std::make_unique<DBGAligner<SuffixSeeder<>>>(graph, aligner_config);

    } else if (aligner_config.max_seed_length == graph.get_k()) {
        assert(aligner_config.min_seed_length == graph.get_k());

        // seeds are single k-mers
        return std::make_unique<DBGAligner<>>(graph, aligner_config);

    } else {
        // seeds are maximal matches within unitigs (uni-MEMs)
        return std::make_unique<DBGAligner<UniMEMSeeder<>>>(graph, aligner_config);
    }
}

void map_sequences_in_file(const std::string &file,
                           const DeBruijnGraph &graph,
                           std::shared_ptr<DBGSuccinct> dbg,
                           const Config &config,
                           const Timer &timer,
                           ThreadPool *thread_pool = nullptr,
                           std::mutex *print_mutex = nullptr) {
    // TODO: multithreaded
    std::ignore = std::tie(thread_pool, print_mutex);

    std::ostream *out = config.outfbase.size()
        ? new std::ofstream(config.outfbase)
        : &std::cout;

    Timer data_reading_timer;

    seq_io::read_fasta_file_critical(file, [&](kseq_t *read_stream) {
        logger->trace("Sequence: {}", read_stream->seq.s);

        if (config.query_presence
                && config.alignment_length == graph.get_k()) {

            bool found = graph.find(read_stream->seq.s,
                                    config.discovery_fraction);

            if (!config.filter_present) {
                *out << found << "\n";

            } else if (found) {
                *out << ">" << read_stream->name.s << "\n"
                            << read_stream->seq.s << "\n";
            }

            return;
        }

        assert(config.alignment_length <= graph.get_k());

        std::vector<DeBruijnGraph::node_index> graphindices;
        if (config.alignment_length == graph.get_k()) {
            graph.map_to_nodes(read_stream->seq.s,
                               [&](const auto &node) {
                                   graphindices.emplace_back(node);
                               });
        } else if (config.query_presence || config.count_kmers) {
            // TODO: make more efficient
            // TODO: canonicalization
            for (size_t i = 0; i + graph.get_k() <= read_stream->seq.l; ++i) {
                dbg->call_nodes_with_suffix(
                    std::string_view(read_stream->seq.s + i, config.alignment_length),
                    [&](auto node, auto) {
                        if (graphindices.empty())
                            graphindices.emplace_back(node);
                    },
                    config.alignment_length
                );
            }
        }

        size_t num_discovered = std::count_if(graphindices.begin(), graphindices.end(),
                                              [](const auto &x) { return x > 0; });

        const size_t num_kmers = graphindices.size();

        if (config.query_presence) {
            const size_t min_kmers_discovered =
                num_kmers - num_kmers * (1 - config.discovery_fraction);
            if (config.filter_present) {
                if (num_discovered >= min_kmers_discovered)
                    *out << ">" << read_stream->name.s << "\n"
                                << read_stream->seq.s << "\n";
            } else {
                *out << (num_discovered >= min_kmers_discovered) << "\n";
            }
            return;
        }

        if (config.count_kmers) {
            std::sort(graphindices.begin(), graphindices.end());
            size_t num_unique_matching_kmers = std::inner_product(
                graphindices.begin() + 1, graphindices.end(),
                graphindices.begin(),
                size_t(graphindices.front() != DeBruijnGraph::npos),
                std::plus<size_t>(),
                [](DeBruijnGraph::node_index next, DeBruijnGraph::node_index prev) {
                    return next != DeBruijnGraph::npos && next != prev;
                }
            );
            *out << read_stream->name.s << "\t"
                 << num_discovered << "/" << num_kmers << "/"
                 << num_unique_matching_kmers << "\n";
            return;
        }

        if (config.alignment_length == graph.get_k()) {
            for (size_t i = 0; i < graphindices.size(); ++i) {
                assert(i + config.alignment_length <= read_stream->seq.l);
                *out << std::string_view(read_stream->seq.s + i, config.alignment_length)
                     << ": " << graphindices[i] << "\n";
            }
        } else {
            // map input subsequences to multiple nodes
            for (size_t i = 0; i + graph.get_k() <= read_stream->seq.l; ++i) {
                // TODO: make more efficient
                std::string_view subseq(read_stream->seq.s + i, config.alignment_length);

                dbg->call_nodes_with_suffix(subseq,
                                            [&](auto node, auto) {
                                                *out << subseq << ": "
                                                     << node
                                                     << "\n";
                                            },
                                            config.alignment_length);
            }
        }

    }, config.forward_and_reverse);

    logger->trace("File '{}' processed in {} sec, current mem usage: {} MiB, total time {} sec",
                  file, data_reading_timer.elapsed(), get_curr_RSS() >> 20, timer.elapsed());

    if (config.outfbase.size())
        delete out;
}


int align_to_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase.size());

    // initialize aligner
    auto graph = load_critical_dbg(config->infbase);
    auto dbg = std::dynamic_pointer_cast<DBGSuccinct>(graph);

    // This speeds up mapping, and allows for node suffix matching
    if (dbg)
        dbg->reset_mask();

    if (config->canonical) {
        logger->trace("Loading as canonical DBG");
        graph.reset(new CanonicalDBG(graph, config->kmers_in_single_form));
    }

    Timer timer;
    ThreadPool thread_pool(get_num_threads());
    std::mutex print_mutex;

    if (config->map_sequences) {
        if (!config->alignment_length) {
            config->alignment_length = graph->get_k();
        } else if (config->alignment_length > graph->get_k()) {
            logger->warn("Mapping to k-mers longer than k is not supported");
            config->alignment_length = graph->get_k();
        }

        if ((!dbg || std::dynamic_pointer_cast<const CanonicalDBG>(graph))
                && config->alignment_length != graph->get_k()) {
            logger->error("Matching k-mers shorter than k only "
                          "supported for DBGSuccinct without --canonical flag");
            exit(1);
        }

        logger->trace("Map sequences against the de Bruijn graph with k={}",
                      graph->get_k());
        logger->trace("Length of mapped k-mers: {}", config->alignment_length);

        for (const auto &file : files) {
            logger->trace("Map sequences from file '{}'", file);

            map_sequences_in_file(file,
                                  *graph,
                                  dbg,
                                  *config,
                                  timer,
                                  &thread_pool,
                                  &print_mutex);
        }

        thread_pool.join();

        return 0;
    }

    auto aligner = build_aligner(*graph, *config);

    if (aligner->get_config().min_seed_length < graph->get_k()
            && std::dynamic_pointer_cast<const CanonicalDBG>(graph)) {
        logger->error("Seeds of length < k not supported with --canonical flag");
        exit(1);
    }

    for (const auto &file : files) {
        logger->trace("Align sequences from file '{}'", file);

        Timer data_reading_timer;

        std::ostream *out = config->outfbase.size()
            ? new std::ofstream(config->outfbase)
            : &std::cout;

        seq_io::read_fasta_file_critical(file, [&](kseq_t *read_stream) {
            thread_pool.enqueue([&](const std::string &query, const std::string &header) {
                auto paths = aligner->align(query);

                std::lock_guard<std::mutex> lock(print_mutex);
                if (!config->output_json) {
                    *out << header << "\t" << paths.get_query();
                    if (paths.empty()) {
                        *out << "\t*\t*\t" << config->alignment_min_path_score
                             << "\t*\t*\t*";
                    } else {
                        for (const auto &path : paths) {
                            *out << "\t" << path;
                        }
                    }

                    *out << "\n";
                } else {
                    Json::StreamWriterBuilder builder;
                    builder["indentation"] = "";

                    bool secondary = false;
                    for (const auto &path : paths) {
                        const auto& path_query = path.get_orientation()
                            ? paths.get_query_reverse_complement()
                            : paths.get_query();

                        *out << Json::writeString(builder,
                                                  path.to_json(path_query,
                                                               *graph,
                                                               secondary,
                                                               header)) << "\n";

                        secondary = true;
                    }

                    if (paths.empty()) {
                        *out << Json::writeString(builder,
                                                  DBGAligner<>::DBGAlignment().to_json(
                                                      query,
                                                      *graph,
                                                      secondary,
                                                      header)
                                                  ) << "\n";
                    }
                }
            }, std::string(read_stream->seq.s),
               config->fasta_anno_comment_delim != Config::UNINITIALIZED_STR
                   && read_stream->comment.l
                       ? utils::join_strings(
                           { read_stream->name.s, read_stream->comment.s },
                           config->fasta_anno_comment_delim,
                           true)
                       : std::string(read_stream->name.s));
        });

        thread_pool.join();

        logger->trace("File '{}' processed in {} sec, "
                      "current mem usage: {} MiB, total time {} sec",
                      file, data_reading_timer.elapsed(),
                      get_curr_RSS() >> 20, timer.elapsed());

        if (config->outfbase.size())
            delete out;
    }

    return 0;
}

} // namespace cli
} // namespace mtg
