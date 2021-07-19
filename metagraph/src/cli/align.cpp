#include "align.hpp"

#include <tsl/ordered_set.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/alignment/aligner_seeder_methods.hpp"
#include "graph/alignment/aligner_extender_methods.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"

namespace mtg {
namespace cli {

using namespace mtg::graph;
using namespace mtg::graph::align;

using mtg::seq_io::kseq_t;
using mtg::common::logger;


DBGAlignerConfig initialize_aligner_config(const Config &config) {
    assert(config.alignment_num_alternative_paths);

    DBGAlignerConfig aligner_config;

    aligner_config.num_alternative_paths = config.alignment_num_alternative_paths;
    aligner_config.min_seed_length = config.alignment_min_seed_length;
    aligner_config.max_seed_length = config.alignment_max_seed_length;
    aligner_config.max_num_seeds_per_locus = config.alignment_max_num_seeds_per_locus;
    aligner_config.max_nodes_per_seq_char = config.alignment_max_nodes_per_seq_char;
    aligner_config.max_ram_per_alignment = config.alignment_max_ram;
    aligner_config.min_path_score = config.alignment_min_path_score;
    aligner_config.xdrop = config.alignment_xdrop;
    aligner_config.min_exact_match = config.alignment_min_exact_match;
    aligner_config.gap_opening_penalty = -config.alignment_gap_opening_penalty;
    aligner_config.gap_extension_penalty = -config.alignment_gap_extension_penalty;
    aligner_config.forward_and_reverse_complement = !config.align_one_strand;
    aligner_config.alignment_edit_distance = config.alignment_edit_distance;
    aligner_config.alignment_match_score = config.alignment_match_score;
    aligner_config.alignment_mm_transition_score = config.alignment_mm_transition_score;
    aligner_config.alignment_mm_transversion_score = config.alignment_mm_transversion_score;
    aligner_config.rel_score_cutoff = config.alignment_rel_score_cutoff;

    logger->trace("Alignment settings:");
    logger->trace("\t Alignments to report: {}", aligner_config.num_alternative_paths);
    logger->trace("\t Min seed length: {}", aligner_config.min_seed_length);
    logger->trace("\t Max seed length: {}", aligner_config.max_seed_length);
    logger->trace("\t Max num seeds per locus: {}", aligner_config.max_num_seeds_per_locus);
    logger->trace("\t Max num nodes per sequence char: {}", aligner_config.max_nodes_per_seq_char);
    logger->trace("\t Max RAM per alignment: {}", aligner_config.max_ram_per_alignment);
    logger->trace("\t Gap opening penalty: {}", int64_t(aligner_config.gap_opening_penalty));
    logger->trace("\t Gap extension penalty: {}", int64_t(aligner_config.gap_extension_penalty));
    logger->trace("\t Min alignment score: {}", aligner_config.min_path_score);
    logger->trace("\t X drop-off: {}", aligner_config.xdrop);
    logger->trace("\t Exact nucleotide match threshold: {}", aligner_config.min_exact_match);

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

template <class Graph>
std::unique_ptr<graph::align::IDBGAligner>
build_aligner(const Graph &graph, const DBGAlignerConfig &aligner_config) {
    return std::make_unique<DBGAligner<>>(graph, aligner_config);
}

template std::unique_ptr<IDBGAligner> build_aligner<DeBruijnGraph>(const DeBruijnGraph &, const DBGAlignerConfig &);

void map_sequences_in_file(const std::string &file,
                           const DeBruijnGraph &graph,
                           const Config &config,
                           const Timer &timer,
                           ThreadPool *thread_pool = nullptr,
                           std::mutex *print_mutex = nullptr) {
    // TODO: multithreaded
    std::ignore = std::tie(thread_pool, print_mutex);

    const DBGSuccinct *dbg = dynamic_cast<const DBGSuccinct*>(&graph);
    if (!dbg) {
        if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph))
            dbg = dynamic_cast<const DBGSuccinct*>(&canonical->get_graph());
    }

    std::unique_ptr<std::ofstream> ofile;
    if (config.outfbase.size())
        ofile = std::make_unique<std::ofstream>(config.outfbase);

    std::ostream *out = ofile ? ofile.get() : &std::cout;

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
            if (dbg->get_mode() == CanonicalDBG::PRIMARY)
                logger->warn("Sub-k-mers will be mapped to unwrapped primary graph");

            for (size_t i = 0; i + graph.get_k() <= read_stream->seq.l; ++i) {
                dbg->call_nodes_with_suffix_matching_longest_prefix(
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
            size_t num_unique_matching_kmers = 0;
            auto prev = DeBruijnGraph::npos;
            for (auto i : graphindices) {
                if (i != DeBruijnGraph::npos && i != prev)
                    ++num_unique_matching_kmers;

                prev = i;
            }
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
                // TODO: canonicalization
                std::string_view subseq(read_stream->seq.s + i, config.alignment_length);

                dbg->call_nodes_with_suffix_matching_longest_prefix(
                    subseq,
                    [&](auto node, auto) {
                        *out << subseq << ": " << node << "\n";
                    },
                    config.alignment_length
                );
            }
        }

    }, config.forward_and_reverse);

    logger->trace("File {} processed in {} sec, current mem usage: {} MB, total time {} sec",
                  file, data_reading_timer.elapsed(), get_curr_RSS() / 1e6, timer.elapsed());
}

std::string sequence_to_gfa_path(const std::string &seq,
                                 const size_t seq_id,
                                 const DeBruijnGraph &graph,
                                 const tsl::ordered_set<uint64_t> &is_unitig_end_node,
                                 const Config *config) {
    auto path_nodes = map_sequence_to_nodes(graph, seq);

    std::string nodes_on_path;
    std::string cigars_on_path;
    const size_t overlap = graph.get_k() - 1;

    for (size_t i = 0; i < path_nodes.size() - 1; ++i) {
        if (config->output_compacted && !is_unitig_end_node.count(path_nodes[i])) {
            continue;
        }
        nodes_on_path += fmt::format("{}+,", path_nodes[i]);
        cigars_on_path += fmt::format("{}M,", overlap);
    }
    uint64_t last_node_to_print = path_nodes.back();
    // We need to print the id of the last unitig even in the case that
    // this unitig is not completely covered by the query sequence.
    while (config->output_compacted && !is_unitig_end_node.count(last_node_to_print)) {
        uint64_t unique_next_node;
        graph.adjacent_outgoing_nodes(
            last_node_to_print,
            [&](uint64_t node) { unique_next_node = node; }
        );
        last_node_to_print = unique_next_node;
    }
    nodes_on_path += fmt::format("{}+,", last_node_to_print);

    // Remove right trailing comma.
    nodes_on_path.pop_back();
    if (cigars_on_path.size()) {
        cigars_on_path.pop_back();
    }
    return fmt::format("P\t{}\t{}\t{}\n", seq_id, nodes_on_path, cigars_on_path);
}

void gfa_map_files(const Config *config,
                   const std::vector<std::string> &files,
                   const DeBruijnGraph &graph) {
    logger->trace("Starting GFA mapping:");

    tsl::ordered_set<uint64_t> is_unitig_end_node;

    graph.call_unitigs(
        [&](const auto &, const auto &path) {
            is_unitig_end_node.insert(path.back());
        },
        get_num_threads()
    );

    std::ofstream gfa_file(utils::remove_suffix(config->outfbase, ".gfa", ".path") + ".path.gfa");

    for (const std::string &file : files) {
        logger->trace("Loading sequences from FASTA file {} to append GFA paths.", file);

        std::vector<string> seq_queries;
        seq_io::FastaParser fasta_parser(file, false);
        for (const seq_io::kseq_t &kseq : fasta_parser) {
            seq_queries.push_back(kseq.seq.s);
        }
        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic) shared(gfa_file)
        for (size_t i = 0; i < seq_queries.size(); ++i) {
            std::string path_string_gfa
                = sequence_to_gfa_path(seq_queries[i], i + 1, graph, is_unitig_end_node, config);
            #pragma omp critical
            gfa_file << path_string_gfa;
        }
    }
}

std::string format_alignment(std::string_view header,
                             const QueryAlignment &paths,
                             const DeBruijnGraph &graph,
                             const Config &config) {
    std::string sout;
    if (!config.output_json) {
        sout += fmt::format("{}\t{}", header, paths.get_query());
        if (paths.empty()) {
            sout += fmt::format("\t*\t*\t{}\t*\t*\t*", config.alignment_min_path_score);
        } else {
            for (const auto &path : paths) {
                sout += fmt::format("\t{}", path);
            }
        }

        sout += "\n";
    } else {
        Json::StreamWriterBuilder builder;
        builder["indentation"] = "";

        bool secondary = false;
        for (size_t i = 0; i < paths.size(); ++i) {
            const auto &path = paths[i];

            Json::Value json_line = path.to_json(paths.get_query(path.get_orientation()),
                                                 graph, secondary, header);

            sout += fmt::format("{}\n", Json::writeString(builder, json_line));
            secondary = true;
        }

        if (paths.empty()) {
            Json::Value json_line
                    = Alignment().to_json(paths.get_query(), graph, secondary, header);

            sout += fmt::format("{}\n", Json::writeString(builder, json_line));
        }
    }

    return sout;
}

int align_to_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase.size());

    // initialize graph
    auto graph = load_critical_dbg(config->infbase);
    auto base_graph = graph;

    if (utils::ends_with(config->outfbase, ".gfa")) {
        gfa_map_files(config, files, *graph);
        return 0;
    }

    Timer timer;
    ThreadPool thread_pool(get_num_threads());
    std::mutex print_mutex;

    if (graph->get_mode() == DeBruijnGraph::PRIMARY) {
        logger->trace("Primary graph wrapped into canonical");
        graph = std::make_shared<CanonicalDBG>(graph);
    }

    if (config->map_sequences) {
        if (!config->alignment_length) {
            config->alignment_length = graph->get_k();
        } else if (config->alignment_length > graph->get_k()) {
            logger->warn("Mapping to k-mers longer than k is not supported");
            config->alignment_length = graph->get_k();
        } else if (config->alignment_length != graph->get_k()
                && !dynamic_cast<const DBGSuccinct*>(base_graph.get())) {
            logger->error("Matching k-mers shorter than k only supported for succinct graphs");
            exit(1);
        }

        logger->trace("Map sequences against the de Bruijn graph with k={}",
                      graph->get_k());
        logger->trace("Length of mapped k-mers: {}", config->alignment_length);

        for (const auto &file : files) {
            logger->trace("Map sequences from file {}", file);

            map_sequences_in_file(file, *graph, *config, timer, &thread_pool, &print_mutex);
        }

        thread_pool.join();

        return 0;
    }

    DBGAlignerConfig aligner_config = initialize_aligner_config(*config);

    for (const auto &file : files) {
        logger->trace("Align sequences from file {}", file);
        seq_io::FastaParser fasta_parser(file, config->forward_and_reverse);
        bool is_reverse_complement = false;

        Timer data_reading_timer;

        std::shared_ptr<std::ostream> out{ std::shared_ptr<std::ostream>{}, &std::cout };
        if (config->outfbase.size())
            out = std::make_shared<std::ofstream>(config->outfbase);

        const uint64_t batch_size = config->query_batch_size_in_bytes;

        auto it = fasta_parser.begin();
        auto end = fasta_parser.end();

        size_t num_batches = 0;

        while (it != end) {
            uint64_t num_bytes_read = 0;

            // Read a batch to pass on to a thread
            typedef std::vector<IDBGAligner::Query> SeqBatch;
            SeqBatch seq_batch;
            num_bytes_read = 0;
            for ( ; it != end && num_bytes_read <= batch_size; ++it) {
                std::string header
                    = config->fasta_anno_comment_delim != Config::UNINITIALIZED_STR
                        && it->comment.l
                            ? utils::join_strings({ it->name.s, it->comment.s },
                                                  config->fasta_anno_comment_delim,
                                                  true)
                            : std::string(it->name.s);
                seq_batch.emplace_back(std::move(header), it->seq.s, is_reverse_complement);
                is_reverse_complement ^= config->forward_and_reverse;
                num_bytes_read += it->seq.l;
            }

            auto process_batch = [&](SeqBatch batch) {
                auto aln_graph = graph;
                if (auto *canonical = dynamic_cast<CanonicalDBG*>(graph.get()))
                    aln_graph = std::make_shared<CanonicalDBG>(*canonical);

                std::unique_ptr<IDBGAligner> aligner;
                aligner = build_aligner(*aln_graph, aligner_config);

                aligner->align_batch(batch, [&](std::string_view header, auto&& paths) {
                    std::string res = format_alignment(header, paths, *aln_graph, *config);
                    std::lock_guard<std::mutex> lock(print_mutex);
                    *out << res;
                });
            };

            ++num_batches;

            uint64_t mbatch_size = it == end && num_batches < get_num_threads()
                ? num_bytes_read / std::max(get_num_threads() - num_batches,
                                            static_cast<size_t>(1))
                : 0;

            if (mbatch_size) {
                // split remaining batch
                logger->trace("Splitting final batch into minibatches");

                auto it = seq_batch.begin();
                auto b_end = seq_batch.end();
                uint64_t num_minibatches = 0;
                while (it != b_end) {
                    uint64_t cur_minibatch_read = 0;
                    auto last_mv_it = std::make_move_iterator(it);
                    for ( ; it != b_end && cur_minibatch_read < mbatch_size; ++it) {
                        cur_minibatch_read += std::get<1>(*it).size();
                    }

                    thread_pool.enqueue(process_batch,
                                        SeqBatch(last_mv_it, std::make_move_iterator(it)));
                    ++num_minibatches;
                }

                logger->trace("Num minibatches: {}, minibatch size: {} KB",
                              num_minibatches, mbatch_size / 1e3);
            } else {
                thread_pool.enqueue(process_batch, std::move(seq_batch));
            }
        };

        thread_pool.join();

        logger->trace("File {} processed in {} sec, "
                      "num batches: {}, batch size: {} KB, "
                      "current mem usage: {} MB, total time {} sec",
                      file, data_reading_timer.elapsed(), num_batches, batch_size / 1e3,
                      get_curr_RSS() / 1e6, timer.elapsed());
    }

    return 0;
}

} // namespace cli
} // namespace mtg
