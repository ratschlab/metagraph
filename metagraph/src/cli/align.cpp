#include "align.hpp"

#include <tsl/ordered_set.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/alignment/aligner_labeled.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/graph_extensions/node_rc.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"
#include "graph/graph_extensions/unitigs.hpp"
#include "graph/graph_extensions/hll_wrapper.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"

namespace mtg {
namespace cli {

using namespace mtg::graph;
using namespace mtg::graph::align;

using mtg::seq_io::kseq_t;
using mtg::common::logger;


DBGAlignerConfig initialize_aligner_config(const Config &config) {
    assert(config.alignment_num_alternative_paths);

    DBGAlignerConfig c = {
        .num_alternative_paths = config.alignment_num_alternative_paths,
        .min_seed_length = config.alignment_min_seed_length,
        .max_seed_length = config.alignment_max_seed_length,
        .max_num_seeds_per_locus = config.alignment_max_num_seeds_per_locus,
        .min_path_score = config.alignment_min_path_score,
        .xdrop = config.alignment_xdrop,
        .min_exact_match = config.alignment_min_exact_match,
        .max_nodes_per_seq_char = config.alignment_max_nodes_per_seq_char,
        .max_ram_per_alignment = config.alignment_max_ram,
        .rel_score_cutoff = config.alignment_rel_score_cutoff,
        .gap_opening_penalty = static_cast<int8_t>(-config.alignment_gap_opening_penalty),
        .gap_extension_penalty = static_cast<int8_t>(-config.alignment_gap_extension_penalty),
        .left_end_bonus = config.alignment_end_bonus,
        .right_end_bonus = config.alignment_end_bonus,
        .label_change_score = config.alignment_label_change_score,
        .forward_and_reverse_complement = !config.align_only_forwards,
        .chain_alignments = config.alignment_chain,
        .post_chain_alignments = config.alignment_post_chain,
        .no_seed_complexity_filter = config.alignment_no_seed_complexity_filter,
        .label_change_union = config.alignment_label_change_union,
        .alignment_edit_distance = config.alignment_edit_distance,
        .alignment_match_score = config.alignment_match_score,
        .alignment_mm_transition_score = config.alignment_mm_transition_score,
        .alignment_mm_transversion_score = config.alignment_mm_transversion_score,
        .score_matrix = DBGAlignerConfig::ScoreMatrix{},
    };

    c.set_scoring_matrix();

    c.print_summary();

    return c;
}

void map_sequences_in_file(const std::string &file,
                           const DeBruijnGraph &graph,
                           const Config &config,
                           const Timer &timer,
                           ThreadPool *thread_pool = nullptr,
                           std::mutex *print_mutex = nullptr) {
    // TODO: multithreaded
    std::ignore = std::tie(thread_pool, print_mutex);

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

        std::vector<DeBruijnGraph::node_index> nodes;
        if (config.alignment_length == graph.get_k()) {
            nodes = map_to_nodes(graph, read_stream->seq.s);

        } else {
            // TODO: make more efficient
            // TODO: canonicalization
            const DBGSuccinct &dbg = static_cast<const DBGSuccinct&>(graph);
            if (dbg.get_mode() == DeBruijnGraph::PRIMARY)
                logger->warn("Sub-k-mers will be mapped to unwrapped primary graph");

            for (size_t i = 0; i + config.alignment_length <= read_stream->seq.l; ++i) {
                nodes.emplace_back(DeBruijnGraph::npos);
                dbg.call_nodes_with_suffix_matching_longest_prefix(
                    std::string_view(read_stream->seq.s + i, config.alignment_length),
                    [&](auto node, auto) {
                        if (nodes.back() == DeBruijnGraph::npos)
                            nodes.back() = node;
                    },
                    config.alignment_length
                );
            }
        }

        size_t num_discovered = std::count_if(nodes.begin(), nodes.end(),
                                              [](const auto &x) { return x > 0; });

        const size_t num_kmers = nodes.size();

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
            std::sort(nodes.begin(), nodes.end());
            size_t num_unique_matching_kmers = 0;
            auto prev = DeBruijnGraph::npos;
            for (auto i : nodes) {
                if (i != DeBruijnGraph::npos && i != prev)
                    ++num_unique_matching_kmers;

                prev = i;
            }
            *out << fmt::format("{}\t{}/{}/{}\n", read_stream->name.s,
                                num_discovered, num_kmers, num_unique_matching_kmers);
            return;
        }

        // mapping of each k-mer to a graph node
        for (size_t i = 0; i < nodes.size(); ++i) {
            assert(i + config.alignment_length <= read_stream->seq.l);
            *out << fmt::format("{}: {}\n", std::string_view(read_stream->seq.s + i,
                                                             config.alignment_length),
                                            nodes[i]);
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
    auto path_nodes = map_to_nodes_sequentially(graph, seq);

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

std::string format_alignment(const std::string &header,
                             const AlignmentResults &paths,
                             const DeBruijnGraph &graph,
                             const Config &config) {
    std::string sout;
    if (config.output_json) {
        Json::StreamWriterBuilder builder;
        builder["indentation"] = "";

        bool secondary = false;
        for (size_t i = 0; i < paths.size(); ++i) {
            const auto &path = paths[i];

            Json::Value json_line = path.to_json(graph.get_k(), secondary, header);
            sout += fmt::format("{}\n", Json::writeString(builder, json_line));
            secondary = true;
        }

        if (paths.empty()) {
            Json::Value json_line = Alignment().to_json(graph.get_k(), secondary, header);
            sout += fmt::format("{}\n", Json::writeString(builder, json_line));
        }
    } else if (config.output_gaf) {
        for (size_t i = 0; i < paths.size(); ++i) {
            sout += fmt::format("{}\n", paths[i].to_gaf(header));
        }
    } else {
        sout += fmt::format("{}\t{}", header, paths.get_query());
        if (paths.empty()) {
            sout += fmt::format("\t*\t*\t{}\t*\t*\t*\n", config.alignment_min_path_score);
        } else {
            for (const auto &path : paths) {
                sout += fmt::format("\t{}", path);
            }
            sout += "\n";
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

    if (utils::ends_with(config->outfbase, ".gfa")) {
        gfa_map_files(config, files, *graph);
        return 0;
    }

    auto hll = std::make_shared<HLLWrapper<>>();
    if (config->alignment_post_chain) {
        if (hll->load(utils::remove_suffix(config->infbase, graph->file_extension()))) {
            logger->trace("Loaded HLL sketch with {} columns", hll->data().num_columns());
        } else {
            logger->warn("HLL sketch not present or failed to load");
            hll.reset();
        }
    } else {
        hll.reset();
    }

    // For graphs which still feature a mask, this speeds up mapping and allows
    // for dummy nodes to be matched by suffix seeding
    auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(graph);
    std::shared_ptr<Unitigs> graph_unitigs;
    if (dbg_succ) {
        dbg_succ->reset_mask();
        if (dbg_succ->get_mode() == DeBruijnGraph::PRIMARY) {
            auto node_rc = std::make_shared<NodeRC>(*dbg_succ);
            if (node_rc->load(config->infbase)) {
                logger->trace("Loaded the adj-rc index (adjacent to reverse-complement nodes)");
                dbg_succ->add_extension(node_rc);
            } else {
                logger->warn("adj-rc index missing or failed to load. "
                             "Alignment speed will be significantly slower. "
                             "Use metagraph transform to generate an adj-rc index.");
            }
        }

        graph_unitigs = std::make_shared<Unitigs>(*dbg_succ);
        if (graph_unitigs->load(config->infbase)) {
            logger->trace("Loaded the unitig index");
        } else {
            graph_unitigs.reset();
            logger->warn("Unitig index missing or failed to load. "
                         "Alignment speed will be significantly slower. "
                         "Use metagraph transform to generate a unitig index.");
        }
    }

    Timer timer;
    ThreadPool thread_pool(get_num_threads());
    std::mutex print_mutex;

    if (config->map_sequences) {
        if (graph->get_mode() == DeBruijnGraph::PRIMARY) {
            graph = std::make_shared<CanonicalDBG>(graph);
            logger->trace("Primary graph wrapped into canonical");
        }

        if (!config->alignment_length) {
            config->alignment_length = graph->get_k();
        } else if (config->alignment_length > graph->get_k()) {
            logger->warn("Mapping to k-mers longer than k is not supported");
            config->alignment_length = graph->get_k();
        } else if (config->alignment_length != graph->get_k() && !dbg_succ) {
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

    std::unique_ptr<AnnotatedDBG> anno_dbg;
    if (config->infbase_annotators.size()) {
        assert(config->infbase_annotators.size() == 1);
        anno_dbg = initialize_annotated_dbg(graph, *config);
    }

    for (const auto &file : files) {
        logger->trace("Align sequences from file {}", file);
        seq_io::FastaParser fasta_parser(file, config->forward_and_reverse);

        Timer data_reading_timer;

        std::unique_ptr<std::ofstream> ofile;
        if (config->outfbase.size())
            ofile = std::make_unique<std::ofstream>(config->outfbase);

        std::ostream *out = ofile ? ofile.get() : &std::cout;

        const uint64_t batch_size = config->query_batch_size_in_bytes;

        auto it = fasta_parser.begin();
        auto end = fasta_parser.end();

        size_t num_batches = 0;
        size_t seq_id = 0;

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
                                                  config->fasta_anno_comment_delim, true)
                            : std::string(it->name.s);
                seq_batch.emplace_back(std::move(header), it->seq.s);
                num_bytes_read += it->seq.l;
                ++seq_id;
            }

            ++num_batches;
            thread_pool.enqueue([&,graph,seq_id,batch=std::move(seq_batch)]() {
                // Make a dummy shared_ptr
                auto aln_graph
                    = std::shared_ptr<DeBruijnGraph>(std::shared_ptr<DeBruijnGraph>{}, graph.get());

                // Wrap it in CanonicalDBG if needed. This way, each thread gets its
                // own wrapper (and more importantly, its own local cache).
                if (aln_graph->get_mode() == DeBruijnGraph::PRIMARY) {
                    aln_graph = std::make_shared<CanonicalDBG>(aln_graph);
                    logger->trace("Primary graph wrapped into canonical");
                    // If backwards traversal on DBGSuccinct will be needed, then
                    // add a cache to speed it up.
                    if (dbg_succ)
                        aln_graph->add_extension(std::make_shared<NodeFirstCache>(*dbg_succ));
                }

                if (hll)
                    aln_graph->add_extension(hll);

                if (graph_unitigs)
                    aln_graph->add_extension(graph_unitigs);

                std::unique_ptr<IDBGAligner> aligner;

                if (anno_dbg) {
                    aligner = std::make_unique<LabeledAligner<>>(*aln_graph, aligner_config,
                                                                 anno_dbg->get_annotator());
                } else {
                    aligner = std::make_unique<DBGAligner<>>(*aln_graph, aligner_config);
                }

                aligner->align_batch(batch,
                    [&](const std::string &header, AlignmentResults&& paths) {
                        const auto &res = format_alignment(header, paths, *graph, *config);
                        std::lock_guard<std::mutex> lock(print_mutex);
                        *out << res;
                    }, seq_id - batch.size()
                );
            });
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
