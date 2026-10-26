#include "align.hpp"
#include "alignment_config.hpp"

#include <tsl/hopscotch_set.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/alignment/alignment.hpp"
#include "graph/alignment_redone/aligner_config.hpp"
#include "graph/alignment_redone/aln_query.hpp"
#include "graph/alignment_redone/aln_seeder.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"

namespace mtg {
namespace cli {

using namespace mtg::graph;

using mtg::seq_io::kseq_t;
using mtg::common::logger;

// after how many threads should per-thread traversal caches be used instead of
// a single shared cache
const size_t SEPARATE_CACHES_SWITCHOVER_THREADS = 8;


align_redone::DBGAlignerConfig initialize_aligner_config(const Config &config,
                                           const DeBruijnGraph &graph) {
    assert(config.alignment_num_alternative_paths);
    if (config.alignment_post_chain
            && (config.alignment_chain || config.alignment_two_sided || !config.alignment_extend_chains))
        throw std::invalid_argument("Post-chaining cannot be combined with seed chaining, two-sided extension, or disabled chain extension");

    auto c = detail::make_alignment_config<align_redone::DBGAlignerConfig>(config,
        config.alignment_min_seed_length
            ? std::min(config.alignment_min_seed_length, graph.get_k()) : graph.get_k());

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
                                 const tsl::hopscotch_set<uint64_t> &is_unitig_end_node,
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

    tsl::hopscotch_set<uint64_t> is_unitig_end_node;

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

namespace {

std::string format_alignment(const std::string &header, const std::string &query,
                             const std::vector<align_redone::Alignment> &paths,
                             const DeBruijnGraph &graph, const Config &config,
                             const align_redone::DBGAlignerConfig &scoring,
                             const align_redone::AnnotationBuffer *annotations,
                             const annot::CoordToHeader *cth) {
    const auto *encoder = annotations
        ? &annotations->get_annotator().get_label_encoder() : nullptr;
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    std::ostringstream output;
    if (!config.output_json)
        output << header << '\t' << query;
    bool secondary = false;
    for (const auto &path : paths) {
        // Reuse the established output schema and coordinate/header mapping.
        graph::align::Alignment displayed(
            path.get_seed(), std::vector(path.get_path()), std::string(path.get_spelling()),
            align_redone::score_match(path, scoring),
            graph::align::Cigar(path.get_cigar().data()), 0, path.get_orientation(), path.get_offset());
        if (annotations && !path.get_label_classes().empty()) {
            const auto &columns = annotations->get_cached_column_set(path.get_label_classes()[0]);
            displayed.label_columns.assign(columns.begin(), columns.end());
            displayed.label_coordinates = path.get_coords();
        }
        if (config.output_json) {
            output << Json::writeString(builder,
                displayed.to_json(graph.get_k(), secondary, header, {}, encoder, cth,
                                  config.align_output_path)) << '\n';
        } else {
            output << '\t' << displayed;
            if (!displayed.label_coordinates.empty()) {
                output << '\t' << (cth ? displayed.format_coords(*cth, graph.get_k())
                                        : displayed.format_coords(*encoder));
            } else if (!displayed.label_columns.empty()) {
                const char *delimiter = "\t";
                for (auto column : displayed.label_columns) {
                    output << delimiter << encoder->decode(column);
                    delimiter = ";";
                }
            }
        }
        secondary = true;
    }
    if (paths.empty()) {
        if (config.output_json) {
            output << Json::writeString(builder,
                graph::align::Alignment().to_json(graph.get_k(), false, header)) << '\n';
        } else {
            output << "\t*\t*\t" << config.alignment_min_path_score << "\t*\t*\t*";
        }
    }
    if (!config.output_json)
        output << '\n';
    return output.str();
}

} // namespace

int align_to_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase.size());

    // For graph-only outputs/queries, avoid starting annotation loading.
    const bool graph_only_path = config->outfbase.ends_with(".gfa")
                                 || config->map_sequences;
    std::shared_ptr<DeBruijnGraph> graph;
    std::future<std::unique_ptr<AnnotatedDBG>> anno_dbg_future;
    if (graph_only_path) {
        graph = load_critical_dbg(config->infbase);
    } else {
        // Load graph and annotation in parallel for alignment mode.
        auto loaded = load_graph_with_async_annotation(*config);
        graph = loaded.first.get();
        anno_dbg_future = std::move(loaded.second);
    }

    if (config->outfbase.ends_with(".gfa")) {
        gfa_map_files(config, files, *graph);
        return 0;
    }

    // For graphs which still feature a mask, this speeds up mapping and allows
    // for dummy nodes to be matched by suffix seeding
    auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(graph);
    if (dbg_succ)
        dbg_succ->reset_mask();

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

    align_redone::DBGAlignerConfig aligner_config = initialize_aligner_config(*config, *graph);

    if (auto sshash = std::dynamic_pointer_cast<DBGSSHash>(graph)) {
        std::visit([&](const auto &dict) {
            if (dict.m() > aligner_config.min_seed_length)
                common::logger->warn("SSHash minimizer size is {}, only seeds of at least this length will be found.", dict.m());
        }, sshash->data());
    }

    assert(config->infbase_annotators.size() <= 1);
    std::unique_ptr<AnnotatedDBG> anno_dbg;
    if (anno_dbg_future.valid())
        anno_dbg = anno_dbg_future.get();

    auto wrap_graph = [&](auto graph) {
        if (graph->get_mode() == DeBruijnGraph::PRIMARY) {
            graph = std::make_shared<CanonicalDBG>(graph);
            logger->trace("Primary graph wrapped into canonical");
        }

        // If backwards traversal on DBGSuccinct will be needed, then add a cache
        // to speed it up.
        if (dbg_succ)
            graph->add_extension(std::make_shared<NodeFirstCache>(*dbg_succ));

        return graph;
    };

    // if using a small number of threads, use a shared cache for all threads
    if (get_num_threads() <= SEPARATE_CACHES_SWITCHOVER_THREADS)
        graph = wrap_graph(graph);

    timer.reset();

    for (const auto &file : files) {
        logger->trace("Align sequences from file {}", file);
        seq_io::FastaParser fasta_parser(file, config->forward_and_reverse);

        Timer data_reading_timer;

        std::unique_ptr<std::ofstream> ofile;
        if (config->outfbase.size())
            ofile = std::make_unique<std::ofstream>(config->outfbase);

        std::ostream *out = ofile ? ofile.get() : &std::cout;

        const uint64_t batch_size = config->query_batch_size;

        auto it = fasta_parser.begin();
        auto end = fasta_parser.end();

        size_t num_bp = 0;
        std::exception_ptr task_error;

        while (it != end) {
            uint64_t num_bytes_read = 0;

            using SeqBatch = std::vector<std::pair<std::string, std::string>>;
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
            }

            auto run_batch = [&,graph,batch=std::move(seq_batch)]() {
                // Make a dummy shared_ptr
                auto aln_graph
                    = std::shared_ptr<DeBruijnGraph>(std::shared_ptr<DeBruijnGraph>{}, graph.get());

                // Give each thread a separate cache if using many threads to avoid
                // contention
                if (aln_graph->get_mode() == DeBruijnGraph::PRIMARY) {
                    assert(get_num_threads() > SEPARATE_CACHES_SWITCHOVER_THREADS);
                    aln_graph = wrap_graph(aln_graph);
                }

                std::optional<align_redone::AnnotationBuffer> anno_buffer;
                if (anno_dbg)
                    anno_buffer.emplace(*aln_graph, anno_dbg->get_annotator());

                for (const auto &[header, query] : batch) {
                    align_redone::Query aln_query(*aln_graph, query);

                    std::unique_ptr<align_redone::ExactSeeder> seeder;
                    if (anno_dbg) {
                        seeder = std::make_unique<align_redone::LabeledSeeder>(*anno_buffer, aln_query, aligner_config);
                    } else {
                        seeder = std::make_unique<align_redone::ExactSeeder>(aln_query, aligner_config);
                    }

                    std::vector<align_redone::Alignment> paths;
                    auto aln_callback = [&](align_redone::Alignment&& aln) {
                        paths.emplace_back(std::move(aln));
                    };

                    if (config->alignment_extend_chains) {
                        std::unique_ptr<align_redone::Extender> extender;
                        if (anno_dbg) {
                            extender = std::make_unique<align_redone::LabeledExtender>(*anno_buffer, aln_query, aligner_config);
                        } else {
                            extender = std::make_unique<align_redone::Extender>(aln_query, aligner_config);
                        }
                        align_redone::align_query(aln_query, *seeder, *extender,
                                                  aln_callback,
                                                  config->alignment_connect_anchors);
                    } else {
                        align_redone::align_query(aln_query, *seeder,
                                                  aln_callback,
                                                  config->alignment_connect_anchors);
                    }

                    auto output = format_alignment(header, aln_query.get_query(), paths, *aln_graph,
                                                   *config, aligner_config,
                                                   anno_buffer ? &*anno_buffer : nullptr,
                                                   anno_dbg ? anno_dbg->get_coord_to_header() : nullptr);
                    std::lock_guard<std::mutex> lock(print_mutex);
                    *out << output;
                }
            };
            // ThreadPool rethrows task exceptions on the worker. Preserve the
            // first one so the caller can report a regular command failure.
            thread_pool.enqueue([&,run_batch=std::move(run_batch)]() {
                try {
                    run_batch();
                } catch (...) {
                    std::lock_guard<std::mutex> lock(print_mutex);
                    if (!task_error)
                        task_error = std::current_exception();
                }
            });

            num_bp += num_bytes_read;
        };

        thread_pool.join();
        if (task_error) {
            try {
                std::rethrow_exception(task_error);
            } catch (const std::exception &e) {
                logger->error("Alignment failed: {}", e.what());
            } catch (...) {
                logger->error("Alignment failed with an unknown error");
            }
            return 1;
        }

        auto time = data_reading_timer.elapsed();
        logger->trace("File {} with {} base pairs was processed in {} sec, throughput: {:.1f} bp/s, "
                      "current mem usage: {} MB, total time {} sec",
                      file, num_bp, time, (double)num_bp / time,
                      get_curr_RSS() / 1e6, timer.elapsed());
    }

    return 0;
}

} // namespace cli
} // namespace mtg
