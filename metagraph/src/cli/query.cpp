#include "query.hpp"

#include <ips4o.hpp>
#include <fmt/format.h>
#include <tsl/ordered_set.h>

#include "common/logger.hpp"
#include "common/hash/hash.hpp"
#include "common/unix_tools.hpp"
#include "common/hash/hash.hpp"
#include "common/utils/template_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "align.hpp"


namespace mtg {
namespace cli {

const size_t kRowBatchSize = 100'000;
const bool kPrefilterWithBloom = true;
const char ALIGNED_SEQ_HEADER_FORMAT[] = "{}:{}:{}:{}";
const double kHLLCounterError = 0.003;

using mtg::common::logger;


QueryExecutor::QueryExecutor(const Config &config,
                             const AnnotatedDBG &anno_graph,
                             const IDBGAligner *aligner,
                             ThreadPool &thread_pool)
      : config_(config),
        anno_graph_(anno_graph),
        aligner_(aligner),
        thread_pool_(thread_pool) {
    if (config_.get_coverage) {
        label_counts_ = anno_graph_.get_annotation().get_label_counts();
        match_counter_.resize(label_counts_.size(), kHLLCounterError);
    }
}


std::string QueryExecutor::execute_query(const std::string &seq_name,
                                         const std::string &sequence,
                                         bool count_labels,
                                         bool print_signature,
                                         bool suppress_unlabeled,
                                         size_t num_top_labels,
                                         double discovery_fraction,
                                         std::string anno_labels_delimiter,
                                         const AnnotatedDBG &anno_graph) {
    std::string output;
    output.reserve(1'000);

    if (print_signature) {
        auto top_labels
            = anno_graph.get_top_label_signatures(sequence,
                                                  num_top_labels,
                                                  discovery_fraction);

        if (!top_labels.size() && suppress_unlabeled)
            return "";

        output += seq_name;

        for (const auto &[label, kmer_presence_mask] : top_labels) {
            output += fmt::format("\t<{}>:{}:{}:{}", label,
                                  sdsl::util::cnt_one_bits(kmer_presence_mask),
                                  sdsl::util::to_string(kmer_presence_mask),
                                  anno_graph.score_kmer_presence_mask(kmer_presence_mask));
        }

        output += '\n';

    } else if (count_labels) {
        auto top_labels = anno_graph.get_top_labels(sequence,
                                                    num_top_labels,
                                                    discovery_fraction);

        if (!top_labels.size() && suppress_unlabeled)
            return "";

        output += seq_name;

        for (const auto &[label, count] : top_labels) {
            output += "\t<";
            output += label;
            output += ">:";
            output += fmt::format_int(count).c_str();
        }

        output += '\n';

    } else {
        auto labels_discovered = anno_graph.get_labels(sequence, discovery_fraction);

        if (!labels_discovered.size() && suppress_unlabeled)
            return "";

        output += seq_name;
        output += '\t';
        output += utils::join_strings(labels_discovered, anno_labels_delimiter);
        output += '\n';
    }

    return output;
}

/**
 * Construct a de Bruijn graph from the query sequences
 * fetched in |call_sequences|.
 *
 *  Algorithm.
 *
 * 1. Index k-mers from the query sequences in a non-canonical query de Bruijn
 *    graph (with pre-filtering by a Bloom filter, if initialized).
 *    This query graph will be rebuilt as a canonical one in step 2.b), if the
 *    full graph is canonical.
 *
 * 2. Extract contigs from this small de Bruijn graph and map them to the full
 *    graph to map each k-mer to its respective annotation row index.
 *    --> here we map each unique k-mer in sequences only once.
 *
 *    (b, canonical) If the full graph is canonical, rebuild the query graph
 *                   in the canonical mode storing all k-mers found in the full
 *                   graph.
 *
 * 3. If |discovery_fraction| is greater than zero, map all query sequences to
 *    the query graph and filter out those having too few k-mer matches.
 *    Then, remove from the query graph those k-mers occurring only in query
 *    sequences that have been filtered out.
 *
 * 4. Extract annotation for the nodes of the query graph and return.
 */
std::pair<std::unique_ptr<AnnotatedDBG>, std::vector<uint64_t>>
construct_query_graph(const AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      double discovery_fraction,
                      size_t num_threads) {
    const auto &full_dbg = anno_graph.get_graph();
    const auto &full_annotation = anno_graph.get_annotation();

    Timer timer;

    // construct graph storing all k-mers in query
    auto graph_init = std::make_shared<DBGHashOrdered>(full_dbg.get_k(), false);

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct *>(&full_dbg);
    if (kPrefilterWithBloom && dbg_succ) {
        if (dbg_succ->get_bloom_filter())
            logger->trace(
                    "[Query graph construction] Started indexing k-mers pre-filtered "
                    "with Bloom filter");

        call_sequences([&](const std::string &sequence) {
            // TODO: implement add_sequence with filter for all graph representations
            graph_init->add_sequence(
                sequence,
                get_missing_kmer_skipper(dbg_succ->get_bloom_filter(), sequence)
            );
        });
    } else {
        call_sequences([&](const std::string &sequence) {
            graph_init->add_sequence(sequence);
        });
    }

    std::shared_ptr<DeBruijnGraph> graph = std::move(graph_init);

    logger->trace("[Query graph construction] k-mer indexing took {} sec", timer.elapsed());
    timer.reset();

    // pull contigs from query graph
    std::vector<std::pair<std::string, std::vector<DeBruijnGraph::node_index>>> contigs;
    graph->call_sequences(
        [&](auto&&... contig_args) { contigs.emplace_back(std::move(contig_args)...); },
        full_dbg.is_canonical_mode()
    );

    logger->trace("[Query graph construction] Contig extraction took {} sec", timer.elapsed());
    timer.reset();

    // map contigs onto the full graph
    std::vector<uint64_t> index_in_full_graph;

    if (full_dbg.is_canonical_mode()) {
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
        for (size_t i = 0; i < contigs.size(); ++i) {
            const std::string &contig = contigs[i].first;
            auto &nodes_in_full = contigs[i].second;

            size_t j = 0;
            full_dbg.map_to_nodes(contig,
                [&](auto node_in_full) { nodes_in_full[j++] = node_in_full; }
            );
            assert(j == nodes_in_full.size());
        }

        logger->trace("[Query graph construction] Contigs mapped to graph in {} sec", timer.elapsed());
        timer.reset();

        // construct canonical graph storing all k-mers found in the full graph
        graph_init = std::make_shared<DBGHashOrdered>(full_dbg.get_k(), true);

        for (size_t i = 0; i < contigs.size(); ++i) {
            const std::string &contig = contigs[i].first;
            const auto &nodes_in_full = contigs[i].second;
            size_t j = 0;
            graph_init->add_sequence(contig, [&]() { return nodes_in_full[j++] == 0; });
        }

        graph = std::move(graph_init);

        logger->trace("[Query graph construction] k-mers reindexed in canonical mode in {} sec",
                      timer.elapsed());
        timer.reset();

        index_in_full_graph.assign(graph->max_index() + 1, 0);

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
        for (size_t i = 0; i < contigs.size(); ++i) {
            const std::string &contig = contigs[i].first;
            const auto &nodes_in_full = contigs[i].second;

            size_t j = 0;
            graph->map_to_nodes(contig, [&](auto node) {
                index_in_full_graph[node] = nodes_in_full[j++];
            });
            assert(j == nodes_in_full.size());
        }

        logger->trace("[Query graph construction] Mapping between graphs constructed in {} sec",
                      timer.elapsed());
        timer.reset();

    } else {
        index_in_full_graph.assign(graph->max_index() + 1, 0);

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
        for (size_t i = 0; i < contigs.size(); ++i) {
            const std::string &contig = contigs[i].first;
            const auto &path = contigs[i].second;

            size_t j = 0;
            full_dbg.map_to_nodes(contig, [&](auto node_in_full) {
                index_in_full_graph[path[j++]] = node_in_full;
            });
            assert(j == path.size());
        }

        logger->trace("[Query graph construction] Contigs mapped to graph in {} sec", timer.elapsed());
        timer.reset();
    }

    contigs = decltype(contigs)();

    assert(!index_in_full_graph.at(0));

    if (discovery_fraction > 0) {
        sdsl::bit_vector mask(graph->max_index() + 1, false);

        call_sequences([&](const std::string &sequence) {
            if (sequence.length() < graph->get_k())
                return;

            const size_t num_kmers = sequence.length() - graph->get_k() + 1;
            const size_t max_kmers_missing = num_kmers * (1 - discovery_fraction);
            const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
            size_t num_kmers_discovered = 0;
            size_t num_kmers_missing = 0;

            std::vector<DeBruijnGraph::node_index> nodes;
            nodes.reserve(num_kmers);

            graph->map_to_nodes(sequence,
                [&](auto node) {
                    if (index_in_full_graph[node]) {
                        num_kmers_discovered++;
                        nodes.push_back(node);
                    } else {
                        num_kmers_missing++;
                    }
                },
                [&]() { return num_kmers_missing > max_kmers_missing
                                || num_kmers_discovered >= min_kmers_discovered; }
            );

            if (num_kmers_missing <= max_kmers_missing) {
                for (auto node : nodes) {
                    mask[node] = true;
                }
            }
        });

        // correcting the mask
        call_zeros(mask, [&](auto i) { index_in_full_graph[i] = 0; });

        graph = std::make_shared<MaskedDeBruijnGraph>(
            graph,
            std::make_unique<bit_vector_stat>(std::move(mask))
        );

        logger->trace("[Query graph construction] Reduced k-mer dictionary in {} sec",
                      timer.elapsed());
        timer.reset();
    }

    std::vector<std::pair<uint64_t, uint64_t>> from_full_to_query;
    from_full_to_query.reserve(index_in_full_graph.size());

    for (uint64_t node = 0; node < index_in_full_graph.size(); ++node) {
        if (index_in_full_graph[node]) {
            from_full_to_query.emplace_back(
                AnnotatedDBG::graph_to_anno_index(index_in_full_graph[node]),
                AnnotatedDBG::graph_to_anno_index(node)
            );
        }
    }

    ips4o::parallel::sort(from_full_to_query.begin(), from_full_to_query.end(),
                          utils::LessFirst(), num_threads);

    logger->trace("[Query graph construction] Prepared row indexes for query {} sec",
                  timer.elapsed());
    timer.reset();

    // initialize fast query annotation
    using RowSet = tsl::ordered_set<SmallVector<uint32_t>,
                                    utils::VectorHash,
                                    std::equal_to<SmallVector<uint32_t>>,
                                    std::allocator<SmallVector<uint32_t>>,
                                    std::vector<SmallVector<uint32_t>>,
                                    uint32_t>;
    RowSet unique_rows { SmallVector<uint32_t>() };
    std::vector<uint32_t> row_rank(graph->max_index(), 0);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint64_t batch_begin = 0;
                        batch_begin < from_full_to_query.size();
                                        batch_begin += kRowBatchSize) {
        const uint64_t batch_end
            = std::min(batch_begin + kRowBatchSize,
                       static_cast<uint64_t>(from_full_to_query.size()));

        std::vector<uint64_t> row_indexes;
        row_indexes.reserve(batch_end - batch_begin);

        for (uint64_t i = batch_begin; i < batch_end; ++i) {
            assert(from_full_to_query[i].first < full_annotation.num_objects());
            row_indexes.push_back(from_full_to_query[i].first);
            if (unique_rows.size() == std::numeric_limits<uint32_t>::max())
                throw std::runtime_error("There must be less than 2^32 unique rows."
                                         " Reduce the query batch size.");
        }

        auto rows = full_annotation.get_matrix().get_rows(row_indexes);

        assert(rows.size() == batch_end - batch_begin);

        #pragma omp critical
        {
            for (uint64_t i = batch_begin; i < batch_end; ++i) {
                const auto &row = rows[i - batch_begin];
                auto it = unique_rows.emplace(row.begin(), row.end()).first;
                row_rank[from_full_to_query[i].second] = it - unique_rows.begin();
            }
        }
    }

    auto annotation_rows = const_cast<std::vector<SmallVector<uint32_t>>&&>(
        unique_rows.values_container()
    );

    // copy annotations from the full graph to the query graph
    auto annotation = std::make_unique<annotate::UniqueRowAnnotator>(
        std::make_unique<UniqueRowBinmat>(std::move(annotation_rows),
                                          std::move(row_rank),
                                          full_annotation.num_labels()),
        full_annotation.get_label_encoder()
    );

    logger->trace("[Query graph construction] Query annotation constructed in {} sec",
                  timer.elapsed());
    timer.reset();

    // build annotated graph from the query graph and copied annotations
    return std::make_pair(
        std::make_unique<AnnotatedDBG>(graph, std::move(annotation)),
        std::move(index_in_full_graph)
    );
}


std::string get_alignment_header_and_swap_query(const std::string &name,
                                                std::string *query_seq,
                                                QueryAlignment<> *matches) {
    std::string header;

    if (matches->size()) {
        // sequence for querying -- the best alignment
        *query_seq = const_cast<std::string&&>((*matches)[0].get_sequence());
        header = fmt::format(ALIGNED_SEQ_HEADER_FORMAT,
                             name, *query_seq, (*matches)[0].get_score(),
                             (*matches)[0].get_cigar().to_string());

    } else {
        // no alignment was found
        // the original sequence `query_seq` will be queried
        header = fmt::format(ALIGNED_SEQ_HEADER_FORMAT,
                             name, *query_seq, 0,
                             fmt::format("{}S", query_seq->length()));
    }

    return header;
}


int query_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase_annotators.size() == 1);

    auto graph = load_critical_dbg(config->infbase);
    auto anno_graph = initialize_annotated_dbg(graph, *config);

    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);

    Timer timer;

    std::unique_ptr<IDBGAligner> aligner;
    if (config->align_sequences) {
        assert(config->alignment_num_alternative_paths == 1u
                && "only the best alignment is used in query");

        aligner = build_aligner(*graph, *config);
    }

    QueryExecutor executor(*config, *anno_graph, aligner.get(), thread_pool);

    // iterate over input files
    for (const auto &file : files) {
        Timer curr_timer;

        executor.query_fasta(file, [](const std::string &result) { std::cout << result; });
        logger->trace("File '{}' was processed in {} sec, total time: {}", file,
                      curr_timer.elapsed(), timer.elapsed());
    }

    return 0;
}

inline std::string query_sequence(size_t id,
                                  const std::string &name,
                                  const std::string &seq,
                                  const AnnotatedDBG &anno_graph,
                                  const Config &config) {
    return QueryExecutor::execute_query(fmt::format_int(id).str() + '\t' + name, seq,
                                        config.count_labels, config.print_signature,
                                        config.suppress_unlabeled, config.num_top_labels,
                                        config.discovery_fraction, config.anno_labels_delimiter,
                                        anno_graph);
}

inline void update_counters(const AnnotatedDBG &anno_graph,
                            HLLCounter &file_kmer_counter,
                            std::vector<HLLCounter> &match_counter,
                            std::mutex &match_counter_mutex,
                            const std::string *seq = nullptr,
                            const std::vector<uint64_t> *index_in_full_graph = nullptr) {
    const auto &graph = anno_graph.get_graph();
    const auto &matrix = anno_graph.get_annotation().get_matrix();

    std::vector<uint64_t> node_hashes;

    if (!seq) {
        // a query graph was passed, so there's no need to process seq
        assert(index_in_full_graph);
        size_t num_rows = matrix.num_rows();
        node_hashes.reserve(num_rows);
        for (size_t i = 0; i < num_rows; ++i) {
            node_hashes.push_back(utils::Hash<DeBruijnGraph::node_index>()(
                index_in_full_graph->at(AnnotatedDBG::anno_to_graph_index(i))
            ));
        }

        {
            std::lock_guard<std::mutex> lock(match_counter_mutex);
            file_kmer_counter.insert(node_hashes.data(), node_hashes.data() + node_hashes.size());
        }

        for (size_t i = 0; i < num_rows; ++i) {
            for (BinaryMatrix::Column j : matrix.get_row(i)) {
                std::lock_guard<std::mutex> lock(match_counter_mutex);
                match_counter[j].insert(node_hashes[i]);
            }
        }

        return;
    }

    node_hashes.reserve(seq->size() - graph.get_k() + 1);
    std::vector<uint64_t> rows;
    rows.reserve(node_hashes.capacity());

    graph.map_to_nodes(*seq, [&](DeBruijnGraph::node_index node) {
        if (node) {
            assert(!index_in_full_graph || index_in_full_graph->at(node));
            rows.push_back(AnnotatedDBG::graph_to_anno_index(node));
            node_hashes.push_back(utils::Hash<DeBruijnGraph::node_index>()(
                index_in_full_graph ? index_in_full_graph->at(node) : node
            ));
        }
    });

    if (node_hashes.empty())
        return;

    {
        std::lock_guard<std::mutex> lock(match_counter_mutex);
        file_kmer_counter.insert(node_hashes.data(), node_hashes.data() + node_hashes.size());
    }

    std::vector<std::vector<uint64_t>> node_hash_map(match_counter.size());
    auto jt = node_hashes.begin();
    for (const auto &row : matrix.get_rows(rows)) {
        assert(jt != node_hashes.end());

        for (uint64_t j : row) {
            node_hash_map.at(j).push_back(*jt);
        }

        ++jt;
    }
    assert(jt == node_hashes.end());

    auto it = node_hash_map.begin();
    for (auto &counter : match_counter) {
        assert(it != node_hash_map.end());
        std::lock_guard<std::mutex> lock(match_counter_mutex);
        counter.insert(it->data(), it->data() + it->size());
        ++it;
    }
    assert(it == node_hash_map.end());
}

inline std::string print_counter(const AnnotatedDBG &anno_graph,
                                 const std::string &file,
                                 const HLLCounter &file_kmer_counter,
                                 const std::vector<HLLCounter> &match_counter,
                                 const std::vector<size_t> &label_counts) {
    assert(match_counter.size() == label_counts.size());
    std::string output;
    output.reserve(1'000);

    output += fmt::format(
        "{}:{}",
        file,
        static_cast<size_t>(std::roundl(file_kmer_counter.estimate_cardinality()))
    );

    const auto &label_encoder = anno_graph.get_annotation().get_label_encoder();
    for (size_t i = 0; i < label_counts.size(); ++i) {
        const auto &counter = match_counter[i];
        size_t estimate = std::min(
            static_cast<size_t>(std::roundl(counter.estimate_cardinality())),
            label_counts[i]
        );

        if (estimate) {
            output += fmt::format("\t<{}>:{}:{}",
                                  label_encoder.decode(i),
                                  estimate,
                                  label_counts[i]);
        }
    }

    output += "\n";

    return output;
}

inline std::string print_coverage(const AnnotatedDBG &anno_graph,
                                  const std::string &file,
                                  const std::vector<size_t> &full_label_counts) {
    std::string output;
    output.reserve(1'000);

    const auto &annotation = anno_graph.get_annotation();

    // TODO: fix k-mer count
    output += fmt::format(
        "{}:{}",
        file,
        annotation.num_objects() / (1 + anno_graph.get_graph().is_canonical_mode())
    );

    const auto &labels = annotation.get_all_labels();
    const auto label_counts = annotation.get_label_counts();
    assert(labels.size() == label_counts.size());
    assert(labels.size() == full_label_counts.size());
    for (size_t i = 0; i < labels.size(); ++i) {
        if (label_counts[i]) {
            output += fmt::format("\t<{}>:{}:{}",
                                  labels[i],
                                  label_counts[i],
                                  full_label_counts[i]);
        }
    }

    output += "\n";

    return output;
}

void QueryExecutor::query_fasta(const string &file,
                                const std::function<void(const std::string &)> &callback) {
    logger->trace("Parsing sequences from file '{}'", file);

    seq_io::FastaParser fasta_parser(file, config_.forward_and_reverse);

    if (config_.fast || config_.get_coverage) {
        // Construct a query graph and query against it
        batched_query_fasta(fasta_parser, callback);
        return;
    }

    // Query sequences independently

    size_t seq_count = 0;

    for (const seq_io::kseq_t &kseq : fasta_parser) {
        thread_pool_.enqueue(
            [&](size_t id, const std::string &name, std::string &seq) {
                if (!aligner_) {
                    callback(query_sequence(id, name, seq, anno_graph_, config_));
                    return;
                }

                // query the alignment matches against the annotator
                auto matches = aligner_->align(seq);

                std::string seq_header
                    = get_alignment_header_and_swap_query(name, &seq, &matches);

                callback(query_sequence(id, seq_header, seq, anno_graph_, config_));
            },
            seq_count++,
            std::string(kseq.name.s),
            std::string(kseq.seq.s)
        );
    }

    // wait while all threads finish processing the current file
    thread_pool_.join();
}

void QueryExecutor
::batched_query_fasta(seq_io::FastaParser &fasta_parser,
                      const std::function<void(const std::string &)> &callback) {
    auto begin = fasta_parser.begin();
    auto end = fasta_parser.end();

    HLLCounter file_kmer_counter(kHLLCounterError);
    std::mutex match_counter_mutex;
    if (config_.get_coverage) {
        for (auto &counter : match_counter_) {
            counter.reset();
        }
    }

    const uint64_t batch_size = config_.query_batch_size_in_bytes;
    seq_io::FastaParser::iterator it;

    size_t seq_count = 0;
    while (begin != end) {
        Timer batch_timer;

        std::vector<std::tuple<size_t, std::string, std::string>> named_alignments;
        uint64_t num_bytes_read = 0;

        StringGenerator generate_batch = [&](auto call_sequence) {
            num_bytes_read = 0;

            if (!aligner_) {
                // basic query regime
                // the query graph is constructed directly from the input sequences
                for (it = begin; it != end && num_bytes_read <= batch_size; ++it) {
                    call_sequence(it->seq.s);
                    num_bytes_read += it->seq.l;
                }
                return;
            }
            // Query with alignment to graph

            // Check if this isn't the first invocation,
            // if the sequences have already been aligned and
            // if the results are stored in |named_alignments|.
            if (named_alignments.size()) {
                for (const auto &[id, name, seq] : named_alignments) {
                    call_sequence(seq);
                }
                return;
            }

            std::mutex sequence_mutex;
            for ( ; begin != end && num_bytes_read <= batch_size; ++begin) {
                thread_pool_.enqueue(
                    [&](size_t id, const std::string &name, std::string &seq) {
                        // Align the sequence, then add the best match
                        // in the graph to the query graph.
                        auto matches = aligner_->align(seq);

                        std::string seq_header
                            = get_alignment_header_and_swap_query(name, &seq, &matches);

                        std::lock_guard<std::mutex> lock(sequence_mutex);

                        named_alignments.emplace_back(id, std::move(seq_header),
                                                      std::move(seq));

                        call_sequence(std::get<2>(named_alignments.back()));
                    },
                    seq_count++,
                    std::string(begin->name.s),
                    std::string(begin->seq.s)
                );

                num_bytes_read += begin->seq.l;
            }
            thread_pool_.join();
        };

        auto [query_graph, index_in_full_graph] = construct_query_graph(
            anno_graph_,
            generate_batch,
            config_.count_labels || config_.get_coverage ? 0 : config_.discovery_fraction,
            get_num_threads()
        );

        logger->trace("Query graph constructed for batch of {} bytes from '{}' in {} sec",
                      num_bytes_read, fasta_parser.get_filename(), batch_timer.elapsed());

        batch_timer.reset();

        if (!aligner_) {
            if (config_.get_coverage) {
                if (it == end) {
                    callback(print_coverage(*query_graph,
                                            fasta_parser.get_filename(),
                                            label_counts_));
                    return;
                }

                update_counters(*query_graph,
                                file_kmer_counter,
                                match_counter_,
                                match_counter_mutex,
                                nullptr,
                                &index_in_full_graph);

                begin = it;
            } else {
                for ( ; begin != it; ++begin) {
                    assert(begin != end);

                    thread_pool_.enqueue(
                        [&](size_t id, const std::string &name, const std::string &seq) {
                            callback(query_sequence(id, name, seq, *query_graph, config_));
                        },
                        seq_count++, std::string(begin->name.s),
                        std::string(begin->seq.s)
                    );
                }
            }
        } else {
            if (config_.get_coverage && begin == end) {
                auto [align_graph, index_in_full_graph] = construct_query_graph(
                    anno_graph_,
                    [&](auto call_sequence) {
                        for (auto&& [id, name, seq] : named_alignments) {
                            call_sequence(std::move(seq));
                        }
                    },
                    0,
                    get_num_threads()
                );

                callback(print_coverage(*align_graph,
                                        fasta_parser.get_filename(),
                                        label_counts_));

                return;
            }

            for (auto&& [id, name, seq] : named_alignments) {
                thread_pool_.enqueue(
                    [&](size_t id, const std::string &name, const std::string &seq) {
                        if (config_.get_coverage) {
                            update_counters(*query_graph,
                                            file_kmer_counter,
                                            match_counter_,
                                            match_counter_mutex,
                                            &seq,
                                            &index_in_full_graph);
                            return;
                        }

                        callback(query_sequence(id, name, seq, *query_graph, config_));
                    },
                    id, std::move(name), std::move(seq)
                );
            }
        }

        thread_pool_.join();
        logger->trace("Batch of {} bytes from '{}' queried in {} sec", num_bytes_read,
                      fasta_parser.get_filename(), batch_timer.elapsed());
    }

    if (config_.get_coverage) {
        callback(print_counter(anno_graph_,
                               fasta_parser.get_filename(),
                               file_kmer_counter,
                               match_counter_,
                               label_counts_));
    }
}

} // namespace cli
} // namespace mtg
