#include "query.hpp"

#include <ips4o.hpp>
#include <tsl/ordered_set.h>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/hashers/hash.hpp"
#include "common/utils/template_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
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

using namespace mtg::annot::binmat;
using namespace mtg::graph;

using mtg::common::logger;
using mtg::graph::boss::BOSS;
using mtg::graph::boss::BOSSConstructor;


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

void call_while_linear(const DeBruijnGraph &graph,
                       DeBruijnGraph::node_index node,
                       const std::function<void(DeBruijnGraph::node_index, char)> &callback,
                       const std::function<bool()> &terminate = []() { return false; }) {
    while (!terminate() && graph.has_single_outgoing(node)) {
        graph.call_outgoing_kmers(node, [&](auto next, char c) {
            callback(next, c);
            node = next;
        });
    }
}

void call_suffix_match_sequences(const DBGSuccinct &dbg_succ,
                                 const std::string &contig,
                                 const std::function<void(std::string&&)> &callback,
                                 size_t sub_k) {
    auto nodes_in_full = map_sequence_to_nodes(dbg_succ, contig);
    size_t k = dbg_succ.get_k();
    assert(sub_k < k);

    size_t last_k = 0;
    for (size_t i = 0; i < nodes_in_full.size(); ++i) {
        if (nodes_in_full[i] != DeBruijnGraph::npos) {
            last_k = 0;
            continue;
        }

        if (last_k >= sub_k) {
            i += last_k - sub_k;
            last_k = 0;
            continue;
        }

        dbg_succ.call_nodes_with_suffix(
            std::string_view(contig.data() + i, k),
            [&](auto node, size_t seed_length) {
                last_k = seed_length;
                std::string next_seq;
                next_seq.reserve(contig.size() * 2);
                next_seq += dbg_succ.get_node_sequence(node);
                callback(std::move(next_seq));
            },
            sub_k,
            1 // max num nodes per suffix
        );
    }
}

// Expand the query graph by traversing around its nodes which are forks in the
// full graph. Take at most max_fork_count forks and traverse a linear path for
// at most max_traversal_distance steps.
void call_halo_sequences(const DeBruijnGraph &full_dbg,
                         const DeBruijnGraph &graph_init,
                         const std::function<void(std::string&&)> &callback,
                         size_t max_fork_count,
                         size_t max_traversal_distance,
                         size_t num_threads = 1,
                         size_t sub_k = std::numeric_limits<size_t>::max()) {
    const auto *dbg_succ = dynamic_cast<const DBGSuccinct *>(&full_dbg);
    if (sub_k >= full_dbg.get_k())
        dbg_succ = nullptr;

    if (!max_fork_count && !dbg_succ)
        return;

    if (!max_traversal_distance)
        max_traversal_distance = full_dbg.get_k() * 2;

    logger->trace("[Query graph expansion] max traversal distance: {}", max_traversal_distance);
    logger->trace("[Query graph expansion] max fork count: {}", max_fork_count);

    std::vector<std::string> contigs;
    graph_init.call_sequences([&](const std::string &contig, const auto &) {
        if (dbg_succ) {
            call_suffix_match_sequences(*dbg_succ, contig, [&](std::string&& contig) {
                if (max_fork_count)
                    contigs.emplace_back(std::string(contig));

                std::string rev = contig;
                callback(std::move(contig));

                reverse_complement(rev.begin(), rev.end());

                if (max_fork_count)
                    contigs.emplace_back(std::string(rev));

                callback(std::move(rev));
            }, sub_k);
        }

        if (!max_fork_count)
            return;

        contigs.push_back(contig);

        if (!graph_init.is_canonical_mode()) {
            contigs.push_back(contig);
            reverse_complement(contigs.back().begin(), contigs.back().end());
        }
    });

    if (!max_fork_count)
        return;

    // when a node which has already been accessed is visited, only continue
    // traversing if the previous access was in a longer path (i.e., it cut
    // off earlier)
    tsl::hopscotch_map<DeBruijnGraph::node_index, size_t> distance_traversed_until_node;
    std::mutex map_mutex;

    auto update_node = [&](DeBruijnGraph::node_index node, size_t distance) {
        std::lock_guard<std::mutex> lock(map_mutex);
        auto emplace = distance_traversed_until_node.emplace(node, distance);
        if (!emplace.second) {
            if (emplace.first->second > distance) {
                emplace.first.value() = distance;
            } else {
                return false;
            }
        }

        return true;
    };

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
    for (size_t i = 0; i < contigs.size(); ++i) {
        auto nodes = map_sequence_to_nodes(full_dbg, contigs[i]);

        for (size_t j = 0; j < nodes.size(); ++j) {
            // if the next starting node is not in the graph, or if it has a single
            // outgoing node which is also in this contig, skip
            if (!nodes[j] || !update_node(nodes[j], 0)
                    || (j + 1 < nodes.size() && nodes[j + 1] != DeBruijnGraph::npos
                        && full_dbg.has_single_outgoing(nodes[j]))) {
                continue;
            }

            // DFS from branching points
            std::string unitig(contigs[i], j + 1, full_dbg.get_k() - 1);
            std::vector<std::tuple<std::string,
                                   DeBruijnGraph::node_index,
                                   size_t>> unitig_traversal;
            full_dbg.call_outgoing_kmers(nodes[j], [&](auto next_node, char c) {
                if (update_node(next_node, unitig.size())) {
                    unitig_traversal.emplace_back(unitig + c, next_node, 0);
                    std::get<0>(unitig_traversal.back()).reserve(max_traversal_distance);
                }
            });

            DeBruijnGraph::node_index node;
            size_t depth;
            while (unitig_traversal.size()) {
                std::tie(unitig, node, depth) = std::move(unitig_traversal.back());
                unitig_traversal.pop_back();

                if (depth >= max_fork_count)
                    continue;

                bool continue_traversal = true;
                call_while_linear(full_dbg, node, [&](auto next_node, char c) {
                    node = next_node;
                    if ((continue_traversal = update_node(node, unitig.size())))
                        unitig += c;

                }, [&]() {
                    return !continue_traversal || unitig.size() >= max_traversal_distance;
                });

                callback(std::string(unitig));

                if (node != DeBruijnGraph::npos
                        && continue_traversal
                        && full_dbg.has_multiple_outgoing(node)) {
                    unitig = unitig.substr(unitig.size() - full_dbg.get_k() + 1);
                    full_dbg.call_outgoing_kmers(node, [&](auto next_node, char c) {
                        if (update_node(next_node, unitig.size())) {
                            unitig_traversal.emplace_back(unitig + c,
                                                          next_node,
                                                          depth + 1);
                            std::get<0>(unitig_traversal.back()).reserve(
                                max_traversal_distance
                            );
                        }
                    });
                }
            }
        }
    }
}

/**
 * @brief      Construct annotation submatrix with a subset of rows extracted
 *             from the full annotation matrix
 *
 * @param[in]  full_annotation  The full annotation matrix.
 * @param[in]  index_in_full    Indexes of rows in the full annotation matrix.
 *                              index_in_full[i] = -1 means that the i-th row in
 *                              the submatrix is empty.
 * @param[in]  num_threads      The number of threads used.
 *
 * @return     Annotation submatrix in the UniqueRowAnnotator representation
 */
std::unique_ptr<annot::UniqueRowAnnotator>
slice_annotation(const AnnotatedDBG::Annotator &full_annotation,
                 const std::vector<uint64_t> &index_in_full,
                 size_t num_threads) {
    const uint64_t npos = -1;

    if (auto *rb = dynamic_cast<const RainbowMatrix *>(&full_annotation.get_matrix())) {
        // shortcut construction for Rainbow<> annotation
        std::vector<uint64_t> row_indexes;
        row_indexes.reserve(index_in_full.size());
        for (uint64_t i : index_in_full) {
            if (i != npos) {
                row_indexes.push_back(i);
            } else {
                row_indexes.push_back(0);
            }
        }

        // get unique rows and set pointers to them in |row_indexes|
        auto unique_rows = rb->get_rows(&row_indexes, num_threads);

        if (unique_rows.size() >= std::numeric_limits<uint32_t>::max()) {
            throw std::runtime_error("There must be less than 2^32 unique rows."
                                     " Reduce the query batch size.");
        }

        // if the 0-th row is not empty, we must insert an empty unique row
        // and reassign those indexes pointing to npos in |index_in_full|.
        if (rb->get_row(0).size()) {
            logger->trace("Add empty row");
            unique_rows.emplace_back();
            for (uint64_t i = 0; i < index_in_full.size(); ++i) {
                if (index_in_full[i] == npos) {
                    row_indexes[i] = unique_rows.size() - 1;
                }
            }
        }

        // copy annotations from the full graph to the query graph
        return std::make_unique<annot::UniqueRowAnnotator>(
            std::make_unique<UniqueRowBinmat>(std::move(unique_rows),
                                              std::vector<uint32_t>(row_indexes.begin(),
                                                                    row_indexes.end()),
                                              full_annotation.num_labels()),
            full_annotation.get_label_encoder()
        );
    }

    std::vector<std::pair<uint64_t, uint64_t>> from_full_to_small;

    for (uint64_t i = 0; i < index_in_full.size(); ++i) {
        if (index_in_full[i] != npos)
            from_full_to_small.emplace_back(index_in_full[i], i);
    }

    ips4o::parallel::sort(from_full_to_small.begin(), from_full_to_small.end(),
                          utils::LessFirst(), num_threads);

    using RowSet = tsl::ordered_set<BinaryMatrix::SetBitPositions,
                                    utils::VectorHash,
                                    std::equal_to<BinaryMatrix::SetBitPositions>,
                                    std::allocator<BinaryMatrix::SetBitPositions>,
                                    std::vector<BinaryMatrix::SetBitPositions>,
                                    uint32_t>;
    RowSet unique_rows { BinaryMatrix::SetBitPositions() };
    std::vector<uint32_t> row_rank(index_in_full.size(), 0);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint64_t batch_begin = 0;
                        batch_begin < from_full_to_small.size();
                                        batch_begin += kRowBatchSize) {
        const uint64_t batch_end
            = std::min(batch_begin + kRowBatchSize,
                       static_cast<uint64_t>(from_full_to_small.size()));

        std::vector<uint64_t> row_indexes;
        row_indexes.reserve(batch_end - batch_begin);
        for (uint64_t i = batch_begin; i < batch_end; ++i) {
            assert(from_full_to_small[i].first < full_annotation.num_objects());

            row_indexes.push_back(from_full_to_small[i].first);
        }

        auto rows = full_annotation.get_matrix().get_rows(row_indexes);

        assert(rows.size() == batch_end - batch_begin);

        #pragma omp critical
        {
            for (uint64_t i = batch_begin; i < batch_end; ++i) {
                const auto &row = rows[i - batch_begin];
                auto it = unique_rows.emplace(row).first;
                row_rank[from_full_to_small[i].second] = it - unique_rows.begin();
                if (unique_rows.size() == std::numeric_limits<uint32_t>::max())
                    throw std::runtime_error("There must be less than 2^32 unique rows."
                                             " Reduce the query batch size.");
            }
        }
    }

    auto &annotation_rows = const_cast<std::vector<BinaryMatrix::SetBitPositions>&>(
        unique_rows.values_container()
    );

    // copy annotations from the full graph to the query graph
    return std::make_unique<annot::UniqueRowAnnotator>(
        std::make_unique<UniqueRowBinmat>(std::move(annotation_rows),
                                          std::move(row_rank),
                                          full_annotation.num_labels()),
        full_annotation.get_label_encoder()
    );
}

template <class Contigs>
std::shared_ptr<DBGSuccinct> convert_to_succinct(const DeBruijnGraph &full_dbg,
                                                 const Contigs &contigs,
                                                 bool canonical = false,
                                                 size_t num_threads = 1) {
    BOSSConstructor constructor(full_dbg.get_k() - 1, canonical, 0, "", num_threads);
    for (size_t i = 0; i < contigs.size(); ++i) {
        const std::string &contig = contigs[i].first;
        const auto &nodes_in_full = contigs[i].second;
        constructor.add_sequences([&](const CallString &callback) {
            auto it = nodes_in_full.begin();
            while ((it = std::find_if(it, nodes_in_full.end(), [&](auto i) { return i; }))
                    < nodes_in_full.end()) {
                auto next = std::find(it, nodes_in_full.end(), 0);
                assert(full_dbg.find(std::string_view(
                    contig.data() + (it - nodes_in_full.begin()),
                    next - it + full_dbg.get_k() - 1
                )));
                callback(std::string(contig.data() + (it - nodes_in_full.begin()),
                                     next - it + full_dbg.get_k() - 1));
                it = next;
            }
        });
    }

    return std::make_shared<DBGSuccinct>(new BOSS(&constructor), canonical);
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
std::unique_ptr<AnnotatedDBG>
construct_query_graph(const AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      double discovery_fraction,
                      size_t num_threads,
                      bool canonical,
                      size_t sub_k,
                      size_t max_fork_count,
                      size_t max_traversal_distance) {
    const auto &full_dbg = anno_graph.get_graph();
    const auto &full_annotation = anno_graph.get_annotation();

    canonical |= full_dbg.is_canonical_mode();

    Timer timer;

    // construct graph storing all k-mers in query
    auto graph_init = std::make_shared<DBGHashOrdered>(full_dbg.get_k(), false);

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct *>(&full_dbg);
    size_t max_sequence_length = 0;
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
            max_sequence_length = std::max(max_sequence_length, sequence.size());
        });
    } else {
        call_sequences([&](const std::string &sequence) {
            graph_init->add_sequence(sequence);
            max_sequence_length = std::max(max_sequence_length, sequence.size());
        });
    }

    logger->trace("[Query graph construction] k-mer indexing took {} sec", timer.elapsed());
    timer.reset();

    logger->trace("[Query graph construction] extracted {} k-mers",
                  graph_init->num_nodes());

    if (max_fork_count || sub_k < full_dbg.get_k()) {
        // extend graph_init with halo from full graph
        std::vector<std::string> new_sequences;
        std::mutex seq_mutex;
        call_halo_sequences(full_dbg, *graph_init,
            [&](std::string&& unitig) {
                std::lock_guard<std::mutex> lock(seq_mutex);
                new_sequences.emplace_back(std::move(unitig));
            },
            max_fork_count,
            std::min(max_traversal_distance, max_sequence_length),
            num_threads,
            sub_k
        );

        for (auto&& seq : new_sequences) {
            graph_init->add_sequence(std::move(seq));
        }

        logger->trace("[Query graph construction] graph expansion took {} sec", timer.elapsed());
        timer.reset();

        logger->trace("[Query graph construction] expanded to {} k-mers",
                      graph_init->num_nodes());
    }

    std::shared_ptr<DeBruijnGraph> graph = std::move(graph_init);

    // pull contigs from query graph
    std::vector<std::pair<std::string, std::vector<DeBruijnGraph::node_index>>> contigs;
    std::mutex seq_mutex;
    graph->call_sequences([&](auto&&... contig_args) {
                              std::lock_guard<std::mutex> lock(seq_mutex);
                              contigs.emplace_back(std::move(contig_args)...);
                          },
                          get_num_threads(),
                          canonical);  // pull only primary contigs when building canonical query graph

    logger->trace("[Query graph construction] Contig extraction took {} sec", timer.elapsed());
    timer.reset();

    // map contigs onto the full graph
    std::vector<uint64_t> index_in_full_graph;

    if (canonical) {
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
        for (size_t i = 0; i < contigs.size(); ++i) {
            std::string &contig = contigs[i].first;
            auto &nodes_in_full = contigs[i].second;

            if (full_dbg.is_canonical_mode()) {
                size_t j = 0;
                full_dbg.map_to_nodes(contig,
                    [&](auto node_in_full) { nodes_in_full[j++] = node_in_full; }
                );
            } else {
                size_t j = 0;
                // TODO: if a k-mer is found, don't search its reverse-complement
                // TODO: add `primary/canonical` mode to DBGSuccinct?
                full_dbg.map_to_nodes_sequentially(contig,
                    [&](auto node_in_full) { nodes_in_full[j++] = node_in_full; }
                );
                reverse_complement(contig.begin(), contig.end());
                full_dbg.map_to_nodes_sequentially(contig,
                    [&](auto node_in_full) {
                        --j;
                        if (node_in_full)
                            nodes_in_full[j] = node_in_full;
                    }
                );
                reverse_complement(contig.begin(), contig.end());
                assert(j == 0);
            }
        }

        logger->trace("[Query graph construction] Contigs mapped to graph in {} sec", timer.elapsed());
        timer.reset();

        // construct canonical graph storing all k-mers found in the full graph
        if (sub_k >= full_dbg.get_k()) {
            graph_init = std::make_shared<DBGHashOrdered>(full_dbg.get_k(), true);

            for (size_t i = 0; i < contigs.size(); ++i) {
                const std::string &contig = contigs[i].first;
                const auto &nodes_in_full = contigs[i].second;
                size_t j = 0;
                graph_init->add_sequence(contig, [&]() { return nodes_in_full[j++] == 0; });
            }

            graph = std::move(graph_init);
        } else {
            graph = convert_to_succinct(full_dbg, contigs, true, num_threads);
        }

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
        if (sub_k >= full_dbg.get_k()) {
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
        } else {
            graph = convert_to_succinct(full_dbg, contigs, false, num_threads);
            contigs = decltype(contigs)();
            index_in_full_graph.assign(graph->max_index() + 1, 0);

            graph->call_sequences([&](const std::string &contig, const auto &path) {
                size_t j = 0;
                full_dbg.map_to_nodes(contig, [&](auto node_in_full) {
                    index_in_full_graph[path[j++]] = node_in_full;
                });
                assert(j == path.size());
            }, get_num_threads());
        }

        logger->trace("[Query graph construction] Contigs mapped to graph in {} sec", timer.elapsed());
        timer.reset();
    }

    contigs = decltype(contigs)();

    assert(!index_in_full_graph.at(0));

    if (discovery_fraction > 0 && sub_k >= full_dbg.get_k()) {
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

    logger->trace("[Query graph construction] Query graph contains {} k-mers",
                  graph->num_nodes());

    // convert to annotation indexes: remove 0 and shift
    for (size_t i = 1; i < index_in_full_graph.size(); ++i) {
        if (index_in_full_graph[i]) {
            index_in_full_graph[i - 1]
                = AnnotatedDBG::graph_to_anno_index(index_in_full_graph[i]);
        } else {
            index_in_full_graph[i - 1] = -1;  // npos
        }
    }
    index_in_full_graph.pop_back();

    // initialize fast query annotation
    // copy annotations from the full graph to the query graph
    auto annotation = slice_annotation(full_annotation,
                                       index_in_full_graph,
                                       num_threads);

    logger->trace("[Query graph construction] Query annotation constructed in {} sec",
                  timer.elapsed());
    timer.reset();

    // build annotated graph from the query graph and copied annotations
    return std::make_unique<AnnotatedDBG>(graph, std::move(annotation));
}


std::string get_alignment_header_and_swap_query(const std::string &name,
                                                std::string *query_seq,
                                                align::QueryAlignment<> *matches) {
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

    std::unique_ptr<align::IDBGAligner> aligner;
    if (config->align_sequences) {
        assert(config->alignment_num_alternative_paths == 1u
                && "only the best alignment is used in query");

        aligner = build_aligner(*graph, *config);

        // the fwd_and_reverse argument in the aligner config returns the best of
        // the forward and reverse complement alignments, rather than both.
        // so, we want to prevent it from doing this
        auto &aligner_config = const_cast<align::DBGAlignerConfig&>(aligner->get_config());
        aligner_config.forward_and_reverse_complement = false;
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

void QueryExecutor::query_fasta(const string &file,
                                const std::function<void(const std::string &)> &callback) {
    logger->trace("Parsing sequences from file '{}'", file);

    seq_io::FastaParser fasta_parser(file, config_.forward_and_reverse);

    if (config_.fast) {
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

    const uint64_t batch_size = config_.query_batch_size_in_bytes;
    seq_io::FastaParser::iterator it;

    size_t seq_count = 0;
    size_t sub_k = aligner_
        ? aligner_->get_config().min_seed_length
        : std::numeric_limits<size_t>::max();

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

        auto query_graph = construct_query_graph(
            anno_graph_,
            generate_batch,
            config_.count_labels || sub_k < anno_graph_.get_graph().get_k()
                ? 0
                : config_.discovery_fraction,
            get_num_threads(),
            anno_graph_.get_graph().is_canonical_mode() || config_.canonical,
            sub_k,
            config_.max_fork_count,
            config_.max_traversal_distance
        );

        logger->trace("Query graph constructed for batch of {} bytes from '{}' in {} sec",
                      num_bytes_read, fasta_parser.get_filename(), batch_timer.elapsed());

        batch_timer.reset();

        if (!aligner_) {
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
        } else {
            for (auto&& [id, name, seq] : named_alignments) {
                thread_pool_.enqueue(
                    [&](size_t id, const std::string &name, const std::string &seq) {
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
}

} // namespace cli
} // namespace mtg
