#include "query.hpp"

#include <ips4o.hpp>
#include <tsl/ordered_set.h>
#include <fmt/format.h>

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
#include "cli/align.hpp"


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

typedef typename mtg::graph::DeBruijnGraph::node_index node_index;


QueryExecutor::QueryExecutor(const Config &config,
                             const graph::AnnotatedDBG &anno_graph,
                             std::unique_ptr<graph::align::DBGAlignerConfig>&& aligner_config,
                             ThreadPool &thread_pool)
      : config_(config),
        anno_graph_(anno_graph),
        aligner_config_(std::move(aligner_config)),
        thread_pool_(thread_pool) {
    if (aligner_config_ && aligner_config_->forward_and_reverse_complement)
        throw std::runtime_error("Error: align_both_strands must be off when querying");
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

void call_suffix_match_sequences(const DBGSuccinct &dbg_succ,
                                 const std::string_view contig,
                                 const std::vector<node_index> &nodes_in_full,
                                 const std::function<void(std::string&&,
                                                          node_index)> &callback,
                                 size_t sub_k,
                                 size_t max_num_nodes_per_suffix) {
    assert(sub_k < dbg_succ.get_k());
    assert(nodes_in_full.size() == contig.size() - dbg_succ.get_k() + 1);

    for (size_t prev_match_len = 0, i = 0; i < nodes_in_full.size(); ++i) {
        if (!nodes_in_full[i]) {
            // if prefix[i:i+prev_match_len] was a match on the previous step, then
            // prefix[i+1:i+prev_match_len] of length prev_match_len-1 must be a match on this step
            size_t cur_match_len = prev_match_len ? prev_match_len - 1 : 0;
            // TODO: call first |max_num_nodes_per_suffix| matches
            //       and, if there are too many of them, discard them here
            // TODO: test if this heuristic works and we need to discard large ranges at all
            dbg_succ.call_nodes_with_suffix_matching_longest_prefix(
                std::string_view(&contig[i], dbg_succ.get_k()),
                [&](node_index node, size_t match_len) {
                    assert(match_len >= cur_match_len);
                    cur_match_len = match_len;
                    callback(dbg_succ.get_node_sequence(node), node);
                },
                std::max(sub_k, prev_match_len), // new match must be at least as long as previous
                max_num_nodes_per_suffix
            );
            prev_match_len = cur_match_len;
        } else {
            prev_match_len = dbg_succ.get_k();
        }
    }
}

struct HullPathContext {
    std::string last_kmer;
    node_index last_node;
    size_t depth; // not equal to path.size() if the path is cut off
    size_t fork_count;
};

// Expand the query graph by traversing around its nodes which are forks in the
// full graph. Take at most |max_hull_forks| forks and traverse paths for
// at most |max_hull_depth| steps.
// |continue_traversal| is given a node and the distrance traversed so far and
// returns whether traversal should continue.
template <class ContigCallback>
void call_hull_sequences(const DeBruijnGraph &full_dbg,
                         node_index node,
                         std::string kmer,
                         const ContigCallback &callback,
                         const std::function<bool(std::string_view seq,
                                                  node_index last_node,
                                                  size_t depth,
                                                  size_t fork_count)> &continue_traversal) {
    // DFS from branching points
    kmer.push_back('$');
    std::vector<HullPathContext> paths_to_extend;
    full_dbg.call_outgoing_kmers(node, [&](node_index next_node, char c) {
        kmer.back() = c;
        assert(full_dbg.kmer_to_node(kmer) == next_node);
        if (continue_traversal(kmer, next_node, 1, 0)) {
            paths_to_extend.emplace_back(HullPathContext{
                .last_kmer = kmer,
                .last_node = next_node,
                .depth = 1,
                .fork_count = 0
            });
        } else {
            callback(kmer, std::vector<node_index>{ next_node });
        }
    });

    while (paths_to_extend.size()) {
        HullPathContext hull_path = std::move(paths_to_extend.back());
        paths_to_extend.pop_back();

        std::string &seq = hull_path.last_kmer;
        std::vector<node_index> path = { hull_path.last_node };
        size_t depth = hull_path.depth;
        size_t fork_count = hull_path.fork_count;

        assert(seq.size() == path.size() + full_dbg.get_k() - 1);

        bool extend = true;
        while (extend && full_dbg.has_single_outgoing(path.back())) {
            full_dbg.call_outgoing_kmers(path.back(), [&](auto node, char c) {
                path.push_back(node);
                seq.push_back(c);
            });
            depth++;
            extend = continue_traversal(seq, path.back(), depth, fork_count);
        }

        assert(path.size() == seq.size() - full_dbg.get_k() + 1);
        assert(path == map_sequence_to_nodes(full_dbg, seq));

        callback(seq, path);

        if (!extend)
            continue;

        // a fork or a sink has been reached before the path has reached max depth
        assert(!full_dbg.has_single_outgoing(path.back()));

        node = path.back();
        path.resize(0);
        seq.erase(seq.begin(), seq.end() - full_dbg.get_k() + 1);
        seq.push_back('$');

        // schedule further traversals
        full_dbg.call_outgoing_kmers(node, [&](node_index next_node, char c) {
            seq.back() = c;
            assert(full_dbg.kmer_to_node(seq) == next_node);
            if (continue_traversal(seq, next_node, depth + 1, fork_count + 1)) {
                paths_to_extend.emplace_back(HullPathContext{
                    .last_kmer = seq,
                    .last_node = next_node,
                    .depth = depth + 1,
                    .fork_count = fork_count + 1
                });
            } else {
                callback(seq, std::vector<node_index>{ next_node });
            }
        });
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
void add_to_succinct(BOSSConstructor &constructor, const Contigs &contigs, size_t k) {
    for (const auto &[contig, nodes_in_full] : contigs) {
        size_t begin = 0;
        size_t end;
        do {
            end = std::find(nodes_in_full.begin() + begin, nodes_in_full.end(),
                            DeBruijnGraph::npos) - nodes_in_full.begin();
            if (begin != end) {
                constructor.add_sequence(std::string_view(
                    contig.data() + begin, end - begin + k - 1
                ));
            }
            begin = end + 1;
        } while (end < nodes_in_full.size());
    }
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
 * 3. Extract annotation for the nodes of the query graph and return.
 */
std::unique_ptr<AnnotatedDBG>
construct_query_graph(const AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      size_t num_threads,
                      bool canonical,
                      const Config *config) {
    const auto &full_dbg = anno_graph.get_graph();
    const auto &full_annotation = anno_graph.get_annotation();
    const auto *dbg_succ = dynamic_cast<const DBGSuccinct *>(&full_dbg);

    canonical |= full_dbg.is_canonical_mode();

    size_t sub_k = full_dbg.get_k();
    size_t max_hull_forks = 0;
    size_t max_hull_depth = 0;
    size_t max_num_nodes_per_suffix = 1;
    double max_hull_depth_per_seq_char = 0.0;
    if (config) {
        if (config->alignment_min_seed_length
                && config->alignment_min_seed_length < full_dbg.get_k()) {
            if (dbg_succ) {
                sub_k = config->alignment_min_seed_length;
            } else {
                logger->warn("Matching suffixes of k-mers only supported for DBGSuccinct");
            }
        }

        max_hull_forks = config->max_hull_forks;
        max_hull_depth = config->max_hull_depth;
        max_hull_depth_per_seq_char = config->alignment_max_nodes_per_seq_char;
        max_num_nodes_per_suffix = config->alignment_max_num_seeds_per_locus;
    }

    Timer timer;

    // construct graph storing all k-mers in query
    auto graph_init = std::make_shared<DBGHashOrdered>(full_dbg.get_k());

    size_t max_input_sequence_length = 0;

    if (kPrefilterWithBloom && dbg_succ && sub_k >= full_dbg.get_k()) {
        if (dbg_succ->get_bloom_filter())
            logger->trace(
                    "[Query graph construction] Started indexing k-mers pre-filtered "
                    "with Bloom filter");

        call_sequences([&](const std::string &seq) {
            // TODO: implement add_sequence with filter for all graph representations
            graph_init->add_sequence(
                seq, get_missing_kmer_skipper(dbg_succ->get_bloom_filter(), seq)
            );
            if (max_input_sequence_length < seq.size())
                max_input_sequence_length = seq.size();
        });
    } else {
        call_sequences([&](const std::string &seq) {
            graph_init->add_sequence(seq);
            if (max_input_sequence_length < seq.size())
                max_input_sequence_length = seq.size();
        });
    }

    max_hull_depth = std::min(
        max_hull_depth,
        static_cast<size_t>(max_hull_depth_per_seq_char * max_input_sequence_length)
    );

    logger->trace("[Query graph construction] k-mer indexing took {} sec", timer.elapsed());
    timer.reset();

    logger->trace("[Query graph construction] extracted {} k-mers",
                  graph_init->num_nodes());

    // pull contigs from query graph
    std::vector<std::pair<std::string, std::vector<node_index>>> contigs;
    std::vector<std::pair<std::string, std::vector<node_index>>> rev_comp_contigs;

    std::mutex seq_mutex;
    graph_init->call_sequences([&](auto&&... contig_args) {
        std::lock_guard<std::mutex> lock(seq_mutex);
        contigs.emplace_back(std::move(contig_args)...);
    }, get_num_threads(), canonical);  // pull only primary contigs when building canonical query graph

    if (canonical) {
        rev_comp_contigs.resize(contigs.size());
        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (size_t i = 0; i < contigs.size(); ++i) {
            rev_comp_contigs[i].first = contigs[i].first;
            reverse_complement(rev_comp_contigs[i].first.begin(),
                               rev_comp_contigs[i].first.end());
        }
    }

    logger->trace("[Query graph construction] Contig extraction took {} sec", timer.elapsed());
    timer.reset();

    if (!full_dbg.is_canonical_mode() || sub_k < full_dbg.get_k() || max_hull_forks) {
        // map from nodes in query graph to full graph
        logger->trace("[Query graph construction] Mapping k-mers back to full graph");
        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (size_t i = 0; i < contigs.size(); ++i) {
            contigs[i].second = map_sequence_to_nodes(full_dbg, contigs[i].first);
            if (canonical) {
                // the reverse complement mapping only needs to be mapping if the
                // full graph is not canonical, or if this will be used for hull
                // computation
                rev_comp_contigs[i].second = map_sequence_to_nodes(
                    full_dbg, rev_comp_contigs[i].first
                );
            }
        }
        logger->trace("[Query graph construction] Contigs mapped to graph in {} sec",
                      timer.elapsed());
        timer.reset();
    }

    size_t original_size = contigs.size();

    // add nodes with suffix matches to the query
    if (sub_k < full_dbg.get_k()) {
        assert(dbg_succ);
        logger->trace("[Query graph construction] Adding k-mers with matching "
                      "suffixes of length {}", sub_k);
        timer.reset();

        size_t num_added = 0;
        std::atomic<size_t> num_explored = 0;

        std::vector<std::pair<std::string, std::vector<node_index>>> contig_buffer;
        std::vector<std::pair<std::string, std::vector<node_index>>> rev_comp_contig_buffer;

        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (size_t i = 0; i < contigs.size(); ++i) {
            const auto &[contig, path] = contigs[i];
            std::vector<std::pair<std::string, node_index>> added_nodes;
            std::vector<std::pair<std::string, node_index>> added_nodes_rc;
            auto callback = [&](std::string&& seq, node_index node) {
                assert(node == full_dbg.kmer_to_node(seq));
                added_nodes.emplace_back(std::move(seq), node);
                ++num_explored;
                if (canonical) {
                    added_nodes_rc.emplace_back(added_nodes.back().first,
                                                DeBruijnGraph::npos);
                    reverse_complement(added_nodes_rc.back().first.begin(),
                                       added_nodes_rc.back().first.end());
                    added_nodes_rc.back().second = full_dbg.kmer_to_node(
                        added_nodes_rc.back().first
                    );
                    ++num_explored;
                }
            };
            call_suffix_match_sequences(*dbg_succ, contig, path,
                                        callback, sub_k, max_num_nodes_per_suffix);
            if (canonical) {
                call_suffix_match_sequences(*dbg_succ,
                                            rev_comp_contigs[i].first,
                                            rev_comp_contigs[i].second,
                                            callback, sub_k, max_num_nodes_per_suffix);
            }

            #pragma omp critical
            {
                for (size_t j = 0; j < added_nodes.size(); ++j) {
                    contig_buffer.emplace_back(
                        std::move(added_nodes[j].first),
                        std::vector<node_index>{ added_nodes[j].second }
                    );
                }
                for (size_t j = 0; j < added_nodes_rc.size(); ++j) {
                    rev_comp_contig_buffer.emplace_back(
                        std::move(added_nodes_rc[j].first),
                        std::vector<node_index>{ added_nodes_rc[j].second }
                    );
                }
            }
        }

        assert((canonical && contig_buffer.size() == rev_comp_contig_buffer.size())
                || (!canonical && rev_comp_contig_buffer.empty()));

        for (auto&& pair : contig_buffer) {
            graph_init->add_sequence(pair.first, [&](node_index) { ++num_added; });
            contigs.emplace_back(std::move(pair));
        }
        for (auto&& pair : rev_comp_contig_buffer) {
            graph_init->add_sequence(pair.first, [&](node_index) { ++num_added; });
            rev_comp_contigs.emplace_back(std::move(pair));
        }

        assert((canonical && contigs.size() == rev_comp_contigs.size())
                || (!canonical && rev_comp_contigs.empty()));

        logger->trace("[Query graph construction] Finding {} and adding {} suffix-matching k-mers "
                      "took {} sec", num_explored, num_added, timer.elapsed());
    }

    if (max_hull_forks) {
        logger->trace("[Query graph extension] Computing query graph hull");
        logger->trace("[Query graph expansion] max traversal distance: {}\tmax fork count: {}",
                      max_hull_depth, max_hull_forks);
        timer.reset();

        tsl::hopscotch_map<node_index, uint32_t> distance_traversed_until_node;

        size_t hull_contig_count = 0;
        size_t num_added = 0;

        std::vector<std::pair<std::string, std::vector<node_index>>> contig_buffer;
        std::vector<std::pair<std::string, std::vector<node_index>>> rev_comp_contig_buffer;

        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (size_t i = 0; i < contigs.size(); ++i) {
            const auto &[contig, path] = contigs[i];
            std::vector<std::pair<std::string, std::vector<node_index>>> added_paths;
            std::vector<std::pair<std::string, std::vector<node_index>>> added_paths_rc;
            // TODO: combine these two callbacks into one
            auto callback = [&](const std::string &sequence,
                                const std::vector<node_index> &path) {
                added_paths.emplace_back(sequence, path);
                if (canonical) {
                    added_paths_rc.emplace_back(added_paths.back().first,
                                                std::vector<node_index>{});
                    reverse_complement(added_paths_rc.back().first.begin(),
                                       added_paths_rc.back().first.end());
                    if (!full_dbg.is_canonical_mode()) {
                        // no need to map here because these will be remapped
                        // below
                        added_paths_rc.back().second = map_sequence_to_nodes(
                            full_dbg, added_paths_rc.back().first
                        );
                    }
                }
            };
            auto continue_traversal = [&](std::string_view seq,
                                          node_index last_node,
                                          size_t depth,
                                          size_t fork_count) {
                if (fork_count >= max_hull_forks)
                    return false;

                // if the last node is already in the graph, cut off traversal
                // since this node will be covered in another traversal
                if (graph_init->find(seq.substr(seq.size() - graph_init->get_k()))) {
                    return false;
                }

                // when a node which has already been accessed is visited,
                // only continue traversing if the previous access was in a
                // longer path (i.e., it cut off earlier)
                bool ret_val;
                // TODO: check the number of forks too (shorter paths may have more forks)
                #pragma omp critical
                {
                    auto [it, inserted]
                        = distance_traversed_until_node.emplace(last_node, depth);

                    ret_val = inserted || depth < it->second;

                    if (!inserted && depth < it->second)
                        it.value() = depth;
                }
                return ret_val;
            };

            for (size_t j = 0; j < path.size(); ++j) {
                // if the next starting node is not in the graph, or if it has a single
                // outgoing node which is also in this contig, skip
                if (path[j] && !(j + 1 < path.size()
                                        && path[j + 1]
                                        && full_dbg.has_single_outgoing(path[j]))) {
                    call_hull_sequences(full_dbg, path[j],
                                        contig.substr(j + 1, full_dbg.get_k() - 1),
                                        callback, continue_traversal);
                }
            }

            if (i < rev_comp_contigs.size()) {
                const std::string &rev_contig = rev_comp_contigs[i].first;
                const std::vector<node_index> &rev_path = rev_comp_contigs[i].second;

                for (size_t j = 0; j < rev_path.size(); ++j) {
                    // if the next starting node is not in the graph, or if it has a single
                    // outgoing node which is also in this contig, skip
                    if (rev_contig[j] && !(j + 1 < rev_contig.size()
                                            && rev_contig[j + 1]
                                            && full_dbg.has_single_outgoing(rev_contig[j]))) {
                        call_hull_sequences(full_dbg, rev_path[j],
                                            rev_contig.substr(j + 1, full_dbg.get_k() - 1),
                                            callback, continue_traversal);
                    }
                }
            }

            #pragma omp critical
            {
                hull_contig_count += added_paths.size() + added_paths_rc.size();
                for (auto&& pair : added_paths) {
                    assert(pair.first.size() == pair.second.size() + full_dbg.get_k() - 1);
                    contig_buffer.emplace_back(std::move(pair));
                }
                for (auto&& pair : added_paths_rc) {
                    assert(full_dbg.is_canonical_mode()
                        || pair.first.size() == pair.second.size() + full_dbg.get_k() - 1);
                    rev_comp_contig_buffer.emplace_back(std::move(pair));
                }
            }
        }

        assert((canonical && contig_buffer.size() == rev_comp_contig_buffer.size())
                || (!canonical && rev_comp_contig_buffer.empty()));

        for (auto&& pair : contig_buffer) {
            graph_init->add_sequence(pair.first, [&](node_index) { ++num_added; });
            contigs.emplace_back(std::move(pair));
        }
        for (auto&& pair : rev_comp_contig_buffer) {
            graph_init->add_sequence(pair.first, [&](node_index) { ++num_added; });
            rev_comp_contigs.emplace_back(std::move(pair));
        }

        assert((canonical && contigs.size() == rev_comp_contigs.size())
                || (!canonical && rev_comp_contigs.empty()));

        logger->trace("[Query graph extension] Added {} nodes from {} contigs in {} sec",
                      num_added, hull_contig_count, timer.elapsed());
    }

    if (full_dbg.is_canonical_mode()) {
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
        for (size_t i = 0; i < contigs.size(); ++i) {
            size_t j = 0;
            full_dbg.map_to_nodes(contigs[i].first, [&](node_index n) {
                contigs[i].second[j++] = n;
            });
            rev_comp_contigs[i].second.resize(contigs[i].second.size());
            std::copy(contigs[i].second.begin(), contigs[i].second.end(),
                      rev_comp_contigs[i].second.rbegin());
        }
    }

    logger->trace("[Query graph construction] Intersecting batch graph with full graph");
    timer.reset();
    std::shared_ptr<DeBruijnGraph> graph;

    // restrict nodes to those in the full graph
    if (sub_k < full_dbg.get_k()) {
        BOSSConstructor constructor(full_dbg.get_k() - 1, canonical, 0, "", num_threads);
        add_to_succinct(constructor, contigs, full_dbg.get_k());

        if (canonical && !full_dbg.is_canonical_mode())
            add_to_succinct(constructor, rev_comp_contigs, full_dbg.get_k());

        graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor), canonical);

    } else {
        auto graph_intersection = std::make_shared<DBGHashOrdered>(
            full_dbg.get_k(), canonical
        );

        for (const auto &pair : contigs) {
            const std::string &contig = pair.first;
            const std::vector<node_index> &nodes_in_full = pair.second;
            assert(contig.size() == nodes_in_full.size() + full_dbg.get_k() - 1);
            size_t j = 0;
            graph_intersection->add_sequence(contig, [&]() { return !nodes_in_full[j++]; });
            assert(j == nodes_in_full.size());
        }

        if (canonical && !full_dbg.is_canonical_mode()) {
            for (const auto &pair : rev_comp_contigs) {
                const std::string &contig = pair.first;
                const std::vector<node_index> &nodes_in_full = pair.second;
                assert(contig.size() == nodes_in_full.size() + full_dbg.get_k() - 1);
                size_t j = 0;
                graph_intersection->add_sequence(contig, [&]() { return !nodes_in_full[j++]; });
                assert(j == nodes_in_full.size());
            }
        }

        graph = std::move(graph_intersection);
    }

    logger->trace("[Query graph construction] Intersected in {} sec", timer.elapsed());

    logger->trace("[Query graph construction] Remapping nodes to full graph");
    timer.reset();

    std::vector<node_index> index_in_full_graph(graph->max_index() + 1);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
    for (size_t i = 0; i < contigs.size(); ++i) {
        std::string &contig = contigs[i].first;
        std::vector<node_index> &nodes_in_full = contigs[i].second;
        assert(contig.size() == nodes_in_full.size() + full_dbg.get_k() - 1);
        size_t j = 0;
        if (original_size == contigs.size()) {
            graph->map_to_nodes_sequentially(contig, [&](node_index node) {
                index_in_full_graph[node] = nodes_in_full[j++];
            });
        } else {
            // nodes in the query graph hull may overlap
            graph->map_to_nodes_sequentially(contig, [&](node_index node) {
                __atomic_store_n(&index_in_full_graph[node], nodes_in_full[j++],
                                 __ATOMIC_RELAXED);
            });
        }

        assert(j == nodes_in_full.size());

        if (rev_comp_contigs.size()) {
            std::string &contig = rev_comp_contigs[i].first;
            std::vector<node_index> &nodes_in_full = rev_comp_contigs[i].second;
            size_t j = 0;
            if (original_size == contigs.size()) {
                graph->map_to_nodes_sequentially(contig, [&](node_index node) {
                    index_in_full_graph[node] = nodes_in_full[j++];
                });
            } else {
                // nodes in the query graph hull may overlap
                graph->map_to_nodes_sequentially(contig, [&](node_index node) {
                    __atomic_store_n(&index_in_full_graph[node], nodes_in_full[j++],
                                     __ATOMIC_RELAXED);
                });
            }
            assert(j == nodes_in_full.size());
        }
    }

    if (original_size != contigs.size())
        __atomic_thread_fence(__ATOMIC_ACQUIRE);

    logger->trace("[Query graph construction] Mapping between graphs constructed in {} sec",
                  timer.elapsed());

    contigs = decltype(contigs)();

    assert(!index_in_full_graph.at(0));

    logger->trace("[Query graph construction] Query graph contains {} k-mers",
                  graph->num_nodes());

    // convert to annotation indexes: remove 0 and shift
    size_t num_objects = 0;
    for (size_t i = 1; i < index_in_full_graph.size(); ++i) {
        if (index_in_full_graph[i]) {
            index_in_full_graph[i - 1]
                = AnnotatedDBG::graph_to_anno_index(index_in_full_graph[i]);
            ++num_objects;
        } else {
            index_in_full_graph[i - 1] = -1;  // npos
        }
    }
    index_in_full_graph.pop_back();

    logger->trace("[Query graph construction] Slicing {} rows out of full annotation",
                  num_objects);

    // initialize fast query annotation
    // copy annotations from the full graph to the query graph
    auto annotation = slice_annotation(full_annotation,
                                       index_in_full_graph,
                                       num_threads);

    logger->trace("[Query graph construction] Query annotation with {} labels"
                  " and {} set bits constructed in {} sec",
                  annotation->num_labels(), annotation->num_relations(), timer.elapsed());
    timer.reset();

    // build annotated graph from the query graph and copied annotations
    return std::make_unique<AnnotatedDBG>(graph, std::move(annotation));
}


int query_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase_annotators.size() == 1);

    auto graph = load_critical_dbg(config->infbase);
    auto anno_graph = initialize_annotated_dbg(graph, *config);

    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);

    Timer timer;

    std::unique_ptr<align::DBGAlignerConfig> aligner_config;
    if (config->align_sequences) {
        assert(config->alignment_num_alternative_paths == 1u
                && "only the best alignment is used in query");

        aligner_config.reset(new align::DBGAlignerConfig(
            initialize_aligner_config(*graph, *config)
        ));
    }

    QueryExecutor executor(*config, *anno_graph, std::move(aligner_config), thread_pool);

    // iterate over input files
    for (const auto &file : files) {
        Timer curr_timer;

        executor.query_fasta(file, [](const std::string &result) { std::cout << result; });
        logger->trace("File '{}' was processed in {} sec, total time: {}", file,
                      curr_timer.elapsed(), timer.elapsed());
    }

    return 0;
}

inline void query_sequence(size_t id, std::string name, std::string seq,
                           const AnnotatedDBG &anno_graph,
                           const Config &config,
                           const align::DBGAlignerConfig *aligner_config,
                           const std::function<void(const std::string&)> &callback) {
    if (aligner_config) {
        auto matches = build_aligner(anno_graph.get_graph(), *aligner_config)->align(seq);
        if (matches.size()) {
            auto &match = matches[0];
            // sequence for querying -- the best alignment
            if (match.get_offset()) {
                seq = anno_graph.get_graph()
                                .get_node_sequence(match[0])
                                .substr(0, match.get_offset())
                        + match.get_sequence();
            } else {
                seq = const_cast<std::string&&>(match.get_sequence());
            }

            name = fmt::format(ALIGNED_SEQ_HEADER_FORMAT, name, seq,
                               match.get_score(), match.get_cigar().to_string());

        } else {
            // no alignment was found
            // the original sequence will be queried
            name = fmt::format(ALIGNED_SEQ_HEADER_FORMAT, name, seq,
                               0, fmt::format("{}S", seq.length()));
        }
    }

    callback(QueryExecutor::execute_query(fmt::format_int(id).str() + '\t' + name, seq,
                                          config.count_labels, config.print_signature,
                                          config.suppress_unlabeled, config.num_top_labels,
                                          config.discovery_fraction, config.anno_labels_delimiter,
                                          anno_graph));
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
        thread_pool_.enqueue([&](size_t id, std::string name, std::string seq) {
            query_sequence(id, name, seq, anno_graph_,
                           config_, aligner_config_.get(), callback);
        }, seq_count++, std::string(kseq.name.s), std::string(kseq.seq.s));
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

    while (begin != end) {
        Timer batch_timer;

        uint64_t num_bytes_read = 0;

        StringGenerator generate_batch = [&](auto call_sequence) {
            num_bytes_read = 0;

            // the query graph is constructed directly from the input sequences
            for (it = begin; it != end && num_bytes_read <= batch_size; ++it) {
                call_sequence(it->seq.s);
                num_bytes_read += it->seq.l;
            }
        };

        auto query_graph = construct_query_graph(
            anno_graph_,
            generate_batch,
            get_num_threads(),
            anno_graph_.get_graph().is_canonical_mode() || config_.canonical,
            aligner_config_ ? &config_ : nullptr
        );

        logger->trace("Query graph constructed for batch of {} bytes from '{}' in {} sec",
                      num_bytes_read, fasta_parser.get_filename(), batch_timer.elapsed());

        batch_timer.reset();

        for ( ; begin != it; ++begin) {
            assert(begin != end);

            thread_pool_.enqueue([&](size_t id, const std::string &name, const std::string &seq) {
                query_sequence(id, name, seq, *query_graph,
                               config_, aligner_config_.get(), callback);
            }, seq_count++, std::string(begin->name.s), std::string(begin->seq.s));
        }

        thread_pool_.join();
        logger->trace("Batch of {} bytes from '{}' queried in {} sec", num_bytes_read,
                      fasta_parser.get_filename(), batch_timer.elapsed());
    }
}

} // namespace cli
} // namespace mtg
