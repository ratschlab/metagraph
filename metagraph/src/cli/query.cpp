#include "query.hpp"

#include <ips4o.hpp>
#include <tsl/ordered_set.h>
#include <mutex>

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
using namespace mtg::annot::matrix;
using namespace mtg::graph;

using mtg::common::logger;
using mtg::graph::boss::BOSS;
using mtg::graph::boss::BOSSConstructor;

typedef typename mtg::graph::DeBruijnGraph::node_index node_index;

// JSON Field Keys
const std::string SEQ_DESCRIPTION_JSON_FIELD = "seq_description";
const std::string KMER_COUNT_FIELD = "kmer_count";
const std::string LABEL_COUNT_FIELD = "label_count";
const std::string KMER_COORDINATE_FIELD = "kmer_coords";
const std::string SIGNATURE_FIELD = "signature";
const std::string KMER_COUNT_QUANTILE_FIELD = "kmer_count_quantile";
const std::string SCORE_JSON_FIELD = "score";
const std::string SEQUENCE_JSON_FIELD = "sequence";
const std::string ALIGNMENT_JSON_FIELD = "alignments";
const std::string CIGAR_JSON_FIELD = "cigar";


/**
 * Given a vector of kmer label matched coordinates, collapse continuous ranges of coordinates
 * to start-end tuples.
 *
 * @param coords    the vector of tuples as stored in a SeqSearchResult
 *                  originally returned from AnnotatedDBG::get_kmer_coordinates
 * @return vector of 'begin-end' range string representations
 */
std::vector<std::string> get_collapsed_coord_ranges(
        const std::vector<SmallVector<uint64_t>> &tuples) {
    // Build output
    std::vector<std::string> ranges;

    // Keep track of the beginning of the range (just use a sentinel)
    uint64_t begin = 0;
    uint64_t current = 0;
    bool in_range = false;

    for (const auto &coords : tuples) {
        if (coords.empty()) continue;

        // TODO: Look for example where row_tuples are multiple values?
        uint64_t coord = coords[0];

        if (in_range) {
            // Check if we should continue the range
            if (coord == current + 1) {
                current = coord;
            } else {
                // End the range and output
                if (begin == current) {
                    ranges.push_back(fmt::format(":{}", begin));
                } else {
                    ranges.push_back(fmt::format(":{}-{}", begin, current));
                }

                // Start a new range
                begin = coord;
                current = coord;
            }
        } else {
            // Start a new range
            begin = coord;
            current = coord;
            in_range = true;
        }
    }

    // Add the last range
    if (begin == current) {
        ranges.push_back(fmt::format("{}", begin));
    } else {
        ranges.push_back(fmt::format("{}-{}", begin, current));
    }

    return ranges;
}


/**
 * Convert string values into proper JSON values with types (i.e., 'nan' -> null,
 * strings representing numbers -> numbers).
 *
 * @param v value string
 * @return  Json::Value representing value string v
 */
Json::Value adjust_for_types(const std::string &v) {
    if (v == "nan")
        return Json::nullValue;

    try {
        return Json::Value(std::stof(v));
    } catch(...) {}

    return Json::Value(v);
}


/**
 * Given a label string in the format '<sample_name>(;<property_name>=<property_value>)*'
 * return a JSON representation of the label with its properties.
 *
 * @param label     label as a string
 * @return JSON object in the format of
 *  {
 *      "sample": <sample_name>,
 *      "properties": {
 *          <property_name>: <property_value>,
 *          ...
 *      }
 *  }
 */
Json::Value get_label_as_json(const std::string &label) {
    Json::Value label_root;

    // Split by ';' TODO: check this and see why substring was necessary
    std::vector<std::string> label_parts = utils::split_string(label, ";");

    // First field is always the sample
    label_root["sample"] = label_parts[0];

    // Fill properties if existant
    Json::Value properties = Json::objectValue;
    for (auto lit = ++label_parts.begin(); lit != label_parts.end(); ++lit) {
        std::vector<std::string> key_value = utils::split_string(*lit, "=");
        properties[key_value[0]] = adjust_for_types(key_value[1]);
    }

    if (properties.size() > 0) {
        label_root["properties"] = properties;
    }

    return label_root;
}


Json::Value SeqSearchResult::to_json(bool count_kmers,
                                     bool expand_coords,
                                     const graph::AnnotatedDBG &anno_graph) const {
    Json::Value root;

    // Add seq information
    root[SEQ_DESCRIPTION_JSON_FIELD] = sequence.name;

    // Add alignment information if there is any
    if (sequence.alignment.has_value()) {
        const auto &alignment = sequence.alignment.value();

        // Recover alignment sequence sequence string
        root[SEQUENCE_JSON_FIELD] = Json::Value(sequence.sequence);

        // Alignment metrics
        root[SCORE_JSON_FIELD] = Json::Value(alignment.score);
        root[CIGAR_JSON_FIELD] = Json::Value(alignment.cigar);
    }

    // Add discovered labels and extra results
    root["results"] = Json::Value(Json::arrayValue);

    // Different action depending on the result type
    if (std::holds_alternative<label_vec>(result)) {
        // Standard labels only
        for (const auto &label : std::get<label_vec>(result)) {
            root["results"].append(get_label_as_json(label));
        }
    } else if (std::holds_alternative<label_count_vec>(result)) {
        // Labels with count data
        for (const auto &[label, count] : std::get<label_count_vec>(result)) {
            Json::Value label_obj = get_label_as_json(label);

            std::string count_key = (count_kmers) ? KMER_COUNT_FIELD : LABEL_COUNT_FIELD;
            label_obj[count_key] = Json::Value((Json::Int64) count);

            root["results"].append(label_obj);
        }
    } else if (std::holds_alternative<label_sig_vec>(result)) {
        // Count signatures
        for (const auto &[label, kmer_presence_mask] : std::get<label_sig_vec>(result)) {
            Json::Value label_obj = get_label_as_json(label);

            Json::Value sig_obj = Json::objectValue;
            sig_obj["num_present"] = Json::Value(sdsl::util::cnt_one_bits(kmer_presence_mask));
            sig_obj["presence_mask"] = Json::Value(sdsl::util::to_string(kmer_presence_mask));
            sig_obj["score"] =
                    Json::Value(anno_graph.score_kmer_presence_mask(kmer_presence_mask));

            label_obj[SIGNATURE_FIELD] = sig_obj;
            root["results"].append(label_obj);
        }
    } else if (std::holds_alternative<label_quantile_vec>(result)) {
        // Count quantiles
        for (const auto &[label, quantiles] : std::get<label_quantile_vec>(result)) {
            Json::Value label_obj = get_label_as_json(label);
            Json::Value quantile_array = Json::arrayValue;

            for (size_t quantile : quantiles) {
                quantile_array.append(Json::Value((Json::Int64) quantile));
            }

            label_obj[KMER_COUNT_QUANTILE_FIELD] = quantile_array;
            root["results"].append(label_obj);
        }
    } else if (std::holds_alternative<label_coord_vec>(result)) {
        // Kmer coordinates
        for (const auto &[label, tuples] : std::get<label_coord_vec>(result)) {
            Json::Value label_obj = get_label_as_json(label);
            Json::Value coord_array = Json::arrayValue;

            // Each tuple is represented as a folly::SmallVector<uint64_t> instance
            if (expand_coords) {
                for (const auto &coords : tuples) {
                    coord_array.append(Json::Value(fmt::format("{}", fmt::join(coords, ","))));
                }
            } else {
                for (const auto &range : get_collapsed_coord_ranges(tuples)) {
                    coord_array.append(range);
                }
            }


            label_obj[KMER_COORDINATE_FIELD] = coord_array;
            root["results"].append(label_obj);
        }
    }

    return root;
}


std::string SeqSearchResult::to_string(const std::string delimiter,
                                       bool suppress_unlabeled,
                                       bool expand_coords,
                                       const graph::AnnotatedDBG &anno_graph) const {
    // Return the size of the result type vector using std::visit (ducktyping variant)
    if (!std::visit([] (auto&& vec) { return vec.size(); }, result) && suppress_unlabeled) {
        // If empty and unlabeled sequences are supressed then return empty string
        return "";
    }

    // Get modified name if sequence has an alignment result
    std::string mod_seq_name = sequence.name;

    if (sequence.alignment.has_value()) {
        const QuerySequence::Alignment &alignment = sequence.alignment.value();
        mod_seq_name = fmt::format(ALIGNED_SEQ_HEADER_FORMAT, sequence.name, sequence.sequence,
                                   alignment.score, alignment.cigar);
    }

    // Print sequence ID and name
    std::string output = fmt::format("{}\t{}", sequence.id, mod_seq_name);

    // ... with diff result output depending on the type of result being stored.
    if (std::holds_alternative<label_vec>(result)) {
        // Standard labels only
        output += fmt::format("\t{}",
                              utils::join_strings(std::get<label_vec>(result),
                                                  delimiter));
    } else if (std::holds_alternative<label_count_vec>(result)) {
        // Labels with count data
        for (const auto &[label, count] : std::get<label_count_vec>(result)) {
            output += fmt::format("\t<{}>:{}", label, count);
        }
    } else if (std::holds_alternative<label_sig_vec>(result)) {
        // Count signatures
        for (const auto &[label, kmer_presence_mask] : std::get<label_sig_vec>(result)) {
            output += fmt::format("\t<{}>:{}:{}:{}", label,
                                  sdsl::util::cnt_one_bits(kmer_presence_mask),
                                  sdsl::util::to_string(kmer_presence_mask),
                                  anno_graph.score_kmer_presence_mask(kmer_presence_mask));
        }
    } else if (std::holds_alternative<label_quantile_vec>(result)) {
        // Count quantiles
        for (const auto &[label, quantiles] : std::get<label_quantile_vec>(result)) {
            output += fmt::format("\t<{}>:{}", label, fmt::join(quantiles, ":"));
        }
    } else if (std::holds_alternative<label_coord_vec>(result)) {
        // Kmer coordinates
        for (const auto &[label, tuples] : std::get<label_coord_vec>(result)) {
            output += "\t<" + label + ">";

            if (expand_coords) {
                for (const auto &coords : tuples) {
                    output += fmt::format(":{}", fmt::join(coords, ","));
                }
            } else {
                output += fmt::format(":{}",
                                      fmt::join(get_collapsed_coord_ranges(tuples), ":"));
            }
        }
    }

    return output;
}


SeqSearchResult QueryExecutor::execute_query(QuerySequence&& sequence,
                                             bool count_labels,
                                             bool print_signature,
                                             size_t num_top_labels,
                                             double discovery_fraction,
                                             double presence_fraction,
                                             const graph::AnnotatedDBG &anno_graph,
                                             bool with_kmer_counts,
                                             const std::vector<double> &count_quantiles,
                                             bool query_coords) {
    // Perform a different action depending on the type (specified by config flags)
    SeqSearchResult::result_type result;

    if (print_signature) {
        // Get labels with presence/absence signatures
        result = anno_graph.get_top_label_signatures(sequence.sequence,
                                                     num_top_labels,
                                                     discovery_fraction,
                                                     presence_fraction);
    } else if (query_coords) {
        // Get labels with query coordinates
        result = anno_graph.get_kmer_coordinates(sequence.sequence,
                                                 num_top_labels,
                                                 discovery_fraction,
                                                 presence_fraction);
    } else if (count_quantiles.size()) {
        // Get labels with count quantiles
        result = anno_graph.get_label_count_quantiles(sequence.sequence,
                                                      num_top_labels,
                                                      discovery_fraction,
                                                      presence_fraction,
                                                      count_quantiles);
    } else if (count_labels || with_kmer_counts) {
        // Get labels with label counts / kmer counts
        result = anno_graph.get_top_labels(sequence.sequence,
                                           num_top_labels,
                                           discovery_fraction,
                                           presence_fraction,
                                           with_kmer_counts);
    } else {
        // Default, just get labels
        result = anno_graph.get_labels(sequence.sequence,
                                       discovery_fraction,
                                       presence_fraction);
    }

    // Clear sequence if we did not align to save space
    if (!sequence.alignment.has_value()) {
        sequence.sequence = "";
    }

    // Create seq search instance and move sequence into it
    return SeqSearchResult(std::move(sequence), std::move(result));
}

void call_suffix_match_sequences(const DBGSuccinct &dbg_succ,
                                 std::string_view contig,
                                 const std::vector<node_index> &nodes_in_full,
                                 const std::function<void(std::string&&,
                                                          node_index)> &callback,
                                 size_t sub_k,
                                 size_t max_num_nodes_per_suffix) {
    assert(sub_k < dbg_succ.get_k());
    assert(nodes_in_full.size() == contig.length() - dbg_succ.get_k() + 1);

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
// full graph.
// |continue_traversal| is given a node and the distrance traversed so far and
// returns whether traversal should continue.
template <class ContigCallback>
void call_hull_sequences(const DeBruijnGraph &full_dbg,
                         std::string kmer,
                         const ContigCallback &callback,
                         const std::function<bool(std::string_view seq,
                                                  node_index last_node,
                                                  size_t depth,
                                                  size_t fork_count)> &continue_traversal) {
    // DFS from branching points
    node_index node = full_dbg.kmer_to_node(kmer);
    if (!node)
        return;

    kmer.erase(kmer.begin());
    kmer.push_back('$');
    std::vector<HullPathContext> paths_to_extend;
    full_dbg.call_outgoing_kmers(node, [&](node_index next_node, char c) {
        if (c == '$')
            return;

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

        assert(path.size() == seq.length() - full_dbg.get_k() + 1);

        bool extend = true;
        while (extend && full_dbg.has_single_outgoing(path.back())) {
            full_dbg.call_outgoing_kmers(path.back(), [&](auto node, char c) {
                if (c == '$')
                    return;

                path.push_back(node);
                seq.push_back(c);
            });
            depth++;
            extend = continue_traversal(seq, path.back(), depth, fork_count);
        }

        assert(path.size() == seq.length() - full_dbg.get_k() + 1);
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
            if (c == '$')
                return;

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

template <typename T>
annot::LabelEncoder<> reencode_labels(const annot::LabelEncoder<> &encoder,
                                      std::vector<T> *rows) {
    assert(rows);
    annot::LabelEncoder<std::string> new_encoder;
    tsl::hopscotch_map<size_t, size_t> old_to_new;
    for (auto &row : *rows) {
        for (auto &v : row) {
            auto &j = utils::get_first(v);
            auto [it, inserted] = old_to_new.emplace(j, new_encoder.size());
            if (inserted)
                new_encoder.insert_and_encode(encoder.decode(j));

            assert(encoder.decode(j) == new_encoder.decode(it->second));
            j = it->second;
        }
    }
    return new_encoder;
}

/**
 * @brief      Construct annotation submatrix with a subset of rows extracted
 *             from the full annotation matrix
 *
 * @param[in]  full_annotation  The full annotation matrix.
 * @param[in]  num_rows         The number of rows in the target submatrix.
 * @param[in]  full_to_small    The mapping between the rows in the full matrix
 *                              and its submatrix.
 * @param[in]  num_threads      The number of threads used.
 *
 * @return     Annotation submatrix
 */
std::unique_ptr<AnnotatedDBG::Annotator>
slice_annotation(const AnnotatedDBG::Annotator &full_annotation,
                 uint64_t num_rows,
                 std::vector<std::pair<uint64_t, uint64_t>>&& full_to_small,
                 size_t num_threads) {
    if (auto *rb = dynamic_cast<const RainbowMatrix *>(&full_annotation.get_matrix())) {
        // shortcut construction for Rainbow<> annotation
        std::vector<uint64_t> row_indexes(full_to_small.size());
        for (size_t i = 0; i < full_to_small.size(); ++i) {
            row_indexes[i] = full_to_small[i].first;
        }

        // get unique rows and set pointers to them in |row_indexes|
        auto unique_rows = rb->get_rows(&row_indexes, num_threads);

        if (unique_rows.size() >= std::numeric_limits<uint32_t>::max()) {
            throw std::runtime_error("There must be less than 2^32 unique rows."
                                     " Reduce the query batch size.");
        }

        // insert one empty row for representing unmatched rows
        unique_rows.emplace_back();
        std::vector<uint32_t> row_ids(num_rows, unique_rows.size() - 1);
        for (size_t i = 0; i < row_indexes.size(); ++i) {
            row_ids[full_to_small[i].second] = row_indexes[i];
        }

        auto label_encoder = reencode_labels(full_annotation.get_label_encoder(), &unique_rows);

        // copy annotations from the full graph to the query graph
        return std::make_unique<annot::UniqueRowAnnotator>(
            std::make_unique<UniqueRowBinmat>(std::move(unique_rows),
                                              std::move(row_ids),
                                              label_encoder.size()),
            std::move(label_encoder)
        );
    }

    // don't break the topological order for row-diff annotation
    if (!dynamic_cast<const IRowDiff *>(&full_annotation.get_matrix())) {
        ips4o::parallel::sort(full_to_small.begin(), full_to_small.end(),
                              utils::LessFirst(), num_threads);
    }

    if (const auto *mat = dynamic_cast<const IntMatrix *>(&full_annotation.get_matrix())) {
        std::vector<uint64_t> row_indexes;
        row_indexes.reserve(full_to_small.size());
        for (const auto &[in_full, _] : full_to_small) {
            assert(in_full < full_annotation.num_objects());
            row_indexes.push_back(in_full);
        }

        auto slice = mat->get_row_values(row_indexes);

        auto label_encoder = reencode_labels(full_annotation.get_label_encoder(), &slice);

        Vector<CSRMatrix::RowValues> rows(num_rows);

        for (uint64_t i = 0; i < slice.size(); ++i) {
            rows[full_to_small[i].second] = std::move(slice[i]);
        }

        // copy annotations from the full graph to the query graph
        return std::make_unique<annot::IntRowAnnotator>(
            std::make_unique<CSRMatrix>(std::move(rows), label_encoder.size()),
            std::move(label_encoder)
        );
    }

    using RowSet = tsl::ordered_set<BinaryMatrix::SetBitPositions,
                                    utils::VectorHash,
                                    std::equal_to<BinaryMatrix::SetBitPositions>,
                                    std::allocator<BinaryMatrix::SetBitPositions>,
                                    std::vector<BinaryMatrix::SetBitPositions>,
                                    uint32_t>;
    RowSet unique_rows { BinaryMatrix::SetBitPositions() };
    std::vector<uint32_t> row_rank(num_rows, 0);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint64_t batch_begin = 0;
                        batch_begin < full_to_small.size();
                                        batch_begin += kRowBatchSize) {
        const uint64_t batch_end
            = std::min(batch_begin + kRowBatchSize,
                       static_cast<uint64_t>(full_to_small.size()));

        std::vector<uint64_t> row_indexes;
        row_indexes.reserve(batch_end - batch_begin);
        for (uint64_t i = batch_begin; i < batch_end; ++i) {
            assert(full_to_small[i].first < full_annotation.num_objects());

            row_indexes.push_back(full_to_small[i].first);
        }

        auto rows = full_annotation.get_matrix().get_rows(row_indexes);

        assert(rows.size() == batch_end - batch_begin);

        #pragma omp critical
        {
            for (uint64_t i = batch_begin; i < batch_end; ++i) {
                const auto &row = rows[i - batch_begin];
                auto it = unique_rows.emplace(row).first;
                row_rank[full_to_small[i].second] = it - unique_rows.begin();
                if (unique_rows.size() == std::numeric_limits<uint32_t>::max())
                    throw std::runtime_error("There must be less than 2^32 unique rows."
                                             " Reduce the query batch size.");
            }
        }
    }

    auto &annotation_rows = const_cast<std::vector<BinaryMatrix::SetBitPositions>&>(
        unique_rows.values_container()
    );

    auto label_encoder = reencode_labels(full_annotation.get_label_encoder(), &annotation_rows);

    // copy annotations from the full graph to the query graph
    return std::make_unique<annot::UniqueRowAnnotator>(
        std::make_unique<UniqueRowBinmat>(std::move(annotation_rows),
                                          std::move(row_rank),
                                          label_encoder.size()),
        std::move(label_encoder)
    );
}

void add_nodes_with_suffix_matches(const DBGSuccinct &full_dbg,
                                   size_t sub_k,
                                   size_t max_num_nodes_per_suffix,
                                   std::vector<std::pair<std::string, std::vector<node_index>>> *contigs,
                                   bool check_reverse_complement) {
    std::vector<std::pair<std::string, std::vector<node_index>>> contig_buffer;

    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t i = 0; i < contigs->size(); ++i) {
        const auto &[contig, path] = (*contigs)[i];
        std::vector<std::pair<std::string, node_index>> added_nodes;
        auto callback = [&](std::string&& kmer, node_index node) {
            assert(node == full_dbg.kmer_to_node(kmer));

            if (full_dbg.get_mode() == DeBruijnGraph::CANONICAL)
                full_dbg.map_to_nodes(kmer, [&](node_index cn) { node = cn; });

            added_nodes.emplace_back(std::move(kmer), node);
        };
        call_suffix_match_sequences(full_dbg, contig, path,
                                    callback, sub_k, max_num_nodes_per_suffix);

        if (check_reverse_complement) {
            std::string rev_contig = contig;
            reverse_complement(rev_contig);
            call_suffix_match_sequences(full_dbg, rev_contig,
                                        map_sequence_to_nodes(full_dbg, rev_contig),
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
        }
    }

    for (auto&& pair : contig_buffer) {
        contigs->emplace_back(std::move(pair));
    }
}

void add_hull_contigs(const DeBruijnGraph &full_dbg,
                      const DeBruijnGraph &batch_graph,
                      size_t max_hull_forks,
                      size_t max_hull_depth,
                      std::vector<std::pair<std::string, std::vector<node_index>>> *contigs) {
    tsl::hopscotch_map<node_index, uint32_t> distance_traversed_until_node;

    std::mutex mu;
    std::vector<std::pair<std::string, std::vector<node_index>>> contig_buffer;

    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t i = 0; i < contigs->size(); ++i) {
        const auto &[contig, path] = (*contigs)[i];
        std::vector<std::pair<std::string, std::vector<node_index>>> added_paths;
        // TODO: combine these two callbacks into one?
        auto callback = [&](const std::string &sequence,
                            const std::vector<node_index> &path) {
            added_paths.emplace_back(sequence, path);
            if (full_dbg.get_mode() == DeBruijnGraph::CANONICAL) {
                added_paths.back().second.resize(0);
                full_dbg.map_to_nodes(sequence,
                    [&](node_index cn) { added_paths.back().second.push_back(cn); }
                );
            }
        };
        auto continue_traversal = [&](std::string_view seq,
                                      node_index last_node,
                                      size_t depth,
                                      size_t fork_count) {
            if (fork_count > max_hull_forks || depth >= max_hull_depth)
                return false;

            // if the last node is already in the graph, cut off traversal
            // since this node will be covered in another traversal
            if (batch_graph.find(seq.substr(seq.length() - batch_graph.get_k()))) {
                return false;
            }

            // when a node which has already been accessed is visited,
            // only continue traversing if the previous access was in a
            // longer path (i.e., it cut off earlier)

            // TODO: check the number of forks too (shorter paths may have more forks)
            std::lock_guard<std::mutex> lock(mu);

            auto [it, inserted]
                = distance_traversed_until_node.emplace(last_node, depth);

            bool extend = inserted || depth < it->second;

            if (!inserted && depth < it->second)
                it.value() = depth;

            return extend;
        };

        for (size_t j = 0; j < path.size(); ++j) {
            if (!path[j]) {
                // If a node is unmatched, start expansion from all the nodes
                // adjacent to it.
                std::string kmer = contig.substr(j, full_dbg.get_k());
                // forward expansion from the incoming nodes
                batch_graph.adjacent_incoming_nodes(batch_graph.kmer_to_node(kmer),
                    [&](node_index next) {
                        call_hull_sequences(full_dbg, batch_graph.get_node_sequence(next),
                                            callback, continue_traversal);
                    }
                );
                if (batch_graph.get_mode() == DeBruijnGraph::CANONICAL) {
                    // In canonical graphs, (incoming to forward) is equivalent
                    // to (outgoing from rev-comp), so the following code invokes
                    // backward expansion.
                    reverse_complement(kmer);
                    batch_graph.adjacent_incoming_nodes(batch_graph.kmer_to_node(kmer),
                        [&](node_index next) {
                            call_hull_sequences(full_dbg, batch_graph.get_node_sequence(next),
                                                callback, continue_traversal);
                        }
                    );
                }
                // TODO: For non-canonical graphs, also expand backwards from
                //       the outgoing nodes?
            }
        }

        std::string last_kmer = contig.substr(contig.length() - full_dbg.get_k(), full_dbg.get_k());
        if (!batch_graph.outdegree(batch_graph.kmer_to_node(last_kmer)))
            call_hull_sequences(full_dbg, last_kmer, callback, continue_traversal);

        if (batch_graph.get_mode() == DeBruijnGraph::CANONICAL) {
            last_kmer = contig.substr(0, full_dbg.get_k());
            reverse_complement(last_kmer);
            if (!batch_graph.outdegree(batch_graph.kmer_to_node(last_kmer)))
                call_hull_sequences(full_dbg, last_kmer, callback, continue_traversal);
        }

        #pragma omp critical
        {
            for (auto&& pair : added_paths) {
                assert(pair.second.size() == pair.first.length() - full_dbg.get_k() + 1);
                contig_buffer.emplace_back(std::move(pair));
            }
        }
    }

    for (auto&& pair : contig_buffer) {
        contigs->emplace_back(std::move(pair));
    }
}

template <class Graph, class Contigs>
void add_to_graph(Graph &graph, const Contigs &contigs, size_t k) {
    for (const auto &[contig, nodes_in_full] : contigs) {
        size_t begin = 0;
        size_t end;
        do {
            end = std::find(nodes_in_full.begin() + begin,
                            nodes_in_full.end(),
                            DeBruijnGraph::npos)
                    - nodes_in_full.begin();

            if (begin != end) {
                graph.add_sequence(std::string_view(
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
                      const Config *config) {
    const auto &full_dbg = anno_graph.get_graph();
    const auto &full_annotation = anno_graph.get_annotation();
    const auto *dbg_succ = dynamic_cast<const DBGSuccinct *>(&full_dbg);

    assert(full_dbg.get_mode() != DeBruijnGraph::PRIMARY
            && "primary graphs must be wrapped into canonical");

    size_t sub_k = full_dbg.get_k();
    size_t max_hull_forks = 0;
    size_t max_hull_depth = 0;
    size_t max_num_nodes_per_suffix = 1;
    double max_hull_depth_per_seq_char = 0.0;
    if (config) {
        if (config->alignment_min_seed_length > full_dbg.get_k()) {
            logger->warn("Can't match suffixes longer than k={}."
                         " The value of k={} will be used.",
                         full_dbg.get_k(), full_dbg.get_k());
        }
        if (config->alignment_min_seed_length
                && config->alignment_min_seed_length < full_dbg.get_k()) {
            if (!dbg_succ) {
                logger->error("Matching suffixes of k-mers only supported for DBGSuccinct");
                exit(1);
            }
            sub_k = config->alignment_min_seed_length;
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

    logger->trace("[Query graph construction] Building the batch graph...");

    if (kPrefilterWithBloom && dbg_succ && sub_k == full_dbg.get_k()) {
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
            if (max_input_sequence_length < sequence.length())
                max_input_sequence_length = sequence.length();
        });
    } else {
        call_sequences([&](const std::string &sequence) {
            graph_init->add_sequence(sequence);
            if (max_input_sequence_length < sequence.length())
                max_input_sequence_length = sequence.length();
        });
    }

    max_hull_depth = std::min(
        max_hull_depth,
        static_cast<size_t>(max_hull_depth_per_seq_char * max_input_sequence_length)
    );

    logger->trace("[Query graph construction] Batch graph contains {} k-mers"
                  " and took {} sec to construct",
                  graph_init->num_nodes(), timer.elapsed());
    timer.reset();

    // pull contigs from query graph
    std::vector<std::pair<std::string, std::vector<node_index>>> contigs;
    std::mutex seq_mutex;
    graph_init->call_sequences([&](const std::string &contig, const auto &) {
                                   std::lock_guard<std::mutex> lock(seq_mutex);
                                   contigs.emplace_back(contig, std::vector<node_index>{});
                               },
                               get_num_threads(),
                               // pull only primary contigs when building canonical query graph
                               full_dbg.get_mode() == DeBruijnGraph::CANONICAL);

    logger->trace("[Query graph construction] Contig extraction took {} sec", timer.elapsed());
    timer.reset();

    logger->trace("[Query graph construction] Mapping k-mers back to full graph...");
    // map from nodes in query graph to full graph
    #pragma omp parallel for num_threads(get_num_threads())
    for (size_t i = 0; i < contigs.size(); ++i) {
        contigs[i].second.reserve(contigs[i].first.length() - graph_init->get_k() + 1);
        full_dbg.map_to_nodes(contigs[i].first,
                              [&](node_index node) { contigs[i].second.push_back(node); });
    }
    logger->trace("[Query graph construction] Contigs mapped to the full graph in {} sec",
                  timer.elapsed());
    timer.reset();

    size_t original_size = contigs.size();

    // add nodes with suffix matches to the query
    if (sub_k < full_dbg.get_k()) {
        assert(dbg_succ);
        logger->trace("[Query graph construction] Adding k-mers with matching "
                      "suffixes of length {}...", sub_k);
        timer.reset();

        add_nodes_with_suffix_matches(*dbg_succ, sub_k, max_num_nodes_per_suffix,
                                      &contigs, full_dbg.get_mode() == DeBruijnGraph::CANONICAL);

        logger->trace("[Query graph construction] Found {} suffix-matching k-mers, took {} sec",
                      contigs.size() - original_size, timer.elapsed());
    }

    if (max_hull_forks) {
        logger->trace("[Query graph augmentation] Computing query graph hull...");
        logger->trace("[Query graph augmentation] max traversal distance: {}, max fork count: {}",
                      max_hull_depth, max_hull_forks);
        timer.reset();

        // add k-mers with sub_k-suffix matches
        for (size_t i = original_size; i < contigs.size(); ++i) {
            graph_init->add_sequence(contigs[i].first);
        }

        size_t hull_contigs_begin = contigs.size();

        add_hull_contigs(full_dbg, *graph_init, max_hull_forks, max_hull_depth, &contigs);

        logger->trace("[Query graph augmentation] Augmented the batch graph with {} contigs in {} sec",
                      contigs.size() - hull_contigs_begin, timer.elapsed());
    }

    graph_init.reset();

    logger->trace("[Query graph construction] Building the query graph...");
    timer.reset();
    std::shared_ptr<DeBruijnGraph> graph;

    // restrict nodes to those in the full graph
    if (sub_k < full_dbg.get_k()) {
        BOSSConstructor constructor(full_dbg.get_k() - 1,
                                    full_dbg.get_mode() == DeBruijnGraph::CANONICAL,
                                    0, "", num_threads);
        add_to_graph(constructor, contigs, full_dbg.get_k());

        graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor), full_dbg.get_mode());

    } else {
        graph = std::make_shared<DBGHashOrdered>(full_dbg.get_k(), full_dbg.get_mode());
        add_to_graph(*graph, contigs, full_dbg.get_k());
    }

    logger->trace("[Query graph construction] Query graph contains {} k-mers"
                  " and took {} sec to construct",
                  graph->num_nodes(), timer.elapsed());
    timer.reset();

    logger->trace("[Query graph construction] Mapping the contigs back to the query graph...");

    std::vector<std::pair<uint64_t, uint64_t>> from_full_to_small;

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < contigs.size(); ++i) {
        const std::string &contig = contigs[i].first;
        const std::vector<node_index> &nodes_in_full = contigs[i].second;

        std::vector<uint64_t> path(nodes_in_full.size());
        size_t j = 0;
        // nodes in the query graph hull may overlap
        graph->map_to_nodes(contig, [&](node_index node) {
            path[j++] = node;
        });
        assert(j == nodes_in_full.size());

        #pragma omp critical
        {
            for (size_t j = 0; j < path.size(); ++j) {
                if (nodes_in_full[j]) {
                    assert(path[j]);
                    from_full_to_small.emplace_back(nodes_in_full[j], path[j]);
                }
            }
        }
    }

    logger->trace("[Query graph construction] Mapping between graphs constructed in {} sec",
                  timer.elapsed());

    contigs = decltype(contigs)();

    for (auto &[first, second] : from_full_to_small) {
        assert(first && second);
        first = AnnotatedDBG::graph_to_anno_index(first);
        second = AnnotatedDBG::graph_to_anno_index(second);
    }

    logger->trace("[Query graph construction] Slicing {} rows out of full annotation...",
                  from_full_to_small.size());

    // initialize fast query annotation
    // copy annotations from the full graph to the query graph
    auto annotation = slice_annotation(full_annotation,
                                       graph->max_index(),
                                       std::move(from_full_to_small),
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

    std::shared_ptr<DeBruijnGraph> graph = load_critical_dbg(config->infbase);

    std::unique_ptr<AnnotatedDBG> anno_graph = initialize_annotated_dbg(graph, *config);

    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);

    Timer timer;

    std::unique_ptr<align::DBGAlignerConfig> aligner_config;
    if (config->align_sequences) {
        assert(config->alignment_num_alternative_paths == 1u
                && "only the best alignment is used in query");

        aligner_config.reset(new align::DBGAlignerConfig(
            initialize_aligner_config(*config)
        ));
    }

    QueryExecutor executor(*config, *anno_graph, std::move(aligner_config), thread_pool);

    // iterate over input files
    for (const auto &file : files) {
        Timer curr_timer;

        // Callback, which captures the config pointer and a const reference to the anno_graph
        // instance pointed to by our unique_ptr...
        executor.query_fasta(file,
                [config, &anno_graph
                        = std::as_const(*anno_graph)](const SeqSearchResult &result) {
            std::cout << result.to_string(config->anno_labels_delimiter,
                                          config->suppress_unlabeled,
                                          config->expand_coords,
                                          anno_graph) << "\n";
        });
        logger->trace("File '{}' was processed in {} sec, total time: {}", file,
                      curr_timer.elapsed(), timer.elapsed());
    }

    return 0;
}


/**
 * Align a sequence to the annotated DBG before querying. Modify the sequence in place and return
 * alignment information (score and cigar string).
 *
 * @param seq               reference to sequence string (which will be modified if alignment found)
 * @param anno_graph        reference to annotated DBG we align to
 * @param aligner_config    alignment config
 *
 * @return Alignment struct instance with score and cigar string
 */
QuerySequence::Alignment align_sequence(std::string &seq,
                                        const AnnotatedDBG &anno_graph,
                                        const align::DBGAlignerConfig &aligner_config) {
    const DeBruijnGraph &graph = anno_graph.get_graph();
    auto alignments = build_aligner(graph, aligner_config)->align(seq);

    assert(alignments.size() <= 1 && "Only the best alignment is needed");

    if (alignments.size()) {
        auto &match = alignments[0];
        // modify sequence for querying with the best alignment
        if (match.get_offset()) {
            seq = graph.get_node_sequence(match.get_nodes()[0]).substr(0, match.get_offset())
                  + match.get_sequence();
        } else {
            seq = const_cast<std::string&&>(match.get_sequence());
        }

        return QuerySequence::Alignment(match.get_score(), match.get_cigar().to_string());
    } else {
        return QuerySequence::Alignment(0, fmt::format("{}S", seq.length()));
    }
}


SeqSearchResult query_sequence(QuerySequence&& sequence,
                               const AnnotatedDBG &anno_graph,
                               const Config &config,
                               const align::DBGAlignerConfig *aligner_config) {
    // Align sequence and modify sequence struct to store alignment results
    if (aligner_config) {
        sequence.alignment = align_sequence(sequence.sequence, anno_graph, *aligner_config);
    };

    // Get sequence search result by executing query
    SeqSearchResult result = QueryExecutor::execute_query(std::move(sequence),
            config.count_labels, config.print_signature,
            config.num_top_labels, config.discovery_fraction,
            config.presence_fraction, anno_graph,
            config.count_kmers, config.count_quantiles,
            config.query_coords);

    return result;
}


void QueryExecutor::query_fasta(const string &file,
                                const std::function<void(const SeqSearchResult &)> &callback) {
    logger->trace("Parsing sequences from file '{}'", file);

    seq_io::FastaParser fasta_parser(file, config_.forward_and_reverse);

    if (config_.fast) {
        // TODO: Implement batch mode for query_coords queries
        if (config_.query_coords) {
            logger->error("Querying coordinates in batch mode is not supported");
            exit(1);
        }
        // Construct a query graph and query against it
        batched_query_fasta(fasta_parser, callback);
        return;
    }

    // Query sequences independently
    size_t seq_count = 0;

    for (const seq_io::kseq_t &kseq : fasta_parser) {
        thread_pool_.enqueue([&](const auto&... args) {
            // Create a QuerySequence instance and pass to query_sequence
            QuerySequence sequence(args...);

            // Callback with the SeqSearchResult
            callback(query_sequence(std::move(sequence),
                                    anno_graph_,
                                    config_,
                                    aligner_config_.get()));
        }, seq_count++, std::string(kseq.name.s), std::string(kseq.seq.s));
    }

    // wait while all threads finish processing the current file
    thread_pool_.join();
}

void QueryExecutor
::batched_query_fasta(seq_io::FastaParser &fasta_parser,
                      const std::function<void(const SeqSearchResult &)> &callback) {
    auto it = fasta_parser.begin();
    auto end = fasta_parser.end();

    const uint64_t batch_size = config_.query_batch_size_in_bytes;

    size_t seq_count = 0;

    while (it != end) {
        Timer batch_timer;

        uint64_t num_bytes_read = 0;

        // A generator that can be called multiple times until all sequences
        // are called
        std::vector<QuerySequence> seq_batch;
        num_bytes_read = 0;

        for ( ; it != end && num_bytes_read <= batch_size; ++it) {
            seq_batch.emplace_back(seq_count++, it->name.s, it->seq.s);
            num_bytes_read += it->seq.l;
        }

        // Align sequences ahead of time on full graph if we don't have batch_align
        if (aligner_config_ && !config_.batch_align) {
            logger->trace("Aligning sequences from batch against the full graph...");
            batch_timer.reset();

            #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
            for (size_t i = 0; i < seq_batch.size(); ++i) {
                // Set alignment for this seq_batch
                seq_batch[i].alignment = align_sequence(seq_batch[i].sequence,
                                                        anno_graph_,
                                                        *aligner_config_);
            }
            logger->trace("Sequences alignment took {} sec", batch_timer.elapsed());
            batch_timer.reset();
        }

        // Construct the query graph for this batch
        auto query_graph = construct_query_graph(
            anno_graph_,
            [&](auto callback) {
                for (const auto &seq : seq_batch) {
                    callback(seq.sequence);
                }
            },
            get_num_threads(),
            aligner_config_ && config_.batch_align ? &config_ : NULL
        );

        logger->trace("Query graph constructed for batch of sequences"
                      " with {} bases from '{}' in {} sec",
                      num_bytes_read, fasta_parser.get_filename(), batch_timer.elapsed());
        batch_timer.reset();

        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (size_t i = 0; i < seq_batch.size(); ++i) {
            callback(query_sequence(std::move(seq_batch[i]),
                                    *query_graph,
                                    config_,
                                    config_.batch_align ? aligner_config_.get() : NULL));
        }

        logger->trace("Batch of {} bytes from '{}' queried in {} sec", num_bytes_read,
                      fasta_parser.get_filename(), batch_timer.elapsed());
    }
}

} // namespace cli
} // namespace mtg
