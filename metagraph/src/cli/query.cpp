#include "query.hpp"

#include <mutex>
#include <sstream>

#include <ips4o.hpp>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/hashers/hash.hpp"
#include "common/utils/template_utils.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/representation/canonical_dbg.hpp"
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

const char ALIGNED_SEQ_HEADER_FORMAT[] = "{}:{}:{}:{}";

using namespace mtg::graph;
using namespace mtg::annot::matrix;

using mtg::common::logger;
using mtg::graph::boss::BOSS;
using mtg::graph::boss::BOSSConstructor;

typedef typename mtg::graph::DeBruijnGraph::node_index node_index;


// Format a range of coordinates (start, first_coord, last_coord) to a string representation
std::string get_range(const std::tuple<uint64_t, uint64_t, uint64_t> &range) {
    const auto &[pos, first, last] = range;
    if (first == last) {
        return fmt::format("{}-{}", pos, first);
    } else {
        return fmt::format("{}-{}-{}", pos, first, last);
    }
}

/**
 * Given a vector of kmer label matched coordinates, collapse continuous ranges of coordinates
 * to start-end tuples.
 *
 * @param coords    the vector of tuples as stored in a SeqSearchResult
 *                  originally returned from AnnotatedDBG::get_kmer_coordinates
 * @return vector of 'begin-end' range string representations
 */
std::vector<std::string>
collapse_coord_ranges(const std::vector<SmallVector<uint64_t>> &tuples) {
    // Build output
    std::vector<std::string> range_strings;

    // Keep track of the ranges: start position, first coord, last coord
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> ranges;
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> next_ranges;

    for (size_t i = 0; i < tuples.size(); ++i) {
        const auto &coords = tuples[i];
        assert(std::is_sorted(coords.begin(), coords.end()));

        size_t j = 0;
        next_ranges.resize(0);

        for (uint64_t coord : coords) {
            while (j < ranges.size() && std::get<2>(ranges[j]) + 1 < coord) {
                // end of range
                range_strings.push_back(get_range(ranges[j++]));
            }

            if (j < ranges.size() && std::get<2>(ranges[j]) + 1 == coord) {
                // extend range
                next_ranges.push_back(ranges[j++]);
                std::get<2>(next_ranges.back())++;
            } else {
                // start new range
                next_ranges.emplace_back(i, coord, coord);
            }
        }
        while (j < ranges.size()) {
            range_strings.push_back(get_range(ranges[j++]));
        }

        ranges.swap(next_ranges);
    }

    for (const auto &range : ranges) {
        range_strings.push_back(get_range(range));
    }

    return range_strings;
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
    for (size_t i = 1; i < label_parts.size(); ++i) {
        std::vector<std::string> key_value = utils::split_string(label_parts[i], "=");
        if (key_value.size() != 2) {
            logger->error("Can't read key-value pair in part {} of label {}", label_parts[i], label);
            throw std::runtime_error("Label formatting error");
        }
        properties[key_value[0]] = adjust_for_types(key_value[1]);
    }

    if (properties.size() > 0) {
        label_root["properties"] = properties;
    }

    return label_root;
}



/**
 * Helper functions for score_kmer_presence_mask
 */

std::array<AlignedVector<size_t>, 2>
tabulate_score(const sdsl::bit_vector &presence, size_t correction = 0) {
    std::array<AlignedVector<size_t>, 2> table;
    table[0].reserve(presence.size());
    table[1].reserve(presence.size());

    if (!presence.size())
        return table;

    bool last_block = presence[0];
    size_t last_size = 1;
    size_t i = 1;

    if (presence.size() >= 64) {
        const auto word = *presence.data();
        if (!word || word == 0xFFFFFFFFFFFFFFFF) {
            last_size = 64;
            i = 64;
        }
    }

    for ( ; i < presence.size(); ++i) {
        if (!(i & 0x3F) && i + 64 <= presence.size()) {
            // if at a word boundary and the next word is either all zeros or
            // all ones
            const uint64_t word = presence.get_int(i);
            if (!word || word == 0xFFFFFFFFFFFFFFFF) {
                if ((!word && last_block) || (word == 0xFFFFFFFFFFFFFFFF && !last_block)) {
                    table[last_block].push_back(last_size + correction);
                    last_block = !last_block;
                    last_size = 0;
                }

                last_size += 64;
                i += 63;
                continue;
            }
        }

        if (last_block == presence[i]) {
            ++last_size;
        } else {
            table[last_block].push_back(last_size + correction);
            last_block = !last_block;
            last_size = 1;
        }
    }

    table[last_block].push_back(last_size);

    assert(std::accumulate(table[0].begin(), table[0].end(), size_t(0))
         + std::accumulate(table[1].begin(), table[1].end(), size_t(0))
         - correction * (table[0].size() + table[1].size() - 1)
        == presence.size());

    assert(std::accumulate(table[1].begin(), table[1].end(), size_t(0))
         - correction * (table[1].size() - last_block)
        == sdsl::util::cnt_one_bits(presence));

#ifndef NDEBUG
    if (correction == 1) {
        std::array<AlignedVector<size_t>, 2> check;

        size_t cnt = 1;
        for (size_t i = 0; i < presence.size(); ++i) {
            if (i < presence.size() - 1) {
                ++cnt;
                if (presence[i] != presence[i + 1]) {
                    check[presence[i]].push_back(cnt);
                    cnt = 1;
                }
            } else {
                check[presence[i]].push_back(cnt);
            }
        }

        assert(check[0] == table[0]);
        assert(check[1] == table[1]);
    }
#endif // NDEBUG

    return table;
}

int32_t score_kmer_presence_mask(size_t k,
                                 const sdsl::bit_vector &kmer_presence_mask,
                                 int32_t match_score = 1,
                                 int32_t mismatch_score = 2) {
    if (!kmer_presence_mask.size())
        return 0;

    const int32_t kmer_adjust = 3;

    const size_t sequence_length = kmer_presence_mask.size() + k - 1;
    const int32_t SNP_t = k + kmer_adjust;

    auto score_counter = tabulate_score(autocorrelate(kmer_presence_mask, kmer_adjust), 1);

    double score = std::accumulate(score_counter[1].begin(),
                                   score_counter[1].end(), 0) * match_score;
    if (score == 0)
        return 0;

    if (score_counter[0].empty())
        return score * sequence_length / kmer_presence_mask.size();

    for (double count : score_counter[0]) {
        // A penalty function used in BIGSI
        double min_N_snps = count / SNP_t;
        double max_N_snps = std::max(count - SNP_t + 1, min_N_snps);
        double mean_N_snps = max_N_snps * 0.05 + min_N_snps;

        assert(count >= mean_N_snps);

        double mean_penalty = mean_N_snps * mismatch_score;
        score += (count - mean_penalty) * match_score - mean_penalty;
    }

    return std::max(score * sequence_length / kmer_presence_mask.size(), 0.);
}


Json::Value SeqSearchResult::to_json(bool verbose_output, size_t k) const {
    Json::Value root;

    // Add seq information
    root[SEQ_DESCRIPTION_JSON_FIELD] = sequence_.name;

    // Add alignment information if there is any
    if (alignment_) {
        // Recover alignment sequence sequence string
        root[SEQUENCE_JSON_FIELD] = Json::Value(sequence_.sequence);

        // Alignment metrics
        root[SCORE_JSON_FIELD] = Json::Value(alignment_->score);
        root[MAX_SCORE_JSON_FIELD] = Json::Value(alignment_->max_score);
        root[CIGAR_JSON_FIELD] = Json::Value(alignment_->cigar);
        root[ORIENTATION_JSON_FIELD] = Json::Value(alignment_->orientation);
    }

    // Add discovered labels and extra results
    root["results"] = Json::Value(Json::arrayValue);

    // Different action depending on the result type
    if (const auto *v = std::get_if<LabelVec>(&result_)) {
        // Standard labels only
        for (const auto &label : *v) {
            root["results"].append(get_label_as_json(label));
        }
    } else if (const auto *v = std::get_if<LabelCountVec>(&result_)) {
        // Labels with count data
        for (const auto &[label, count] : *v) {
            Json::Value &label_obj = root["results"].append(get_label_as_json(label));
            label_obj[KMER_COUNT_FIELD] = static_cast<Json::Int64>(count);
        }
    } else if (const auto *v = std::get_if<LabelSigVec>(&result_)) {
        // Count signatures
        for (const auto &[label, kmer_presence_mask] : *v) {
            Json::Value &label_obj = root["results"].append(get_label_as_json(label));
            // Store the presence mask and score in a separate object
            Json::Value &sig_obj = (label_obj[SIGNATURE_FIELD] = Json::objectValue);
            sig_obj["presence_mask"] = Json::Value(sdsl::util::to_string(kmer_presence_mask));
            sig_obj["score"] = Json::Value(score_kmer_presence_mask(k, kmer_presence_mask));
            // Add kmer_counts calculated using bitmask
            label_obj[KMER_COUNT_FIELD] = Json::Value(sdsl::util::cnt_one_bits(kmer_presence_mask));
        }
    } else if (const auto *v = std::get_if<LabelCountAbundancesVec>(&result_)) {
        // k-mer counts (or quantiles)
        for (const auto &[label, count, abundances] : *v) {
            assert(abundances.size());
            Json::Value &label_obj = root["results"].append(get_label_as_json(label));
            label_obj[KMER_COUNT_FIELD] = static_cast<Json::Int64>(count);
            Json::Value &counts_array = (label_obj[KMER_ABUNDANCE_FIELD] = Json::arrayValue);
            if (verbose_output) {
                for (size_t c : abundances) {
                    counts_array.append(static_cast<Json::Int64>(c));
                }
            } else {
                std::pair<size_t, size_t> last(0, abundances.at(0));
                for (size_t i = 1; i <= abundances.size(); ++i) {
                    if (i < abundances.size() && abundances[i] == last.second)
                        continue; // extend run

                    // end of run
                    if (last.second) {
                        if (i == last.first + 1) {
                            counts_array.append(fmt::format("{}={}", last.first, last.second));
                        } else {
                            counts_array.append(fmt::format("{}-{}={}", last.first, i - 1, last.second));
                        }
                    }

                    // start new run
                    if (i < abundances.size()) {
                        last.first = i;
                        last.second = abundances[i];
                    }
                }
            }
        }
    } else {
        // Kmer coordinates
        for (const auto &[label, count, tuples] : std::get<LabelCountCoordsVec>(result_)) {
            Json::Value &label_obj = root["results"].append(get_label_as_json(label));
            label_obj[KMER_COUNT_FIELD] = static_cast<Json::Int64>(count);
            Json::Value &coord_array = (label_obj[KMER_COORDINATE_FIELD] = Json::arrayValue);

            // Each tuple is represented as a folly::SmallVector<uint64_t> instance
            if (verbose_output) {
                for (const auto &coords : tuples) {
                    coord_array.append(fmt::format("{}", fmt::join(coords, ",")));
                }
            } else {
                for (const auto &range : collapse_coord_ranges(tuples)) {
                    coord_array.append(range);
                }
            }
        }
    }

    return root;
}


std::string SeqSearchResult::to_string(const std::string delimiter,
                                       bool suppress_unlabeled,
                                       bool verbose_output,
                                       size_t k) const {
    // Return the size of the result type vector using std::visit (ducktyping variant)
    // If empty and unlabeled sequences are supressed then return empty string
    if (suppress_unlabeled && !std::visit([] (auto&& vec) { return vec.size(); }, result_))
        return "";

    std::string output;

    // Print sequence ID and name
    if (alignment_) {
        // Get modified name if sequence has an alignment result
        const std::string &mod_seq_name
            = fmt::format(ALIGNED_SEQ_HEADER_FORMAT, sequence_.name, sequence_.sequence,
                                   alignment_->score, alignment_->cigar);
        output = fmt::format("{}\t{}", sequence_.id, mod_seq_name);
    } else {
        output = fmt::format("{}\t{}", sequence_.id, sequence_.name);
    }

    // ... with diff result output depending on the type of result being stored.
    if (const auto *v = std::get_if<LabelVec>(&result_)) {
        // Standard labels only
        output += fmt::format("\t{}", utils::join_strings(*v, delimiter));

    } else if (const auto *v = std::get_if<LabelCountVec>(&result_)) {
        // Labels with count data
        for (const auto &[label, count] : *v) {
            output += fmt::format("\t<{}>:{}", label, count);
        }
    } else if (const auto *v = std::get_if<LabelSigVec>(&result_)) {
        // Count signatures
        for (const auto &[label, kmer_presence_mask] : *v) {
            output += fmt::format("\t<{}>:{}:{}:{}", label,
                                  sdsl::util::cnt_one_bits(kmer_presence_mask),
                                  sdsl::util::to_string(kmer_presence_mask),
                                  score_kmer_presence_mask(k, kmer_presence_mask));
        }
    } else if (const auto *v = std::get_if<LabelCountAbundancesVec>(&result_)) {
        // k-mer counts (or quantiles)
        for (const auto &[label, count, abundances] : *v) {
            output += "\t<" + label + ">";
            if (verbose_output) {
                output += fmt::format(":{}", fmt::join(abundances, ":"));
            } else {
                std::pair<size_t, size_t> last(0, abundances.at(0));
                for (size_t i = 1; i <= abundances.size(); ++i) {
                    if (i < abundances.size() && abundances[i] == last.second)
                        continue; // extend run

                    // end of run
                    if (last.second) {
                        if (i == last.first + 1) {
                            output += fmt::format(":{}={}", last.first, last.second);
                        } else {
                            output += fmt::format(":{}-{}={}", last.first, i - 1, last.second);
                        }
                    }

                    // start new run
                    if (i < abundances.size()) {
                        last.first = i;
                        last.second = abundances[i];
                    }
                }
            }
        }
    } else {
        // Kmer coordinates
        for (const auto &[label, count, tuples] : std::get<LabelCountCoordsVec>(result_)) {
            output += "\t<" + label + ">";

            if (verbose_output) {
                for (const auto &coords : tuples) {
                    output += fmt::format(":{}", fmt::join(coords, ","));
                }
            } else {
                output += fmt::format(":{}", fmt::join(collapse_coord_ranges(tuples), ":"));
            }
        }
    }

    return output;
}


SeqSearchResult execute_query(QuerySequence&& sequence,
                              QueryMode query_mode,
                              size_t num_top_labels,
                              double discovery_fraction,
                              double presence_fraction,
                              const graph::AnnotatedDBG &anno_graph) {
    // Perform a different action depending on the type (specified by config flags)
    SeqSearchResult::result_type result;

    switch (query_mode) {
        case SIGNATURE: {
            // Get labels with presence/absence signatures
            result = anno_graph.get_top_label_signatures(sequence.sequence,
                                                         num_top_labels,
                                                         discovery_fraction,
                                                         presence_fraction);
            break;
        }
        case COORDS: {
            // Get labels with k-mer coordinates
            result = anno_graph.get_kmer_coordinates(sequence.sequence,
                                                     num_top_labels,
                                                     discovery_fraction,
                                                     presence_fraction);
            break;
        }
        case COUNTS: {
            // Get labels with k-mer counts
            result = anno_graph.get_kmer_counts(sequence.sequence,
                                                num_top_labels,
                                                discovery_fraction,
                                                presence_fraction);
            break;
        }
        case MATCHES: {
            // Get labels with label counts
            result = anno_graph.get_top_labels(sequence.sequence,
                                               num_top_labels,
                                               discovery_fraction,
                                               presence_fraction,
                                               false);
            break;
        }
        case COUNTS_SUM: {
            // Get labels with label counts / kmer counts
            result = anno_graph.get_top_labels(sequence.sequence,
                                               num_top_labels,
                                               discovery_fraction,
                                               presence_fraction,
                                               true);
            break;
        }
        case LABELS: {
            // Default, just get labels
            result = anno_graph.get_labels(sequence.sequence,
                                           discovery_fraction,
                                           presence_fraction);
            break;
        }
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
        assert(path == map_to_nodes_sequentially(full_dbg, seq));

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
    if (const auto *mat = dynamic_cast<const IntMatrix *>(&full_annotation.get_matrix())) {
        // don't break the topological order for row-diff annotation
        if (!dynamic_cast<const IRowDiff *>(&full_annotation.get_matrix())) {
            ips4o::parallel::sort(full_to_small.begin(), full_to_small.end(),
                                  utils::LessFirst(), num_threads);
        }

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

    // shortcut construction for Rainbow<> annotation
    std::vector<uint64_t> row_indexes(full_to_small.size());
    for (size_t i = 0; i < full_to_small.size(); ++i) {
        row_indexes[i] = full_to_small[i].first;
    }

    // get unique rows and set pointers to them in |row_indexes|
    auto unique_rows = full_annotation.get_matrix().get_rows_dict(&row_indexes, num_threads);

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

void add_nodes_with_suffix_matches(const DBGSuccinct &full_dbg,
                                   size_t sub_k,
                                   size_t max_num_nodes_per_suffix,
                                   std::vector<std::pair<std::string, std::vector<node_index>>> *contigs,
                                   bool check_reverse_complement) {
    std::vector<std::pair<std::string, std::vector<node_index>>> contig_buffer;

    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic, 10)
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
                                        map_to_nodes_sequentially(full_dbg, rev_contig),
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

    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic, 10)
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

            if (begin != end)
                graph.add_sequence(std::string_view(contig.data() + begin, end - begin + k - 1));

            begin = end + 1;
        } while (end < nodes_in_full.size());
    }
}

/**
 * Construct contigs from the query sequences.
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
 */
std::pair<std::vector<std::pair<std::string, std::vector<node_index>>>, size_t>
construct_contigs(const DeBruijnGraph &full_dbg,
                  const std::vector<QuerySequence> &seq_batch,
                  size_t num_threads,
                  const Config *config,
                  size_t *sub_k_ptr) {
    const auto *dbg_succ = dynamic_cast<const DBGSuccinct *>(&full_dbg);

    assert(full_dbg.get_mode() != DeBruijnGraph::PRIMARY
            && "primary graphs must be wrapped into canonical");

    assert(sub_k_ptr);
    auto &sub_k = *sub_k_ptr;
    sub_k = full_dbg.get_k();
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
    size_t k = full_dbg.get_k();
    auto graph_init = std::make_shared<DBGHashOrdered>(k, full_dbg.get_mode());
    size_t max_input_sequence_length = 0;

    logger->trace("[Query graph construction] Building the batch graph...");

    std::vector<std::pair<std::string, std::vector<node_index>>> contigs(seq_batch.size());

    // map from nodes in query graph to full graph
    std::atomic<uint64_t> num_found_kmers = 0;

    for (size_t i = 0; i < seq_batch.size(); ++i) {
        const auto &str = seq_batch[i].sequence;
        contigs[i].first = str;
        auto &path = contigs[i].second;
        path.assign(str.size() - k + 1, 0);
#ifndef NDEBUG
        uint64_t offset = graph_init->num_nodes();
#endif
        graph_init->add_sequence(str, [](){ return false; }, [&](size_t i, node_index node) {
            if (i < path.size()) {
                assert(node > offset);
                path[i] = node;
            }
        });
        if (max_input_sequence_length < str.length())
            max_input_sequence_length = str.length();
    }

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
    for (size_t i = 0; i < contigs.size(); ++i) {
        const auto &str = contigs[i].first;
        auto &path = contigs[i].second;
        if (static_cast<size_t>(std::count(path.begin(), path.end(), 0)) == path.size()) {
            path = std::vector<node_index>{};
            continue;  // no new k-mers
        }

        size_t begin = 0;
        size_t end;
        do {
            for (end = begin; end < path.size() && path[end]; ++end) {
            }

            if (begin != end) {
                std::string_view segment(str.data() + begin, end - begin + k - 1);
                // nodes in the query graph hull may overlap
                full_dbg.map_to_nodes(segment, [&](node_index node) {
                    path[begin++] = node;
                });
            }
            begin = end + 1;
        } while (end < path.size());

        num_found_kmers += (path.size() - std::count(path.begin(), path.end(), 0));
    }
    uint64_t num_kmers = graph_init->num_nodes() / (graph_init->get_mode() == DeBruijnGraph::CANONICAL ? 2 : 1);

    max_hull_depth = std::min(
        max_hull_depth,
        static_cast<size_t>(max_hull_depth_per_seq_char * max_input_sequence_length)
    );

    logger->trace("[Query graph construction] Batch graph contains {} k-mers (found in full {} / {} k-mers)"
                  " and took {} sec to construct",
                  graph_init->num_nodes(), num_found_kmers, num_kmers, timer.elapsed());
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

    return std::make_pair(std::move(contigs), (size_t)num_found_kmers);
}

auto construct_query_graph(std::vector<std::pair<std::string, std::vector<node_index>>>&& contigs,
                           size_t num_threads,
                           size_t k,
                           DeBruijnGraph::Mode mode,
                           size_t sub_k) {
    assert(mode != DeBruijnGraph::PRIMARY && "primary graphs must be wrapped into canonical");

    logger->trace("[Query graph construction] Building the query graph...");
    Timer timer;
    std::shared_ptr<DeBruijnGraph> graph;

    // restrict nodes to those in the full graph
    if (sub_k < k) {
        BOSSConstructor constructor(k - 1, mode == DeBruijnGraph::CANONICAL, 0, "", num_threads);
        add_to_graph(constructor, contigs, k);

        graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor), mode);

    } else {
        graph = std::make_shared<DBGHashOrdered>(k, mode);
        add_to_graph(*graph, contigs, k);
    }

    logger->trace("[Query graph construction] Query graph contains {} k-mers"
                  " and took {} sec to construct",
                  graph->num_nodes(), timer.elapsed());
    timer.reset();

    VectorMap<uint64_t, uint64_t> from_full_to_small_map;
    from_full_to_small_map.reserve(graph->num_nodes());

    for (size_t i = 0; i < contigs.size(); ++i) {
        const std::string &contig = contigs[i].first;
        const std::vector<node_index> &nodes_in_full = contigs[i].second;

        std::vector<uint64_t> path(nodes_in_full.size(), 0);
        size_t begin = 0;
        size_t end;
        do {
            end = std::find(nodes_in_full.begin() + begin,
                            nodes_in_full.end(),
                            DeBruijnGraph::npos)
                    - nodes_in_full.begin();

            if (begin != end) {
                std::string_view segment(contig.data() + begin, end - begin + k - 1);
                // nodes in the query graph hull may overlap
                graph->map_to_nodes(segment, [&](node_index node) {
                    path[begin++] = node;
                });
            }
            begin = end + 1;
        } while (end < nodes_in_full.size());

        for (size_t j = 0; j < path.size(); ++j) {
            if (nodes_in_full[j]) {
                assert(path[j]);
                from_full_to_small_map.try_emplace(nodes_in_full[j], path[j]);
            }
        }
    }

    auto &from_full_to_small
        = const_cast<std::vector<std::pair<uint64_t, uint64_t>>&>(from_full_to_small_map.values_container());

    logger->trace("[Query graph construction] Mapping between graphs constructed in {} sec",
                  timer.elapsed());

    contigs.clear();

    for (auto &[first, second] : from_full_to_small) {
        assert(first && second);
        first = AnnotatedDBG::graph_to_anno_index(first);
        second = AnnotatedDBG::graph_to_anno_index(second);
    }

    return std::make_pair(std::move(graph), std::move(from_full_to_small));
}

std::unique_ptr<AnnotatedDBG> construct_query_graph(const AnnotatedDBG::Annotator &full_annotation,
                                                    std::vector<std::pair<uint64_t, uint64_t>>&& from_full_to_small,
                                                    std::shared_ptr<DeBruijnGraph> graph) {
    Timer timer;
    // initialize fast query annotation
    // copy annotations from the full graph to the query graph
    auto annotation = slice_annotation(full_annotation,
                                       graph->max_index(),
                                       std::move(from_full_to_small),
                                       get_num_threads());

    logger->trace("[Query graph construction] Query annotation with {} rows, {} labels,"
                  " and {} set bits constructed in {} sec", annotation->num_objects(),
                  annotation->num_labels(), annotation->num_relations(), timer.elapsed());
    timer.reset();

    // build annotated graph from the query graph and copied annotations
    return std::make_unique<AnnotatedDBG>(graph, std::move(annotation));
}

std::unique_ptr<graph::AnnotatedDBG>
construct_query_graph(const graph::AnnotatedDBG &anno_graph,
                      StringGenerator call_sequences,
                      size_t num_threads,
                      const Config *config) {
    size_t sub_k;
    std::mutex mu;
    std::vector<QuerySequence> seq_batch;
    call_sequences([&](const auto &seq) { seq_batch.push_back(QuerySequence{ 0, "", seq }); });
    auto contigs = construct_contigs(anno_graph.get_graph(), seq_batch, num_threads, config, &sub_k).first;
    auto [query_graph, from_full_to_small] = construct_query_graph(
        std::move(contigs), num_threads,
        anno_graph.get_graph().get_k(), anno_graph.get_graph().get_mode(), sub_k
    );
    return construct_query_graph(anno_graph.get_annotator(), std::move(from_full_to_small), std::move(query_graph));
}


size_t batched_query_fasta(const std::string &file,
                           const std::function<void(const SeqSearchResult &)> &callback,
                           const Config &config,
                           const graph::align::DBGAlignerConfig *aligner_config,
                           const DeBruijnGraph &graph,
                           const std::function<const AnnotatedDBG::Annotator *()> &get_annotation);

void check_annotation(const Config &config, const annot::matrix::BinaryMatrix &anno_matrix);

int query_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase_annotators.size() == 1);

    std::shared_ptr<DeBruijnGraph> graph = load_critical_dbg(config->infbase);

    ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1, 1000);

    std::unique_ptr<align::DBGAlignerConfig> aligner_config;
    if (config->align_sequences) {
        assert(config->alignment_num_alternative_paths == 1u
                && "only the best alignment is used in query");

        aligner_config.reset(new align::DBGAlignerConfig(
            initialize_aligner_config(*config, *graph)
        ));
    }

    auto query_callback = [config, k=graph->get_k()](const SeqSearchResult &result) {
        if (config->output_json) {
            std::ostringstream ss;
            ss << result.to_json(config->verbose_output
                                     || !(config->query_mode == COUNTS || config->query_mode == COORDS),
                                 k) << "\n";
            std::cout << ss.str();
        } else {
            std::cout << result.to_string(config->anno_labels_delimiter,
                                          config->suppress_unlabeled,
                                          config->verbose_output
                                            || !(config->query_mode == COUNTS || config->query_mode == COORDS),
                                          k) + "\n";
        }
    };

    std::unique_ptr<graph::AnnotatedDBG> anno_graph;

    // iterate over input files
    for (const auto &file : files) {
        logger->trace("Parsing sequences from file '{}'", file);
        Timer curr_timer;
        size_t num_bp = 0;

        // If there is a single file and it's a batch query, run a special mode that
        // loads the index in in parts, only when it's needed for query.
        if (files.size() == 1 && config->query_batch_size && config->query_mode != COORDS) {
            if (graph->get_mode() == DeBruijnGraph::PRIMARY) {
                graph = std::make_shared<graph::CanonicalDBG>(graph);
                logger->trace("Primary graph wrapped into canonical");
            }

            ThreadPool annotation_loader(1, 1);
            std::shared_future<std::unique_ptr<AnnotatedDBG::Annotator>> annotation = annotation_loader.enqueue([&]() {
                auto anno = load_annotation(graph, *config);
                check_annotation(*config, anno->get_matrix());
                return anno;
            });

            num_bp = batched_query_fasta(file, query_callback, *config, aligner_config.get(),
                                         *graph, [&](){ return annotation.get().get(); }
            );
        } else {
            if (!anno_graph)
                anno_graph = initialize_annotated_dbg(graph, *config);
            num_bp = query_fasta(file, query_callback, *config, *anno_graph, aligner_config.get(), thread_pool);
        }

        auto time = curr_timer.elapsed();
        logger->trace("File '{}' with {} base pairs was processed in {} sec, throughput: {:.1f} bp/s",
                      file, num_bp, time, (double)num_bp / time);
    }

    return 0;
}


/**
 * Align a sequence to the annotated DBG before querying. Modify the sequence in place and return
 * alignment information (score and cigar string).
 *
 * @param seq               pointer to sequence string (which will be modified if alignment found)
 * @param graph             reference to DBG we align to
 * @param aligner_config    alignment config
 *
 * @return Alignment struct instance with score and cigar string
 */
Alignment align_sequence(std::string *seq,
                         const DeBruijnGraph &graph,
                         const align::DBGAlignerConfig &aligner_config) {
    align::DBGAligner aligner(graph, aligner_config);
    const align::DBGAlignerConfig &revised_config = aligner.get_config();
    align::DBGAlignerConfig::score_t max_score = revised_config.match_score(*seq)
        + revised_config.left_end_bonus + revised_config.right_end_bonus;
    auto alignments = aligner.align(*seq);

    assert(alignments.size() <= 1 && "Only the best alignment is needed");

    if (alignments.size()) {
        auto &match = alignments[0];
        // modify sequence for querying with the best alignment
        if (match.get_offset()) {
            *seq = graph.get_node_sequence(match.get_nodes()[0]).substr(0, match.get_offset())
                  + std::string(match.get_sequence());
        } else {
            *seq = match.get_sequence();
        }

        return { match.get_score(), max_score, match.get_cigar().to_string(), match.get_orientation() };
    } else {
        return { 0, max_score, fmt::format("{}S", seq->length()), false };
    }
}


SeqSearchResult query_sequence(QuerySequence&& sequence,
                               const AnnotatedDBG &anno_graph,
                               const Config &config,
                               const align::DBGAlignerConfig *aligner_config) {
    // Align sequence and modify sequence struct to store alignment results
    std::optional<Alignment> alignment;
    if (aligner_config)
        alignment = align_sequence(&sequence.sequence, anno_graph.get_graph(), *aligner_config);

    // Get sequence search result by executing query
    SeqSearchResult result = execute_query(std::move(sequence),
            config.query_mode,
            config.num_top_labels, config.discovery_fraction,
            config.presence_fraction, anno_graph);

    if (aligner_config)
        result.get_alignment() = alignment;

    return result;
}


// Check the the annotation matrix stores the data required for the query mode
void check_annotation(const Config &config, const annot::matrix::BinaryMatrix &anno_matrix) {
    // Only query_coords/count_kmers if using coord/count aware index.
    if (config.query_mode == COORDS
            && !dynamic_cast<const annot::matrix::MultiIntMatrix *>(&anno_matrix)) {
        logger->error("Annotation does not support k-mer coordinate queries. "
                      "First transform this annotation to include coordinate data "
                      "(e.g., {}, {}, {}, {}, {}).",
                      Config::annotype_to_string(Config::ColumnCoord),
                      Config::annotype_to_string(Config::BRWTCoord),
                      Config::annotype_to_string(Config::RowDiffCoord),
                      Config::annotype_to_string(Config::RowDiffDiskCoord),
                      Config::annotype_to_string(Config::RowDiffBRWTCoord));
        exit(1);
    }

    if ((config.query_mode == COUNTS || config.query_mode == COUNTS_SUM)
            && !dynamic_cast<const annot::matrix::IntMatrix *>(&anno_matrix)) {
        logger->error("Annotation does not support k-mer count queries. "
                      "First transform this annotation to include count data "
                      "(e.g., {}, {}, {}).",
                      Config::annotype_to_string(Config::IntBRWT),
                      Config::annotype_to_string(Config::IntRowDiffDisk),
                      Config::annotype_to_string(Config::IntRowDiffBRWT));
        exit(1);
    }
}

size_t query_fasta(const std::string &file,
                   const std::function<void(const SeqSearchResult &)> &callback,
                   const Config &config,
                   const graph::AnnotatedDBG &anno_graph,
                   const graph::align::DBGAlignerConfig *aligner_config,
                   ThreadPool &thread_pool) {
    check_annotation(config, anno_graph.get_annotator().get_matrix());

    if (config.query_batch_size) {
        if (config.query_mode != COORDS) {
            // Construct a query graph and query against it
            return batched_query_fasta(file, callback, config, aligner_config,
                                       anno_graph.get_graph(),
                                       [&]() { return &anno_graph.get_annotator(); });
        } else {
            // TODO: Implement batch mode for query_coords queries
            logger->warn("Querying coordinates in batch mode is currently not supported. Querying sequentially...");
        }
    }

    seq_io::FastaParser fasta_parser(file, config.forward_and_reverse);

    // Query sequences independently
    size_t seq_count = 0;
    size_t num_bp = 0;

    for (const seq_io::kseq_t &kseq : fasta_parser) {
        thread_pool.enqueue([&](QuerySequence &sequence) {
            // Callback with the SeqSearchResult
            callback(query_sequence(std::move(sequence), anno_graph,
                                    config, aligner_config));
        }, QuerySequence { seq_count++, std::string(kseq.name.s), std::string(kseq.seq.s) });
        num_bp += kseq.seq.l;
    }

    // wait while all threads finish processing the current file
    thread_pool.join();

    return num_bp;
}

size_t batched_query_fasta(const std::string &file,
                           const std::function<void(const SeqSearchResult &)> &callback,
                           const Config &config,
                           const graph::align::DBGAlignerConfig *aligner_config,
                           const DeBruijnGraph &graph,
                           const std::function<const AnnotatedDBG::Annotator *()> &get_annotation) {
    assert(graph.get_mode() != DeBruijnGraph::PRIMARY
            && "Primary graphs must be wrapped into canonical");

    seq_io::FastaParser fasta_parser(file, config.forward_and_reverse);

    auto it = fasta_parser.begin();
    auto end = fasta_parser.end();

    const uint64_t batch_size = config.query_batch_size;

    size_t seq_count = 0;
    size_t num_bp = 0;

    ThreadPool thread_pool(config.parallel_each);
    size_t threads_per_batch = get_num_threads() / config.parallel_each;

    std::mutex mu;
    std::vector<QuerySequence> seq_chunk;
    std::vector<std::pair<std::string, std::vector<node_index>>> contigs_chunk;
    size_t num_kmer_matches_chunk = 0;
    uint64_t num_bytes_read_chunk = 0;

    std::mutex query_mu;

    while (it != end) {
        uint64_t num_bytes_read = 0;

        // A generator that can be called multiple times until all sequences
        // are called
        std::vector<QuerySequence> seq_batch;

        for ( ; it != end && num_bytes_read <= batch_size; ++it) {
            seq_batch.push_back(QuerySequence { seq_count++, it->name.s, it->seq.s });
            num_bytes_read += it->seq.l;
        }

        thread_pool.enqueue([&](std::vector<QuerySequence> seq_batch, uint64_t num_bytes_read) {
            Timer batch_timer;
            std::vector<Alignment> alignments_batch;
            // Align sequences ahead of time on full graph if we don't have batch_align
            if (aligner_config && !config.batch_align) {
                alignments_batch.resize(seq_batch.size());
                logger->trace("Aligning sequences from batch against the full graph...");
                batch_timer.reset();

                #pragma omp parallel for num_threads(threads_per_batch) schedule(dynamic, 10)
                for (size_t i = 0; i < seq_batch.size(); ++i) {
                    // Set alignment for this seq_batch
                    alignments_batch[i] = align_sequence(&seq_batch[i].sequence, graph, *aligner_config);
                }
                logger->trace("Sequences alignment took {} sec", batch_timer.elapsed());
                batch_timer.reset();
            }

            // Construct the query graph for this batch
            size_t sub_k;
            auto [contigs, num_kmer_matches] = construct_contigs(
                graph, seq_batch, threads_per_batch,
                aligner_config && config.batch_align ? &config : NULL,
                &sub_k
            );

            // we accumulate batches until a certain number of k-mer matches to make query graph large enough
            {
                std::lock_guard<std::mutex> lock(mu);
                contigs_chunk.insert(contigs_chunk.end(),
                                     std::make_move_iterator(contigs.begin()),
                                     std::make_move_iterator(contigs.end()));
                contigs.clear();
                seq_chunk.insert(seq_chunk.end(),
                             std::make_move_iterator(seq_batch.begin()),
                             std::make_move_iterator(seq_batch.end()));
                seq_batch.clear();
                num_kmer_matches_chunk += num_kmer_matches;
                num_kmer_matches = 0;
                num_bytes_read_chunk += num_bytes_read;
                num_bytes_read = 0;
                // if the current batch is too small and there are more batches in the queue, go to next
                if (num_kmer_matches_chunk < batch_size * config.batch_min_matches
                        && thread_pool.num_waiting_tasks())
                    return;
                std::swap(contigs, contigs_chunk);
                std::swap(seq_batch, seq_chunk);
                std::swap(num_kmer_matches, num_kmer_matches_chunk);
                std::swap(num_bytes_read, num_bytes_read_chunk);
                // seq_chunk, contigs_chunk, and num_kmer_matches_chunk are now reset
            }
            size_t contigs_total_bp = 0;
            size_t contigs_total_kmers = 0;
            for (size_t i = 0; i < contigs.size(); ++i) {
                contigs_total_bp += contigs[i].first.length();
                contigs_total_kmers += contigs[i].second.size();
            }

            // work with seq_batch, contigs
            auto [small_graph, from_full_to_small] = construct_query_graph(std::move(contigs), threads_per_batch,
                                                                           graph.get_k(), graph.get_mode(), sub_k);

            std::lock_guard<std::mutex> query_lock(query_mu);

            auto query_graph = construct_query_graph(*get_annotation(), std::move(from_full_to_small), std::move(small_graph));

            auto query_graph_construction = batch_timer.elapsed();
            batch_timer.reset();

            #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic, 10)
            for (size_t i = 0; i < seq_batch.size(); ++i) {
                SeqSearchResult search_result
                    = query_sequence(std::move(seq_batch[i]), *query_graph, config,
                                     config.batch_align ? aligner_config : NULL);

                if (alignments_batch.size())
                    search_result.get_alignment() = std::move(alignments_batch[i]);

                callback(search_result);
            }

            logger->trace("Query graph constructed for {} bp from '{}' in {:.5f} sec, "
                          "query redundancy: {:.2f} bp/kmer, total bp/k-mers in contigs: {}/{}, queried in {:.5f} sec",
                          num_bytes_read, fasta_parser.get_filename(), query_graph_construction,
                          (double)num_bytes_read / query_graph->get_graph().num_nodes(),
                          contigs_total_bp, contigs_total_kmers,
                          batch_timer.elapsed());
        }, std::move(seq_batch), num_bytes_read);

        num_bp += num_bytes_read;
    }
    thread_pool.join();

    return num_bp;
}

} // namespace cli
} // namespace mtg
