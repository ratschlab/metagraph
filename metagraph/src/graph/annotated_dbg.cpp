#include "annotated_dbg.hpp"

#include <array>
#include <cstdlib>

#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "common/utils/simd_utils.hpp"
#include "common/aligned_vector.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vector_map.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;
using mtg::annot::matrix::IntMatrix;
using mtg::annot::matrix::MultiIntMatrix;

typedef AnnotatedDBG::Label Label;
typedef std::pair<Label, size_t> StringCountPair;


AnnotatedSequenceGraph
::AnnotatedSequenceGraph(std::shared_ptr<SequenceGraph> graph,
                         std::shared_ptr<Annotator> annotation,
                         bool force_fast)
      : graph_(graph), annotator_(annotation),
        force_fast_(force_fast) {
    assert(graph_.get());
    assert(annotator_.get());
}

void AnnotatedSequenceGraph
::annotate_sequence(std::string_view sequence,
                    const std::vector<Label> &labels) {
    assert(check_compatibility());

    std::vector<row_index> indices;
    indices.reserve(sequence.size());

    graph_->map_to_nodes(sequence, [&](node_index i) {
        if (i > 0)
            indices.push_back(graph_to_anno_index(i));
    });

    if (!indices.size())
        return;

    std::lock_guard<std::mutex> lock(mutex_);

    if (force_fast_) {
        auto row_major = dynamic_cast<annot::RowCompressed<Label>*>(annotator_.get());
        if (row_major) {
            row_major->add_labels_fast(indices, labels);
            return;
        }
    }

    annotator_->add_labels(indices, labels);
}

void AnnotatedSequenceGraph
::annotate_sequences(const std::vector<std::pair<std::string, std::vector<Label>>> &data) {
    assert(check_compatibility());

    std::vector<std::vector<row_index>> ids(data.size());
    size_t last = 0;
    for (size_t t = 0; t < data.size(); ++t) {
        // if the labels are the same, write indexes to the same array
        auto &indices = data[t].second == data[last].second ? ids[last] : ids[t];
        indices.reserve(data[t].first.size());

        graph_->map_to_nodes(data[t].first, [&](node_index i) {
            if (i > 0)
                indices.push_back(graph_to_anno_index(i));
        });
    }

    std::lock_guard<std::mutex> lock(mutex_);

    for (size_t t = 0; t < data.size(); ++t) {
        if (!ids[t].size())
            continue;

        if (force_fast_) {
            auto row_major = dynamic_cast<annot::RowCompressed<Label>*>(annotator_.get());
            if (row_major) {
                row_major->add_labels_fast(ids[t], data[t].second);
                continue;
            }
        }

        annotator_->add_labels(ids[t], data[t].second);
    }
}

void AnnotatedDBG::add_kmer_counts(std::string_view sequence,
                                   const std::vector<Label> &labels,
                                   std::vector<uint64_t>&& kmer_counts) {
    assert(check_compatibility());
    assert(kmer_counts.size() == sequence.size() - dbg_.get_k() + 1);

    if (sequence.size() < dbg_.get_k())
        return;

    std::vector<row_index> indices;
    indices.reserve(sequence.size() - dbg_.get_k() + 1);
    size_t end = 0;

    graph_->map_to_nodes(sequence, [&](node_index i) {
        // only insert indexes for matched k-mers and shift counts accordingly
        if (i > 0) {
            indices.push_back(graph_to_anno_index(i));
            kmer_counts[indices.size() - 1] = kmer_counts[end++];
        }
    });

    kmer_counts.resize(end);

    if (indices.size())
        annotator_->add_label_counts(indices, labels, kmer_counts);
}

void AnnotatedDBG::add_kmer_coord(std::string_view sequence,
                                  const std::vector<Label> &labels,
                                  uint64_t coord) {
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return;

    std::vector<row_index> indices;
    indices.reserve(sequence.size() - dbg_.get_k() + 1);

    graph_->map_to_nodes(sequence, [&](node_index i) { indices.push_back(i); });

    if (!indices.size())
        return;

    std::lock_guard<std::mutex> lock(mutex_);

    for (node_index i : indices) {
        // only insert coordinates for matched k-mers and increment the coordinates
        if (i > 0)
            annotator_->add_label_coord(graph_to_anno_index(i), labels, coord);
        coord++;
    }
}

void AnnotatedDBG::add_kmer_coords(
        const std::vector<std::tuple<std::string, std::vector<Label>, uint64_t>> &data) {
    assert(check_compatibility());

    std::vector<std::vector<row_index>> ids(data.size());
    for (size_t t = 0; t < data.size(); ++t) {
        const auto &[sequence, labels, _] = data[t];
        if (sequence.size() < dbg_.get_k())
            continue;

        auto &indices = ids[t];
        indices.reserve(sequence.size() - dbg_.get_k() + 1);

        graph_->map_to_nodes(sequence, [&](node_index i) { indices.push_back(i); });
    }

    std::lock_guard<std::mutex> lock(mutex_);

    for (size_t t = 0; t < data.size(); ++t) {
        const auto &labels = std::get<1>(data[t]);
        uint64_t coord = std::get<2>(data[t]);

        for (node_index i : ids[t]) {
            // only insert coordinates for matched k-mers and increment the coordinates
            if (i > 0)
                annotator_->add_label_coord(graph_to_anno_index(i), labels, coord);
            coord++;
        }
    }
}

std::vector<Label> AnnotatedDBG::get_labels(std::string_view sequence,
                                            double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    VectorMap<row_index, size_t> index_counts;
    index_counts.reserve(sequence.size() - dbg_.get_k() + 1);

    size_t num_present_kmers = 0;
    size_t num_missing_kmers = 0;

    graph_->map_to_nodes(sequence, [&](node_index i) {
        if (i > 0) {
            index_counts[graph_to_anno_index(i)]++;
            num_present_kmers++;
        } else {
            num_missing_kmers++;
        }
    });

    size_t min_count = std::max(1.0, std::ceil(presence_ratio
                                                 * (num_present_kmers
                                                     + num_missing_kmers)));

    if (num_present_kmers < min_count)
        return {};

    return get_labels(index_counts.values_container(), min_count);
}

std::vector<Label>
AnnotatedDBG::get_labels(const std::vector<std::pair<row_index, size_t>> &index_counts,
                         size_t min_count) const {
    assert(check_compatibility());

    auto code_counts = annotator_->get_matrix().sum_rows(
        index_counts,
        min_count,
        std::max(min_count, size_t(1))
    );

    std::vector<Label> labels;
    labels.reserve(code_counts.size());

    const auto &label_encoder = annotator_->get_label_encoder();

    for (const auto &pair : code_counts) {
        assert(pair.second >= min_count);
        labels.push_back(label_encoder.decode(pair.first));
    }

    return labels;
}

std::vector<Label>
AnnotatedSequenceGraph::get_labels(node_index index) const {
    assert(check_compatibility());
    assert(index != SequenceGraph::npos);

    return annotator_->get(graph_to_anno_index(index));
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(std::string_view sequence,
                             size_t num_top_labels,
                             double presence_ratio,
                             bool with_kmer_counts) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    VectorMap<row_index, size_t> index_counts;
    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;
    index_counts.reserve(num_kmers);

    size_t num_present_kmers = 0;

    graph_->map_to_nodes(sequence, [&](node_index i) {
        if (i > 0) {
            index_counts[graph_to_anno_index(i)]++;
            num_present_kmers++;
        }
    });

    uint64_t min_count = std::max(1.0, std::ceil(presence_ratio * num_kmers));
    if (num_present_kmers < min_count)
        return {};

    auto top_labels = get_top_labels(index_counts.values_container(),
                                     num_top_labels, min_count, with_kmer_counts);

    assert(std::all_of(top_labels.begin(), top_labels.end(),
                       [&](const auto &pair) { return pair.second <= num_kmers; })
            || with_kmer_counts);

    return top_labels;
}

std::vector<std::pair<std::string, std::vector<size_t>>>
AnnotatedDBG::get_label_count_quantiles(std::string_view sequence,
                                        size_t num_top_labels,
                                        double presence_ratio,
                                        const std::vector<double> &count_quantiles) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());
    if (!std::is_sorted(count_quantiles.begin(), count_quantiles.end()))
        throw std::runtime_error("Quantiles must be sorted");
    if (count_quantiles.at(0) < 0. || count_quantiles.back() > 1.)
        throw std::runtime_error("Quantiles must be in range [0, 1]");

    if (sequence.size() < dbg_.get_k())
        return {};

    std::vector<row_index> rows;
    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;
    rows.reserve(num_kmers);

    graph_->map_to_nodes(sequence, [&](node_index i) {
        if (i > 0)
            rows.push_back(graph_to_anno_index(i));
    });

    uint64_t min_count = std::max(1.0, std::ceil(presence_ratio * num_kmers));
    if (rows.size() < min_count)
        return {};

    std::vector<size_t> q_low(count_quantiles.size());
    for (size_t i = 0; i < count_quantiles.size(); ++i) {
        q_low[i] = (num_kmers - 1) * count_quantiles[i];
    }

    const auto *int_matrix = dynamic_cast<const IntMatrix *>(&annotator_->get_matrix());
    if (!int_matrix) {
        logger->error("k-mer counts are not indexed in this annotator");
        exit(1);
    }

    VectorMap<size_t, std::vector<uint64_t>> code_to_counts;
    for (const auto &row_values : int_matrix->get_row_values(rows)) {
        for (const auto &[column, count] : row_values) {
            code_to_counts[column].push_back(count);
        }
    }

    std::vector<std::pair<size_t, std::vector<uint64_t>>> code_counts;
    code_counts.reserve(code_to_counts.size());
    for (auto &[j, counts] : code_to_counts.values_container()) {
        // filter by the number of matched k-mers
        if (counts.size() >= min_count)
            code_counts.emplace_back(j, std::move(counts));
    }
    // sort by the number of matched k-mers
    std::sort(code_counts.begin(), code_counts.end(),
              [](const auto &x, const auto &y) {
                  return x.second.size() > y.second.size()
                      || (x.second.size() == y.second.size() && x.first < y.first);
              });
    // keep only the first |num_top_labels| top labels
    if (code_counts.size() > num_top_labels)
        code_counts.resize(num_top_labels);

    std::vector<std::pair<Label, std::vector<size_t>>> label_quantiles;
    label_quantiles.reserve(code_counts.size());
    // Quantiles are defined as `count[i]` where `i < q * N <= i + 1`
    for (auto &[j, counts] : code_counts) {
        std::sort(counts.begin(), counts.end());
        const size_t num_zeros = num_kmers - counts.size();

        label_quantiles.emplace_back(annotator_->get_label_encoder().decode(j),
                                     std::vector<size_t>(q_low.size()));

        std::vector<size_t> &quantiles = label_quantiles.back().second;
        for (size_t q = 0; q < q_low.size(); ++q) {
            if (q_low[q] < num_zeros) {
                quantiles[q] = 0;
            } else {
                quantiles[q] = counts[q_low[q] - num_zeros];
            }
        }
    }

    return label_quantiles;
}

std::vector<std::pair<std::string, std::vector<SmallVector<uint64_t>>>>
AnnotatedDBG::get_kmer_coordinates(std::string_view sequence,
                                   size_t num_top_labels,
                                   double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    std::vector<node_index> path;
    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;
    path.reserve(num_kmers);

    graph_->map_to_nodes(sequence, [&](node_index i) {
        path.push_back(i);
    });

    return get_kmer_coordinates(path, num_top_labels, presence_ratio);
}

std::vector<std::pair<std::string, std::vector<SmallVector<uint64_t>>>>
AnnotatedDBG::get_kmer_coordinates(const std::vector<node_index> &path,
                                   size_t num_top_labels,
                                   double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    if (!path.size())
        return {};

    std::vector<row_index> rows;
    rows.reserve(path.size());

    std::vector<size_t> ids;
    ids.reserve(path.size());

    for (node_index i : path) {
        if (i > 0) {
            ids.push_back(rows.size());
            rows.push_back(graph_to_anno_index(i));
        } else {
            ids.push_back(-1);
        }
    }

    uint64_t min_count = std::max(1.0, std::ceil(presence_ratio * path.size()));
    if (rows.size() < min_count)
        return {};

    const auto *tuple_matrix = dynamic_cast<const MultiIntMatrix *>(&annotator_->get_matrix());
    if (!tuple_matrix) {
        logger->error("k-mer coordinates are not indexed in this annotator");
        exit(1);
    }

    auto rows_tuples = tuple_matrix->get_row_tuples(rows);

    VectorMap<size_t, size_t> code_to_count;
    for (const auto &row_tuples : rows_tuples) {
        for (const auto &[column, tuple] : row_tuples) {
            code_to_count[column] += 1;
        }
    }

    auto code_counts = code_to_count.values_container();
    // sort by the number of matched k-mers
    std::sort(code_counts.begin(), code_counts.end(),
              [](const auto &x, const auto &y) {
                  return x.second > y.second || (x.second == y.second && x.first < y.first);
              });

    // keep only the first |num_top_labels| top labels
    if (code_counts.size() > num_top_labels)
        code_counts.resize(num_top_labels);

    // filter by the number of matched k-mers
    code_counts.erase(
        std::upper_bound(code_counts.begin(), code_counts.end(), min_count,
                         [](uint64_t min, const auto &x) { return x.second < min; }),
        code_counts.end()
    );

    code_to_count = VectorMap<size_t, size_t>(code_counts.begin(), code_counts.end());

    std::vector<std::pair<std::string, std::vector<SmallVector<uint64_t>>>> result(code_to_count.size());

    for (size_t j = 0; j < result.size(); ++j) {
        result[j].first = annotator_->get_label_encoder().decode(code_counts[j].first);
    }

    for (size_t i : ids) {
        for (size_t j = 0; j < result.size(); ++j) {
            // append empty tuple
            result[j].second.emplace_back();
        }

        // leave all tuples empty if the k-mer is missing
        if (i == (size_t)-1)
            continue;

        // set the non-empty tuples
        for (auto &[j, tuple] : rows_tuples[i]) {
            auto it = code_to_count.find(j);
            if (it != code_to_count.end())
                result[it - code_to_count.begin()].second.back() = std::move(tuple);
        }
    }

    return result;
}

std::vector<std::pair<Label, sdsl::bit_vector>>
AnnotatedDBG::get_top_label_signatures(std::string_view sequence,
                                       size_t num_top_labels,
                                       double presence_ratio) const {
    assert(presence_ratio >= 0.);
    assert(presence_ratio <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;

    if (presence_ratio == 1.) {
        std::vector<std::pair<Label, sdsl::bit_vector>> presence_vectors;

        auto label_counts = get_top_labels(sequence, num_top_labels, presence_ratio);
        presence_vectors.reserve(label_counts.size());
        for (auto&& [label, count] : label_counts) {
            presence_vectors.emplace_back(
                std::move(label),
                sdsl::bit_vector(num_kmers, true)
            );
        }
        return presence_vectors;
    }

    // kmers and their positions in the query sequence
    std::vector<row_index> row_indices;
    row_indices.reserve(num_kmers);

    std::vector<size_t> kmer_positions;
    kmer_positions.reserve(num_kmers);

    size_t j = 0;
    graph_->map_to_nodes(sequence, [&](node_index i) {
        if (i > 0) {
            kmer_positions.push_back(j);
            row_indices.push_back(graph_to_anno_index(i));
        }
        j++;
    });
    assert(j == num_kmers);

    const uint64_t min_count = std::max(1.0, std::ceil(presence_ratio * num_kmers));

    if (kmer_positions.size() < min_count)
        return {};

    // k-mer presence mask with the k-mer count
    using SignatureCount = std::pair<std::vector<uint8_t>, size_t>;

    typedef uint64_t LabelCode;
    // map each label code to a k-mer presence mask and its popcount
    VectorMap<LabelCode, SignatureCount> label_codes_to_presence;

    auto label_codes = annotator_->get_matrix().get_rows(row_indices);

    assert(label_codes.size() == row_indices.size());

    for (size_t i = 0; i < row_indices.size(); ++i) {
        for (size_t label_code : label_codes[i]) {
            auto& [mask, label_count] = label_codes_to_presence[label_code];

            if (mask.empty())
                mask.resize(num_kmers, 0);

            mask[kmer_positions[i]] = true;
            label_count++;
        }
    }

    // get label codes with k-mer match signatures and label counts
    auto &vector = const_cast<std::vector<std::pair<LabelCode, SignatureCount>>&>(
        label_codes_to_presence.values_container()
    );

    // sort to get top |num_top_labels| labels
    if (vector.size() > num_top_labels) {
        std::sort(vector.begin(), vector.end(),
                  [](const auto &a, const auto &b) {
                      return a.second.second > b.second.second
                          || (a.second.second == b.second.second && a.first < b.first);
                  });
        vector.resize(num_top_labels);
    }

    std::vector<std::pair<Label, sdsl::bit_vector>> results;
    results.reserve(vector.size());

    for (const auto &[code, mask_count] : vector) {
        // TODO: check this before sorting
        if (mask_count.second < min_count)
            continue;

        // TODO: remove the decoding step?
        results.emplace_back(annotator_->get_label_encoder().decode(code),
                             to_sdsl(mask_count.first));

        assert(sdsl::util::cnt_one_bits(results.back().second)
                == mask_count.second);
        assert(results.back().second.size() == num_kmers);
    }

#ifndef NDEBUG
    // sanity check, make sure that the same matches are output by get_top_labels
    auto top_labels = get_top_labels(sequence, num_top_labels, presence_ratio);
    assert(top_labels.size() == results.size());

    std::unordered_map<Label, uint64_t> check(top_labels.begin(), top_labels.end());
    for (const auto &[label, mask] : results) {
        auto find = check.find(label);
        assert(find != check.end());
        assert(find->second == sdsl::util::cnt_one_bits(mask));
        check.erase(find);
    }
    assert(check.empty());
#endif // NDEBUG

    return results;
}

template <class Container>
std::vector<StringCountPair> top_labels(Container&& code_counts,
                                        const annot::LabelEncoder<> &label_encoder,
                                        size_t num_top_labels,
                                        size_t min_count) {
    assert(std::all_of(
        code_counts.begin(), code_counts.end(),
        [&](const auto &code_count) { return code_count.second >= min_count; }
    ));
    std::ignore = min_count;

    if (code_counts.size() > num_top_labels) {
        // sort labels by counts to get the top |num_top_labels|
        std::sort(code_counts.begin(), code_counts.end(),
                  [](const auto &first, const auto &second) {
                      return first.second > second.second
                          || (first.second == second.second && first.first < second.first);
                  });
        // leave only the first |num_top_labels| top labels
        code_counts.resize(num_top_labels);
    }

    // TODO: remove this step?
    std::vector<StringCountPair> label_counts(code_counts.size());
    for (size_t i = 0; i < code_counts.size(); ++i) {
        label_counts[i].first = label_encoder.decode(code_counts[i].first);
        label_counts[i].second = code_counts[i].second;
    }

    return label_counts;
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(const std::vector<std::pair<row_index, size_t>> &index_counts,
                             size_t num_top_labels,
                             size_t min_count,
                             bool with_kmer_counts) const {
    assert(check_compatibility());

    if (with_kmer_counts) {
        // TODO: Don't take into account counts when comparing to min_count.
        //       It should be compared to sum_rows and not sum_row_values.
        return top_labels(dynamic_cast<const IntMatrix &>(annotator_->get_matrix())
                                                  .sum_row_values(index_counts, min_count),
                          annotator_->get_label_encoder(),
                          num_top_labels, min_count);
    } else {
        return top_labels(annotator_->get_matrix().sum_rows(index_counts, min_count),
                          annotator_->get_label_encoder(),
                          num_top_labels, min_count);
    }
}

bool AnnotatedSequenceGraph::label_exists(const Label &label) const {
    return annotator_->label_exists(label);
}

bool AnnotatedSequenceGraph::has_label(node_index index, const Label &label) const {
    assert(check_compatibility());
    assert(index != SequenceGraph::npos);

    return annotator_->has_label(graph_to_anno_index(index), label);
}

void AnnotatedSequenceGraph
::call_annotated_nodes(const Label &label,
                       std::function<void(node_index)> callback) const {
    assert(check_compatibility());

    annotator_->call_objects(
        label,
        [&](row_index index) { callback(anno_to_graph_index(index)); }
    );
}

bool AnnotatedSequenceGraph::check_compatibility() const {
    return graph_->max_index() == annotator_->num_objects();
}

bool AnnotatedDBG::check_compatibility() const {
    return dbg_.get_base_graph().max_index() == annotator_->num_objects();
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

#ifdef __AVX2__
__m256d get_penalty_bigsi_avx2(__m256d counts,
                               __m256d match_score,
                               __m256d mismatch_score,
                               __m256d SNP_t) {
    __m256d neg_ones = _mm256_set1_pd(-1.0);

    __m256d min_N_snps = _mm256_div_pd(counts, SNP_t);

    // subtract -1 instead of add 1 to reuse a register
    __m256d max_N_snps = _mm256_max_pd(
        _mm256_sub_pd(_mm256_sub_pd(counts, SNP_t), neg_ones),
        min_N_snps
    );

    // separate mul and add is faster than using FMA instructions
    __m256d mean_N_snps = fmafast_pd(max_N_snps, _mm256_set1_pd(0.05), min_N_snps);

    assert(!_mm256_movemask_pd(_mm256_cmp_pd(mean_N_snps, counts, 14)));

    __m256d mean_penalty = _mm256_mul_pd(mean_N_snps, mismatch_score);

    return _mm256_mul_pd(fmafast_pd(_mm256_sub_pd(mean_penalty, counts),
                                    match_score,
                                    mean_penalty),
                         neg_ones);
}
#endif // __AVX2__


/**
 * A penalty function used in BIGSI
 */
double penalty_bigsi(double count,
                     double match_score,
                     double mismatch_score,
                     double SNP_t) {
    double min_N_snps = count / SNP_t;
    double max_N_snps = std::max(count - SNP_t + 1, min_N_snps);
    double mean_N_snps = max_N_snps * 0.05 + min_N_snps;

    assert(count >= mean_N_snps);

    double mean_penalty = mean_N_snps * mismatch_score;
    return (count - mean_penalty) * match_score - mean_penalty;
}

int32_t AnnotatedDBG
::score_kmer_presence_mask(const sdsl::bit_vector &kmer_presence_mask,
                           int32_t match_score,
                           int32_t mismatch_score) const {
    if (!kmer_presence_mask.size())
        return 0;

    const size_t k = dbg_.get_k();
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

    const auto *it = score_counter[0].data();
    const auto *end = it + score_counter[0].size();

#ifdef __AVX2__
    __m256d match_score_packed = _mm256_set1_pd(match_score);
    __m256d mismatch_score_packed = _mm256_set1_pd(mismatch_score);
    __m256d SNP_t_packed = _mm256_set1_pd(SNP_t);
    __m256d penalties = _mm256_setzero_pd();

    for ( ; it + 4 <= end; it += 4) {
        __m256d penalty_add = get_penalty_bigsi_avx2(
            uint64_to_double(_mm256_load_si256(reinterpret_cast<const __m256i*>(it))),
            match_score_packed,
            mismatch_score_packed,
            SNP_t_packed
        );

        // TODO: the least-significant bits are usually off by one, is there
        //       a way to fix this?
        assert(float(penalty_bigsi(*it, match_score, mismatch_score, SNP_t)
                    + penalty_bigsi(*(it + 1), match_score, mismatch_score, SNP_t)
                    + penalty_bigsi(*(it + 2), match_score, mismatch_score, SNP_t)
                    + penalty_bigsi(*(it + 3), match_score, mismatch_score, SNP_t))
            == float(haddall_pd(penalty_add)));

        penalties = _mm256_add_pd(penalties, penalty_add);
    }

    // reduce
    score += haddall_pd(penalties);
#endif // __AVX2__

    for ( ; it != end; ++it) {
        score += penalty_bigsi(*it, match_score, mismatch_score, SNP_t);
    }

    return std::max(score * sequence_length / kmer_presence_mask.size(), 0.);
}

} // namespace graph
} // namespace mtg
