#include "annotated_dbg.hpp"

#include <array>
#include <cstdlib>

#include <tsl/hopscotch_set.h>
#include <tsl/hopscotch_map.h>

#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/coord_to_header.hpp"
#include "common/algorithms.hpp"
#include "common/aligned_vector.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vector_map.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;
using mtg::annot::matrix::IntMatrix;
using mtg::annot::matrix::MultiIntMatrix;
using Column = mtg::annot::matrix::BinaryMatrix::Column;
using Tuple = mtg::annot::matrix::MultiIntMatrix::Tuple;

typedef AnnotatedDBG::Label Label;
typedef std::pair<Label, size_t> StringCountPair;


AnnotatedSequenceGraph
::AnnotatedSequenceGraph(std::shared_ptr<SequenceGraph> graph,
                         std::unique_ptr<Annotator>&& annotation,
                         bool force_fast)
      : graph_(graph), annotator_(std::move(annotation)), force_fast_(force_fast) {
    assert(graph_.get());
    assert(annotator_.get());
    assert(check_compatibility());
}

AnnotatedDBG::AnnotatedDBG(std::shared_ptr<DeBruijnGraph> dbg,
                           std::unique_ptr<Annotator>&& annotation,
                           bool force_fast,
                           std::unique_ptr<annot::CoordToHeader> coord_to_header)
      : AnnotatedSequenceGraph(dbg, std::move(annotation), force_fast), dbg_(*dbg),
        coord_to_header_(std::move(coord_to_header)) {}

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
        if (!indices.capacity())
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

    std::vector<row_index> indices = map_to_nodes(dbg_, sequence);

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

    std::vector<std::vector<row_index>> ids;
    ids.reserve(data.size());
    for (const auto &[sequence, labels, _] : data) {
        if (sequence.size() >= dbg_.get_k())
            ids.push_back(map_to_nodes(dbg_, sequence));
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

void AnnotatedDBG::annotate_kmer_coords(
        const std::vector<std::tuple<std::string, std::vector<Label>, uint64_t>> &data) {
    assert(check_compatibility());

    std::vector<std::vector<row_index>> ids(data.size());
    std::vector<std::vector<std::pair<row_index, uint64_t>>> coords(data.size());
    size_t last = 0;

    for (size_t t = 0; t < data.size(); ++t) {
        const auto &[sequence, labels, _] = data[t];
        uint64_t coord = std::get<2>(data[t]);
        if (sequence.size() < dbg_.get_k())
            continue;

        // if the labels are the same, write indexes to the same array
        if (labels != std::get<1>(data[last]))
            last = t;

        if (!ids[last].size()) {
            ids[last].reserve(sequence.size() - dbg_.get_k() + 1);
            coords[last].reserve(sequence.size() - dbg_.get_k() + 1);
        }

        graph_->map_to_nodes(sequence, [&](node_index i) {
            if (i > 0) {
                ids[last].push_back(graph_to_anno_index(i));
                coords[last].emplace_back(graph_to_anno_index(i), coord);
            }
            coord++;
        });
    }

    std::lock_guard<std::mutex> lock(mutex_);

    for (size_t t = 0; t < data.size(); ++t) {
        if (ids[t].size()) {
            const auto &labels = std::get<1>(data[t]);
            annotator_->add_labels(ids[t], labels);
            annotator_->add_label_coords(coords[t], labels);
        }
    }
}

size_t get_min_count(double discovery_fraction, double presence_fraction,
                     size_t num_kmers, size_t num_present) {
    if (num_present < std::max(1.0, std::ceil(presence_fraction * num_kmers)))
        return num_kmers + 1;  // too few k-mers are in the graph -> set the threshold too high
    return std::max(1.0, std::ceil(discovery_fraction * num_kmers));
}

std::vector<Label> AnnotatedDBG::get_labels(std::string_view sequence,
                                            double discovery_fraction,
                                            double presence_fraction) const {
    assert(discovery_fraction >= 0.);
    assert(discovery_fraction <= 1.);
    assert(presence_fraction >= 0.);
    assert(presence_fraction <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    if (coord_to_header_) {
        auto kmer_coord_res = get_kmer_coordinates(sequence, std::numeric_limits<size_t>::max(),
                                                   discovery_fraction, presence_fraction);
        std::vector<Label> result;
        result.reserve(kmer_coord_res.size());
        for (auto &[label, count, coords] : kmer_coord_res) {
            result.emplace_back(std::move(label));
        }
        return result;
    }

    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;

    VectorMap<row_index, size_t> index_counts;
    index_counts.reserve(sequence.size() - dbg_.get_k() + 1);

    size_t num_present_kmers = 0;

    graph_->map_to_nodes(sequence, [&](node_index i) {
        if (i > 0) {
            index_counts[graph_to_anno_index(i)]++;
            num_present_kmers++;
        }
    });

    size_t min_count = get_min_count(discovery_fraction, presence_fraction,
                                     num_kmers, num_present_kmers);
    if (num_present_kmers < min_count)
        return {};

    auto code_counts = annotator_->get_matrix().sum_rows(index_counts.values_container(), min_count);

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

    return annotator_->get_labels(graph_to_anno_index(index));
}

// Container `code_counts` must have pairs (index, count)
template <class Container>
void top_n_sorted(Container& code_counts, size_t num_top_labels) {
    auto comp = [](const auto &x, const auto &y) {
        return std::make_pair(y.second, x.first)
              < std::make_pair(x.second, y.first);
    };
    std::nth_element(code_counts.begin(),
                     code_counts.begin() + num_top_labels,
                     code_counts.end(),
                     comp);
    // leave only the first |num_top_labels| top labels
    code_counts.resize(num_top_labels);
    // sort the top labels by counts
    std::sort(code_counts.begin(), code_counts.end(), comp);
}

template <class Container>
std::vector<StringCountPair> filter_and_decode(Container&& code_counts,
                                               const annot::LabelEncoder<> &label_encoder,
                                               size_t num_top_labels) {
    if (code_counts.size() > num_top_labels)
        top_n_sorted(code_counts, num_top_labels);

    // TODO: remove this step?
    std::vector<StringCountPair> label_counts;
    label_counts.reserve(code_counts.size());
    for (const auto &[j, count] : code_counts) {
        label_counts.emplace_back(label_encoder.decode(j), count);
    }

    return label_counts;
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(std::string_view sequence,
                             size_t num_top_labels,
                             double discovery_fraction,
                             double presence_fraction,
                             bool with_kmer_counts) const {
    assert(discovery_fraction >= 0.);
    assert(discovery_fraction <= 1.);
    assert(presence_fraction >= 0.);
    assert(presence_fraction <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    if (coord_to_header_) {
        auto kmer_coord_res = get_kmer_coordinates(sequence, num_top_labels,
                                                   discovery_fraction, presence_fraction);
        std::vector<StringCountPair> result;
        result.reserve(kmer_coord_res.size());
        for (auto &[label, count, coords] : kmer_coord_res) {
            if (with_kmer_counts) {
                count = 0;
                for (const auto &tuple : coords) {
                    count += tuple.size();
                }
            }
            result.emplace_back(std::move(label), count);
        }
        return result;
    }

    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;

    VectorMap<row_index, size_t> index_counts;
    index_counts.reserve(num_kmers);

    size_t num_present_kmers = 0;

    graph_->map_to_nodes(sequence, [&](node_index i) {
        if (i > 0) {
            index_counts[graph_to_anno_index(i)]++;
            num_present_kmers++;
        }
    });

    size_t min_count = get_min_count(discovery_fraction, presence_fraction,
                                     num_kmers, num_present_kmers);
    if (num_present_kmers < min_count)
        return {};

    std::vector<StringCountPair> result;
    if (with_kmer_counts) {
        result = filter_and_decode(dynamic_cast<const IntMatrix &>(annotator_->get_matrix())
                                        .sum_row_values(index_counts.values_container(), min_count),
                                   annotator_->get_label_encoder(), num_top_labels);
    } else {
        result = filter_and_decode(annotator_->get_matrix()
                                        .sum_rows(index_counts.values_container(), min_count),
                                   annotator_->get_label_encoder(), num_top_labels);
    }

    assert(std::all_of(result.begin(), result.end(),
                       [&](const auto &pair) { return pair.second <= num_kmers; })
            || with_kmer_counts);

    return result;
}


template <class Result, class Container>
Result filter_and_aggregate(Container&& rows,
                            const annot::LabelEncoder<Label> &label_encoder,
                            size_t min_count,
                            size_t num_top_labels,
                            const std::vector<size_t> &kmer_positions,
                            size_t num_kmers) {
    static_assert(utils::is_one_of_v<Result,
                    std::vector<std::tuple<Label, size_t, sdsl::bit_vector>>,
                    std::vector<std::tuple<Label, size_t, std::vector<size_t>>>,
                    std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>>>);
    using ValueType = std::tuple_element_t<2, typename Result::value_type>;

    auto call_bits = [&](const std::function<void(Column, size_t)> &callback) {
        for (const auto &row : rows) {
            for (const auto &j : row) {
                callback(utils::get_first(j), 1);
            }
        }
    };
    std::vector<std::pair<Column, size_t>> counts
            = utils::accumulate_counts(call_bits, label_encoder.size(), min_count);

    if (counts.size() > num_top_labels)
        top_n_sorted(counts, num_top_labels);

    // Aggregate results (group by labels)
    tsl::hopscotch_map<Column, ValueType> values;
    values.reserve(counts.size());
    for (const auto &[j, count] : counts) {
        values.emplace(j, num_kmers);
    }

    for (size_t i = 0; i < rows.size(); ++i) {
        for (auto &item : rows[i]) {
            auto j = utils::get_first(item);
            auto it = values.find(j);
            if (it != values.end()) {
                if constexpr(std::is_same_v<ValueType, sdsl::bit_vector>) {
                    it.value()[kmer_positions[i]] = true;
                } else {
                    it.value()[kmer_positions[i]] = std::move(item.second);
                }
            }
        }
        rows[i].clear();
    }
    rows.clear();

    Result result;
    result.reserve(counts.size());
    for (const auto &[j, count] : counts) {
        result.emplace_back(label_encoder.decode(j), count, std::move(values[j]));
    }

    return result;
}

void map_to_rows(const DeBruijnGraph &graph, std::string_view sequence,
                 std::vector<AnnotatedDBG::row_index> *row_indices,
                 std::vector<size_t> *kmer_positions) {
    assert(sequence.size() >= graph.get_k());

    size_t num_kmers = sequence.size() - graph.get_k() + 1;
    row_indices->reserve(num_kmers);
    kmer_positions->reserve(num_kmers);
    size_t j = 0;
    graph.map_to_nodes(sequence, [&](AnnotatedDBG::node_index i) {
        if (i > 0) {
            kmer_positions->push_back(j);
            row_indices->push_back(AnnotatedDBG::graph_to_anno_index(i));
        }
        j++;
    });
    assert(j == num_kmers);
}

std::vector<std::tuple<std::string, size_t, std::vector<size_t>>>
AnnotatedDBG::get_kmer_counts(std::string_view sequence,
                              size_t num_top_labels,
                              double discovery_fraction,
                              double presence_fraction) const {
    assert(discovery_fraction >= 0.);
    assert(discovery_fraction <= 1.);
    assert(presence_fraction >= 0.);
    assert(presence_fraction <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    if (coord_to_header_) {
        auto kmer_coord_res = get_kmer_coordinates(sequence, num_top_labels,
                                                   discovery_fraction, presence_fraction);
        std::vector<std::tuple<std::string, size_t, std::vector<size_t>>> result;
        result.reserve(kmer_coord_res.size());
        for (auto &[label, count, coords] : kmer_coord_res) {
            result.emplace_back(std::move(label), count, coords.size());
            auto &counts = std::get<2>(result.back());
            for (size_t i = 0; i < coords.size(); ++i) {
                counts[i] = coords[i].size();
            }
        }
        return result;
    }

    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;

    // kmers in the graph and their positions in the query sequence
    std::vector<row_index> row_indices;
    std::vector<size_t> kmer_positions;  // positions of k-mers in the query sequence
    map_to_rows(dbg_, sequence, &row_indices, &kmer_positions);

    size_t min_count = get_min_count(discovery_fraction, presence_fraction,
                                     num_kmers, kmer_positions.size());
    if (kmer_positions.size() < min_count)
        return {};

    const auto *int_matrix = dynamic_cast<const IntMatrix *>(&annotator_->get_matrix());
    if (!int_matrix) {
        logger->error("k-mer counts are not indexed in this annotator");
        exit(1);
    }

    auto row_values = int_matrix->get_row_values(row_indices);

    return filter_and_aggregate<std::vector<std::tuple<std::string, size_t, std::vector<size_t>>>>(
        row_values, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers);
}

std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>>
AnnotatedDBG::get_kmer_coordinates(std::string_view sequence,
                                   size_t num_top_labels,
                                   double discovery_fraction,
                                   double presence_fraction) const {
    assert(discovery_fraction >= 0.);
    assert(discovery_fraction <= 1.);
    assert(presence_fraction >= 0.);
    assert(presence_fraction <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;

    // kmers in the graph and their positions in the query sequence
    std::vector<row_index> row_indices;
    std::vector<size_t> kmer_positions;  // positions of k-mers in the query sequence
    map_to_rows(dbg_, sequence, &row_indices, &kmer_positions);

    size_t min_count = get_min_count(discovery_fraction, presence_fraction,
                                     num_kmers, kmer_positions.size());
    if (kmer_positions.size() < min_count)
        return {};

    const auto *tuple_matrix = dynamic_cast<const MultiIntMatrix *>(&annotator_->get_matrix());
    if (!tuple_matrix) {
        logger->error("k-mer coordinates are not indexed in this annotator");
        exit(1);
    }

    auto rows_tuples = tuple_matrix->get_row_tuples(row_indices);

    if (coord_to_header_) {
        if (annotator_->num_labels() != coord_to_header_->num_columns()) {
            logger->error("Incompatible number of columns in CoordToHeader and annotation matrix: "
                          "{} != {}", coord_to_header_->num_columns(), annotator_->num_labels());
            exit(1);
        }

        coord_to_header_->map_to_local_coords(&rows_tuples);
        // Now each coord in `rows_tuples` is `coord = local_coord * num_headers[j] + header_id`

        // Before splitting coords across headers, filter them by the number of k-mer matches
        using Header = std::pair<Column, size_t>; // global header id
        using Count = size_t;
        std::vector<std::pair<Header, Count>> counts;
        {
            VectorMap<Header, Count> counts_map;
            tsl::hopscotch_set<Header> matches;
            for (const auto &row : rows_tuples) {
                matches.clear();
                for (auto &[col, coords] : row) {
                    for (uint64_t coord : coords) {
                        size_t header = coord % coord_to_header_->num_headers(col);
                        matches.emplace(col, header);
                    }
                }
                for (const Header &h : matches) {
                    counts_map[h]++;
                }
            }
            counts = to_vector(std::move(counts_map));
        }

        // remove headers with insufficient number of k-mer matches
        counts.erase(std::remove_if(counts.begin(), counts.end(),
                                    [&](const auto &pair) { return pair.second < min_count; }),
                     counts.end());

        // keep only top-n headers
        if (counts.size() > num_top_labels)
            top_n_sorted(counts, num_top_labels);

        // split coordinates across headers
        tsl::hopscotch_map<Header, std::vector<Tuple>> coords_map;
        coords_map.reserve(counts.size());
        for (const auto &[h, count] : counts) {
            coords_map.emplace(h, num_kmers);
        }
        for (size_t i = 0; i < rows_tuples.size(); ++i) {
            for (auto &[col, coords] : rows_tuples[i]) {
                for (uint64_t coord : coords) {
                    size_t header = coord % coord_to_header_->num_headers(col);
                    uint64_t local_coord = coord / coord_to_header_->num_headers(col);
                    auto it = coords_map.find({col, header});
                    if (it != coords_map.end())
                        it.value()[kmer_positions[i]].push_back(local_coord);
                }
                coords.clear();
            }
            rows_tuples[i].clear();
        }
        rows_tuples.clear();

        std::vector<std::tuple<std::string, Count, std::vector<Tuple>>> result;
        result.reserve(counts.size());
        for (const auto &[h, count] : counts) {
            const auto &[col, header] = h;
            auto &coords = coords_map[h];
            result.emplace_back(coord_to_header_->get_headers(col)[header],
                                count, std::move(coords));
        }

        return result;
    }

    return filter_and_aggregate<std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>>>(
        rows_tuples, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers);
}

// A sanity check to ensure the same matches are for get_top_label_signatures and get_top_labels
bool same_results(const std::vector<std::tuple<Label, size_t, sdsl::bit_vector>> &result,
                  const std::vector<StringCountPair> &top_labels) {
    assert(top_labels.size() == result.size());
    std::unordered_map<Label, uint64_t> check(top_labels.begin(), top_labels.end());
    for (const auto &[label, count, mask] : result) {
        assert(count == sdsl::util::cnt_one_bits(mask));
        auto find = check.find(label);
        assert(find != check.end());
        assert(find->second == sdsl::util::cnt_one_bits(mask));
        check.erase(find);
    }
    return check.empty();
}

std::vector<std::tuple<Label, size_t, sdsl::bit_vector>>
AnnotatedDBG::get_top_label_signatures(std::string_view sequence,
                                       size_t num_top_labels,
                                       double discovery_fraction,
                                       double presence_fraction) const {
    assert(discovery_fraction >= 0.);
    assert(discovery_fraction <= 1.);
    assert(presence_fraction >= 0.);
    assert(presence_fraction <= 1.);
    assert(check_compatibility());

    if (sequence.size() < dbg_.get_k())
        return {};

    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;

    if (discovery_fraction == 1.) {
        auto label_counts = get_top_labels(sequence, num_top_labels, discovery_fraction, presence_fraction);
        std::vector<std::tuple<Label, size_t, sdsl::bit_vector>> presence_vectors;
        presence_vectors.reserve(label_counts.size());
        for (auto&& [label, count] : label_counts) {
            presence_vectors.emplace_back(
                std::move(label),
                count,
                sdsl::bit_vector(num_kmers, true)
            );
        }
        return presence_vectors;
    }

    if (coord_to_header_) {
        auto kmer_coord_res = get_kmer_coordinates(sequence, num_top_labels,
                                                   discovery_fraction, presence_fraction);
        std::vector<std::tuple<Label, size_t, sdsl::bit_vector>> result;
        result.reserve(kmer_coord_res.size());
        for (auto &[label, count, coords] : kmer_coord_res) {
            assert(coords.size() == sequence.size() - dbg_.get_k() + 1);
            result.emplace_back(std::move(label), count, coords.size());
            auto &mask = std::get<2>(result.back());
            for (size_t i = 0; i < coords.size(); ++i) {
                if (coords[i].size())
                    mask[i] = true;
            }
        }
        assert(same_results(result, get_top_labels(sequence, num_top_labels, discovery_fraction, presence_fraction)));
        return result;
    }

    // kmers in the graph and their positions in the query sequence
    std::vector<row_index> row_indices;
    std::vector<size_t> kmer_positions;  // positions of k-mers in the query sequence
    map_to_rows(dbg_, sequence, &row_indices, &kmer_positions);

    size_t min_count = get_min_count(discovery_fraction, presence_fraction,
                                     num_kmers, kmer_positions.size());
    if (row_indices.size() < min_count)
        return {};

    auto rows = annotator_->get_matrix().get_rows(row_indices);

    auto result = filter_and_aggregate<std::vector<std::tuple<std::string, size_t, sdsl::bit_vector>>>(
        rows, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers);

    assert(same_results(result, get_top_labels(sequence, num_top_labels, discovery_fraction, presence_fraction)));

    return result;
}

bool AnnotatedSequenceGraph::label_exists(const Label &label) const {
    return annotator_->get_label_encoder().label_exists(label);
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
    // TODO: what if CanonicalDBG is not the highest level? find a better way to do this
    if (const auto *canonical = dynamic_cast<const CanonicalDBG *>(graph_.get()))
        return canonical->get_graph().max_index() == annotator_->num_objects();

    return graph_->max_index() == annotator_->num_objects();
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

} // namespace graph
} // namespace mtg
