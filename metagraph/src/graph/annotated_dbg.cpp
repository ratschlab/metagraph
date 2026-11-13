#include "annotated_dbg.hpp"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <span>

#include <fmt/format.h>

#include <tsl/hopscotch_set.h>
#include <tsl/hopscotch_map.h>

#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/binary_matrix/row_vector/unique_row_binmat.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "annotation/int_matrix/csr_matrix/csr_matrix.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/coord_to_header.hpp"
#include "common/algorithms.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vector_map.hpp"


namespace mtg {
namespace graph {

using mtg::annot::matrix::BinaryMatrix;
using mtg::annot::matrix::IntMatrix;
using mtg::annot::matrix::MultiIntMatrix;
using mtg::annot::matrix::CSRMatrix;
using mtg::annot::matrix::RowMajor;
using mtg::annot::matrix::UniqueRowBinmat;
using Column = mtg::annot::matrix::BinaryMatrix::Column;

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

    for (node_index i : graph_->node_range(sequence)) {
        if (i > 0)
            indices.push_back(graph_to_anno_index(i));
    }

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

        for (node_index i : graph_->node_range(data[t].first)) {
            if (i > 0)
                indices.push_back(graph_to_anno_index(i));
        }
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

    for (node_index i : graph_->node_range(sequence)) {
        // only insert indexes for matched k-mers and shift counts accordingly
        if (i > 0) {
            indices.push_back(graph_to_anno_index(i));
            kmer_counts[indices.size() - 1] = kmer_counts[end++];
        }
    }

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

        for (node_index i : graph_->node_range(sequence)) {
            if (i > 0) {
                ids[last].push_back(graph_to_anno_index(i));
                coords[last].emplace_back(graph_to_anno_index(i), coord);
            }
            coord++;
        }
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
        return num_kmers + 1;  // not enough k-mers in the graph -> return an unattainable threshold
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
        auto counts = get_top_labels(sequence, std::numeric_limits<size_t>::max(),
                                     discovery_fraction, presence_fraction);
        std::vector<Label> result;
        result.reserve(counts.size());
        for (auto &[label, count] : counts)
            result.emplace_back(std::move(label));
        return result;
    }

    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;

    VectorMap<row_index, size_t> index_counts;
    index_counts.reserve(sequence.size() - dbg_.get_k() + 1);

    size_t num_present_kmers = 0;

    for (node_index i : graph_->node_range(sequence)) {
        if (i > 0) {
            index_counts[graph_to_anno_index(i)]++;
            num_present_kmers++;
        }
    }

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
    assert(num_top_labels <= code_counts.size());
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

    // TODO: remove this step? (return (code, count) pairs and defer label decoding to the caller)
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
        auto counts = get_header_results<size_t>(sequence, num_top_labels,
                                                discovery_fraction, presence_fraction,
                                                with_kmer_counts);
        std::vector<StringCountPair> result;
        result.reserve(counts.size());
        for (auto &[label, matches, count] : counts)
            result.emplace_back(std::move(label), count);
        return result;
    }

    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;

    VectorMap<row_index, size_t> index_counts;
    index_counts.reserve(num_kmers);

    size_t num_present_kmers = 0;

    for (node_index i : graph_->node_range(sequence)) {
        if (i > 0) {
            index_counts[graph_to_anno_index(i)]++;
            num_present_kmers++;
        }
    }

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

    assert(with_kmer_counts || std::all_of(result.begin(), result.end(),
                                    [&](const auto &pair) { return pair.second <= num_kmers; }));

    return result;
}


namespace {
using Header = std::pair<Column, size_t>;
using HeaderCounts = std::vector<std::pair<Header, size_t>>;

// a query position whose k-mer is not in the graph
constexpr size_t kNoRow = std::numeric_limits<size_t>::max();

// Input coordinates are packed by CoordToHeader::map_to_local_coords. Preserve
// the existing per-row set iteration order, including untruncated result order.
template <bool WithAbundances = false>
HeaderCounts count_header_matches(const std::vector<MultiIntMatrix::RowTuples> &rows,
                                  const annot::CoordToHeader &mapping,
                                  size_t min_count, size_t num_top_labels,
                                  tsl::hopscotch_map<Header, size_t> *abundances = nullptr,
                                  std::span<const size_t> weights = {},
                                  size_t *num_passed = nullptr) {
    VectorMap<Header, size_t> counts_map;
    tsl::hopscotch_set<Header> matches;
    size_t row_index = 0;
    for (const auto &row : rows) {
        const size_t weight = weights.empty() ? 1 : weights[row_index++];
        matches.clear();
        for (const auto &[column, coords] : row) {
            for (uint64_t coord : coords) {
                Header h{column, coord % mapping.num_sequences(column)};
                matches.emplace(h);
                if constexpr (WithAbundances)
                    (*abundances)[h] += weight;
            }
        }
        for (const auto &h : matches)
            counts_map[h] += weight;
    }
    auto counts = to_vector(std::move(counts_map));
    std::erase_if(counts, [&](const auto &p) { return p.second < min_count; });
    if (num_passed)
        *num_passed = counts.size();
    if (counts.size() > num_top_labels)
        top_n_sorted(counts, num_top_labels);
    return counts;
}
} // namespace


template <class Result, class RowEnumerator>
Result filter_and_aggregate(const RowEnumerator &enumerate_rows,
                            const annot::LabelEncoder<Label> &label_encoder,
                            size_t min_count,
                            size_t num_top_labels,
                            const std::vector<size_t> &kmer_positions,
                            size_t num_kmers,
                            size_t *num_passed = nullptr) {
    static_assert(utils::is_one_of_v<Result,
                    std::vector<std::tuple<Label, size_t, PresenceRuns>>,
                    std::vector<std::tuple<Label, size_t, CountRuns>>,
                    std::vector<std::tuple<Label, size_t, PositionCoords>>>);

    auto call_bits = [&](const auto &callback) {
        enumerate_rows([&](const auto &row) {
            for (const auto &j : row) {
                callback(utils::get_first(j), 1);
            }
        });
    };
    std::vector<std::pair<Column, size_t>> counts
            = utils::accumulate_counts(call_bits, label_encoder.size(), min_count);

    if (num_passed)
        *num_passed = counts.size();
    if (counts.size() > num_top_labels)
        top_n_sorted(counts, num_top_labels);

    if (counts.empty())
        return {};

    // the rows come in the order of their query positions, in which the values are built
    using ValueType = std::tuple_element_t<2, typename Result::value_type>;
    utils::ValueStore<Column, ValueType> value_store(label_encoder.size(), counts.size());
    for (const auto &[j, count] : counts) {
        value_store.emplace(j, num_kmers);
    }

    size_t i = 0;
    enumerate_rows([&](const auto &row) {
        for (const auto &item : row) {
            auto j = utils::get_first(item);
            ValueType *value_p = value_store.get(j);
            if (value_p) {
                if constexpr(std::is_same_v<ValueType, PresenceRuns>) {
                    value_p->add(kmer_positions[i]);
                } else if constexpr(std::is_same_v<ValueType, CountRuns>) {
                    value_p->add(kmer_positions[i], item.second);
                } else {
                    value_p->add(kmer_positions[i], PositionCoords::Coords(item.second));
                }
            }
        }
        i++;
    });

    Result result;
    result.reserve(counts.size());
    for (const auto &[j, count] : counts) {
        ValueType *value = value_store.get(j);
        value->shrink();
        result.emplace_back(label_encoder.decode(j), count, std::move(*value));
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
    for (AnnotatedDBG::node_index i : graph.node_range(sequence)) {
        if (i > 0) {
            kmer_positions->push_back(j);
            row_indices->push_back(AnnotatedDBG::graph_to_anno_index(i));
        }
        j++;
    }
    assert(j == num_kmers);
}

template <class Values>
std::vector<std::tuple<Label, size_t, Values>>
AnnotatedDBG::get_header_results(std::string_view sequence, size_t num_top_labels,
                                 double discovery_fraction, double presence_fraction,
                                 bool with_kmer_counts, size_t *num_passed) const {
    static_assert(utils::is_one_of_v<Values, size_t, CountRuns, PresenceRuns>);
    const size_t num_kmers = sequence.size() - dbg_.get_k() + 1;
    // Repeated query positions still count separately, but their annotation
    // row only needs decoding once. Only positional results retain positions.
    using Positions = std::conditional_t<std::is_same_v<Values, size_t>,
                                         size_t, std::vector<size_t>>;
    std::vector<std::pair<row_index, Positions>> index_counts;
    size_t num_present = 0;
    {
        VectorMap<row_index, Positions> counts;
        counts.reserve(num_kmers);
        size_t position = 0;
        for (node_index node : dbg_.node_range(sequence)) {
            if (node) {
                auto &positions = counts[graph_to_anno_index(node)];
                if constexpr (std::is_same_v<Values, size_t>)
                    ++positions;
                else
                    positions.push_back(position);
                ++num_present;
            }
            ++position;
        }
        index_counts = to_vector(std::move(counts));
    } // release the hash table before decoding annotation rows
    const size_t min_count = get_min_count(discovery_fraction, presence_fraction,
                                           num_kmers, num_present);
    if (num_present < min_count)
        return {};
    std::vector<row_index> rows;
    std::vector<size_t> weights;
    rows.reserve(index_counts.size());
    weights.reserve(index_counts.size());
    for (const auto &[row, positions] : index_counts) {
        rows.push_back(row);
        if constexpr (std::is_same_v<Values, size_t>)
            weights.push_back(positions);
        else
            weights.push_back(positions.size());
    }
    // the index in |rows| of each query position (kNoRow where the k-mer is not in the graph),
    // to build the per-position values in the order of positions
    std::vector<size_t> row_at;
    if constexpr (!std::is_same_v<Values, size_t>) {
        row_at.assign(num_kmers, kNoRow);
        for (size_t i = 0; i < index_counts.size(); ++i) {
            for (size_t position : index_counts[i].second) {
                row_at[position] = i;
            }
        }
        index_counts = decltype(index_counts)();
    }
    const auto *matrix = dynamic_cast<const MultiIntMatrix *>(&annotator_->get_matrix());
    if (!matrix)
        throw std::runtime_error("k-mer coordinates are not indexed in this annotator");
    if (annotator_->num_labels() != coord_to_header_->num_columns()) {
        throw std::runtime_error(fmt::format("Incompatible number of columns in CoordToHeader "
                                             "and annotation matrix: {} != {}",
                                             coord_to_header_->num_columns(),
                                             annotator_->num_labels()));
    }
    auto tuples = matrix->get_row_tuples(rows);
    coord_to_header_->map_to_local_coords(&tuples);
    std::vector<std::tuple<Label, size_t, Values>> result;
    if constexpr (std::is_same_v<Values, size_t>) {
        tsl::hopscotch_map<Header, size_t> abundances;
        auto counts = with_kmer_counts
            ? count_header_matches<true>(tuples, *coord_to_header_, min_count, num_top_labels,
                                         &abundances, weights)
            : count_header_matches(tuples, *coord_to_header_, min_count, num_top_labels,
                                    nullptr, weights);
        result.reserve(counts.size());
        for (const auto &[h, count] : counts)
            result.emplace_back(coord_to_header_->get_headers(h.first)[h.second], count,
                                with_kmer_counts ? abundances.at(h) : count);
    } else {
        auto counts = count_header_matches(tuples, *coord_to_header_, min_count,
                                           num_top_labels, nullptr, weights, num_passed);
        tsl::hopscotch_map<Header, size_t> returned; // index in |counts|
        returned.reserve(counts.size());
        for (size_t r = 0; r < counts.size(); ++r)
            returned.emplace(counts[r].first, r);
        std::vector<Values> values(counts.size(), Values(num_kmers));

        // Collapse target occurrences to one count per header and distinct
        // query k-mer, then add it at that k-mer's query positions. The positions
        // are visited in ascending order, a stretch of one row in one step. A row's
        // headers are found at its first position and kept until its last one;
        // its decoded tuples are released then.
        std::vector<std::vector<std::pair<size_t, size_t>>> hits(tuples.size());
        std::vector<bool> found(tuples.size(), false);
        std::vector<size_t> &due = weights; // positions of each row not visited yet
        tsl::hopscotch_map<Header, size_t> occurrences;
        for (size_t position = 0; position < num_kmers; ) {
            const size_t i = row_at[position];
            if (i == kNoRow) {
                ++position;
                continue;
            }
            size_t end = position + 1;
            while (end < num_kmers && row_at[end] == i) {
                ++end;
            }
            if (!found[i]) {
                found[i] = true;
                occurrences.clear();
                for (const auto &[column, coords] : tuples[i])
                    for (uint64_t coord : coords)
                        ++occurrences[{column, coord % coord_to_header_->num_sequences(column)}];
                for (const auto &[h, count] : occurrences) {
                    auto it = returned.find(h);
                    if (it != returned.end())
                        hits[i].emplace_back(it->second, count);
                }
                tuples[i] = MultiIntMatrix::RowTuples();
            }
            for (const auto &[r, count] : hits[i]) {
                if constexpr (std::is_same_v<Values, PresenceRuns>)
                    values[r].add_run(position, end);
                else
                    values[r].add_run(position, end, count);
            }
            due[i] -= end - position;
            if (!due[i])
                hits[i] = std::vector<std::pair<size_t, size_t>>();
            position = end;
        }
        result.reserve(counts.size());
        for (size_t r = 0; r < counts.size(); ++r) {
            const auto &[h, count] = counts[r];
            values[r].shrink();
            result.emplace_back(coord_to_header_->get_headers(h.first)[h.second],
                                count, std::move(values[r]));
        }
    }
    return result;
}

#ifndef NDEBUG
// A sanity check to ensure the same matches are for get_top_label_signatures and get_top_labels
template <typename ValueType>
bool same_results(const std::vector<std::tuple<Label, size_t, ValueType>> &result,
                  const std::vector<StringCountPair> &top_labels) {
    assert(top_labels.size() == result.size());
    std::unordered_map<Label, uint64_t> check(top_labels.begin(), top_labels.end());
    for (const auto &[label, count, value] : result) {
        if constexpr(utils::is_one_of_v<ValueType, PresenceRuns, CountRuns, PositionCoords>) {
            assert(count == value.count());
        } else {
            static_assert(utils::dependent_false<ValueType>::value, "Unsupported value type");
        }
        auto find = check.find(label);
        assert(find != check.end());
        assert(find->second == count);
        check.erase(find);
    }
    return check.empty();
}
#endif // NDEBUG

std::vector<std::tuple<std::string, size_t, CountRuns>>
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
        return get_header_results<CountRuns>(sequence, num_top_labels,
                                             discovery_fraction, presence_fraction);
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
    if (!int_matrix)
        throw std::runtime_error("k-mer counts are not indexed in this annotator");

    auto get_results = [&](const auto &enumerate_rows) {
        auto results = filter_and_aggregate<std::vector<std::tuple<std::string, size_t, CountRuns>>>(
            enumerate_rows, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers);
        assert(same_results(results, get_top_labels(sequence, num_top_labels, discovery_fraction, presence_fraction)));
        return results;
    };

    // Per-index get_row_values is cheap for CSRMatrix; other matrices get rows in a batch.
    if (const auto *mat = dynamic_cast<const CSRMatrix *>(int_matrix)) {
        auto enumerate_row_values = [mat, &row_indices](const auto &callback) {
            for (row_index row : row_indices) {
                callback(mat->get_row_values(row));
            }
        };
        return get_results(enumerate_row_values);
    }

    // Get rows in a batch for other matrices.
    auto row_values = int_matrix->get_row_values(row_indices);
    auto enumerate_row_values = [rows = std::move(row_values)](const auto &callback) {
        std::for_each(rows.begin(), rows.end(), callback);
    };
    return get_results(enumerate_row_values);
}

std::vector<std::tuple<Label, size_t, PositionCoords, size_t>>
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
    if (!tuple_matrix)
        throw std::runtime_error("k-mer coordinates are not indexed in this annotator");

    auto rows_tuples = tuple_matrix->get_row_tuples(row_indices);

    if (coord_to_header_) {
        if (annotator_->num_labels() != coord_to_header_->num_columns()) {
            throw std::runtime_error(fmt::format("Incompatible number of columns in CoordToHeader "
                                                 "and annotation matrix: {} != {}",
                                                 coord_to_header_->num_columns(),
                                                 annotator_->num_labels()));
        }

        coord_to_header_->map_to_local_coords(&rows_tuples);
        // Now each coord in `rows_tuples` is `coord = local_coord * num_sequences[j] + seq_id`

        using Count = size_t;
        auto counts = count_header_matches(rows_tuples, *coord_to_header_,
                                           min_count, num_top_labels);

        // split coordinates across sequences, each held at the positions it matches
        // (the rows come in the order of their query positions)
        tsl::hopscotch_map<Header, size_t> returned; // index in |counts|
        returned.reserve(counts.size());
        for (size_t r = 0; r < counts.size(); ++r) {
            returned.emplace(counts[r].first, r);
        }
        std::vector<PositionCoords> coords_of(counts.size(), PositionCoords(num_kmers));
        for (size_t i = 0; i < rows_tuples.size(); ++i) {
            for (auto &[col, coords] : rows_tuples[i]) {
                for (uint64_t coord : coords) {
                    size_t header = coord % coord_to_header_->num_sequences(col);
                    uint64_t local_coord = coord / coord_to_header_->num_sequences(col);
                    auto it = returned.find({col, header});
                    if (it != returned.end())
                        coords_of[it->second].add(kmer_positions[i], local_coord);
                }
                coords.clear();
            }
            rows_tuples[i].clear();
        }
        rows_tuples.clear();

        std::vector<std::tuple<std::string, Count, PositionCoords, size_t>> result;
        result.reserve(counts.size());
        for (size_t r = 0; r < counts.size(); ++r) {
            const auto &[h, count] = counts[r];
            const auto &[col, header] = h;
            coords_of[r].shrink();
            result.emplace_back(coord_to_header_->get_headers(col)[header],
                                count, std::move(coords_of[r]),
                                coord_to_header_->num_kmers_in_sequence(col, header));
        }

        return result;
    }

    auto call_rows_tuples = [&](const auto &callback) {
        std::for_each(rows_tuples.begin(), rows_tuples.end(), callback);
    };

    auto results = filter_and_aggregate<std::vector<std::tuple<Label, size_t, PositionCoords>>>(
        call_rows_tuples, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers);
    assert(same_results(results, get_top_labels(sequence, num_top_labels, discovery_fraction, presence_fraction)));

    // Without a CoordToHeader mapping, the label refers to a whole annotation
    // column rather than a single indexed sequence, so num_kmers_in_target is
    // undefined.
    std::vector<std::tuple<Label, size_t, PositionCoords, size_t>> with_lengths;
    with_lengths.reserve(results.size());
    for (auto &[label, count, coords] : results) {
        with_lengths.emplace_back(std::move(label), count, std::move(coords), 0);
    }
    return with_lengths;
}

std::vector<std::tuple<Label, size_t, PresenceRuns>>
AnnotatedDBG::get_top_label_signatures(std::string_view sequence,
                                       size_t num_top_labels,
                                       double discovery_fraction,
                                       double presence_fraction,
                                       size_t *num_passed,
                                       bool header_labels) const {
    assert(discovery_fraction >= 0.);
    assert(discovery_fraction <= 1.);
    assert(presence_fraction >= 0.);
    assert(presence_fraction <= 1.);
    assert(check_compatibility());

    if (num_passed)
        *num_passed = 0;

    if (sequence.size() < dbg_.get_k())
        return {};

    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;

    // the shortcut ranks with get_top_labels, which reports no count and ranks headers when
    // a CoordToHeader mapping is loaded
    if (discovery_fraction == 1. && !num_passed && (header_labels || !coord_to_header_)) {
        auto label_counts = get_top_labels(sequence, num_top_labels, discovery_fraction, presence_fraction);
        std::vector<std::tuple<Label, size_t, PresenceRuns>> presence_vectors;
        presence_vectors.reserve(label_counts.size());
        // every position marked: one run
        PresenceRuns all(num_kmers);
        all.add_run(0, num_kmers);
        all.shrink();
        for (auto&& [label, count] : label_counts) {
            presence_vectors.emplace_back(std::move(label), count, all);
        }
        return presence_vectors;
    }

    if (coord_to_header_ && header_labels) {
        return get_header_results<PresenceRuns>(sequence, num_top_labels,
                                                discovery_fraction, presence_fraction,
                                                false, num_passed);
    }

    // kmers in the graph and their positions in the query sequence
    std::vector<row_index> row_indices;
    std::vector<size_t> kmer_positions;  // positions of k-mers in the query sequence
    map_to_rows(dbg_, sequence, &row_indices, &kmer_positions);

    size_t min_count = get_min_count(discovery_fraction, presence_fraction,
                                     num_kmers, kmer_positions.size());
    if (row_indices.size() < min_count)
        return {};

    const auto &matrix = annotator_->get_matrix();
    auto get_results = [&](const auto &enumerate_rows) {
        auto results = filter_and_aggregate<std::vector<std::tuple<std::string, size_t, PresenceRuns>>>(
            enumerate_rows, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers,
            num_passed);
        // get_top_labels ranks headers when a CoordToHeader mapping is loaded
        assert(coord_to_header_ || same_results(results, get_top_labels(sequence, num_top_labels,
                                                                         discovery_fraction, presence_fraction)));
        return results;
    };

    // Per-index get_row is cheap for these RowMajor types; other matrices get rows in a batch.
    if (const auto *mat = dynamic_cast<const UniqueRowBinmat *>(&matrix)) {
        auto enumerate_rows = [mat, &row_indices](const auto &callback) {
            for (row_index row : row_indices) {
                callback(mat->get_row_ref(row));
            }
        };
        return get_results(enumerate_rows);
    } else if (const auto *mat = dynamic_cast<const CSRMatrix *>(&matrix)) {
        auto enumerate_rows = [mat, &row_indices](const auto &callback) {
            for (row_index row : row_indices) {
                callback(mat->get_row_values(row));
            }
        };
        return get_results(enumerate_rows);
    }
    auto rows = matrix.get_rows(row_indices);
    auto enumerate_rows = [rows = std::move(rows)](const auto &callback) {
        std::for_each(rows.begin(), rows.end(), callback);
    };
    return get_results(enumerate_rows);
}

bool AnnotatedSequenceGraph::label_exists(const Label &label) const {
    return annotator_->get_label_encoder().label_exists(label);
}

bool AnnotatedSequenceGraph::check_compatibility() const {
    // TODO: what if CanonicalDBG is not the highest level? find a better way to do this
    if (const auto *canonical = dynamic_cast<const CanonicalDBG *>(graph_.get()))
        return canonical->get_graph().max_index() == annotator_->num_objects();

    return graph_->max_index() == annotator_->num_objects();
}

std::string encode_presence_mask(const sdsl::bit_vector &mask) {
    std::string result;
    size_t pos = 0;
    while (pos < mask.size()) {
        bool bit = mask[pos];
        size_t run_start = pos;
        while (pos < mask.size() && (bool)mask[pos] == bit) ++pos;
        result += (bit ? 'x' : 'o');
        result += std::to_string(pos - run_start);
    }
    return result;
}

std::string encode_presence_mask(const PresenceRuns &mask) {
    std::string result;
    size_t pos = 0;
    mask.call_runs([&](size_t begin, size_t end) {
        assert(begin < end && end <= mask.size());
        assert(begin > pos || !begin);
        if (begin > pos) {
            result += 'o';
            result += std::to_string(begin - pos);
        }
        result += 'x';
        result += std::to_string(end - begin);
        pos = end;
    });
    if (pos < mask.size()) {
        result += 'o';
        result += std::to_string(mask.size() - pos);
    }
    return result;
}


namespace {

// a run of 32-bit positions ends at most here: a larger query is held densely from the start
constexpr size_t kMaxRunsSize = std::numeric_limits<uint32_t>::max();

// Make room for one more element in |v| within |budget| bytes: double its capacity,
// up to what |budget| holds. Return false if |budget| holds no more elements.
template <typename T>
bool grow_within(std::vector<T> &v, size_t budget) {
    const size_t max_capacity = budget / sizeof(T);
    if (v.size() >= max_capacity)
        return false;
    v.reserve(std::min(max_capacity, std::max<size_t>(2, 2 * v.capacity())));
    return true;
}

} // namespace

PresenceRuns::PresenceRuns(size_t size) : size_(size) {
    if (size_ > kMaxRunsSize)
        to_dense_form();
}

PresenceRuns::PresenceRuns(const sdsl::bit_vector &mask) : PresenceRuns(mask.size()) {
    for (size_t i = 0; i < mask.size(); ++i) {
        if (mask[i])
            add(i);
    }
}

bool PresenceRuns::grow() {
    return grow_within(runs_, dense_bytes(size_));
}

void PresenceRuns::to_dense_form() {
    dense_ = sdsl::bit_vector(size_, false);
    dense_form_ = true;
    for (const Run &run : runs_) {
        mark(run.begin, run.end);
    }
    std::vector<Run>().swap(runs_);
}

void PresenceRuns::mark(size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) {
        dense_[i] = true;
    }
}

size_t PresenceRuns::count() const {
    if (dense_form_)
        return sdsl::util::cnt_one_bits(dense_);
    size_t count = 0;
    for (const Run &run : runs_) {
        count += run.end - run.begin;
    }
    return count;
}

sdsl::bit_vector PresenceRuns::to_dense() const {
    if (dense_form_)
        return dense_;
    sdsl::bit_vector mask(size_, false);
    for (const Run &run : runs_) {
        for (size_t i = run.begin; i < run.end; ++i) {
            mask[i] = true;
        }
    }
    return mask;
}

size_t PresenceRuns::heap_bytes() const {
    return dense_form_ ? dense_bytes(size_) : runs_.capacity() * sizeof(Run);
}

bool PresenceRuns::operator==(const PresenceRuns &other) const {
    return size_ == other.size_ && to_dense() == other.to_dense();
}

CountRuns::CountRuns(size_t size) : size_(size) {
    if (size_ > kMaxRunsSize)
        to_dense_form();
}

CountRuns::CountRuns(const std::vector<size_t> &values) : CountRuns(values.size()) {
    for (size_t i = 0; i < values.size(); ++i) {
        add(i, values[i]);
    }
}

bool CountRuns::grow() {
    return grow_within(runs_, dense_bytes(size_));
}

void CountRuns::to_dense_form() {
    dense_.assign(size_, 0);
    dense_form_ = true;
    for (const Run &run : runs_) {
        std::fill(dense_.begin() + run.begin, dense_.begin() + run.end, run.value);
    }
    std::vector<Run>().swap(runs_);
}

size_t CountRuns::count() const {
    if (dense_form_)
        return size_ - std::count(dense_.begin(), dense_.end(), size_t(0));
    size_t count = 0;
    for (const Run &run : runs_) {
        count += run.end - run.begin;
    }
    return count;
}

std::vector<size_t> CountRuns::to_dense() const {
    if (dense_form_)
        return dense_;
    std::vector<size_t> values(size_, 0);
    for (const Run &run : runs_) {
        std::fill(values.begin() + run.begin, values.begin() + run.end, run.value);
    }
    return values;
}

size_t CountRuns::heap_bytes() const {
    return dense_form_ ? dense_.capacity() * sizeof(size_t) : runs_.capacity() * sizeof(Run);
}

bool CountRuns::operator==(const CountRuns &other) const {
    return size_ == other.size_ && to_dense() == other.to_dense();
}

PositionCoords::PositionCoords(size_t size) : size_(size) {
    if (size_ > kMaxRunsSize)
        to_dense_form();
}

PositionCoords::PositionCoords(const std::vector<Coords> &values) : PositionCoords(values.size()) {
    for (size_t i = 0; i < values.size(); ++i) {
        add(i, Coords(values[i]));
    }
}

bool PositionCoords::open(size_t position) {
    const bool new_run = runs_.empty() || runs_.back().end != position;
    // the runs and the coordinate lists share the budget: each grows within what
    // the other's capacity leaves
    const size_t budget = dense_bytes(size_);
    if ((coords_.size() == coords_.capacity()
            && !grow_within(coords_, budget - runs_.capacity() * sizeof(Run)))
        || (new_run && runs_.size() == runs_.capacity()
            && !grow_within(runs_, budget - coords_.capacity() * sizeof(Coords)))) {
        to_dense_form();
        return false;
    }
    if (new_run) {
        runs_.push_back({ static_cast<uint32_t>(position), static_cast<uint32_t>(position + 1) });
    } else {
        ++runs_.back().end;
    }
    return true;
}

void PositionCoords::to_dense_form() {
    dense_.resize(size_);
    dense_form_ = true;
    size_t c = 0;
    for (const Run &run : runs_) {
        for (size_t i = run.begin; i < run.end; ++i) {
            dense_[i] = std::move(coords_[c++]);
        }
    }
    assert(c == coords_.size());
    std::vector<Run>().swap(runs_);
    std::vector<Coords>().swap(coords_);
}

size_t PositionCoords::count() const {
    if (!dense_form_)
        return coords_.size();
    return std::count_if(dense_.begin(), dense_.end(),
                         [](const Coords &coords) { return !coords.empty(); });
}

std::vector<PositionCoords::Coords> PositionCoords::to_dense() const {
    if (dense_form_)
        return dense_;
    std::vector<Coords> values(size_);
    call_positions([&](size_t i, const Coords &coords) { values[i] = coords; });
    return values;
}

size_t PositionCoords::heap_bytes() const {
    return dense_form_ ? dense_.capacity() * sizeof(Coords)
                       : runs_.capacity() * sizeof(Run) + coords_.capacity() * sizeof(Coords);
}

bool PositionCoords::operator==(const PositionCoords &other) const {
    return size_ == other.size_ && to_dense() == other.to_dense();
}

} // namespace graph
} // namespace mtg
