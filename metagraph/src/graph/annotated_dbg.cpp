#include "annotated_dbg.hpp"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>

#include <tsl/hopscotch_set.h>
#include <tsl/hopscotch_map.h>

#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/binary_matrix/row_vector/unique_row_binmat.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "annotation/int_matrix/csr_matrix/csr_matrix.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "annotation/coord_to_header.hpp"
#include "common/algorithms.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vector_map.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;
using mtg::annot::matrix::BinaryMatrix;
using mtg::annot::matrix::IntMatrix;
using mtg::annot::matrix::MultiIntMatrix;
using mtg::annot::matrix::CSRMatrix;
using mtg::annot::matrix::RowMajor;
using mtg::annot::matrix::UniqueRowBinmat;
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
                           std::shared_ptr<const annot::CoordToHeader> coord_to_header)
      : AnnotatedSequenceGraph(dbg, std::move(annotation), force_fast), dbg_(*dbg),
        coord_to_header_(std::move(coord_to_header)),
        kmer_abundances_(holds_kmer_abundances(*annotator_)) {
    // a mapping of another column count maps no coordinate right: refused here, so that no
    // query has to (the loaders refuse it first, naming the files: check_record_mapping)
    if (coord_to_header_ && coord_to_header_->num_columns() != annotator_->num_labels()) {
        throw std::invalid_argument(fmt::format(
                "A record mapping (.seqs) of {} columns for an annotation of {} labels",
                coord_to_header_->num_columns(), annotator_->num_labels()));
    }
}

bool AnnotatedDBG::holds_kmer_abundances(const Annotator &annotation) {
    const auto &matrix = annotation.get_matrix();
    return dynamic_cast<const IntMatrix *>(&matrix)
            && !dynamic_cast<const MultiIntMatrix *>(&matrix);
}

const IntMatrix* AnnotatedDBG::abundance_matrix() const {
    return kmer_abundances_ ? dynamic_cast<const IntMatrix *>(&annotator_->get_matrix())
                            : nullptr;
}

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

uint64_t lower_median(std::vector<uint64_t> values) {
    values.erase(std::remove(values.begin(), values.end(), 0u), values.end());
    if (values.empty())
        return 0;
    auto median = values.begin() + (values.size() - 1) / 2;
    std::nth_element(values.begin(), median, values.end());
    return *median;
}

bool search_row_before(const SearchRow &a, const SearchRow &b) {
    if (a.score != b.score)
        return a.score > b.score;
    if (a.median != b.median)
        return a.median > b.median;
    return a.label < b.label;
}

std::vector<std::pair<size_t, std::string>>
search_top_n(const std::vector<size_t> &scores,
             size_t num_top_labels,
             const std::function<const std::string &(size_t)> &label,
             const std::function<std::vector<uint64_t>(const std::vector<size_t> &)> &medians,
             SearchCut *cut) {
    SearchCut c;
    c.abundances = static_cast<bool>(medians);

    // (a label is the caller's, never copied but for the rows returned: a tie group at a cut
    // can hold millions of rows, all of them ordered by label)
    struct Row {
        size_t index;
        size_t score;
        uint64_t median;
        std::string_view label;
    };
    std::vector<size_t> order(scores.size());
    std::iota(order.begin(), order.end(), 0);

    // [0, above): the rows above the N-th score, every one returned; [above, rows.size()):
    // the rows of the N-th score (the tie group at the cut), of which |slots| are returned
    size_t above = order.size();
    size_t slots = 0;
    if (order.size() > num_top_labels) {
        c.complete = false;
        if (!num_top_labels) {
            // no row returned, so no N-th row
            order.clear();
            above = 0;
        } else {
            std::nth_element(order.begin(), order.begin() + num_top_labels - 1, order.end(),
                             [&](size_t i, size_t j) { return scores[i] > scores[j]; });
            c.nth_score = scores[order[num_top_labels - 1]];
            above = std::partition(order.begin(), order.end(),
                                   [&](size_t i) { return scores[i] > c.nth_score; })
                    - order.begin();
            order.erase(std::partition(order.begin() + above, order.end(),
                                       [&](size_t i) { return scores[i] == c.nth_score; }),
                        order.end());
            slots = num_top_labels - above;
            c.tied = order.size() - above;
            c.tied_returned = slots;
        }
    }

    // the rows above the N-th score (every row of a complete result), each returned; the tie
    // group at a cut, [above, order.size()) of |order|, is held as its candidates' indexes
    // (and medians where they are read) until its |slots| places are taken, so that a tie
    // group of millions costs a few words per row, as the selection before release 10 did
    std::vector<Row> rows;
    rows.reserve(above + slots);
    for (size_t p = 0; p < above; ++p) {
        rows.push_back(Row{ order[p], scores[order[p]], 0, std::string_view() });
    }
    struct Tied {
        uint64_t median;
        size_t index;
    };
    // (filled where the medians of the tie group are read; without them the group stays in
    // |order|)
    std::vector<Tied> tie;
    const size_t tied = order.size() - above;

    // The medians break ties only, so a row's median is read only where it is compared: the
    // rows that share their score with another one, the tie group at a cut (whose N-th row's
    // median the cut then states) and the rows above it of one score. A row with a score of
    // its own keeps median 0, which no comparison reads, and a cut whose N-th score has one
    // row states no median. Reading the abundances costs a pass over the query's rows and
    // holds values of the columns asked, so it is asked for no other row, and once.
    if (medians && (rows.size() || tied)) {
        std::vector<size_t> asked;
        // the rows above the N-th score in score order: the runs of one score
        std::vector<size_t> by_score(above);
        std::iota(by_score.begin(), by_score.end(), 0);
        std::sort(by_score.begin(), by_score.end(), [&](size_t r, size_t s) {
            return rows[r].score < rows[s].score;
        });
        for (size_t b = 0; b < by_score.size(); ) {
            size_t e = b + 1;
            while (e < by_score.size() && rows[by_score[e]].score == rows[by_score[b]].score) {
                ++e;
            }
            if (e - b > 1)
                asked.insert(asked.end(), by_score.begin() + b, by_score.begin() + e);
            b = e;
        }
        std::vector<size_t> indexes(asked.size());
        std::transform(asked.begin(), asked.end(), indexes.begin(),
                       [&](size_t r) { return rows[r].index; });
        // the tie group at the cut, if the N-th score is shared
        if (tied > 1) {
            indexes.insert(indexes.end(), order.begin() + above, order.end());
            c.nth_median_read = true;
        }
        if (indexes.size()) {
            std::vector<uint64_t> values = medians(indexes);
            assert(values.size() == indexes.size());
            for (size_t a = 0; a < asked.size(); ++a) {
                rows[asked[a]].median = values[a];
            }
            if (c.nth_median_read) {
                tie.reserve(tied);
                for (size_t t = 0; t < tied; ++t) {
                    tie.push_back(Tied{ values[asked.size() + t], order[above + t] });
                }
            }
        }
    }

    // search_row_before's order, with one comparison of the labels, and candidates of one label
    // (a header in two columns) in their order in |scores|
    auto before = [](const Row &a, const Row &b) {
        if (a.score != b.score)
            return a.score > b.score;
        if (a.median != b.median)
            return a.median > b.median;
        if (const int labels = a.label.compare(b.label))
            return labels < 0;
        return a.index < b.index;
    };

    if (!c.complete && slots) {
        // The tie group competes for |slots| places: by median where its medians were read,
        // then by label. Only the rows of the N-th row's median compete by label, so only
        // they are decoded beside the rows returned (a large tie group, a common k-mer set at
        // the cut, is not decoded entirely where the medians separate it), and the places
        // are taken by a selection, not a sort of the whole group: O(T + slots log slots),
        // the returned rows sorted below.
        // the candidates of the N-th row's median: the whole group, or (medians read) those
        // left once the |slots|-th largest median is found and the group split around it
        std::vector<std::pair<const std::string *, size_t>> same;
        if (c.nth_median_read) {
            std::nth_element(tie.begin(), tie.begin() + slots - 1, tie.end(),
                             [](const Tied &a, const Tied &b) { return a.median > b.median; });
            c.nth_median = tie[slots - 1].median;
            auto same_begin = std::partition(tie.begin(), tie.end(),
                    [&](const Tied &t) { return t.median > c.nth_median; });
            auto same_end = std::partition(same_begin, tie.end(),
                    [&](const Tied &t) { return t.median == c.nth_median; });
            // every row of a larger median is returned
            for (auto it = tie.begin(); it != same_begin; ++it) {
                rows.push_back(Row{ it->index, c.nth_score, it->median, label(it->index) });
            }
            same.reserve(same_end - same_begin);
            for (auto it = same_begin; it != same_end; ++it) {
                same.emplace_back(&label(it->index), it->index);
            }
        } else {
            // no median: the whole group competes by label, selected where it is (|order|)
            c.tied_same_median = tied;
            std::nth_element(order.begin() + above, order.begin() + above + slots - 1,
                             order.end(), [&](size_t a, size_t b) {
                                 if (const int labels = label(a).compare(label(b)))
                                     return labels < 0;
                                 return a < b;
                             });
            for (size_t p = above; p < above + slots; ++p) {
                rows.push_back(Row{ order[p], c.nth_score, 0, label(order[p]) });
            }
        }
        if (c.nth_median_read) {
            c.tied_same_median = same.size();
            tie = std::vector<Tied>();
            // the first rows of the N-th row's median in the C5 order take the places left
            const size_t left = slots - (rows.size() - above);
            assert(left && left <= same.size());
            std::nth_element(same.begin(), same.begin() + left - 1, same.end(),
                             [](const auto &a, const auto &b) {
                                 if (const int labels = a.first->compare(*b.first))
                                     return labels < 0;
                                 return a.second < b.second;
                             });
            for (size_t r = 0; r < left; ++r) {
                rows.push_back(Row{ same[r].second, c.nth_score, c.nth_median, *same[r].first });
            }
        }
    }
    if (!c.abundances) {
        c.nth_median = c.tied_same_median = 0;
        c.nth_median_read = false;
    }

    for (size_t r = 0; r < above; ++r) {
        rows[r].label = label(rows[r].index);
    }
    std::sort(rows.begin(), rows.end(), before);

    c.returned = rows.size();
    if (cut)
        *cut = c;

    std::vector<std::pair<size_t, std::string>> result;
    result.reserve(rows.size());
    for (const Row &row : rows) {
        result.emplace_back(row.index, std::string(row.label));
    }
    return result;
}

namespace {

// the largest position a run of 32-bit positions takes (its end): a query of more k-mers is
// held densely from the start
constexpr size_t kMaxRunsSize = std::numeric_limits<uint32_t>::max();

// Room for one more element in |v| within |budget| bytes of buffer: its capacity doubled, at
// most what |budget| holds; false if |budget| holds no more elements than |v| has
template <typename T>
bool grow_within(std::vector<T> &v, size_t budget) {
    const size_t max_capacity = budget / sizeof(T);
    if (v.size() >= max_capacity)
        return false;
    v.reserve(std::min(max_capacity, std::max<size_t>(2, 2 * v.capacity())));
    return true;
}

} // namespace

PresenceRuns::PresenceRuns(size_t size) : size(size) {
    if (size > kMaxRunsSize)
        to_dense_form();
}

PresenceRuns::PresenceRuns(const sdsl::bit_vector &mask) : PresenceRuns(mask.size()) {
    for (size_t i = 0; i < mask.size(); ++i) {
        if (mask[i])
            add(i);
    }
}

bool PresenceRuns::grow() {
    return grow_within(runs_, dense_bytes(size));
}

void PresenceRuns::to_dense_form() {
    dense_ = sdsl::bit_vector(size, false);
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
    sdsl::bit_vector mask(size, false);
    for (const Run &run : runs_) {
        for (size_t i = run.begin; i < run.end; ++i) {
            mask[i] = true;
        }
    }
    return mask;
}

size_t PresenceRuns::heap_bytes() const {
    return dense_form_ ? dense_bytes(size) : runs_.capacity() * sizeof(Run);
}

bool PresenceRuns::operator==(const PresenceRuns &other) const {
    return size == other.size && to_dense() == other.to_dense();
}

CountRuns::CountRuns(size_t size) : size(size) {
    if (size > kMaxRunsSize)
        to_dense_form();
}

CountRuns::CountRuns(const std::vector<size_t> &values) : CountRuns(values.size()) {
    for (size_t i = 0; i < values.size(); ++i) {
        add(i, values[i]);
    }
}

bool CountRuns::grow() {
    return grow_within(runs_, dense_bytes(size));
}

void CountRuns::to_dense_form() {
    dense_.assign(size, 0);
    dense_form_ = true;
    for (const Run &run : runs_) {
        std::fill(dense_.begin() + run.begin, dense_.begin() + run.end, run.value);
    }
    std::vector<Run>().swap(runs_);
}

size_t CountRuns::count() const {
    if (dense_form_)
        return size - std::count(dense_.begin(), dense_.end(), size_t(0));
    size_t count = 0;
    for (const Run &run : runs_) {
        count += run.end - run.begin;
    }
    return count;
}

std::vector<size_t> CountRuns::to_dense() const {
    if (dense_form_)
        return dense_;
    std::vector<size_t> values(size, 0);
    for (const Run &run : runs_) {
        std::fill(values.begin() + run.begin, values.begin() + run.end, run.value);
    }
    return values;
}

size_t CountRuns::heap_bytes() const {
    return dense_form_ ? dense_.capacity() * sizeof(size_t) : runs_.capacity() * sizeof(Run);
}

bool CountRuns::operator==(const CountRuns &other) const {
    return size == other.size && to_dense() == other.to_dense();
}

PositionCoords::PositionCoords(size_t size) : size(size) {
    if (size > kMaxRunsSize)
        to_dense_form();
}

PositionCoords::PositionCoords(const std::vector<Coords> &values) : PositionCoords(values.size()) {
    for (size_t i = 0; i < values.size(); ++i) {
        add(i, Coords(values[i]));
    }
}

bool PositionCoords::open(size_t position) {
    const bool new_run = runs_.empty() || runs_.back().end != position;
    // the runs and the coordinates within the bytes of a coordinate list per position: each
    // grows within what the other's capacity leaves
    const size_t budget = dense_bytes(size);
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
    dense_.resize(size);
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
    std::vector<Coords> values(size);
    call_positions([&](size_t i, const Coords &coords) { values[i] = coords; });
    return values;
}

size_t PositionCoords::heap_bytes() const {
    return dense_form_ ? dense_.capacity() * sizeof(Coords)
                       : runs_.capacity() * sizeof(Run) + coords_.capacity() * sizeof(Coords);
}

bool PositionCoords::operator==(const PositionCoords &other) const {
    return size == other.size && to_dense() == other.to_dense();
}

using ColumnMedians = std::function<std::vector<uint64_t>(const std::vector<Column> &)>;

void LowerMedianCounter::add(uint64_t value, uint64_t positions) {
    if (!value || !positions)
        return;
    if (counts_.size() && counts_.back().first == value) {
        // a run of one value along the query
        counts_.back().second += positions;
        return;
    }
    counts_.emplace_back(value, positions);
    if (counts_.size() >= 2 * merged_ + 16)
        merge();
}

void LowerMedianCounter::merge() {
    std::sort(counts_.begin(), counts_.end(), utils::LessFirst());
    size_t last = 0;
    for (size_t i = 1; i < counts_.size(); ++i) {
        if (counts_[i].first == counts_[last].first) {
            counts_[last].second += counts_[i].second;
        } else {
            counts_[++last] = counts_[i];
        }
    }
    if (counts_.size())
        counts_.resize(last + 1);
    merged_ = counts_.size();
}

uint64_t LowerMedianCounter::lower_median() {
    merge();
    uint64_t n = 0;
    for (const auto &[value, count] : counts_) {
        n += count;
    }
    if (!n)
        return 0;
    // the element at index floor((n - 1) / 2) of the values sorted ascending
    uint64_t rank = (n - 1) / 2;
    for (const auto &[value, count] : counts_) {
        if (rank < count)
            return value;
        rank -= count;
    }
    assert(false);
    return 0;
}

// The lower medians (lower_median) of the abundances of |columns| over the query's k-mer
// positions: |call_rows| calls its callback with the row of (column, abundance) pairs of each
// k-mer of the query found in the graph and the number of query positions that row stands for
// (a k-mer repeated in the query is a position each time). Each column asked holds its
// distinct values with their counts (LowerMedianCounter), not a value per position.
template <class CallRows>
std::vector<uint64_t> abundance_medians(const CallRows &call_rows,
                                        const std::vector<Column> &columns,
                                        size_t num_columns) {
    utils::ValueStore<Column, LowerMedianCounter> values(num_columns, columns.size());
    for (Column j : columns) {
        values.emplace(j, LowerMedianCounter());
    }
    call_rows([&](const auto &row, size_t positions) {
        for (const auto &[j, value] : row) {
            if (LowerMedianCounter *counter = values.get(j))
                counter->add(value, positions);
        }
    });
    std::vector<uint64_t> medians;
    medians.reserve(columns.size());
    for (Column j : columns) {
        medians.push_back(values.get(j)->lower_median());
    }
    return medians;
}

std::vector<uint64_t>
abundance_medians(const IntMatrix &matrix,
                  const std::vector<std::pair<uint64_t, size_t>> &row_counts,
                  const std::vector<uint32_t> &columns,
                  size_t num_columns,
                  size_t rows_per_read) {
    assert(rows_per_read);
    // per-row access is cheap for a CSRMatrix (a batched query's slice); other matrices
    // decode their rows in batches of |rows_per_read|, so that the values decoded at a time
    // stay bounded whatever the query's length
    if (const auto *csr = dynamic_cast<const CSRMatrix *>(&matrix)) {
        return abundance_medians([&](const auto &callback) {
            for (const auto &[row, positions] : row_counts) {
                callback(csr->get_row_values(row), positions);
            }
        }, columns, num_columns);
    }
    // in the order of the rows (the medians are of the values, whatever their order): the
    // matrix is read as a batched query's slice reads it, a row-diff matrix in the query's
    // order (along its paths), every other sorted (its rank structures read in order)
    const std::vector<std::pair<uint64_t, size_t>> *rows_in_order = &row_counts;
    std::vector<std::pair<uint64_t, size_t>> sorted;
    if (!dynamic_cast<const annot::matrix::IRowDiff *>(&matrix)) {
        sorted = row_counts;
        std::sort(sorted.begin(), sorted.end(), utils::LessFirst());
        rows_in_order = &sorted;
    }
    return abundance_medians([&](const auto &callback) {
        const auto &counts = *rows_in_order;
        std::vector<AnnotatedDBG::row_index> rows;
        for (size_t begin = 0; begin < counts.size(); begin += rows_per_read) {
            const size_t end = std::min(begin + rows_per_read, counts.size());
            rows.clear();
            for (size_t i = begin; i < end; ++i) {
                rows.push_back(counts[i].first);
            }
            auto row_values = matrix.get_row_values(rows);
            for (size_t i = begin; i < end; ++i) {
                callback(row_values[i - begin], counts[i].second);
            }
        }
    }, columns, num_columns);
}

// The rows of |row_indices| (one per query position), each with the number of positions it
// stands for, in the order of their first position
std::vector<std::pair<AnnotatedDBG::row_index, size_t>>
count_rows(const std::vector<AnnotatedDBG::row_index> &row_indices) {
    VectorMap<AnnotatedDBG::row_index, size_t> row_counts;
    row_counts.reserve(row_indices.size());
    for (AnnotatedDBG::row_index row : row_indices) {
        row_counts[row]++;
    }
    return to_vector(std::move(row_counts));
}

// The rows read at a time where the medians of a search's ties are read from a matrix other
// than a CSRMatrix: integer values take four times the bytes of the set bits a binary row
// holds, so a quarter of the rows BinaryMatrix::sum_rows decodes at a time
constexpr size_t kMedianRowsPerRead = 2'500;

ColumnMedians
AnnotatedDBG::tie_medians(std::string_view sequence,
                          const std::vector<std::pair<row_index, size_t>> &row_counts) const {
    if (const IntMatrix *abundances = abundance_matrix()) {
        // the matrix is captured by value, its pointer ends with the condition's scope
        return [this, &row_counts, abundances](const std::vector<Column> &columns) {
            return abundance_medians(*abundances, row_counts, columns,
                                     annotator_->num_labels(), kMedianRowsPerRead);
        };
    }
    return source_tie_medians(sequence);
}

ColumnMedians
AnnotatedDBG::tie_medians(std::string_view sequence,
                          const std::vector<row_index> &row_indices) const {
    if (const IntMatrix *abundances = abundance_matrix()) {
        // the rows collapsed with their multiplicity (a k-mer repeated in the query is read
        // once), only when a tie asks
        return [this, &row_indices, abundances](const std::vector<Column> &columns) {
            return abundance_medians(*abundances, count_rows(row_indices), columns,
                                     annotator_->num_labels(), kMedianRowsPerRead);
        };
    }
    return source_tie_medians(sequence);
}

ColumnMedians AnnotatedDBG::source_tie_medians(std::string_view sequence) const {
    if (!kmer_abundances_ || !abundance_source_ || !abundance_source_->abundance_matrix())
        return {};
    // a binary slice of a counts annotation: the rows of the full index the query maps to,
    // each with its multiplicity in the query (the rows of this index's graph stand for the
    // same k-mers; a k-mer the full graph lacks has an empty row here)
    return [this, sequence](const std::vector<Column> &columns) {
        const AnnotatedDBG &source = *abundance_source_;
        VectorMap<row_index, size_t> source_counts;
        source.get_graph().map_to_nodes(sequence, [&](node_index i) {
            if (i > 0)
                source_counts[graph_to_anno_index(i)]++;
        });
        return abundance_medians(*source.abundance_matrix(), source_counts.values_container(),
                                 columns, source.get_annotator().num_labels(),
                                 kMedianRowsPerRead);
    };
}

// The candidates' indexes mapped to the columns |medians| reads: what search_top_n asks of
// a list of (column, score) pairs
template <class Container>
std::function<std::vector<uint64_t>(const std::vector<size_t> &)>
medians_of_candidates(const Container &code_counts, const ColumnMedians &medians) {
    if (!medians)
        return {};
    return [&code_counts, medians](const std::vector<size_t> &indexes) {
        std::vector<Column> columns;
        columns.reserve(indexes.size());
        for (size_t i : indexes) {
            columns.push_back(code_counts[i].first);
        }
        return medians(columns);
    };
}

// The C5 top N of (column, count) pairs (search_top_n), the labels decoded
template <class Container>
std::vector<StringCountPair> filter_and_decode(Container&& code_counts,
                                               const annot::LabelEncoder<> &label_encoder,
                                               size_t num_top_labels,
                                               const ColumnMedians &medians,
                                               SearchCut *cut) {
    std::vector<size_t> scores;
    scores.reserve(code_counts.size());
    for (const auto &[j, count] : code_counts) {
        scores.push_back(count);
    }
    auto top = search_top_n(scores, num_top_labels,
                            [&](size_t i) -> const std::string & {
                                return label_encoder.decode(code_counts[i].first);
                            },
                            medians_of_candidates(code_counts, medians), cut);

    std::vector<StringCountPair> label_counts;
    label_counts.reserve(top.size());
    for (auto &[i, label] : top) {
        label_counts.emplace_back(std::move(label), code_counts[i].second);
    }

    return label_counts;
}

std::vector<StringCountPair>
AnnotatedDBG::get_top_labels(std::string_view sequence,
                             size_t num_top_labels,
                             double discovery_fraction,
                             double presence_fraction,
                             bool with_kmer_counts,
                             SearchCut *cut) const {
    assert(discovery_fraction >= 0.);
    assert(discovery_fraction <= 1.);
    assert(presence_fraction >= 0.);
    assert(presence_fraction <= 1.);
    assert(check_compatibility());

    // what a search that returns no row states (complete), replaced by its cut otherwise
    if (cut) {
        *cut = SearchCut();
        cut->abundances = kmer_abundances_;
    }

    if (sequence.size() < dbg_.get_k())
        return {};

    if (coord_to_header_) {
        // the records ranked by the score reported: with |with_kmer_counts| the number of
        // their coordinates matched
        auto counts = get_header_results<size_t>(sequence, num_top_labels,
                                                 discovery_fraction, presence_fraction,
                                                 with_kmer_counts, cut);
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

    // the abundances of the ties, read from the rows matched (each with its multiplicity in
    // the query) only where a tie needs them
    ColumnMedians medians = tie_medians(sequence, index_counts.values_container());

    std::vector<StringCountPair> result;
    if (with_kmer_counts) {
        result = filter_and_decode(dynamic_cast<const IntMatrix &>(annotator_->get_matrix())
                                        .sum_row_values(index_counts.values_container(), min_count),
                                   annotator_->get_label_encoder(), num_top_labels, medians, cut);
    } else {
        result = filter_and_decode(annotator_->get_matrix()
                                        .sum_rows(index_counts.values_container(), min_count),
                                   annotator_->get_label_encoder(), num_top_labels, medians, cut);
    }

    assert(with_kmer_counts || std::all_of(result.begin(), result.end(),
                                    [&](const auto &pair) { return pair.second <= num_kmers; }));

    return result;
}


namespace {
using Header = std::pair<Column, size_t>; // (column, seq_id)

// a query position without a row (its k-mer not in the graph)
constexpr size_t kNoRow = std::numeric_limits<size_t>::max();

// A record of a CoordToHeader mapping the query's k-mers match: the query positions they
// cover, and the number of the record's coordinates they match
struct HeaderMatch {
    Header header;
    size_t kmers;
    uint64_t coordinates;
};

/**
 * The records the rows |rows| match (their coordinates packed by map_to_local_coords), each
 * matched at |min_count| query positions or more, the top |num_top_labels| of them in the C5
 * order (search_top_n) of the score reported: the number of their coordinates matched if
 * |by_coordinates|, of the query positions otherwise; ties by the header's bytes (coordinates
 * are positions, not abundances), a header in two columns in (column, seq_id) order. |weights|:
 * the number of query positions each row stands for (one each if empty), so that a row the
 * query repeats is decoded once. |cut| states the cut. Returns the records with their headers.
 */
std::vector<std::pair<HeaderMatch, std::string>>
top_header_matches(const std::vector<MultiIntMatrix::RowTuples> &rows,
                   const annot::CoordToHeader &mapping,
                   size_t min_count, size_t num_top_labels, bool by_coordinates,
                   std::span<const size_t> weights, SearchCut *cut) {
    assert(weights.empty() || weights.size() == rows.size());
    std::vector<HeaderMatch> candidates;
    {
        VectorMap<Header, std::pair<size_t, uint64_t>> counts_map;
        tsl::hopscotch_set<Header> matches;
        for (size_t r = 0; r < rows.size(); ++r) {
            const size_t weight = weights.empty() ? 1 : weights[r];
            matches.clear();
            for (const auto &[column, coords] : rows[r]) {
                for (uint64_t coord : coords) {
                    Header h(column, coord % mapping.num_sequences(column));
                    matches.emplace(h);
                    if (by_coordinates)
                        counts_map[h].second += weight;
                }
            }
            for (const Header &h : matches) {
                counts_map[h].first += weight;
            }
        }
        for (auto &[h, c] : to_vector(std::move(counts_map))) {
            // the records with too few k-mer matches are no candidates
            if (c.first >= min_count)
                candidates.push_back({ h, c.first, c.second });
        }
    }
    // the records in (column, seq_id) order, so that a header in two columns is returned in
    // the same order by every request
    std::sort(candidates.begin(), candidates.end(),
              [](const HeaderMatch &a, const HeaderMatch &b) { return a.header < b.header; });

    std::vector<size_t> scores;
    scores.reserve(candidates.size());
    for (const HeaderMatch &c : candidates) {
        scores.push_back(by_coordinates ? c.coordinates : c.kmers);
    }
    auto top = search_top_n(scores, num_top_labels, [&](size_t i) -> const std::string & {
        const auto &[column, seq_id] = candidates[i].header;
        return mapping.get_headers(column)[seq_id];
    }, {}, cut);

    std::vector<std::pair<HeaderMatch, std::string>> result;
    result.reserve(top.size());
    for (auto &[i, label] : top) {
        result.emplace_back(candidates[i], std::move(label));
    }
    return result;
}
} // namespace


template <class Values>
std::vector<std::tuple<Label, size_t, Values>>
AnnotatedDBG::get_header_results(std::string_view sequence, size_t num_top_labels,
                                 double discovery_fraction, double presence_fraction,
                                 bool with_kmer_counts, SearchCut *cut) const {
    static_assert(utils::is_one_of_v<Values, size_t, CountRuns, PresenceRuns>);
    assert((!with_kmer_counts || std::is_same_v<Values, size_t>));
    assert(coord_to_header_);
    const size_t num_kmers = sequence.size() - dbg_.get_k() + 1;
    // Repeated query positions still count separately, but their annotation row is decoded
    // once. Only the results per position keep the positions.
    using Positions = std::conditional_t<std::is_same_v<Values, size_t>,
                                         size_t, std::vector<size_t>>;
    std::vector<std::pair<row_index, Positions>> index_counts;
    size_t num_present = 0;
    {
        VectorMap<row_index, Positions> counts;
        counts.reserve(num_kmers);
        size_t position = 0;
        dbg_.map_to_nodes(sequence, [&](node_index node) {
            if (node) {
                auto &positions = counts[graph_to_anno_index(node)];
                if constexpr (std::is_same_v<Values, size_t>) {
                    ++positions;
                } else {
                    positions.push_back(position);
                }
                ++num_present;
            }
            ++position;
        });
        index_counts = to_vector(std::move(counts));
    } // the hash table is released before the rows are decoded
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
        if constexpr (std::is_same_v<Values, size_t>) {
            weights.push_back(positions);
        } else {
            weights.push_back(positions.size());
        }
    }
    // The values per position below visit the query's positions in ascending order: the row
    // (its index in |rows|) of each position, |kNoRow| where none; the positions per row are
    // released
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
    // the constructor refused a mapping of another column count: a request never exits here
    if (annotator_->num_labels() != coord_to_header_->num_columns()) {
        throw std::logic_error(fmt::format(
                "Incompatible number of columns in CoordToHeader and annotation matrix: {} != {}",
                coord_to_header_->num_columns(), annotator_->num_labels()));
    }
    auto tuples = matrix->get_row_tuples(rows);
    coord_to_header_->map_to_local_coords(&tuples);
    // Now each coord in `tuples` is `coord = local_coord * num_sequences[j] + seq_id`

    auto top = top_header_matches(tuples, *coord_to_header_, min_count, num_top_labels,
                                  with_kmer_counts, weights, cut);

    std::vector<std::tuple<Label, size_t, Values>> result;
    result.reserve(top.size());
    if constexpr (std::is_same_v<Values, size_t>) {
        for (auto &[match, label] : top) {
            result.emplace_back(std::move(label), match.kmers,
                                with_kmer_counts ? match.coordinates : match.kmers);
        }
    } else {
        // Each record's values are held only at the positions it matches while that takes
        // fewer bytes than a value of num_kmers (which made the memory the records returned
        // times the query's k-mers), and as that value otherwise: the query's positions are
        // visited in ascending order, so that each record's runs are appended in order, a
        // block of consecutive positions of one row at once. A row's records and their
        // occurrences in it are found at its first position and kept while its other
        // positions are due; its decoded tuples are released then.
        tsl::hopscotch_map<Header, size_t> returned; // a record's index in |top|
        returned.reserve(top.size());
        for (size_t r = 0; r < top.size(); ++r) {
            returned.emplace(top[r].first.header, r);
        }
        std::vector<Values> values(top.size(), Values(num_kmers));
        // (record index, its occurrences in the row) per row found, while the row is due
        std::vector<std::vector<std::pair<size_t, size_t>>> hits(tuples.size());
        std::vector<bool> found(tuples.size(), false);
        // the positions of each row not visited yet
        std::vector<size_t> &due = weights;
        tsl::hopscotch_map<Header, size_t> occurrences;
        for (size_t position = 0; position < num_kmers; ) {
            const size_t i = row_at[position];
            if (i == kNoRow) {
                ++position;
                continue;
            }
            // the positions [position, end) of row |i| (a k-mer repeated along the query)
            size_t end = position + 1;
            while (end < num_kmers && row_at[end] == i) {
                ++end;
            }
            if (!found[i]) {
                found[i] = true;
                occurrences.clear();
                for (const auto &[column, coords] : tuples[i]) {
                    for (uint64_t coord : coords) {
                        ++occurrences[{ column, coord % coord_to_header_->num_sequences(column) }];
                    }
                }
                for (const auto &[h, count] : occurrences) {
                    auto it = returned.find(h);
                    if (it != returned.end())
                        hits[i].emplace_back(it->second, count);
                }
                tuples[i] = MultiIntMatrix::RowTuples();
            }
            for (const auto &[r, count] : hits[i]) {
                if constexpr (std::is_same_v<Values, PresenceRuns>) {
                    values[r].add_run(position, end);
                } else {
                    values[r].add_run(position, end, count);
                }
            }
            due[i] -= end - position;
            if (!due[i])
                hits[i] = std::vector<std::pair<size_t, size_t>>();
            position = end;
        }
        for (size_t r = 0; r < top.size(); ++r) {
            values[r].shrink();
            result.emplace_back(std::move(top[r].second), top[r].first.kmers,
                                std::move(values[r]));
        }
    }
    return result;
}


// The C5 top N of the labels of the rows |enumerate_rows| calls (the query's k-mers found in
// the graph, at |kmer_positions|) with at least |min_count| of them, each with its value per
// query position; |medians| and |cut| as filter_and_decode's
template <class Result, class RowEnumerator>
Result filter_and_aggregate(const RowEnumerator &enumerate_rows,
                            const annot::LabelEncoder<Label> &label_encoder,
                            size_t min_count,
                            size_t num_top_labels,
                            const std::vector<size_t> &kmer_positions,
                            size_t num_kmers,
                            const ColumnMedians &medians,
                            SearchCut *cut) {
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
    std::vector<std::pair<Column, size_t>> candidates
            = utils::accumulate_counts(call_bits, label_encoder.size(), min_count);

    std::vector<size_t> scores;
    scores.reserve(candidates.size());
    for (const auto &[j, count] : candidates) {
        scores.push_back(count);
    }
    auto top = search_top_n(scores, num_top_labels,
                            [&](size_t i) -> const std::string & {
                                return label_encoder.decode(candidates[i].first);
                            },
                            medians_of_candidates(candidates, medians), cut);
    if (top.empty())
        return {};

    // the rows returned, in the C5 order
    std::vector<std::pair<Column, size_t>> counts;
    counts.reserve(top.size());
    for (const auto &[i, label] : top) {
        counts.push_back(candidates[i]);
    }
    candidates = decltype(candidates)();

    // each label's values held at the positions it matches (or as a value per position where
    // that takes fewer bytes); the rows are enumerated in the order of their positions, which
    // the values are built in
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
    for (size_t r = 0; r < counts.size(); ++r) {
        const auto &[j, count] = counts[r];
        ValueType *value = value_store.get(j);
        value->shrink();
        result.emplace_back(std::move(top[r].second), count, std::move(*value));
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

#ifndef NDEBUG
// A sanity check to ensure the same matches are for get_top_label_signatures and get_top_labels
template <typename ValueType>
bool same_results(const std::vector<std::tuple<Label, size_t, ValueType>> &result,
                  const std::vector<StringCountPair> &top_labels) {
    assert(top_labels.size() == result.size());
    std::unordered_map<Label, uint64_t> check(top_labels.begin(), top_labels.end());
    for (const auto &[label, count, value] : result) {
        if constexpr(utils::is_one_of_v<ValueType, PresenceRuns, CountRuns, PositionCoords>) {
            // the positions marked, with a nonzero count or with coordinates
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
                              double presence_fraction,
                              SearchCut *cut) const {
    assert(discovery_fraction >= 0.);
    assert(discovery_fraction <= 1.);
    assert(presence_fraction >= 0.);
    assert(presence_fraction <= 1.);
    assert(check_compatibility());

    // what a search that returns no row states (complete), replaced by its cut otherwise
    if (cut) {
        *cut = SearchCut();
        cut->abundances = kmer_abundances_;
    }

    if (sequence.size() < dbg_.get_k())
        return {};

    if (coord_to_header_) {
        return get_header_results<CountRuns>(sequence, num_top_labels,
                                             discovery_fraction, presence_fraction,
                                             false, cut);
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
        // the abundances of the ties, from the rows read for the counts (one per position)
        ColumnMedians medians;
        if (kmer_abundances_) {
            medians = [&](const std::vector<Column> &columns) {
                return abundance_medians([&](const auto &callback) {
                    enumerate_rows([&](const auto &row) { callback(row, 1); });
                }, columns, annotator_->num_labels());
            };
        }
        auto results = filter_and_aggregate<std::vector<std::tuple<std::string, size_t, CountRuns>>>(
            enumerate_rows, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers,
            medians, cut);
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
                                   double presence_fraction,
                                   SearchCut *cut) const {
    assert(discovery_fraction >= 0.);
    assert(discovery_fraction <= 1.);
    assert(presence_fraction >= 0.);
    assert(presence_fraction <= 1.);
    assert(check_compatibility());

    // what a search that returns no row states (complete), replaced by its cut otherwise; a
    // coordinate annotation holds no abundances
    if (cut) {
        *cut = SearchCut();
        cut->abundances = kmer_abundances_;
    }

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
        // the constructor refused a mapping of another column count: a request never exits
        if (annotator_->num_labels() != coord_to_header_->num_columns()) {
            throw std::logic_error(fmt::format(
                    "Incompatible number of columns in CoordToHeader and annotation matrix: "
                    "{} != {}", coord_to_header_->num_columns(), annotator_->num_labels()));
        }

        coord_to_header_->map_to_local_coords(&rows_tuples);
        // Now each coord in `rows_tuples` is `coord = local_coord * num_sequences[j] + seq_id`

        // the records ranked by the number of query positions they match
        auto top = top_header_matches(rows_tuples, *coord_to_header_, min_count,
                                      num_top_labels, false, {}, cut);
        tsl::hopscotch_map<Header, size_t> returned; // a record's index in |top|
        returned.reserve(top.size());
        for (size_t r = 0; r < top.size(); ++r) {
            returned.emplace(top[r].first.header, r);
        }

        // split coordinates across sequences, each record's held at the positions it matches
        // (or as a list per position where that takes fewer bytes): the rows are in the order
        // of the query's positions, which they are added in
        std::vector<PositionCoords> coords_of(top.size(), PositionCoords(num_kmers));
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

        std::vector<std::tuple<std::string, size_t, PositionCoords, size_t>> result;
        result.reserve(top.size());
        for (size_t r = 0; r < top.size(); ++r) {
            const auto &[col, header] = top[r].first.header;
            coords_of[r].shrink();
            result.emplace_back(std::move(top[r].second),
                                top[r].first.kmers, std::move(coords_of[r]),
                                coord_to_header_->num_kmers_in_sequence(col, header));
        }

        return result;
    }

    auto call_rows_tuples = [&](const auto &callback) {
        std::for_each(rows_tuples.begin(), rows_tuples.end(), callback);
    };

    // ties ordered by label: the coordinates are positions, not abundances
    auto results = filter_and_aggregate<std::vector<std::tuple<Label, size_t, PositionCoords>>>(
        call_rows_tuples, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers,
        {}, cut);
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
                                       SearchCut *cut) const {
    assert(discovery_fraction >= 0.);
    assert(discovery_fraction <= 1.);
    assert(presence_fraction >= 0.);
    assert(presence_fraction <= 1.);
    assert(check_compatibility());

    // what a search that returns no row states (complete), replaced by its cut otherwise
    if (cut) {
        *cut = SearchCut();
        cut->abundances = kmer_abundances_;
    }

    if (sequence.size() < dbg_.get_k())
        return {};

    size_t num_kmers = sequence.size() - dbg_.get_k() + 1;

    if (discovery_fraction == 1.) {
        auto label_counts = get_top_labels(sequence, num_top_labels, discovery_fraction,
                                           presence_fraction, false, cut);
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

    if (coord_to_header_) {
        return get_header_results<PresenceRuns>(sequence, num_top_labels,
                                                discovery_fraction, presence_fraction,
                                                false, cut);
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
    // the abundances of the ties: the rows enumerated below may be binary, so they are read
    // from the integer matrix (one value per position), only where a tie needs them
    ColumnMedians medians = tie_medians(sequence, row_indices);
    auto get_results = [&](const auto &enumerate_rows) {
        auto results = filter_and_aggregate<std::vector<std::tuple<std::string, size_t, PresenceRuns>>>(
            enumerate_rows, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers,
            medians, cut);
        assert(same_results(results, get_top_labels(sequence, num_top_labels, discovery_fraction, presence_fraction)));
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
    return compatible(*graph_, *annotator_);
}

bool AnnotatedSequenceGraph::compatible(const SequenceGraph &graph,
                                        const Annotator &annotation) {
    // TODO: what if CanonicalDBG is not the highest level? find a better way to do this
    if (const auto *canonical = dynamic_cast<const CanonicalDBG *>(&graph)) {
        return canonical->get_graph().max_index() == annotation.num_objects();
    }
    return graph.max_index() == annotation.num_objects();
}

} // namespace graph
} // namespace mtg
