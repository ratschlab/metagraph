#include "annotated_dbg.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
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
#include "annotation/binary_matrix/base/decode_budget.hpp"
#include "annotation/binary_matrix/base/hopscotch_layout.hpp"
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
using mtg::annot::matrix::alloc_bytes;
using mtg::annot::matrix::buffer_bytes;

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

const char* refused_account(SearchCut::Refused refused) {
    switch (refused) {
        case SearchCut::Refused::VALUES:
            return "values";
        case SearchCut::Refused::ANSWER:
            return "answer";
        case SearchCut::Refused::NO:
            break;
    }
    return nullptr;
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
        for (node_index i : source.get_graph().node_range(sequence)) {
            if (i > 0)
                source_counts[graph_to_anno_index(i)]++;
        }
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
using HeaderCounts = std::vector<std::pair<Header, size_t>>;

// Annotation rows decoded at a time for the queries by sequence record. A row holds a
// coordinate for every record with its k-mer, so the rows of a whole query can take far
// more memory than the result. Row-diff rows share their paths within a block, which a
// small block loses.
constexpr size_t kRowBlock = 256;

// Bytes of per-position values that the fold of a query by sequence record holds before
// it frees those of the records that can no longer be returned. If the records that
// still can be returned hold more than half of it and are more than the query returns,
// the fold keeps only the counts, and the rows are read again after the cut for the
// returned records. Debug builds have no budget, so that the tests check both ways.
// (Staging: one fold per query record on one graph, a task of the unit's values account;
// what the fold holds before the cut is not charged to that account, so that the cut is
// made before the first charge. Lane T7-M5.)
#ifdef NDEBUG
constexpr size_t kFoldValuesBudget = size_t(16) << 20;
#else
constexpr size_t kFoldValuesBudget = 0;
#endif
// what the folds use (tests set it: set_record_fold_values_budget, review M2, lane T7-M5),
// and what they did (record_fold_counts)
std::atomic<size_t> fold_values_budget = kFoldValuesBudget;
std::atomic<uint64_t> folds_freed = 0;
std::atomic<uint64_t> folds_met_late = 0;
std::atomic<uint64_t> folds_counts_only = 0;

/**
 * The buffers that hold the values of the returned records after the cut, until they are
 * those records' values: the runs of the rows read again for them, the values at the
 * repeated rows, and a record's values before they are sorted. Each is charged to the
 * unit's values account (null: none, the CLI's) before it grows, as the scratch of the
 * task (SearchAccount::charge_scratch: the most one task holds at once, times the tasks that
 * run at once), and released when it is freed. The values themselves charge as they grow
 * (PresenceRuns, CountRuns, PositionCoords). What the fold holds before the cut is not
 * charged here (its budget: the fold's budget_, set_record_fold_values_budget, by default
 * kFoldValuesBudget). Lane T7-M5: train 7 P1's R8 (a).
 */
class CutScratch {
  public:
    explicit CutScratch(SearchAccount *account) : account_(account) {}

    // before |bytes| more are allocated
    void grow(uint64_t bytes) {
        charge_search_scratch(account_, held_ + bytes);
        held_ += bytes;
    }
    void release(uint64_t bytes) { held_ -= std::min(held_, bytes); }

    // room for one more element in |v|, its capacity doubled (charged first)
    template <typename T>
    void make_room(std::vector<T> &v) {
        if (v.size() < v.capacity())
            return;
        reserve(v, std::max<size_t>(1, 2 * v.capacity()));
    }
    // |v| reserved for |n| elements (the new buffer less the old one charged first)
    template <typename T>
    void reserve(std::vector<T> &v, size_t n) {
        if (n <= v.capacity())
            return;
        grow(buffer_bytes(n, sizeof(T)) - buffer_bytes(v.capacity(), sizeof(T)));
        v.reserve(n);
    }
    // |v| freed
    template <typename T>
    void free(std::vector<T> &v) {
        release(buffer_bytes(v.capacity(), sizeof(T)));
        std::vector<T>().swap(v);
    }

  private:
    SearchAccount *account_;
    uint64_t held_ = 0;
};

// Append the rows of |rows| not seen yet to |distinct|, in the order of their first
// occurrence, and the index in |distinct| of each row to |index|
void index_distinct_rows(std::span<const AnnotatedDBG::row_index> rows,
                         tsl::hopscotch_map<AnnotatedDBG::row_index, uint32_t> &seen,
                         std::vector<AnnotatedDBG::row_index> &distinct,
                         std::vector<uint32_t> &index) {
    for (AnnotatedDBG::row_index row : rows) {
        auto [it, inserted] = seen.try_emplace(row, distinct.size());
        if (inserted)
            distinct.push_back(row);
        index.push_back(it->second);
    }
}

// Call |callback|(i, row of rows[i]) for i in order, decoding |kRowBlock| rows at a time
// and each distinct row of a block once
template <class Callback>
void call_row_tuples(const MultiIntMatrix &matrix, std::span<const AnnotatedDBG::row_index> rows,
                     const Callback &callback) {
    std::vector<AnnotatedDBG::row_index> distinct;
    std::vector<uint32_t> index;
    tsl::hopscotch_map<AnnotatedDBG::row_index, uint32_t> seen;
    for (size_t begin = 0; begin < rows.size(); begin += kRowBlock) {
        const size_t end = std::min(rows.size(), begin + kRowBlock);
        distinct.clear();
        index.clear();
        seen.clear();
        index_distinct_rows(rows.subspan(begin, end - begin), seen, distinct, index);
        const auto decoded = matrix.get_row_tuples(distinct);
        for (size_t i = begin; i < end; ++i) {
            callback(i, decoded[index[i - begin]]);
        }
    }
}

// The rows of a query, each distinct row once. A k-mer repeated m times (a poly-A
// query) has one row, which can hold millions of coordinates.
struct QueryRows {
    // the rows in the order of their first position
    std::vector<AnnotatedDBG::row_index> distinct;
    // the index in |distinct| of the row at each position
    std::vector<uint32_t> row_of;
    // of each distinct row: its first query position and its number of positions
    std::vector<size_t> first_position;
    std::vector<uint32_t> weight;
    // the distinct rows of more than one position, and the index in |repeated| of each
    std::vector<AnnotatedDBG::row_index> repeated;
    std::vector<uint32_t> repeated_of;

    QueryRows(std::span<const AnnotatedDBG::row_index> rows, std::span<const size_t> positions) {
        {
            tsl::hopscotch_map<AnnotatedDBG::row_index, uint32_t> seen;
            index_distinct_rows(rows, seen, distinct, row_of);
        }
        first_position.resize(distinct.size());
        weight.assign(distinct.size(), 0);
        for (size_t i = 0; i < row_of.size(); ++i) {
            if (!weight[row_of[i]]++)
                first_position[row_of[i]] = positions[i];
        }
        repeated_of.resize(distinct.size());
        for (size_t i = 0; i < distinct.size(); ++i) {
            if (weight[i] > 1) {
                repeated_of[i] = repeated.size();
                repeated.push_back(distinct[i]);
            }
        }
    }
};

// Call |callback|(j, r, coordinate) for each coordinate of rows[j] in the record
// counts[r], for j in order. A column without sequence records is skipped, as master's
// (#764); a column outside the mapping is refused by map_coords (std::out_of_range, staging's)
template <class Callback>
void call_selected_coords(const MultiIntMatrix &matrix, const annot::CoordToHeader &mapping,
                          std::span<const AnnotatedDBG::row_index> rows,
                          const HeaderCounts &counts, const Callback &callback) {
    if (rows.empty() || counts.empty())
        return;
    tsl::hopscotch_map<Header, uint32_t> selected; // index in |counts|
    selected.reserve(counts.size());
    for (size_t r = 0; r < counts.size(); ++r) {
        selected.emplace(counts[r].first, r);
    }
    call_row_tuples(matrix, rows, [&](size_t j, const MultiIntMatrix::RowTuples &row) {
        for (const auto &[col, coords] : row) {
            if (col < mapping.num_columns() && !mapping.num_sequences(col))
                continue;
            mapping.map_coords(col, coords, [&](size_t seq, uint64_t coord) {
                auto it = selected.find(Header{ col, seq });
                if (it != selected.end())
                    callback(j, it->second, coord);
            });
        }
    });
}

// The values of the records of |counts| at the positions of the repeated rows, whose
// values the fold does not hold. The value is each coordinate of the record there or,
// with |count|, their number. A value is held once per row, and found again for each
// position of the row as the positions are visited in order (call_in_order): it is not
// spread over the positions first, so nothing grows with a row's positions times its
// values, and nothing is sorted (review M1/m3 of lane T7-M5, Q-M16). A record visits the
// positions of the rows it has values at only (review A-m1). The values held, and the
// rows a record visits, are charged to |scratch| as they grow. The positions of the
// repeated rows (4 B a position and 4 B a row, made only when records are returned) are
// not: they are the query's, as its rows and positions are.
class RepeatedValues {
  public:
    RepeatedValues(const MultiIntMatrix &matrix, const annot::CoordToHeader &mapping,
                   const QueryRows &q, std::span<const size_t> positions,
                   const HeaderCounts &counts, bool count, CutScratch *scratch)
          : hits_(counts.size()), scratch_(scratch) {
        if (counts.empty())
            return;
        // the positions of each repeated row, ascending: those of row j are
        // positions_[offsets_[j], offsets_[j + 1])
        offsets_.assign(q.repeated.size() + 1, 0);
        for (size_t d = 0; d < q.distinct.size(); ++d) {
            if (q.weight[d] > 1)
                offsets_[q.repeated_of[d] + 1] = q.weight[d];
        }
        std::partial_sum(offsets_.begin(), offsets_.end(), offsets_.begin());
        positions_.resize(offsets_.back());
        for (size_t i = 0; i < q.row_of.size(); ++i) {
            const uint32_t d = q.row_of[i];
            if (q.weight[d] > 1)
                positions_[offsets_[q.repeated_of[d]]++] = positions[i];
        }
        // each offset was moved to the end of its row: the start of the next one
        std::copy_backward(offsets_.begin(), offsets_.end() - 1, offsets_.end());
        offsets_[0] = 0;
        call_selected_coords(matrix, mapping, q.repeated, counts,
                             [&](size_t j, uint32_t r, uint64_t coord) {
            auto &hits = hits_[r];
            if (count && hits.size() && hits.back().first == j) {
                ++hits.back().second;
            } else {
                scratch_->make_room(hits);
                hits.emplace_back(j, count ? 1 : coord);
            }
        });
    }

    // Call |callback|(position, value) for the values of counts[r], in ascending order of
    // position and, at a position, of value (a row's coordinates are in ascending order).
    // Once: the values are freed after. The positions of a record's one row are walked
    // directly, those of several rows merged by position (a heap of the rows).
    template <class Callback>
    void call_in_order(size_t r, const Callback &callback) {
        auto &hits = hits_[r];
        // by row, as call_selected_coords adds them, and at a row in ascending order
        assert(std::is_sorted(hits.begin(), hits.end()));
        if (hits.empty()) {
            // nothing to visit
        } else if (hits.front().first == hits.back().first) {
            // one row: its positions in order
            const uint32_t j = hits.front().first;
            for (uint32_t i = offsets_[j]; i < offsets_[j + 1]; ++i) {
                for (const auto &hit : hits) {
                    callback(positions_[i], hit.second);
                }
            }
        } else {
            merge_rows(hits, callback);
        }
        scratch_->free(hits);
    }

  private:
    // Call |callback| for the values in |hits| of several rows, merged by position: a heap
    // of the rows by their next position (no position is in two rows), charged to the
    // scratch as it grows
    template <class Callback>
    void merge_rows(const std::vector<std::pair<uint32_t, uint64_t>> &hits,
                    const Callback &callback) {
        // the next of a row's positions, and its values in |hits|
        struct Row {
            uint32_t next;
            uint32_t end;
            size_t hits_begin;
            size_t hits_end;
        };
        std::vector<Row> rows;
        for (size_t h = 0; h < hits.size(); ) {
            const uint32_t j = hits[h].first;
            size_t e = h + 1;
            while (e < hits.size() && hits[e].first == j) {
                ++e;
            }
            scratch_->make_room(rows);
            rows.push_back({ offsets_[j], offsets_[j + 1], h, e });
            h = e;
        }
        auto later = [&](const Row &a, const Row &b) {
            return positions_[a.next] > positions_[b.next];
        };
        std::make_heap(rows.begin(), rows.end(), later);
        for (size_t left = rows.size(); left; ) {
            std::pop_heap(rows.begin(), rows.begin() + left, later);
            Row &row = rows[left - 1];
            for (size_t h = row.hits_begin; h < row.hits_end; ++h) {
                callback(positions_[row.next], hits[h].second);
            }
            if (++row.next < row.end) {
                std::push_heap(rows.begin(), rows.begin() + left, later);
            } else {
                --left;
            }
        }
        scratch_->free(rows);
    }

    // the query positions of the repeated rows, by row, and the start of each row's in it
    std::vector<uint32_t> positions_;
    std::vector<uint32_t> offsets_;
    // (index of the repeated row, value) of each record, by row
    std::vector<std::vector<std::pair<uint32_t, uint64_t>>> hits_;
    CutScratch *scratch_;
};

/**
 * The sequence records (headers) that the coordinates of annotation rows fall in, with
 * the number of rows matching each record and, on request, a value per query position.
 * Rows are added one by one, so they can be decoded in blocks and dropped. The records
 * are kept in the order of their first match.
 */
class RecordFold {
  public:
    enum Values {
        kNone,
        kAbundance, // the coordinates of the record over all rows
        kPresence,  // the positions matching the record
        kCounts,    // the number of the record's coordinates at each matching position
        kCoords,    // the record's local coordinates at each matching position
    };

    // a run of query positions [pos, pos + len) with one value (kPresence, kCounts), or
    // with the coordinates value, value + 1, ... (kCoords)
    struct Run {
        uint32_t pos;
        uint32_t len;
        uint64_t value;
    };

    // |min_count|, |num_top|: the cut that will be made (top_records: staging's
    // search_top_n, which ranks the records with values by their count, as master's select()
    // did). |total_weight|: the weight of all rows to be added. With them, the fold knows the
    // records that can no longer be returned, whose per-position values it need not hold.
    RecordFold(const annot::CoordToHeader &mapping, Values values, size_t min_count = 0,
               size_t num_top = std::numeric_limits<size_t>::max(), size_t total_weight = 0)
          : mapping_(mapping), values_(values), min_count_(min_count), num_top_(num_top),
            total_weight_(total_weight) {}

    ~RecordFold() {
        folds_freed += freed_;
        folds_met_late += met_late_;
        folds_counts_only += dropped_;
    }

    // Add the row of query position |pos|, counted |weight| times, and without |values|
    // only counted. With per-position values, the positions must come in ascending
    // order, each once.
    void add_row(const MultiIntMatrix::RowTuples &row, size_t pos, size_t weight = 1,
                 bool values = true) {
        const Values values_to_add = values ? values_ : kNone;
        // runs hold 32-bit positions, and the end of a run is one past its last position
        if (pos >= std::numeric_limits<uint32_t>::max()) {
            throw std::runtime_error("A sequence queried by sequence record must have fewer "
                                     "than 2^32 k-mers");
        }
        ++row_;
        for (const auto &entry : row) {
            const Column column = entry.first;
            const auto &coords = entry.second;
            // a column without sequence records has no header labels (master's #764); a column
            // outside the mapping is refused by map_coords (std::out_of_range, staging's)
            if (column < mapping_.num_columns() && !mapping_.num_sequences(column))
                continue;
            Slots &slots = slots_[column];
            mapping_.map_coords(column, coords, [&](size_t seq_id, uint64_t local_coord) {
                uint32_t &slot = slots.at(seq_id, mapping_.num_sequences(column));
                if (!slot) {
                    if (records_.size() == std::numeric_limits<uint32_t>::max())
                        throw std::runtime_error("A query matches 2^32 or more sequence records");
                    // a record first met too late to be returned holds no values
                    const bool keeps_values = !dropped_ && remaining_weight() >= out_below_;
                    met_late_ += !dropped_ && !keeps_values;
                    records_.push_back({ Header{column, seq_id}, 0, 0, 0, {}, 0, 0, 0,
                                         keeps_values });
                    slot = records_.size();
                }
                Record &record = records_[slot - 1];
                const bool first = record.last_row != row_;
                if (first) {
                    record.last_row = row_;
                    if (!reading_again_)
                        record.count += weight;
                }
                const size_t runs_before = record.runs.size();
                switch (record.keeps_values ? values_to_add : kNone) {
                    case kNone:
                        break;
                    case kAbundance:
                        record.abundance += weight;
                        break;
                    case kPresence:
                    case kCounts:
                        if (first) {
                            close_last_run(record.runs);
                            push_run(record.runs, { static_cast<uint32_t>(pos), 1, 1 });
                        } else if (values_ == kCounts) {
                            ++record.runs.back().value;
                        }
                        break;
                    case kCoords:
                        if (first)
                            open_runs(record, pos);
                        add_coord(record, pos, local_coord);
                        break;
                }
                // a run is only merged into the one before when one is added
                num_runs_ += record.runs.size() - runs_before;
            });
        }
        seen_weight_ += weight;
        if (!reading_again_ && num_runs_ * sizeof(Run) > check_bytes_)
            limit_values();
    }

    // After select(): true if the fold kept only the counts. The values of |selected|
    // are then held by adding the rows again after hold_only(selected).
    bool dropped(const HeaderCounts &selected) const {
        // a record that lost its values alone could not be returned: one that is would be
        // returned without them (a cut that ranks by anything but the count)
        if (!dropped_ && !std::all_of(selected.begin(), selected.end(), [&](const auto &h) {
                return record(h.first).keeps_values; })) {
            throw std::logic_error("The fold of a query by record freed the values of a "
                                   "record it returns");
        }
        return dropped_;
    }

    // Hold the values of the records of |selected| only, for the rows added from now
    // on: the rows added before, added again. The counts stay. The runs added from now on
    // are charged to |scratch| before they grow (they are the values of the returned
    // records, built after the cut).
    void hold_only(const HeaderCounts &selected, CutScratch *scratch) {
        assert(!num_runs_);
        for (const auto &[header, count] : selected) {
            record(header).keeps_values = true;
        }
        reading_again_ = true;
        scratch_ = scratch;
    }

    // The records matched by at least |min_count| rows, in the order of their first match
    // (staging cuts them with search_top_n: top_records)
    HeaderCounts select(size_t min_count) const {
        HeaderCounts counts;
        for (const Record &record : records_) {
            if (record.count >= min_count)
                counts.emplace_back(record.header, record.count);
        }
        return counts;
    }

    size_t abundance(const Header &header) const { return record(header).abundance; }

    // The runs of a record: in ascending order of position, or (kCoords) in any order.
    // Runs read again after the cut were charged to the scratch: free them with it.
    std::vector<Run> take_runs(const Header &header) {
        std::vector<Run> runs = std::move(record(header).runs);
        if (values_ != kCoords)
            close_last_run(runs);
        return runs;
    }
    // the scratch the runs of take_runs() were charged to, null if none
    CutScratch* runs_scratch() const { return scratch_; }

  private:
    struct Record {
        Header header;
        size_t count;
        size_t abundance;
        size_t last_row;
        std::vector<Run> runs;
        // kCoords: the runs [open_begin, open_end) can continue at the position being
        // added, in the order of their next coordinate. Those before |open_next| were
        // continued or passed.
        uint32_t open_begin;
        uint32_t open_next;
        uint32_t open_end;
        // false once the record can no longer be returned, or for all records once the
        // fold keeps only the counts
        bool keeps_values;
    };

    // The index + 1 in |records_| of each sequence of a column, 0 for the ones not met:
    // a hash map until an eighth of the column's sequences are met, an array after that
    struct Slots {
        tsl::hopscotch_map<size_t, uint32_t> sparse;
        std::vector<uint32_t> dense;

        uint32_t& at(size_t seq_id, size_t num_sequences) {
            if (dense.empty()) {
                if (sparse.size() * 8 < num_sequences)
                    return sparse[seq_id];
                dense.assign(num_sequences, 0);
                for (const auto &[seq, slot] : sparse) {
                    dense[seq] = slot;
                }
                sparse = {};
            }
            return dense[seq_id];
        }
        uint32_t get(size_t seq_id) const {
            if (dense.size())
                return dense[seq_id];
            auto it = sparse.find(seq_id);
            return it == sparse.end() ? 0 : it->second;
        }
    };

    Record& record(const Header &header) {
        return records_[slots_.at(header.first).get(header.second) - 1];
    }
    const Record& record(const Header &header) const {
        return records_[slots_.at(header.first).get(header.second) - 1];
    }

    // Append |run| to |runs|: charged to the scratch before the runs grow when they are
    // read again after the cut (the values of a returned record)
    void push_run(std::vector<Run> &runs, const Run &run) {
        if (scratch_)
            scratch_->make_room(runs);
        runs.push_back(run);
    }

    // Join the last run with the one before it if it continues it with the same value.
    // The last run of a record is the position being added, whose value can still grow.
    static void close_last_run(std::vector<Run> &runs) {
        if (runs.size() < 2)
            return;
        Run &previous = runs[runs.size() - 2];
        const Run &last = runs.back();
        if (previous.pos + previous.len == last.pos && previous.value == last.value) {
            previous.len += last.len;
            runs.pop_back();
        }
    }

    // the coordinate that continues |run| at the position after its last one
    static uint64_t next_coord(const Run &run) { return run.value + run.len; }

    // Before the first coordinate of |pos|: find the runs that end at the position
    // before, which the coordinates of |pos| can continue. A record met at several
    // coordinates per position (a repeat) has one such run per coordinate.
    static void open_runs(Record &record, size_t pos) {
        std::vector<Run> &runs = record.runs;
        const auto open = std::partition(runs.begin() + record.open_begin, runs.end(),
                                         [&](const Run &run) {
            return static_cast<size_t>(run.pos) + run.len != pos;
        });
        std::sort(open, runs.end(), [](const Run &a, const Run &b) {
            return next_coord(a) < next_coord(b);
        });
        record.open_begin = open - runs.begin();
        record.open_next = record.open_begin;
        record.open_end = runs.size();
    }

    // Continue the open run that |coord| at |pos| follows, or start a run. The
    // coordinates of a position come in ascending order, as the open runs do.
    void add_coord(Record &record, size_t pos, uint64_t coord) {
        std::vector<Run> &runs = record.runs;
        while (record.open_next < record.open_end
                && next_coord(runs[record.open_next]) < coord) {
            ++record.open_next;
        }
        if (record.open_next < record.open_end && next_coord(runs[record.open_next]) == coord) {
            ++runs[record.open_next++].len;
        } else {
            push_run(runs, { static_cast<uint32_t>(pos), 1, coord });
        }
    }

    // the weight of the rows still to be added
    size_t remaining_weight() const {
        return total_weight_ > seen_weight_ ? total_weight_ - seen_weight_ : 0;
    }

    void drop_values(Record &record) {
        num_runs_ -= record.runs.size();
        // free the buffer: assigning {} would keep its capacity
        std::vector<Run>().swap(record.runs);
        record.open_begin = record.open_next = record.open_end = 0;
        record.keeps_values = false;
    }

    // Past the budget: free the values of the records that can no longer be returned.
    // Such a record has fewer matches, even if all rows to come match it, than
    // |min_count_| or than |num_top_| records have now (counts only grow). The values
    // of a record seen at scattered positions take a run each, for every record until
    // the cut. If the records that still can be returned are more than |num_top_| and
    // hold more than half the budget, keep only the counts: they and the order of the
    // records don't depend on the values, which are read again after the cut.
    // (Staging's cut, search_top_n, ranks by the same count, so a record freed here scores
    // below the N-th score whatever the ties.)
    void limit_values() {
        size_t out_below = min_count_;
        if (num_top_ && records_.size() >= num_top_) {
            std::vector<size_t> counts;
            counts.reserve(records_.size());
            for (const Record &record : records_) {
                counts.push_back(record.count);
            }
            std::nth_element(counts.begin(), counts.begin() + (num_top_ - 1), counts.end(),
                             std::greater<size_t>());
            out_below = std::max(out_below, counts[num_top_ - 1]);
        }
        out_below_ = std::max(out_below_, out_below);
        const size_t remaining = remaining_weight();
        size_t num_kept = 0;
        for (Record &record : records_) {
            if (record.keeps_values && record.count + remaining < out_below_) {
                drop_values(record);
                ++freed_;
            }
            num_kept += record.keeps_values;
        }
        const size_t bytes = num_runs_ * sizeof(Run);
        if (num_kept > num_top_ && bytes > budget_ / 2) {
            for (Record &record : records_) {
                drop_values(record);
            }
            dropped_ = true;
        } else {
            // check again when the values have doubled
            check_bytes_ = std::max(budget_, 2 * bytes);
        }
    }

    const annot::CoordToHeader &mapping_;
    const Values values_;
    const size_t min_count_;
    const size_t num_top_;
    const size_t total_weight_;
    size_t seen_weight_ = 0;
    // the budget (fold_values_budget when the fold is made), the runs held, and the bytes
    // of runs past which limit_values() is called
    const size_t budget_ = fold_values_budget;
    size_t num_runs_ = 0;
    size_t check_bytes_ = budget_;
    // a record with fewer matches than this at the end is not returned
    size_t out_below_ = 0;
    // the fold kept only the counts
    bool dropped_ = false;
    // records whose values limit_values() freed, and records met too late (the counts)
    uint64_t freed_ = 0;
    uint64_t met_late_ = 0;
    // the rows are added a second time, for the values of the returned records
    bool reading_again_ = false;
    // what the runs added a second time are charged to (hold_only)
    CutScratch *scratch_ = nullptr;
    size_t row_ = 0;
    std::vector<Record> records_;
    tsl::hopscotch_map<Column, Slots> slots_;
};

// If |fold| kept only the counts, read the rows of one position again to hold the
// values of the records of |counts|, charged to |scratch| as they grow
void read_again_if_dropped(const MultiIntMatrix &matrix, const QueryRows &q,
                           const HeaderCounts &counts, RecordFold *fold, CutScratch *scratch) {
    if (counts.empty() || !fold->dropped(counts))
        return;
    std::vector<AnnotatedDBG::row_index> rows;
    std::vector<size_t> positions;
    for (size_t i = 0; i < q.distinct.size(); ++i) {
        if (q.weight[i] == 1) {
            rows.push_back(q.distinct[i]);
            positions.push_back(q.first_position[i]);
        }
    }
    fold->hold_only(counts, scratch);
    call_row_tuples(matrix, rows, [&](size_t i, const MultiIntMatrix::RowTuples &row) {
        fold->add_row(row, positions[i]);
    });
}

// A record of the fold returned by top_records, with its label
struct TopRecord {
    Header header;
    size_t count;
    size_t abundance;
    std::string label;
};

/**
 * The records |fold| met at |min_count| query positions or more, the top |num_top_labels|
 * of them in the C5 order (search_top_n) of the score reported: the number of their
 * coordinates matched if |by_coordinates| (the fold's abundance), of the query positions
 * otherwise; ties by the header's bytes (coordinates are positions, not abundances), a
 * header in two columns in (column, seq_id) order. |cut| states the cut. Nothing is
 * charged: the cut is made before the first charge (train 7 P1's R8 (c)).
 */
std::vector<TopRecord> top_records(const RecordFold &fold, const annot::CoordToHeader &mapping,
                                   size_t min_count, size_t num_top_labels,
                                   bool by_coordinates, SearchCut *cut) {
    // the records in (column, seq_id) order, so that a header in two columns is returned in
    // the same order by every request
    HeaderCounts candidates = fold.select(min_count);
    std::sort(candidates.begin(), candidates.end());

    std::vector<size_t> scores;
    scores.reserve(candidates.size());
    for (const auto &[header, count] : candidates) {
        scores.push_back(by_coordinates ? fold.abundance(header) : count);
    }
    auto top = search_top_n(scores, num_top_labels, [&](size_t i) -> const std::string & {
        const auto &[column, seq_id] = candidates[i].first;
        return mapping.get_headers(column)[seq_id];
    }, {}, cut);

    std::vector<TopRecord> result;
    result.reserve(top.size());
    for (auto &[i, label] : top) {
        const auto &[header, count] = candidates[i];
        result.push_back({ header, count, by_coordinates ? fold.abundance(header) : count,
                           std::move(label) });
    }
    return result;
}

// The records of |top| with their counts, as the fold's helpers take them
HeaderCounts header_counts(const std::vector<TopRecord> &top) {
    HeaderCounts counts;
    counts.reserve(top.size());
    for (const TopRecord &record : top) {
        counts.emplace_back(record.header, record.count);
    }
    return counts;
}
} // namespace

size_t record_fold_values_budget() { return fold_values_budget; }
void set_record_fold_values_budget(size_t bytes) { fold_values_budget = bytes; }
RecordFoldCounts record_fold_counts() {
    return { folds_freed, folds_met_late, folds_counts_only };
}
void reset_record_fold_counts() {
    folds_freed = 0;
    folds_met_late = 0;
    folds_counts_only = 0;
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
                                 bool with_kmer_counts, SearchCut *cut,
                                 SearchAccount *values_account) const {
    static_assert(utils::is_one_of_v<Values, size_t, CountRuns, PresenceRuns>);
    assert((!with_kmer_counts || std::is_same_v<Values, size_t>));
    assert(coord_to_header_);
    const size_t num_kmers = sequence.size() - dbg_.get_k() + 1;
    // The rows to decode and the query position of each. Without per-position values,
    // a repeated k-mer is one row counted once per position.
    std::vector<row_index> rows;
    std::vector<size_t> positions;
    std::vector<size_t> weights;
    size_t num_present = 0;
    if constexpr (std::is_same_v<Values, size_t>) {
        VectorMap<row_index, size_t> counts;
        counts.reserve(num_kmers);
        for (node_index node : dbg_.node_range(sequence)) {
            if (node) {
                ++counts[graph_to_anno_index(node)];
                ++num_present;
            }
        }
        rows.reserve(counts.size());
        weights.reserve(counts.size());
        for (const auto &[row, count] : counts.values_container()) {
            rows.push_back(row);
            weights.push_back(count);
        }
    } else {
        map_to_rows(dbg_, sequence, &rows, &positions);
        num_present = rows.size();
    }
    const size_t min_count = get_min_count(discovery_fraction, presence_fraction,
                                           num_kmers, num_present);
    if (num_present < min_count)
        return {};

    const auto *matrix = dynamic_cast<const MultiIntMatrix *>(&annotator_->get_matrix());
    if (!matrix)
        throw std::runtime_error("k-mer coordinates are not indexed in this annotator");
    // the constructor refused a mapping of another column count: a request never exits here
    if (annotator_->num_labels() != coord_to_header_->num_columns()) {
        throw std::logic_error(fmt::format(
                "Incompatible number of columns in CoordToHeader and annotation matrix: {} != {}",
                coord_to_header_->num_columns(), annotator_->num_labels()));
    }
    RecordFold fold(*coord_to_header_,
                    std::is_same_v<Values, PresenceRuns> ? RecordFold::kPresence
                        : std::is_same_v<Values, CountRuns> ? RecordFold::kCounts
                        : with_kmer_counts ? RecordFold::kAbundance : RecordFold::kNone,
                    min_count, num_top_labels, positions.size());
    if constexpr (std::is_same_v<Values, size_t>) {
        call_row_tuples(*matrix, rows, [&](size_t i, const MultiIntMatrix::RowTuples &row) {
            fold.add_row(row, 0, weights[i]);
        });
        // the records ranked by the score reported: with |with_kmer_counts| the number of
        // their coordinates matched
        auto top = top_records(fold, *coord_to_header_, min_count, num_top_labels,
                               with_kmer_counts, cut);
        std::vector<std::tuple<Label, size_t, Values>> result;
        result.reserve(top.size());
        for (TopRecord &record : top) {
            result.emplace_back(std::move(record.label), record.count, record.abundance);
        }
        return result;

    } else {
        // the positions of the repeated rows are kept as 32-bit, as the runs of the fold
        if (num_kmers >= std::numeric_limits<uint32_t>::max()) {
            throw std::runtime_error("A sequence queried by sequence record must have fewer "
                                     "than 2^32 k-mers");
        }
        // As in get_kmer_coordinates, each row is added once, at its first position,
        // counted for all its positions. A repeated row is only counted here, and its
        // values are added below for the returned records.
        const QueryRows q(rows, positions);
        call_row_tuples(*matrix, q.distinct, [&](size_t i, const MultiIntMatrix::RowTuples &row) {
            fold.add_row(row, q.first_position[i], q.weight[i], q.weight[i] == 1);
        });
        // the cut, before anything is charged
        auto top = top_records(fold, *coord_to_header_, min_count, num_top_labels, false, cut);
        const HeaderCounts counts = header_counts(top);

        // the records' values, charged to the unit's values account (if any) as they grow,
        // and what they are built from after the cut, as its scratch
        charge_search(values_account, buffer_bytes(top.size(), sizeof(Values)));
        CutScratch scratch(values_account);
        read_again_if_dropped(*matrix, q, counts, &fold, &scratch);
        // the values of the returned records at the repeated rows
        RepeatedValues unheld_values(*matrix, *coord_to_header_, q, positions, counts, true,
                                     &scratch);

        std::vector<std::tuple<Label, size_t, Values>> result;
        result.reserve(top.size());
        for (size_t r = 0; r < top.size(); ++r) {
            Values values(num_kmers, values_account);
            // merge the runs of the fold and the values it does not hold (at the positions
            // of the repeated rows, which no run holds), both in ascending order of position
            auto runs = fold.take_runs(top[r].header);
            auto next = runs.begin();
            auto add_runs_before = [&](size_t pos) {
                for ( ; next != runs.end() && next->pos < pos; ++next) {
                    if constexpr (std::is_same_v<Values, PresenceRuns>) {
                        values.add_run(next->pos, next->pos + next->len);
                    } else {
                        values.add_run(next->pos, next->pos + next->len, next->value);
                    }
                }
            };
            unheld_values.call_in_order(r, [&](uint32_t pos, uint64_t value) {
                add_runs_before(pos);
                if constexpr (std::is_same_v<Values, PresenceRuns>) {
                    values.add(pos);
                    (void)value;
                } else {
                    values.add(pos, value);
                }
            });
            add_runs_before(num_kmers);
            if (fold.runs_scratch())
                scratch.free(runs);
            values.shrink();
            result.emplace_back(std::move(top[r].label), top[r].count, std::move(values));
        }
        return result;
    }
}


template <class Values>
uint64_t search_store_bytes(size_t num_labels, size_t num_rows) {
    using Store = utils::ValueStore<Column, Values>;
    if (Store::dense(num_labels)) {
        // ids_: a size_t per label; dense_values_: reserved for the rows
        return buffer_bytes(num_labels, sizeof(size_t)) + buffer_bytes(num_rows, sizeof(Values));
    }
    if (!num_rows)
        return 0;
    // tsl's reserve: the buckets for |num_rows| at the table's maximum load, rounded up by its
    // growth policy (the table holds them and its padding, HopscotchLayout)
    using Map = typename Store::SparseMap;
    size_t buckets = size_t(std::ceil(float(num_rows) / Map().max_load_factor()));
    tsl::hh::power_of_two_growth_policy<2> policy(buckets);
    return annot::matrix::HopscotchLayout<Map>::array_bytes(buckets);
}
template uint64_t search_store_bytes<PresenceRuns>(size_t, size_t);
template uint64_t search_store_bytes<CountRuns>(size_t, size_t);
template uint64_t search_store_bytes<PositionCoords>(size_t, size_t);

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
                            SearchCut *cut,
                            SearchAccount *values_account = nullptr) {
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
    // the store, this record's scratch (freed when its rows are returned; the values, which
    // charge their buffers as they grow, move into them), charged to the unit's values account
    // (if any) as what it allocates before it is built: the dense index over the labels and the
    // values' buffer, or the hash table reserved for the labels kept (REVIEW 9726cb420 [P2],
    // lane T6-W2: the smaller of the two was charged, but the store picks its form by the
    // labels alone)
    charge_search_scratch(values_account,
                          search_store_bytes<ValueType>(label_encoder.size(), counts.size()));
    utils::ValueStore<Column, ValueType> value_store(label_encoder.size(), counts.size());
    for (const auto &[j, count] : counts) {
        value_store.emplace(j, ValueType(num_kmers, values_account));
    }
    // keys that collide in the table's hash (label ids chosen for it) can make it grow past its
    // reserve: what it holds then, charged once its keys are in, before any value is added
    if (values_account && !value_store.is_dense())
        values_account->charge_scratch(
                annot::matrix::hopscotch_held_bytes(value_store.sparse_values()));

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
                              SearchCut *cut,
                              SearchAccount *values) const {
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
                                             false, cut, values);
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
            medians, cut, values);
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
                                   SearchCut *cut,
                                   SearchAccount *values) const {
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

    if (coord_to_header_) {
        // the constructor refused a mapping of another column count: a request never exits
        if (annotator_->num_labels() != coord_to_header_->num_columns()) {
            throw std::logic_error(fmt::format(
                    "Incompatible number of columns in CoordToHeader and annotation matrix: "
                    "{} != {}", coord_to_header_->num_columns(), annotator_->num_labels()));
        }
        // the rows and positions of a query are indexed as 32-bit
        if (num_kmers >= std::numeric_limits<uint32_t>::max()) {
            throw std::runtime_error("A sequence queried for coordinates must have fewer "
                                     "than 2^32 k-mers");
        }
        // each distinct row is decoded once, not once per position
        const QueryRows q(row_indices, kmer_positions);

        // Split the coordinates across sequences, a block of rows at a time. Each row is
        // added once, at its first position, counted for all its positions. A row of one
        // position adds its coordinates. A repeated row is only counted here: its
        // coordinates at all positions, for all records, would take time and memory in
        // proportion to both. They are added below for the returned records.
        RecordFold fold(*coord_to_header_, RecordFold::kCoords, min_count, num_top_labels,
                        kmer_positions.size());
        call_row_tuples(*tuple_matrix, q.distinct,
                        [&](size_t i, const MultiIntMatrix::RowTuples &row) {
            fold.add_row(row, q.first_position[i], q.weight[i], q.weight[i] == 1);
        });
        // the records ranked by the number of query positions they match, cut before
        // anything is charged
        auto top = top_records(fold, *coord_to_header_, min_count, num_top_labels, false, cut);
        const HeaderCounts counts = header_counts(top);

        // the records' coordinates, charged to the unit's values account (if any) as they
        // grow, and what they are built from after the cut, as its scratch
        charge_search(values, buffer_bytes(top.size(), sizeof(PositionCoords)));
        CutScratch scratch(values);
        read_again_if_dropped(*tuple_matrix, q, counts, &fold, &scratch);
        // (position, coordinate) of each returned record at the repeated rows
        RepeatedValues unheld_coords(*tuple_matrix, *coord_to_header_, q, kmer_positions,
                                     counts, false, &scratch);

        std::vector<std::tuple<std::string, size_t, PositionCoords, size_t>> result;
        result.reserve(top.size());
        for (size_t r = 0; r < top.size(); ++r) {
            const auto &[col, header] = top[r].header;
            // (position, coordinate) of the fold's runs, sorted: a record at several
            // coordinates of a position has a run for each, and the runs come in any order
            std::vector<std::pair<uint32_t, uint64_t>> held;
            {
                auto runs = fold.take_runs(top[r].header);
                size_t num_held = 0;
                for (const RecordFold::Run &run : runs) {
                    num_held += run.len;
                }
                scratch.reserve(held, num_held);
                for (const RecordFold::Run &run : runs) {
                    for (uint32_t i = 0; i < run.len; ++i) {
                        held.emplace_back(run.pos + i, run.value + i);
                    }
                }
                if (fold.runs_scratch())
                    scratch.free(runs);
            }
            std::sort(held.begin(), held.end());
            // merged with the coordinates at the repeated rows, whose positions no run
            // holds: each position's coordinates in ascending order, as the row holds a
            // record's
            PositionCoords coords(num_kmers, values);
            auto next = held.begin();
            auto add_held_before = [&](size_t pos) {
                for ( ; next != held.end() && next->first < pos; ++next) {
                    coords.add(next->first, next->second);
                }
            };
            unheld_coords.call_in_order(r, [&](uint32_t pos, uint64_t coord) {
                add_held_before(pos);
                coords.add(pos, coord);
            });
            add_held_before(num_kmers);
            scratch.free(held);
            coords.shrink();
            result.emplace_back(std::move(top[r].label), top[r].count, std::move(coords),
                                coord_to_header_->num_kmers_in_sequence(col, header));
        }
        return result;
    }

    // Without a record mapping: each distinct row is decoded once, its positions still counted
    // separately: a row was decoded and held once per query position, so a k-mer the query
    // repeats m times held m copies of its coordinates (8 B x m x the row's coordinates: GBs
    // for a homopolymer query, T5b-Y3). |rows|: the distinct rows in the order of their first
    // position; |row_of|: the index in |rows| of each present position's row
    std::vector<row_index> rows;
    std::vector<size_t> row_of(row_indices.size());
    {
        tsl::hopscotch_map<row_index, size_t> index_of;
        index_of.reserve(row_indices.size());
        for (size_t i = 0; i < row_indices.size(); ++i) {
            auto [it, inserted] = index_of.emplace(row_indices[i], rows.size());
            if (inserted)
                rows.push_back(row_indices[i]);
            row_of[i] = it->second;
        }
    } // the hash table is released before the rows are decoded
    row_indices = std::vector<row_index>();

    // a row's content does not depend on the rows read with it, so that each position gets the
    // row it got when every position's row was read
    auto rows_tuples = tuple_matrix->get_row_tuples(rows);

    // each present position's row, in the order of the positions
    auto call_rows_tuples = [&](const auto &callback) {
        for (size_t d : row_of) {
            callback(rows_tuples[d]);
        }
    };

    // ties ordered by label: the coordinates are positions, not abundances
    auto results = filter_and_aggregate<std::vector<std::tuple<Label, size_t, PositionCoords>>>(
        call_rows_tuples, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers,
        {}, cut, values);
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
                                       SearchCut *cut,
                                       SearchAccount *values) const {
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
        // every position marked: one run, a copy per label (each charged to the unit's values
        // account, if any)
        PresenceRuns all(num_kmers, values);
        all.add_run(0, num_kmers);
        all.shrink();
        charge_search(values, label_counts.size() * alloc_bytes(all.heap_bytes()));
        for (auto&& [label, count] : label_counts) {
            presence_vectors.emplace_back(std::move(label), count, all);
        }
        return presence_vectors;
    }

    if (coord_to_header_) {
        return get_header_results<PresenceRuns>(sequence, num_top_labels,
                                                discovery_fraction, presence_fraction,
                                                false, cut, values);
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
            medians, cut, values);
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
