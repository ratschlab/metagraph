#include "annotated_dbg.hpp"

#include <cstdlib>
#include <numeric>

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
        kmer_abundances_(holds_kmer_abundances(*annotator_)) {}

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
        auto kmer_coord_res = get_kmer_coordinates(sequence, std::numeric_limits<size_t>::max(),
                                                   discovery_fraction, presence_fraction);
        std::vector<Label> result;
        result.reserve(kmer_coord_res.size());
        for (auto &[label, count, coords, num_kmers_in_target] : kmer_coord_res) {
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
             const std::function<std::string(size_t)> &label,
             const std::function<std::vector<uint64_t>(const std::vector<size_t> &)> &medians,
             SearchCut *cut) {
    SearchCut c;
    c.abundances = static_cast<bool>(medians);

    struct Row {
        size_t index;
        size_t score;
        uint64_t median;
        std::string label;
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

    std::vector<Row> rows;
    rows.reserve(order.size());
    for (size_t i : order) {
        rows.push_back(Row{ i, scores[i], 0, std::string() });
    }

    // The medians break ties only, so a row's median is read only where it is compared: the
    // rows that share their score with another one, the tie group at a cut (whose N-th row's
    // median the cut then states) and the rows above it of one score. A row with a score of
    // its own keeps median 0, which no comparison reads, and a cut whose N-th score has one
    // row states no median. Reading the abundances costs a pass over the query's rows and
    // holds values of the columns asked, so it is asked for no other row.
    if (medians && rows.size()) {
        std::vector<size_t> asked;
        // the rows above the N-th score (every row of a complete result) in score order: the
        // runs of one score
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
        // the tie group at the cut, if the N-th score is shared
        if (rows.size() - above > 1) {
            for (size_t r = above; r < rows.size(); ++r) {
                asked.push_back(r);
            }
            c.nth_median_read = true;
        }
        if (asked.size()) {
            std::vector<size_t> indexes(asked.size());
            std::transform(asked.begin(), asked.end(), indexes.begin(),
                           [&](size_t r) { return rows[r].index; });
            std::vector<uint64_t> values = medians(indexes);
            assert(values.size() == asked.size());
            for (size_t a = 0; a < asked.size(); ++a) {
                rows[asked[a]].median = values[a];
            }
        }
    }

    auto before = [](const Row &a, const Row &b) {
        if (search_row_before({ a.score, a.median, a.label }, { b.score, b.median, b.label }))
            return true;
        if (search_row_before({ b.score, b.median, b.label }, { a.score, a.median, a.label }))
            return false;
        // one label in two candidates (a header in two columns): their order in |scores|
        return a.index < b.index;
    };

    if (!c.complete && slots) {
        // The tie group competes for |slots| places: by median where its medians were read,
        // then by label. Only the rows of the N-th row's median compete by label, so only
        // they are decoded beside the rows returned (a large tie group, a common k-mer set at
        // the cut, is not decoded entirely where the medians separate it), and the places
        // are taken by a selection, not a sort of the whole group: O(T + slots log slots),
        // the returned rows sorted below.
        auto tie_begin = rows.begin() + above;
        auto same_begin = tie_begin;
        auto same_end = rows.end();
        if (c.nth_median_read) {
            // the |slots|-th largest median, and the group split around it
            std::nth_element(tie_begin, tie_begin + slots - 1, rows.end(),
                             [](const Row &a, const Row &b) { return a.median > b.median; });
            c.nth_median = tie_begin[slots - 1].median;
            same_begin = std::partition(tie_begin, rows.end(),
                    [&](const Row &row) { return row.median > c.nth_median; });
            same_end = std::partition(same_begin, rows.end(),
                    [&](const Row &row) { return row.median == c.nth_median; });
        }
        c.tied_same_median = same_end - same_begin;
        for (auto it = tie_begin; it != same_end; ++it) {
            it->label = label(it->index);
        }
        // the first rows of the N-th row's median in the C5 order take the places left
        auto end_returned = tie_begin + slots;
        assert(same_begin < end_returned && end_returned <= same_end);
        std::nth_element(same_begin, end_returned - 1, same_end, before);
        rows.erase(end_returned, rows.end());
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
    for (auto &row : rows) {
        result.emplace_back(row.index, std::move(row.label));
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
    return abundance_medians([&](const auto &callback) {
        std::vector<AnnotatedDBG::row_index> rows;
        for (size_t begin = 0; begin < row_counts.size(); begin += rows_per_read) {
            const size_t end = std::min(begin + rows_per_read, row_counts.size());
            rows.clear();
            for (size_t i = begin; i < end; ++i) {
                rows.push_back(row_counts[i].first);
            }
            auto row_values = matrix.get_row_values(rows);
            for (size_t i = begin; i < end; ++i) {
                callback(row_values[i - begin], row_counts[i].second);
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
                            [&](size_t i) { return label_encoder.decode(code_counts[i].first); },
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
        // their coordinates matched, which the counts below sum up
        auto kmer_coord_res = query_kmer_coordinates(sequence, num_top_labels,
                                                     discovery_fraction, presence_fraction,
                                                     cut, with_kmer_counts);
        std::vector<StringCountPair> result;
        result.reserve(kmer_coord_res.size());
        for (auto &[label, count, coords, num_kmers_in_target] : kmer_coord_res) {
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
                    std::vector<std::tuple<Label, size_t, sdsl::bit_vector>>,
                    std::vector<std::tuple<Label, size_t, std::vector<size_t>>>,
                    std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>>>);

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
                            [&](size_t i) { return label_encoder.decode(candidates[i].first); },
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
                if constexpr(std::is_same_v<ValueType, sdsl::bit_vector>) {
                    (*value_p)[kmer_positions[i]] = true;
                } else {
                    (*value_p)[kmer_positions[i]] = std::move(item.second);
                }
            }
        }
        i++;
    });

    Result result;
    result.reserve(counts.size());
    for (size_t r = 0; r < counts.size(); ++r) {
        const auto &[j, count] = counts[r];
        result.emplace_back(std::move(top[r].second), count, std::move(*value_store.get(j)));
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

#ifndef NDEBUG
// A sanity check to ensure the same matches are for get_top_label_signatures and get_top_labels
template <typename ValueType>
bool same_results(const std::vector<std::tuple<Label, size_t, ValueType>> &result,
                  const std::vector<StringCountPair> &top_labels) {
    assert(top_labels.size() == result.size());
    std::unordered_map<Label, uint64_t> check(top_labels.begin(), top_labels.end());
    for (const auto &[label, count, value] : result) {
        if constexpr(std::is_same_v<ValueType, sdsl::bit_vector>) {
            assert(count == sdsl::util::cnt_one_bits(value));
        } else if constexpr(std::is_same_v<ValueType, std::vector<size_t>>) {
            assert(count == (size_t)std::count_if(value.begin(), value.end(), [](const auto &v) { return v > 0; }));
        } else if constexpr(std::is_same_v<ValueType, std::vector<SmallVector<uint64_t>>>) {
            assert(count == (size_t)std::count_if(value.begin(), value.end(), [](const auto &v) { return v.size() > 0; }));
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

std::vector<std::tuple<std::string, size_t, std::vector<size_t>>>
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
        auto kmer_coord_res = get_kmer_coordinates(sequence, num_top_labels,
                                                   discovery_fraction, presence_fraction, cut);
        std::vector<std::tuple<std::string, size_t, std::vector<size_t>>> result;
        result.reserve(kmer_coord_res.size());
        for (auto &[label, count, coords, num_kmers_in_target] : kmer_coord_res) {
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
        auto results = filter_and_aggregate<std::vector<std::tuple<std::string, size_t, std::vector<size_t>>>>(
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

std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>, size_t>>
AnnotatedDBG::get_kmer_coordinates(std::string_view sequence,
                                   size_t num_top_labels,
                                   double discovery_fraction,
                                   double presence_fraction,
                                   SearchCut *cut) const {
    return query_kmer_coordinates(sequence, num_top_labels, discovery_fraction,
                                  presence_fraction, cut, false);
}

std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>, size_t>>
AnnotatedDBG::query_kmer_coordinates(std::string_view sequence,
                                     size_t num_top_labels,
                                     double discovery_fraction,
                                     double presence_fraction,
                                     SearchCut *cut,
                                     bool score_by_coordinates) const {
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
        // Now each coord in `rows_tuples` is `coord = local_coord * num_sequences[j] + seq_id`

        // Before splitting coords across sequences, filter them by the number of k-mer matches
        using Header = std::pair<Column, size_t>; // (column, seq_id)
        using Count = size_t;
        // per record: the query's k-mers it matches, and (|score_by_coordinates|) the number
        // of its coordinates they match
        std::vector<std::pair<Header, std::pair<Count, uint64_t>>> candidates;
        {
            VectorMap<Header, std::pair<Count, uint64_t>> counts_map;
            tsl::hopscotch_set<Header> matches;
            for (const auto &row : rows_tuples) {
                matches.clear();
                for (auto &[col, coords] : row) {
                    for (uint64_t coord : coords) {
                        size_t header = coord % coord_to_header_->num_sequences(col);
                        matches.emplace(col, header);
                        if (score_by_coordinates)
                            counts_map[Header(col, header)].second++;
                    }
                }
                for (const Header &h : matches) {
                    counts_map[h].first++;
                }
            }
            candidates = to_vector(std::move(counts_map));
        }

        // remove sequences with insufficient number of k-mer matches
        candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                                        [&](const auto &c) { return c.second.first < min_count; }),
                         candidates.end());
        // the records in (column, seq_id) order, so that a header in two columns is
        // returned in the same order by every request
        std::sort(candidates.begin(), candidates.end(), utils::LessFirst());

        // keep only the top-n sequences, by the score reported, ties by their header (C5:
        // coordinates are positions, not abundances)
        std::vector<size_t> scores;
        scores.reserve(candidates.size());
        for (const auto &[h, c] : candidates) {
            scores.push_back(score_by_coordinates ? c.second : c.first);
        }
        auto top = search_top_n(scores, num_top_labels, [&](size_t i) {
            const auto &[col, header] = candidates[i].first;
            return coord_to_header_->get_headers(col)[header];
        }, {}, cut);
        std::vector<std::pair<Header, Count>> counts;
        counts.reserve(top.size());
        for (const auto &[i, label] : top) {
            counts.emplace_back(candidates[i].first, candidates[i].second.first);
        }
        candidates = decltype(candidates)();

        // split coordinates across sequences
        tsl::hopscotch_map<Header, std::vector<Tuple>> coords_map;
        coords_map.reserve(counts.size());
        for (const auto &[h, count] : counts) {
            coords_map.emplace(h, num_kmers);
        }
        for (size_t i = 0; i < rows_tuples.size(); ++i) {
            for (auto &[col, coords] : rows_tuples[i]) {
                for (uint64_t coord : coords) {
                    size_t header = coord % coord_to_header_->num_sequences(col);
                    uint64_t local_coord = coord / coord_to_header_->num_sequences(col);
                    auto it = coords_map.find({col, header});
                    if (it != coords_map.end())
                        it.value()[kmer_positions[i]].push_back(local_coord);
                }
                coords.clear();
            }
            rows_tuples[i].clear();
        }
        rows_tuples.clear();

        std::vector<std::tuple<std::string, Count, std::vector<Tuple>, size_t>> result;
        result.reserve(counts.size());
        for (size_t r = 0; r < counts.size(); ++r) {
            const auto &[h, count] = counts[r];
            const auto &[col, header] = h;
            auto &coords = coords_map[h];
            result.emplace_back(std::move(top[r].second),
                                count, std::move(coords),
                                coord_to_header_->num_kmers_in_sequence(col, header));
        }

        return result;
    }

    auto call_rows_tuples = [&](const auto &callback) {
        std::for_each(rows_tuples.begin(), rows_tuples.end(), callback);
    };

    // ties ordered by label: the coordinates are positions, not abundances
    auto results = filter_and_aggregate<std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>>>>(
        call_rows_tuples, annotator_->get_label_encoder(), min_count, num_top_labels, kmer_positions, num_kmers,
        {}, cut);
    assert(same_results(results, get_top_labels(sequence, num_top_labels, discovery_fraction, presence_fraction)));

    // Without a CoordToHeader mapping, the label refers to a whole annotation
    // column rather than a single indexed sequence, so num_kmers_in_target is
    // undefined.
    std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>, size_t>> with_lengths;
    with_lengths.reserve(results.size());
    for (auto &[label, count, coords] : results) {
        with_lengths.emplace_back(std::move(label), count, std::move(coords), 0);
    }
    return with_lengths;
}

std::vector<std::tuple<Label, size_t, sdsl::bit_vector>>
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
                                                   discovery_fraction, presence_fraction, cut);
        std::vector<std::tuple<Label, size_t, sdsl::bit_vector>> result;
        result.reserve(kmer_coord_res.size());
        for (auto &[label, count, coords, num_kmers_in_target] : kmer_coord_res) {
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

    const auto &matrix = annotator_->get_matrix();
    // the abundances of the ties: the rows enumerated below may be binary, so they are read
    // from the integer matrix (one value per position), only where a tie needs them
    ColumnMedians medians = tie_medians(sequence, row_indices);
    auto get_results = [&](const auto &enumerate_rows) {
        auto results = filter_and_aggregate<std::vector<std::tuple<std::string, size_t, sdsl::bit_vector>>>(
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
