#include "annotated_dbg.hpp"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <span>

#include <fmt/format.h>

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

// Annotation rows decoded at a time for the queries by sequence record. A row holds a
// coordinate for every record with its k-mer, so the rows of a whole query can take far
// more memory than the result. Row-diff rows share their paths within a block, which a
// small block loses.
constexpr size_t kRowBlock = 256;

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

    RecordFold(const annot::CoordToHeader &mapping, Values values)
          : mapping_(mapping), values_(values) {}

    // Add the row of query position |pos|, counted |weight| times. With per-position
    // values, the positions must come in ascending order, each once.
    void add_row(const MultiIntMatrix::RowTuples &row, size_t pos, size_t weight = 1) {
        // runs hold 32-bit positions, and the end of a run is one past its last position
        if (pos >= std::numeric_limits<uint32_t>::max()) {
            throw std::runtime_error("A sequence queried by sequence record must have fewer "
                                     "than 2^32 k-mers");
        }
        ++row_;
        for (const auto &entry : row) {
            const Column column = entry.first;
            const auto &coords = entry.second;
            if (column >= mapping_.num_columns()) {
                throw std::runtime_error(fmt::format("Column {} is out of range: CoordToHeader "
                                                     "has {} columns", column, mapping_.num_columns()));
            }
            // a column without sequence records has no header labels
            if (!mapping_.num_sequences(column))
                continue;
            Slots &slots = slots_[column];
            mapping_.map_coords(column, coords, [&](size_t seq_id, uint64_t local_coord) {
                uint32_t &slot = slots.at(seq_id, mapping_.num_sequences(column));
                if (!slot) {
                    if (records_.size() == std::numeric_limits<uint32_t>::max())
                        throw std::runtime_error("A query matches 2^32 or more sequence records");
                    records_.push_back({ Header{column, seq_id}, 0, 0, 0, {}, 0, 0, 0 });
                    slot = records_.size();
                }
                Record &record = records_[slot - 1];
                const bool first = record.last_row != row_;
                if (first) {
                    record.last_row = row_;
                    record.count += weight;
                }
                switch (values_) {
                    case kNone:
                        break;
                    case kAbundance:
                        record.abundance += weight;
                        break;
                    case kPresence:
                    case kCounts:
                        if (first) {
                            close_last_run(record.runs);
                            record.runs.push_back({ static_cast<uint32_t>(pos), 1, 1 });
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
            });
        }
    }

    // The records matched by at least |min_count| rows: all of them in the order of
    // their first match, or the |num_top| with the most matches, sorted
    HeaderCounts select(size_t min_count, size_t num_top, size_t *num_passed = nullptr) const {
        HeaderCounts counts;
        for (const Record &record : records_) {
            if (record.count >= min_count)
                counts.emplace_back(record.header, record.count);
        }
        if (num_passed)
            *num_passed = counts.size();
        if (counts.size() > num_top)
            top_n_sorted(counts, num_top);
        return counts;
    }

    size_t abundance(const Header &header) const { return record(header).abundance; }

    // The runs of a record: in ascending order of position, or (kCoords) in any order
    std::vector<Run> take_runs(const Header &header) {
        std::vector<Run> runs = std::move(record(header).runs);
        if (values_ != kCoords)
            close_last_run(runs);
        return runs;
    }

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
    static void add_coord(Record &record, size_t pos, uint64_t coord) {
        std::vector<Run> &runs = record.runs;
        while (record.open_next < record.open_end
                && next_coord(runs[record.open_next]) < coord) {
            ++record.open_next;
        }
        if (record.open_next < record.open_end && next_coord(runs[record.open_next]) == coord) {
            ++runs[record.open_next++].len;
        } else {
            runs.push_back({ static_cast<uint32_t>(pos), 1, coord });
        }
    }

    const annot::CoordToHeader &mapping_;
    const Values values_;
    size_t row_ = 0;
    std::vector<Record> records_;
    tsl::hopscotch_map<Column, Slots> slots_;
};
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
    if (annotator_->num_labels() != coord_to_header_->num_columns()) {
        throw std::runtime_error(fmt::format("Incompatible number of columns in CoordToHeader "
                                             "and annotation matrix: {} != {}",
                                             coord_to_header_->num_columns(),
                                             annotator_->num_labels()));
    }
    RecordFold fold(*coord_to_header_,
                    std::is_same_v<Values, PresenceRuns> ? RecordFold::kPresence
                        : std::is_same_v<Values, CountRuns> ? RecordFold::kCounts
                        : with_kmer_counts ? RecordFold::kAbundance : RecordFold::kNone);
    call_row_tuples(*matrix, rows, [&](size_t i, const MultiIntMatrix::RowTuples &row) {
        if constexpr (std::is_same_v<Values, size_t>) {
            fold.add_row(row, 0, weights[i]);
        } else {
            fold.add_row(row, positions[i]);
        }
    });
    const auto counts = fold.select(min_count, num_top_labels, num_passed);

    std::vector<std::tuple<Label, size_t, Values>> result;
    result.reserve(counts.size());
    for (const auto &[h, count] : counts) {
        const Label &label = coord_to_header_->get_headers(h.first)[h.second];
        if constexpr (std::is_same_v<Values, size_t>) {
            result.emplace_back(label, count, with_kmer_counts ? fold.abundance(h) : count);
        } else {
            Values values(num_kmers);
            for (const RecordFold::Run &run : fold.take_runs(h)) {
                if constexpr (std::is_same_v<Values, PresenceRuns>) {
                    values.add_run(run.pos, run.pos + run.len);
                } else {
                    values.add_run(run.pos, run.pos + run.len, run.value);
                }
            }
            values.shrink();
            result.emplace_back(label, count, std::move(values));
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

    if (coord_to_header_) {
        if (annotator_->num_labels() != coord_to_header_->num_columns()) {
            throw std::runtime_error(fmt::format("Incompatible number of columns in CoordToHeader "
                                                 "and annotation matrix: {} != {}",
                                                 coord_to_header_->num_columns(),
                                                 annotator_->num_labels()));
        }

        // split the coordinates across sequences, a block of rows at a time
        RecordFold fold(*coord_to_header_, RecordFold::kCoords);
        call_row_tuples(*tuple_matrix, row_indices,
                        [&](size_t i, const MultiIntMatrix::RowTuples &row) {
            fold.add_row(row, kmer_positions[i]);
        });
        const auto counts = fold.select(min_count, num_top_labels);

        std::vector<std::tuple<std::string, size_t, PositionCoords, size_t>> result;
        result.reserve(counts.size());
        std::vector<std::pair<uint32_t, uint64_t>> coords_at; // (position, coordinate)
        for (const auto &[h, count] : counts) {
            const auto &[col, header] = h;
            coords_at.clear();
            for (const RecordFold::Run &run : fold.take_runs(h)) {
                for (uint32_t i = 0; i < run.len; ++i) {
                    coords_at.emplace_back(run.pos + i, run.value + i);
                }
            }
            std::sort(coords_at.begin(), coords_at.end());
            PositionCoords coords(num_kmers);
            for (const auto &[position, coord] : coords_at) {
                coords.add(position, coord);
            }
            coords.shrink();
            result.emplace_back(coord_to_header_->get_headers(col)[header],
                                count, std::move(coords),
                                coord_to_header_->num_kmers_in_sequence(col, header));
        }

        return result;
    }

    // the rows are read twice (to count, then to fill the values), so each distinct row is
    // decoded once and kept: a k-mer repeated m times would otherwise be decoded m times
    std::vector<row_index> distinct;
    std::vector<uint32_t> row_of;
    {
        tsl::hopscotch_map<row_index, uint32_t> seen;
        index_distinct_rows(row_indices, seen, distinct, row_of);
    }
    const auto rows_tuples = tuple_matrix->get_row_tuples(distinct);

    auto call_rows_tuples = [&](const auto &callback) {
        for (uint32_t i : row_of) {
            callback(rows_tuples[i]);
        }
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
    return encode_presence_mask(PresenceRuns(mask));
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

} // namespace graph
} // namespace mtg
