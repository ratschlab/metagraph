#ifndef __ANNOTATED_DBG_HPP__
#define __ANNOTATED_DBG_HPP__

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <mutex>
#include <vector>

#include <sdsl/int_vector.hpp>

#include "representation/base/sequence_graph.hpp"
#include "annotation/representation/base/annotation.hpp"
#include "common/vector.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "annotation/coord_to_header.hpp"


namespace mtg {
namespace graph {

/**
 * The per-position values of one label returned by get_top_label_signatures,
 * get_kmer_counts and get_kmer_coordinates. Each stands for a vector with one
 * value per k-mer of the query (to_dense() returns it) and is built by adding
 * positions in ascending order.
 *
 * A value holds only the positions its label matches, as long as that takes fewer
 * bytes than the vector, and the vector itself otherwise: its buffers never grow
 * beyond dense_bytes(size()), so it never holds more than the vector would. Runs
 * store 32-bit positions; a query of 2^32 or more k-mers is held as the vector from
 * the start.
 */

// k-mer presence mask: the maximal runs [begin, end) of matched positions, or the mask
class PresenceRuns {
  public:
    struct Run {
        uint32_t begin;
        uint32_t end;
    };

    PresenceRuns() = default;
    explicit PresenceRuns(size_t size);
    explicit PresenceRuns(const sdsl::bit_vector &mask);

    size_t size() const { return size_; }

    // mark |position|, at or after the last position marked
    void add(size_t position) { add_run(position, position + 1); }
    // mark the positions [begin, end), at or after the last position marked
    void add_run(size_t begin, size_t end) {
        assert(begin <= end && end <= size_);
        if (begin == end)
            return;
        if (dense_form_) {
            mark(begin, end);
            return;
        }
        assert(runs_.empty() || begin >= runs_.back().end);
        if (!runs_.empty() && runs_.back().end == begin) {
            runs_.back().end = end;
        } else if (runs_.size() < runs_.capacity() || grow()) {
            runs_.push_back({ static_cast<uint32_t>(begin), static_cast<uint32_t>(end) });
        } else {
            to_dense_form();
            mark(begin, end);
        }
    }

    // call |callback(begin, end)| for the maximal runs of marked positions, in order
    template <class Callback>
    void call_runs(const Callback &callback) const;

    // number of marked positions
    size_t count() const;
    sdsl::bit_vector to_dense() const;
    bool is_dense() const { return dense_form_; }
    // bytes of the buffers held, at most dense_bytes(size())
    size_t heap_bytes() const;
    static size_t dense_bytes(size_t size) { return (size + 63) / 64 * sizeof(uint64_t); }
    // release spare capacity once all positions are added
    void shrink() { runs_.shrink_to_fit(); }

    // same size and same positions marked, in either form
    bool operator==(const PresenceRuns &other) const;

  private:
    size_t size_ = 0;
    std::vector<Run> runs_;
    sdsl::bit_vector dense_;
    bool dense_form_ = false;

    // room for one more run within dense_bytes(size_), false if there is none
    bool grow();
    void to_dense_form();
    void mark(size_t begin, size_t end);
};

// k-mer counts: the maximal runs [begin, end) of one nonzero value, or the counts
class CountRuns {
  public:
    struct Run {
        uint32_t begin;
        uint32_t end;
        uint64_t value;
    };

    CountRuns() = default;
    explicit CountRuns(size_t size);
    explicit CountRuns(const std::vector<size_t> &values);

    size_t size() const { return size_; }

    // set |position| to |value|, at or after the last position set (zero is no value)
    void add(size_t position, size_t value) { add_run(position, position + 1, value); }
    // set the positions [begin, end) to |value|, at or after the last position set
    void add_run(size_t begin, size_t end, size_t value) {
        assert(begin <= end && end <= size_);
        if (begin == end || !value)
            return;
        if (dense_form_) {
            std::fill(dense_.begin() + begin, dense_.begin() + end, value);
            return;
        }
        assert(runs_.empty() || begin >= runs_.back().end);
        if (!runs_.empty() && runs_.back().end == begin && runs_.back().value == value) {
            runs_.back().end = end;
        } else if (runs_.size() < runs_.capacity() || grow()) {
            runs_.push_back({ static_cast<uint32_t>(begin), static_cast<uint32_t>(end), value });
        } else {
            to_dense_form();
            std::fill(dense_.begin() + begin, dense_.begin() + end, value);
        }
    }

    // call |callback(begin, end, value)| for the maximal runs of one nonzero value, in order
    template <class Callback>
    void call_runs(const Callback &callback) const;

    // number of positions with a nonzero value
    size_t count() const;
    std::vector<size_t> to_dense() const;
    bool is_dense() const { return dense_form_; }
    // bytes of the buffers held, at most dense_bytes(size())
    size_t heap_bytes() const;
    static size_t dense_bytes(size_t size) { return size * sizeof(size_t); }
    // release spare capacity once all positions are added
    void shrink() { runs_.shrink_to_fit(); }

    // same size and same values, in either form
    bool operator==(const CountRuns &other) const;

  private:
    size_t size_ = 0;
    std::vector<Run> runs_;
    std::vector<size_t> dense_;
    bool dense_form_ = false;

    bool grow();
    void to_dense_form();
};

// k-mer coordinates: the maximal runs [begin, end) of positions with coordinates and
// the coordinates of each of these positions, in order, or a coordinate list per position
class PositionCoords {
  public:
    using Coords = SmallVector<uint64_t>;
    struct Run {
        uint32_t begin;
        uint32_t end;
    };

    PositionCoords() = default;
    explicit PositionCoords(size_t size);
    explicit PositionCoords(const std::vector<Coords> &values);

    size_t size() const { return size_; }

    // append |coord| to the coordinates of |position|, at or after the last position
    // with coordinates
    void add(size_t position, uint64_t coord) {
        assert(position < size_);
        if (dense_form_) {
            dense_[position].push_back(coord);
            return;
        }
        assert(runs_.empty() || position + 1 >= runs_.back().end);
        if (!runs_.empty() && runs_.back().end == position + 1) {
            coords_.back().push_back(coord);
        } else if (open(position)) {
            coords_.emplace_back();
            coords_.back().push_back(coord);
        } else {
            dense_[position].push_back(coord);
        }
    }
    // set the coordinates of |position|, after the last position with coordinates
    // (no coordinates is no value)
    void add(size_t position, Coords&& position_coords) {
        assert(position < size_);
        if (position_coords.empty())
            return;
        if (dense_form_) {
            dense_[position] = std::move(position_coords);
            return;
        }
        assert(runs_.empty() || position >= runs_.back().end);
        if (open(position)) {
            coords_.push_back(std::move(position_coords));
        } else {
            dense_[position] = std::move(position_coords);
        }
    }

    // call |callback(position, coords)| for the positions with coordinates, in order
    template <class Callback>
    void call_positions(const Callback &callback) const;

    // number of positions with coordinates
    size_t count() const;
    std::vector<Coords> to_dense() const;
    bool is_dense() const { return dense_form_; }
    // bytes of the buffers held, not counting the coordinate lists' own (the same in
    // both forms), at most dense_bytes(size())
    size_t heap_bytes() const;
    static size_t dense_bytes(size_t size) { return size * sizeof(Coords); }
    // release spare capacity once all positions are added
    void shrink() {
        runs_.shrink_to_fit();
        coords_.shrink_to_fit();
    }

    // same size and same coordinates, in either form
    bool operator==(const PositionCoords &other) const;

  private:
    size_t size_ = 0;
    std::vector<Run> runs_;
    std::vector<Coords> coords_;
    std::vector<Coords> dense_;
    bool dense_form_ = false;

    // start |position| (after the last position with coordinates) with room for its
    // coordinates; if the runs and coordinates would exceed dense_bytes(size_), switch
    // to a coordinate list per position and return false
    bool open(size_t position);
    void to_dense_form();
};

template <class Callback>
void PresenceRuns::call_runs(const Callback &callback) const {
    if (!dense_form_) {
        for (const Run &run : runs_) {
            callback(run.begin, run.end);
        }
        return;
    }
    // the ones of the mask, joined into maximal runs
    size_t begin = 0;
    size_t end = 0;
    call_ones(dense_, [&](size_t i) {
        if (i != end) {
            if (begin < end)
                callback(begin, end);
            begin = i;
        }
        end = i + 1;
    });
    if (begin < end)
        callback(begin, end);
}

template <class Callback>
void CountRuns::call_runs(const Callback &callback) const {
    if (!dense_form_) {
        for (const Run &run : runs_) {
            callback(run.begin, run.end, run.value);
        }
        return;
    }
    for (size_t i = 0; i < size_; ) {
        const size_t value = dense_[i];
        size_t end = i + 1;
        while (end < size_ && dense_[end] == value) {
            ++end;
        }
        if (value)
            callback(i, end, value);
        i = end;
    }
}

template <class Callback>
void PositionCoords::call_positions(const Callback &callback) const {
    if (!dense_form_) {
        size_t c = 0;
        for (const Run &run : runs_) {
            for (size_t i = run.begin; i < run.end; ++i) {
                callback(i, coords_[c++]);
            }
        }
        assert(c == coords_.size());
        return;
    }
    for (size_t i = 0; i < size_; ++i) {
        if (!dense_[i].empty())
            callback(i, dense_[i]);
    }
}


class AnnotatedSequenceGraph {
  public:
    typedef std::string Label;
    typedef annot::MultiLabelAnnotation<Label> Annotator;
    using node_index = SequenceGraph::node_index;
    using row_index = Annotator::Index;

    AnnotatedSequenceGraph(std::shared_ptr<SequenceGraph> graph,
                           std::unique_ptr<Annotator>&& annotation,
                           bool force_fast = false);

    virtual ~AnnotatedSequenceGraph() {}

    virtual std::vector<Label> get_labels(node_index index) const;

    // thread-safe, can be called from multiple threads concurrently
    virtual void annotate_sequence(std::string_view sequence,
                                   const std::vector<Label> &labels);
    // thread-safe, can be called from multiple threads concurrently
    virtual void annotate_sequences(
        const std::vector<std::pair<std::string, std::vector<Label>>> &data);

    virtual bool label_exists(const Label &label) const;

    virtual bool check_compatibility() const;

    virtual const SequenceGraph& get_graph() const { return *graph_; }
    std::shared_ptr<const SequenceGraph> get_graph_ptr() const { return graph_; }

    virtual const Annotator& get_annotator() const { return *annotator_; }

    static row_index graph_to_anno_index(node_index kmer_index) {
        assert(kmer_index);
        return kmer_index - 1;
    }
    static node_index anno_to_graph_index(row_index anno_index) {
        return anno_index + 1;
    }

  protected:
    std::shared_ptr<SequenceGraph> graph_;
    std::unique_ptr<Annotator> annotator_;

    std::mutex mutex_;
    bool force_fast_;
};


class AnnotatedDBG : public AnnotatedSequenceGraph {
  public:
    AnnotatedDBG(std::shared_ptr<DeBruijnGraph> dbg,
                 std::unique_ptr<Annotator>&& annotation,
                 bool force_fast = false,
                 std::unique_ptr<annot::CoordToHeader> coord_to_header = {});

    using AnnotatedSequenceGraph::get_labels;

    const DeBruijnGraph& get_graph() const { return dbg_; }

    // Returns pointer to the coordinate-to-header mapping if available, nullptr otherwise.
    // This mapping transforms file-based coordinates to sequence-header-based coordinates
    // for query results. The returned pointer is valid for the lifetime of this object.
    const annot::CoordToHeader* get_coord_to_header() const { return coord_to_header_.get(); }

    // add k-mer counts to the annotation, thread-safe for concurrent calls
    void add_kmer_counts(std::string_view sequence,
                         const std::vector<Label> &labels,
                         std::vector<uint64_t>&& kmer_counts);

    // add k-mer coordinates to the annotation, the binary annotation must exist
    void add_kmer_coord(std::string_view sequence,
                        const std::vector<Label> &labels,
                        uint64_t start);

    // add k-mer coordinates to the annotation, the binary annotation must exist
    void add_kmer_coords(
        const std::vector<std::tuple<std::string, std::vector<Label>, uint64_t>> &data);

    // annotate k-mer and their coordinates (combines annotate_sequences and add_kmer_coords)
    void annotate_kmer_coords(
        const std::vector<std::tuple<std::string, std::vector<Label>, uint64_t>> &data);

    /*********************** Special queries **********************/

    // return labels that occur at least in |discovery_fraction| k-mers
    // but skip the sequence if fewer than |presence_fraction| k-mers are matched against the graph
    std::vector<Label> get_labels(std::string_view sequence,
                                  double discovery_fraction = 0.0,
                                  double presence_fraction = 0.0) const;

    // Return top |num_top_labels| labels with their counts.
    // The returned counts are weighted by the annotated relation counts if
    // |with_kmer_counts| is true.
    std::vector<std::pair<Label, size_t>>
    get_top_labels(std::string_view sequence,
                   size_t num_top_labels,
                   double discovery_fraction = 0.0,
                   double presence_fraction = 0.0,
                   bool with_kmer_counts = false) const;

    // The per-position queries below hold each label's values only at the positions
    // it matches (PresenceRuns, CountRuns, PositionCoords), so their memory grows with
    // the matches, not with the labels returned times the k-mers of the query.

    // returns tuples (label, num_kmer_matches, kmer_abundances)
    std::vector<std::tuple<Label, size_t, CountRuns>>
    get_kmer_counts(std::string_view sequence,
                    size_t num_top_labels,
                    double discovery_fraction,
                    double presence_fraction) const;

    // returns tuples (label, num_kmer_matches, kmer_coordinates, num_kmers_in_target).
    // num_kmers_in_target is 0 when no CoordToHeader mapping is loaded (the
    // label refers to a whole annotation column rather than a single indexed
    // sequence).
    std::vector<std::tuple<Label, size_t, PositionCoords, size_t>>
    get_kmer_coordinates(std::string_view sequence,
                         size_t num_top_labels,
                         double discovery_fraction,
                         double presence_fraction) const;

    // With |num_passed|, also report how many labels passed the thresholds, before the cut
    // to |num_top_labels|. With |header_labels| false, rank annotation columns even when a
    // CoordToHeader mapping is loaded.
    std::vector<std::tuple<Label, size_t, PresenceRuns>>
    get_top_label_signatures(std::string_view sequence,
                             size_t num_top_labels,
                             double discovery_fraction = 0.0,
                             double presence_fraction = 0.0,
                             size_t *num_passed = nullptr,
                             bool header_labels = true) const;

  private:
    template <class Values>
    std::vector<std::tuple<Label, size_t, Values>>
    get_header_results(std::string_view sequence, size_t num_top_labels,
                       double discovery_fraction, double presence_fraction,
                       bool with_kmer_counts = false, size_t *num_passed = nullptr) const;

    DeBruijnGraph &dbg_;
    std::unique_ptr<annot::CoordToHeader> coord_to_header_;
};

// Encode a k-mer presence bitmask as alternating x<len>o<len>... runs (x = ones,
// o = zeros), trailing zeros included. Example: 11100110 -> x3o2x2o1
std::string encode_presence_mask(const sdsl::bit_vector &mask);
// The same string, from the runs of the mask
std::string encode_presence_mask(const PresenceRuns &mask);

} // namespace graph
} // namespace mtg

#endif // __ANNOTATED_DBG_HPP__
