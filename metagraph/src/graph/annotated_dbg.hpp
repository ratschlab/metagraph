#ifndef __ANNOTATED_DBG_HPP__
#define __ANNOTATED_DBG_HPP__

#include <cassert>
#include <functional>
#include <memory>
#include <mutex>
#include <string_view>

#include <sdsl/int_vector.hpp>

#include "representation/base/sequence_graph.hpp"
#include "annotation/representation/base/annotation.hpp"
#include "common/vector.hpp"
#include "annotation/coord_to_header.hpp"


namespace mtg {
namespace graph {

/**
 * What the top-N cut of one query sequence's search took (SPEC-multigraph-envelope.md section
 * 3): the rows returned; whether every match above the request's thresholds was returned (no
 * cut, |complete|); otherwise the score of the N-th returned row, the labels with that score,
 * returned or not, and how many of them were returned. On an annotation with k-mer counts
 * (|abundances|) also the N-th returned row's median abundance and the labels with the N-th
 * score and that median. The median is read only where the N-th score is shared (|tied| > 1,
 * |nth_median_read|): a cut whose N-th score has one label orders no tie by median and states
 * none. A complete cut leaves the members after |complete| at 0, which a writer states as null
 * (the score and the median) or 0 (the counts); so does a cut with no row returned (N = 0),
 * which has no N-th row.
 */
struct SearchCut {
    size_t returned = 0;
    bool complete = true;
    size_t nth_score = 0;
    size_t tied = 0;
    size_t tied_returned = 0;
    bool abundances = false;
    uint64_t nth_median = 0;
    size_t tied_same_median = 0;
    bool nth_median_read = false;
};

/**
 * The lower median of a label's abundances over the query's k-mer positions where it is
 * present, one value per position: the zeros of |values| excluded, the others sorted
 * ascending, the element at index floor((n - 1) / 2); 0 when no value is left.
 */
uint64_t lower_median(std::vector<uint64_t> values);

/**
 * lower_median of the abundances added, held as their distinct values with their counts
 * instead of one value per position: a run of one value along the query is one entry, and
 * the entries are merged by value whenever they have doubled since the last merge, so a
 * label's counter holds O(its distinct values), not O(its positions).
 */
class LowerMedianCounter {
  public:
    // |positions| query positions with abundance |value|; a zero (the label absent) is no
    // abundance and not counted
    void add(uint64_t value, uint64_t positions = 1);
    // the lower median of the values added, 0 if none
    uint64_t lower_median();

  private:
    void merge();

    std::vector<std::pair<uint64_t, uint64_t>> counts_;
    size_t merged_ = 0;
};

/**
 * The lower medians of the abundances of |columns| (in that order) over a query's k-mer
 * positions, read from |matrix|: |row_counts| are the rows the query's k-mers map to, each
 * with the number of query positions it stands for. A matrix other than a CSRMatrix is read
 * |rows_per_read| rows at a time; each column holds its distinct values (LowerMedianCounter).
 */
std::vector<uint64_t>
abundance_medians(const annot::matrix::IntMatrix &matrix,
                  const std::vector<std::pair<uint64_t, size_t>> &row_counts,
                  const std::vector<uint32_t> &columns,
                  size_t num_columns,
                  size_t rows_per_read);

// One row of a search in the C5 order (search_row_before), |median| 0 on an annotation
// without k-mer counts
struct SearchRow {
    size_t score;
    uint64_t median;
    std::string_view label;
};

// The C5 order of a search's rows (SPEC-multigraph-envelope.md section 3): score descending;
// then the median abundance descending (an annotation with k-mer counts; every median is 0 on
// another annotation); then the label ascending, by the bytes of the decoded label
bool search_row_before(const SearchRow &a, const SearchRow &b);

/**
 * The top |num_top_labels| of a search's candidate rows in the C5 order, which every
 * engine-local cut of a search takes (a row cut away here cannot be recovered by a client).
 * |scores|: the candidates' scores, every candidate above the request's thresholds. |label|:
 * a candidate's label, decoded only for the rows returned and for the rows of the N-th score
 * that compete for the last places (the tie group at the cut, of the N-th row's median where
 * medians are given), never for every candidate. |medians|: on an annotation with k-mer counts
 * the lower medians (lower_median) of the candidates it is passed, asked only for the rows a
 * median is compared for (the tie group at a cut, and the rows above it that share a score);
 * empty otherwise. Candidates with one label (a header in two
 * columns of one pair) keep their order in |scores|.
 * Returns the returned candidates (their indexes in |scores|) with their labels in the C5
 * order, and states the cut in |cut|.
 */
std::vector<std::pair<size_t, std::string>>
search_top_n(const std::vector<size_t> &scores,
             size_t num_top_labels,
             const std::function<std::string(size_t)> &label,
             const std::function<std::vector<uint64_t>(const std::vector<size_t> &)> &medians,
             SearchCut *cut);


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
    // Whether |annotation| has one row per node of |graph| (of the graph a CanonicalDBG wraps):
    // check_compatibility of the pair, for a loader to refuse a pair before it builds the index
    // of it, whose constructor asserts it
    static bool compatible(const SequenceGraph &graph, const Annotator &annotation);

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
    // |coord_to_header|: the mapping loaded with the annotation (a unique_ptr converts), or
    // the one another index of the same annotation already holds (share_coord_to_header)
    AnnotatedDBG(std::shared_ptr<DeBruijnGraph> dbg,
                 std::unique_ptr<Annotator>&& annotation,
                 bool force_fast = false,
                 std::shared_ptr<const annot::CoordToHeader> coord_to_header = {});

    using AnnotatedSequenceGraph::get_labels;

    const DeBruijnGraph& get_graph() const { return dbg_; }

    // Returns pointer to the coordinate-to-header mapping if available, nullptr otherwise.
    // This mapping transforms file-based coordinates to sequence-header-based coordinates
    // for query results. The returned pointer is valid for the lifetime of this object.
    const annot::CoordToHeader* get_coord_to_header() const { return coord_to_header_.get(); }

    // The mapping to build another index of the same annotation with (a server's copy of a
    // pair loaded into RAM for one request): the mapping depends on the .seqs alone, so the
    // copy shares the object and its reverse header index instead of loading and building
    // identical ones. Null without a mapping.
    std::shared_ptr<const annot::CoordToHeader> share_coord_to_header() const {
        return coord_to_header_;
    }

    // Whether |annotation| holds k-mer counts (abundances): its matrix is an integer matrix
    // and not a coordinate matrix (row_diff_int_brwt, int_row_diff_disk, the integer CSC
    // forms). The coordinates of a coordinate matrix are positions, not abundances, so a
    // search's ties there are ordered by label alone (SPEC-multigraph-envelope.md section 3)
    static bool holds_kmer_abundances(const Annotator &annotation);

    // Whether this index's annotation holds k-mer counts, which order the ties of a search's
    // top N by median abundance (search_top_n): holds_kmer_abundances of the annotation, or
    // what set_kmer_abundances stated
    bool has_kmer_abundances() const { return kmer_abundances_; }

    // A query graph's annotation is sliced from a full one into an integer matrix whenever a
    // query reads integers, the coordinate counts of a coordinate annotation included: it
    // states the full annotation's kind here, so that a batched search orders its ties as the
    // full index does
    void set_kmer_abundances(bool abundances) { kmer_abundances_ = abundances; }

    // A query graph sliced from an annotation with k-mer counts into a binary matrix (the
    // query modes that read no counts keep the binary slice, which costs a fraction of the
    // integer one) reads the abundances of a search's ties from the full index |source|: the
    // query sequence mapped onto its graph, the values read from its integer matrix, only
    // where a tie needs a median (search_top_n). |source| must outlive this index.
    void set_abundance_source(const AnnotatedDBG *source) { abundance_source_ = source; }

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

    // The top-N queries below return their rows in the C5 order (search_top_n): score
    // descending, then the median abundance descending on an annotation with k-mer counts,
    // then the label ascending, the top N taken in that order; |cut|, if given, states the
    // cut (SearchCut).

    // Return top |num_top_labels| labels with their counts.
    // The returned counts are weighted by the annotated relation counts if
    // |with_kmer_counts| is true.
    std::vector<std::pair<Label, size_t>>
    get_top_labels(std::string_view sequence,
                   size_t num_top_labels,
                   double discovery_fraction = 0.0,
                   double presence_fraction = 0.0,
                   bool with_kmer_counts = false,
                   SearchCut *cut = nullptr) const;

    // returns tuples (label, num_kmer_matches, kmer_abundances)
    std::vector<std::tuple<Label, size_t, std::vector<size_t>>>
    get_kmer_counts(std::string_view sequence,
                    size_t num_top_labels,
                    double discovery_fraction,
                    double presence_fraction,
                    SearchCut *cut = nullptr) const;

    // returns tuples (label, num_kmer_matches, kmer_coordinates, num_kmers_in_target).
    // num_kmers_in_target is 0 when no CoordToHeader mapping is loaded (the
    // label refers to a whole annotation column rather than a single indexed
    // sequence).
    std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>, size_t>>
    get_kmer_coordinates(std::string_view sequence,
                         size_t num_top_labels,
                         double discovery_fraction,
                         double presence_fraction,
                         SearchCut *cut = nullptr) const;

    std::vector<std::tuple<Label, size_t, sdsl::bit_vector>>
    get_top_label_signatures(std::string_view sequence,
                             size_t num_top_labels,
                             double discovery_fraction = 0.0,
                             double presence_fraction = 0.0,
                             SearchCut *cut = nullptr) const;

  private:
    // get_kmer_coordinates, with the records of a CoordToHeader mapping ranked by the number
    // of their coordinates the query's k-mers match if |score_by_coordinates| (the score
    // get_top_labels reports with |with_kmer_counts|), by the number of k-mers otherwise
    std::vector<std::tuple<Label, size_t, std::vector<SmallVector<uint64_t>>, size_t>>
    query_kmer_coordinates(std::string_view sequence,
                           size_t num_top_labels,
                           double discovery_fraction,
                           double presence_fraction,
                           SearchCut *cut,
                           bool score_by_coordinates) const;

    // The integer matrix to read the abundances of a search's ties from: the annotation's
    // matrix when it holds k-mer counts (has_kmer_abundances), null otherwise
    const annot::matrix::IntMatrix* abundance_matrix() const;

    // The lower medians of the abundances of the columns it is passed over the query
    // |sequence|'s k-mer positions where they are present, which search_top_n asks for the
    // ties it orders; empty where this index holds no abundances. |row_counts|: the rows of
    // this index's graph the query's k-mers map to, each with the number of query positions
    // it stands for. Read from this index's integer matrix, or from the abundance source's
    // (set_abundance_source), the query mapped onto the source's graph. |sequence| and
    // |row_counts| must outlive the function returned.
    std::function<std::vector<uint64_t>(const std::vector<uint32_t> &)>
    tie_medians(std::string_view sequence,
                const std::vector<std::pair<row_index, size_t>> &row_counts) const;
    // tie_medians over |row_indices|, the rows of the query's k-mers found, one per position:
    // collapsed into rows with their multiplicity only when a tie asks
    std::function<std::vector<uint64_t>(const std::vector<uint32_t> &)>
    tie_medians(std::string_view sequence, const std::vector<row_index> &row_indices) const;
    // tie_medians read from the abundance source (set_abundance_source), empty without one
    std::function<std::vector<uint64_t>(const std::vector<uint32_t> &)>
    source_tie_medians(std::string_view sequence) const;

    DeBruijnGraph &dbg_;
    std::shared_ptr<const annot::CoordToHeader> coord_to_header_;
    bool kmer_abundances_;
    const AnnotatedDBG *abundance_source_ = nullptr;
};

} // namespace graph
} // namespace mtg

#endif // __ANNOTATED_DBG_HPP__
