#ifndef __SEARCH_CUT_HPP__
#define __SEARCH_CUT_HPP__

#include <cstddef>
#include <cstdint>

namespace mtg::graph {

/**
 * What the top-N cut of one query sequence's search took (SPEC-multigraph-envelope.md section
 * 3): the labels above the request's thresholds before the cut (|matched|, master's member:
 * every candidate of the cut; a unit refused alone keeps it); the rows returned; whether every match above the request's thresholds was returned (no
 * cut, |complete|); otherwise the score of the N-th returned row, the labels with that score,
 * returned or not, and how many of them were returned. On an annotation with k-mer counts
 * (|abundances|) also the N-th returned row's median abundance and the labels with the N-th
 * score and that median. The median is read only where the N-th score is shared (|tied| > 1,
 * |nth_median_read|): a cut whose N-th score has one label orders no tie by median and states
 * none. A complete cut leaves the score, the counts and the median at 0, which a writer states
 * as null (the score and the median) or 0 (the counts); so does a cut with no row returned (N = 0),
 * which has no N-th row. A unit refused alone (|refused|) states no row: its cut keeps
 * |abundances| (the member set its entry states) and is written as refused.
 */
struct SearchCut {
    // the account that refused the unit (this query sequence on this pair) alone, with a
    // /search's field `partial_refusal` (SPEC-multigraph-envelope.md section 3, lane T7-P1)
    enum class Refused : uint8_t { NO, VALUES, ANSWER };

    size_t matched = 0;
    size_t returned = 0;
    size_t nth_score = 0;
    size_t tied = 0;
    size_t tied_returned = 0;
    uint64_t nth_median = 0;
    size_t tied_same_median = 0;
    // a refused unit: the bound (--search-max-answer-memory-mb) and, refused by the answer
    // account, its own object's price, both in MiB (rounded up)
    uint64_t refused_max_mb = 0;
    uint64_t refused_priced_mb = 0;
    // (the flags last, so that the cut takes the bytes it took before |refused|)
    bool complete = true;
    bool abundances = false;
    bool nth_median_read = false;
    Refused refused = Refused::NO;
};

} // namespace mtg::graph
#endif
