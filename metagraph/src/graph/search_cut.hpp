#ifndef __SEARCH_CUT_HPP__
#define __SEARCH_CUT_HPP__

#include <algorithm>
#include <cstddef>
#include <optional>
#include <utility>

namespace mtg::graph {

// Facts about the count-ranked cut, before per-position values are materialized.
// A tie is on count only. This does not certify a median-abundance ordering.
struct SearchCut {
    size_t matched = 0;
    size_t returned = 0;
    std::optional<size_t> nth_score;
    size_t tied = 0;
    size_t tied_returned = 0;
};

// Container `code_counts` must have pairs (index, count)
template <class Container>
void top_n_sorted(Container& code_counts, size_t num_top_labels, SearchCut *cut = nullptr) {
    if (cut) {
        *cut = {};
        cut->matched = code_counts.size();
        cut->returned = std::min(num_top_labels, code_counts.size());
    }
    if (code_counts.size() <= num_top_labels)
        return;
    auto comp = [](const auto &x, const auto &y) {
        return std::make_pair(y.second, x.first)
              < std::make_pair(x.second, y.first);
    };
    std::nth_element(code_counts.begin(),
                     code_counts.begin() + num_top_labels,
                     code_counts.end(),
                     comp);
    if (cut && num_top_labels) {
        auto end = code_counts.begin() + num_top_labels;
        const size_t nth = std::min_element(code_counts.begin(), end,
            [](const auto &a, const auto &b) { return a.second < b.second; })->second;
        cut->nth_score = nth;
        auto tied = [nth](const auto &value) { return value.second == nth; };
        cut->tied_returned = std::count_if(code_counts.begin(), end, tied);
        cut->tied = cut->tied_returned + std::count_if(end, code_counts.end(), tied);
    }
    // leave only the first |num_top_labels| top labels
    code_counts.resize(num_top_labels);
    // sort the top labels by counts
    std::sort(code_counts.begin(), code_counts.end(), comp);
}

} // namespace mtg::graph
#endif
