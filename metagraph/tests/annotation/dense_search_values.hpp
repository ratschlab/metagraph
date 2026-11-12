#ifndef __DENSE_SEARCH_VALUES_HPP__
#define __DENSE_SEARCH_VALUES_HPP__

#include <string>
#include <tuple>
#include <vector>

#include <sdsl/int_vector.hpp>

#include "common/vector.hpp"
#include "graph/annotated_dbg.hpp"


namespace mtg {
namespace test {

// The rows of AnnotatedDBG's per-position queries (get_top_label_signatures, get_kmer_counts,
// get_kmer_coordinates) with each row's values as the vector of the query's k-mers they stand
// for (PresenceRuns::to_dense and the others), to compare with expected vectors

inline std::vector<std::tuple<std::string, size_t, sdsl::bit_vector>>
dense(const std::vector<std::tuple<std::string, size_t, graph::PresenceRuns>> &rows) {
    std::vector<std::tuple<std::string, size_t, sdsl::bit_vector>> result;
    for (const auto &[label, count, values] : rows) {
        result.emplace_back(label, count, values.to_dense());
    }
    return result;
}

inline std::vector<std::tuple<std::string, size_t, std::vector<size_t>>>
dense(const std::vector<std::tuple<std::string, size_t, graph::CountRuns>> &rows) {
    std::vector<std::tuple<std::string, size_t, std::vector<size_t>>> result;
    for (const auto &[label, count, values] : rows) {
        result.emplace_back(label, count, values.to_dense());
    }
    return result;
}

inline std::vector<std::tuple<std::string, size_t, std::vector<SmallVector<uint64_t>>, size_t>>
dense(const std::vector<std::tuple<std::string, size_t, graph::PositionCoords, size_t>> &rows) {
    std::vector<std::tuple<std::string, size_t, std::vector<SmallVector<uint64_t>>, size_t>> result;
    for (const auto &[label, count, values, length] : rows) {
        result.emplace_back(label, count, values.to_dense(), length);
    }
    return result;
}

} // namespace test
} // namespace mtg

#endif // __DENSE_SEARCH_VALUES_HPP__
