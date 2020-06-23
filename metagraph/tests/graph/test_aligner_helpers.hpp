#ifndef __TEST_DBG_ALIGNER_HELPERS_HPP__
#define __TEST_DBG_ALIGNER_HELPERS_HPP__

#include <string>

#include "gtest/gtest.h"

#include "graph/alignment/aligner_helper.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace test {

template <typename NodeType>
void check_json_dump_load(const DeBruijnGraph &graph,
                          const Alignment<NodeType> &alignment,
                          const std::string &query,
                          const std::string &rc_query = "") {
    ASSERT_TRUE(!rc_query.size() || query.size() == rc_query.size());

    const auto& path_query = alignment.get_orientation() && rc_query.size()
        ? rc_query
        : query;

    ASSERT_EQ(std::string(path_query.c_str() + alignment.get_clipping(),
                          alignment.get_query().size()),
              alignment.get_query());

    Alignment<NodeType> load_alignment;
    auto load_sequence = load_alignment.load_from_json(
        alignment.to_json(path_query, graph),
        graph
    );

    EXPECT_EQ(path_query, *load_sequence);

    EXPECT_EQ(alignment, load_alignment)
        << alignment.get_orientation() << " "
        << load_alignment.get_orientation() << "\n"
        << alignment.get_score() << " "
        << load_alignment.get_score() << "\n"
        << alignment.get_num_matches() << " "
        << load_alignment.get_num_matches() << "\n"
        << alignment.get_sequence() << " "
        << load_alignment.get_sequence() << "\n"
        << alignment.get_cigar().to_string() << " "
        << load_alignment.get_cigar().to_string() << "\n"
        << alignment.get_query() << " "
        << load_alignment.get_query() << "\n";
}

} // namespace test
} // namespace mtg

#endif // __TEST_DBG_ALIGNER_HELPERS_HPP__
