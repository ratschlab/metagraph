#ifndef __TEST_DBG_ALIGNER_HELPERS_HPP__
#define __TEST_DBG_ALIGNER_HELPERS_HPP__

#include <string>

#include "gtest/gtest.h"

#include "../test_helpers.hpp"

#include "graph/alignment/aligner_helper.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/alignment/dbg_aligner.hpp"


namespace mtg {
namespace test {

using namespace mtg::graph;
using namespace mtg::graph::align;


template <typename NodeType>
inline void check_json_dump_load(const DeBruijnGraph &graph,
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

DBGAligner<>::DBGQueryAlignment
inline get_extend(std::shared_ptr<const DeBruijnGraph> graph,
                  const DBGAlignerConfig &config,
                  const DBGAligner<>::DBGQueryAlignment &paths,
                  const std::string &query) {
    assert(graph.get());
    EXPECT_EQ(query, paths.get_query());
    auto uniconfig = config;
    uniconfig.max_seed_length = std::numeric_limits<size_t>::max();

    return std::dynamic_pointer_cast<const DBGSuccinct>(graph)
        ? DBGAligner<SuffixSeeder<>>(*graph, uniconfig).align(query)
        : DBGAligner<UniMEMSeeder<>>(*graph, uniconfig).align(query);
}

inline void check_extend(std::shared_ptr<const DeBruijnGraph> graph,
                         const DBGAlignerConfig &config,
                         const DBGAligner<>::DBGQueryAlignment &paths,
                         const std::string &query) {
    auto unimem_paths = get_extend(graph, config, paths, query);

    ASSERT_EQ(paths.size(), unimem_paths.size());

    for (size_t i = 0; i < paths.size(); ++i) {
        EXPECT_EQ(paths[i], unimem_paths[i]) << paths[i] << "\n" << unimem_paths[i];
    }
}

inline void check_chain(const DBGAligner<>::DBGQueryAlignment &paths,
                        const DeBruijnGraph &graph,
                        const DBGAlignerConfig &config,
                        size_t expected_size) {
    for (const auto &path : paths) {
        EXPECT_TRUE(path.is_valid(graph, &config)) << path;
        check_json_dump_load(graph, path,
                             paths.get_query(),
                             paths.get_query_reverse_complement());
    }

    if (!std::is_sorted(paths.begin(), paths.end(),
                        [](const auto &a, const auto &b) {
                            return a.get_query().end() < b.get_query().begin();
                        })
            || paths.size() != expected_size) {
        for (const auto &path : paths) {
            TEST_COUT << uintptr_t(path.get_query().data()) << " "
                      << uintptr_t(path.get_query().data() + path.get_query().size()) << " "
                      << path << "\n";
        }
        ASSERT_EQ(expected_size, paths.size());
    }
}

} // namespace test
} // namespace mtg

#endif // __TEST_DBG_ALIGNER_HELPERS_HPP__
