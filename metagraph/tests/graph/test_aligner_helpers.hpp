#ifndef __TEST_DBG_ALIGNER_HELPERS_HPP__
#define __TEST_DBG_ALIGNER_HELPERS_HPP__

#include <string>

#include "gtest/gtest.h"

#include "graph/alignment/alignment.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::align;
using namespace mtg::test;
using namespace mtg::kmer;


void check_json_dump_load(const DeBruijnGraph &graph,
                          const Alignment &alignment,
                          const std::string &query,
                          const std::string &rc_query = "") {
    ASSERT_TRUE(!rc_query.size() || query.size() == rc_query.size());

    const auto& path_query = alignment.get_orientation() && rc_query.size()
        ? rc_query
        : query;

    ASSERT_EQ(std::string(path_query.c_str() + alignment.get_clipping(),
                          alignment.get_query_view().size()),
              alignment.get_query_view());

    Alignment load_alignment;
    std::string load_sequence;
    load_alignment.load_from_json(alignment.to_json(graph.get_k()), graph, &load_sequence);

    EXPECT_EQ(path_query, load_sequence);

    EXPECT_EQ(alignment, load_alignment)
        << alignment.get_orientation() << " "
        << load_alignment.get_orientation() << "\n"
        << alignment.get_score() << " "
        << load_alignment.get_score() << "\n"
        << alignment.get_cigar().get_num_matches() << " "
        << load_alignment.get_cigar().get_num_matches() << "\n"
        << alignment.get_sequence() << " "
        << load_alignment.get_sequence() << "\n"
        << alignment.get_cigar().to_string() << " "
        << load_alignment.get_cigar().to_string() << "\n"
        << alignment.get_query_view() << " "
        << load_alignment.get_query_view() << "\n";
}

AlignmentResults get_extend(std::shared_ptr<const DeBruijnGraph> graph,
                            const DBGAlignerConfig &config,
                            const AlignmentResults &paths,
                            const std::string &query) {
    assert(graph.get());
    EXPECT_EQ(query, paths.get_query());
    auto uniconfig = config;
    uniconfig.max_seed_length = std::numeric_limits<size_t>::max();
    return DBGAligner<>(*graph, uniconfig).align(query);
}

inline void check_extend(std::shared_ptr<const DeBruijnGraph> graph,
                         const DBGAlignerConfig &config,
                         const AlignmentResults &paths,
                         const std::string &query) {
    auto unimem_paths = get_extend(graph, config, paths, query);

    ASSERT_EQ(paths.size(), unimem_paths.size());

    for (size_t i = 0; i < paths.size(); ++i) {
        if (graph->get_mode() == DeBruijnGraph::CANONICAL) {
            auto alt = unimem_paths[i];
            alt.reverse_complement(*graph, paths.get_query(!alt.get_orientation()));
            if (alt.size()) {
                EXPECT_TRUE(paths[i] == unimem_paths[i] || paths[i] == alt)
                     << i << "\n" << paths[i] << "\n" << unimem_paths[i];
            } else {
                EXPECT_EQ(paths[i], unimem_paths[i]) << i << "\n" << paths[i] << "\n" << unimem_paths[i];
            }
        } else {
            EXPECT_EQ(paths[i], unimem_paths[i]) << i << "\n" << paths[i] << "\n" << unimem_paths[i];
        }
    }
}

} // namespace

#endif // __TEST_DBG_ALIGNER_HELPERS_HPP__
