#include <gtest/gtest.h>

#include "all/test_dbg_helpers.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/annotation_key_mapper.hpp"
#include "graph/representation/canonical_dbg.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::test;

const std::string kSequence = "ATGCGATCGATTAGCCATGCGATTCAGGCAT";

TEST(AnnotationKeyMapper, BasicKeysAreNodes) {
    auto graph = build_graph_batch<DBGSuccinct>(5, { kSequence }, DeBruijnGraph::BASIC);
    AnnotationKeyMapper mapper(*graph);
    ASSERT_FALSE(mapper.needs_sequence());

    auto path = map_to_nodes_sequentially(*graph, kSequence);
    EXPECT_EQ(path, mapper.keys(path));
}

// The test helper returns primary graphs wrapped in CanonicalDBG
const CanonicalDBG& wrapped_primary(const std::shared_ptr<DeBruijnGraph> &graph) {
    const auto *canonical = dynamic_cast<const CanonicalDBG*>(graph.get());
    EXPECT_TRUE(canonical);
    return *canonical;
}

TEST(AnnotationKeyMapper, PrimaryGraphMustBeWrapped) {
    auto graph = build_graph_batch<DBGSuccinct>(5, { kSequence }, DeBruijnGraph::PRIMARY);
    const DeBruijnGraph &primary = wrapped_primary(graph).get_graph();
    ASSERT_EQ(DeBruijnGraph::PRIMARY, primary.get_mode());
    EXPECT_THROW(AnnotationKeyMapper mapper(primary), std::invalid_argument);
}

#if ! _PROTEIN_GRAPH
// A k-mer and its reverse complement carry the same annotation row, so their
// keys must coincide, in reverse order along the path.
void check_reverse_complement_keys(const DeBruijnGraph &graph) {
    AnnotationKeyMapper mapper(graph);

    std::string rc = kSequence;
    reverse_complement(rc.begin(), rc.end());

    auto path = map_to_nodes_sequentially(graph, kSequence);
    auto rc_path = map_to_nodes_sequentially(graph, rc);
    ASSERT_EQ(path.size(), rc_path.size());

    auto keys = mapper.keys(path, mapper.needs_sequence() ? kSequence : "");
    auto rc_keys = mapper.keys(rc_path, mapper.needs_sequence() ? rc : "");
    std::reverse(rc_keys.begin(), rc_keys.end());
    EXPECT_EQ(keys, rc_keys);
    EXPECT_EQ(keys.end(), std::find(keys.begin(), keys.end(), DeBruijnGraph::npos));

    for (size_t i = 0; i < path.size(); ++i) {
        EXPECT_EQ(keys[i], mapper.key(path[i], std::string_view(kSequence).substr(i, graph.get_k())));
    }
}

TEST(AnnotationKeyMapper, CanonicalGraphMapsBothStrandsToOneKey) {
    auto graph = build_graph_batch<DBGSuccinct>(5, { kSequence }, DeBruijnGraph::CANONICAL);
    ASSERT_TRUE(AnnotationKeyMapper(*graph).needs_sequence());
    check_reverse_complement_keys(*graph);
}

TEST(AnnotationKeyMapper, WrappedPrimaryGraphMapsBothStrandsToOneKey) {
    auto primary = build_graph_batch<DBGSuccinct>(5, { kSequence }, DeBruijnGraph::PRIMARY);
    const CanonicalDBG &graph = wrapped_primary(primary);
    ASSERT_FALSE(AnnotationKeyMapper(graph).needs_sequence());
    check_reverse_complement_keys(graph);
}
#endif

} // namespace
