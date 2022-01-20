#ifndef __ALIGNER_CHAINER_HPP__
#define __ALIGNER_CHAINER_HPP__

#include "aligner_alignment.hpp"

namespace mtg {
namespace graph {
namespace align {

typedef std::vector<Alignment> Chain;
typedef Alignment::score_t score_t;
typedef std::tuple<Alignment::Column /* label */,
                   ssize_t /* coordinate */,
                   size_t /* seed clipping */,
                   Alignment::node_index /* first node of seed */,
                   ssize_t /* seed length */,
                   score_t /* chain score */,
                   uint32_t /* previous seed index */,
                   uint32_t /* current seed index */> TableElem;
typedef std::vector<TableElem> ChainDPTable;

// Given sets of forward and reverse-complement seeds, construct and call chains of seeds
// until terminate() returns true.
std::pair<size_t, size_t>
call_seed_chains_both_strands(std::string_view forward,
                              std::string_view reverse,
                              const DeBruijnGraph &graph,
                              const DBGAlignerConfig &config,
                              std::vector<Alignment>&& fwd_seeds,
                              std::vector<Alignment>&& bwd_seeds,
                              const std::function<void(Chain&&, score_t)> &callback,
                              const std::function<bool()> &terminate = []() { return false; });

// Given a set of local alignments, use sparse dynamic programming to construct
// longer alignments, potentially with gaps.
template <class AlignmentCompare>
std::vector<Alignment> chain_alignments(std::vector<Alignment>&& alignments,
                                        std::string_view query,
                                        std::string_view rc_query,
                                        const DBGAlignerConfig &config,
                                        const DeBruijnGraph &graph);

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_CHAINER_HPP__
