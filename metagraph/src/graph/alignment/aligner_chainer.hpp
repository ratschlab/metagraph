#ifndef __ALIGNER_CHAINER_HPP__
#define __ALIGNER_CHAINER_HPP__

#include "alignment.hpp"


namespace mtg {
namespace graph {
namespace align {

typedef std::vector<std::pair<Alignment, int64_t>> Chain;
typedef Alignment::score_t score_t;

// Given forward and reverse-complement seeds, construct and call chains of seeds
std::pair<size_t, size_t>
call_seed_chains_both_strands(std::string_view forward,
                              std::string_view reverse,
                              size_t node_overlap,
                              const DBGAlignerConfig &config,
                              std::vector<Alignment>&& fwd_seeds,
                              std::vector<Alignment>&& bwd_seeds,
                              const std::function<void(Chain&&, score_t)> &callback);

// Given a set of local alignments, use sparse dynamic programming to construct
// longer alignments, potentially with gaps.
template <class AlignmentCompare>
std::vector<Alignment> chain_alignments(std::vector<Alignment>&& alignments,
                                        std::string_view query,
                                        std::string_view rc_query,
                                        const DBGAlignerConfig &config,
                                        size_t node_overlap);

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_CHAINER_HPP__
