#ifndef __ALIGNER_CHAINER_HPP__
#define __ALIGNER_CHAINER_HPP__

#include "alignment.hpp"


namespace mtg {
namespace graph {
namespace align {

class IDBGAligner;
class SeedFilteringExtender;
class ISeeder;

typedef std::vector<std::pair<Alignment, int64_t>> Chain;
typedef Alignment::score_t score_t;


// Given forward and reverse-complement seeds, construct and call chains of seeds
std::pair<size_t, size_t>
call_seed_chains_both_strands(const IDBGAligner &aligner,
                              std::string_view forward,
                              std::string_view reverse,
                              const DBGAlignerConfig &config,
                              std::vector<Seed>&& fwd_seeds,
                              std::vector<Seed>&& bwd_seeds,
                              const std::function<void(Chain&&, score_t)> &callback,
                              const std::function<bool(Alignment::Column)> &skip_column
                                  = [](Alignment::Column) { return false; },
                              const std::function<bool()> &terminate
                                  = []() { return false; });

std::tuple<size_t, size_t, size_t>
chain_and_filter_seeds(const IDBGAligner &aligner,
                       std::shared_ptr<ISeeder> &seeder,
                       SeedFilteringExtender&& extender);

void chain_alignments(const IDBGAligner &aligner,
                      const std::vector<Alignment> &alignments,
                      const std::function<void(Alignment&&)> &callback);

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_CHAINER_HPP__
