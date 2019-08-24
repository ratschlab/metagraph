#ifndef __DBG_ALIGNER_METHODS_HPP__
#define __DBG_ALIGNER_METHODS_HPP__


#include "aligner_helper.hpp"


template <typename NodeType>
using Seeder = std::function<std::vector<Alignment<NodeType>>(
    const DeBruijnGraph&,
    const DBGAlignerConfig&,
    const char*, // begin
    const char*, // end
    size_t, // clipping
    bool // orientation
)>;

template <typename NodeType, class... Args>
using SeederBuilder = std::function<Seeder<NodeType>(Args&&... args)>;

template <typename NodeType>
using MapExtendSeederBuilder = SeederBuilder<NodeType,
                                             const std::vector<NodeType>&,
                                             const DeBruijnGraph&>;


template <typename NodeType, class... Args>
using Extender = std::function<void(
    const DeBruijnGraph&,
    const Alignment<NodeType>&,
    std::vector<Alignment<NodeType>>*, // output vector
    const char*, // end iterator
    const DBGAlignerConfig&,
    bool, // orientation
    typename Alignment<NodeType>::score_t // min path score
)>;

template <typename NodeType>
using PriorityFunction = std::function<bool(const Alignment<NodeType>&,
                                            const Alignment<NodeType>&)>;


template <typename NodeType>
Seeder<NodeType>
build_mem_seeder(const std::vector<NodeType> &nodes,
                 std::function<bool(NodeType,
                                    const DeBruijnGraph &)> stop_matching,
                 const DeBruijnGraph &graph);

template <typename NodeType>
Seeder<NodeType> build_unimem_seeder(const std::vector<NodeType> &nodes,
                                     const DeBruijnGraph &graph);

template <typename NodeType>
std::vector<Alignment<NodeType>>
exact_seeder(const DeBruijnGraph &graph,
             const DBGAlignerConfig &config,
             const char *seed_begin,
             const char *seed_end,
             size_t clipping = 0,
             bool orientation = false);

template <typename NodeType>
std::vector<Alignment<NodeType>>
suffix_seeder(const DeBruijnGraph &graph,
               const DBGAlignerConfig &config,
               const char *seed_begin,
               const char *seed_end,
               size_t clipping = 0,
               bool orientation = false);

template <typename NodeType>
void default_extender(const DeBruijnGraph &graph,
                      const Alignment<NodeType> &path,
                      std::vector<Alignment<NodeType>> *next_paths, // output vector
                      const char *sequence_end, // end iterator
                      const DBGAlignerConfig &config,
                      bool orientation, // orientation
                      typename Alignment<NodeType>::score_t min_path_score
                          = std::numeric_limits<score_t>::min());

#endif // __DBG_ALIGNER_METHODS_HPP__
