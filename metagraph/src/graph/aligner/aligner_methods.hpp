#ifndef __DBG_ALIGNER_METHODS_HPP__
#define __DBG_ALIGNER_METHODS_HPP__


#include "aligner_helper.hpp"


typedef Alignment<DeBruijnGraph::node_index> DBGAlignment;

typedef std::function<std::vector<DBGAlignment>(const DeBruijnGraph&,
                                                const DBGAlignerConfig&,
                                                const char*, // begin
                                                const char*, // end
                                                size_t, // clipping
                                                bool // orientation
                                                )> Seeder;

typedef std::function<void(const DeBruijnGraph&,
                           const DBGAlignment&,
                           std::vector<DBGAlignment>*, // output vector
                           const char*, // end iterator
                           const DBGAlignerConfig&,
                           bool, // orientation
                           typename DBGAlignment::score_t // min path score
                           )> Extender;

typedef std::function<bool(const DBGAlignment&, const DBGAlignment&)> PriorityFunction;


Seeder make_unimem_seeder(const std::vector<DeBruijnGraph::node_index> &nodes);

extern const Seeder default_seeder;
extern const Extender default_extender;

#endif // __DBG_ALIGNER_METHODS_HPP__
