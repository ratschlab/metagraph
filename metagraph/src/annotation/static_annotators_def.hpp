#ifndef __STATIC_ANNOTATOR_DEFS_HPP__
#define __STATIC_ANNOTATOR_DEFS_HPP__

#include "annotate_static.hpp"
#include "bin_rel_wt.hpp"
#include "bin_rel_wt_sdsl.hpp"
#include "rainbowfish.hpp"
#include "BRWT.hpp"


namespace annotate {

const char kRowPackedExtension[] = ".flat.annodbg";
const char kRainbowfishExtension[] = ".rbfish.annodbg";
const char kBRWTExtension[] = ".brwt.annodbg";
const char kBinRelWT_sdslExtension[] = ".bin_rel_wt_sdsl.annodbg";
const char kBinRelWTExtension[] = ".bin_rel_wt.annodbg";


typedef StaticBinRelAnnotator<RowConcatenated<>, std::string> RowFlatAnnotator;

typedef StaticBinRelAnnotator<Rainbowfish, std::string> RainbowfishAnnotator;

template <typename Label = std::string>
using BRWTCompressed = StaticBinRelAnnotator<BRWT, Label>;

typedef StaticBinRelAnnotator<BinRelWT_sdsl, std::string> BinRelWT_sdslAnnotator;

typedef StaticBinRelAnnotator<BinRelWT, std::string> BinRelWTAnnotator;

} // namespace annotate

#endif // __STATIC_ANNOTATOR_DEFS_HPP__
