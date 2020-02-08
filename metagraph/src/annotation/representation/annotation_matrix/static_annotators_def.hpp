#ifndef __STATIC_ANNOTATOR_DEFS_HPP__
#define __STATIC_ANNOTATOR_DEFS_HPP__

#include "annotation_matrix.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt_sdsl.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbowfish.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"


namespace annotate {

typedef StaticBinRelAnnotator<RowConcatenated<>, std::string> RowFlatAnnotator;

typedef StaticBinRelAnnotator<Rainbowfish, std::string> RainbowfishAnnotator;

typedef StaticBinRelAnnotator<BRWT, std::string> MultiBRWTAnnotator;

typedef StaticBinRelAnnotator<BinRelWT_sdsl, std::string> BinRelWT_sdslAnnotator;

typedef StaticBinRelAnnotator<BinRelWT, std::string> BinRelWTAnnotator;

template <>
inline const std::string RowFlatAnnotator::kExtension = ".flat.annodbg";
template <>
inline const std::string RainbowfishAnnotator::kExtension = ".rbfish.annodbg";
template <>
inline const std::string MultiBRWTAnnotator::kExtension = ".brwt.annodbg";
template <>
inline const std::string BinRelWT_sdslAnnotator::kExtension = ".bin_rel_wt_sdsl.annodbg";
template <>
inline const std::string BinRelWTAnnotator::kExtension = ".bin_rel_wt.annodbg";

} // namespace annotate

#endif // __STATIC_ANNOTATOR_DEFS_HPP__
