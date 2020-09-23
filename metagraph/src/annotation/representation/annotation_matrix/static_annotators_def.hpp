#ifndef __STATIC_ANNOTATOR_DEFS_HPP__
#define __STATIC_ANNOTATOR_DEFS_HPP__

#include "annotation_matrix.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt_sdsl.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbowfish.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbow.hpp"
#include "annotation/binary_matrix/row_diff/column_diff.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/binary_matrix/row_vector/unique_row_binmat.hpp"


namespace mtg {
namespace annot {

typedef StaticBinRelAnnotator<binmat::RowConcatenated<>, std::string> RowFlatAnnotator;

typedef StaticBinRelAnnotator<binmat::Rainbowfish, std::string> RainbowfishAnnotator;

typedef StaticBinRelAnnotator<binmat::BRWT, std::string> MultiBRWTAnnotator;

typedef StaticBinRelAnnotator<binmat::BinRelWT_sdsl, std::string> BinRelWT_sdslAnnotator;

typedef StaticBinRelAnnotator<binmat::BinRelWT, std::string> BinRelWTAnnotator;

typedef StaticBinRelAnnotator<binmat::UniqueRowBinmat, std::string> UniqueRowAnnotator;

typedef StaticBinRelAnnotator<binmat::Rainbow<binmat::BRWT>, std::string> RbBRWTAnnotator;

typedef StaticBinRelAnnotator<binmat::RowDiff, std::string> RowDiffAnnotator;

typedef StaticBinRelAnnotator<binmat::ColumnDiff<binmat::ColumnMajor>, std::string> ColumnDiffAnnotator;


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
template <>
inline const std::string UniqueRowAnnotator::kExtension = ".unique_row.annodbg";
template <>
inline const std::string RbBRWTAnnotator::kExtension = ".rb_brwt.annodbg";
template <>
inline const std::string RowDiffAnnotator::kExtension = ".row_diff.annodbg";
template <>
inline const std::string ColumnDiffAnnotator::kExtension = ".column_diff.annodbg";
} // namespace annot
} // namespace mtg

#endif // __STATIC_ANNOTATOR_DEFS_HPP__
