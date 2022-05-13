#ifndef __STATIC_ANNOTATOR_DEFS_HPP__
#define __STATIC_ANNOTATOR_DEFS_HPP__

#include <sdsl/dac_vector.hpp>

#include "annotation_matrix.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt_sdsl.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt_disk.hpp"
#include "annotation/binary_matrix/multi_brwt/row_sparse_brwt_disk.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbowfish.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbow.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "annotation/binary_matrix/row_sparse_disk/row_sparse_disk.hpp"
#include "annotation/binary_matrix/row_vector/unique_row_binmat.hpp"
#include "annotation/int_matrix/rank_extended/csc_matrix.hpp"
#include "annotation/int_matrix/row_diff/int_row_diff.hpp"
#include "annotation/int_matrix/row_diff/tuple_row_diff.hpp"
#include "annotation/int_matrix/csr_matrix/csr_matrix.hpp"
#include "annotation/int_matrix/rank_extended/tuple_csc_matrix.hpp"


namespace mtg {
namespace annot {

typedef StaticBinRelAnnotator<binmat::RowConcatenated<>, std::string> RowFlatAnnotator;

typedef StaticBinRelAnnotator<binmat::RowSparse, std::string> RowSparseAnnotator;

typedef StaticBinRelAnnotator<binmat::RowSparseBRWT_Disk, std::string> RowSparseBRWTDiskAnnotator;

typedef StaticBinRelAnnotator<binmat::RowSparseDisk, std::string> RowSparseDiskAnnotator;

typedef StaticBinRelAnnotator<binmat::Rainbowfish, std::string> RainbowfishAnnotator;

typedef StaticBinRelAnnotator<binmat::BRWT, std::string> MultiBRWTAnnotator;

typedef StaticBinRelAnnotator<binmat::BRWT_Disk, std::string> MultiBRWTDiskAnnotator;

typedef StaticBinRelAnnotator<binmat::BinRelWT_sdsl, std::string> BinRelWT_sdslAnnotator;

typedef StaticBinRelAnnotator<binmat::BinRelWT, std::string> BinRelWTAnnotator;

typedef StaticBinRelAnnotator<binmat::UniqueRowBinmat, std::string> UniqueRowAnnotator;

typedef StaticBinRelAnnotator<binmat::Rainbow<binmat::BRWT>, std::string> RbBRWTAnnotator;

typedef StaticBinRelAnnotator<binmat::RowDiff<binmat::ColumnMajor>, std::string> RowDiffColumnAnnotator;

typedef StaticBinRelAnnotator<binmat::RowDiff<binmat::BRWT>, std::string> RowDiffBRWTAnnotator;

typedef StaticBinRelAnnotator<binmat::RowDiff<binmat::BRWT_Disk>, std::string> RowDiffBRWTDiskAnnotator;

typedef StaticBinRelAnnotator<binmat::RowDiff<binmat::RowSparseBRWT_Disk>, std::string> RowDiffRowSparseBRWTDiskAnnotator;

typedef StaticBinRelAnnotator<binmat::RowDiff<binmat::RowSparse>, std::string> RowDiffRowSparseAnnotator;

typedef StaticBinRelAnnotator<binmat::RowDiff<binmat::RowSparseDisk>, std::string> RowDiffRowSparseDiskAnnotator;

typedef sdsl::dac_vector_dp<> CountsVector;

typedef StaticBinRelAnnotator<matrix::CSCMatrix<binmat::BRWT, CountsVector>, std::string> IntMultiBRWTAnnotator;

typedef StaticBinRelAnnotator<matrix::IntRowDiff<matrix::CSCMatrix<binmat::BRWT, CountsVector>>, std::string> IntRowDiffBRWTAnnotator;

typedef StaticBinRelAnnotator<matrix::CSRMatrix, std::string> IntRowAnnotator;

typedef StaticBinRelAnnotator<matrix::TupleCSCMatrix<binmat::ColumnMajor>, std::string> ColumnCoordAnnotator;

typedef StaticBinRelAnnotator<matrix::TupleCSCMatrix<binmat::BRWT>, std::string> MultiBRWTCoordAnnotator;

typedef StaticBinRelAnnotator<matrix::TupleRowDiff<matrix::TupleCSCMatrix<binmat::ColumnMajor>>, std::string> RowDiffCoordAnnotator;

typedef StaticBinRelAnnotator<matrix::TupleRowDiff<matrix::TupleCSCMatrix<binmat::BRWT>>, std::string> RowDiffBRWTCoordAnnotator;


template <>
inline const std::string RowFlatAnnotator::kExtension = ".flat.annodbg";
template <>
inline const std::string RowSparseAnnotator::kExtension = ".row_sparse.annodbg";
template <>
inline const std::string RowSparseDiskAnnotator::kExtension = ".row_sparse_disk.annodbg";
template <>
inline const std::string RainbowfishAnnotator::kExtension = ".rbfish.annodbg";
template <>
inline const std::string MultiBRWTAnnotator::kExtension = ".brwt.annodbg";
template <>
inline const std::string MultiBRWTDiskAnnotator::kExtension = ".brwt_disk.annodbg";
template <>
inline const std::string BinRelWT_sdslAnnotator::kExtension = ".bin_rel_wt_sdsl.annodbg";
template <>
inline const std::string BinRelWTAnnotator::kExtension = ".bin_rel_wt.annodbg";
template <>
inline const std::string UniqueRowAnnotator::kExtension = ".unique_row.annodbg";
template <>
inline const std::string RbBRWTAnnotator::kExtension = ".rb_brwt.annodbg";
template <>
inline const std::string RowDiffColumnAnnotator::kExtension = ".row_diff.annodbg";
template <>
inline const std::string RowDiffBRWTAnnotator::kExtension = ".row_diff_brwt.annodbg";
template <>
inline const std::string RowDiffBRWTDiskAnnotator::kExtension = ".row_diff_brwt_disk.annodbg";
template <>
inline const std::string RowSparseBRWTDiskAnnotator::kExtension = ".row_sparse_brwt_disk.annodbg";
template <>
inline const std::string RowDiffRowSparseBRWTDiskAnnotator::kExtension = ".row_diff_row_sparse_brwt_disk.annodbg";
template <>
inline const std::string RowDiffRowSparseAnnotator::kExtension = ".row_diff_sparse.annodbg";
template <>
inline const std::string RowDiffRowSparseDiskAnnotator::kExtension = ".row_diff_sparse_disk.annodbg";
template <>
inline const std::string IntMultiBRWTAnnotator::kExtension = ".int_brwt.annodbg";
template <>
inline const std::string IntRowDiffBRWTAnnotator::kExtension = ".row_diff_int_brwt.annodbg";
template <>
inline const std::string IntRowAnnotator::kExtension = ".int_csr.annodbg";
template <>
inline const std::string ColumnCoordAnnotator::kExtension = ".column_coord.annodbg";
template <>
inline const std::string MultiBRWTCoordAnnotator::kExtension = ".brwt_coord.annodbg";
template <>
inline const std::string RowDiffCoordAnnotator::kExtension = ".row_diff_coord.annodbg";
template <>
inline const std::string RowDiffBRWTCoordAnnotator::kExtension = ".row_diff_brwt_coord.annodbg";

} // namespace annot
} // namespace mtg

#endif // __STATIC_ANNOTATOR_DEFS_HPP__
