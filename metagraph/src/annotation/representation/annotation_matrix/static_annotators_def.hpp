#ifndef __STATIC_ANNOTATOR_DEFS_HPP__
#define __STATIC_ANNOTATOR_DEFS_HPP__

#include <sdsl/dac_vector.hpp>

#include "annotation_matrix.hpp"
#include "annotation/binary_matrix/bin_rel_wt/bin_rel_wt.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbowfish.hpp"
#include "annotation/binary_matrix/rainbowfish/rainbow.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "annotation/binary_matrix/row_disk/row_disk.hpp"
#include "annotation/binary_matrix/row_vector/unique_row_binmat.hpp"
#include "annotation/int_matrix/rank_extended/csc_matrix.hpp"
#include "annotation/int_matrix/row_diff/int_row_diff.hpp"
#include "annotation/int_matrix/row_diff/tuple_row_diff.hpp"
#include "annotation/int_matrix/csr_matrix/csr_matrix.hpp"
#include "annotation/int_matrix/row_disk/int_row_disk.hpp"
#include "annotation/int_matrix/row_disk/coord_row_disk.hpp"
#include "annotation/int_matrix/rank_extended/tuple_csc_matrix.hpp"


namespace mtg {
namespace annot {

typedef StaticBinRelAnnotator<matrix::RowFlat<>, std::string> RowFlatAnnotator;

typedef StaticBinRelAnnotator<matrix::RowSparse, std::string> RowSparseAnnotator;

typedef StaticBinRelAnnotator<matrix::Rainbowfish, std::string> RainbowfishAnnotator;

typedef StaticBinRelAnnotator<matrix::BRWT, std::string> MultiBRWTAnnotator;

typedef StaticBinRelAnnotator<matrix::BinRelWT, std::string> BinRelWTAnnotator;

typedef StaticBinRelAnnotator<matrix::UniqueRowBinmat, std::string> UniqueRowAnnotator;

typedef StaticBinRelAnnotator<matrix::Rainbow<matrix::BRWT>, std::string> RbBRWTAnnotator;

typedef StaticBinRelAnnotator<matrix::RowDiff<matrix::ColumnMajor>, std::string> RowDiffColumnAnnotator;

typedef StaticBinRelAnnotator<matrix::RowDiff<matrix::BRWT>, std::string> RowDiffBRWTAnnotator;

typedef StaticBinRelAnnotator<matrix::RowDiff<matrix::RowFlat<>>, std::string> RowDiffRowFlatAnnotator;

typedef StaticBinRelAnnotator<matrix::RowDiff<matrix::RowSparse>, std::string> RowDiffRowSparseAnnotator;

typedef StaticBinRelAnnotator<matrix::RowDiff<matrix::RowDisk>, std::string> RowDiffDiskAnnotator;

typedef sdsl::dac_vector_dp<> CountsVector;

typedef StaticBinRelAnnotator<matrix::CSCMatrix<matrix::BRWT, CountsVector>, std::string> IntMultiBRWTAnnotator;

typedef StaticBinRelAnnotator<matrix::IntRowDiff<matrix::CSCMatrix<matrix::BRWT, CountsVector>>, std::string> IntRowDiffBRWTAnnotator;

typedef StaticBinRelAnnotator<matrix::CSRMatrix, std::string> IntRowAnnotator;

typedef StaticBinRelAnnotator<matrix::IntRowDiff<matrix::IntRowDisk>, std::string> IntRowDiffDiskAnnotator;

typedef StaticBinRelAnnotator<matrix::TupleRowDiff<matrix::CoordRowDisk>, std::string> RowDiffDiskCoordAnnotator;

typedef StaticBinRelAnnotator<matrix::TupleCSCMatrix<matrix::ColumnMajor>, std::string> ColumnCoordAnnotator;

typedef StaticBinRelAnnotator<matrix::TupleCSCMatrix<matrix::BRWT>, std::string> MultiBRWTCoordAnnotator;

typedef StaticBinRelAnnotator<matrix::TupleRowDiff<matrix::TupleCSCMatrix<matrix::ColumnMajor>>, std::string> RowDiffCoordAnnotator;

typedef StaticBinRelAnnotator<matrix::TupleRowDiff<matrix::TupleCSCMatrix<matrix::BRWT>>, std::string> RowDiffBRWTCoordAnnotator;


template <>
inline const std::string RowFlatAnnotator::kExtension = ".flat.annodbg";
template <>
inline const std::string RowSparseAnnotator::kExtension = ".row_sparse.annodbg";
template <>
inline const std::string RainbowfishAnnotator::kExtension = ".rbfish.annodbg";
template <>
inline const std::string MultiBRWTAnnotator::kExtension = ".brwt.annodbg";
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
inline const std::string RowDiffRowFlatAnnotator::kExtension = ".row_diff_flat.annodbg";
template <>
inline const std::string RowDiffRowSparseAnnotator::kExtension = ".row_diff_sparse.annodbg";
template <>
inline const std::string RowDiffDiskAnnotator::kExtension = ".row_diff_disk.annodbg";
template <>
inline const std::string IntMultiBRWTAnnotator::kExtension = ".int_brwt.annodbg";
template <>
inline const std::string IntRowDiffBRWTAnnotator::kExtension = ".row_diff_int_brwt.annodbg";
template <>
inline const std::string IntRowAnnotator::kExtension = ".int_csr.annodbg";
template <>
inline const std::string IntRowDiffDiskAnnotator::kExtension = ".row_diff_int_disk.annodbg";
template <>
inline const std::string RowDiffDiskCoordAnnotator::kExtension = ".row_diff_disk_coord.annodbg";
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
