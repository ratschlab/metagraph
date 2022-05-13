#include "row_sparse_brwt_disk.hpp"

#include <queue>
#include <numeric>

#include "common/algorithms.hpp"
#include "common/serialization.hpp"
#include "common/utils/template_utils.hpp"


namespace mtg {
namespace annot {
namespace binmat {

bool RowSparseBRWT_Disk::get(Row row, Column column) const {
    logger->error("RowSparseBRWT_Disk::get not implemented {} {}", row, column);
    exit(1);
    // mkokot, TODO: implement
}

RowSparseBRWT_Disk::SetBitPositions RowSparseBRWT_Disk::get_row(Row row) const {
    if (stored_as_row_sparse[row]) {
        ++n_row_sparse_queried_rows;
        Row row_in_row_sparse = stored_as_row_sparse.rank1(row) - 1;
        return row_sparse_disk.get_row(row_in_row_sparse);
    } else {
        ++n_brwt_queried_rows;
        Row row_in_brwt = stored_as_row_sparse.rank0(row) - 1;
        return brwt_disk.get_row(row_in_brwt);
    }    
}

// mkokot, TODO: przepisac to zeby wywolywac get_rows z brwt_disk i z row_sparse_disk
std::vector<RowSparseBRWT_Disk::SetBitPositions>
RowSparseBRWT_Disk::get_rows(const std::vector<Row> &row_ids) const {   
    std::vector<Row> rows_row_sparse;
    std::vector<Row> rows_brwt;

    for (size_t i = 0 ; i < row_ids.size(); ++i) {
        if (stored_as_row_sparse[row_ids[i]]) {
            ++n_row_sparse_queried_rows;
            rows_row_sparse.push_back(stored_as_row_sparse.rank1(row_ids[i]) - 1);            
        } else {
            ++n_brwt_queried_rows;            
            rows_brwt.push_back(stored_as_row_sparse.rank0(row_ids[i]) - 1);            
        }
    }
    auto res_row_sparse = row_sparse_disk.get_rows(rows_row_sparse);
    auto res_brwt_disk = brwt_disk.get_rows(rows_brwt); // mkokot, TODO: just checking if it still crashes without this
    std::vector<RowSparseBRWT_Disk::SetBitPositions> res(row_ids.size());
    size_t pos_row_sprase{};
    size_t pos_brwt{};
    for (size_t i = 0 ; i < row_ids.size(); ++i) {
        if (stored_as_row_sparse[row_ids[i]])
            res[i] = std::move(res_row_sparse[pos_row_sprase++]);
        else        
            res[i] = std::move(res_brwt_disk[pos_brwt++]);
    }

    return res;
}

std::vector<RowSparseBRWT_Disk::Row> RowSparseBRWT_Disk::get_column(Column column) const {
    assert(column < num_columns());
    logger->error("RowSparseBRWT_Disk::get_column not implemented {}", column);
    exit(1);
    // mkokot, TODO: implement
}

bool RowSparseBRWT_Disk::load(std::istream &in) {
    return  stored_as_row_sparse.load(in) &&
            row_sparse_disk.load(in) &&
            brwt_disk.load(in);    
}

void RowSparseBRWT_Disk::serialize(std::ostream &out) const {
    logger->error("RowSparseBRWT_Disk::serialize not implemented {}", !!out);
    exit(1);
    // mkokot, TODO: implement
}

uint64_t RowSparseBRWT_Disk::num_relations() const {
    return brwt_disk.num_relations() + row_sparse_disk.num_relations();    
}

} // namespace binmat
} // namespace annot
} // namespace mtg
