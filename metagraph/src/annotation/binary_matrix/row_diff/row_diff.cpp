#include "row_diff.hpp"

#include <fstream>

#include "annotation/binary_matrix/column_sparse/column_major.hpp"

namespace mtg {
namespace annot {
namespace binmat {

template <class BaseMatrix>
void RowDiff<BaseMatrix>::serialize(const std::string &filename) const {
    std::ofstream f(filename, ios::binary);
    serialize(f);
    f.close();
}

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::load(const std::string &filename) {
    std::ifstream f(filename, ios::binary);
    bool result = load(f);
    f.close();

    load_terminal(anchors_filename_, &terminal_);

    return result;
}

template
class RowDiff<ColumnMajor>;

} // namespace binmat
} // namespace annot
} // namespace mtg
