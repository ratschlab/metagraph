#include "row_diff.hpp"

#include "annotation/binary_matrix/column_sparse/column_major.hpp"

namespace mtg {
namespace annot {
namespace binmat {

template <class BaseMatrix>
void RowDiff<BaseMatrix>::serialize(const std::string &name) const {
    std::ofstream f(name, ios::binary);
    serialize(f);
    f.close();
}

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::load(const std::string &name) {
    std::ifstream f(name, ios::binary);
    bool result = load(f);
    f.close();

    load_terminal(terminal_file_, &terminal_);

    return result;
}

template
class RowDiff<ColumnMajor>;

} // namespace binmat
} // namespace annot
} // namespace mtg
