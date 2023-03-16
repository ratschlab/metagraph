#ifndef __MST_MATRIX__
#define __MST_MATRIX__

#include "annotation/binary_matrix/base/binary_matrix.hpp"

namespace mtg {
namespace annot {
namespace binmat {

using mtg::common::logger;

template <class BaseMatrix>
class MSTMatrix : public BinaryMatrix {
public:
    typedef sdsl::int_vector<64> parents_type;

    uint64_t num_columns() const override { return base_matrix_.num_columns(); }
    uint64_t num_rows() const override { return base_matrix_.num_rows(); }

    bool get(Row /*row*/, Column /*column*/) const override {
        logger->error("MSTMatrix<BaseMatrix>::get not implemented");
        exit(1);
        // mkokot_TODO: implement
    }

    SetBitPositions get_row(Row row) const override {
        return decode_row(row);
    }

    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override {
        std::vector<SetBitPositions> res(rows.size());
        for (size_t i = 0 ; i < rows.size() ; ++i)
            res[i] = get_row(rows[i]);
        return res;
    }
    std::vector<Row> get_column(Column /*column*/) const override {
        logger->error("MSTMatrix<BaseMatrix>::get_column not implemented");
        exit(1);
        // mkokot_TODO: implement
    }

    bool load(std::istream &in) override  {
        num_relations_ = load_number(in);
        auto pos = in.tellg();
        if (!base_matrix_.load(in))
            return false;

        //mkokot_TODO: remove
        logger->trace("mst matrix, size of base matrix: {}", in.tellg() - pos);
        pos = in.tellg();

        parents_.load(in);

        //mkokot_TODO: remove
        logger->trace("mst matrix, size of base parents: {}", in.tellg() - pos);

        //mkokot_TODO: remove
        logger->trace("mst matrix, underlying matrix num_relations: {}", base_matrix_.num_relations());

        return true;
    }

    void serialize(std::ostream &/*out*/) const override {
        logger->error("MSTMatrix<BaseMatrix>::serialize  not implemented"); //mkokot_TODO: remove
        exit(1);
        // mkokot_TODO: implement
    }

    uint64_t num_relations() const override { return num_relations_; }

private:
    uint64_t num_relations_ = 0;
    BaseMatrix base_matrix_;
    parents_type parents_;

    SetBitPositions get_xor(const SetBitPositions& v1, const SetBitPositions& v2) const {
        SetBitPositions res;

        size_t i1 = 0;
        size_t i2 = 0;
        while (i1 < v1.size() && i2 < v2.size()) {
            if(v1[i1] == v2[i2]) {
                ++i1, ++i2;
            } else {
                if (v1[i1] < v2[i2])
                    res.push_back(v1[i1++]);
                else
                    res.push_back(v2[i2++]);
            }
        }
        while (i1 < v1.size())
            res.push_back(v1[i1++]);
        while (i2 < v2.size())
            res.push_back(v2[i2++]);

        return res;
    }
    SetBitPositions decode_row(Row row) const {
         size_t recu_depth{};
        auto res = decode_row_count_depth(row, recu_depth);
        //n_sum_decode_depth += recu_depth; //mkokot_TODO: may be used to count the recu depth
        return res;
    }
    SetBitPositions decode_row_count_depth(Row row, size_t& recu_depth) const {
        if (parents_[row] == row) {
            return base_matrix_.get_row(row);
        }
        ++recu_depth;
        return get_xor(decode_row_count_depth(parents_[row], recu_depth), base_matrix_.get_row(row));
    }
};

} // namespace binmat
} // namespace annot
} // namespace mtg


#endif // __MST_MATRIX__
