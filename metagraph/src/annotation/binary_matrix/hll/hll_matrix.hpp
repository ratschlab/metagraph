#ifndef __HLL_MATRIX_HPP__
#define __HLL_MATRIX_HPP__

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/hashers/hll_counter.hpp"
#include "common/serialization.hpp"

namespace mtg {
namespace annot {
namespace binmat {

template <typename T>
struct IntHash {
    inline uint64_t operator()(const T &x) const;
};

template<>
inline uint64_t IntHash<uint64_t>::operator()(const uint64_t &in) const {
    uint64_t x = in;
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
}

template <template <class T> class Hash = std::hash>
class HLLMatrix : public BinaryMatrix {
  public:
    typedef HLLCounter ColumnSketch;

    HLLMatrix() {}
    HLLMatrix(const std::vector<std::unique_ptr<bit_vector>> &columns, double precision)
          : precision_(precision) {
        if (columns.empty()) {
            num_rows_ = 0;
            num_relations_ = 0;
            return;
        }

        num_rows_ = columns[0]->size();
        columns_.resize(columns.size(), precision_);
        for (size_t j = 0; j < columns_.size(); ++j) {
            std::cerr << "before: " << columns_[j].estimate_cardinality() << " + " << columns[j]->num_set_bits() << "\tafter: ";
            num_relations_ += columns[j]->num_set_bits();
            columns[j]->call_ones([&](size_t i) { columns_[j].insert(i); });
            std::cerr << columns_[j].estimate_cardinality() << "\n";
        }
    }

    inline const ColumnSketch& get_column_sketch(Column column) const {
        assert(column < columns_.size());
        return columns_[column];
    }

    template <class Columns>
    ColumnSketch merge_columns(const Columns &columns) const {
        if (columns.empty())
            return {};

        ColumnSketch result = get_column_sketch(columns[0]);
        for (size_t i = 1; i < columns.size(); ++i) {
            result.merge(get_column_sketch(columns[i]));
        }

        return result;
    }

    uint64_t num_columns() const override { return columns_.size(); }
    uint64_t num_rows() const override { return num_rows_; }

    bool get(Row row, Column column) const override {
        return get_column_sketch(column).check(hasher_(row));
    }

    SetBitPositions get_row(Row row) const override {
        SetBitPositions result;
        uint64_t row_hash = hasher_(row);
        for (Column i = 0; i < columns_.size(); ++i) {
            if (columns_[i].check(row_hash))
                result.emplace_back(i);
        }

        return result;
    }

    std::vector<Row> get_column(Column column) const override {
        std::vector<uint64_t> hashes;
        hashes.reserve(num_rows_);
        for (size_t i = 0; i < num_rows_; ++i) {
            hashes.push_back(hasher_(i));
        }

        std::vector<Row> result;
        get_column_sketch(column).check(hashes.data(), hashes.data() + hashes.size(),
                                        [&](Row row) { result.push_back(row); });

        return result;
    }

    bool load(std::istream &in) override {
        if (!in.good())
            return false;

        size_t num_columns;
        try {
            Deserializer ds(in);
            precision_ = ds.operator()<double>();
            num_rows_ = load_number(in);
            num_relations_ = load_number(in);
            num_columns = load_number(in);
        } catch (...) {
            return false;
        }

        columns_.resize(num_columns);
        for (size_t i = 0; i < columns_.size(); ++i) {
            if (!columns_[i].load(in))
                return false;
        }

        return true;
    }

    void serialize(std::ostream &out) const override {
        Serializer s(out);
        s(precision_);
        serialize_number(out, num_rows_);
        serialize_number(out, num_relations_);
        serialize_number(out, columns_.size());
        for (size_t i = 0; i < columns_.size(); ++i) {
            columns_[i].serialize(out);
        }
    }

    // number of ones in the matrix
    uint64_t num_relations() const override { return num_relations_; }

    auto& data() { return columns_; }
    const auto& data() const { return columns_; }

  private:
    double precision_;
    std::vector<ColumnSketch> columns_;
    size_t num_rows_;
    size_t num_relations_;
    Hash<Row> hasher_;

    inline uint64_t hash(Row row) { return hasher_(row); }
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __HLL_MATRIX_HPP__
