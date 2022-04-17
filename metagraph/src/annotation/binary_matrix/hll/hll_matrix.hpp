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

template <template <class T> class Hash = IntHash>
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
        num_set_bits_.reserve(columns.size());
        for (size_t j = 0; j < columns_.size(); ++j) {
            num_set_bits_.push_back(columns[j]->num_set_bits());
            num_relations_ += num_set_bits_.back();
            columns[j]->call_ones([&](size_t i) { columns_[j].insert(hasher_(i)); });
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

    template <class Columns>
    std::pair<double, double>
    estimate_column_union_intersection_cardinality(const Columns &a,
                                                   const Columns &b) const {
        if (a.empty() && b.empty())
            return {};

        double a_cardinality_est;
        std::shared_ptr<ColumnSketch> a_sketch;
        if (a.empty()) {
            a_sketch = std::make_shared<ColumnSketch>();
            a_cardinality_est = 0.0;
        } else {
            a_sketch = std::make_shared<ColumnSketch>(merge_columns(a));
            a_cardinality_est = a_sketch->estimate_cardinality();
        }

        double b_cardinality_est;
        std::shared_ptr<const ColumnSketch> b_sketch;
        if (b.empty()) {
            b_sketch = std::make_shared<const ColumnSketch>();
            b_cardinality_est = 0.0;
        } else if (b.size() == 1) {
            b_sketch = std::shared_ptr<const ColumnSketch>{
                std::shared_ptr<const ColumnSketch>{}, &get_column_sketch(b[0])
            };
            b_cardinality_est = num_set_bits_[b[0]];
        } else {
            b_sketch = std::make_shared<const ColumnSketch>(merge_columns(b));
            b_cardinality_est = b_sketch->estimate_cardinality();
        }

        a_sketch->merge(*b_sketch);
        double union_est = a_sketch->estimate_cardinality();
        return { union_est, a_cardinality_est + b_cardinality_est - union_est };
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

        try {
            Deserializer ds(in);
            precision_ = ds.operator()<double>();
            num_rows_ = load_number(in);
            num_relations_ = load_number(in);
        } catch (...) {
            return false;
        }

        if (!load_number_vector(in, &num_set_bits_))
            return false;

        columns_.resize(num_set_bits_.size());
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
        serialize_number_vector(out, num_set_bits_);
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
    std::vector<uint64_t> num_set_bits_;
    size_t num_rows_;
    size_t num_relations_;
    Hash<Row> hasher_;

    inline uint64_t hash(Row row) { return hasher_(row); }
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __HLL_MATRIX_HPP__
