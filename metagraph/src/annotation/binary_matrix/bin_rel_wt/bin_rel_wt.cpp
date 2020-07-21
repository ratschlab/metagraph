#include "bin_rel_wt.hpp"

#include <cassert>
#include <iterator>

#include "common/serialization.hpp"
#include "common/vectors/bit_vector.hpp"


namespace mtg {
namespace annot {
namespace binmat {

typedef brwt::binary_relation::object_id object_id;
typedef brwt::binary_relation::label_id label_id;


label_id to_label_id(uint64_t x) {
    return static_cast<label_id>(x + 1);
}

object_id to_object_id(uint64_t x) {
    return static_cast<object_id>(x);
}

uint64_t to_index(object_id y) {
    return static_cast<uint64_t>(y);
}

uint64_t to_index(label_id y) {
    return static_cast<uint64_t>(y) - 1;
}

bool BinRelWT::is_zero_row(Row row) const {
    auto first_label = to_label_id(0);
    auto last_label = to_label_id(max_used_label);
    auto cur_object = to_object_id(row);
    return (binary_relation_.size() == 0 || row > max_used_object
        || binary_relation_.obj_rank(cur_object, first_label, last_label) == 0);
}

bool BinRelWT::is_zero_column(Column column) const {
    auto first_object = to_object_id(0);
    auto last_object = to_object_id(max_used_object);
    auto cur_label = to_label_id(column);
    return (binary_relation_.size() == 0 || column > max_used_label
        || binary_relation_.rank(first_object, last_object, cur_label) == 0);
}

BinRelWT::BinRelWT(std::vector<std::unique_ptr<bit_vector>>&& columns)
      : num_labels(columns.size()) {

    uint64_t num_relations = 0;
    for (const auto &col_ptr : columns) {
        num_relations += col_ptr->num_set_bits();
    }

    std::vector<brwt::binary_relation::pair_type> relation_pairs;
    relation_pairs.reserve(num_relations);

    num_objects = columns.size() ? columns[0]->size() : 0;

    for (size_t j = 0; j < columns.size(); ++j) {
        const auto &col = *columns[j];

        col.call_ones([&](auto i) {
            brwt::binary_relation::pair_type relation = {
                .object = to_object_id(i),
                .label = to_label_id(j)
            };

            relation_pairs.push_back(std::move(relation));

            max_used_label = std::max(max_used_label, static_cast<uint64_t>(j));
            max_used_object = std::max(max_used_object, static_cast<uint64_t>(i));
        });

        columns[j].reset();
    }

    binary_relation_ = brwt::binary_relation(std::move(relation_pairs));
    assert(max_used_label <= num_labels);
}

BinRelWT::BinRelWT(const std::function<void(const RowCallback &)> &generate_rows,
                   uint64_t num_relations,
                   uint64_t num_columns)
      : num_labels(num_columns) {
    std::vector<brwt::binary_relation::pair_type> relation_pairs;
    relation_pairs.reserve(num_relations);

    uint64_t row_counter = 0;
    generate_rows([&](const SetBitPositions &row_set_bits) {
        for (const auto &col_index : row_set_bits) {
            brwt::binary_relation::pair_type relation = {
                .object = to_object_id(row_counter),
                .label = to_label_id(col_index)};
            relation_pairs.push_back(relation);
            max_used_label = std::max(max_used_label, col_index);
            max_used_object = std::max(max_used_object, row_counter);
        }
        row_counter++;
    });
    num_objects = row_counter;
    binary_relation_ = brwt::binary_relation(std::move(relation_pairs));
    assert(max_used_label <= num_labels);
}

uint64_t BinRelWT::num_columns() const {
    return num_labels;
}

uint64_t BinRelWT::num_rows() const {
    return num_objects;
}

BinRelWT::SetBitPositions BinRelWT::get_row(Row row) const {
    assert(row < num_objects);
    if (is_zero_row(row)) {
        return {};
    }
    auto first_label = to_label_id(0);
    auto last_label = to_label_id(max_used_label);
    auto cur_object = to_object_id(row);

    auto num_relations_in_row = static_cast<uint64_t>(
        binary_relation_.count_distinct_labels(cur_object,
                                               cur_object,
                                               first_label,
                                               last_label)
    );

    SetBitPositions relations_in_row;
    relations_in_row.reserve(num_relations_in_row);
    for (uint64_t relation_it = 1; relation_it <= num_relations_in_row; ++relation_it) {
        auto element =
            binary_relation_.nth_element(cur_object,
                                         cur_object,
                                         first_label,
                                         relation_it,
                                         brwt::lab_major);
        assert(element);
        relations_in_row.push_back(to_index(element->label));
    }
    return relations_in_row;
}

std::vector<BinRelWT::Row> BinRelWT::get_column(Column column) const {
    assert(column < num_labels);
    if (is_zero_column(column)) {
        return {};
    }
    auto first_object = to_object_id(0);
    auto last_object = to_object_id(max_used_object);
    auto cur_label = to_label_id(column);
    auto num_relations_in_column = (binary_relation_.size() == 0 ? 0 :
        static_cast<uint64_t>(binary_relation_.obj_rank(last_object, cur_label)));

    std::vector<Row> relations_in_column;
    relations_in_column.reserve(num_relations_in_column);
    for (uint64_t relation_it = 1; relation_it <= num_relations_in_column; ++relation_it) {
        auto element = binary_relation_.nth_element(first_object,
                                                    last_object,
                                                    cur_label,
                                                    relation_it,
                                                    brwt::lab_major);
        relations_in_column.push_back(to_index(element->object));
    }
    return relations_in_column;
}

bool BinRelWT::get(Row row, Column column) const {
    assert(column < num_labels);
    assert(row < num_objects);
    if (is_zero_column(column) || is_zero_row(row)) {
        return 0;
    }
    auto cur_object = to_object_id(row);
    auto cur_label = to_label_id(column);
    auto first_relation_in_range = binary_relation_.nth_element(cur_object,
                                                                cur_object,
                                                                cur_label,
                                                                1,
                                                                brwt::lab_major);
    return (first_relation_in_range &&
            first_relation_in_range->object == cur_object &&
            first_relation_in_range->label == cur_label);
}

bool BinRelWT::load(std::istream &in) {
    if (!in.good())
        return false;
    try {
        num_labels = load_number(in);
        num_objects = load_number(in);
        max_used_label = load_number(in);
        max_used_object = load_number(in);
        return binary_relation_.load(in);
    } catch (...) {
        return false;
    }
}

void BinRelWT::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Bad stream");
    serialize_number(out, num_labels);
    serialize_number(out, num_objects);
    serialize_number(out, max_used_label);
    serialize_number(out, max_used_object);
    binary_relation_.serialize(out);
}

uint64_t BinRelWT::num_relations() const {
    return binary_relation_.size();
}

} // namespace binmat
} // namespace annot
} // namespace mtg
