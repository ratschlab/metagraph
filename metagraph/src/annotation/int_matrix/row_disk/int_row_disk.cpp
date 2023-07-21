#include "int_row_disk.hpp"

#include "common/logger.hpp"
#include "common/utils/file_utils.hpp"

namespace mtg {
namespace annot {
namespace matrix {

using mtg::common::logger;

std::vector<BinaryMatrix::Row> IntRowDisk::View::get_column(Column column) const {
    logger->warn("get_column is extremely inefficient for IntRowDisk, consider"
                 " using a column-major format");
    const uint64_t num_rows = boundary_.num_set_bits();
    std::vector<Row> result;
    for (Row row = 0; row < num_rows; ++row) {
        SetBitPositions set_bits = get_row(row);
        if (std::binary_search(set_bits.begin(), set_bits.end(), column))
            result.push_back(row);
    }
    return result;
}

BinaryMatrix::SetBitPositions IntRowDisk::View::get_row(Row row) const {
    assert(boundary_[boundary_.size() - 1] == 1);
    uint64_t begin = row == 0 ? 0 : boundary_.select1(row) + 1 - row;
    uint64_t end = boundary_.select1(row + 1) - row;

    SetBitPositions result;
    result.reserve(end - begin);
    // In each row, the first value in `set_bits_` stores the first set bit,
    // and all next ones store deltas pos[i] - pos[i-1].
    uint64_t last = 0;
    set_bits_.start_reading_at(begin * (bits_for_col_id_ + bits_for_value_));
    Column col_id;
    for (uint64_t i = begin; i != end; ++i) {
        col_id = set_bits_.get(bits_for_col_id_);
        set_bits_.get(bits_for_value_);
        result.push_back(last += col_id);
    }

    return result;
}

IntMatrix::RowValues IntRowDisk::View::get_row_values(Row row) const {
    assert(boundary_[boundary_.size() - 1] == 1);
    uint64_t begin = row == 0 ? 0 : boundary_.select1(row) + 1 - row;
    uint64_t end = boundary_.select1(row + 1) - row;

    RowValues result;
    result.reserve(end - begin);
    // In each row, the first value in `set_bits_` stores the first set bit,
    // and all next ones store deltas pos[i] - pos[i-1].
    uint64_t last = 0;
    set_bits_.start_reading_at(begin * (bits_for_col_id_ + bits_for_value_));
    Column col_id;
    Value val;
    for (uint64_t i = begin; i != end; ++i) {
        col_id = set_bits_.get(bits_for_col_id_);
        val = set_bits_.get(bits_for_value_);
        result.emplace_back(last += col_id, val);
    }

    return result;
}


std::vector<IntMatrix::RowValues>
IntRowDisk::View::get_row_values(const std::vector<Row> &row_ids) const {
    std::vector<RowValues> rows_with_values(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        rows_with_values[i] = get_row_values(row_ids[i]);
    }

    return rows_with_values;
}

bool IntRowDisk::load(std::istream &f) {
    auto _f = dynamic_cast<sdsl::mmap_ifstream *>(&f);
    assert(_f);
    try {
        num_columns_ = load_number(f);
        auto boundary_start = load_number(f);

        num_set_bits_ = load_number(f);
        bits_for_col_id_ = load_number(f);
        bits_for_value_ = load_number(f);

        buffer_params_.filename = _f->get_filename();
        buffer_params_.offset = f.tellg();

        assert(boundary_start >= buffer_params_.offset);

        f.seekg(boundary_start, std::ios_base::beg);

        boundary_.load(f);

        num_rows_ = boundary_.num_set_bits();

    } catch (...) {
        return false;
    }

    return true;
}

void IntRowDisk::serialize(
            const std::string &filename,
            const std::function<void(std::function<void(const RowValues &)>)> &call_rows,
            uint64_t num_cols,
            uint64_t num_rows,
            uint64_t num_set_bits,
            uint8_t max_val) {
    // std::ios::ate needed because labels are serialized before
    // std::ios::in needed for tellp to return absolute file position
    std::ofstream outstream(filename, std::ios::binary | std::ios::ate | std::ios::in);

    if (!outstream.good())
        throw std::ofstream::failure("Cannot write to file " + filename);

    serialize_number(outstream, num_cols);

    auto boundary_start_pos = outstream.tellp();

    // boundary is built in memory
    // the boundary must be loaded in the IntRowDisk::load so its position must be known
    // binary format for this representation is
    // [boundary pos] (8B)
    // [serialized rows with values]
    // [boundary] (that are stored in the file at boundary pos)

    // write "empty" boundary start
    std::streampos boundary_start = 0;
    serialize_number(outstream, boundary_start);

    serialize_number(outstream, num_set_bits);

    uint8_t bits_for_col_id = sdsl::bits::hi(num_cols) + 1;
    uint8_t bits_for_value = sdsl::bits::hi(max_val) + 1;

    serialize_number(outstream, bits_for_col_id);
    serialize_number(outstream, bits_for_value);

    common::DiskWriter writer(outstream, 1024 * 1024);

    //  call_bits should be called for each argument that should be set (row boundaries)
    auto call_bits = [&](const std::function<void(uint64_t)> &call_bit) {
        uint64_t t = 0;
        call_rows([&](const auto &row_with_values) {
            uint64_t last = 0;
            for (auto [col_id, val] : row_with_values) {
                writer.add(col_id - last, bits_for_col_id);
                writer.add(val, bits_for_value);
                last = col_id;
            }
            call_bit((t += (row_with_values.size() + 1)) - 1);
        });
        writer.flush();
    };

    bit_vector_small boundary(call_bits, num_rows + num_set_bits, num_rows);

    boundary_start = outstream.tellp();
    boundary.serialize(outstream);
    outstream.seekp(boundary_start_pos, std::ios_base::beg);
    serialize_number(outstream, boundary_start);
}

} // namespace matrix
} // namespace annot
} // namespace mtg
