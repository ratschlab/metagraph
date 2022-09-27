#include "row_sparse_disk.hpp"

#include "common/logger.hpp"
#include "common/utils/file_utils.hpp"
#include "common/ifstream_with_name.hpp"
#include "graph/representation/succinct/boss.hpp"


namespace mtg {
namespace annot {
namespace binmat {

using mtg::common::logger;

bool RowSparseDisk::Impl::get(Row row, Column column) const {
    SetBitPositions set_bits = get_row(row);
    SetBitPositions::iterator v = std::lower_bound(set_bits.begin(), set_bits.end(), column);
    return v != set_bits.end() && *v == column;
}

std::vector<BinaryMatrix::Row> RowSparseDisk::Impl::get_column(uint64_t num_rows,
                                                               Column column) const {
    std::vector<Row> result;
    for (Row row = 0; row < num_rows; ++row) {
        if (get(row, column))
            result.push_back(row);
    }
    return result;
}

BinaryMatrix::SetBitPositions RowSparseDisk::Impl::get_row(Row row) const {
    assert(boundary_[boundary_.size() - 1] == 1);
    uint64_t start_idx = row == 0 ? 0 : boundary_.select1(row) + 1;
    uint64_t end_idx = boundary_.next1(start_idx);
    SetBitPositions result;
    result.reserve(end_idx - start_idx);
    // In each row, the first value in `set_bits_` stores the first set bit,
    // and all next ones store deltas pos[i] - pos[i-1].
    uint64_t last = 0;
    for (uint64_t i = start_idx; i != end_idx; ++i) {
        result.push_back(last += set_bits_[i - row]);
    }
    return result;
}

std::vector<BinaryMatrix::SetBitPositions>
RowSparseDisk::Impl::get_rows(const std::vector<Row> &row_ids) const {
    std::vector<SetBitPositions> rows(row_ids.size());

    for (size_t i = 0; i < row_ids.size(); ++i) {
        rows[i] = get_row(row_ids[i]);
    }

    return rows;
}

bool RowSparseDisk::load(std::istream &f) {
    auto _f = dynamic_cast<mtg::common::IfstreamWithName *>(&f);
    assert(_f);
    try {
        num_columns_ = load_number(f);

        auto boundary_start = load_number(f);

        buffer_params_.filename = _f->get_name();
        buffer_params_.offset = f.tellg();

        assert(boundary_start >= buffer_params_.offset);
        iv_size_on_disk_ = boundary_start - buffer_params_.offset;

        f.seekg(boundary_start, ios_base::beg);

        boundary_.load(f);

        num_rows_ = boundary_.num_set_bits();

    } catch (...) {
        return false;
    }

    return true;
}

void RowSparseDisk::serialize(std::ostream & f) const {
    serialize_number(f, num_columns_);
    auto boundary_start = iv_size_on_disk_ + f.tellp() + sizeof(size_t);
    serialize_number(f, boundary_start);

    size_t chunk_size = 1024 * 1024;

    std::ifstream iv_in(buffer_params_.filename, std::ios::binary);
    if (!iv_in.good())
        throw std::ofstream::failure("Cannot read from file " + buffer_params_.filename);

    iv_in.seekg(buffer_params_.offset);
    std::vector<char> chunk(chunk_size);

    // Perform binary rewrite of int_vector_buffer
    auto left_to_read = iv_size_on_disk_;
    while (left_to_read) {
        auto to_read = left_to_read < chunk_size ? left_to_read : chunk_size;
        iv_in.read(chunk.data(), to_read);
        f.write(chunk.data(), to_read);
        left_to_read -= to_read;
    }
    iv_in.close();

    assert(f.tellp() == boundary_start);
    boundary_.serialize(f);
}


void RowSparseDisk::serialize(const std::function<void(binmat::BinaryMatrix::RowCallback)> &call_rows,
                              const std::string &filename,
                              uint64_t num_cols,
                              uint64_t num_set_bits,
                              uint64_t num_rows) {
    // std::ios::ate needed because labels are serialized before
    // std::ios::in needed for tellp to return absolute file position
    std::ofstream outstream(filename,
                            std::ios::binary | std::ios::ate
                                    | std::ios::in);
    if (!outstream.good())
        throw std::ofstream::failure("Cannot write to file " + filename);

    serialize_number(outstream, num_cols);
    auto boundary_start_pos = outstream.tellp();

    // boundary is built in memory while writing rows data to disk via int_vector_buffer
    // the boundary must be loaded in the RowSparseDisk::load so its position must be known
    // binary format for this representation is
    // [boundary pos] (8B)
    // [serialized rows]
    // [boundary] (that are stored in the file at boundary pos)

    // write "empty" boundary start
    std::streampos boundary_start = 0;
    serialize_number(outstream, boundary_start);

    const uint64_t iv_offs = outstream.tellp();
    outstream.close();

    auto outbuf = sdsl::int_vector_buffer<>(filename, std::ios::out | std::ios::binary,
                                            1024 * 1024, sdsl::bits::hi(num_cols) + 1,
                                            false, iv_offs);
    // call_bits should be called for each argument that should be set (row boundaries)
    auto call_bits = [&](const std::function<void(uint64_t)> &call_bit) {
        uint64_t t = 0;
        call_rows([&](const auto &row) {
            uint64_t last = 0;
            for (auto val : row) {
                outbuf.push_back(val - last);
                last = val;
            }
            call_bit((t += row.size() + 1) - 1);
        });
        outbuf.close();
    };

    bit_vector_small boundary(call_bits, num_rows + num_set_bits, num_rows);

    outstream.open(filename, std::ios::in | std::ios::out | std::ios::binary | std::ios::ate);
    boundary_start = outstream.tellp();
    boundary.serialize(outstream);
    outstream.seekp(boundary_start_pos, ios_base::beg);
    serialize_number(outstream, boundary_start);
}

} // namespace binmat
} // namespace annot
} // namespace mtg
