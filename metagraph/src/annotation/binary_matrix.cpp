#include "binary_matrix.hpp"

#include "utils.hpp"


template <typename RowType>
StreamRows<RowType>::StreamRows(const std::string &filename, size_t offset) {
    std::ifstream instream(filename, std::ios::binary);

    if (!instream.good() || !instream.seekg(offset).good())
        throw std::ifstream::failure("Cannot read rows from file " + filename);

    (void)load_number(instream);
    (void)load_number(instream);

    inbuf_ = sdsl::int_vector_buffer<>(filename,
                                       std::ios::in | std::ios::binary,
                                       1024 * 1024,
                                       0,
                                       false,
                                       instream.tellg());
}

template <typename RowType>
RowType* StreamRows<RowType>::next_row() {
    row_.clear();

    while (i_ < inbuf_.size()) {
        auto value = inbuf_[i_++];
        if (value - 1 > std::numeric_limits<typename RowType::value_type>::max())
            throw std::ifstream::failure("Integer overflow: trying to read too"
                                         " large column index: " + std::to_string(value - 1));
        if (value) {
            row_.push_back(value - 1);
        } else {
            return &row_;
        }
    }
    return nullptr;
}

template class StreamRows<BinaryMatrix::SetBitPositions>;


void append_row_major(const std::string &filename,
                      const std::function<void(BinaryMatrix::RowCallback)> &call_rows,
                      uint64_t num_cols) {
    std::ofstream outstream(filename, std::ios::binary | std::ios::app);

    uint64_t num_rows = 0;

    // write dummy num_rows value to fill in later
    const uint64_t header_offs = outstream.tellp();
    serialize_number(outstream, 0);
    serialize_number(outstream, num_cols);
    const uint64_t iv_offs = outstream.tellp();
    outstream.close();

    {
        auto outbuf = sdsl::int_vector_buffer<>(filename,
                                                std::ios::out | std::ios::binary,
                                                1024 * 1024,
                                                utils::code_length(num_cols),
                                                false,
                                                iv_offs);

        call_rows([&](const auto &row) {
            for (auto val : row) {
                outbuf.push_back(val + 1);
            }
            outbuf.push_back(0);
            num_rows++;
        });

        outbuf.close();
    }

    outstream.open(filename, std::ios::in | std::ios::out | std::ios::binary);
    outstream.seekp(header_offs);
    serialize_number(outstream, num_rows);
}
