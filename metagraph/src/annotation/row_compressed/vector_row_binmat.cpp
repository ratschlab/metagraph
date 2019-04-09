#include "vector_row_binmat.hpp"

#include <cstdio>
#include <sdsl/int_vector_buffer.hpp>

#include "threading.hpp"


bool VectorRowBinMat::get(Row row, Column column) const {
    assert(row < vector_.size());
    return std::find(vector_[row].begin(), vector_[row].end(), column)
                != vector_[row].end();
}

void VectorRowBinMat::set(Row row, Column column) {
    assert(row < vector_.size());

    if (!get(row, column))
        vector_[row].push_back(column);

    if (column >= num_columns_)
        num_columns_ = column + 1;
}

void VectorRowBinMat::force_set(Row row, Column column) {
    assert(row < vector_.size());

    vector_[row].push_back(column);

    if (column >= num_columns_)
        num_columns_ = column + 1;
}

void VectorRowBinMat::standardize_rows() {
    #pragma omp parallel for num_threads(get_num_threads())
    for (size_t i = 0; i < vector_.size(); ++i) {
        std::sort(vector_[i].begin(), vector_[i].end());
        vector_[i].erase(std::unique(vector_[i].begin(), vector_[i].end()),
                         vector_[i].end());
    }
}

std::vector<VectorRowBinMat::Column>
VectorRowBinMat::get_row(Row row) const {
    assert(row < vector_.size());
    const auto &v = vector_[row];
    return std::vector<Column>(v.begin(), v.end());
}

void VectorRowBinMat::clear_row(Row row) {
    assert(row < vector_.size());
    vector_[row].clear();
}

std::vector<VectorRowBinMat::Row>
VectorRowBinMat::get_column(Column column) const {
    std::vector<Row> result;
    for (uint64_t i = 0; i < vector_.size(); ++i) {
        if (get(i, column))
            result.push_back(i);
    }
    return result;
}

void VectorRowBinMat::insert_rows(const std::vector<Row> &rows) {
    assert(std::is_sorted(rows.begin(), rows.end()));
    utils::insert_default_values(rows, &vector_);
}

bool VectorRowBinMat::load(std::istream &instream) {
    if (!instream.good())
        return false;

    try {
        size_t num_rows = load_number(instream);
        num_columns_ = load_number(instream);
        vector_.clear();
        vector_.shrink_to_fit();
        vector_.resize(num_rows);

        sdsl::int_vector<> full_vector;
        full_vector.load(instream);

        for (size_t k = 0, i = 0; k < full_vector.size(); ++k) {
            if (full_vector[k]) {
                vector_[i].push_back(full_vector[k] - 1);
                if (vector_[i].back() >= num_columns_)
                    return false;
            } else {
                i++;
            }
        }

        return true;
    } catch (...) {
        return false;
    }
}

void VectorRowBinMat::serialize(std::ostream &outstream) const {
    serialize_number(outstream, num_rows());
    serialize_number(outstream, num_columns());

    sdsl::int_vector<> full_vector(num_relations() + num_rows(),
                                   0,
                                   utils::code_length(num_columns()));

    for (uint64_t i = 0, p = 0; i < vector_.size(); ++i) {
        for (uint64_t value : vector_[i]) {
            full_vector[p++] = value + 1;
        }
        full_vector[p++] = 0;
    }

    full_vector.serialize(outstream);
}

// number of ones in the matrix
uint64_t VectorRowBinMat::num_relations() const {
    return std::accumulate(
        vector_.begin(), vector_.end(), uint64_t(0),
        [](uint64_t sum, const auto &v) { return sum + v.size(); }
    );
}

// matrix density
double VectorRowBinMat::density() const {
    return static_cast<double>(num_relations()) / num_columns() / num_rows();
}

VectorRowBinMat::StreamRows::StreamRows(std::ifstream &instream, const std::string &filename) {
    if (!instream.good())
        throw std::ifstream::failure("Bad stream");

    (void)load_number(instream);
    (void)load_number(instream);

    uint64_t offs = instream.tellg();
    inbuf_ = new sdsl::int_vector_buffer<>(filename,
                                           std::ios::in,
                                           1024*1024,
                                           0,
                                           false,
                                           offs);
}

std::unique_ptr<std::vector<VectorRowBinMat::Row> > VectorRowBinMat::StreamRows::next_row() {

    auto row_vector = std::make_unique<std::vector<VectorRowBinMat::Row> >();

    while (i_ < inbuf_->size()) {
        auto value = (*inbuf_)[i_];
        i_++;
        if (value) {
            row_vector->push_back(value - 1);
        } else {
            return row_vector;
        }
    }
    return nullptr;
}

uint64_t VectorRowBinMat::write_rows(std::ofstream &outstream,
                const std::string &filename,
                const std::function<void (const std::function<void (const std::vector<uint64_t> &)>&)> &callback,
                uint64_t num_cols) {

    uint8_t width = utils::code_length(num_cols);
    uint64_t num_rows = 0;

    // write dummy num_rows value to fill in later
    uint64_t header_offs = outstream.tellp();
    serialize_number(outstream, 0);
    serialize_number(outstream, num_cols);
    outstream.flush();

    uint64_t iv_offs = outstream.tellp();
    auto outbuf = sdsl::int_vector_buffer<>(filename,
                                            std::ios::out,
                                            1024*1024,
                                            width,
                                            false,
                                            iv_offs);

    callback([&](const std::vector<uint64_t> &row) {
        for(auto val : row)
            outbuf.push_back(val + 1);
        outbuf.push_back(0);
        num_rows++;
    });

    outbuf.close();

    outstream.seekp(header_offs);
    serialize_number(outstream, num_rows);

    return num_rows;
}
