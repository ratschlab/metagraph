#include "vector_row_binmat.hpp"
#include "binary_matrix.hpp"

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

VectorRowBinMat::StreamRows::StreamRows(std::ifstream &instream) {
    if (!instream.good())
        throw std::ifstream::failure("Bad stream");

    (void)load_number(instream);
    (void)load_number(instream);

    sdsl::int_vector<>::read_header(size_, width_, instream);
    //std::cout << size << std::endl;
    //std::cout << (int)width << std::endl;
    in_ = std::move(instream);
}

std::unique_ptr<std::vector<VectorRowBinMat::Row> > VectorRowBinMat::StreamRows::next_row() {
    auto row_vector = std::make_unique<std::vector<VectorRowBinMat::Row> >();

    VectorRowBinMat::Column value = 0;
    const size_t bits_per = 8 * sizeof(value);

    while(++i_ <= size_) {

        if (0 == (i_ - 1) % 8) {
            if (!in_.good())
                throw std::ifstream::failure("Bad stream");

            b_ = in_.get();
            if (in_.eof())
                throw std::ifstream::failure("Unexpected end of stream");
        } else {
            b_ >>= 1;
        }

        value >>= 1;
        if (b_ & 1)
            value += 1ul << (bits_per - 1);

        if (0 == i_ % width_) {
            value >>= bits_per - width_;
            if (value) {
                row_vector->push_back(value - 1);
            } else {
                return std::move(row_vector);
            }

            value = 0;
        }
    }
    return nullptr;
}

void VectorRowBinMat::write_rows(std::ofstream &outstream,
                const std::function<void (const std::function<void(void)>&, const std::function<void (const std::vector<uint64_t> &)>&, const std::function<void(void)>&)> &callback,
                uint64_t num_rows,
                uint64_t num_cols) {
    uint8_t width = utils::code_length(num_cols);
    //std::cout << num_rows << std::endl;
    //std::cout << num_cols << std::endl;
    serialize_number(outstream, num_rows);
    serialize_number(outstream, num_cols);

    auto size_pos = outstream.tellp();
    sdsl::int_vector<>::write_header(0, 0, outstream);

    if (num_rows <= 0)
        return;

    uint8_t b = 0;
    uint64_t i = 0;
    uint64_t value = 0;
    //uint64_t row_count = 0;
    uint64_t num_relations = 0;
    std::vector<uint64_t>::const_iterator iter;

    callback([&]() {
    }, [&](const std::vector<uint64_t> &row) {
        // std::for_each(row.begin(), row.end(), [](uint64_t s) { std::cout << s << std::endl; });
        // std::cout << "," << std::endl;
        iter = row.begin();

        if (iter != row.end()) {
            value = *iter + 1;
            iter++;
            num_relations++;
        }

        bool flag = false;
        while(++i) {

            if (value & 1)
                b += 1u<<7;
            value >>= 1;

            if (0 == i % 8) {
                outstream.put(b);
                b = 0;
            }

            b >>= 1;
            if (0 == i % width) {
                if (iter == row.end()) {
                    if(!flag) {
                        value = 0;
                        flag = true;
                    } else {
                        return;
                    }
                } else {
                    value = *iter + 1;
                    iter++;
                    num_relations++;
                }
            }
        }
    }, [&]() {
        b >>= 8 - (i % 8);
        outstream.put(b);
        outstream.put(0);

        outstream.seekp(size_pos);
        std::cout << "num_relations: " << num_relations << ", num_rows: " << num_rows << std::endl;
        sdsl::int_vector<>::write_header((num_relations + num_rows) * width, width, outstream);
    });
}
