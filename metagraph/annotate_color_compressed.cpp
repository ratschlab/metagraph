#include "annotate_color_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"


namespace annotate {

// // Merge constructor
// ColorCompressed::ColorCompressed(const std::vector<ColorCompressed> &categories,
//                                  const std::vector<std::vector<size_t>> &merge_plan)
//       : annotation_curr_(NULL) {
//     std::cout << categories.size();
//     std::cout << merge_plan.size();
// }

ColorCompressed::SetStr ColorCompressed::get(Index i) const {
    assert(i < graph_size_);

    const_cast<ColorCompressed*>(this)->flush();

    SetStr label;
    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        assert(bitmatrix_[j]);
        if (bitmatrix_[j]->operator[](i)) {
            label.insert(id_to_label_[j]);
        }
    }
    return label;
}

std::vector<uint64_t> ColorCompressed::get_row(Index i) const {
    assert(i < graph_size_);

    std::vector<uint64_t> row((bitmatrix_.size() + 63) >> 6);

    for (size_t j = 0; j < bitmatrix_.size(); ++j) {
        assert(bitmatrix_[j]);
        if (bitmatrix_[j]->operator[](i)) {
            row[j >> 6] |= (1llu << (j % 64));
        }
    }
    return row;
}

void ColorCompressed::serialize_uncompressed_rows(const std::string &filename) const {
    const_cast<ColorCompressed*>(this)->flush();
    std::ofstream outstream(filename);
    libmaus2::util::NumberSerialisation::serialiseNumber(outstream, graph_size_);

    std::cout << graph_size_ << " " << bitmatrix_.size() << "\n";

    for (size_t i = 0; i < graph_size_; ++i) {
        auto row = get_row(i);
        libmaus2::util::NumberSerialisation::serialiseNumberVector(outstream, row);
    }
    outstream.close();
}

bool ColorCompressed::has_label(Index i, const SetStr &label) const {
    assert(i < graph_size_);

    for (auto &column_name : label) {
        if (!has_label(i, column_name))
            return false;
    }
    return true;
}

bool ColorCompressed::has_label(Index i, const std::string &label) const {
    assert(i < graph_size_);

    const_cast<ColorCompressed*>(this)->flush();

    auto it = label_to_id_.find(label);
    if (it == label_to_id_.end())
        return false;

    return bitmatrix_[it->second]->operator[](i);
}

void ColorCompressed::set_label(Index i, const SetStr &label) {
    assert(i < graph_size_);

    for (const auto &value : label) {
        add_label(i, value);
    }
}

void ColorCompressed::add_labels(const std::string &sequence, const SetStr &labels) {
    if (!graph_) {
        throw std::runtime_error("Please initialize with a graph\n");
    }
    graph_->align(sequence, [&](uint64_t i) {
        if (i > 0) {
            for (const auto &label : labels) {
                add_label(i, label);
            }
        }
    });
}

void ColorCompressed::add_label(Index i, const std::string &label) {
    assert(i < graph_size_);

    auto id_it = label_to_id_.find(label);

    if (id_it == label_to_id_.end()) {
        // current label does not exist yet -> assign new column
        if (annotation_curr_) {
            flush();
            delete annotation_curr_;
        }

        curr_id_ = id_to_label_.size();
        annotation_curr_ = new sdsl::bit_vector(graph_size_, 0);

        label_to_id_[label] = curr_id_;
        id_to_label_.push_back(label);
    } else if (!annotation_curr_ || curr_id_ != id_it->second) {
        if (annotation_curr_) {
            flush();
            delete annotation_curr_;
        }

        // decompress existing column
        curr_id_ = id_it->second;
        assert(bitmatrix_[curr_id_]);
        annotation_curr_ = inflate_column(curr_id_);
    }

    (*annotation_curr_)[i] = 1;
}

// write annotation to disk
void ColorCompressed::serialize(const std::string &filename) const {
    const_cast<ColorCompressed*>(this)->flush();

    std::ofstream outstream(filename);
    libmaus2::util::NumberSerialisation::serialiseNumber(outstream, graph_size_);
    serialize_string_number_map(outstream, label_to_id_);
    libmaus2::util::StringSerialisation::serialiseStringVector(outstream, id_to_label_);
    for (auto *column : bitmatrix_) {
        assert(column);
        column->serialize(outstream);
    }
}

// read annotation from disk
bool ColorCompressed::load(const std::string &filename) {
    release();

    std::ifstream instream(filename);
    if (!instream.good())
        return false;

    graph_size_ = libmaus2::util::NumberSerialisation::deserialiseNumber(instream);
    label_to_id_ = load_string_number_map(instream);
    id_to_label_ = libmaus2::util::StringSerialisation::deserialiseStringVector(instream);
    bitmatrix_.resize(id_to_label_.size(), NULL);
    for (auto &column : bitmatrix_) {
        column = new sdsl::sd_vector<>();
        column->load(instream);
    }
    return true;
}

void ColorCompressed::release() {
    label_to_id_.clear();
    id_to_label_.clear();
    for (auto *column : bitmatrix_) {
        if (column)
            delete column;
    }
    bitmatrix_.clear();
    if (annotation_curr_) {
        delete annotation_curr_;
        annotation_curr_ = NULL;
    }
}

void ColorCompressed::flush() {
    if (!annotation_curr_)
        return;

    if (bitmatrix_.size() <= curr_id_) {
        bitmatrix_.push_back(NULL);
    }
    assert(curr_id_ < bitmatrix_.size());

    if (bitmatrix_[curr_id_]) {
        sdsl::bit_vector *tmp = inflate_column(curr_id_);
        *annotation_curr_ |= *tmp;
        delete tmp;
        delete bitmatrix_[curr_id_];
    }

    bitmatrix_[curr_id_] = new sdsl::sd_vector<>(*annotation_curr_);
}

sdsl::bit_vector* ColorCompressed::inflate_column(uint32_t id) const {
    auto *col = bitmatrix_.at(id);

    assert(col);

    sdsl::bit_vector *result = new sdsl::bit_vector(col->size(), 0);

    auto slct = sdsl::select_support_sd<>(col);
    auto rank = sdsl::rank_support_sd<>(col);
    size_t num_set_bits = rank(col->size());

    for (size_t i = 1; i <= num_set_bits; ++i) {
        size_t idx = slct(i);
        if (idx >= col->size())
            break;
        (*result)[idx] = 1;
    }

    return result;
}

} // namespace annotate
