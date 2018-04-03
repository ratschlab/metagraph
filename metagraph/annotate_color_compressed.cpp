#include "annotate_color_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"

using libmaus2::util::NumberSerialisation;
using libmaus2::util::StringSerialisation;


namespace annotate {

ColorCompressed::ColorCompressed(const DBG_succ &graph, size_t cache_size)
      : ColorCompressed(graph.num_edges() + 1, cache_size) {
    graph_ = &graph;
}

ColorCompressed::ColorCompressed(uint64_t graph_size, size_t cache_size)
      : graph_(NULL),
        graph_size_(graph_size),
        cached_colors_(
            cache_size,
            caches::LRUCachePolicy<uint32_t>(),
            [this](uint32_t curr_id, sdsl::bit_vector *annotation_curr) {
                this->flush(curr_id, annotation_curr);
                delete annotation_curr;
            }
        ) {}

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

    for (uint64_t i = 0; i < graph_size_; ++i) {
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

    auto it = label_to_id_.find(label);
    if (it == label_to_id_.end())
        return false;

    if (cached_colors_.Cached(it->second))
        const_cast<ColorCompressed*>(this)->flush(it->second,
                                                  cached_colors_.Get(it->second));

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
    size_t curr_id;

    if (id_it == label_to_id_.end()) {
        // current label does not exist yet -> assign new column
        curr_id = id_to_label_.size();
        cached_colors_.Put(curr_id, new sdsl::bit_vector(graph_size_, 0));

        label_to_id_[label] = curr_id;
        id_to_label_.push_back(label);

    } else {
        curr_id = id_it->second;

        if (!cached_colors_.Cached(id_it->second))
            cached_colors_.Put(curr_id, new sdsl::bit_vector(graph_size_, 0));
    }

    assert(cached_colors_.Cached(curr_id));
    (*cached_colors_.Get(curr_id))[i] = 1;
}

// write annotation to disk
void ColorCompressed::serialize(const std::string &filename) const {
    const_cast<ColorCompressed*>(this)->flush();

    std::ofstream outstream(filename);
    NumberSerialisation::serialiseNumber(outstream, graph_size_);
    serialize_string_number_map(outstream, label_to_id_);
    StringSerialisation::serialiseStringVector(outstream, id_to_label_);
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

    graph_size_ = NumberSerialisation::deserialiseNumber(instream);
    label_to_id_ = load_string_number_map(instream);
    id_to_label_ = StringSerialisation::deserialiseStringVector(instream);
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
    cached_colors_.Clear();
    for (auto *column : bitmatrix_) {
        if (column)
            delete column;
    }
    bitmatrix_.clear();
}

void ColorCompressed::flush() {
    for (auto it = cached_colors_.Begin(); it != cached_colors_.End(); ++it) {
        flush(it->first, it->second);
    }
}

void ColorCompressed::flush(uint32_t curr_id,
                            sdsl::bit_vector *annotation_curr) {
    assert(annotation_curr);

    if (bitmatrix_.size() <= curr_id) {
        bitmatrix_.resize(curr_id + 1, NULL);
    }
    assert(curr_id < bitmatrix_.size());

    if (bitmatrix_[curr_id]) {
        // decompress existing column
        std::unique_ptr<sdsl::bit_vector> tmp(inflate_column(curr_id));
        *annotation_curr |= *tmp;
        delete bitmatrix_[curr_id];
    }

    bitmatrix_[curr_id] = new sdsl::sd_vector<>(*annotation_curr);
}

sdsl::bit_vector* ColorCompressed::inflate_column(uint32_t id) const {
    auto *col = bitmatrix_.at(id);

    assert(col);

    sdsl::bit_vector *result = new sdsl::bit_vector(col->size(), 0);

    auto slct = sdsl::select_support_sd<>(col);
    auto rank = sdsl::rank_support_sd<>(col);
    auto num_set_bits = rank(col->size());

    for (uint64_t i = 1; i <= num_set_bits; ++i) {
        uint64_t idx = slct(i);
        if (idx >= col->size())
            break;
        (*result)[idx] = 1;
    }

    return result;
}

} // namespace annotate
