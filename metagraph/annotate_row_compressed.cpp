#include "annotate_row_compressed.hpp"

#include <string>
#include <algorithm>
#include <stdexcept>

#include "serialization.hpp"

using libmaus2::util::NumberSerialisation;
using libmaus2::util::StringSerialisation;


namespace annotate {

RowCompressed::RowCompressed(const DBG_succ &graph)
      : RowCompressed(graph.num_edges() + 1) {
    graph_ = &graph;
}

RowCompressed::RowCompressed(uint64_t graph_size)
      : graph_(NULL),
        index_to_ids_(graph_size) {}

RowCompressed::SetStr RowCompressed::get(Index i) const {
    assert(i < index_to_ids_.size());

    SetStr label;
    for (auto &j : index_to_ids_[i]) {
        label.insert(id_to_label_[j]);
    }

    return label;
}

std::vector<uint64_t> RowCompressed::get_row(Index i) const {
    assert(i < index_to_ids_.size());

    std::vector<uint64_t> row((id_to_label_.size() + 63) >> 6);

    for (auto &j : index_to_ids_[i]) {
        row[j >> 6] |= (1llu << (j & 0x3F));
    }

    return row;
}

bool RowCompressed::has_label(Index i, const SetStr &label) const {
    assert(i < index_to_ids_.size());

    for (auto &column_name : label) {
        if (!has_label(i, column_name))
            return false;
    }
    return true;
}

bool RowCompressed::has_label(Index i, const std::string &label) const {
    assert(i < index_to_ids_.size());

    auto it = label_to_id_.find(label);
    if (it == label_to_id_.end())
        return false;

    return find(i, it->second).second;
}

void RowCompressed::set_label(Index i, const SetStr &label) {
    assert(i < index_to_ids_.size());

    for (const auto &value : label) {
        add_label(i, value);
    }
}

void RowCompressed::add_labels(const std::string &sequence, const SetStr &labels) {
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

void RowCompressed::add_label(Index i, const std::string &label) {
    assert(i < index_to_ids_.size());

    auto id_it = label_to_id_.find(label);
    size_t curr_id;

    if (id_it == label_to_id_.end()) {
        // current label does not exist yet -> assign new column
        curr_id = id_to_label_.size();

        label_to_id_[label] = curr_id;
        id_to_label_.push_back(label);
    } else {
        curr_id = id_it->second;
    }

    auto row_it = find(i, curr_id);
    // ensures the ids are in order
    if (!row_it.second) {
        index_to_ids_[i].insert(row_it.first, curr_id);
    }
}

// write annotation to disk
void RowCompressed::serialize(const std::string &filename) const {

    std::ofstream outstream(filename + ".row.annodbg");
    NumberSerialisation::serialiseNumber(outstream, index_to_ids_.size());
    serialize_string_number_map(outstream, label_to_id_);
    StringSerialisation::serialiseStringVector(outstream, id_to_label_);
    std::vector<uint64_t> full_vector;
    for (auto &indices : index_to_ids_) {
        for (auto &j : indices) {
            full_vector.push_back(j + 1);
        }
        full_vector.push_back(0);
    }
    serialize_number_vector(outstream,
                            full_vector,
                            std::log2(label_to_id_.size() + 1) + 1);
}

// read annotation from disk
bool RowCompressed::load(const std::string &filename) {

    std::ifstream instream(filename + ".row.annodbg");
    if (!instream.good())
        return false;

    size_t graph_size_ = NumberSerialisation::deserialiseNumber(instream);
    label_to_id_ = load_string_number_map(instream);
    id_to_label_ = StringSerialisation::deserialiseStringVector(instream);
    index_to_ids_.clear();
    index_to_ids_.resize(graph_size_);
    size_t j = 0;
    auto full_vector = load_number_vector<uint32_t>(instream);
    for (size_t i = 0; i < full_vector.size(); ++i) {
        if (!full_vector[i]) {
            j++;
            continue;
        }
        index_to_ids_[j].insert(index_to_ids_[j].end(), full_vector[i] - 1);
    }
    return true;
}

// assumes the indices are in order
std::pair<RowCompressed::IndexContainer::const_iterator, bool> RowCompressed::find(Index i, uint32_t id) const {
    IndexContainer::const_iterator it = index_to_ids_[i].begin();
    while (it != index_to_ids_[i].end() && *it < id)
        ++it;
    return std::make_pair(it, it != index_to_ids_[i].end() ? *it == id : false);
}

} // namespace annotate
