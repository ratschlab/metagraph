#include "annotation.hpp"

#include "common/serialization.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace annot {

using mtg::common::logger;


template <typename Label>
size_t LabelEncoder<Label>::insert_and_encode(const Label &label) {
    auto [it, inserted] = encode_label_.emplace(label, decode_label_.size());
    if (inserted)
        decode_label_.push_back(label);

    return it->second;
}

template<>
void LabelEncoder<std::string>::serialize(std::ostream &outstream) const {
    serialize_string_number_map(outstream, encode_label_);
    serialize_string_vector(outstream, decode_label_);
}

template<typename Label>
void LabelEncoder<Label>::merge(const LabelEncoder<Label> &other) {
    for (size_t i = 0; i < other.size(); ++i) {
        insert_and_encode(other.decode(i));
    }
}

template<>
bool LabelEncoder<std::string>::load(std::istream &instream) {
    if (!instream.good())
        return false;

    try {
        if (!load_string_number_map(instream, &encode_label_))
            return false;

        if (!load_string_vector(instream, &decode_label_))
            return false;

        return instream.good();
    } catch (...) {
        return false;
    }
}


// For each pair (first, second) in the dictionary, renames column |first|
// to |second| and merges columns with matching names, if supported.
template <typename LabelType>
void MultiLabelEncoded<LabelType>
::rename_labels(const tsl::hopscotch_map<Label, Label> &dict) {
    // old labels
    std::vector<Label> index_to_label = label_encoder_.get_labels();

    // new labels
    for (const auto &pair : dict) {
        try {
            index_to_label[label_encoder_.encode(pair.first)] = pair.second;
        } catch (const std::out_of_range &) {
            logger->warn("Label '{}' not found, instruction '{} -> {}' skipped",
                         pair.first, pair.first, pair.second);
        }
    }

    label_encoder_.clear();

    // insert new column labels
    for (const auto &label : index_to_label) {
        if (label_encoder_.label_exists(label)) {
            // no exception -> there already exists a column with this name
            logger->error("Detected multiple labels renamed to '{}'"
                          ". Annotation merge is not implemented"
                          " for this annotation type.", label);
            exit(1);
        }

        // this is the first column with this name
        label_encoder_.insert_and_encode(label);
    }
}

template <typename LabelType>
void MultiLabelEncoded<LabelType>
::call_objects(const Label &label, std::function<void(Index)> callback) const {
    if (!label_exists(label))
        return;

    for (Index index : get_matrix().get_column(label_encoder_.encode(label))) {
        callback(index);
    }
}

// calls get_row(i)
template <typename LabelType>
typename MultiLabelEncoded<LabelType>::VLabels
MultiLabelEncoded<LabelType>::get(Index i) const {
    assert(i < this->num_objects());

    const auto &label_codes = get_matrix().get_row(i);

    VLabels labels(label_codes.size());

    for (size_t j = 0; j < label_codes.size(); ++j) {
        labels[j] = label_encoder_.decode(label_codes[j]);
    }

    return labels;
}

template <typename IndexType, typename LabelType>
void MultiLabelAnnotation<IndexType, LabelType>
::add_label_counts(const std::vector<Index> &,
                   const VLabels &,
                   const std::vector<uint32_t> &) {
    logger->error("Adding label counts is not implemented for this annotator");
    exit(1);
}

template class MultiLabelEncoded<std::string>;

template class LabelEncoder<std::string>;

} // namespace annot
} // namespace mtg
