#include "annotation.hpp"

#include "common/serialization.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace annot {

using mtg::common::logger;

template <typename Label>
size_t LabelEncoder<Label>::insert_and_encode(const Label &label) {
    auto it = encode_label_.emplace(label).first;
    return it - encode_label_.begin();
}

template <typename Label>
size_t LabelEncoder<Label>::encode(const Label &label) const {
    auto it = encode_label_.find(label);
    if (it == encode_label_.end())
        throw std::out_of_range("Label not found");
    return it - encode_label_.begin();
}

template<>
void LabelEncoder<std::string>::serialize(std::ostream &outstream) const {
    outstream.write("LE-v2.0", 7);
    Serializer serializer(outstream);
    encode_label_.serialize(serializer);
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
        auto pos = instream.tellg();
        std::string version(7, '\0');
        if (instream.read(version.data(), 7) && version == "LE-v2.0") {
            Deserializer deserializer(instream);
            encode_label_ = VectorSet<std::string>::deserialize(deserializer, true);
            return instream.good();
        }
        // backward compatibility
        instream.seekg(pos);
        {
            std::vector<std::string> labels;
            if (!load_string_vector(instream, &labels))
                return false;
            std::vector<uint64_t> values;
            if (!load_number_vector(instream, &values))
                return false;
        }
        std::vector<std::string> labels;
        if (!load_string_vector(instream, &labels))
            return false;
        encode_label_ = VectorSet<std::string>(labels.begin(), labels.end());
        return instream.good();
    } catch (...) {
        return false;
    }
}


// For each pair (first, second) in the dictionary, renames column |first|
// to |second| and merges columns with matching names, if supported.
template <typename LabelType>
void MultiLabelAnnotation<LabelType>
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
void MultiLabelAnnotation<LabelType>
::call_objects(const Label &label, std::function<void(Index)> callback) const {
    if (!get_label_encoder().label_exists(label))
        return;

    for (Index index : get_matrix().get_column(label_encoder_.encode(label))) {
        callback(index);
    }
}

// TODO: remove?
template <typename LabelType>
typename MultiLabelAnnotation<LabelType>::VLabels
MultiLabelAnnotation<LabelType>::get_labels(Index i) const {
    assert(i < this->num_objects());

    auto label_codes = get_matrix().get_rows({ i })[0];

    VLabels labels(label_codes.size());

    for (size_t j = 0; j < label_codes.size(); ++j) {
        labels[j] = label_encoder_.decode(label_codes[j]);
    }

    return labels;
}

template <typename LabelType>
void MultiLabelAnnotation<LabelType>::add_label_counts(const std::vector<Index> &,
                                                       const VLabels &,
                                                       const std::vector<uint64_t> &) {
    logger->error("Adding label counts is not implemented for this annotator");
    exit(1);
}

template <typename LabelType>
void MultiLabelAnnotation<LabelType>::add_label_coord(Index, const VLabels &, uint64_t) {
    logger->error("Adding relation attributes is not implemented for this annotator");
    exit(1);
}

template <typename LabelType>
void MultiLabelAnnotation<LabelType>::add_label_coords(const std::vector<std::pair<Index, uint64_t>> &,
                                                       const VLabels &) {
    logger->error("Adding relation attributes is not implemented for this annotator");
    exit(1);
}

template class MultiLabelAnnotation<std::string>;

template class LabelEncoder<std::string>;

} // namespace annot
} // namespace mtg
