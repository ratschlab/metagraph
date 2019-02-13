#include "annotate.hpp"

#include <libmaus2/util/StringSerialisation.hpp>

#include "serialization.hpp"

using libmaus2::util::StringSerialisation;


namespace annotate {

template <typename Label>
size_t LabelEncoder<Label>::encode(const Label &label) const {
    auto it = encode_label_.find(label);
    if (it != encode_label_.end()) {
        return it->second;
    } else {
        throw std::runtime_error("ERROR: No such label");
    }
}

template <typename Label>
size_t LabelEncoder<Label>::insert_and_encode(const Label &label) {
    try {
        return encode(label);
    } catch (const std::runtime_error &e) {
        encode_label_[label] = decode_label_.size();
        decode_label_.push_back(label);
        return decode_label_.size() - 1;
    }
}

template<>
void LabelEncoder<std::string>::serialize(std::ostream &outstream) const {
    serialize_string_number_map(outstream, encode_label_);
    StringSerialisation::serialiseStringVector(outstream, decode_label_);
}

template<>
bool LabelEncoder<std::string>::load(std::istream &instream) {
    if (!instream.good())
        return false;

    try {
        load_string_number_map(instream, &encode_label_);
        decode_label_ = StringSerialisation::deserialiseStringVector(instream);

        return true;
    } catch (...) {
        return false;
    }
}


// For each pair (first, second) in the dictionary, renames column |first|
// to |second| and merges columns with matching names, if supported.
template <typename IndexType, typename LabelType>
void MultiLabelEncoded<IndexType, LabelType>
::rename_columns(const std::unordered_map<Label, Label> &dict) {
    std::vector<Label> index_to_label(label_encoder_.size());
    for (size_t i = 0; i < index_to_label.size(); ++i) {
        index_to_label[i] = label_encoder_.decode(i);
    }
    for (const auto &pair : dict) {
        index_to_label[label_encoder_.encode(pair.first)] = pair.second;
    }

    label_encoder_.clear();

    // insert new column labels
    for (const auto &label : index_to_label) {
        try {
            label_encoder_.encode(label);
            // no exception -> there already exists a column with this name
            std::cerr << "Error: detected more than one column with"
                      << " target name " << label
                      << ". Merging columns is not implemented"
                      << " for this annotation type."
                      << std::endl;
            exit(1);
        } catch (...) {
            // this is the first column with this name
            label_encoder_.insert_and_encode(label);
        }
    }
}


// Count all labels collected from the given rows
// and return top |num_top| with the their counts.
template <typename IndexType, typename LabelType>
auto MultiLabelEncoded<IndexType, LabelType>
::get_top_labels(const std::vector<Index> &indices, size_t num_top) const
-> std::vector<std::pair<Label, size_t>> {
    auto counter = count_labels(indices);

    std::vector<std::pair<size_t, size_t>> counts;
    for (size_t j = 0; j < counter.size(); ++j) {
        if (counter[j])
            counts.emplace_back(j, counter[j]);
    }
    // sort in decreasing order
    std::sort(counts.begin(), counts.end(),
              [](const auto &first, const auto &second) {
                  return first.second > second.second;
              });

    counts.resize(std::min(counts.size(), num_top));

    std::vector<std::pair<Label, size_t>> top_counts;
    for (const auto &encoded_pair : counts) {
        top_counts.emplace_back(label_encoder_.decode(encoded_pair.first),
                                encoded_pair.second);
    }

    return top_counts;
}

template class MultiLabelEncoded<uint64_t, std::string>;

template class LabelEncoder<std::string>;

} // namespace annotate
