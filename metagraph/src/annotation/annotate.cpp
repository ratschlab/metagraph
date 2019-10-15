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
::rename_labels(const std::unordered_map<Label, Label> &dict) {
    std::vector<Label> index_to_label(label_encoder_.size());
    // old labels
    for (size_t i = 0; i < index_to_label.size(); ++i) {
        index_to_label[i] = label_encoder_.decode(i);
    }
    // new labels
    for (const auto &pair : dict) {
        try {
            index_to_label[label_encoder_.encode(pair.first)] = pair.second;
        } catch (const std::runtime_error&) {
            std::cerr << "Warning: label '" << pair.first << "' not"
                      << " found in annotation. Skipping instruction"
                      << " '" << pair.first << " -> " << pair.second << "'."
                      << std::endl;
        }
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

template <typename IndexType, typename LabelType>
void MultiLabelEncoded<IndexType, LabelType>
::call_rows(const std::vector<Index> &indices,
            const std::function<void(std::vector<uint64_t>&&)> &row_callback,
            const std::function<bool()> &terminate) const {
    for (Index i : indices) {
        if (terminate())
            break;

        row_callback(get_label_codes(i));
    }
}

template <typename Annotator>
class IterateRowsByIndex : public IterateRows {
  public:
    IterateRowsByIndex(const Annotator &annotator)
          : annotator_(annotator) {};

    std::vector<uint64_t> next_row() override final {
        return annotator_.get_label_codes(i_++);
    };

  private:
    typename Annotator::Index i_ = 0;
    const Annotator &annotator_;
};

template <typename IndexType, typename LabelType>
std::unique_ptr<IterateRows>
MultiLabelEncoded<IndexType, LabelType>::iterator() const {
    return std::make_unique<IterateRowsByIndex<MultiLabelEncoded<IndexType, LabelType>>>(*this);
}

// calls get_label_codes(i)
template <typename IndexType, typename LabelType>
typename MultiLabelEncoded<IndexType, LabelType>::VLabels
MultiLabelEncoded<IndexType, LabelType>::get_labels(Index i) const {
    assert(i < this->num_objects());

    const auto &label_codes = get_label_codes(i);

    VLabels labels(label_codes.size());

    for (size_t j = 0; j < label_codes.size(); ++j) {
        labels[j] = label_encoder_.decode(label_codes[j]);
    }

    return labels;
}

// calls get_label_codes(indices)
template <typename IndexType, typename LabelType>
std::vector<typename MultiLabelEncoded<IndexType, LabelType>::VLabels>
MultiLabelEncoded<IndexType, LabelType>
::get_labels(const std::vector<Index> &indices) const {
    auto rows = get_label_codes(indices);

    std::vector<VLabels> annotation(rows.size());

    for (size_t i = 0; i < rows.size(); ++i) {
        auto row = std::move(rows[i]);

        annotation[i].reserve(row.size());
        for (auto label_code : row) {
            annotation[i].push_back(label_encoder_.decode(label_code));
        }
    }

    return annotation;
}

// calls get_label_codes(i)
template <typename IndexType, typename LabelType>
std::vector<std::vector<uint64_t>>
MultiLabelEncoded<IndexType, LabelType>
::get_label_codes(const std::vector<Index> &indices) const {
    std::vector<std::vector<uint64_t>> rows(indices.size());

    for (size_t i = 0; i < indices.size(); ++i) {
        assert(indices[i] < this->num_objects());

        rows[i] = get_label_codes(indices[i]);
    }

    return rows;
}

template class MultiLabelEncoded<uint64_t, std::string>;

template class LabelEncoder<std::string>;

} // namespace annotate
