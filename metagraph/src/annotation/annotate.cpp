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

template <typename IndexType, typename LabelType>
void MultiLabelEncoded<IndexType, LabelType>
::call_rows(const std::vector<Index> &indices,
            std::function<void(std::vector<uint64_t>&&)> row_callback,
            std::function<bool()> terminate) const {
    for (Index i : indices) {
        if (terminate())
            break;

        row_callback(get_label_indices(i));
    }
}

// For each index i in indices, check of i has the label. Return
// true if the finished callback evaluates true during execution.
template <typename IndexType, typename LabelType>
bool MultiLabelEncoded<IndexType, LabelType>
::call_indices_until(const std::vector<Index> &indices,
                     const Label &label,
                     std::function<void(Index)> index_callback,
                     std::function<bool()> finished) const {
    for (Index i : indices) {
        if (finished())
            return true;

        if (this->has_label(i, label))
            index_callback(i);
    }

    return finished();
}

template <typename Annotator>
class IterateRowsByIndex : public IterateRows {
  public:
    IterateRowsByIndex(const Annotator &annotator)
          : annotator_(annotator) {};

    std::vector<uint64_t> next_row() override final {
        return annotator_.get_label_indices(i_++);
    };

  private:
    typename Annotator::Index i_ = 0;
    const Annotator &annotator_;
};

template <typename IndexType, typename LabelType>
std::unique_ptr<IterateRows>
MultiLabelEncoded<IndexType, LabelType>::iterator() const {
    return std::make_unique<IterateRowsByIndex<MultiLabelEncoded<IndexType, LabelType>>>(*this);
};

template class MultiLabelEncoded<uint64_t, std::string>;

template class LabelEncoder<std::string>;

} // namespace annotate
