#include "annotate.hpp"

#include <libmaus2/util/StringSerialisation.hpp>

#include "common/serialization.hpp"
#include "common/logger.hpp"

using libmaus2::util::StringSerialisation;


namespace annotate {

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

        return instream.good();
    } catch (...) {
        return false;
    }
}


// For each pair (first, second) in the dictionary, renames column |first|
// to |second| and merges columns with matching names, if supported.
template <typename IndexType, typename LabelType>
void MultiLabelEncoded<IndexType, LabelType>
::rename_labels(const std::unordered_map<Label, Label> &dict) {
    // old labels
    std::vector<Label> index_to_label = label_encoder_.get_labels();

    // new labels
    for (const auto &pair : dict) {
        try {
            index_to_label[label_encoder_.encode(pair.first)] = pair.second;
        } catch (const std::out_of_range &) {
            mg::common::logger->warn("Label '{}' not found"
                                     ", instruction '{} -> {}' skipped",
                                     pair.first, pair.first, pair.second);
        }
    }

    label_encoder_.clear();

    // insert new column labels
    for (const auto &label : index_to_label) {
        if (label_encoder_.label_exists(label)) {
            // no exception -> there already exists a column with this name
            mg::common::logger->error("Detected multiple labels renamed to '{}'"
                                      ". Annotation merge is not implemented"
                                      " for this annotation type.", label);
            exit(1);
        }

        // this is the first column with this name
        label_encoder_.insert_and_encode(label);
    }
}

template <class Annotator>
class IterateRowsByIndex : public Annotator::IterateRows {
  public:
    IterateRowsByIndex(const Annotator &annotator) : annotator_(annotator) {};

    typename Annotator::SetBitPositions next_row() override final {
        return annotator_.get_label_codes(i_++);
    };

  private:
    typename Annotator::Index i_ = 0;
    const Annotator &annotator_;
};

template <typename IndexType, typename LabelType>
std::unique_ptr<typename MultiLabelEncoded<IndexType, LabelType>::IterateRows>
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
    std::vector<VLabels> annotation;
    annotation.reserve(indices.size());

    for (const auto &label_codes : get_label_codes(indices)) {
        annotation.emplace_back(label_codes.size());

        for (size_t k = 0; k < label_codes.size(); ++k) {
            annotation.back()[k] = label_encoder_.decode(label_codes[k]);
        }
    }

    return annotation;
}

template <typename IndexType, typename LabelType>
std::vector<std::pair<uint64_t /* label_code */, size_t /* count */>>
MultiLabelEncoded<IndexType, LabelType>
::count_labels(const std::unordered_map<Index, size_t> &index_counts,
               size_t min_count,
               size_t count_cap) const {

    min_count = std::max(min_count, size_t(1));

    assert(count_cap >= min_count);

    size_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    std::vector<uint64_t> indices(index_counts.size());
    std::transform(index_counts.begin(), index_counts.end(), indices.begin(),
                   [](const auto &pair) { return pair.first; });

    std::vector<size_t> code_counts(this->num_labels(), 0);
    size_t max_matched = 0;
    size_t total_checked = 0;

    auto it = index_counts.begin();
    for (Index i : indices) {
        assert(it != index_counts.end());

        if (max_matched + (total_sum_count - total_checked) < min_count)
            break;

        for (size_t label_code : get_label_codes(i)) {
            assert(label_code < code_counts.size());

            code_counts[label_code] += it->second;
            max_matched = std::max(max_matched, code_counts[label_code]);
        }

        total_checked += it->second;

        ++it;
    }

    if (max_matched < min_count)
        return {};

    std::vector<std::pair<uint64_t, size_t>> label_counts;
    label_counts.reserve(code_counts.size());

    for (size_t label_code = 0; label_code < code_counts.size(); ++label_code) {
        if (code_counts[label_code] >= min_count) {
            label_counts.emplace_back(label_code,
                                      std::min(code_counts[label_code], count_cap));
        }
    }

    return label_counts;
}

// calls get_label_codes(i)
template <typename IndexType, typename LabelType>
std::vector<typename MultiLabelEncoded<IndexType, LabelType>::SetBitPositions>
MultiLabelEncoded<IndexType, LabelType>
::get_label_codes(const std::vector<Index> &indices) const {
    std::vector<SetBitPositions> rows(indices.size());

    for (size_t i = 0; i < indices.size(); ++i) {
        assert(indices[i] < this->num_objects());

        rows[i] = get_label_codes(indices[i]);
    }

    return rows;
}

template class MultiLabelEncoded<uint64_t, std::string>;

template class LabelEncoder<std::string>;

} // namespace annotate
