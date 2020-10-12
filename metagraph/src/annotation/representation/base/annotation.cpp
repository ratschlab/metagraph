#include "annotation.hpp"

#include <libmaus2/util/StringSerialisation.hpp>

#include "common/serialization.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace annot {

using libmaus2::util::StringSerialisation;
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

template <typename LabelType>
std::vector<std::pair<uint64_t /* label_code */, size_t /* count */>>
MultiLabelEncoded<LabelType>
::count_labels(const std::vector<std::pair<Index, size_t>> &index_counts,
               size_t min_count,
               size_t count_cap) const {
    assert(count_cap >= min_count);

    if (!count_cap)
        return {};

    min_count = std::max(min_count, size_t(1));

    size_t total_sum_count = 0;
    for (const auto &pair : index_counts) {
        total_sum_count += pair.second;
    }

    if (total_sum_count < min_count)
        return {};

    std::vector<size_t> code_counts(this->num_labels(), 0);
    size_t max_matched = 0;
    size_t total_checked = 0;

    for (auto [i, count] : index_counts) {
        if (max_matched + (total_sum_count - total_checked) < min_count)
            break;

        for (size_t label_code : get_matrix().get_row(i)) {
            assert(label_code < code_counts.size());

            code_counts[label_code] += count;
            max_matched = std::max(max_matched, code_counts[label_code]);
        }

        total_checked += count;
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

template <typename IndexType, typename LabelType>
void MultiLabelAnnotation<IndexType, LabelType>
::add_label_counts(const std::vector<Index> &indices,
                   const VLabels &labels,
                   const std::vector<uint32_t> &counts) {
    std::ignore = indices;
    std::ignore = labels;
    std::ignore = counts;
    logger->error("Adding label counts is not implemented for this annotator");
    exit(1);
}

template class MultiLabelEncoded<std::string>;

template class LabelEncoder<std::string>;

} // namespace annot
} // namespace mtg
