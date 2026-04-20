#include "annotation.hpp"

#include "common/serialization.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace annot {

using mtg::common::logger;

template <typename Label>
size_t LabelEncoder<Label>::insert_and_encode(const Label &label) {
    // if it's a static encoder, copy the data to a local mutable object
    if (remote_data_) {
        encode_label_ = *remote_data_;
        remote_data_ = nullptr;
    }
    auto it = encode_label_.emplace(label).first;
    return it - encode_label_.begin();
}

template <typename Label>
size_t LabelEncoder<Label>::encode(const Label &label) const {
    auto it = data().find(label);
    if (it == data().end())
        throw std::out_of_range("Label not found");
    return it - data().begin();
}

template<>
void LabelEncoder<std::string>::serialize(std::ostream &outstream) const {
    outstream.write("LE-v2.0", 7);
    Serializer serializer(outstream);
    data().serialize(serializer);
}

template<typename Label>
void LabelEncoder<Label>::merge(const LabelEncoder<Label> &other) {
    for (size_t i = 0; i < other.size(); ++i) {
        insert_and_encode(other.decode(i));
    }
}

template<>
bool LabelEncoder<std::string>::load(std::istream &instream) {
    clear();

    // DEBUG: trace load path — remove once corrupt-annotation test flake is diagnosed.
    logger->warn("[DEBUG LabelEncoder::load] enter, good={}", instream.good());
    if (!instream.good())
        return false;

    try {
        auto pos = instream.tellg();
        std::string version(7, '\0');
        bool read_ok = (bool)instream.read(version.data(), 7);
        std::string hex;
        for (unsigned char c : version) {
            char buf[4]; std::snprintf(buf, sizeof(buf), "%02x ", c);
            hex += buf;
        }
        logger->warn("[DEBUG LabelEncoder::load] read_header read_ok={} gcount={} hex={}",
                     read_ok, (long)instream.gcount(), hex);
        if (read_ok && version == "LE-v2.0") {
            logger->warn("[DEBUG LabelEncoder::load] v2 branch");
            Deserializer deserializer(instream);
            encode_label_ = VectorSet<std::string>::deserialize(deserializer, true);
            logger->warn("[DEBUG LabelEncoder::load] v2 done good={} size={}",
                         instream.good(), encode_label_.size());
            return instream.good();
        }
        logger->warn("[DEBUG LabelEncoder::load] legacy branch");
        // backward compatibility
        if (!instream.seekg(pos)) {
            logger->error("Couldn't seek in the input stream when reading label encoder");
            return false;
        }
        {
            std::vector<std::string> labels;
            if (!load_string_vector(instream, &labels)) {
                logger->warn("[DEBUG LabelEncoder::load] legacy first load_string_vector failed");
                return false;
            }
            logger->warn("[DEBUG LabelEncoder::load] legacy 1st string_vector size={}", labels.size());
            std::vector<uint64_t> values;
            if (!load_number_vector(instream, &values)) {
                logger->warn("[DEBUG LabelEncoder::load] legacy load_number_vector failed");
                return false;
            }
            logger->warn("[DEBUG LabelEncoder::load] legacy number_vector size={}", values.size());
        }
        std::vector<std::string> labels;
        if (!load_string_vector(instream, &labels)) {
            logger->warn("[DEBUG LabelEncoder::load] legacy 2nd load_string_vector failed");
            return false;
        }
        logger->warn("[DEBUG LabelEncoder::load] legacy 2nd string_vector size={} final good={}",
                     labels.size(), instream.good());
        encode_label_ = VectorSet<std::string>(std::make_move_iterator(labels.begin()),
                                               std::make_move_iterator(labels.end()));
        return instream.good();
    } catch (const std::exception &e) {
        logger->warn("[DEBUG LabelEncoder::load] std::exception: {}", e.what());
        return false;
    } catch (...) {
        logger->warn("[DEBUG LabelEncoder::load] unknown exception");
        return false;
    }
}

template<typename Label>
LabelEncoder<Label> LabelEncoder<Label>::make_static_copy() const {
    LabelEncoder<Label> static_copy;
    static_copy.remote_data_ = &data();
    return static_copy;
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
