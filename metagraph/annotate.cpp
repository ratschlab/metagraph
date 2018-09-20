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
        encode_label_ = load_string_number_map(instream);
        decode_label_ = StringSerialisation::deserialiseStringVector(instream);

        return true;
    } catch (...) {
        return false;
    }
}

template class LabelEncoder<std::string>;

} // namespace annotate
