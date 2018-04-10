#include "annotate.hpp"

#include <libmaus2/util/StringSerialisation.hpp>

#include "serialization.hpp"

using libmaus2::util::StringSerialisation;


namespace annotate {

size_t StringEncoder::encode(const Color &color,
                             bool insert_if_not_exists) {
    auto it = encode_color_.find(color);
    if (it != encode_color_.end())
        return it->second;

    if (!insert_if_not_exists)
        throw std::runtime_error("ERROR: No such color");

    encode_color_[color] = decode_color_.size();
    decode_color_.push_back(color);

    return decode_color_.size() - 1;
}

const StringEncoder::Color& StringEncoder::decode(size_t code) const {
    return decode_color_.at(code);
}

void StringEncoder::serialize(std::ostream &outstream) const {
    serialize_string_number_map(outstream, encode_color_);
    StringSerialisation::serialiseStringVector(outstream, decode_color_);
}

bool StringEncoder::load(std::istream &instream) {
    if (!instream.good())
        return false;

    try {
        encode_color_ = load_string_number_map(instream);
        decode_color_ = StringSerialisation::deserialiseStringVector(instream);

        return true;
    } catch (...) {
        return false;
    }
}

} // namespace annotate
