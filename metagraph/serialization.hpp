#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

#include <string>
#include <set>
#include <unordered_map>

#include <libmaus2/util/NumberSerialisation.hpp>
#include <libmaus2/util/StringSerialisation.hpp>

void serialize_annotation_map(std::ostream & out, std::unordered_map<std::uint16_t, std::set<std::string> > const & M ) {
    libmaus2::util::NumberSerialisation::serialiseNumber(out, M.size()); 

    for (std::unordered_map<std::uint16_t, std::set<std::string> >::const_iterator it = M.begin(); it != M.end(); ++it) {
        libmaus2::util::NumberSerialisation::serialiseNumberShort(out, it->first); 
        libmaus2::util::StringSerialisation::serialiseStringSet(out, it->second);
    }
}


std::unordered_map<std::uint16_t, std::set<std::string> > deserialize_annotation_map(std::istream & in) {
    size_t const num = libmaus2::util::NumberSerialisation::deserialiseNumber(in); 

    std::unordered_map<std::uint16_t, std::set<std::string> > M;
    for (size_t i = 0; i < num; ++i) {
        uint16_t key = libmaus2::util::NumberSerialisation::deserialiseNumberShort(in); 
        std::set<std::string> value = libmaus2::util::StringSerialisation::deserialiseStringSet(in);
        M[key] = value;
    }
    return M;
}

#endif
