#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

#include <assert.h>
#include <string>
#include <set>
#include <unordered_map>
#include <map>
#include <cmath>

#include <libmaus2/util/NumberSerialisation.hpp>
#include <libmaus2/util/StringSerialisation.hpp>


void serialize_number_number_map(std::ostream &out, const std::map<uint32_t, uint32_t> &M) {
    libmaus2::util::NumberSerialisation::serialiseNumber(out, M.size());
    for (auto it = M.begin(); it != M.end(); ++it) {
        libmaus2::util::NumberSerialisation::serialiseNumber32(out, it->first);
        libmaus2::util::NumberSerialisation::serialiseNumber32(out, it->second);
    }
}

std::map<std::uint32_t, uint32_t> load_number_number_map(std::istream &in) {
    std::map<std::uint32_t, uint32_t> M;
    size_t const num = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
    for (size_t i = 0; i < num; ++i) {
        uint32_t key = libmaus2::util::NumberSerialisation::deserialiseNumber32(in);
        uint32_t value = libmaus2::util::NumberSerialisation::deserialiseNumber32(in);
        M.insert(std::make_pair(key, value));
    }
    return M;
};

void serialize_number_vector(std::ostream &out, const std::vector<std::uint32_t> &V) {
    libmaus2::util::NumberSerialisation::serialiseNumber32Vector(out, V);
}

std::vector<uint32_t> load_number_vector(std::istream &in) {
    return libmaus2::util::NumberSerialisation::deserialiseNumber32Vector<uint32_t>(in);
}

void serialize_string_number_map(std::ostream &out, const std::unordered_map<std::string, uint32_t> &M) {
    libmaus2::util::NumberSerialisation::serialiseNumber(out, M.size());
    for (auto it = M.begin(); it != M.end(); ++it) {
        libmaus2::util::StringSerialisation::serialiseString(out, it->first);
        libmaus2::util::NumberSerialisation::serialiseNumber32(out, it->second);
    }
}

std::unordered_map<std::string, uint32_t> load_string_number_map(std::istream &in) {
    std::unordered_map<std::string, uint32_t> M;
    size_t const num = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
    for (size_t i = 0; i < num; ++i) {
        std::string key = libmaus2::util::StringSerialisation::deserialiseString(in);
        uint32_t value = libmaus2::util::NumberSerialisation::deserialiseNumber32(in);
        M.insert(std::make_pair(key, value));
    }
    return M;
}

void serialize_strings_vector(std::ostream &out, const std::vector<std::string> &V) {
    libmaus2::util::StringSerialisation::serialiseStringVector(out, V);
}

std::vector<std::string> load_strings_vector(std::istream &in) {
    return libmaus2::util::StringSerialisation::deserialiseStringVector(in);
}

#endif // __SERIALIZATION_HPP__
