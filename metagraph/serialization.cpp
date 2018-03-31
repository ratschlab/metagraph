#include "serialization.hpp"

#include <cassert>

#include <libmaus2/util/NumberSerialisation.hpp>
#include <libmaus2/util/StringSerialisation.hpp>
#include <sdsl/int_vector.hpp>

using libmaus2::util::NumberSerialisation;


template <typename T>
void serialize_number_vector(std::ostream &out,
                             const std::vector<T> &vector,
                             size_t bits_per_number) {
    std::ignore = bits_per_number;
    NumberSerialisation::serialiseNumberVector(out, vector);
}

template void serialize_number_vector<uint64_t>(std::ostream &out,
                                                const std::vector<uint64_t> &vector,
                                                size_t bits_per_number = 64);
template void serialize_number_vector<uint8_t>(std::ostream &out,
                                               const std::vector<uint8_t> &vector,
                                               size_t bits_per_number = 8);
template void serialize_number_vector<bool>(std::ostream &out,
                                            const std::vector<bool> &vector,
                                            size_t bits_per_number = 1);


template <typename T>
std::vector<T> load_number_vector(std::istream &in) {
    return NumberSerialisation::deserialiseNumberVector<T>(in);
}

template std::vector<uint64_t> load_number_vector<uint64_t>(std::istream &in);
template std::vector<uint8_t> load_number_vector<uint8_t>(std::istream &in);
template std::vector<bool> load_number_vector<bool>(std::istream &in);


void serialize_number_number_map(std::ostream &out, const std::map<uint32_t, uint32_t> &M) {
    NumberSerialisation::serialiseNumber(out, M.size());
    for (auto it = M.begin(); it != M.end(); ++it) {
        NumberSerialisation::serialiseNumber32(out, it->first);
        NumberSerialisation::serialiseNumber32(out, it->second);
    }
}

std::map<std::uint32_t, uint32_t> load_number_number_map(std::istream &in) {
    std::map<std::uint32_t, uint32_t> M;
    size_t const num = NumberSerialisation::deserialiseNumber(in);
    for (size_t i = 0; i < num; ++i) {
        uint32_t key = NumberSerialisation::deserialiseNumber32(in);
        uint32_t value = NumberSerialisation::deserialiseNumber32(in);
        M.insert(std::make_pair(key, value));
    }
    return M;
};


void serialize_string_number_map(std::ostream &out, const std::unordered_map<std::string, uint32_t> &M) {
    NumberSerialisation::serialiseNumber(out, M.size());
    for (auto it = M.begin(); it != M.end(); ++it) {
        libmaus2::util::StringSerialisation::serialiseString(out, it->first);
        NumberSerialisation::serialiseNumber32(out, it->second);
    }
}

std::unordered_map<std::string, uint32_t> load_string_number_map(std::istream &in) {
    std::unordered_map<std::string, uint32_t> M;
    size_t const num = NumberSerialisation::deserialiseNumber(in);
    for (size_t i = 0; i < num; ++i) {
        std::string key = libmaus2::util::StringSerialisation::deserialiseString(in);
        uint32_t value = NumberSerialisation::deserialiseNumber32(in);
        M.insert(std::make_pair(key, value));
    }
    return M;
}
