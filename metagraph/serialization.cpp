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
    if (!out.good()) {
        throw std::ofstream::failure("Bad stream");
    }
    sdsl::int_vector<> int_vector(vector.size(), 0, bits_per_number);
    for (size_t i = 0; i < vector.size(); ++i) {
        int_vector[i] = vector[i];
    }
    int_vector.serialize(out);
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

uint64_t load_number_vector_size(std::istream &in) {
    if (!in.good()) {
        throw std::ifstream::failure("Bad stream");
    }
    typename sdsl::int_vector<>::size_type size;
    typename sdsl::int_vector<>::int_width_type int_width;
    sdsl::int_vector<>::read_header(size, int_width, in);
    return size / int_width;
}

template <typename T>
std::vector<T> load_number_vector(std::istream &in) {
    if (!in.good()) {
        throw std::ifstream::failure("Bad stream");
    }
    try {
        sdsl::int_vector<> int_vector;
        int_vector.load(in);
        return std::vector<T>(int_vector.begin(), int_vector.end());
    } catch (...) {
        throw std::ifstream::failure("Bad stream");
    }
}

template std::vector<uint64_t> load_number_vector<uint64_t>(std::istream &in);
template std::vector<uint32_t> load_number_vector<uint32_t>(std::istream &in);
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
