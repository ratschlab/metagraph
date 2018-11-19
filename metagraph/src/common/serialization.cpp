#include "serialization.hpp"

#include <cassert>
#include <map>
#include <unordered_map>

#include <sparsepp/spp.h>
#include <tsl/hopscotch_map.h>
#include <tsl/ordered_set.h>
#include <libmaus2/util/NumberSerialisation.hpp>
#include <libmaus2/util/StringSerialisation.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>

using libmaus2::util::NumberSerialisation;
using libmaus2::util::StringSerialisation;


template <typename Key>
using OrderedSet = tsl::ordered_set<Key,
                                    std::hash<Key>,
                                    std::equal_to<Key>,
                                    std::allocator<Key>,
                                    std::deque<Key, std::allocator<Key>>,
                                    std::uint64_t>;


void serialize_number(std::ostream &out, uint64_t number) {
    NumberSerialisation::serialiseNumber(out, number);
}

uint64_t load_number(std::istream &in) {
    if (!in.good()) {
        throw std::ifstream::failure("Bad stream");
    }
    return NumberSerialisation::deserialiseNumber(in);
}

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
template void serialize_number_vector<uint32_t>(std::ostream &out,
                                                const std::vector<uint32_t> &vector,
                                                size_t bits_per_number = 32);
template void serialize_number_vector<uint8_t>(std::ostream &out,
                                               const std::vector<uint8_t> &vector,
                                               size_t bits_per_number = 8);
template void serialize_number_vector<bool>(std::ostream &out,
                                            const std::vector<bool> &vector,
                                            size_t bits_per_number = 1);

uint64_t get_number_vector_size(std::istream &in) {
    if (!in.good())
        throw std::ifstream::failure("Bad stream");

    // save the position
    auto position = in.tellg();

    typename sdsl::int_vector<>::size_type size;
    typename sdsl::int_vector<>::int_width_type int_width;
    sdsl::int_vector<>::read_header(size, int_width, in);
    if (!int_width || size % int_width)
        throw std::ifstream::failure("Error when loading size and int_width in vector");

    // restore the initial position in stream
    in.seekg(position, in.beg);

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


template <class Map>
void serialize_number_number_map(std::ostream &out, const Map &map) {
    std::vector<typename Map::key_type> keys;
    keys.reserve(map.size());

    std::vector<typename Map::mapped_type> values;
    values.reserve(map.size());

    for (const auto &pair : map) {
        keys.push_back(pair.first);
        values.push_back(pair.second);
    }

    serialize_number_vector(out, keys);
    serialize_number_vector(out, values);
}

template
void serialize_number_number_map(std::ostream &out,
                                 const std::map<uint64_t, uint64_t> &map);
template
void serialize_number_number_map(std::ostream &out,
                                 const std::map<uint64_t, uint32_t> &map);
template
void serialize_number_number_map(std::ostream &out,
                                 const std::map<uint32_t, uint32_t> &map);
template
void serialize_number_number_map(std::ostream &out,
                                 const std::unordered_map<uint64_t, uint64_t> &map);
template
void serialize_number_number_map(std::ostream &out,
                                 const std::unordered_map<uint64_t, uint32_t> &map);
template
void serialize_number_number_map(std::ostream &out,
                                 const std::unordered_map<uint32_t, uint32_t> &map);

template <class Map>
void load_number_number_map(std::istream &in, Map *map) {
    assert(map);
    map->clear();

    auto keys = load_number_vector<typename Map::key_type>(in);
    auto values = load_number_vector<typename Map::mapped_type>(in);

    if (keys.size() != values.size())
        throw std::ifstream::failure("Different number of keys and values");

    for (size_t i = 0; i < keys.size(); ++i) {
        map->emplace(keys[i], values[i]);
    }
}

template
void load_number_number_map(std::istream &in,
                            std::map<uint64_t, uint64_t> *map);
template
void load_number_number_map(std::istream &in,
                            std::map<uint64_t, uint32_t> *map);
template
void load_number_number_map(std::istream &in,
                            std::map<uint32_t, uint32_t> *map);
template
void load_number_number_map(std::istream &in,
                            std::unordered_map<uint64_t, uint64_t> *map);
template
void load_number_number_map(std::istream &in,
                            std::unordered_map<uint64_t, uint32_t> *map);
template
void load_number_number_map(std::istream &in,
                            std::unordered_map<uint32_t, uint32_t> *map);


template <class Map>
void serialize_string_number_map(std::ostream &out, const Map &map) {
    serialize_number(out, map.size());

    std::vector<typename Map::mapped_type> values;
    values.reserve(map.size());

    for (const auto &pair : map) {
        StringSerialisation::serialiseString(out, pair.first);
        values.push_back(pair.second);
    }

    serialize_number_vector(out, values);
}

template
void serialize_string_number_map(std::ostream &out,
                                 const std::unordered_map<std::string, uint64_t> &map);
template
void serialize_string_number_map(std::ostream &out,
                                 const std::unordered_map<std::string, uint32_t> &map);
template
void serialize_string_number_map(std::ostream &out,
                                 const tsl::hopscotch_map<std::string, uint64_t> &map);
template
void serialize_string_number_map(std::ostream &out,
                                 const tsl::hopscotch_map<std::string, uint32_t> &map);
template
void serialize_string_number_map(std::ostream &out,
                                 const spp::sparse_hash_map<std::string, uint64_t> &map);
template
void serialize_string_number_map(std::ostream &out,
                                 const spp::sparse_hash_map<std::string, uint32_t> &map);

template <class Map>
void load_string_number_map(std::istream &in, Map *map) {
    assert(map);
    map->clear();

    auto pos = in.tellg();

    try {
        auto const size = load_number(in);

        std::vector<std::string> keys;
        keys.reserve(size);

        for (size_t i = 0; i < size; ++i) {
            keys.push_back(StringSerialisation::deserialiseString(in));
        }

        if (get_number_vector_size(in) != size)
            throw std::ifstream::failure("Different number of keys and values");

        auto values = load_number_vector<typename Map::mapped_type>(in);

        if (keys.size() != values.size())
            throw std::ifstream::failure("Different number of keys and values");

        for (size_t i = 0; i < keys.size(); ++i) {
            map->emplace(std::move(keys[i]), std::move(values[i]));
        }
    } catch (...) {
        // backward compatibility
        in.seekg(pos, in.beg);

        size_t const num = NumberSerialisation::deserialiseNumber(in);
        for (size_t i = 0; i < num; ++i) {
            std::string key = StringSerialisation::deserialiseString(in);
            auto value = NumberSerialisation::deserialiseNumber32(in);
            map->emplace(std::move(key), std::move(value));
        }
    }
}

template
void load_string_number_map(std::istream &in,
                            std::unordered_map<std::string, uint64_t> *map);
template
void load_string_number_map(std::istream &in,
                            std::unordered_map<std::string, uint32_t> *map);
template
void load_string_number_map(std::istream &in,
                            tsl::hopscotch_map<std::string, uint64_t> *map);
template
void load_string_number_map(std::istream &in,
                            tsl::hopscotch_map<std::string, uint32_t> *map);
template
void load_string_number_map(std::istream &in,
                            spp::sparse_hash_map<std::string, uint64_t> *map);
template
void load_string_number_map(std::istream &in,
                            spp::sparse_hash_map<std::string, uint32_t> *map);


template <class Set>
void serialize_set(std::ostream &out, const Set &set) {
    serialize_number(out, set.size());

    for (const auto &value : set) {
        StringSerialisation::serialiseString(out, value);
    }
}

template void serialize_set(std::ostream &out,
                            const OrderedSet<std::string> &set);

template <class Set>
void load_set(std::istream &in, Set *set) {
    assert(set);
    set->clear();

    auto const size = load_number(in);
    set->reserve(size + 1);
    set->rehash(size + 1);

    for (size_t i = 0; i < size; ++i) {
        set->insert(StringSerialisation::deserialiseString(in));
    }
}

template void load_set(std::istream &in, OrderedSet<std::string> *set);


VectorFileStream::VectorFileStream(const std::string &file)
      : istream_(MappedFile(file)) {
    if (!istream_.good())
        throw std::ifstream::failure(std::string("Bad stream file ") + file);

    values_left_ = load_number(istream_);
}

uint64_t VectorFileStream::next_value() {
    assert(values_left_ > 0);
    values_left_--;
    return load_number(istream_);
}

VectorBitStream::VectorBitStream(const bit_vector &vec,
                                 uint64_t begin,
                                 uint64_t end)
      : vector_(vec),
        begin_(begin),
        current_rank_(begin ? vector_.rank1(begin - 1) : 0) {
    assert(begin <= end);
    assert(current_rank_ <= begin);
    max_rank_ = (end == static_cast<uint64_t>(-1)
                    ? vector_.num_set_bits()
                    : (end ? vector_.rank1(end - 1) : 0));
}

uint64_t VectorBitStream::next_value() {
    assert(current_rank_ < max_rank_);
    return vector_.select1(++current_rank_) - begin_;
}
