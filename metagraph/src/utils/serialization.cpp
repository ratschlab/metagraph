#include "serialization.hpp"

#include <cassert>
#include <map>
#include <unordered_map>

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
    if (!in.good())
        throw std::ifstream::failure("Bad stream");

    return NumberSerialisation::deserialiseNumber(in);
}

template <typename T>
void serialize_number_vector(std::ostream &out,
                             const std::vector<T> &vector,
                             size_t bits_per_number) {
    if (!out.good())
        throw std::ofstream::failure("Bad stream");

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
bool load_number_vector(std::istream &in, std::vector<T> *vector) {
    assert(vector);

    if (!in.good()) {
        std::cerr << "Bad stream" << std::endl;
        return false;
    }

    try {
        sdsl::int_vector<> int_vector;
        int_vector.load(in);
        vector->assign(int_vector.begin(), int_vector.end());
        return true;
    } catch (...) {
        return false;
    }
}

template bool load_number_vector<uint64_t>(std::istream &in, std::vector<uint64_t> *);
template bool load_number_vector<uint32_t>(std::istream &in, std::vector<uint32_t> *);
template bool load_number_vector<uint8_t>(std::istream &in, std::vector<uint8_t> *);
template bool load_number_vector<bool>(std::istream &in, std::vector<bool> *);


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
bool load_number_number_map(std::istream &in, Map *map) {
    assert(map);
    map->clear();

    std::vector<typename Map::key_type> keys;
    if (!load_number_vector(in, &keys)) {
        std::cerr << "Failed to load map keys" << std::endl;
        return false;
    }

    std::vector<typename Map::mapped_type> values;
    if (!load_number_vector(in, &values)) {
        std::cerr << "Failed to load map values" << std::endl;
        return false;
    }

    if (keys.size() != values.size()) {
        std::cerr << "Different number of keys and values" << std::endl;
        return false;
    }

    for (size_t i = 0; i < keys.size(); ++i) {
        map->emplace(keys[i], values[i]);
    }

    return true;
}

template
bool load_number_number_map(std::istream &in,
                            std::map<uint64_t, uint64_t> *map);
template
bool load_number_number_map(std::istream &in,
                            std::map<uint64_t, uint32_t> *map);
template
bool load_number_number_map(std::istream &in,
                            std::map<uint32_t, uint32_t> *map);
template
bool load_number_number_map(std::istream &in,
                            std::unordered_map<uint64_t, uint64_t> *map);
template
bool load_number_number_map(std::istream &in,
                            std::unordered_map<uint64_t, uint32_t> *map);
template
bool load_number_number_map(std::istream &in,
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

template <class Map>
bool load_string_number_map(std::istream &in, Map *map) {
    assert(map);
    map->clear();

    // save for the backwards compatibility case
    auto pos = in.tellg();

    try {
        auto const size = load_number(in);

        std::vector<std::string> keys;
        keys.reserve(size);

        try {
            for (size_t i = 0; i < size; ++i) {
                keys.push_back(StringSerialisation::deserialiseString(in));
            }
        } catch (...) {
            return false;
        }

        if (get_number_vector_size(in) != size)
            throw std::ifstream::failure("Different number of keys and values");

        std::vector<typename Map::mapped_type> values;
        if (!load_number_vector(in, &values))
            throw std::ifstream::failure("Failed to read map values");

        if (keys.size() != values.size())
            throw std::ifstream::failure("Different number of keys and values");

        for (size_t i = 0; i < keys.size(); ++i) {
            map->emplace(std::move(keys[i]), std::move(values[i]));
        }
    } catch (...) {
        // backward compatibility
        in.seekg(pos, in.beg);

        try {
            size_t const num = NumberSerialisation::deserialiseNumber(in);
            for (size_t i = 0; i < num; ++i) {
                std::string key = StringSerialisation::deserialiseString(in);
                auto value = NumberSerialisation::deserialiseNumber32(in);
                map->emplace(std::move(key), std::move(value));
            }
        } catch (...) {
            return false;
        }
    }

    return true;
}

template
bool load_string_number_map(std::istream &in,
                            std::unordered_map<std::string, uint64_t> *map);
template
bool load_string_number_map(std::istream &in,
                            std::unordered_map<std::string, uint32_t> *map);
template
bool load_string_number_map(std::istream &in,
                            tsl::hopscotch_map<std::string, uint64_t> *map);
template
bool load_string_number_map(std::istream &in,
                            tsl::hopscotch_map<std::string, uint32_t> *map);

template <class Map>
void serialize_number_string_map(std::ostream &out, const Map &map) {
    serialize_number(out, map.size());

    std::vector<typename Map::key_type> keys;
    keys.reserve(map.size());

    for (const auto &pair : map) {
        keys.push_back(pair.first);
        StringSerialisation::serialiseString(out, pair.second);
    }

    serialize_number_vector(out, keys);
}

template
void serialize_number_string_map(std::ostream &out,
                                 const std::unordered_map<uint64_t, std::string> &map);
template
void serialize_number_string_map(std::ostream &out,
                                 const std::unordered_map<uint32_t, std::string> &map);
template
void serialize_number_string_map(std::ostream &out,
                                 const tsl::hopscotch_map<uint64_t, std::string> &map);
template
void serialize_number_string_map(std::ostream &out,
                                 const tsl::hopscotch_map<uint32_t, std::string> &map);

template <class Map>
bool load_number_string_map(std::istream &in, Map *map) {
    assert(map);
    map->clear();

    // save for the backwards compatibility case
    auto pos = in.tellg();

    try {
        auto const size = load_number(in);

        std::vector<std::string> values;
        values.reserve(size);

        for (size_t i = 0; i < size; ++i) {
            values.push_back(StringSerialisation::deserialiseString(in));
        }

        if (get_number_vector_size(in) != size)
            throw std::ifstream::failure("Different number of keys and values");

        std::vector<typename Map::key_type> keys;
        if (!load_number_vector(in, &keys))
            throw std::ifstream::failure("Failed to read map keys");

        if (keys.size() != values.size())
            throw std::ifstream::failure("Different number of keys and values");

        for (size_t i = 0; i < values.size(); ++i) {
            map->emplace(std::move(keys[i]), std::move(values[i]));
        }
    } catch (...) {
        // backward compatibility
        in.seekg(pos, in.beg);

        try {
            size_t const num = NumberSerialisation::deserialiseNumber(in);
            for (size_t i = 0; i < num; ++i) {
                auto key = NumberSerialisation::deserialiseNumber32(in);
                std::string value = StringSerialisation::deserialiseString(in);
                map->emplace(std::move(key), std::move(value));
            }
        } catch (...) {
            return false;
        }
    }

    return true;
}

template
bool load_number_string_map(std::istream &in,
                            std::unordered_map<uint64_t, std::string> *map);
template
bool load_number_string_map(std::istream &in,
                            std::unordered_map<uint32_t, std::string> *map);
template
bool load_number_string_map(std::istream &in,
                            tsl::hopscotch_map<uint64_t, std::string> *map);
template
bool load_number_string_map(std::istream &in,
                            tsl::hopscotch_map<uint32_t, std::string> *map);


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
bool load_set(std::istream &in, Set *set) {
    assert(set);
    set->clear();

    auto const size = load_number(in);
    set->reserve(size + 1);
    set->rehash(size + 1);

    try {
        for (size_t i = 0; i < size; ++i) {
            set->insert(StringSerialisation::deserialiseString(in));
        }
    } catch (...) {
        return false;
    }

    return true;
}

template bool load_set(std::istream &in, OrderedSet<std::string> *set);


BitVectorFileInStream::BitVectorFileInStream(const std::string &file)
      : istream_(file, std::ios::binary) {
    if (!istream_.good())
        throw std::ifstream::failure(std::string("Bad stream file ") + file);

    istream_ >> length_;
    istream_ >> values_left_;
}

uint64_t BitVectorFileInStream::next_value() {
    assert(values_left_ > 0);
    --values_left_;

    uint64_t next_val;
    istream_ >> next_val;

    if (next_val >= length_)
        throw std::runtime_error("Error: index >= length: "
            + std::to_string(next_val)
            + " >= " + std::to_string(length_));

    return next_val;
}

VectorBitInStream::VectorBitInStream(const bit_vector &vec,
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

uint64_t VectorBitInStream::next_value() {
    assert(current_rank_ < max_rank_);
    return vector_.select1(++current_rank_) - begin_;
}


VectorFileOutStream::VectorFileOutStream(const std::string &file)
      : ostream_(file) {
    if (!ostream_.good())
        throw std::ofstream::failure(std::string("Bad stream file ") + file);
}

void VectorFileOutStream::write_value(uint64_t value) {
    ostream_ << value << std::endl;
}
