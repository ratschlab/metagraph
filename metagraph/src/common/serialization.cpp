#include "serialization.hpp"

#include <cassert>
#include <map>
#include <string>
#include <unordered_map>

#include <tsl/hopscotch_map.h>
#include <tsl/ordered_set.h>
#include <sdsl/int_vector.hpp>


template <typename Key>
using OrderedSet = tsl::ordered_set<Key,
                                    std::hash<Key>,
                                    std::equal_to<Key>,
                                    std::allocator<Key>,
                                    std::deque<Key, std::allocator<Key>>,
                                    std::uint64_t>;


void serialize_number(std::ostream &out, uint64_t number) {
    out.put((number >> (7 * 8)) & 0xFF);
    out.put((number >> (6 * 8)) & 0xFF);
    out.put((number >> (5 * 8)) & 0xFF);
    out.put((number >> (4 * 8)) & 0xFF);
    out.put((number >> (3 * 8)) & 0xFF);
    out.put((number >> (2 * 8)) & 0xFF);
    out.put((number >> (1 * 8)) & 0xFF);
    out.put((number >> (0 * 8)) & 0xFF);

    if (!out)
        throw std::ostream::failure("Bad stream");
}

uint64_t load_number(std::istream &in) {
    if (!in.good())
        throw std::istream::failure("Bad stream");

    int const c0 = in.get();
    int const c1 = in.get();
    int const c2 = in.get();
    int const c3 = in.get();
    int const c4 = in.get();
    int const c5 = in.get();
    int const c6 = in.get();
    int const c7 = in.get();

    if (  c0 < 0 || c1 < 0 || c2 < 0 || c3 < 0 || c4 < 0 || c5 < 0 || c6 < 0 || c7 < 0 )
        throw std::istream::failure("Bad stream");

    uint64_t const u0 = c0;
    uint64_t const u1 = c1;
    uint64_t const u2 = c2;
    uint64_t const u3 = c3;
    uint64_t const u4 = c4;
    uint64_t const u5 = c5;
    uint64_t const u6 = c6;
    uint64_t const u7 = c7;

    uint64_t const u =
          (u0 << (7*8))
        | (u1 << (6*8))
        | (u2 << (5*8))
        | (u3 << (4*8))
        | (u4 << (3*8))
        | (u5 << (2*8))
        | (u6 << (1*8))
        | (u7 << (0*8))
        ;

    return u;
}

uint32_t load_number32(std::istream &in) {
    if (!in.good())
        throw std::istream::failure("Bad stream");

    int const c0 = in.get();
    int const c1 = in.get();
    int const c2 = in.get();
    int const c3 = in.get();

    if (  c0 < 0 || c1 < 0 || c2 < 0 || c3 < 0 )
        throw std::istream::failure("Bad stream");

    uint32_t const u0 = c0;
    uint32_t const u1 = c1;
    uint32_t const u2 = c2;
    uint32_t const u3 = c3;

    uint32_t const u =
          (u0 << (3*8))
        | (u1 << (2*8))
        | (u2 << (1*8))
        | (u3 << (0*8))
        ;

    return u;
}

template <typename T>
void serialize_number_vector_raw(std::ostream &out, const std::vector<T> &vector) {
    serialize_number(out, vector.size());

    for (uint64_t n : vector) {
        serialize_number(out, n);
    }
}

template void serialize_number_vector_raw<uint64_t>(std::ostream &out,
                                                    const std::vector<uint64_t> &vector);

template <typename T>
std::vector<T> load_number_vector_raw(std::istream &in) {
    assert(vector);

    const uint64_t n = load_number(in);
    std::vector<T> vector;

    for (uint64_t i = 0; i < n; ++i) {
        vector.push_back(load_number(in));
    }

    return vector;
}

template std::vector<uint64_t> load_number_vector_raw<uint64_t>(std::istream &in);

template <typename T>
void serialize_number_vector(std::ostream &out,
                             const std::vector<T> &vector,
                             size_t bits_per_number) {
    if (!out.good())
        throw std::ostream::failure("Bad stream");

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
        throw std::istream::failure("Bad stream");

    // save the position
    auto position = in.tellg();

    typename sdsl::int_vector<>::size_type size;
    typename sdsl::int_vector<>::int_width_type int_width;
    sdsl::int_vector<>::read_header(size, int_width, in);
    if (!int_width || size % int_width)
        throw std::istream::failure("Error when loading size and int_width in vector");

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


void encode_utf8(uint32_t num, std::ostream &out) {
    if (num <= 0x7F) {
        out.put(static_cast<uint8_t>(num));
    } else if (num <= 0x7FF) {
        out.put(static_cast<uint8_t>(128 | 64 | (((1<<5)-1) & (num >> 6))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num))));
    } else if (num <= 0x0000FFFF) {
        out.put(static_cast<uint8_t>(128 | 64 | 32 | (((1<<4)-1) & (num >> 12))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num >> 6))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num))));
    } else if (num <= 0x001FFFFF) {
        out.put(static_cast<uint8_t>(128 | 64 | 32 | 16 | (((1<<3)-1) & (num >> 18))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num >> 12))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num >> 6))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num))));
    } else if (num <= 0x03FFFFFF) {
        out.put(static_cast<uint8_t>(128 | 64 | 32 | 16 | 8 | (((1<<2)-1) & (num >> 24))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num >> 18))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num >> 12))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num >> 6))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num))));
    } else if (num <= 0x7FFFFFFF) {
        out.put(static_cast<uint8_t>(128 | 64 | 32 | 16 | 8 | 4 | (((1<<1)-1) & (num >> 30))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num >> 24))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num >> 18))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num >> 12))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num >> 6))));
        out.put(static_cast<uint8_t>(128 | (((1<<6)-1) & (num))));
    } else {
        throw std::runtime_error("Encoding value out of range for code.");
    }
}

uint32_t decode_utf8(std::istream &istr) {
    int len = 0;
    unsigned char mask = 0x80u;

    // compute length of utf8 repr.
    int const str0 = istr.get();

    if (str0 < 0)
        throw std::istream::failure("Bad stream");

    if ((str0 & 0xc0) == 0x80)
        throw std::istream::failure("Bad stream");

    while (str0 & mask) {
        len++;
        mask >>= 1;
    }

    // get useable bits from first byte
    unsigned int bitsinfirstbyte = 8 - len - 1;
    uint32_t number = str0 & (bitsinfirstbyte < 64 ? (static_cast<uint64_t>(1ULL)<<bitsinfirstbyte)-1 : 0xFFFFFFFFFFFFFFFFULL);

    // every additional byte provides 6 bits
    // of information
    while (--len > 0) {
        number <<= 6;
        int const strn = istr.get();

        if (strn < 0)
            throw std::istream::failure("Bad stream");

        if ((strn & 0xc0) != 0x80)
            throw std::istream::failure("Bad stream");

        number |= (strn) & 0x3f;
    }

    return number;
}


void serialize_string(std::ostream &out, const std::string_view str) {
    if (!out.good())
        throw std::ostream::failure("Bad stream");

    encode_utf8(str.size(), out);
    out.write(str.data(), str.size());

    if (!out)
        throw std::ostream::failure("Bad stream");
}

bool load_string(std::istream &in, std::string *str) {
    assert(str);

    if (!in.good()) {
        std::cerr << "Bad stream" << std::endl;
        return false;
    }

    try {
        uint64_t const u = decode_utf8(in);
        str->resize(u);
        in.read(str->data(), u);

        if (!in)
            throw std::istream::failure("Bad stream");

        return true;
    } catch (...) {
        return false;
    }
}

void serialize_string_vector(std::ostream &out, const std::vector<std::string> &vector) {
    serialize_number(out, vector.size());

    for (const std::string &s : vector) {
        serialize_string(out, s);
    }
}

bool load_string_vector(std::istream &in, std::vector<std::string> *vector) {
    assert(vector);

    try {
        uint64_t s = load_number(in);

        vector->resize(s);
        for (uint64_t i = 0; i < s; ++i) {
            load_string(in, &(*vector)[i]);
        }

        return true;
    } catch (...) {
        return false;
    }
}


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
        serialize_string(out, pair.first);
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
                keys.emplace_back();
                load_string(in, &keys.back());
            }
        } catch (...) {
            return false;
        }

        if (get_number_vector_size(in) != size)
            throw std::istream::failure("Different number of keys and values");

        std::vector<typename Map::mapped_type> values;
        if (!load_number_vector(in, &values))
            throw std::istream::failure("Failed to read map values");

        if (keys.size() != values.size())
            throw std::istream::failure("Different number of keys and values");

        for (size_t i = 0; i < keys.size(); ++i) {
            map->emplace(std::move(keys[i]), std::move(values[i]));
        }
    } catch (...) {
        // backward compatibility
        in.seekg(pos, in.beg);

        try {
            size_t const num = load_number(in);
            for (size_t i = 0; i < num; ++i) {
                std::string key;
                load_string(in, &key);
                auto value = load_number32(in);
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
        serialize_string(out, pair.second);
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
            values.emplace_back();
            load_string(in, &values.back());
        }

        if (get_number_vector_size(in) != size)
            throw std::istream::failure("Different number of keys and values");

        std::vector<typename Map::key_type> keys;
        if (!load_number_vector(in, &keys))
            throw std::istream::failure("Failed to read map keys");

        if (keys.size() != values.size())
            throw std::istream::failure("Different number of keys and values");

        for (size_t i = 0; i < values.size(); ++i) {
            map->emplace(std::move(keys[i]), std::move(values[i]));
        }
    } catch (...) {
        // backward compatibility
        in.seekg(pos, in.beg);

        try {
            size_t const num = load_number(in);
            for (size_t i = 0; i < num; ++i) {
                auto key = load_number32(in);
                std::string value;
                load_string(in, &value);
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
        serialize_string(out, value);
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
            std::string val;
            load_string(in, &val);
            set->insert(std::move(val));
        }
    } catch (...) {
        return false;
    }

    return true;
}

template bool load_set(std::istream &in, OrderedSet<std::string> *set);
