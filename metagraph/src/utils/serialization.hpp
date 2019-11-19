#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

#include <fstream>
#include <string>
#include <vector>

#include "utils/bit_vectors/bit_vector.hpp"


void serialize_number(std::ostream &out, uint64_t number);

uint64_t load_number(std::istream &in);


template <typename T>
void serialize_number_vector(std::ostream &out,
                             const std::vector<T> &vector,
                             size_t bits_per_number = sizeof(T) * 8);

uint64_t get_number_vector_size(std::istream &in);

template <typename T>
bool load_number_vector(std::istream &in, std::vector<T> *vector);


template <class Map>
void serialize_number_number_map(std::ostream &out, const Map &map);

template <class Map>
bool load_number_number_map(std::istream &in, Map *map);


template <class Map>
void serialize_string_number_map(std::ostream &out, const Map &M);

template <class Map>
bool load_string_number_map(std::istream &in, Map *map);

template <class Map>
void serialize_number_string_map(std::ostream &out, const Map &M);

template <class Map>
bool load_number_string_map(std::istream &in, Map *map);


template <class Set>
void serialize_set(std::ostream &out, const Set &set);

template <class Set>
bool load_set(std::istream &in, Set *set);


/** Vector streaming
 *
 * Used to extract the indices of set bits from either a file, or compressed data structures
 *
 */

// Read vector from a source
class VectorStream {
  public:
    virtual ~VectorStream() {}
    virtual uint64_t next_value() = 0;
    virtual uint64_t values_left() const = 0;
};

// Read numbers from file
class VectorFileStream : public VectorStream {
  public:
    VectorFileStream(const std::string &file);

    uint64_t next_value();
    uint64_t values_left() const { return values_left_; }

  private:
    std::ifstream istream_;
    uint64_t length_;
    uint64_t values_left_;
};

// Return set bits from a bit vector
class VectorBitStream : public VectorStream {
  public:
    VectorBitStream(const bit_vector &vector,
                    uint64_t begin = 0,
                    uint64_t end = static_cast<uint64_t>(-1));

    uint64_t next_value();
    uint64_t values_left() const { return max_rank_ - current_rank_; }

  private:
    const bit_vector &vector_;
    uint64_t begin_;
    uint64_t current_rank_;
    uint64_t max_rank_;
};

#endif // __SERIALIZATION_HPP__
