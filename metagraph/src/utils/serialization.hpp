#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

#include <fstream>
#include <vector>


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


#endif // __SERIALIZATION_HPP__
