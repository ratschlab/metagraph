#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

#include <fstream>
#include <string>
#include <string_view>
#include <vector>


void serialize_number(std::ostream &out, uint64_t number);

uint64_t load_number(std::istream &in);

uint32_t load_number32(std::istream &in);

template <typename T>
void serialize_number_vector_raw(std::ostream &out, const std::vector<T> &vector);

template <typename T>
std::vector<T> load_number_vector_raw(std::istream &in);


template <typename T>
void serialize_number_vector(std::ostream &out,
                             const std::vector<T> &vector,
                             size_t bits_per_number = sizeof(T) * 8);

uint64_t get_number_vector_size(std::istream &in);

template <typename T>
bool load_number_vector(std::istream &in, std::vector<T> *vector);


void serialize_string(std::ostream &out, const std::string_view str);
bool load_string(std::istream &in, std::string *str);

void serialize_string_vector(std::ostream &out, const std::vector<std::string> &vector);
bool load_string_vector(std::istream &in, std::vector<std::string> *vector);


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

class Serializer {
  public:
    explicit Serializer(std::ostream &os) : os_(os) {}

    template <class T>
    void operator()(const T &value) {
        os_.write(reinterpret_cast<const char *>(&value), sizeof(T));
    }

  private:
    std::ostream &os_;
};

class Deserializer {
  public:
    explicit Deserializer(std::istream &is) : is_(is) {}

    template <class T>
    T operator()() {
        T value;
        is_.read(reinterpret_cast<char *>(&value), sizeof(T));
        return value;
    }

  private:
    std::istream &is_;
};

#endif // __SERIALIZATION_HPP__
