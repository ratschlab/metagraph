#ifndef __SERIALIZATION_HPP__
#define __SERIALIZATION_HPP__

#include <string>
#include <vector>
#include <map>
#include <unordered_map>


void serialize_number_number_map(std::ostream &out,
                                 const std::map<uint32_t, uint32_t> &M);

std::map<std::uint32_t, uint32_t> load_number_number_map(std::istream &in);

void serialize_string_number_map(std::ostream &out,
                                 const std::unordered_map<std::string, uint32_t> &M);

std::unordered_map<std::string, uint32_t> load_string_number_map(std::istream &in);


#endif // __SERIALIZATION_HPP__
