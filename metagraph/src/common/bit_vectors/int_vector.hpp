#ifndef __INT_VECTOR_HPP__
#define __INT_VECTOR_HPP__

#include <functional>

#include <sdsl/int_vector.hpp>
#include "bit_vector.hpp"


void call_nonzeros(const sdsl::int_vector<> &vector,
                   uint64_t begin, uint64_t end,
                   std::function<void(uint64_t /* index */,
                                      uint64_t /* value */)> callback);

void call_nonzeros(const sdsl::int_vector<> &vector,
                   std::function<void(uint64_t /* index */,
                                      uint64_t /* value */)> callback);

template <class Vector>
void insert_new_indexes(Vector &vector, bit_vector_dyn *new_indexes);

#endif // __INT_VECTOR_HPP__
