//
// Created by Jan Studený on 2019-06-01.
//

#ifndef METAGRAPH_DENSE_HASHMAP_HPP
#define METAGRAPH_DENSE_HASHMAP_HPP

#include <vector>
#include <sdsl/bit_vectors.hpp>



template <typename T,typename bit_vector=sdsl::bit_vector, typename rank_support=sdsl::rank_support_v<1>>
class DenseHashMap {
public:
    DenseHashMap() = default;


    DenseHashMap(bit_vector *isElement,rank_support *rank) {
        init(isElement,rank);
    }

    void init(bit_vector *isElement, rank_support *rank) {
        this->is_element = isElement;
        this->rank = rank;
        elements = decltype(elements)(rank->rank(rank->size()));
    }

    T& operator[](int64_t n) {
        assert(n < is_element->size());
        assert((*is_element)[n]);
        return elements[rank->rank(n)];
    }

    const T& at(int64_t n) const {
        assert(n < is_element->size());
        assert((*is_element)[n]);
        return elements.at(rank->rank(n));
    }

    int count(int64_t n) const {
        assert(n < is_element->size());
        return (*is_element)[n];
    }

    bit_vector * is_element;
    rank_support * rank;
    std::vector<T> elements;
};

#endif //METAGRAPH_DENSE_HASHMAP_HPP
