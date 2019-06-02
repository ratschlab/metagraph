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

//    void operator=(DenseHashMap<T>&& other) {
//        this->is_element = other;
//        this->
//    }
    DenseHashMap(bit_vector *isElement,rank_support *rank,T default_element = T()) {
        int total_num_elements = rank->rank(rank->size());
        this->is_element = isElement;
        this->rank = rank;
        elements = decltype(elements)(total_num_elements,default_element);
    }

    int64_t underlying_position(int64_t n) const {
        assert(n >= 0);
        assert(n < is_element->size());
        assert((*is_element)[n]);
        return rank->rank(n);
    }

    T& operator[](int64_t n) {
        return elements[underlying_position(n)];
    }

    const T& at(int64_t n) const {
        return elements.at(underlying_position(n));
    }



    bool count(int64_t n) const {
        assert(n < is_element->size());
        return (*is_element)[n];
    }

    T* ptr_to(int64_t n) {
        return &elements[underlying_position(n)];
    }


    bit_vector * is_element;
    rank_support * rank;
    std::vector<T> elements;
};

template <typename T,typename bit_vector=sdsl::bit_vector, typename rank_support=sdsl::rank_support_v<1>>
class ChunkedDenseHashMap : public DenseHashMap<T,bit_vector,rank_support> {
public:
    ChunkedDenseHashMap() = default;
    
    ChunkedDenseHashMap(bit_vector *isElement,rank_support *rank,T default_element = T(),int chunks=omp_get_num_threads()*100) {
        int total_num_elements = rank->rank(rank->size());
        this->is_element = isElement;
        this->rank = rank;
        divisor = get_divisor(total_num_elements,chunks);
        this->elements = decltype(this->elements)(chunks,default_element);
    }

    static int64_t get_divisor(int64_t num_elements,int64_t chunks) {
        return (num_elements+chunks-1)/chunks;
    }


    int64_t unchunked_position(int64_t n) {
        return this->rank->rank(n);
    }

    int64_t chunk(int64_t n) {
        unchunked_position(n)/divisor;
    }

    int64_t position_in_chunk(int64_t n) {
        return unchunked_position(n) % divisor;
    }

    int64_t underlying_position(int64_t n) const {
        assert(n >= 0);
        assert(n < this->is_element->size());
        assert((*this->is_element)[n]);
        return chunk(n);
    }

    T& operator[](int64_t n) {
        return this->elements[underlying_position(n)];
    }

    const T& at(int64_t n) const {
        return this->elements.at(underlying_position(n));
    }

    bool count(int64_t n) const {
        assert(n < this->is_element->size());
        return (*this->is_element)[n];
    }

    T* ptr_to(int64_t n) {
        return &(this->elements[underlying_position(n)]);
    }
    
    int64_t divisor = 0;
};

#endif //METAGRAPH_DENSE_HASHMAP_HPP
