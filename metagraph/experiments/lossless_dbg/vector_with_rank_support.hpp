//
// Created by Jan Studený on 2019-06-02.
//

#ifndef METAGRAPH_VECTOR_WITH_RANK_SUPPORT_HPP
#define METAGRAPH_VECTOR_WITH_RANK_SUPPORT_HPP

#include <vector>
#include <iostream>

using namespace std;

template <typename T>
class VectorWithRankSupport : public vector<T> {
public:
    using vector<T>::vector;
    // rank [0..position]
    int64_t rank(T symbol,int64_t position) const {
        int64_t result = 0;
        for (size_t i = 0; i <= position; ++i) {
            result += vector<T>::operator[](i) == symbol;
        }
        return result;
    }

    int64_t select(T symbol,int64_t rank) const {
        int64_t crank = 1;
        int64_t i = 0;
        for(auto& c : *this) {
            if (c == symbol) {
                if (crank == rank) {
                    return i;
                }
                else {
                    crank++;
                }
            }
            i++;
        }
        return INT_MAX;
    }

    char get(int64_t position) const {
        return this->at(position);
    }

    void insert(int64_t position,T symbol) {
        assert(position <= this->size());
        vector<T>::insert(this->begin() + position, symbol);
    }
};

#include <utility>
#include <iostream>
#include <set>
#include <map>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <boost/range/size_type.hpp>
#include <tsl/hopscotch_map.h>
#include "path_database.hpp"
#include "graph_patch.hpp"
#include "path_database_dynamic.hpp"
#include "utils.hpp"
#include "utilities.hpp"
#include "alphabets.hpp"
#include "graph_preprocessor.hpp"
#include "routing_table_transformation.hpp"
#include "dense_hashmap.hpp"
#include "solid_dynamic_routing_table.hpp"

#endif //METAGRAPH_VECTOR_WITH_RANK_SUPPORT_HPP
