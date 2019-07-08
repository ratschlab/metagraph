//
// Created by studenyj on 7/8/19.
//
#include <iostream>
#include <vector>
#include <omp.h>
#include "utilities.hpp"

using node_index = uint64_t;
using namespace std;


class Graph {
public:
    static const ll bits = 1'000'000;
    static ll num_nodes() {
        return bits;
    }
} graph;

sdsl::bit_vector is_join;
sdsl::bit_vector is_split;
sdsl::bit_vector is_bifurcation;
int main(int argc, char** argv) {
    ll bits_to_set = graph.num_nodes()/10;
    auto additional_splits_t = VerboseTimer("computing additional splits and joins");
    is_split = decltype(is_split)(graph.num_nodes() + 1); // bit
    is_join = decltype(is_join)(graph.num_nodes() + 1);
    is_bifurcation = decltype(is_join)(graph.num_nodes() + 1);
    //
    // uint64_t chunk_size = nodes / chunks + 64ull
//    // chunk_size &= ~63ull
//    #pragma omp parallel for
//    for (size_t i = 0; i < bits_to_set; ++i) {
//        // add additional bifurcation
//        ll node =  graph.kmer_to_node(
//
//        );
//        assert(additional_joins_vec[i]); // node has to be in graph
//        additional_splits_vec[i] = graph.kmer_to_node(
//                transformed_sequence.substr(transformed_sequence.length() - graph.get_k())
//        );
//        assert(additional_splits_vec[i]); // node has to be in graph
//    }
//    auto bifurcation_timer = VerboseTimer("construction of bifurcation bit_vectors");
//
//#pragma omp parallel for num_threads(get_num_threads())
//    for (uint64_t node = 0; node <= graph.num_nodes(); node += 64) {
//        for (int64_t i = node; i < node + 64 && i <= graph.num_nodes(); ++i) {
//            if (!i)
//                continue;
//            is_split[i] = node_is_split_raw(i);
//            is_join[i] = node_is_join_raw(i);
//            is_bifurcation[i] = is_split[i] || is_join[i];
//        }
//    }


    cout << get_used_memory();
}