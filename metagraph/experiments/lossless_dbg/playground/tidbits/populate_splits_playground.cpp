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
    static const ll bits = 1'000;
     ll num_nodes() {
        return bits;
    }
    node_index kmer_to_node(string dummy="A") {
        return rand()%bits + 1;
    }
    int outdegree(node_index node) {
         return rand()%5;
    }
    int indegree(node_index node) {
         return rand()%5;
    }
} graph;

int chunks = 1000;

sdsl::bit_vector is_join;
sdsl::bit_vector is_split;
sdsl::bit_vector is_bifurcation;
int main(int argc, char** argv) {
    omp_set_num_threads(10);

    ll bits_to_set = graph.num_nodes()/10;
    vector<node_index> debug_join_ids(bits_to_set);

    auto additional_splits_t = VerboseTimer("computing additional splits and joins");
    is_split = decltype(is_split)(graph.num_nodes() + 1); // bit
    is_join = decltype(is_join)(graph.num_nodes() + 1);
    is_bifurcation = decltype(is_join)(graph.num_nodes() + 1);
    uint64_t chunk_size = (graph.num_nodes() + 1)/ chunks + 64ull;
    chunk_size &= ~63ull;
    vector<omp_lock_t> node_locks(chunks);
    for(int i = 0; i < node_locks.size(); i++) {
        omp_init_lock(&node_locks[i]);
    }
    #pragma omp parallel for
    for (size_t i = 0; i < bits_to_set; ++i) {
        ll start_node =  graph.kmer_to_node(

        );
        assert(start_node);
        auto lock_ptr = &node_locks[start_node / chunk_size];
        omp_set_lock(lock_ptr);
        is_join[start_node] = true;
        omp_unset_lock(lock_ptr);
        debug_join_ids[i] = start_node;

        //ll end_node = graph.kmer_to_node(
        //
        //);
        //assert(end_node);
    }
    auto bifurcation_timer = VerboseTimer("construction of bifurcation bit_vectors");

#pragma omp parallel for reduction(append : debug_join_ids)
    for (uint64_t id = 0; id <= graph.num_nodes(); id += 64) {
        for (int64_t node = id; node < id + 64 && node <= graph.num_nodes(); ++node) {
            if (!node)
                continue;
//            is_split[i] = is_split[i] or graph.outdegree(i) > 1;
            auto indegree = graph.indegree(node);
            is_join[node] = is_join[node] or indegree > 1;
            if (indegree > 1) {
                debug_join_ids.push_back(node);
            }
//            is_bifurcation[i] = is_split[i] || is_join[i];
        }
    }
    sort(all(debug_join_ids));
    int debug_i = 0;
    for(ll node=0; node <= graph.num_nodes(); node++) {
        while(debug_join_ids[debug_i] < node && debug_i < debug_join_ids.size()) debug_i++;
        if (debug_join_ids[debug_i] == node) {
            assert(is_join[node]);
        }
        else {
            assert(!is_join[node]);
        }

    }


//cout << get_used_memory();
}