//
// Created by Jan Studený on 2019-06-13.
//
#include <vector>
#include <iostream>
#include <omp.h>

#include "dense_hashmap.hpp"
using namespace std;

#ifndef METAGRAPH_EXIT_BARRIER_HPP
#define METAGRAPH_EXIT_BARRIER_HPP

struct exit_barrier_element_t {
    int64_t relative_position = 0;
    uint64_t node = 0;
    char traversed_edge = 0;
    friend ostream& operator<<(ostream& os, const exit_barrier_element_t& dt);
};

ostream& operator<<(ostream& os, const exit_barrier_element_t& wt)
{
    os << ",rp=" << wt.relative_position << ",node=" << wt.node << ",tedg=" << wt.traversed_edge << "]";
    return os;
}
using exit_barrier_t = vector<exit_barrier_element_t>;

template <typename BitVector,typename RankSupport>
class ExitBarrier {
public:
    ExitBarrier(BitVector* is_element,RankSupport* rank_is_element,int chunk_size=DefaultChunks) :
        exit_barriers(is_element,rank_is_element,chunk_size), exit_barrier_locks(is_element, rank_is_element, chunk_size) {
        for(int64_t i=0;i<exit_barriers.elements.size();i++) {
            using BarrierT = typename decltype(exit_barriers)::element_type;
            exit_barriers.elements[i] = BarrierT(omp_get_num_threads(),{0});
            omp_init_lock(&exit_barrier_locks.elements[i]);
        }
    }
    ~ExitBarrier(){
        for(int64_t i=0;i<exit_barriers.elements.size();i++) {
            omp_destroy_lock(&exit_barrier_locks.elements[i]);
        }
    }

    ChunkedDenseHashMap<exit_barrier_t,BitVector, RankSupport,false> exit_barriers;
    ChunkedDenseHashMap<omp_lock_t,BitVector,RankSupport,false> exit_barrier_locks;

    int64_t exit(int64_t node, int tid) {
        auto exit_barrier_lock = exit_barrier_locks.ptr_to(node);
        omp_set_lock(exit_barrier_lock);
        auto& exit_barrier = exit_barriers[node];
        auto& myself = exit_barrier[tid];
        auto result = update_myself_after_wakeup(exit_barrier,tid);
        // remove myself after wakeup
        myself.node = 0;
        omp_unset_lock(exit_barrier_lock);
        return result;
    }

    void enter(int64_t node, char traversed_edge, int64_t relative_position, int tid) {
        auto exit_barrier_lock = exit_barrier_locks.ptr_to(node);
        omp_set_lock(exit_barrier_lock);
        // add myself before sleep
        auto& exit_barrier = exit_barriers[node];
        auto& myself = exit_barrier[tid];
        myself.node = node;
        myself.traversed_edge = traversed_edge;
        myself.relative_position = relative_position;
        update_others_before_sleep(exit_barrier,tid);
        omp_unset_lock(exit_barrier_lock);
    }

    void update_others_before_sleep(exit_barrier_t &exit_barrier, int tid) const {
        auto& myself = exit_barrier[tid];
        for(int i=0;i<exit_barrier.size();i++) {
            auto& other = exit_barrier[i];
            if (i == tid) continue;
            if (other.node == myself.node and
                other.traversed_edge == myself.traversed_edge and
                other.relative_position >= myself.relative_position) {
                other.relative_position++;// I am below other
            }
        }
    }

    int64_t update_myself_after_wakeup(exit_barrier_t &exit_barrier,int tid) const {
        auto& myself = exit_barrier[tid];
        int offset_change = 0;
        for(int i=0;i<exit_barrier.size();i++) {
            auto& other = exit_barrier[i];
            if (i == tid) continue;
            if (other.node == myself.node and
                other.traversed_edge == myself.traversed_edge and
                other.relative_position < myself.relative_position) {
                offset_change--;// Other is below me
            }
        }
        assert(myself.relative_position + offset_change >= 0);
        return myself.relative_position + offset_change;
    }

};

template <typename BitVector,typename RankSupport>
class ReferenceExitBarrier {
ReferenceExitBarrier(BitVector* is_element,RankSupport* rank_is_element,int chunk_size=DefaultChunks) :
        exit_barriers(is_element,rank_is_element,chunk_size), exit_barrier_locks(is_element, rank_is_element, chunk_size) {
    for(int64_t i=0;i<exit_barriers.elements.size();i++) {
        using BarrierT = typename decltype(exit_barriers)::element_type;
        exit_barriers.elements[i] = BarrierT(omp_get_num_threads(),{0});
        omp_init_lock(&exit_barrier_locks.elements[i]);
    }
}

~ReferenceExitBarrier(){
    for(int64_t i=0;i<exit_barriers.elements.size();i++) {
        omp_destroy_lock(&exit_barrier_locks.elements[i]);
    }
}

    ChunkedDenseHashMap<exit_barrier_t,BitVector, RankSupport,false> exit_barriers;
    ChunkedDenseHashMap<omp_lock_t,BitVector,RankSupport,false> exit_barrier_locks;
};

#endif //METAGRAPH_EXIT_BARRIER_HPP
