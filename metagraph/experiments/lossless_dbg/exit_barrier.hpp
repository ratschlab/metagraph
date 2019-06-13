//
// Created by Jan Studený on 2019-06-13.
//
#include <vector>
#include <iostream>
#include <omp.h>

#include "dense_hashmap.hpp"
using namespace std;

#ifndef METAGRAPH_waiting_queue_HPP
#define METAGRAPH_waiting_queue_HPP

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

template <typename BitVector=sdsl::bit_vector,typename RankSupport=sdsl::rank_support_v<1>>
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

    void update_others_before_sleep(exit_barrier_t &exit_barrier, int tid) {
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

    int64_t update_myself_after_wakeup(exit_barrier_t &exit_barrier,int tid) {
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

struct waiting_thread_info_t {
    int16_t thread_id;
    int64_t relative_position;
    uint64_t node;
    char traversed_edge;
    friend ostream& operator<<(ostream& os, const waiting_thread_info_t& dt);
};
ostream& operator<<(ostream& os, const waiting_thread_info_t& wt)
{
    os << "[tid=" << wt.thread_id << ",rp=" << wt.relative_position << ",node=" << wt.node << ",tedg=" << wt.traversed_edge << "]";
    return os;
}

template <typename BitVector=sdsl::bit_vector,typename RankSupport=sdsl::rank_support_v<1>>
class ReferenceExitBarrier {
public:
ReferenceExitBarrier(BitVector* is_element,RankSupport* rank_is_element,int chunk_size=DefaultChunks) :
        waiting_threads(is_element,rank_is_element,chunk_size), waiting_queue_locks(is_element, rank_is_element, chunk_size) {
    for(int64_t i=0;i<waiting_threads.elements.size();i++) {
        using BarrierT = typename decltype(waiting_threads)::element_type;
        waiting_threads.elements[i] = BarrierT(omp_get_num_threads(),{0});
        omp_init_lock(&waiting_queue_locks.elements[i]);
    }
}

~ReferenceExitBarrier(){
    for(int64_t i=0;i<waiting_threads.elements.size();i++) {
        omp_destroy_lock(&waiting_queue_locks.elements[i]);
    }
}

    int64_t exit(int64_t node, int tid) {
        tid++;
        auto waiting_queue_lock = waiting_queue_locks.ptr_to(node);
        omp_set_lock(waiting_queue_lock);
        auto& waiting_queue = waiting_threads[node];
        waiting_thread_info_t myself;
        for(auto&other : waiting_queue) {
            if (other.thread_id == tid) {
                myself = other;
            }
        }
        auto result = update_myself_after_wakeup(tid,myself.node,myself.traversed_edge,waiting_queue,myself.relative_position);
        omp_unset_lock(waiting_queue_lock);
        return result;
    }

    void enter(int64_t node, char traversed_edge, int64_t relative_position, int tid) {
        tid++;
        auto waiting_queue_lock = waiting_queue_locks.ptr_to(node);
        waiting_threads[node].push_back({(int16_t)tid, relative_position, node, traversed_edge});
        omp_unset_lock(waiting_queue_lock);
    }


    int64_t update_myself_after_wakeup(int tid,
                                       int64_t previous_node, char traversed_edge,
                                       deque<waiting_thread_info_t> &waiting_queue,
                                       int64_t relative_position) {
        assert(previous_node);
        bool first_it = 1;
        bool me_first = 0;
        bool was_me = 0;
        int64_t past_offset = 0;
        // todo_wrap in define or better in macro
        int64_t debug_my_id = 0;
        int64_t debug_idx = 0;
        for(auto&other : waiting_queue) {
            if (was_me) {
                // future
                if (other.node == previous_node and
                    other.traversed_edge == traversed_edge and
                    other.thread_id < 0 and // job has already finished
                    other.relative_position <= relative_position) {
#ifdef DEBUG_ORDER_CORRECTION
                    cerr << waiting_queue << endl;
                        PRINT_VAR(tid, node, other.relative_position);
#endif
                    relative_position++;
                }
            }
            else if (tid == other.thread_id) {
                debug_my_id = debug_idx;
                assert(!was_me); // can appear only once
                // present
                was_me = 1;
                if (first_it) {
                    me_first = 1;
                }
                other.thread_id = -tid;
            }  else {
                // past
                if (other.node == previous_node and
                    other.traversed_edge == traversed_edge and
                    other.thread_id > 0 and// job is waiting
                    other.relative_position < relative_position) {
                    past_offset--;
                }
            }
            first_it = 0;
            debug_idx++;
        }
        assert(was_me);
        relative_position += past_offset;
        assert(relative_position >= 0 || [&](){
            cerr << waiting_queue << endl;
            PRINT_VAR(tid,traversed_edge,debug_my_id,past_offset,relative_position);
            return false;
        }());

        if (me_first) {
            while (!waiting_queue.empty() and waiting_queue.front().thread_id < 0) {
                waiting_queue.pop_front();
            }
        }
        return relative_position;
    }

    ChunkedDenseHashMap<deque<waiting_thread_info_t>,BitVector, RankSupport,false> waiting_threads;
    ChunkedDenseHashMap<omp_lock_t,BitVector,RankSupport,false> waiting_queue_locks;
};

#endif //METAGRAPH_waiting_queue_HPP
