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

//struct exit_barrier_element_t {
//    int64_t relative_position;
//    uint64_t edge;// node + last 8 bits is traversed edge
//    friend ostream& operator<<(ostream& os, const exit_barrier_element_t& dt);
//};
#define edges first
#define relative_positions second
constexpr int max_threads = 32;
using exit_barrier_t = pair<array<uint64_t,max_threads>,array<uint32_t,max_threads>>;

template <typename BitVector=sdsl::bit_vector,typename RankSupport=sdsl::rank_support_v<1>>
class ExitBarrier {
public:
    ExitBarrier(BitVector* is_element,RankSupport* rank_is_element,int chunk_size=DefaultChunks) :
        exit_barriers(is_element,rank_is_element,chunk_size)
#ifndef DISABLE_PARALELIZATION
        , exit_barrier_locks(is_element, rank_is_element, chunk_size)
#endif
        {
        if (get_num_threads() > max_threads) {
            throw "Implementation doesn't work with more than "s + to_string(max_threads) + " threads"s;
        }
        for(int64_t i=0;i<exit_barriers.elements.size();i++) {
            using BarrierT = typename decltype(exit_barriers)::element_type;
            //exit_barriers.elements[i] = BarrierT(get_num_threads(),{0});
            exit_barriers.elements[i] = {{0},{0}};
#ifndef DISABLE_PARALELIZATION
            omp_init_lock(&exit_barrier_locks.elements[i]);
#endif
        }
    }
    ~ExitBarrier(){
        for(int64_t i=0;i<exit_barriers.elements.size();i++) {
#ifndef DISABLE_PARALELIZATION
            omp_destroy_lock(&exit_barrier_locks.elements[i]);
#endif
        }
    }

    ChunkedDenseHashMap<exit_barrier_t,BitVector, RankSupport,false> exit_barriers;
#ifndef DISABLE_PARALELIZATION
    ChunkedDenseHashMap<omp_lock_t,BitVector,RankSupport,false> exit_barrier_locks;
#endif

    int64_t exit(uint64_t node, int tid) {
#ifndef DISABLE_PARALELIZATION
        auto exit_barrier_lock = exit_barrier_locks.ptr_to(node);
        omp_set_lock(exit_barrier_lock);
#endif
        auto& exit_barrier = exit_barriers[node];
        auto& my_edge = exit_barrier.edges[tid];
        assert(my_edge);
        auto result = update_myself_after_wakeup(exit_barrier,tid);
        // remove myself after wakeup
        my_edge = 0;
        #ifndef DISABLE_PARALELIZATION
        omp_unset_lock(exit_barrier_lock);
        #endif
        return result;
    }

    void enter(uint64_t node, char traversed_edge, int64_t relative_position, int tid) {
        #ifndef DISABLE_PARALELIZATION
        auto exit_barrier_lock = exit_barrier_locks.ptr_to(node);
        omp_set_lock(exit_barrier_lock);
#endif
        // add myself before sleep
        auto& exit_barrier = exit_barriers[node];
        auto& stored_edge = exit_barrier.edges[tid];
        auto& stored_relative_position = exit_barrier.relative_positions[tid];
        stored_edge = node | (((uint64_t)traversed_edge)<<(7*8));
        stored_relative_position = relative_position;
        update_others_before_sleep(exit_barrier,tid);
        #ifndef DISABLE_PARALELIZATION
        omp_unset_lock(exit_barrier_lock);
#endif
    }

    void update_others_before_sleep(exit_barrier_t &exit_barrier, int tid) {
        auto my_edge = exit_barrier.edges[tid];
        auto my_relative_position = exit_barrier.relative_positions[tid];

        for(int i=0;i<max_threads;i++) {
            auto& other_edge = exit_barrier.edges[i];
            auto& other_relative_position = exit_barrier.relative_positions[i];
//            if (i == tid) continue;
//            if (other.node == myself.node and
//                other.traversed_edge == myself.traversed_edge and
//                other.relative_position >= myself.relative_position) {
//                other.relative_position++;// I am below other
//            }
            // vectorized version
            other_relative_position += other_edge == my_edge and
                    other_relative_position >= my_relative_position and
                                       i != tid;
        }
    }

    int64_t update_myself_after_wakeup(exit_barrier_t &exit_barrier,int tid) {
        auto my_edge = exit_barrier.edges[tid];
        auto my_relative_position = exit_barrier.relative_positions[tid];
        int offset_change = 0;

        for(int i=0;i<max_threads;i++) {
            auto& other_edge = exit_barrier.edges[i];
            auto& other_relative_position = exit_barrier.relative_positions[i];
//            if (i == tid) continue;
//            if (other.node == myself.node and
//                other.traversed_edge == myself.traversed_edge and
//                other.relative_position < myself.relative_position) {
//                offset_change--;// Other is below me
//            }
            offset_change -= other_edge == my_edge and
                    other_relative_position < my_relative_position and
                i != tid;
        }
        assert(my_relative_position + offset_change >= 0);
        return my_relative_position + offset_change;
    }

    void print_content(uint64_t node, uint64_t tid=0) {
        cout << exit_barriers[node] << endl;
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
        omp_set_lock(waiting_queue_lock);
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

    void print_content(uint64_t node, uint64_t tid=0) {
        cout << waiting_threads[node] << endl;
    }

    ChunkedDenseHashMap<deque<waiting_thread_info_t>,BitVector, RankSupport,false> waiting_threads;
    ChunkedDenseHashMap<omp_lock_t,BitVector,RankSupport,false> waiting_queue_locks;
};

#endif //METAGRAPH_waiting_queue_HPP
