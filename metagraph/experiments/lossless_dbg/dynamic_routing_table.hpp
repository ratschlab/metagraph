#include <utility>

//
// Created by Jan Studený on 2019-04-26.
//

#ifndef METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP
#define METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP

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

using node_index = SequenceGraph::node_index;

template <typename... Args>
void doPrint(std::ostream& out, Args&&... args)
{
    ((out << ',' << std::forward<Args>(args)), ...);
}
#undef protected
template <typename Test,typename Reference>
class IdentityComparator : public Reference {
public:
    using Reference::Reference;
    template<typename ...Args>
    int rank(Args... args) const {
        auto target = Reference::rank(args...);
        auto value = t.rank(args...);
        if (target != value) {
            doPrint(cout,args...);
        }
        assert(target==value);
        return target;
    }

    template<typename ...Args>
    int size(Args... args) const {
        auto target = Reference::size(args...);
        auto value = t.size(args...);
        if (target != value) {
            doPrint(cout,args...);
        }
        assert(target==value);
        return target;
    }

    template<typename ...Args>
    int select(Args... args) const {
        auto target = Reference::select(args...);
        auto value = t.select(args...);
        if (target != value) {
            doPrint(cout,args...);
        }
        assert(target==value);
        return target;
    }

    template<typename ...Args>
    int operator[](Args... args) const {
        auto target = Reference::operator[](args...);
        auto value = t.operator[](args...);
        if (target != value) {
            doPrint(cout,args...);
        }
        assert(target==value);
        return target;
    }

    template<typename ...Args>
    void insert(Args... args) {
        Reference::insert(args...);
        t.insert(args...);
    }

    Test t;
};

template <typename T>
class VectorWithRankSupport : public vector<T> {
public:
    using vector<T>::vector;
    // rank [0..position]
    int rank(T symbol,int position) const {
        int result = 0;
        for (size_t i = 0; i <= position; ++i) {
            result += vector<T>::operator[](i) == symbol;
        }
        return result;
    }

    int select(T symbol,int rank) const {
        int crank = 1;
        int i = 0;
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

    char get(int position) const {
        return this->at(position);
    }

    void insert(int position,T symbol) {
        assert(position <= this->size());
        vector<T>::insert(this->begin() + position, symbol);
    }
};

//template<typename EntryT=wavelet_tree_dyn>
template <typename EntryT=VectorWithRankSupport<int64_t>>
class DynamicRoutingTableCore {
public:
    DynamicRoutingTableCore(const DBGSuccinct& graph) : graph(graph) {}

//    int select(node_index node, int occurrence, char symbol) const {
//    }

    // rank [0..position)
    int rank(node_index node, int position, char symbol) const {
        int encoded = graph.encode(symbol);
        auto &node_entry = routing_table.at(node);
        int result = 0;
        assert(position <= node_entry.size());
        if (position) {
            result = node_entry.rank(encoded, position - 1);
        }
        return result;
    }
public:
    int select(node_index node,int rank,char symbol) const {
        int encoded = graph.encode(symbol);
        auto &node_entry = routing_table.at(node);
        return node_entry.select(encoded,rank);
    }
    char get(node_index node, int position) const {
        auto &node_entry = routing_table.at(node);
        return graph.decode(node_entry[position]);
    }

    string print_content(node_index node) const {
        stringstream out;
        auto table_size = size(node);
        for (int i=0;i<table_size;i++) {
            out << get(node, i);
        }
        out << endl;
        cerr << out.str();
        return out.str();
    }

    char traversed_base(node_index node, int position) const {
        assert(position >= 0 && "traversing o");
        return get(node,position);
    }

    int new_relative_position(node_index node, int position) const {
        auto base = get(node,position);
        auto base_rank = rank(node,position,base);
        return base_rank;
    }

    int size(node_index node) const {
        if (routing_table.count(node))
            return routing_table.at(node).size();
        else
            return 0;
    }

    void insert(node_index node, int position, char symbol) {
        int encoded = graph.encode(symbol);
        assert(position <= size(node));
        routing_table[node].insert(position,encoded);
    }

    DenseHashMap<EntryT> routing_table;
    const DBGSuccinct& graph;
};



//template <typename EntryT=wavelet_tree_dyn>
template <typename EntryT=VectorWithRankSupport<int64_t>>
class DynamicRoutingTable : public TransformationsEnabler<DynamicRoutingTableCore<EntryT>> {
    using TransformationsEnabler<DynamicRoutingTableCore<EntryT>>::TransformationsEnabler;
};

#endif //METAGRAPH_DYNAMIC_ROUTING_TABLE_HPP
