#include <utility>

//
// Created by Jan Studený on 2019-05-06.
//

#ifndef METAGRAPH_INCOMING_TABLE_HPP
#define METAGRAPH_INCOMING_TABLE_HPP

#include <sdsl/enc_vector.hpp>

using default_bit_vector = bit_vector_small;

using MultiplicityT = sdsl::int_vector<>;
template <typename BitVector_=default_bit_vector,typename GraphT=DBGSuccinct>
class IncomingTable {
public:
    using edge_identifier_t = node_index;
    using BitVector = BitVector_;
    IncomingTable() = default;
    explicit IncomingTable(shared_ptr<const GraphT> graph, BitVector joins,
                           MultiplicityT edge_multiplicity_table) :
                                            graph(graph), joins(joins), edge_multiplicity_table(std::move(edge_multiplicity_table)) {}
    using bit_vector_t = BitVector;
    BitVector joins;
    MultiplicityT edge_multiplicity_table;
    int64_t branch_offset(node_index node,node_index prev_node) const {
        int64_t branch_offset = relative_offset(node,prev_node);
        int64_t result = 0;
        for(int64_t i=0;i<branch_offset;i++) {
            result += branch_size_rank(node,i);
        }
        return result;
    }



    bool is_join(node_index node) const {
        return size(node); // as 1
    }

    int64_t branch_size_rank(node_index node,int64_t offset) const {
        if (offset < 0) { return 0; }
        int64_t joins_position = joins.select1(node);
        int64_t table_offset = joins.rank0(joins_position);
        if (offset == size(node)) {
            return INT_MIN;
        }
        assert(offset < size(node));
        return edge_multiplicity_table[table_offset+offset];
    }

    int64_t branch_size(node_index node,node_index prev_node) const {
        return branch_size_rank(node,relative_offset(node,prev_node));
    }

    bool has_size(node_index node,node_index prev_node) const {
#ifndef FULL_INCOMING_TABLE
        int64_t branch_offset = relative_offset(node,prev_node);
        return branch_offset < size(node);
#else
        return true;
#endif
    }

    int64_t size(node_index node) const {
        return joins.select1(node+1)-joins.select1(node)-1;
    }

    bool has_new_reads(node_index node) const {
        // todo : remove or after the indegree bug is fixed and use assertion
        // assert(!(graph->indegree(node) < 2 and size(node)) or size(node) > graph->indegree(node));
        // indegree smaller than two with nonempty incoming table implies that node has new reads
#ifndef FULL_INCOMING_TABLE
        return size(node) == graph->indegree(node) or (graph->indegree(node) < 2 and size(node));
#else
        return (size(node) > graph->indegree(node)) or (graph->indegree(node) < 2 and size(node));
#endif
    }

    int64_t relative_offset(node_index node,node_index prev_node) const {
        bool increment = has_new_reads(node);
        int64_t result;
        if (prev_node) {
            result = graph->branch_id(node,prev_node);
        }
        else {
            result = -1; // starting branch is before all other branches
        }
        return result + increment;
    }

    string print_content(node_index node) const {
        stringstream out;
        auto table_size = size(node);
//        for(char c : "$ACGTN") {
//            out << c << ":" << branch_size(node, c) << endl;
//        }
        for (int i=0;i<table_size;i++) {
            cout << i << ":" << branch_size_rank(node,i) << endl;
        }
        cerr << out.str();
        return out.str();
    }

    shared_ptr<const GraphT> graph;
};

#endif //METAGRAPH_INCOMING_TABLE_HPP
