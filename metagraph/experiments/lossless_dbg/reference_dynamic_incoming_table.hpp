//
// Created by Jan Studený on 2019-06-04.
//

#ifndef __REFERENCE_DYNAMIC_INCOMING_TABLE_HPP__
#define __REFERENCE_DYNAMIC_INCOMING_TABLE_HPP__

#include "utilities.hpp"

template <typename GraphT = DBGSuccinct, typename _edge_identifier_t = char>
class ReferenceDynamicIncomingTable {
  public:
    using edge_identifier_t = _edge_identifier_t;
    explicit ReferenceDynamicIncomingTable(shared_ptr<const GraphT> /* graph */) {}

    template <typename BitVector, typename RankSupport>
    ReferenceDynamicIncomingTable(shared_ptr<const DBGSuccinct> /* graph */,
                                  BitVector *is_element,
                                  RankSupport *rank_element,
                                  ll /* chunks = DefaultChunks */)
        : incoming_table(is_element, rank_element) {}


    int branch_offset(node_index node, edge_identifier_t incoming) const {
        int result = 0;
        for (char base : "$ACGTN") {
            if (base < incoming) {
                result += branch_size(node, base);
            }
        }
        return result;
    }

    bool is_join(node_index node) const { return incoming_table.count(node); }

    int branch_size(node_index node, edge_identifier_t incoming) const {
        assert(node);
        int result = 0;
        int encoded = encode(incoming);
        if (incoming_table.count(node)) {
            result = incoming_table.at(node)[encoded];
        }
        return result;
    }

    int size(node_index node) const {
        int result = 0;
        for (char base : "$ACGTN") {
            result += branch_size(node, base);
        }
        return result;
    }

    int branch_offset_and_increment(node_index node, edge_identifier_t incoming) {
        assert(node);
        assert(incoming == '$' or incoming == 'A' or incoming == 'C' or incoming == 'G'
               or incoming == 'T' or incoming == 'N');
        int result = 0;
        int encoded = encode(incoming);
        auto &array = incoming_table[node];
        for (auto i = 0; i < 6; i++) {
            if (i >= encoded)
                break;
            result += array[i];
        }
        array[encoded]++;
        return result;
    }

    string print_content(node_index node) const {
        stringstream out;
        for (char c : "$ACGTN") {
            out << c << ":" << branch_size(node, c) << endl;
        }
        mg::common::logger->debug(out.str());
        return out.str();
    }

    bool has_new_reads(node_index node) const { return branch_size(node, '$'); }

    DenseHashMap<array<int, 6>> incoming_table;
};
#endif // __REFERENCE_DYNAMIC_INCOMING_TABLE_HPP__
