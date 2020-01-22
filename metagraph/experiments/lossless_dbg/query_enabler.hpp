//
// Created by Jan Studený on 2019-05-21.
//

#ifndef __QUERY_ENABLER_HPP__
#define __QUERY_ENABLER_HPP__

#include <iostream>
#include <set>
#include <map>
#include <tsl/hopscotch_set.h>
#include <optional>

#include "path_database.hpp"
#include "dynamic_routing_table.hpp"
#include "dynamic_incoming_table.hpp"
#include "unix_tools.hpp"
#include "threading.hpp"

template <typename Database>
class QueryEnabler : public Database {
  public:
    using Database::Database;

    using range_t = pair<int64_t, int64_t>;
    using history_t = vector<range_t>;
    using score_t = int64_t;
    using next_nodes_with_extended_info_t = map<node_index, pair<score_t, history_t>>;

    history_t get_initial_history(node_index node) const {
        history_t history = { { 0, get_coverage(node) } };
        return history;
    }

    static bool range_is_empty(range_t range) { return range.first >= range.second; }

    history_t gather_history(const string &str_history) {
        // todo: add function to gather only consistent history (with strong aka first class, 0-th order support),
        //       now we are gathering reads that are inconsistent with longest path #future-work
        node_index node
                = this->graph.kmer_to_node(str_history.substr(0, this->graph.get_k()));
        assert(node);
        auto history = get_initial_history(node);
        for (auto &c : str_history.substr(this->graph.get_k())) {
            auto nodes_with_support = get_next_nodes_with_support(node, history);
            node = this->graph.traverse(node, c);
            history = nodes_with_support[node].second;
            if (history.empty()) {
                break;
            }
        }
        return history;
    }

    node_index get_next_consistent_node(node_index node) {
        return get_next_consistent_node(node, this->graph.kmer);
    }

    node_index get_next_consistent_node(node_index node, const string &str_history) {
        // consistent node should have score 0 and be the only node to go to
        assert(node
               == this->graph.kmer_to_node(
                       str_history.substr(str_history.size() - this->graph.get_k())));
        auto history = gather_history(str_history);
        auto support = get_next_nodes_with_support(node, history);
        node_index consistent_node = 0;
        for (auto &[next_node, info] : support) {
            auto &[score, history] = info;
            if (score == 0) {
                if (!consistent_node) {
                    consistent_node = next_node;
                } else {
                    consistent_node = 0;
                    break;
                }
            } else {
                // histories are sorted
                break;
            }
        }
        return consistent_node;
    }


    next_nodes_with_extended_info_t get_next_nodes_with_support(node_index node,
                                                                history_t &history) {
        // todo: merge ranges for successive joining reads #future-work
        next_nodes_with_extended_info_t result;
        if (this->node_is_split(node)) {
            int64_t range_score = 0;
            for (auto &range : history) {
                for (auto &c : "ACGT"s) {
                    range_t new_range;
                    new_range.first = this->routing_table.rank(node, range.first, c);
                    new_range.second = this->routing_table.rank(node, range.second, c);
                    if (!range_is_empty(new_range)) {
                        auto new_node = this->graph.traverse(node, c);
                        if (!result.count(new_node)) {
                            result[new_node].first = range_score;
                        }
                        result[new_node].second.push_back(new_range);
                    }
                }
                range_score++;
            }
        } else {
            auto base = get_outgoing_base(this->graph, node);
            auto new_node = this->graph.traverse(node, base);
            result[new_node].first = 0;
            result[new_node].second = history;
        }
        for (auto &[new_node, additional_info] : result) {
            auto &[score, result_history] = additional_info;
            if (this->node_is_join(new_node)) {
                auto offset = this->incoming_table.branch_offset(new_node, node);
                for (auto &range : result_history) {
                    range.first += offset;
                    range.second += offset;
                }
                if (this->number_of_reads_starting_at_node(new_node)) {
                    result_history.push_back(
                            { 0, this->number_of_reads_starting_at_node(new_node) });
                }
            }
        }
        return result;
    }

    vector<path_id> get_paths_going_through(node_index node) const {
        // catch: relative indices in node can be from the same sequence if the read is going there multiple times
        // todo: be smarter and stop when arrived to the starting node again #future-work
        auto coverage = get_coverage(node);
        set<path_id> out;
        for (int64_t position = 0; position < coverage; position++) {
            out.insert(this->get_global_path_id(node, position));
        }
        return vector<path_id>(all(out));
    }

    int64_t get_coverage(node_index node) const {
        // one can also traverse backwards
        if (this->node_is_split(node)) {
            return this->routing_table.size(node);
        }

        char base;
        this->graph.call_outgoing_kmers(node, [&base](node_index, char edge_label) {
            base = edge_label;
        });

        auto new_node = this->graph.traverse(node, base);

        if (this->node_is_join(new_node)) {
            if (this->incoming_table.has_size(new_node, node)) {
                return this->incoming_table.branch_size(new_node, node);
            } else {
                return get_coverage(new_node)
                        - this->incoming_table.branch_offset(new_node, node);
            }
        }

        return get_coverage(new_node);
    }
};
#endif // __QUERY_ENABLER_HPP__
