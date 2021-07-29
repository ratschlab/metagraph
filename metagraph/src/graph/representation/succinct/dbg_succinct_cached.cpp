#include "dbg_succinct_cached.hpp"

#include "common/seq_tools/reverse_complement.hpp"

namespace mtg {
namespace graph {


auto DBGSuccinctCached::get_rev_comp_boss_next_node(node_index node) const -> edge_index {
    // 78% effective
    edge_index ret_val = 0;
    size_t chars_unmatched = 0;
    if (auto fetch = rev_comp_next_cache_.TryGet(node)) {
        std::tie(ret_val, chars_unmatched) = *fetch;
        if (!ret_val && !chars_unmatched)
            return 0;
    } else {
        //   AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //->  GAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //->  CTCAAGCAGAAGACGGCATACGAGATCCTC
        std::string rev_seq = get_node_sequence(node).substr(1, get_k() - 1);
        if (rev_seq[0] == boss::BOSS::kSentinel) {
            rev_comp_next_cache_.Put(node, std::make_pair(0, 0));
            return 0;
        }

        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        auto encoded = boss_->encode(rev_seq);
        auto [edge, edge_2, end] = boss_->index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end() || edge == edge_2);

        if (end == encoded.end()) {
            ret_val = edge;
        } else {
            chars_unmatched = encoded.end() - end;
        }

        auto stored_val = std::make_pair(ret_val, chars_unmatched);
        edge_index start = boss_->pred_W(graph_->kmer_to_boss_index(node),
                                         complement(encoded.front()));
        boss_->call_incoming_to_target(start, complement(encoded.front()),
                                       [&](edge_index adj) {
            if (node_index a = graph_->boss_to_kmer_index(adj))
                rev_comp_next_cache_.Put(a, stored_val);
        });
    }

    if (!ret_val && chars_unmatched > 1) {
        // we only matched a suffix
        // $$CTCAAGCAGAAGACGGCATACGAGATCC
        // so if we try to match the next k-mer
        // GAGGATCTCGTATGCCGTCTTCTGCTTGAGX
        // corresponding to
        //  xCTCAAGCAGAAGACGGCATACGAGATCCT
        // we may find the subrange
        // $xCTCAAGCAGAAGACGGCATACGAGATCC
        // or something shorter
        size_t next_chars_unmatched = chars_unmatched - 1;
        auto stored_val = std::make_pair(0, next_chars_unmatched);
        adjacent_outgoing_nodes(node, [&](node_index next) {
            rev_comp_next_cache_.Put(next, stored_val);
#ifndef NDEBUG
            std::string test_seq = graph_->get_node_sequence(next).substr(1, get_k() - 1);
            ::reverse_complement(test_seq.begin(), test_seq.end());
            auto encoded = boss_->encode(test_seq);
            auto [edge, edge_2, end] = boss_->index_range(encoded.begin(), encoded.end());
            assert(end != encoded.end());
            if (next_chars_unmatched > static_cast<size_t>(encoded.end() - end)) {
                common::logger->error("Unmatched chars incorrectly stored: {} vs. {}",
                                      next_chars_unmatched, encoded.end() - end);
                assert(false);
            }
#endif
        });
    } else if (ret_val) {
        // traverse node forward by X, then ret_val back by comp(X)
        std::vector<std::pair<node_index, edge_index>> parents(alphabet().size());
        size_t parents_count = 0;
        edge_index start = 0;
        call_outgoing_kmers(node, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel) {
                parents[boss_->encode(c)].first = next;
                edge_index next_edge = graph_->kmer_to_boss_index(next);
                if (boss_->get_last(next_edge))
                    start = next_edge;

                ++parents_count;
            }
        });

        if (parents_count) {
            edge_index ret_prev = boss_->bwd(ret_val);
            TAlphabet w = boss_->get_node_last_value(ret_val);
            boss_->call_incoming_to_target(ret_prev, w, [&](edge_index prev_edge) {
                if (boss_->get_last(prev_edge)) {
                    TAlphabet s = complement(get_first_value(prev_edge));
                    if (parents[s].first)
                        rev_comp_next_cache_.Put(parents[s].first, std::make_pair(prev_edge, 0));
                }
            });
        }
    }

#ifndef NDEBUG
    std::string test_seq = graph_->get_node_sequence(node).substr(1, get_k() - 1);
    ::reverse_complement(test_seq.begin(), test_seq.end());
    if (ret_val) {
        assert(test_seq == boss_->get_node_str(ret_val));
    } else {
        auto encoded = boss_->encode(test_seq);
        auto [edge, edge_2, end] = boss_->index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end());
        if (chars_unmatched > static_cast<size_t>(encoded.end() - end)) {
            common::logger->error("Unmatched chars incorrectly stored: {} vs. {}",
                                  chars_unmatched, encoded.end() - end);
            assert(false);
        }
    }
#endif

    return ret_val;
}

auto DBGSuccinctCached::get_rev_comp_boss_prev_node(node_index node) const -> edge_index {
    // 8% effective
    edge_index ret_val = 0;
    size_t chars_unmatched = 0;
    edge_index edge_1 = 0;
    edge_index edge_2 = 0;
    std::vector<TAlphabet> encoded;
    std::vector<TAlphabet>::iterator end = encoded.begin();
    if (auto fetch = rev_comp_prev_cache_.TryGet(node)) {
        std::tie(ret_val, chars_unmatched) = *fetch;
        if (!ret_val && !chars_unmatched)
            return 0;
    } else {
        //   AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //-> AGAGGATCTCGTATGCCGTCTTCTGCTTGA
        //-> TCAAGCAGAAGACGGCATACGAGATCCTCT
        std::string rev_seq = get_node_sequence(node);
        if (rev_seq[0] == boss::BOSS::kSentinel) {
            rev_comp_prev_cache_.Put(node, std::make_pair(0, 0));
            return 0;
        }

        rev_seq.pop_back();
        ::reverse_complement(rev_seq.begin(), rev_seq.end());
        auto encoded = boss_->encode(rev_seq);
        std::tie(edge_1, edge_2, end) = boss_->index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end() || edge_1 == edge_2);
        if (end == encoded.end()) {
            ret_val = edge_1;
        } else {
            chars_unmatched = encoded.end() - end;
        }

        auto stored_val = std::make_pair(ret_val, chars_unmatched);
        boss_->call_outgoing(boss_->succ_last(graph_->kmer_to_boss_index(node)),
                             [&](edge_index adj) {
            if (node_index a = graph_->boss_to_kmer_index(adj))
                rev_comp_prev_cache_.Put(a, stored_val);
        });
    }

    if (!ret_val && chars_unmatched > 1 && encoded.size()) {
        // TODO: given an LCS array, we can do the same check as get_rev_comp_boss_next_node
        //       to find neighbouring nodes with rev comp matches
        //       e.g., we looked for
        //       AGAGGATCTCGTATGCCGTCTTCTGCTTGAG
        //       TCAAGCAGAAGACGGCATACGAGATCCTCT
        //       and matched a suffix
        //       XXXXXXXXTCAAGCAGAAGACGGCATACGA
        //       so if we try to match the previous k-mer
        //       XAGAGGATCTCGTATGCCGTCTTCTGCTTGA
        //       corresponding to
        //        CAAGCAGAAGACGGCATACGAGATCCTCTx
        //       we can use the LCS array to delete the first character from the match
        //       and work from there
        // expand the range up
        bool expand_up_exists = false;
        bool expand_down_exists = false;
        auto cur_seq = decode(encoded);
        if (edge_1 > 1) {
            edge_index edge_prev = boss_->pred_last(edge_1 - 1);
            if (node_index prev = graph_->boss_to_kmer_index(edge_prev)) {
                auto prev_seq = get_node_sequence(prev).substr(0, get_k() - 1);
                if (std::equal(prev_seq.begin() + chars_unmatched + 1,
                               prev_seq.end(), cur_seq.begin() + 1)) {
                    expand_up_exists = true;
                }
            }
        }

        // expand the range down
        if (edge_2 < boss_->num_edges()) {
            edge_index edge_next = boss_->succ_last(edge_2 + 1);
            if (node_index next = graph_->boss_to_kmer_index(edge_next)) {
                auto next_seq = get_node_sequence(next).substr(0, get_k() - 1);
                if (std::equal(next_seq.begin() + chars_unmatched + 1,
                               next_seq.end(), cur_seq.begin() + 1)) {
                    expand_down_exists = true;
                }
            }
        }

        if (!expand_up_exists && !expand_down_exists) {
            size_t next_chars_unmatched = chars_unmatched - 1;
            auto stored_val = std::make_pair(0, next_chars_unmatched);
            adjacent_incoming_nodes(node, [&](node_index prev) {
                rev_comp_prev_cache_.Put(prev, stored_val);
#ifndef NDEBUG
                std::string test_seq = graph_->get_node_sequence(node);
                test_seq.pop_back();
                ::reverse_complement(test_seq.begin(), test_seq.end());
                auto encoded = boss_->encode(test_seq);
                auto [edge, edge_2, end] = boss_->index_range(encoded.begin(), encoded.end());
                assert(end != encoded.end());
                if (next_chars_unmatched > static_cast<size_t>(encoded.end() - end)) {
                    common::logger->error("Unmatched chars incorrectly stored: {} vs. {}",
                                          next_chars_unmatched, encoded.end() - end);
                    assert(false);
                }
#endif
            });
        }
    } else if (ret_val) {
        // traverse ret_val forward by X, then node back by comp(X)
        std::vector<std::pair<node_index, edge_index>> parents(alphabet().size());
        size_t parents_count = 0;
        if (graph_->boss_to_kmer_index(ret_val)) {
            boss_->call_outgoing(ret_val, [&](edge_index next_edge) {
                if (TAlphabet w = boss_->get_W(next_edge) % boss_->alph_size) {
                    parents[w].second = boss_->fwd(boss_->pred_W(ret_val, w), w);
                    ++parents_count;
                }
            });
        }

        if (parents_count) {
            call_incoming_kmers(node, [&](node_index prev, char c) {
                TAlphabet s = complement(boss_->encode(c));
                if (parents[s].second)
                    rev_comp_prev_cache_.Put(prev, std::make_pair(parents[s].second, 0));
            });
        }
    }

#ifndef NDEBUG
    std::string test_seq = graph_->get_node_sequence(node);
    test_seq.pop_back();
    ::reverse_complement(test_seq.begin(), test_seq.end());
    if (ret_val) {
        assert(test_seq == boss_->get_node_str(ret_val));
    } else {
        auto encoded = boss_->encode(test_seq);
        auto [edge, edge_2, end] = boss_->index_range(encoded.begin(), encoded.end());
        assert(end != encoded.end());
        assert(chars_unmatched <= static_cast<size_t>(encoded.end() - end));
    }
#endif

    return ret_val;
}


} // namespace graph
} // namespace mtg
