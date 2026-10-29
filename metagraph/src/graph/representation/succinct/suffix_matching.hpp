#pragma once

#include <cassert>
#include <tuple>
#include <type_traits>
#include <vector>

#include "dbg_succinct.hpp"

namespace mtg::graph {

struct NonSentinelSymbols {
    template <class Range, class TrySymbol>
    void operator()(const boss::BOSS &boss, const Range &, const TrySymbol &try_symbol) const {
        for (boss::BOSS::TAlphabet s = 1; s < boss.alph_size; ++s)
            try_symbol(s);
    }
};

// Enumerate valid k-mers having the range's suffix as a prefix. The range is
// (first edge, last edge, suffix length), with length in [1, k - 1]. Returning
// false from callback stops the entire enumeration.
// Validity follows the graph's mask; an unmasked graph can return sink edges.
//
// symbols(boss, range, try_symbol) selects which symbols to append at each
// depth. Its range has its length already incremented for the children. The
// default tries every non-sentinel symbol, including N in a DNA5 graph.
template <class Callback, class SymbolSet = NonSentinelSymbols>
void suffix_to_prefix(const DBGSuccinct &graph,
                      const std::tuple<boss::BOSS::edge_index,
                                       boss::BOSS::edge_index, size_t> &index_range,
                      const Callback &callback,
                      const SymbolSet &symbols = {}) {
    using Range = std::decay_t<decltype(index_range)>;
    const auto &boss = graph.get_boss();
    assert(std::get<2>(index_range));
    assert(std::get<2>(index_range) < graph.get_k());

    auto call_nodes_in_range = [&](const Range &range) {
        const auto &[first, last, length] = range;
        assert(length == boss.get_k());
        for (auto edge = first; edge <= last; ++edge) {
            auto node = graph.validate_edge(edge);
            if (node && !callback(node))
                return false;
        }
        return true;
    };

    if (std::get<2>(index_range) == boss.get_k()) {
        call_nodes_in_range(index_range);
        return;
    }

    std::vector<Range> stack { index_range };
    while (!stack.empty()) {
        auto range = stack.back();
        stack.pop_back();
        ++std::get<2>(range);

        bool keep_going = true;
        symbols(boss, range, [&](boss::BOSS::TAlphabet symbol) {
            if (!keep_going)
                return;
            auto next = range;
            auto &[first, last, length] = next;
            if (boss.tighten_range(&first, &last, symbol)) {
                if (length == boss.get_k())
                    keep_going = call_nodes_in_range(next);
                else
                    stack.push_back(next);
            }
        });
        if (!keep_going)
            return;
    }
}

} // namespace mtg::graph
