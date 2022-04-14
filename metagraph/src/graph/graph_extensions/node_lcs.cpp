#include "node_lcs.hpp"

#include "common/vectors/vector_algorithm.hpp"

#include <filesystem>
#include <algorithm>
#include <mutex>

#include <sdsl/int_vector.hpp>
#include <sdsl/bits.hpp>
#include <progress_bar.hpp>

namespace mtg {
namespace graph {

using common::logger;
static const uint64_t kBlockSize = 10'000;

NodeLCS::NodeLCS(const DeBruijnGraph &graph) {
    const DBGSuccinct *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);

    if (!dbg_succ) {
        common::logger->error("Only implemented for succinct graphs");
        exit(1);
    }

    graph_ = dbg_succ;

    const boss::BOSS &boss = dbg_succ->get_boss();
    uint8_t bits_per_lcs = sdsl::bits::hi(boss.get_k()) + 1;
    auto lcs = aligned_int_vector(graph.max_index() + 1, 0, bits_per_lcs, 16);
    std::mutex backup_mu;

    ProgressBar progress_bar(lcs.size(), "Nodes traversed",
                             std::cerr, !common::get_verbose());

    // TODO: make this more efficient
    progress_bar += 2;
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t j = 2; j < lcs.size(); j += kBlockSize) {
        auto first = boss.get_node_seq(j - 1);

        size_t end = std::min(j + kBlockSize, lcs.size());
        for (size_t i = j; i < end; ++i) {
            auto second = boss.get_node_seq(i);
            uint64_t old = atomic_fetch_and_add(
                lcs, i,
                std::mismatch(first.rbegin(), first.rend(),
                              second.rbegin()).first - first.rbegin(),
                backup_mu, std::memory_order_relaxed
            );
            std::ignore = old;
            assert(!old);
            std::swap(first, second);
        }

        progress_bar += end - j;
    }

    lcs_ = std::make_unique<wavelet_tree_t>(bits_per_lcs, std::move(lcs));
}

bool NodeLCS::load(const std::string &filename_base) {
    const auto lcs_filename = utils::make_suffix(filename_base, kLCSExtension);
    try {
        std::ifstream instream(lcs_filename, std::ios::binary);
        if (!instream.good())
            return false;

        lcs_ = std::make_unique<wavelet_tree_t>(0);
        lcs_->load(instream);
        return true;

    } catch (...) {
        std::cerr << "ERROR: Cannot load graph LCS from file "
                  << lcs_filename << std::endl;
        return false;
    }
}

auto NodeLCS::expand(edge_index first, edge_index last, size_t width, size_t step) const
        -> std::pair<edge_index, edge_index> {
    assert(lcs_);
    assert(width >= step);

    if (width == step)
        return std::make_pair(0, lcs_->size());

    width -= step;

    while (first && (*lcs_)[first] >= width) {
        --first;
    }

    while (last + 1 < lcs_->size() && (*lcs_)[last + 1] >= width) {
        ++last;
    }

    assert(operator[](first) < width);
    assert(last + 1 == lcs_->size() || operator[](last + 1) < width);

#ifndef NDEBUG
    if (first < last) {
        for (uint8_t s = 0; s < width; ++s) {
            assert(lcs_->rank(s, first + 1) == lcs_->rank(s, last));
        }
    }
#endif

    return std::make_pair(first, last);
}

void NodeLCS
::call_contractions(edge_index first,
                    edge_index last,
                    size_t width,
                    const std::function<void(edge_index, edge_index)> &callback) const {
    assert(lcs_);
    edge_index cur = first;
    for (edge_index i = first + 1; i <= last; ++i) {
        if ((*lcs_)[i] <= width) {
            callback(cur, i - 1);
            cur = i;
        }
    }
    callback(cur, last);
}

auto NodeLCS::call_edges(const std::pair<edge_index, edge_index> &node_range,
                         char c,
                         const std::function<void(node_index, node_index)> &callback) const
        -> std::pair<edge_index, edge_index> {
    if (c == boss::BOSS::kSentinel)
        return {};

    const boss::BOSS &boss = graph_->get_boss();
    auto next_range = node_range;
    boss::BOSS::TAlphabet last_c = boss.encode(c);

    if (!boss.tighten_range(&next_range.first, &next_range.second, last_c))
        return {};

    assert(boss.rank_W(node_range.second, last_c)
            - boss.rank_W(node_range.first - 1, last_c)
                == boss.rank_last(next_range.second) - boss.rank_last(next_range.first - 1));

    edge_index node_wp = boss.succ_W(node_range.first, last_c);
    edge_index node_wm = boss.succ_W(node_range.first, last_c + boss.alph_size);
    edge_index i = boss.succ_last(next_range.first);
    edge_index i_start = boss.pred_last(i - 1) + 1;

    while (i <= next_range.second) {
        assert(node_wp != node_wm);
        edge_index node_w = std::min(node_wp, node_wm);
        assert(boss.fwd(node_w, last_c) == i);

        if (node_index n = graph_->boss_to_kmer_index(node_w)) {
            for (edge_index next_edge = i_start; next_edge <= i; ++next_edge) {
                if (node_index m = graph_->boss_to_kmer_index(next_edge))
                    callback(n, m);
            }
        }

        if (node_wp < node_wm) {
            node_wp = boss.succ_W(node_wp + 1, last_c);
        } else {
            node_wm = boss.succ_W(node_wm + 1, last_c + boss.alph_size);
        }

        if (node_wp <= node_wm) {
            i_start = i + 1;
            i = boss.succ_last(i_start);
        }
    }

    return next_range;
}

auto NodeLCS::call_edges(const std::pair<edge_index, edge_index> &node_range,
                         const std::function<void(node_index, const SmallVector<node_index>&)> &callback) const
        -> std::vector<std::pair<edge_index, edge_index>> {
    const boss::BOSS &boss = graph_->get_boss();

    std::vector<std::pair<edge_index, edge_index>> next_ranges(boss.alph_size, node_range);
    std::vector<std::pair<edge_index, edge_index>> next_is(boss.alph_size);
    for (boss::BOSS::TAlphabet s = 1; s < boss.alph_size; ++s) {
        auto &[first, second] = next_ranges[s];
        if (!boss.tighten_range(&first, &second, s)) {
            first = 0;
            second = 0;
        } else {
            assert(boss.rank_W(node_range.second, s)
                - boss.rank_W(node_range.first - 1, s)
                    == boss.rank_last(second) - boss.rank_last(first - 1));
            next_is[s].second = boss.succ_last(first);
            next_is[s].first = boss.pred_last(next_is[s].second - 1) + 1;
        }
    }

    for (edge_index node_w = node_range.first; node_w <= node_range.second; ++node_w) {
        boss::BOSS::TAlphabet last_c = boss.get_W(node_w) % boss.alph_size;
        if (!last_c)
            continue;

        auto &[i_start, i] = next_is[last_c];
        assert(boss.fwd(node_w, last_c) == i);
        assert(i <= next_ranges[last_c].second);

        if (node_index n = graph_->boss_to_kmer_index(node_w)) {
            SmallVector<node_index> edges(boss.alph_size);
            for (edge_index next_edge = i_start; next_edge <= i; ++next_edge) {
                if (node_index m = graph_->boss_to_kmer_index(next_edge))
                    edges[boss.get_W(next_edge) % boss.alph_size] = m;
            }
            callback(n, edges);
        }

        assert(node_w + 1 < boss.get_W().size() || node_w == node_range.second);

        if (node_w + 1 < boss.get_W().size()
                && boss.succ_W(node_w + 1, last_c + boss.alph_size)
                    >= boss.succ_W(node_w + 1, last_c)) {
            i_start = i + 1;
            i = boss.succ_last(i_start);
        }
    }

    return next_ranges;
}

auto NodeLCS
::call_parallel_edges(node_index node,
                      size_t suffix_length,
                      const std::function<void(node_index, node_index)> &callback) const
        -> std::pair<edge_index, edge_index> {
    const boss::BOSS &boss = graph_->get_boss();
    assert(boss.get_k() >= suffix_length);

    edge_index edge = graph_->kmer_to_boss_index(node);
    edge_index last_edge = boss.succ_last(edge);

    std::pair<edge_index, edge_index> node_range {
        boss.pred_last(last_edge - 1) + 1, last_edge
    };
    assert(node_range.first);
    assert(node_range.second);
    node_range = expand(node_range.first, node_range.second, boss.get_k(),
                        boss.get_k() - suffix_length);
    assert(operator[](node_range.first) < suffix_length);
    assert(!lcs_ || node_range.second + 1 == lcs_->size()
        || operator[](node_range.second + 1) < suffix_length);

    return call_edges(node_range, boss.decode(boss.get_W(edge) % boss.alph_size), callback);
}

void NodeLCS::serialize(const std::string &filename_base) const {
    if (!lcs_) {
        logger->error("Nothing to serialize");
        return;
    }

    const auto fname = utils::make_suffix(filename_base, kLCSExtension);

    std::ofstream outstream(fname, std::ios::binary);
    lcs_->serialize(outstream);
}

bool NodeLCS::is_compatible(const SequenceGraph &graph, bool) const {
    if (graph_ != dynamic_cast<const DBGSuccinct*>(&graph)) {
        logger->error("Graph does not match stored pointer");
        return false;
    }

    if (!lcs_)
        return true;

    // nodes plus dummy npos
    if (graph.max_index() + 1 == lcs_->size())
        return true;

    logger->error("LCS file does not match number of nodes in graph");
    return false;
}

void NodeLCS::set_graph(const DBGSuccinct &graph) {
    graph_ = &graph;
    assert(is_compatible(graph));
}

} // namespace graph
} // namespace mtg
