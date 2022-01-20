#include "node_lcs.hpp"

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "common/vectors/vector_algorithm.hpp"

#include <filesystem>
#include <algorithm>
#include <mutex>

#include <sdsl/int_vector.hpp>
#include <sdsl/bits.hpp>
#include <progress_bar.hpp>

namespace mtg {
namespace graph {

static const uint64_t kBlockSize = 10'000;

NodeLCS::NodeLCS(const DeBruijnGraph &graph) {
    const DBGSuccinct *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);

    if (!dbg_succ) {
        common::logger->error("Only implemented for succinct graphs");
        exit(1);
    }

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

    lcs_ = std::make_unique<wavelet_tree_stat>(bits_per_lcs, std::move(lcs));
}

bool NodeLCS::load(const std::string &filename_base) {
    const auto lcs_filename = utils::make_suffix(filename_base, kLCSExtension);
    try {
        std::ifstream instream(lcs_filename, std::ios::binary);
        if (!instream.good())
            return false;

        lcs_ = std::make_unique<wavelet_tree_stat>(0);
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
    assert(width >= step);

    if (width == step)
        return std::make_pair(0, lcs_->size());

    width -= step;

    if (const auto *wt_sdsl = dynamic_cast<const wavelet_tree_stat*>(lcs_.get())) {
        const auto &wt = wt_sdsl->data();
        using size_type = typename std::remove_reference<decltype(wt)>::type::size_type;
        using value_type = typename std::remove_reference<decltype(wt)>::type::value_type;

        size_type num_cs = 0;
        std::vector<value_type> cs(wt.sigma);
        std::vector<size_type> rank_c_0(wt.sigma);
        std::vector<size_type> rank_c_i(wt.sigma);

        wt.interval_symbols(0, first + 1, num_cs, cs, rank_c_0, rank_c_i);
        first = 1;
        for (size_t i = 0; i < num_cs; ++i) {
            if (cs[i] < width)
                first = std::max(first, wt.select(rank_c_i[i], cs[i]));
        }

        wt.interval_symbols(last + 1, lcs_->size(), num_cs, cs, rank_c_0, rank_c_i);
        last = lcs_->size() - 1;
        for (size_t i = 0; i < num_cs; ++i) {
            if (cs[i] < width && rank_c_i[i] > rank_c_0[i])
                last = std::min(last, wt.select(rank_c_0[i] + 1, cs[i]) - 1);
        }

    } else {
        while (first && (*lcs_)[first] >= width) {
            --first;
        }

        while (last + 1 < lcs_->size() && (*lcs_)[last + 1] >= width) {
            ++last;
        }
    }

    assert(operator[](first) < width);
    assert(operator[](last + 1) < width);

#ifndef NDEBUG
    if (first < last) {
        for (uint8_t s = 0; s < width; ++s) {
            assert(lcs_->rank(s, first + 1) == lcs_->rank(s, last));
        }
    }
#endif

    return std::make_pair(first, last);
}

void NodeLCS::serialize(const std::string &filename_base) const {
    const auto fname = utils::make_suffix(filename_base, kLCSExtension);

    std::ofstream outstream(fname, std::ios::binary);
    lcs_->serialize(outstream);
}

void NodeLCS::serialize(sdsl::int_vector_buffer<>&& lcs,
                        const std::string &filename_base) {
    const auto fname = utils::make_suffix(filename_base, kLCSExtension);
    const std::string old_fname = lcs.filename();
    lcs.close(false); // close without removing the file
    std::filesystem::rename(old_fname, fname);
}

bool NodeLCS::is_compatible(const SequenceGraph &graph, bool verbose) const {
    // nodes plus dummy npos
    if (graph.max_index() + 1 == lcs_->size())
        return true;

    if (verbose)
        std::cerr << "ERROR: LCS file does not match number of nodes in graph"
                  << std::endl;
    return false;
}

} // namespace graph
} // namespace mtg
