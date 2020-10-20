#include "row_diff.hpp"

#include <filesystem>

#include <progress_bar.hpp>

#include "annotation/binary_matrix/column_sparse/column_major.hpp"

namespace mtg {
namespace annot {
namespace binmat {

template <class BaseMatrix>
void RowDiff<BaseMatrix>::serialize(const std::string &filename) const {
    std::ofstream f(filename, ios::binary);
    serialize(f);
    f.close();
}

template <class BaseMatrix>
bool RowDiff<BaseMatrix>::load(const std::string &filename) {
    std::ifstream f(filename, ios::binary);
    bool result = load(f);
    f.close();

    load_terminal(anchors_filename_, &terminal_);

    return result;
}

template
class RowDiff<ColumnMajor>;

void build_successor(const std::string &graph_fname,
                     const std::string &outfbase,
                     uint32_t max_length,
                     uint32_t num_threads) {
    bool must_build = false;
    for (const auto &suffix : { ".succ", ".pred", ".pred_boundary", ".terminal" }) {
        if (!std::filesystem::exists(outfbase + suffix)) {
            common::logger->trace(
                    "Building and writing successor, predecessor and terminal files to {}.*",
                    outfbase);
            must_build = true;
            break;
        }
    }
    if (!must_build) {
        common::logger->trace("Using existing pred/succ/terminal files in {}.*", outfbase);
        return;
    }

    graph::DBGSuccinct graph(2);
    common::logger->trace("Loading graph...");
    if (!graph.load(graph_fname)) {
        common::logger->error("Cannot load graph from {}", graph_fname);
        std::exit(1);
    }

    using graph::boss::BOSS;
    const BOSS &boss = graph.get_boss();
    sdsl::bit_vector terminal;
    sdsl::bit_vector dummy;
    boss.call_sequences_row_diff([&](const vector<uint64_t> &, std::optional<uint64_t>) {},
                                 num_threads, max_length, &terminal, &dummy);

    // terminal uses BOSS edges as indices, so we need to map it to annotation indices
    sdsl::bit_vector term(graph.num_nodes(), 0);
    for (uint64_t i = 1; i < terminal.size(); ++i) {
        if (terminal[i]) {
            uint64_t graph_idx = graph.boss_to_kmer_index(i);
            uint64_t anno_index = graph::AnnotatedSequenceGraph::graph_to_anno_index(
                    graph_idx);
            assert(anno_index < graph.num_nodes());
            term[anno_index] = 1;
        }
    }

    std::ofstream fterm(outfbase + ".terminal.unopt", ios::binary);
    term.serialize(fterm);
    fterm.close();
    common::logger->trace("Anchor nodes written to {}.terminal", outfbase);

    // create the succ file, indexed using annotation indices
    uint32_t width = sdsl::bits::hi(graph.num_nodes()) + 1;
    sdsl::int_vector_buffer succ(outfbase + ".succ", std::ios::out, 1024 * 1024, width);
    sdsl::int_vector_buffer<1> pred_boundary(outfbase + ".pred_boundary", std::ios::out, 1024 * 1024);
    sdsl::int_vector_buffer pred(outfbase + ".pred", std::ios::out, 1024 * 1024, width);

    // traverse BOSS table in parallel processing num_threads chunks of size 1'000'000
    constexpr uint64_t chunk_size = 1'000'000;
    const uint64_t batch_size = chunk_size * num_threads;
    const uint64_t batch_count = (graph.num_nodes() - 1) / batch_size + 1;

    ProgressBar progress_bar(batch_count, "Compute successors", std::cerr,
                             !common::get_verbose());

    for (uint64_t batch = 0; batch < batch_count; ++batch) {
        std::vector<std::vector<uint64_t>> pred_buf(num_threads);
        std::vector<std::vector<uint64_t>> succ_buf(num_threads);
        std::vector<std::vector<bool>> pred_boundary_buf(num_threads);

        #pragma omp parallel for num_threads(num_threads) schedule(static, 1)
        for (uint32_t chunk = 0; chunk < std::max((uint32_t)1, num_threads); ++chunk) {
            uint64_t start = batch * batch_size + chunk * chunk_size + 1;
            for (uint64_t i = start;
                 i < std::min(graph.num_nodes() + 1, start + chunk_size); ++i) {
                uint64_t boss_idx = graph.kmer_to_boss_index(i);
                if (dummy[boss_idx] || terminal[boss_idx]) {
                    succ_buf[chunk].push_back(0);
                } else {
                    const graph::boss::BOSS::TAlphabet w
                            = boss.get_W(boss_idx) % boss.alph_size;
                    uint64_t next = w ? graph.boss_to_kmer_index(boss.fwd(boss_idx, w)) : 0;
                    succ_buf[chunk].push_back(next);
                }
                uint64_t back_idx = boss.bwd(boss_idx);
                boss.call_incoming_to_target(
                        back_idx, boss.get_W(back_idx), [&](BOSS::edge_index incoming_edge) {
                          // terminal and dummy predecessors are ignored; also ignore
                          // predecessors for which boss_idx is not the last outgoing
                          // edge (bc. we only traverse the last outgoing at a bifurcation)
                          if (terminal[incoming_edge] || dummy[incoming_edge]
                                  || boss.fwd(incoming_edge,
                                              boss.get_W(incoming_edge) % boss.alph_size)
                                          != boss_idx)
                              return;

                          pred_buf[chunk].push_back(graph.boss_to_kmer_index(incoming_edge));
                          pred_boundary_buf[chunk].push_back(0);
                        });
                pred_boundary_buf[chunk].push_back(1);
            }
        }
        for (uint32_t i = 0; i < num_threads; ++i) {
            for (uint32_t j = 0; j < succ_buf[i].size(); ++j) {
                succ.push_back(succ_buf[i][j]);
            }
            for (uint32_t j = 0; j < pred_buf[i].size(); ++j) {
                pred.push_back(pred_buf[i][j]);
            }
            for (uint32_t j = 0; j < pred_boundary_buf[i].size(); ++j) {
                pred_boundary.push_back(pred_boundary_buf[i][j]);
            }
        }
        ++progress_bar;
    }
    succ.close();
    pred.close();
    pred_boundary.close();

    common::logger->trace("Pred/succ nodes written to {}.pred/succ", outfbase);
}

} // namespace binmat
} // namespace annot
} // namespace mtg
