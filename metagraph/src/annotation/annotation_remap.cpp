#include "annotation_remap.hpp"

#include <atomic>
#include <mutex>
#include <stdexcept>

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/logger.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace annot {

using graph::DeBruijnGraph;
using node_index = DeBruijnGraph::node_index;
using mtg::common::logger;

static uint64_t width_for(uint64_t max_value) {
    return max_value ? sdsl::bits::hi(max_value) + 1 : 1;
}

sdsl::int_vector<> map_rows_to_graph(const DeBruijnGraph &source,
                                     const DeBruijnGraph &target,
                                     size_t num_threads) {
    if (source.get_mode() != DeBruijnGraph::BASIC
            || target.get_mode() != DeBruijnGraph::BASIC)
        throw std::invalid_argument("Row remapping supports only graphs in basic mode");

    if (source.get_k() != target.get_k())
        throw std::invalid_argument("Row remapping needs graphs with the same k");

    const uint64_t num_target_rows = target.max_index();
    sdsl::int_vector<> row_map(source.max_index(), 0, width_for(num_target_rows));
    sdsl::bit_vector found(num_target_rows, false);

    // Entries of row_map are packed, so concurrent writes are guarded by a lock
    // per block. A block of 2^12 entries spans a whole number of 64-bit words.
    constexpr uint64_t kBlockBits = 12;
    constexpr uint64_t kNumLocks = 1 << 12;
    std::unique_ptr<std::mutex[]> locks(new std::mutex[kNumLocks]);

    std::atomic<uint64_t> num_mapped(0);
    std::atomic<uint64_t> num_missing(0);
    std::atomic<uint64_t> num_repeated(0);

    source.call_sequences([&](const std::string &sequence, const std::vector<node_index> &path) {
        uint64_t mapped = 0;
        uint64_t missing = 0;
        uint64_t repeated = 0;
        size_t i = 0;
        target.map_to_nodes_sequentially(sequence, [&](node_index target_node) {
            node_index source_node = path.at(i++);
            if (!target_node) {
                ++missing;
                return;
            }
            if (fetch_and_set_bit(found.data(), target_node - 1, true,
                                  std::memory_order_relaxed)) {
                ++repeated;
                return;
            }
            uint64_t row = source_node - 1;
            std::lock_guard<std::mutex> lock(locks[(row >> kBlockBits) % kNumLocks]);
            row_map[row] = target_node;
            ++mapped;
        });
        num_mapped += mapped;
        num_missing += missing;
        num_repeated += repeated;
    }, num_threads);

    if (num_missing || num_repeated || num_mapped != target.num_nodes()) {
        throw std::runtime_error(fmt::format(
            "The graphs have different k-mers: {} source k-mers are missing from the"
            " target graph, {} target nodes were matched more than once, and {} of"
            " {} target nodes were matched",
            num_missing.load(), num_repeated.load(), num_mapped.load(), target.num_nodes()));
    }

    logger->trace("Matched {} of {} source rows to target rows",
                  num_mapped.load(), source.max_index());
    return row_map;
}

std::unique_ptr<MultiLabelAnnotation<std::string>>
remap_annotation_rows(const MultiLabelAnnotation<std::string> &annotation,
                      const DeBruijnGraph &source,
                      const sdsl::int_vector<> &row_map,
                      uint64_t target_num_rows,
                      size_t num_threads) {
    using matrix::BRWT;
    using matrix::RowDiff;

    if (annotation.get_matrix().num_rows() != row_map.size())
        throw std::invalid_argument("The row map does not match the annotation's rows");

    uint64_t max_target_row = 0;
    #pragma omp parallel for num_threads(std::max<size_t>(1, num_threads)) reduction(max:max_target_row)
    for (uint64_t row = 0; row < row_map.size(); ++row) {
        max_target_row = std::max<uint64_t>(max_target_row, row_map[row]);
    }
    if (max_target_row > target_num_rows)
        throw std::invalid_argument("The row map points past the end of the target graph");

    const auto &label_encoder = annotation.get_label_encoder();

    if (const auto *brwt = dynamic_cast<const MultiBRWTAnnotator *>(&annotation)) {
        auto matrix = brwt->get_matrix().remap_rows(row_map, target_num_rows, num_threads);
        return std::make_unique<MultiBRWTAnnotator>(std::move(matrix), label_encoder);
    }

    const auto *row_diff = dynamic_cast<const RowDiffBRWTAnnotator *>(&annotation);
    if (!row_diff)
        throw std::invalid_argument("Row remapping supports only brwt and row_diff_brwt annotations");

    const RowDiff<BRWT> &matrix = row_diff->get_matrix();
    const auto &anchor = matrix.anchor();
    const auto &fork_succ = matrix.fork_succ();

    if (anchor.size() != row_map.size())
        throw std::invalid_argument("The annotation has no anchors for every row");

    if (!fork_succ.size() && !dynamic_cast<const graph::DBGSuccinct *>(&source))
        throw std::invalid_argument("The annotation has no successor bitmap, which is"
                                    " needed for graphs other than succinct");

    // Keep each row's successor from the source graph. Successors are marked by
    // node index, so the bitmap has one extra entry.
    sdsl::bit_vector new_anchor(target_num_rows, false);
    sdsl::bit_vector new_succ(target_num_rows + 1, false);
    std::atomic<uint64_t> num_lost_successors(0);

    #pragma omp parallel for num_threads(std::max<size_t>(1, num_threads)) schedule(static, 1 << 20)
    for (uint64_t row = 0; row < row_map.size(); ++row) {
        uint64_t target_row = row_map[row];
        if (!target_row)
            continue;

        if (anchor[row]) {
            set_bit(new_anchor.data(), target_row - 1, true, std::memory_order_relaxed);
            continue;
        }

        node_index succ = row_diff_successor(source, row + 1, fork_succ);
        uint64_t target_succ = succ && succ <= row_map.size() ? row_map[succ - 1] : 0;
        if (!target_succ) {
            ++num_lost_successors;
            continue;
        }
        set_bit(new_succ.data(), target_succ, true, std::memory_order_relaxed);
    }

    if (num_lost_successors) {
        throw std::runtime_error(fmt::format(
            "{} rows that are not anchors have a dropped row-diff successor",
            num_lost_successors.load()));
    }

    auto diffs = matrix.diffs().remap_rows(row_map, target_num_rows, num_threads);
    logger->trace("Remapped row-diff matrix: {} of {} relations kept",
                  diffs->num_relations(), matrix.diffs().num_relations());

    auto result = std::make_unique<RowDiff<BRWT>>(nullptr, std::move(*diffs));
    result->set_anchor(matrix::IRowDiff::anchor_bv_type(std::move(new_anchor)));
    result->set_fork_succ(matrix::IRowDiff::fork_succ_bv_type(std::move(new_succ)));

    return std::make_unique<RowDiffBRWTAnnotator>(std::move(result), label_encoder);
}

} // namespace annot
} // namespace mtg
