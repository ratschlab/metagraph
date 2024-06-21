#ifndef __ANNOTATED_GRAPH_ALGORITHM_HPP__
#define __ANNOTATED_GRAPH_ALGORITHM_HPP__

#include <tsl/hopscotch_set.h>

#include "graph/annotated_dbg.hpp"
#include "common/vectors/bit_vector.hpp"
#include "common/utils/file_utils.hpp"


namespace mtg {
namespace graph {


class MaskedDeBruijnGraph;

/**
 * A container for differential assembly parameters.
 */
struct DifferentialAssemblyConfig {
    bool count_kmers = false;
    bool clean = false;
    double family_wise_error_rate = 0.05;
    bool test_by_unitig = false;
    std::string test_type = "nbinom_exact";
    bool assemble_shared = false;
    uint64_t min_count = 1;
    uint64_t min_recurrence = 1;
    uint64_t min_in_recurrence = 0;
    uint64_t min_out_recurrence = 0;
    uint64_t max_in_recurrence = std::numeric_limits<uint64_t>::max();
    uint64_t max_out_recurrence = std::numeric_limits<uint64_t>::max();
    std::string outfbase;
};

template <class PValStorage>
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, PValStorage, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                         const std::vector<std::string> &files,
                         const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &labels_in,
                         const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &labels_out,
                         const DifferentialAssemblyConfig &config,
                         size_t num_threads = 1,
                         std::filesystem::path = "",
                         size_t num_parallel_files = std::numeric_limits<size_t>::max());

template <class ValuesContainer, class PValStorage>
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, PValStorage, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                         std::vector<std::unique_ptr<const bit_vector>> &columns_all,
                         std::vector<std::unique_ptr<const ValuesContainer>> &column_values_all,
                         const std::function<typename ValuesContainer::value_type(uint64_t /* row_i */, uint64_t /* col_j */)> &get_value,
                         const std::vector<bool> &groups,
                         const DifferentialAssemblyConfig &config,
                         size_t num_threads = 1,
                         std::filesystem::path = "",
                         size_t num_parallel_files = std::numeric_limits<size_t>::max(),
                         bool deallocate = false,
                         uint8_t max_width = 64);

} // namespace graph
} // namespace mtg

#endif // __ANNOTATED_GRAPH_ALGORITHM_HPP__
