#ifndef __UNITIGS_HPP__
#define __UNITIGS_HPP__

#include <progress_bar.hpp>

#include "cli/config/config.hpp"
#include "cli/load/load_graph.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/alignment/aligner_seeder_methods.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/int_matrix/rank_extended/csc_matrix.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/int_matrix/row_diff/int_row_diff.hpp"
#include "annotation/annotation_converters.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/unix_tools.hpp"


namespace mtg {
namespace graph {
namespace align {

class Unitigs : public SequenceGraph::GraphExtension {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef std::function<void(std::function<void(const std::vector<node_index>&)>)> PathGenerator;

    Unitigs(const DBGSuccinct &graph) : graph_(std::shared_ptr<const DeBruijnGraph>{}, &graph) {}

    Unitigs(const cli::Config &config, bool keep_graph = true) : graph_(load_graph_impl(config.fnames[0])) {
        std::vector<size_t> unitigs(graph_->num_nodes(), 0);
        sdsl::bit_vector indicator(unitigs.size(), false);

        std::atomic<size_t> i { 1 };
        std::atomic_thread_fence(std::memory_order_release);
        graph_->call_unitigs([&](const std::string&, const auto &path) {
            if (path.size() == 1)
                return;

            size_t idx = i.fetch_add(1, std::memory_order_relaxed);
            for (node_index node : path) {
                assert(node);

                unitigs[node - 1] = idx;
                set_bit(indicator.data(), node - 1, true, std::memory_order_relaxed);
            }
        }, get_num_threads());
        std::atomic_thread_fence(std::memory_order_acquire);

        common::logger->trace("Indexed {} unitigs from {} nodes", i - 1, graph_->num_nodes());
        graph_.reset();

        std::filesystem::path tmp_dir = utils::create_temp_dir(config.tmp_dir, "unitigs");
        std::string out_path = tmp_dir/"unitigs";
        std::vector<std::string> files;
        {
            annot::ColumnCompressed<> colcomp(
                std::move(indicator), "", 1, tmp_dir, 1e9, sdsl::bits::hi(i - 1) + 1
            );

            std::vector<std::string> labels;
            labels.push_back("");
            std::vector<annot::ColumnCompressed<>::Index> indices;
            indices.push_back(0);
            std::vector<uint64_t> counts;
            counts.push_back(0);
            ProgressBar progress_bar(unitigs.size(), "Packing unitig ids",
                                     std::cerr, !common::get_verbose());
            for (size_t i = 0; i < unitigs.size(); ++i) {
                if (unitigs[i]) {
                    indices.back() = i;
                    counts.back() = unitigs[i];
                    colcomp.add_label_counts(indices, labels, counts);
                }
                ++progress_bar;
            }
            files.push_back(out_path + colcomp.file_extension());
            colcomp.serialize(files[0]);
        }
        unitigs = std::vector<size_t>();
        {
            common::logger->trace("Compressing unitig index");
            common::logger->trace("Step 0");
            convert_to_row_diff(files,
                                config.fnames[0],
                                config.memory_available * 1e9,
                                config.max_path_length,
                                tmp_dir,
                                tmp_dir,
                                static_cast<annot::RowDiffStage>(0),
                                out_path + ".row_count", true);
            common::logger->trace("Step 1");
            convert_to_row_diff(files,
                                config.fnames[0],
                                config.memory_available * 1e9,
                                config.max_path_length,
                                tmp_dir,
                                tmp_dir,
                                static_cast<annot::RowDiffStage>(1),
                                out_path + ".row_reduction", true);
            common::logger->trace("Step 2");
            convert_to_row_diff(files,
                                config.fnames[0],
                                config.memory_available * 1e9,
                                config.max_path_length,
                                tmp_dir,
                                tmp_dir,
                                static_cast<annot::RowDiffStage>(2),
                                out_path + ".row_reduction", true);
            common::logger->trace("done");
            const std::string anchors_file = config.fnames[0] + annot::binmat::kRowDiffAnchorExt;
            if (!std::filesystem::exists(anchors_file)) {
                common::logger->error("Anchor bitmap {} does not exist.", anchors_file);
                std::exit(1);
            }
            const std::string fork_succ_file = config.fnames[0] + annot::binmat::kRowDiffForkSuccExt;
            if (!std::filesystem::exists(fork_succ_file)) {
                common::logger->error("Fork successor bitmap {} does not exist", fork_succ_file);
                std::exit(1);
            }

            common::logger->trace("Wrapping as RD");
            annot::ColumnCompressed<> colcomp;
            if (!colcomp.merge_load(files) || !colcomp.num_labels()) {
                common::logger->error("Failed to reload matrix from {}", files[0]);
                std::exit(1);
            }

            std::vector<sdsl::int_vector<>> column_values;
            annot::ColumnCompressed<>::load_column_values(files,
                [&](size_t, const std::string &, sdsl::int_vector<>&& values) {
                    column_values.emplace_back(std::move(values));
                }
            );
            if (column_values.size() != 1) {
                common::logger->error("Failed to reload values from {}", files[0]);
                std::exit(1);
            }

            auto colmap = colcomp.release_matrix();

            auto colmat = std::make_unique<RawUnitigs>(std::move(*colmap), std::move(column_values));
            unitigs_ = CompUnitigs(nullptr, std::move(*colmat));
            unitigs_.load_anchor(anchors_file);
            unitigs_.load_fork_succ(fork_succ_file);
        }

        if (keep_graph)
            load_graph(config.fnames[0]);
    }

    bool load(const std::string &filename_base) {
        std::string fname = utils::make_suffix(filename_base, kUnitigsExtension);
        std::ifstream fin(fname, std::ios::binary);
        if (!fin.good())
            return false;

        try {
            unitigs_.load(fin);
            unitigs_.set_graph(graph_.get());
            return true;
        } catch (...) {
            return false;
        }
    }

    void serialize(const std::string &filename_base) const {
        std::string fname = utils::make_suffix(filename_base, kUnitigsExtension);
        std::ofstream fout(fname, std::ios::binary);
        unitigs_.serialize(fout);
    }

    bool is_compatible(const SequenceGraph &, bool = true) const { return true; }

    std::vector<Seed> cluster_and_filter_seeds(const std::vector<Seed> &seeds,
                                               size_t min_seed_length) const {
        tsl::hopscotch_map<size_t, std::vector<size_t>> unitig_to_bucket;
        for (size_t i = 0; i < seeds.size(); ++i) {
            node_index node = seeds[i].get_nodes().back();
            if (node > graph_->max_index()) {
                node_index base_node = node - graph_->max_index();
                if (seeds[i].size() == 1 && graph_->get_boss().get_W(graph_->kmer_to_boss_index(base_node)) == boss::BOSS::kSentinel)
                    continue;

                node = base_node;
            }

            auto values = unitigs_.get_row_values(node - 1);
            if (values.size())
                unitig_to_bucket[values[0].second].emplace_back(i);
        }

        std::vector<Seed> filtered_seeds;
        for (auto it = unitig_to_bucket.begin(); it != unitig_to_bucket.end(); ++it) {
            if (it->second.size() > 1 || seeds[it->second[0]].get_query_view().size() > min_seed_length) {
                filtered_seeds.emplace_back(std::move(seeds[*std::min_element(
                    it.value().begin(), it.value().end(), [&seeds](size_t a, size_t b) {
                        return seeds[a].get_query_view().data() < seeds[b].get_query_view().data();
                    }
                )]));
            }
        }
        return filtered_seeds;
    }

    void load_graph(const std::string &fname) {
        graph_ = load_graph_impl(fname);
        unitigs_.set_graph(graph_.get());
    }

  private:
    typedef annot::matrix::CSCMatrix<annot::binmat::ColumnMajor> RawUnitigs;
    typedef annot::matrix::IntRowDiff<RawUnitigs> CompUnitigs;
    std::shared_ptr<const DBGSuccinct> graph_;
    CompUnitigs unitigs_;
    static constexpr auto kUnitigsExtension = ".unitigs";

    static std::shared_ptr<const DBGSuccinct> load_graph_impl(const std::string &fname) {
        common::logger->trace("Graph loading...");
        Timer timer;
        return std::dynamic_pointer_cast<const DBGSuccinct>(cli::load_critical_dbg(fname));
        common::logger->trace("Graph loaded in {} sec", timer.elapsed());
    }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __MER_DISTANCES_HPP__
