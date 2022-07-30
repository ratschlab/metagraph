#ifndef __UNITIGS_HPP__
#define __UNITIGS_HPP__

#include <progress_bar.hpp>

#include "cli/config/config.hpp"
#include "cli/load/load_graph.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/int_matrix/rank_extended/csc_matrix.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/int_matrix/row_diff/int_row_diff.hpp"
#include "annotation/annotation_converters.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/unix_tools.hpp"
#include "common/vectors/bit_vector_dyn.hpp"


namespace mtg {
namespace graph {
namespace align {

class Unitigs : public SequenceGraph::GraphExtension {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef bit_vector_small Indicator;
    typedef annot::CountsVector IDVector;
    typedef annot::matrix::CSCMatrix<annot::binmat::ColumnMajor, IDVector> RawUnitigs;
    typedef annot::matrix::IntRowDiff<RawUnitigs> CompUnitigs;

    Unitigs(const DBGSuccinct &graph) : graph_(std::shared_ptr<const DeBruijnGraph>{}, &graph) {}

    Unitigs(const cli::Config &config) : graph_(load_graph_impl(config.fnames[0])) {
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
        common::logger->trace("Marking dummy k-mers");
        {
            DBGSuccinct &ncgraph = const_cast<DBGSuccinct&>(*graph_);
            ncgraph.mask_dummy_kmers(get_num_threads(), false);
            valid_edges_.reset(ncgraph.release_mask());
        }
        graph_.reset();

        std::filesystem::path tmp_dir = utils::create_temp_dir(config.tmp_dir, "unitigs");
        std::string out_path = tmp_dir/"unitigs";
        std::vector<std::string> files;
        {
            annot::LabelEncoder<> encoder;
            encoder.insert_and_encode("");
            std::vector<std::unique_ptr<bit_vector>> cols;
            cols.emplace_back(std::make_unique<Indicator>(std::move(indicator)));
            annot::ColumnCompressed<> colcomp(
                std::move(cols), std::move(encoder), 1, tmp_dir, 1e9, sdsl::bits::hi(i - 1) + 1
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

            common::logger->trace("Wrapping as IntRowDiff");
            annot::ColumnCompressed<> colcomp;
            if (!colcomp.merge_load(files) || !colcomp.num_labels()) {
                common::logger->error("Failed to reload matrix from {}", files[0]);
                std::exit(1);
            }

            std::vector<IDVector> column_values;
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
            load_graph(config.fnames[0]);
        }
    }

    bool load(const std::string &filename_base) {
        std::string fname = utils::make_suffix(filename_base, kUnitigsExtension);
        std::ifstream fin(fname, std::ios::binary);
        if (!fin.good())
            return false;

        if (!unitigs_.load(fin))
            return false;

        unitigs_.set_graph(graph_.get());
        switch (graph_->get_state()) {
            case boss::BOSS::State::STAT: {
                valid_edges_.reset(new bit_vector_small());
                break;
            }
            case boss::BOSS::State::FAST: {
                valid_edges_.reset(new bit_vector_stat());
                break;
            }
            case boss::BOSS::State::DYN: {
                valid_edges_.reset(new bit_vector_dyn());
                break;
            }
            case boss::BOSS::State::SMALL: {
                valid_edges_.reset(new bit_vector_small());
                break;
            }
        }
        return valid_edges_->load(fin);
    }

    void serialize(const std::string &filename_base) const {
        std::string fname = utils::make_suffix(filename_base, kUnitigsExtension);
        std::ofstream fout(fname, std::ios::binary);
        unitigs_.serialize(fout);
        valid_edges_->serialize(fout);
    }

    bool is_compatible(const SequenceGraph &, bool = true) const { return true; }

    size_t cluster_and_filter_seeds(IDBGAligner::BatchSeeders &batch_seeders,
                                    const DBGAlignerConfig &config) const {
        std::vector<bool> indicator;
        std::vector<node_index> nodes;
        size_t new_seed_count = 0;
        {
            ProgressBar progress_bar(batch_seeders.size(), "Extracting seed nodes",
                                     std::cerr, !common::get_verbose());
            for (const auto &[seeder, seeder_rc] : batch_seeders) {
                auto parse_seeder = [&](const auto &cur_seeder) {
                    for (const auto &seed : cur_seeder.get_seeds()) {
                        node_index node = seed.get_nodes().back();
                        if (node > graph_->max_index())
                            node = node - graph_->max_index();

                        if ((*valid_edges_)[graph_->kmer_to_boss_index(node)]) {
                            indicator.emplace_back(true);
                            nodes.emplace_back(node - 1);
                        } else {
                            indicator.emplace_back(false);
                        }
                    }
                };

                if (seeder)
                    parse_seeder(*seeder);

                if (seeder_rc)
                    parse_seeder(*seeder_rc);

                ++progress_bar;
            }
        }

        common::logger->trace("Fetching unitig IDs");
        auto seed_values = unitigs_.get_row_values(nodes);

        {
            ProgressBar progress_bar(batch_seeders.size(), "Clustering and filtering seeds",
                                     std::cerr, !common::get_verbose());
            size_t i = 0;
            for (auto &[seeder, seeder_rc] : batch_seeders) {
                size_t j = 0;
                auto parse_seeder = [&](auto &cur_seeder) {
                    tsl::hopscotch_map<size_t, std::vector<size_t>> unitig_to_bucket;
                    const auto &seeds = cur_seeder->get_seeds();
                    for (size_t k = 0; k < seeds.size(); ++k) {
                        if (!indicator[j++])
                            continue;

                        const auto &unitig_ids = seed_values[i++];
                        if (unitig_ids.size())
                            unitig_to_bucket[unitig_ids[0].second].emplace_back(k);
                    }

                    if (unitig_to_bucket.empty()) {
                        cur_seeder = std::make_shared<ManualMatchingSeeder>(std::vector<Seed>{}, 0, config);
                        return;
                    }

                    std::vector<Seed> filtered_seeds;
                    for (auto it = unitig_to_bucket.begin(); it != unitig_to_bucket.end(); ++it) {
                        if (it->second.size() > 1 || seeds[it->second[0]].get_query_view().size() > config.min_seed_length) {
                            filtered_seeds.emplace_back(seeds[*std::max_element(
                                it.value().begin(), it.value().end(), [&seeds](size_t a, size_t b) {
                                    return std::make_pair(seeds[a].get_query_view().size(),
                                                          seeds[b].get_query_view().data())
                                         < std::make_pair(seeds[b].get_query_view().size(),
                                                          seeds[a].get_query_view().data());
                                }
                            )]);
                        }
                    }

                    size_t num_matches = get_num_char_matches_in_seeds(filtered_seeds.begin(), filtered_seeds.end());
                    new_seed_count += filtered_seeds.size();
                    cur_seeder = std::make_shared<ManualMatchingSeeder>(
                        std::move(filtered_seeds), num_matches, config
                    );
                };

                if (seeder)
                    parse_seeder(seeder);

                if (seeder_rc)
                    parse_seeder(seeder_rc);

                ++progress_bar;
            }
        }

        return new_seed_count;
    }

    void load_graph(const std::string &fname) {
        graph_ = load_graph_impl(fname);
        unitigs_.set_graph(graph_.get());
    }

  private:
    std::shared_ptr<const DBGSuccinct> graph_;
    CompUnitigs unitigs_;
    std::unique_ptr<bit_vector> valid_edges_;
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
