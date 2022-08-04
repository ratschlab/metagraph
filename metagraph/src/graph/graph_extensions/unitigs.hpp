#ifndef __UNITIGS_HPP__
#define __UNITIGS_HPP__

#include <mutex>

#include <progress_bar.hpp>

#include "cli/config/config.hpp"
#include "cli/load/load_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/int_matrix/row_diff/tuple_row_diff.hpp"
#include "annotation/annotation_converters.hpp"
#include "common/utils/file_utils.hpp"
#include "common/unix_tools.hpp"
#include "common/vectors/bit_vector_dyn.hpp"
#include "graph/alignment/aligner_labeled.hpp"


namespace mtg {
namespace graph {
namespace align {

class Unitigs : public SequenceGraph::GraphExtension {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef bit_vector_small Indicator;
    typedef annot::CountsVector IDVector;
    typedef annot::ColumnCoordAnnotator::binary_matrix_type RawUnitigs;
    typedef annot::matrix::TupleRowDiff<RawUnitigs> CompUnitigs;

    Unitigs(const DBGSuccinct &graph) : graph_(std::shared_ptr<const DeBruijnGraph>{}, &graph) {}

    Unitigs(const cli::Config &config) : graph_(load_graph_impl(config.fnames[0])) {
        std::filesystem::path tmp_dir = utils::create_temp_dir(config.tmp_dir, "unitigs");
        std::string out_path = tmp_dir/"unitigs";
        std::vector<std::string> files;
        size_t width = sdsl::bits::hi(graph_->num_nodes()) + 1;
        std::vector<size_t> unitigs;
        std::vector<std::unique_ptr<bit_vector>> cols;

        common::logger->trace("Marking dummy k-mers");
        {
            DBGSuccinct &ncgraph = const_cast<DBGSuccinct&>(*graph_);
            ncgraph.mask_dummy_kmers(get_num_threads(), false);
            valid_edges_.reset(ncgraph.release_mask());
            ProgressBar progress_bar(valid_edges_->num_set_bits(), "Transforming mask",
                                     std::cerr, !common::get_verbose());
            cols.emplace_back(std::make_unique<bit_vector_smart>([&](const auto &callback) {
                valid_edges_->call_ones([&](node_index i) {
                    assert(i);
                    callback(i - 1);
                    ++progress_bar;
                });
            }, ncgraph.num_nodes(), valid_edges_->num_set_bits()));
        }

        annot::LabelEncoder<> encoder;
        encoder.insert_and_encode("");
        auto colcomp = std::make_unique<annot::ColumnCompressed<>>(
            std::move(cols), std::move(encoder), 1, tmp_dir, 1e9, width
        );

        std::vector<std::string> labels;
        labels.push_back("");

        common::logger->trace("Annotating unitigs");
        std::mutex mu;
        size_t counter = 0;
        size_t max_unitig = 0;
        graph_->call_unitigs([&](const std::string&, const auto &path) {
            if (path.size() == 1)
                return;

#ifndef NDEBUG
            if (graph_->has_single_outgoing(path.back())
                    && graph_->has_single_incoming(path.front())) {
                graph_->adjacent_outgoing_nodes(path.back(), [&](node_index next) {
                    assert(next == path.front() || !(*valid_edges_)[next]
                            || graph_->indegree(next) > 1);
                });
                graph_->adjacent_incoming_nodes(path.front(), [&](node_index prev) {
                    assert(prev == path.back() || !(*valid_edges_)[prev]
                            || graph_->outdegree(prev) > 1);
                });
            }
#endif

            std::lock_guard<std::mutex> lock(mu);
            unitigs.emplace_back(path.front());
            unitigs.emplace_back(path.back());
            unitigs.emplace_back(counter);
            max_unitig = std::max({ path.front(), path.back(), max_unitig });
            counter += path.size();
            std::vector<std::pair<annot::ColumnCompressed<>::Index, uint64_t>> coords;
            for (size_t i = 0; i < path.size(); ++i) {
                coords.emplace_back(path[i] - 1, i + unitigs.back());
            }

            colcomp->add_label_coords(coords, labels);
        });

        common::logger->trace("Initializing unitig vector");
        size_t num_unitigs = unitigs.size() / 3;
        // TODO: replace with int_vector_buffer
        sdsl::int_vector<> boundaries(num_unitigs * 2, 0, sdsl::bits::hi(max_unitig + 1) + 1);
        indicator_ = Indicator([&](const auto &callback) {
            ProgressBar progress_bar(num_unitigs, "Packing unitigs",
                                     std::cerr, !common::get_verbose());
            for (size_t i = 0, j = 0; i < unitigs.size(); i += 3, j += 2) {
                boundaries[j] = unitigs[i];
                boundaries[j + 1] = unitigs[i + 1];
                callback(unitigs[i + 2]);
                ++progress_bar;
            }
        }, counter, num_unitigs);
        boundaries_ = IDVector(std::move(boundaries));
        unitigs = std::vector<size_t>();

        common::logger->trace("Serializing initial annotation");
        files.push_back(out_path + colcomp->file_extension());
        colcomp->serialize(files[0]);
        colcomp.reset();
        graph_.reset();
        common::logger->trace("Compressing unitig index");
        common::logger->trace("Step 0");
        convert_to_row_diff(files,
                            config.fnames[0],
                            config.memory_available * 1e9,
                            config.max_path_length,
                            tmp_dir,
                            tmp_dir,
                            static_cast<annot::RowDiffStage>(0),
                            out_path + ".row_count", false, true);
        common::logger->trace("Step 1");
        convert_to_row_diff(files,
                            config.fnames[0],
                            config.memory_available * 1e9,
                            config.max_path_length,
                            tmp_dir,
                            tmp_dir,
                            static_cast<annot::RowDiffStage>(1),
                            out_path + ".row_reduction", false, true);
        common::logger->trace("Step 2");
        convert_to_row_diff(files,
                            config.fnames[0],
                            config.memory_available * 1e9,
                            config.max_path_length,
                            tmp_dir,
                            tmp_dir,
                            static_cast<annot::RowDiffStage>(2),
                            out_path + ".row_reduction", false, true);
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

        common::logger->trace("Loading column");
        auto annotator = std::make_unique<annot::ColumnCompressed<>>(0);
        if (!annotator->merge_load(files)) {
            common::logger->error("Cannot load annotations");
            exit(1);
        }

        common::logger->trace("Wrapping as TupleRowDiff");

        std::vector<bit_vector_smart> delimiters(1);
        std::vector<sdsl::int_vector<>> column_values(1);

        auto coords_fname = utils::remove_suffix(files[0], annot::ColumnCompressed<>::kExtension)
                                                        + annot::ColumnCompressed<>::kCoordExtension;
        std::ifstream in(coords_fname, std::ios::binary);
        try {
            RawUnitigs::load_tuples(in, 1, [&](auto&& delims, auto&& values) {
                delimiters[0] = std::move(delims);
                column_values[0] = std::move(values);
            });
        } catch (const std::exception &e) {
            common::logger->error("Couldn't load coordinates from {}\nException: {}", coords_fname, e.what());
            exit(1);
        } catch (...) {
            common::logger->error("Couldn't load coordinates from {}", coords_fname);
            exit(1);
        }

        unitigs_ = CompUnitigs(nullptr,
                               RawUnitigs(std::move(*annotator->release_matrix()),
                                          std::move(delimiters),
                                          std::move(column_values)));

        unitigs_.load_anchor(anchors_file);
        unitigs_.load_fork_succ(fork_succ_file);
        common::logger->trace("RowDiff support bitmaps loaded");
        load_graph(config.fnames[0]);
        unitigs_.set_graph(graph_.get());
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

        if (!valid_edges_->load(fin))
            return false;

        try {
            boundaries_.load(fin);
            return true;
        } catch (...) {
            return false;
        }

        return indicator_.load(fin);
    }

    void serialize(const std::string &filename_base) const {
        std::string fname = utils::make_suffix(filename_base, kUnitigsExtension);
        std::ofstream fout(fname, std::ios::binary);
        unitigs_.serialize(fout);
        valid_edges_->serialize(fout);
        boundaries_.serialize(fout);
        indicator_.serialize(fout);
    }

    bool is_compatible(const SequenceGraph &, bool = true) const { return true; }

    void adjacent_outgoing_nodes(size_t unitig_id,
                                 const std::function<void(size_t)> &callback) const {
        graph_->adjacent_outgoing_nodes(get_unitig(unitig_id).second, [&](node_index next) {
            auto next_unitig_ids = get_unitig_ids({ next });
            callback(next_unitig_ids.size() ? next_unitig_ids[0] : next);
        });
    }

    void adjacent_incoming_nodes(size_t unitig_id,
                                 const std::function<void(size_t)> &callback) const {
        graph_->adjacent_incoming_nodes(get_unitig(unitig_id).second, [&](node_index prev) {
            auto prev_unitig_ids = get_unitig_ids({ prev });
            callback(prev_unitig_ids.size() ? prev_unitig_ids[0] : prev);
        });
    }

    std::pair<node_index, node_index> get_unitig(size_t unitig_id) const {
        size_t unitig_id_offset = get_unitig_id_offset();
        if (unitig_id <= unitig_id_offset)
            return std::make_pair(unitig_id, unitig_id);

        unitig_id -= unitig_id_offset;
        return std::make_pair(boundaries_[(unitig_id - 1) * 2],
                              boundaries_[(unitig_id - 1) * 2 + 1]);
    }

    std::pair<std::pair<node_index, node_index>, std::pair<size_t, size_t>>
    get_unitig_bounds(size_t unitig_id) const {
        size_t unitig_id_offset = get_unitig_id_offset();
        if (unitig_id <= unitig_id_offset) {
            return std::make_pair(std::make_pair(unitig_id, unitig_id),
                                  std::make_pair(0, 1));
        }

        unitig_id -= unitig_id_offset;
        return std::make_pair(std::make_pair(boundaries_[(unitig_id - 1) * 2],
                                             boundaries_[(unitig_id - 1) * 2 + 1]),
                              std::make_pair(indicator_.select1(unitig_id),
                                             indicator_.select1(unitig_id + 1))
        );
    }

    std::vector<size_t> get_unitig_ids(const std::vector<node_index> &nodes) const {
        auto [indicator, rows] = nodes_to_rows(nodes);
        common::logger->trace("Fetching unitig IDs");
        auto seed_tuples = unitigs_.get_row_tuples(rows);
        std::vector<size_t> results;
        results.reserve(nodes.size());

        size_t j = 0;
        size_t unitig_id_offset = get_unitig_id_offset();
        for (size_t i = 0; i < nodes.size(); ++i) {
            size_t unitig_id = std::numeric_limits<size_t>::max();
            if (indicator[i]) {
                const auto &coordinates = seed_tuples[j++];
                if (coordinates.size()) {
                    assert(coordinates.size() == 1);
                    assert(!coordinates[0].first);
                    if (coordinates[0].second.size()) {
                        assert(coordinates[0].second.size() == 1);
                        unitig_id = indicator_.rank1(coordinates[0].second[0]);
                        results.emplace_back(unitig_id + unitig_id_offset);
                    }
                }
            }

            if (unitig_id == std::numeric_limits<size_t>::max()) {
                results.emplace_back(nodes[i] > graph_->max_index()
                                         ? nodes[i] - graph_->max_index()
                                         : nodes[i]);
            }
        }

        return results;
    }

    std::vector<std::pair<size_t, size_t>>
    get_unitig_ids_and_coordinates(const std::vector<node_index> &nodes) const {
        auto [indicator, rows] = nodes_to_rows(nodes);
        common::logger->trace("Fetching unitig IDs");
        auto seed_tuples = unitigs_.get_row_tuples(rows);
        std::vector<std::pair<size_t, size_t>> results;
        results.reserve(nodes.size());

        size_t j = 0;
        size_t unitig_id_offset = graph_->max_index() + 1;
        for (size_t i = 0; i < nodes.size(); ++i) {
            size_t unitig_id = std::numeric_limits<size_t>::max();
            if (indicator[i]) {
                const auto &coordinates = seed_tuples[j++];
                if (coordinates.size()) {
                    assert(coordinates.size() == 1);
                    assert(!coordinates[0].first);
                    if (coordinates[0].second.size()) {
                        assert(coordinates[0].second.size() == 1);
                        unitig_id = indicator_.rank1(coordinates[0].second[0]);
                        results.emplace_back(
                            unitig_id + unitig_id_offset,
                            indicator_.select1(unitig_id) - coordinates[0].second[0]
                        );
                    }
                }
            }

            if (unitig_id == std::numeric_limits<size_t>::max()) {
                results.emplace_back(nodes[i] > graph_->max_index()
                                         ? nodes[i] - graph_->max_index()
                                         : nodes[i],
                                     0);
            }
        }

        return results;
    }

    size_t cluster_and_filter_seeds(const IDBGAligner &aligner,
                                    IDBGAligner::BatchSeeders &batch_seeders) const {
        const DBGAlignerConfig &config = aligner.get_config();
        std::vector<bool> indicator;
        std::vector<node_index> nodes;
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

        size_t new_seed_count = 0;
        common::logger->trace("Fetching unitig IDs");
        auto seed_tuples = unitigs_.get_row_tuples(nodes);

        {
            ProgressBar progress_bar(batch_seeders.size(), "Clustering and filtering seeds",
                                     std::cerr, !common::get_verbose());
            size_t i = 0;
            for (auto &[seeder, seeder_rc] : batch_seeders) {
                size_t j = 0;
                auto parse_seeder = [&](auto &cur_seeder) {
                    tsl::hopscotch_map<size_t, std::vector<std::pair<size_t, size_t>>> unitig_to_bucket;
                    const auto &seeds = cur_seeder->get_seeds();
                    for (size_t k = 0; k < seeds.size(); ++k) {
                        if (!indicator[j++])
                            continue;

                        const auto &coordinates = seed_tuples[i];
                        size_t unitig_id = std::numeric_limits<size_t>::max();
                        if (coordinates.size()) {
                            assert(coordinates.size() == 1);
                            assert(!coordinates[0].first);
                            if (coordinates[0].second.size()) {
                                assert(coordinates[0].second.size() == 1);
                                unitig_id = indicator_.rank1(coordinates[0].second[0]);
                                unitig_to_bucket[unitig_id].emplace_back(k, coordinates[0].second[0]);
                            }
                        }
                        ++i;
                    }

                    if (unitig_to_bucket.empty()) {
                        cur_seeder = std::make_shared<ManualMatchingSeeder>(std::vector<Seed>{}, 0, config);
                        return;
                    }

                    std::vector<Seed> filtered_seeds;
                    DEBUG_LOG("Found {} unitigs", unitig_to_bucket.size());
                    for (auto it = unitig_to_bucket.begin(); it != unitig_to_bucket.end(); ++it) {
                        auto &bucket = it.value();
                        if (bucket.size() == 1 && seeds[bucket[0].first].get_query_view().size() > config.min_seed_length) {
                            filtered_seeds.emplace_back(seeds[bucket[0].first]);
                        } else if (bucket.size() >= 2) {
                            // look for co-linear seeds
                            DEBUG_LOG("Chaining bucket of size {}", bucket.size());
                            std::string_view forward;
                            std::string_view reverse;
                            if (seeds[bucket[0].first].get_orientation()) {
                                reverse = std::string_view(seeds[bucket[0].first].get_query_view().data() - seeds[bucket[0].first].get_clipping(),
                                                           seeds[bucket[0].first].get_query_view().size() + seeds[bucket[0].first].get_clipping() + seeds[bucket[0].first].get_end_clipping());
                            } else {
                                forward = std::string_view(seeds[bucket[0].first].get_query_view().data() - seeds[bucket[0].first].get_clipping(),
                                                           seeds[bucket[0].first].get_query_view().size() + seeds[bucket[0].first].get_clipping() + seeds[bucket[0].first].get_end_clipping());
                            }
                            std::vector<Seed> seed_bucket_fwd;
                            std::vector<Seed> seed_bucket_bwd;
                            for (const auto &[k, coord] : bucket) {
                                seed_bucket_fwd.emplace_back(seeds[k]);
                                const auto &columns = seeds[k].get_columns();
                                seed_bucket_fwd.back().label_coordinates.resize(
                                    columns.size(),
                                    Alignment::Tuple(1, coord - seeds[k].size() + 1 + seeds[k].get_offset())
                                );
                                DEBUG_LOG("\t{}", Alignment(seed_bucket_fwd.back(), config));
                            }
                            if (reverse.data())
                                std::swap(seed_bucket_fwd, seed_bucket_bwd);

                            tsl::hopscotch_set<Alignment::Column> used_labels;
                            size_t old_filtered_size = filtered_seeds.size();
                            call_seed_chains_both_strands(aligner,
                                                          forward, reverse,
                                                          config,
                                                          std::move(seed_bucket_fwd),
                                                          std::move(seed_bucket_bwd),
                                                          [&](Chain&& chain, score_t) {
                                if (chain.size() > 1 || chain[0].first.get_query_view().size() > config.min_seed_length) {
                                    filtered_seeds.emplace_back(
                                        chain[0].first.get_query_view(),
                                        std::vector<node_index>(chain[0].first.get_nodes()),
                                        chain[0].first.get_orientation(),
                                        chain[0].first.get_offset(),
                                        chain[0].first.get_clipping(),
                                        chain[0].first.get_end_clipping()
                                    );
                                    filtered_seeds.back().label_columns = chain[0].first.label_columns;
                                    filtered_seeds.back().label_encoder = chain[0].first.label_encoder;
                                    const auto &columns = filtered_seeds.back().get_columns();
                                    used_labels.insert(columns.begin(), columns.end());
                                }

                            }, [&](Alignment::Column c) -> bool { return used_labels.count(c); });
                            std::sort(filtered_seeds.begin() + old_filtered_size, filtered_seeds.end(),
                                      [](const auto &a, const auto &b) {
                                          return std::make_pair(a.get_query_view().data(), a.get_query_view().size())
                                            < std::make_pair(b.get_query_view().data(), b.get_query_view().size());
                                      });

                            for (size_t i = old_filtered_size + 1; i < filtered_seeds.size(); ++i) {
                                if (filtered_seeds[i] == filtered_seeds[i - 1]) {
                                    if (filtered_seeds[i].label_encoder) {
                                        const auto &a = filtered_seeds[i - 1].get_columns();
                                        const auto &b = filtered_seeds[i].get_columns();
                                        assert(filtered_seeds[i - 1].label_encoder->check_node_labels_is_superset(a, filtered_seeds[i - 1].get_nodes()));
                                        assert(filtered_seeds[i].label_encoder->check_node_labels_is_superset(b, filtered_seeds[i].get_nodes()));
                                        Vector<Alignment::Column> merged;
                                        merged.reserve(a.size() + b.size());
                                        std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                                                       std::back_inserter(merged));
                                        assert(merged.empty()
                                            || filtered_seeds[i].label_encoder->check_node_labels_is_superset(merged, filtered_seeds[i].get_nodes()));

                                        filtered_seeds[i].label_columns
                                            = filtered_seeds[i].label_encoder->cache_column_set(std::move(merged));
                                    }
                                    filtered_seeds[i - 1] = Seed();
                                }
                            }
                            auto it = std::remove_if(filtered_seeds.begin() + old_filtered_size,
                                                     filtered_seeds.end(),
                                                     [](const auto &a) { return a.empty(); });
                            filtered_seeds.erase(it, filtered_seeds.end());

                        } else {
                            // TODO: check if there are any neighboring buckets
                        }
                    }

                    DEBUG_LOG("Added {} seeds", filtered_seeds.size());
#ifndef NDEBUG
                    for (const auto &seed : filtered_seeds) {
                        DEBUG_LOG("\t{}", Alignment(seed, config));
                    }
#endif
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
    IDVector boundaries_;
    Indicator indicator_;
    std::unique_ptr<bit_vector> valid_edges_;
    static constexpr auto kUnitigsExtension = ".unitigs";

    static std::shared_ptr<const DBGSuccinct> load_graph_impl(const std::string &fname) {
        common::logger->trace("Graph loading...");
        Timer timer;
        return std::dynamic_pointer_cast<const DBGSuccinct>(cli::load_critical_dbg(fname));
        common::logger->trace("Graph loaded in {} sec", timer.elapsed());
    }

    std::pair<sdsl::bit_vector, std::vector<uint64_t>>
    nodes_to_rows(const std::vector<node_index> &nodes) const {
        sdsl::bit_vector indicator(nodes.size(), false);
        std::vector<uint64_t> rows;
        for (size_t i = 0; i < nodes.size(); ++i) {
            node_index node = nodes[i];
            if (node > graph_->max_index())
                node -= graph_->max_index();

            if ((*valid_edges_)[graph_->kmer_to_boss_index(node)]) {
                indicator[i] = true;
                rows.emplace_back(node - 1);
            }
        }

        return std::make_pair(std::move(indicator), std::move(rows));
    }

    size_t get_unitig_id_offset() const { return graph_->max_index() + 1; }
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __MER_DISTANCES_HPP__
