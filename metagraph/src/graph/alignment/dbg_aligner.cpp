#include "dbg_aligner.hpp"

#include "common/algorithms.hpp"
#include "graph/representation/rc_dbg.hpp"

namespace mtg {
namespace graph {
namespace align {


QueryAlignment IDBGAligner::align(std::string_view query,
                                  bool is_reverse_complement) const {
    QueryAlignment result(query);
    align_batch({ Query{ std::string{}, query, is_reverse_complement} },
        [&](std::string_view, QueryAlignment&& alignment) {
            result = std::move(alignment);
        }
    );

    return result;
}

template <class AlignmentCompare>
ISeedAndExtendAligner<AlignmentCompare>
::ISeedAndExtendAligner(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
      : graph_(graph), config_(config) {
    if (!config_.min_seed_length)
        config_.min_seed_length = graph_.get_k();

    if (!config_.max_seed_length)
        config_.max_seed_length = graph_.get_k();

    assert(config_.max_seed_length >= config_.min_seed_length);
    assert(config_.num_alternative_paths);
    assert(graph_.get_mode() != DeBruijnGraph::PRIMARY
        && "primary graphs must be wrapped into canonical");

    if (!config_.check_config_scores())
        throw std::runtime_error("Error: sum of min_cell_score and lowest penalty too low.");
}

template <class AlignmentCompare>
void ISeedAndExtendAligner<AlignmentCompare>
::align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
              const AlignmentCallback &callback) const {
    for (const auto &[header, query, is_reverse_complement] : seq_batch) {
        QueryAlignment paths(query, is_reverse_complement);
        Aggregator aggregator(paths.get_query(false), paths.get_query(true), config_);

        auto add_alignment = [&](Alignment&& alignment) {
            assert(alignment.is_valid(graph_, &config_));
            aggregator.add_alignment(std::move(alignment));
        };

        auto get_min_path_score = [&](const Alignment &seed) {
            return aggregator.get_min_path_score(seed);
        };

        std::string_view this_query = paths.get_query(is_reverse_complement);
        assert(this_query == query);

        std::vector<node_index> nodes = map_sequence_to_nodes(graph_, query);

        auto seeder = build_seeder(this_query, is_reverse_complement, nodes);
        auto extender = build_extender(this_query, aggregator);

        size_t num_explored_nodes = 0;

#if ! _PROTEIN_GRAPH
        if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                || config_.forward_and_reverse_complement) {
            std::vector<node_index> nodes_rc(nodes);
            std::string dummy(query);
            reverse_complement_seq_path(graph_, dummy, nodes_rc);
            assert(dummy == paths.get_query(!is_reverse_complement));
            assert(nodes_rc.size() == nodes.size());

            std::string_view reverse = paths.get_query(!is_reverse_complement);

            auto seeder_rc = build_seeder(reverse, !is_reverse_complement, nodes_rc);
            auto extender_rc = build_extender(reverse, aggregator);

            align_both_directions(paths.get_query(false), paths.get_query(true),
                                  *seeder, *seeder_rc, *extender, *extender_rc,
                                  add_alignment, get_min_path_score);

            num_explored_nodes += extender_rc->num_explored_nodes();

        } else {
            align_core(this_query, *seeder, *extender, add_alignment, get_min_path_score);
        }
#else
        align_core(this_query, *seeder, *extender, add_alignment, get_min_path_score);
#endif

        num_explored_nodes += extender->num_explored_nodes();

        aggregator.call_alignments([&](Alignment&& alignment) {
            assert(alignment.is_valid(graph_, &config_));
            paths.emplace_back(std::move(alignment));
        });

        common::logger->trace(
            "{}\tlength: {}\texplored nodes: {}\texplored nodes/k-mer: {}",
            header, query.size(), num_explored_nodes,
            static_cast<double>(num_explored_nodes) / nodes.size()
        );

        callback(header, std::move(paths));
    };
}

template <class AlignmentCompare>
void ISeedAndExtendAligner<AlignmentCompare>
::align_core(std::string_view query,
             const ISeeder &seeder,
             IExtender &extender,
             const std::function<void(Alignment&&)> &callback,
             const std::function<score_t(const Alignment&)> &get_min_path_score) const {
    for (Alignment &seed : seeder.get_seeds()) {
        if (seed.empty())
            continue;

        score_t min_path_score = get_min_path_score(seed);

        DEBUG_LOG("Min path score: {}\tSeed: {}", min_path_score, seed);

        auto extensions = extender.get_extensions(seed, min_path_score);

        if (extensions.empty() && seed.get_score() >= min_path_score) {
            seed.extend_query_end(query.data() + query.size());
            seed.trim_offset();
            DEBUG_LOG("Alignment (seed): {}", seed);
            callback(std::move(seed));
        }

        for (Alignment &extension : extensions) {
            DEBUG_LOG("Alignment (extension): {}", extension);
            callback(std::move(extension));
        }
    }
}

template <class AlignmentCompare>
void ISeedAndExtendAligner<AlignmentCompare>
::align_both_directions(std::string_view forward,
                        std::string_view reverse,
                        const ISeeder &forward_seeder,
                        const ISeeder &reverse_seeder,
                        IExtender &forward_extender,
                        IExtender &reverse_extender,
                        const std::function<void(Alignment&&)> &callback,
                        const std::function<score_t(const Alignment&)> &get_min_path_score) const {
#if _PROTEIN_GRAPH
    assert(false && "Only alignment in one direction supported for Protein graphs");
#endif

    RCDBG rc_dbg(graph_);
    bool use_rcdbg = graph_.get_mode() != DeBruijnGraph::CANONICAL
                        && config_.forward_and_reverse_complement;

    const DeBruijnGraph &rc_graph = use_rcdbg ? rc_dbg : graph_;

    auto is_reversible = [this](const Alignment &alignment) {
        return graph_.get_mode() == DeBruijnGraph::CANONICAL
            && alignment.get_orientation()
            && !alignment.get_offset();
    };

    auto get_forward_alignments = [&](std::string_view query,
                                      std::string_view query_rc,
                                      const ISeeder &seeder,
                                      IExtender &extender) {
        size_t farthest_reach = 0;
        score_t max_score = config_.min_cell_score;
        std::vector<Alignment> rc_of_alignments;

        DEBUG_LOG("Extending in forwards direction");
        align_core(query, seeder, extender,
            [&](Alignment&& path) {
                score_t min_path_score = get_min_path_score(path);

                farthest_reach = std::max(farthest_reach,
                                          path.get_query().size() + path.get_clipping());
                max_score = std::max(max_score, path.get_score());

                if (path.get_score() >= min_path_score) {
                    if (is_reversible(path)) {
                        Alignment out_path = path;
                        out_path.reverse_complement(graph_, query_rc);
                        assert(out_path.size());
                        callback(std::move(out_path));
                    } else {
                        callback(Alignment(path));
                    }
                }

                if (!path.get_clipping() || path.get_offset())
                    return;

                Alignment rev = path;
                rev.reverse_complement(rc_graph, query_rc);

                if (rev.empty()) {
                    DEBUG_LOG("This local alignment cannot be reversed, skipping");
                    return;
                }

                // Remove any character skipping from the end so that the
                // alignment can proceed
                assert(rev.get_end_clipping());
                rev.trim_end_clipping();

                assert(rev.is_valid(rc_graph, &config_));

                // Pass the reverse complement of the forward alignment
                // as a seed for extension
                rc_of_alignments.emplace_back(std::move(rev));
            },
            [&](const Alignment &seed) {
                return seed.get_clipping() <= farthest_reach
                    && config_.rel_score_cutoff > 0
                        ? max_score * config_.rel_score_cutoff
                        : config_.min_cell_score;
            }
        );

        std::sort(rc_of_alignments.begin(), rc_of_alignments.end(),
                  LocalAlignmentGreater());

        return rc_of_alignments;
    };

    ManualSeeder rc_of_reverse(get_forward_alignments(
        reverse, forward, reverse_seeder, reverse_extender
    ));

    ManualSeeder rc_of_forward(get_forward_alignments(
        forward, reverse, forward_seeder, forward_extender
    ));

    auto finish_alignment = [&](std::string_view query,
                                std::string_view query_rc,
                                const ManualSeeder &seeder,
                                IExtender &extender) {
        if (use_rcdbg)
            extender.set_graph(rc_dbg);

        align_core(query_rc, seeder, extender,
            [&](Alignment&& path) {
                if (use_rcdbg || is_reversible(path)) {
                    path.reverse_complement(rc_graph, query);
                    if (path.empty())
                        return;
                }

                assert(path.is_valid(graph_, &config_));
                callback(std::move(path));
            },
            get_min_path_score
        );
    };

    if (rc_of_forward.data().size() && (rc_of_reverse.data().empty()
            || rc_of_forward.data()[0].get_score() >= rc_of_reverse.data()[0].get_score())) {
        finish_alignment(forward, reverse, rc_of_forward, reverse_extender);
        finish_alignment(reverse, forward, rc_of_reverse, forward_extender);
    } else {
        finish_alignment(reverse, forward, rc_of_reverse, forward_extender);
        finish_alignment(forward, reverse, rc_of_forward, reverse_extender);
    }
}

template class ISeedAndExtendAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg
