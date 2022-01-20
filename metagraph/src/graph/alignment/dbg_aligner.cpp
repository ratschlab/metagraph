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

/**
 * Partition the alignment at the last k-mer. Return a pair containing the
 * alignment of all but the last k-mers, and the alignment of the last k-mer.
 */
std::pair<Alignment, Alignment> split_seed(const DeBruijnGraph &graph,
                                           const DBGAlignerConfig &config,
                                           const Alignment &alignment) {
    if (alignment.get_sequence().size() < graph.get_k() * 2)
        return std::make_pair(Alignment(), alignment);

    auto ret_val = std::make_pair(alignment, alignment);
    ret_val.first.trim_reference_suffix(graph.get_k(), graph, config, false);
    ret_val.second.trim_reference_prefix(alignment.get_sequence().size() - graph.get_k(),
                                         graph, config, true);

    return ret_val;
}

/**
 * filter_seed: clear the seed a if it has no unexplored labels or coordinates
 *              relative to the seed prev
 */
struct CoordDiff {
    CoordDiff(size_t offset = 0) : offset_(offset) {}

    template <typename It1, typename It2, typename Out>
    void operator()(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end, Out out) const {
        while (a_begin != a_end) {
            if (b_begin == b_end || *a_begin + offset_ < *b_begin) {
                *out = *a_begin;
                ++out;
                ++a_begin;
            } else if (*a_begin + offset_ > *b_begin) {
                ++b_begin;
            } else {
                ++a_begin;
                ++b_begin;
            }
        }
    }

    size_t offset_;
};

void filter_seed(const Alignment &prev, Alignment &a) {
    if (prev.label_columns.empty()) {
        a = Alignment();
    } else if (prev.label_coordinates.empty()) {
        Vector<Alignment::Column> diff;
        std::set_difference(a.label_columns.begin(),
                            a.label_columns.end(),
                            prev.label_columns.begin(),
                            prev.label_columns.end(),
                            std::back_inserter(diff));
        if (diff.empty()) {
            a = Alignment();
        } else {
            std::swap(a.label_columns, diff);
        }
    } else {
        Vector<Alignment::Column> diff;
        Vector<Alignment::Tuple> diff_coords;
        utils::indexed_set_op<Alignment::Tuple, CoordDiff>(
            a.label_columns.begin(),
            a.label_columns.end(),
            a.label_coordinates.begin(),
            prev.label_columns.begin(),
            prev.label_columns.end(),
            prev.label_coordinates.begin(),
            std::back_inserter(diff), std::back_inserter(diff_coords)
        );
        if (diff.empty()) {
            a = Alignment();
        } else {
            std::swap(a.label_columns, diff);
            std::swap(a.label_coordinates, diff_coords);
        }
    }

}

// Compute the maximum coordinate distance between two alignments
struct AlignmentPairedCoordinatesDist {
    size_t operator()(const Alignment &a, const Alignment &b) const {
        size_t max_dist = 0;
        if (a.label_coordinates.empty() || b.label_coordinates.empty())
            return b.get_query().data() + b.get_query().size() - a.get_query().data();

        auto a_begin = a.label_columns.begin();
        auto a_end = a.label_columns.end();
        auto b_begin = b.label_columns.begin();
        auto b_end = b.label_columns.end();
        auto a_cs_begin = a.label_coordinates.begin();
        auto b_cs_begin = b.label_coordinates.begin();
        size_t len = b.get_sequence().size();
        while (a_begin != a_end && b_begin != b_end) {
            if (*a_begin < *b_begin) {
                ++a_begin;
                ++a_cs_begin;
            } else if (*a_begin > *b_begin) {
                ++b_begin;
                ++b_cs_begin;
            } else {
                assert(a_cs_begin->size() == b_cs_begin->size());
                for (size_t i = 0; i < a_cs_begin->size(); ++i) {
                    max_dist = std::max(
                        max_dist,
                        static_cast<size_t>((*b_cs_begin)[i] + len - (*a_cs_begin)[i])
                    );
                }
                ++a_begin;
                ++b_begin;
                ++a_cs_begin;
                ++b_cs_begin;
            }
        }

        return max_dist;
    }
};

// Extend the alignment first until it reaches the end of the alignment second
template <class BuildExtender>
size_t align_connect(std::string_view query,
                     const DeBruijnGraph &graph,
                     DBGAlignerConfig &config,
                     Alignment &first,
                     Alignment &second,
                     const BuildExtender &build_extender) {
    auto [left, next] = split_seed(graph, config, first);
    config.xdrop += first.get_score() - next.get_score();
    config.target_distance = AlignmentPairedCoordinatesDist()(next, second);
    std::string_view query_window(
        query.data(),
        second.get_query().data() + second.get_query().size() - query.data()
    );
    assert(next.get_query()
        == query_window.substr(next.get_clipping(), next.get_query().size()));
    auto extender = build_extender(query_window, config);
    auto extensions = extender->get_extensions(next, config.ninf, true);
    if (extensions.size() && extensions[0].get_query().data() + extensions[0].get_query().size()
            > first.get_query().data() + first.get_query().size()) {
        assert(extensions[0].get_nodes().front() == next.get_nodes().front());
        extensions[0].extend_query_end(query.data() + query.size());
        if (left.empty()) {
            std::swap(first, extensions[0]);
        } else {
            left.trim_end_clipping();
            extensions[0].trim_clipping();
            left.append(std::move(extensions[0]));
            std::swap(left, first);
        }

        if (second.get_score() > first.get_score()) {
            std::swap(first, second);
        }

    } else if (!second.get_offset()
            && first.get_clipping() + first.get_query().size() > second.get_clipping()) {
        ssize_t overlap = first.get_clipping() + first.get_query().size()
                            - second.get_clipping();
        second.trim_query_prefix(overlap, graph, config, false);
        second.trim_clipping();
        second.insert_gap_prefix(-overlap, graph, config);
        first.trim_end_clipping();
        first.append(std::move(second));
        assert(first.get_nodes().front());
        assert(first.get_nodes().back());
    } else {
        std::swap(first, second);
    }

    return extender->num_explored_nodes();
}

template <class AlignmentCompare>
void ISeedAndExtendAligner<AlignmentCompare>
::align_batch(const std::vector<IDBGAligner::Query> &seq_batch,
              const AlignmentCallback &callback) const {
    for (const auto &[header, query, is_reverse_complement] : seq_batch) {
        size_t num_seeds = 0;
        size_t num_explored_nodes = 0;
        size_t num_extensions = 0;

        QueryAlignment paths(query, is_reverse_complement);
        Aggregator aggregator(graph_, paths.get_query(false), paths.get_query(true),
                              config_);

        auto add_alignment = [&](Alignment&& alignment) {
            assert(alignment.is_valid(graph_, &config_));
            aggregator.add_alignment(std::move(alignment));
        };

        auto get_min_path_score = [&](const Alignment &seed) {
            return std::max(config_.min_path_score, aggregator.get_min_path_score(seed));
        };

        std::string_view this_query = paths.get_query(is_reverse_complement);
        assert(this_query == query);

        std::vector<node_index> nodes;
        if (config_.max_seed_length >= graph_.get_k()) {
            nodes = map_sequence_to_nodes(graph_, query);
        } else if (this_query.size() >= graph_.get_k()) {
            nodes.resize(this_query.size() - graph_.get_k() + 1);
        }

        auto seeder = build_seeder(this_query, is_reverse_complement, nodes);
        auto extender = build_extender(this_query, aggregator, config_);

#if ! _PROTEIN_GRAPH
        if (graph_.get_mode() == DeBruijnGraph::CANONICAL
                || config_.forward_and_reverse_complement) {
            std::vector<node_index> nodes_rc(nodes);
            std::string dummy(query);
            if (config_.max_seed_length >= graph_.get_k()) {
                reverse_complement_seq_path(graph_, dummy, nodes_rc);
                assert(dummy == paths.get_query(!is_reverse_complement));
            }

            assert(nodes_rc.size() == nodes.size());

            std::string_view reverse = paths.get_query(!is_reverse_complement);

            auto seeder_rc = build_seeder(reverse, !is_reverse_complement, nodes_rc);
            auto extender_rc = build_extender(reverse, aggregator, config_);

            auto [seeds, extensions, explored_nodes] =
                align_both_directions(this_query, reverse, *seeder, *seeder_rc,
                                      *extender, *extender_rc,
                                      add_alignment, get_min_path_score);

            num_seeds += seeds;
            num_extensions += extensions;
            num_explored_nodes += explored_nodes;

        } else {
            align_core(*seeder, *extender, add_alignment, get_min_path_score, false);
        }
#else
        align_core(*seeder, *extender, add_alignment, get_min_path_score, false);
#endif

        num_explored_nodes += extender->num_explored_nodes();
        num_extensions += extender->num_extensions();
        size_t num_labels = std::max(size_t { 1 }, aggregator.num_labels() - 1);
        score_t best_score = std::numeric_limits<score_t>::min();
        size_t query_coverage = 0;



        for (auto&& alignment : chain_alignments<AlignmentCompare>(aggregator.get_alignments(),
                                                                   paths.get_query(false),
                                                                   paths.get_query(true),
                                                                   config_, graph_)) {
            assert(alignment.is_valid(graph_, &config_));
            if (alignment.get_score() > best_score) {
                best_score = alignment.get_score();
                query_coverage = alignment.get_query().size();
            }
            paths.emplace_back(std::forward<decltype(alignment)>(alignment));
        }

        common::logger->trace(
            "{}\tlength: {}\tcovered: {}\tbest score: {}\tseeds: {}\textensions: {}\t"
            "explored nodes: {}\texplored nodes/extension: {:.2f}\texplored nodes/k-mer: {:.2f}\t"
            "labels: {}\tnodes/k-mer/label: {:.2f}",
            header, query.size(), query_coverage, best_score, num_seeds,
            num_extensions, num_explored_nodes,
            static_cast<double>(num_explored_nodes) / num_extensions,
            static_cast<double>(num_explored_nodes) / nodes.size(),
            num_labels,
            static_cast<double>(num_explored_nodes) / nodes.size() / num_labels
        );

        callback(header, std::move(paths));
    };
}

template <class AlignmentCompare>
void ISeedAndExtendAligner<AlignmentCompare>
::align_core(const ISeeder &seeder,
             IExtender &extender,
             const std::function<void(Alignment&&)> &callback,
             const std::function<score_t(const Alignment&)> &get_min_path_score,
             bool force_fixed_seed) const {
    auto seeds = seeder.get_seeds();

    for (size_t i = 0; i < seeds.size(); ++i) {
        if (seeds[i].empty())
            continue;

        score_t min_path_score = get_min_path_score(seeds[i]);

        DEBUG_LOG("Min path score: {}\tSeed: {}", min_path_score, seeds[i]);

        for (auto&& extension : extender.get_extensions(seeds[i], min_path_score,
                                                        force_fixed_seed)) {
            DEBUG_LOG("\tExtension: {}", extension);
            callback(std::move(extension));
        }

        for (size_t j = i + 1; j < seeds.size(); ++j) {
            if (seeds[j].size() && !extender.check_seed(seeds[j]))
                filter_seed(seeds[i], seeds[j]);
        }
    }
}

// Construct a full alignment from a chain by aligning the query agaisnt
// the graph in the regions of the query in between the chain seeds.
template <class AlignmentCompare>
void ISeedAndExtendAligner<AlignmentCompare>
::extend_chain(std::string_view query,
               std::string_view query_rc,
               Chain&& chain,
               score_t score,
               size_t &num_extensions,
               size_t &num_explored_nodes,
               const std::function<void(Alignment&&)> &callback) const {
    DBGAlignerConfig gap_fill_config = config_;
    gap_fill_config.semiglobal = true;
    gap_fill_config.allow_left_trim = false;
    gap_fill_config.trim_offset_after_extend = false;

    std::string_view fw = chain[0].get_orientation() ? query_rc : query;
    std::string_view bw = chain[0].get_orientation() ? query : query_rc;
    AlignmentAggregator<AlignmentCompare> dummy(graph_, fw, bw, config_);
    const char *query_end = query.data() + query.size();
    assert(chain.size());

    std::ignore = score;
    DEBUG_LOG("Chain size: {}, score: {}", chain.size(), score);
#ifndef NDEBUG
    if (common::get_verbose()) {
        for (const Alignment &c : chain) {
            DEBUG_LOG("\t{}", c);
        }
    }
#endif

    Alignment cur = std::move(chain[0]);
    Alignment best = cur;

    for (size_t i = 1; i < chain.size(); ++i) {
        gap_fill_config.xdrop = config_.xdrop;
        gap_fill_config.terminal_node = chain[i].get_nodes().back();
        gap_fill_config.right_end_bonus = chain[i].get_query().data()
                                            + chain[i].get_query().size() == query_end
            ? config_.right_end_bonus : 0;
        num_explored_nodes += align_connect(query, graph_, gap_fill_config, cur, chain[i],
                                            [&](std::string_view query_window,
                                                const DBGAlignerConfig &config) {
            return build_extender(query_window, dummy, config);
        });
        if (cur.empty())
            return;
    }
    num_extensions += chain.size() - 1;


    if (AlignmentCompare()(best, cur))
        std::swap(best, cur);

    best.extend_query_begin(query.data());
    best.extend_query_end(query_end);

    if (best.get_query().data() + best.get_query().size() < query_end) {
        DBGAlignerConfig end_cfg = config_;
        end_cfg.trim_offset_after_extend = false;
        end_cfg.allow_left_trim = false;

        auto extender = build_extender(query, dummy, end_cfg);
        auto extensions = extender->get_extensions(best, config_.ninf, true);
        ++num_extensions;
        num_explored_nodes += extender->num_explored_nodes();
        if (extensions.size() && extensions[0].get_query().data()
                                    + extensions[0].get_query().size()
                > best.get_query().data() + best.get_query().size()
                && extensions[0].get_score() > best.get_score()) {
            std::swap(best, extensions[0]);
            assert(best.is_valid(graph_, &config_));
        }
    }

    best.extend_query_begin(query.data());
    best.extend_query_end(query_end);
    best.trim_offset();

    assert(best.is_valid(graph_, &config_));

    if (best.get_clipping()) {
        RCDBG rc_dbg(graph_);
        const DeBruijnGraph &rc_graph = graph_.get_mode() != DeBruijnGraph::CANONICAL
            ? rc_dbg : graph_;

        auto rev = best;
        rev.reverse_complement(rc_graph, query_rc);
        if (rev.size() && rev.get_nodes().back()) {
            assert(rev.get_end_clipping());
            DBGAlignerConfig no_trim_config = config_;
            no_trim_config.allow_left_trim = false;
            auto extender = build_extender(query_rc, dummy, no_trim_config);
            extender->set_graph(rc_graph);
            auto extensions = extender->get_extensions(rev, config_.ninf, true);
            if (extensions.size() && extensions[0].get_query().data()
                                        + extensions[0].get_query().size()
                    > rev.get_query().data() + rev.get_query().size()
                    && extensions[0].get_query().data() == rev.get_query().data()
                    && extensions[0].get_nodes()[0] == rev.get_nodes()[0]) {
                extensions[0].extend_query_begin(query_rc.data());
                extensions[0].extend_query_end(query_rc.data() + query_rc.size());
                extensions[0].reverse_complement(rc_graph, query);
                if (extensions[0].size()) {
                    std::swap(best, extensions[0]);
                    best.extend_query_begin(query.data());
                    best.extend_query_end(query_end);
                }
            }

            ++num_extensions;
            num_explored_nodes += extender->num_explored_nodes() - rev.size();
        }
    }

    callback(std::move(best));
}

template <class AlignmentCompare>
std::tuple<size_t, size_t, size_t> ISeedAndExtendAligner<AlignmentCompare>
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

    size_t num_seeds = 0;
    size_t num_extensions = 0;
    size_t num_explored_nodes = 0;

    auto fwd_seeds = forward_seeder.get_seeds();
    auto bwd_seeds = reverse_seeder.get_seeds();

    if (config_.chain_alignments) {
        bool can_chain = false;
        for (const Alignment &seed : fwd_seeds) {
            if (seed.label_coordinates.size()) {
                can_chain = true;
                break;
            }
        }

        if (!can_chain) {
            for (const Alignment &seed : bwd_seeds) {
                if (seed.label_coordinates.size()) {
                    can_chain = true;
                    break;
                }
            }
        }

        if (!can_chain) {
            common::logger->error("Chaining only supported for seeds with coordinates. Skipping seed chaining.");
            exit(1);
        }

        Aggregator aggregator(graph_, forward, reverse, config_);

        bool terminate = false;
        size_t this_num_explored;
        std::tie(num_seeds, this_num_explored) = call_seed_chains_both_strands(
            forward, reverse, graph_, config_, std::move(fwd_seeds), std::move(bwd_seeds),
            [&](Chain&& chain, score_t score) {
                extend_chain(chain[0].get_orientation() ? reverse : forward,
                             chain[0].get_orientation() ? forward : reverse,
                             std::move(chain), score, num_extensions, num_explored_nodes,
                             [&](Alignment&& aln) {
                                 terminate |= aggregator.add_alignment(std::move(aln));
                             });
            },
            [&]() { return terminate; }
        );
        num_explored_nodes += this_num_explored;

        for (Alignment &alignment : aggregator.get_alignments()) {
            if (alignment.get_score() < get_min_path_score(alignment))
                continue;

            if (graph_.get_mode() == DeBruijnGraph::CANONICAL && alignment.get_orientation()) {
                Alignment rev(alignment);
                rev.reverse_complement(graph_, forward);
                if (rev.size())
                    std::swap(rev, alignment);
            }

            callback(std::move(alignment));
        }

        return std::make_tuple(num_seeds, num_extensions, num_explored_nodes);
    }

    RCDBG rc_dbg(graph_);
    bool use_rcdbg = graph_.get_mode() != DeBruijnGraph::CANONICAL
                        && config_.forward_and_reverse_complement;

    const DeBruijnGraph &rc_graph = use_rcdbg ? rc_dbg : graph_;

    auto is_reversible = [this](const Alignment &alignment) {
        return graph_.get_mode() == DeBruijnGraph::CANONICAL
            && alignment.get_orientation()
            && !alignment.get_offset();
    };

    auto aln_both = [&](std::string_view query,
                        std::string_view query_rc,
                        std::vector<Alignment>&& seeds,
                        IExtender &fwd_extender,
                        IExtender &bwd_extender,
                        const std::function<void(Alignment&&)> &callback) {
        fwd_extender.set_graph(graph_);
        bwd_extender.set_graph(rc_graph);
        num_seeds = seeds.size();

        if (seeds.empty())
            return;

        for (size_t i = 0; i < seeds.size(); ++i) {
            if (seeds[i].empty())
                continue;

            score_t min_path_score = config_.min_cell_score;

            DEBUG_LOG("Min path score: {}\tSeed: {}", min_path_score, seeds[i]);

            auto extensions = fwd_extender.get_extensions(seeds[i], min_path_score, false);

            std::vector<Alignment> rc_of_alignments;

            for (Alignment &path : extensions) {
                if (path.get_score() >= get_min_path_score(path)) {
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
                    continue;

                path.reverse_complement(rc_graph, query_rc);

                if (path.empty()) {
                    DEBUG_LOG("This local alignment cannot be reversed, skipping");
                    continue;
                }

                // Remove any character skipping from the end so that the
                // alignment can proceed
                assert(path.get_end_clipping());
                assert(path.is_valid(rc_graph, &config_));

                rc_of_alignments.emplace_back(std::move(path));
            }

            align_core(ManualSeeder(std::move(rc_of_alignments)), bwd_extender,
                [&](Alignment&& path) {
                    if (use_rcdbg || is_reversible(path)) {
                        path.reverse_complement(rc_graph, query);
                        if (path.empty())
                            return;

                        if (auto *filter = dynamic_cast<SeedFilteringExtender*>(&fwd_extender)) {
                            for (node_index node : path.get_nodes()) {
                                filter->filter_nodes(
                                    node, path.get_clipping(),
                                    query.size() - path.get_end_clipping()
                                );
                            }
                        }
                    }

                    assert(path.is_valid(graph_, &config_));

                    callback(std::move(path));
                },
                get_min_path_score,
                true /* alignments must have the seed as a prefix */
            );

            for (size_t j = i + 1; j < seeds.size(); ++j) {
                if (seeds[j].size() && !fwd_extender.check_seed(seeds[j]))
                    filter_seed(seeds[i], seeds[j]);
            }
        }
    };

    size_t fwd_num_matches = forward_seeder.get_num_matches();
    size_t bwd_num_matches = reverse_seeder.get_num_matches();

    if (fwd_num_matches >= bwd_num_matches) {
        aln_both(forward, reverse, std::move(fwd_seeds),
                 forward_extender, reverse_extender, callback);
        if (bwd_num_matches >= fwd_num_matches * config_.rel_score_cutoff) {
            aln_both(reverse, forward, std::move(bwd_seeds),
                     reverse_extender, forward_extender, callback);
        }
    } else {
        aln_both(reverse, forward, std::move(bwd_seeds),
                 reverse_extender, forward_extender, callback);
        if (fwd_num_matches >= bwd_num_matches * config_.rel_score_cutoff) {
            aln_both(forward, reverse, std::move(fwd_seeds),
                     forward_extender, reverse_extender, callback);
        }
    }

    return std::make_tuple(num_seeds, num_extensions, num_explored_nodes);
}

template class ISeedAndExtendAligner<>;

} // namespace align
} // namespace graph
} // namespace mtg
