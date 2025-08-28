#include "aln_seeder.hpp"

#include <queue>

#include <tsl/hopscotch_set.h>

#include "aln_chainer.hpp"
#include "common/logger.hpp"
#include "common/vector_map.hpp"
#include "common/algorithms.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/representation/canonical_dbg.hpp"

namespace mtg::graph::align_redone {

void align_query(const Query &query,
                 const Seeder &seeder,
                 const Extender &extender,
                 const std::function<void(Alignment&&)> &callback,
                 bool connect_anchors_in_chain) {
    align_query(query, seeder,
                [&](Alignment&& base_path) { extender.extend(base_path, callback); });
}

void align_query(const Query &query,
                 const Seeder &seeder,
                 const std::function<void(Alignment&&)> &callback,
                 bool connect_anchors_in_chain) {
    for (Alignment &path : seeder.get_inexact_anchors(connect_anchors_in_chain)) {
        callback(std::move(path));
    }
}

using TraverseCallback = std::function<void(DeBruijnGraph::node_index, Anchor::label_class_t)>;
using EdgeCallback = std::function<void(DeBruijnGraph::node_index, char, Anchor::label_class_t)>;
using Terminator = std::function<bool()>;

class AlignmentGraph {
  public:
    using node_index = DeBruijnGraph::node_index;
    using label_class_t = Anchor::label_class_t;

    AlignmentGraph(const DeBruijnGraph &graph,
                   AnnotationBuffer *anno_buffer = nullptr,
                   label_class_t target = Anchor::nannot)
          : graph_(graph), anno_buffer_(anno_buffer), target_(target) {
        assert((target == Anchor::nannot) == !anno_buffer_);
    }

    AnnotationBuffer* get_anno_buffer() const { return anno_buffer_; }

    std::pair<node_index, label_class_t> traverse(node_index node, char c) const {
        assert(c != boss::BOSS::kSentinel);
        node_index next = graph_.traverse(node, c);
        if (target_ == Anchor::nannot)
            return std::make_pair(next, target_);

        if (!next)
            return std::make_pair(next, 0);

        label_class_t next_labels = get_intersect_labels(next, target_);
        assert(next_labels != Anchor::nannot);
        if (!next_labels)
            next = DeBruijnGraph::npos;

        return std::make_pair(next, next_labels);
    }

    std::pair<node_index, label_class_t> traverse_back(node_index node, char c) const {
        assert(c != boss::BOSS::kSentinel);
        node_index prev = graph_.traverse_back(node, c);
        if (target_ == Anchor::nannot)
            return std::make_pair(prev, target_);

        if (!prev)
            return std::make_pair(prev, 0);

        label_class_t prev_labels = get_intersect_labels(prev, target_);
        assert(prev_labels != Anchor::nannot);
        if (!prev_labels)
            prev = DeBruijnGraph::npos;

        return std::make_pair(prev, prev_labels);
    }

    void traverse(node_index node, std::string_view seq,
                  const TraverseCallback &callback,
                  const Terminator &terminate = [](){ return false; }) const {
        if (target_ == Anchor::nannot) {
            graph_.traverse(node, seq,
                            [&](node_index next) {
                                assert(has_labels(next, target_));
                                callback(next, target_);
                            }, terminate);
            return;
        }

        std::vector<node_index> forward_nodes;
        graph_.traverse(node, seq,
                        [&](node_index next) { forward_nodes.emplace_back(next); },
                        terminate);
        anno_buffer_->queue_path(forward_nodes);
        anno_buffer_->fetch_queued_annotations();

        label_class_t cur_target = target_;
        for (node_index next : forward_nodes) {
            cur_target = get_intersect_labels(next, cur_target);
            assert(cur_target != Anchor::nannot);
            if (cur_target) {
                assert(has_labels(next, cur_target));
                callback(next, cur_target);
            } else {
                return;
            }
        }
    }

    void traverse_back(node_index node, std::string_view prefix_seq,
                       const TraverseCallback &callback,
                       const Terminator &terminate = [](){ return false; }) const {
        if (target_ == Anchor::nannot) {
            graph_.traverse_back(node, prefix_seq,
                                 [&](node_index prev) {
                                     assert(has_labels(prev, target_));
                                     callback(prev, target_);
                                 },
                                 terminate);
            return;
        }

        std::vector<node_index> backward_nodes;
        graph_.traverse_back(node, prefix_seq,
                             [&](node_index prev) { backward_nodes.emplace_back(prev); },
                             terminate);
        anno_buffer_->queue_path(std::vector<node_index>(backward_nodes.rbegin(),
                                                         backward_nodes.rend()));
        anno_buffer_->fetch_queued_annotations();

        label_class_t cur_target = target_;
        for (node_index prev : backward_nodes) {
            cur_target = get_intersect_labels(prev, cur_target);
            assert(cur_target != Anchor::nannot);
            if (cur_target) {
                assert(has_labels(prev, cur_target));
                callback(prev, cur_target);
            } else {
                return;
            }
        }
    }

    void call_outgoing_kmers(node_index node, const EdgeCallback &callback) const {
        if (target_ == Anchor::nannot) {
            graph_.call_outgoing_kmers(node, [&](node_index next, char c) { callback(next, c, target_); });
        } else {
            std::vector<node_index> forward_n;
            std::vector<char> forward_c;
            graph_.call_outgoing_kmers(node, [&](node_index next, char c) {
                if (c != boss::BOSS::kSentinel) {
                    forward_n.emplace_back(next);
                    forward_c.emplace_back(c);
                    anno_buffer_->queue_path(std::vector<node_index>{ next });
                }
            });
            anno_buffer_->fetch_queued_annotations();

            for (size_t i = 0; i < forward_n.size(); ++i) {
                node_index next = forward_n[i];
                char c = forward_c[i];
                label_class_t next_target = get_intersect_labels(next, target_);
                assert(next_target != Anchor::nannot);
                if (next_target) {
                    assert(has_labels(next, next_target));
                    callback(next, c, next_target);
                }
            }
        }
    }

    void call_incoming_kmers(node_index node, const EdgeCallback &callback) const {
        // std::cerr << "cik\t" << node << "\t" << target_ << std::endl;
        if (target_ == Anchor::nannot) {
            graph_.call_incoming_kmers(node, [&](node_index prev, char c) { callback(prev, c, target_); });
        } else {
            std::vector<node_index> backward_n;
            std::vector<char> backward_c;
            graph_.call_incoming_kmers(node, [&](node_index prev, char c) {
                if (c != boss::BOSS::kSentinel) {
                    backward_n.emplace_back(prev);
                    backward_c.emplace_back(c);
                    anno_buffer_->queue_path(std::vector<node_index>{ prev });
                }
            });
            anno_buffer_->fetch_queued_annotations();

            for (size_t i = 0; i < backward_n.size(); ++i) {
                node_index prev = backward_n[i];
                char c = backward_c[i];
                label_class_t prev_target = get_intersect_labels(prev, target_);
                assert(prev_target != Anchor::nannot);
                if (prev_target) {
                    assert(has_labels(prev, prev_target));
                    callback(prev, c, prev_target);
                }
            }
        }
    }

    bool has_single_outgoing(node_index node) const {
        if (target_ == Anchor::nannot)
            return graph_.has_single_outgoing(node);

        size_t outdegree = 0;
        adjacent_outgoing_nodes(node, [&](node_index, label_class_t) { ++outdegree; });
        return outdegree == 1;
    }

    bool has_single_incoming(node_index node) const {
        // std::cerr << "hsi\t" << node << "\t" << target_ << std::endl;
        if (target_ == Anchor::nannot) {
            return graph_.has_single_incoming(node);
        } else {
            size_t indegree = 0;
            adjacent_incoming_nodes(node, [&](node_index, label_class_t) { ++indegree; });
            return indegree == 1;
        }
    }

    void adjacent_outgoing_nodes(node_index node, const TraverseCallback &callback) const {
        if (target_ == Anchor::nannot) {
            graph_.adjacent_outgoing_nodes(node, [&](node_index next) { callback(next, target_); });
        } else {
            // TODO: replace this once we can guarantee that we won't get
            //       dummy nodes
            call_outgoing_kmers(node, [&](node_index next, char, label_class_t next_target) {
                callback(next, next_target);
            });
            // std::vector<node_index> forward_n;
            // graph_.adjacent_outgoing_nodes(node, [&](node_index next) {
            //     forward_n.emplace_back(next);
            //     anno_buffer_->queue_path(std::vector<node_index>{ next });
            // });
            // anno_buffer_->fetch_queued_annotations();

            // for (node_index next : forward_n) {
            //     if (has_labels(next))
            //         callback(next);
            // }
        }
    }

    void adjacent_incoming_nodes(node_index node, const TraverseCallback &callback) const {
        if (target_ == Anchor::nannot) {
            graph_.adjacent_incoming_nodes(node, [&](node_index prev) { callback(prev, target_); });
        } else {
            // TODO: replace this once we can guarantee that we won't get
            //       dummy nodes
            call_incoming_kmers(node, [&](node_index prev, char, label_class_t prev_target) {
                callback(prev, prev_target);
            });
            // std::vector<node_index> backward_n;
            // graph_.adjacent_incoming_nodes(node, [&](node_index prev) {
            //     backward_n.emplace_back(prev);
            //     anno_buffer_->queue_path(std::vector<node_index>{ prev });
            // });
            // anno_buffer_->fetch_queued_annotations();

            // for (node_index prev : backward_n) {
            //     if (has_labels(prev))
            //         callback(prev);
            // }
        }
    }

    bool has_labels(node_index node, label_class_t target_labels_id) const {
        if (!anno_buffer_) {
            assert(target_labels_id == Anchor::nannot);
            return true;
        }

        assert(target_labels_id != Anchor::nannot);
        label_class_t node_labels_id = get_label_class(node);
        if (node_labels_id == Anchor::nannot)
            return false;

        if (node_labels_id == target_labels_id)
            return true;

        const auto &node_labels = anno_buffer_->get_cached_column_set(node_labels_id);
        const auto &target_labels = anno_buffer_->get_cached_column_set(target_labels_id);

        return utils::count_intersection(node_labels.begin(), node_labels.end(),
                                         target_labels.begin(), target_labels.end()) == target_labels.size();
    }

  private:
    const DeBruijnGraph &graph_;
    AnnotationBuffer* const anno_buffer_;
    label_class_t target_;

    label_class_t get_label_class(node_index node) const {
        assert(anno_buffer_);

        label_class_t node_labels_id = anno_buffer_->get_labels_id(node);
        if (node_labels_id == Anchor::nannot) {
            anno_buffer_->queue_path(std::vector<node_index>{ node });
            anno_buffer_->fetch_queued_annotations();
            node_labels_id = anno_buffer_->get_labels_id(node);
        }
        assert(node_labels_id != Anchor::nannot);
        return node_labels_id;
    }

    label_class_t get_intersect_labels(node_index node, label_class_t target_labels_id) const {
        assert(anno_buffer_);
        assert(target_labels_id != Anchor::nannot);

        label_class_t node_labels_id = get_label_class(node);
        assert(node_labels_id != Anchor::nannot);

        if (node_labels_id == target_labels_id)
            return node_labels_id;

        const auto &target_labels = anno_buffer_->get_cached_column_set(target_labels_id);
        assert(std::find(target_labels.begin(), target_labels.end(), AnnotationBuffer::ncolumn) == target_labels.end());
        assert(std::all_of(target_labels.begin(), target_labels.end(),
                           [&](auto c) { return c < anno_buffer_->get_annotator().get_label_encoder().size(); }));

        const auto &node_labels = anno_buffer_->get_cached_column_set(node_labels_id);
        assert(std::find(node_labels.begin(), node_labels.end(), AnnotationBuffer::ncolumn) == node_labels.end());
        assert(std::all_of(node_labels.begin(), node_labels.end(),
                           [&](auto c) { return c < anno_buffer_->get_annotator().get_label_encoder().size(); }));

        AnnotationBuffer::Columns columns;
        std::set_intersection(node_labels.begin(), node_labels.end(),
                              target_labels.begin(), target_labels.end(),
                              std::back_inserter(columns));

        return columns.size() ? anno_buffer_->cache_column_set(std::move(columns)) : 0;
    }
};

template <typename T>
struct template_parameter;

template <template <typename...> class C, typename T>
struct template_parameter<C<T>> {
    using type = T;
};

template <typename T>
using get_kmer_t = typename template_parameter<std::decay_t<T>>::type;

std::vector<Anchor> ExactSeeder::get_anchors() const {
    const DeBruijnGraph &graph = query_.get_graph();
    std::vector<Anchor> anchors;

    if (query_.get_query().size() < config_.min_seed_length)
        return anchors;

    std::vector<std::vector<DeBruijnGraph::node_index>> mapped_nodes(2);
    std::vector<std::vector<tsl::hopscotch_map<DeBruijnGraph::node_index, Anchor>>> max_seeds(2);

    for (bool orientation : { false, true }) {
        std::string_view this_query = query_.get_query(orientation);
        auto &nodes = mapped_nodes[orientation];
        nodes = map_to_nodes_sequentially(graph, this_query);

        auto &cur_max_seeds = max_seeds[orientation];
        cur_max_seeds.resize(this_query.size() - config_.min_seed_length + 1);

        for (size_t begin = 0; begin < nodes.size(); ++begin) {
            Match::node_index node = nodes[begin];
            if (node != DeBruijnGraph::npos) {
                cur_max_seeds[begin][node] = Anchor(
                    this_query,
                    begin, begin + graph.get_k(),
                    orientation,
                    std::vector<Match::node_index>{ node }
                );
            }
        }
    }

    for (bool orientation : { false, true }) {
        std::string_view this_query = query_.get_query(orientation);
        auto &nodes = mapped_nodes[orientation];
        auto &cur_max_seeds = max_seeds[orientation];

        if (config_.min_seed_length < graph.get_k()) {
            const DBGSuccinct *dbg_succ = nullptr;
            const DBGSSHash *sshash = nullptr;
            const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph);
            if (canonical) {
                dbg_succ = dynamic_cast<const DBGSuccinct*>(&canonical->get_graph());
                assert(!dynamic_cast<const DBGSSHash*>(&canonical->get_graph()));
            } else {
                dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph);
                sshash = dynamic_cast<const DBGSSHash*>(&graph);
            }

            if (dbg_succ) {
                const auto &boss = dbg_succ->get_boss();
                auto encoded = boss.encode(this_query);

                size_t k = graph.get_k();
                for (size_t i = 0; i + config_.min_seed_length <= this_query.size(); ++i) {
                    if (std::find(encoded.begin() + i, encoded.begin() + i + config_.min_seed_length,
                                boss.alph_size) != encoded.begin() + i + config_.min_seed_length) {
                        continue;
                    }

                    using edge_index = boss::BOSS::edge_index;
                    edge_index first;
                    edge_index last;
                    auto it = encoded.begin();
                    std::tie(first, last, it) = boss.index_range(
                        encoded.begin() + i,
                        encoded.begin() + i + config_.min_seed_length
                    );

                    size_t match_size = it - (encoded.begin() + i);
                    if (match_size < config_.min_seed_length)
                        continue;

                    first = first == boss.select_last(1) ? 1 : boss.pred_last(first - 1) + 1;
                    // std::cerr << i << "\t" << std::string_view(this_query.begin() + i, match_size) << "\t" << first << "," << last << "\t" << std::flush << graph.get_node_sequence(first) << "\t" << graph.get_node_sequence(last) << "\n";

                    if (canonical) {
                        if (dbg_succ->in_graph(last)) {
                            size_t i_rc = this_query.size() - (i + match_size);
                            dbg_succ->adjacent_incoming_nodes(last, [&](DeBruijnGraph::node_index node) {
                                node = canonical->reverse_complement(node);
                                bool inserted = true;
                                auto it = max_seeds[!orientation][i_rc].find(node);
                                if (it != max_seeds[!orientation][i_rc].end() && match_size <= it->second.get_seed().size())
                                    return;

                                std::string spelling = graph.get_node_sequence(node);
                                if (spelling.find(boss::BOSS::kSentinel) != std::string::npos)
                                    return;

                                assert(std::equal(query_.get_query(!orientation).begin() + i_rc,
                                                  query_.get_query(!orientation).begin() + i_rc + match_size,
                                                  spelling.begin(), spelling.begin() + match_size));

                                std::tie(it, inserted) = max_seeds[!orientation][i_rc].try_emplace(node, Anchor());
                                it.value() = Anchor(query_.get_query(!orientation),
                                                    i_rc, i_rc + match_size,
                                                    !orientation,
                                                    std::vector<Match::node_index>{ node },
                                                    spelling.substr(match_size));
                                assert(it->second.get_path_spelling().find(boss::BOSS::kSentinel) == std::string::npos);
                                assert(it->second.get_path_spelling().size() == graph.get_k());
                                assert(it->second.get_spelling().size() + it->second.get_end_trim() == graph.get_k());
                                assert(it->second.is_spelling_valid(graph));
                            });
                        }
                    }

                    std::vector<std::tuple<edge_index, edge_index, std::string, bool>> traverse;
                    std::string suffix;
                    bool is_exact_match = true;
                    traverse.emplace_back(first, last, suffix, is_exact_match);
                    while (traverse.size()) {
                        std::tie(first, last, suffix, is_exact_match) = traverse.back();
                        assert(boss.get_last(last));
                        traverse.pop_back();

                        size_t cur_match_size = match_size + suffix.size();

                        assert(cur_match_size < k);
                        if (cur_match_size == k - 1) {
                            assert(first == last || !boss.get_last(first));
                            assert(boss.succ_last(first) == last);
                            // std::cerr << "foo\t" << i << "\t" << std::string_view(this_query.begin() + i, match_size) << "\n";
                            for (edge_index node = first; node <= last; ++node) {
                                if (i < nodes.size() && node == nodes[i])
                                    continue;

                                if (dbg_succ->in_graph(node)) {
                                    char c = boss.decode(boss.get_W(node) % boss.alph_size);
                                    if (c != boss::BOSS::kSentinel) {
                                        // std::cerr << "\t" << suffix << c << "\n";
                                        auto first_mm = std::mismatch(this_query.begin() + i + match_size, this_query.end(),
                                                                    suffix.begin(), suffix.end()).second;
                                        bool last_char_matches = i + k <= this_query.size() && first_mm == suffix.end() && c == this_query[i + k - 1];
                                        size_t full_match_size = match_size + (first_mm - suffix.begin()) + last_char_matches;

                                        is_exact_match &= last_char_matches;

                                        if (!is_exact_match) {
                                            auto [it, inserted] = cur_max_seeds[i].try_emplace(node, Anchor());
                                            if (full_match_size > it->second.get_seed().size()) {
                                                it.value() = Anchor(this_query,
                                                                    i, i + full_match_size,
                                                                    orientation,
                                                                    std::vector<Match::node_index>{ node },
                                                                    (suffix + c).substr(full_match_size - match_size));
                                                assert(it->second.get_path_spelling().find(boss::BOSS::kSentinel) == std::string::npos);
                                                assert(it->second.get_path_spelling().size() == graph.get_k());
                                                assert(it->second.get_spelling().size() + it->second.get_end_trim() == graph.get_k());
                                                assert(it->second.is_spelling_valid(graph));
                                            }
                                        }
                                    }
                                }
                            }

                            continue;
                        }

                        boss.call_tighter_ranges(first, last,
                            [&](edge_index next_first, edge_index next_last, boss::BOSS::TAlphabet s) {
                                char c = boss.decode(s);
                                bool next_is_exact_match = is_exact_match && i + cur_match_size < this_query.size() && this_query[i + cur_match_size] == c;
                                if (!canonical || next_is_exact_match || nodes[i] == DeBruijnGraph::npos) {
                                    traverse.emplace_back(next_first, next_last, suffix + c,
                                                        next_is_exact_match);
                                }
                                if (canonical && next_is_exact_match) {
                                    size_t next_match_size = cur_match_size + 1;
                                    size_t i_rc = this_query.size() - (i + next_match_size);

                                    std::vector<edge_index> succs;
                                    for (edge_index e = boss.succ_W(first, s); e <= last; e = e + 1 < boss.num_edges() ? boss.succ_W(e + 1, s) : boss.num_edges() + 1) {
                                        succs.emplace_back(e);
                                    }
                                    boss::BOSS::TAlphabet splus = s + boss.alph_size;
                                    for (edge_index e = boss.succ_W(first, splus); e <= last; e = e + 1 < boss.num_edges() ? boss.succ_W(e + 1, splus) : boss.num_edges() + 1) {
                                        succs.emplace_back(e);
                                    }

                                    for (edge_index e : succs) {
                                        auto node = canonical->reverse_complement(e);
                                        bool inserted = true;
                                        auto it = max_seeds[!orientation][i_rc].find(node);
                                        if (it != max_seeds[!orientation][i_rc].end() && next_match_size <= it->second.get_seed().size())
                                            continue;

                                        std::string spelling = graph.get_node_sequence(node);
                                        if (spelling.find(boss::BOSS::kSentinel) != std::string::npos)
                                            continue;

                                        assert(std::equal(query_.get_query(!orientation).begin() + i_rc,
                                                          query_.get_query(!orientation).begin() + i_rc + next_match_size,
                                                          spelling.begin(), spelling.begin() + next_match_size));

                                        std::tie(it, inserted) = max_seeds[!orientation][i_rc].try_emplace(node, Anchor());
                                        it.value() = Anchor(query_.get_query(!orientation),
                                                            i_rc, i_rc + next_match_size,
                                                            !orientation,
                                                            std::vector<Match::node_index>{ node },
                                                            spelling.substr(next_match_size));
                                        assert(it->second.get_path_spelling().find(boss::BOSS::kSentinel) == std::string::npos);
                                        assert(it->second.get_path_spelling().size() == graph.get_k());
                                        assert(it->second.get_spelling().size() + it->second.get_end_trim() == graph.get_k());
                                        assert(it->second.is_spelling_valid(graph));
                                    }
                                }
                            }
                        );
                    }
                }
            } else if (sshash) {
                assert(!canonical);
                auto invalid_char = sshash->get_invalid_char_mask(this_query);

                std::visit([&](const auto &dict) {
                    using kmer_t = get_kmer_t<decltype(dict)>;
                    const auto &buckets = dict.data();
                    const auto &minimizers = dict.get_minimizers();
                    const auto &offsets = buckets.offsets;
                    size_t m = dict.m();
                    size_t k = dict.k();
                    size_t min_match_length = std::max(m, config_.min_seed_length);
                    auto invalid_mmer = utils::drag_and_mark_segments(invalid_char, true, m);
                    for (size_t i = 0; i + m <= this_query.size(); ++i) {
                        assert(i < invalid_mmer.size());
                        if (invalid_mmer[i])
                            continue;

                        kmer_t mmer = sshash::util::string_to_uint_kmer<kmer_t>(
                            this_query.data() + i,
                            m
                        );

                        uint64_t bucket_id = minimizers.lookup(uint64_t(mmer));
                        auto [begin, end] = buckets.locate_bucket(bucket_id);

                        for (size_t super_kmer_id = begin; super_kmer_id < end; ++super_kmer_id) {
                            uint64_t offset = offsets.access(super_kmer_id);

                            auto [res, contig_end] = buckets.offset_to_id(offset, k);
                            assert(res.contig_size);
                            assert(res.kmer_id_in_contig < res.contig_size);
                            assert(contig_end > offset);
                            auto contig_begin = offset - res.kmer_id_in_contig;
                            assert(contig_begin < contig_end);
                            assert(contig_begin <= offset);

                            size_t contig_length = contig_end - contig_begin;
                            std::ignore = contig_length;
                            assert(contig_length == res.contig_size + k - 1);

                            sshash::bit_vector_iterator<kmer_t> start_bv_it(dict.strings(), kmer_t::bits_per_char * offset);
                            assert(offset + m <= contig_end);
                            uint64_t window_size = std::min<uint64_t>({
                                k,
                                contig_end - offset - m
                            });
                            assert(offset + window_size + m <= contig_end);
                            // std::cerr << "foo\t" << i << "\t" << super_kmer_id << "\t" << sshash::util::uint_kmer_to_string<kmer_t>(mmer, m) << "\t"
                            //                      << contig_begin << "," << offset << "," << contig_end << "\t" << window_size << "\n";
                            for (uint64_t j = 0; j < window_size; ++j) {
                                kmer_t read_kmer = start_bv_it.read_and_advance_by_char(kmer_t::bits_per_char * m);
                                // std::cerr << "\tbar\t" << i + j << "\t" << sshash::util::uint_kmer_to_string<kmer_t>(read_kmer, m) << "\n";
                                if (read_kmer != mmer)
                                    continue;

                                // std::cerr << "\t\tgo!\n";

                                auto bv_it = start_bv_it;
                                for (size_t w = j; w < window_size; ++w) {
                                    size_t i_shift = i + w - j;
                                    if (i_shift + m > this_query.size() || invalid_mmer[i_shift])
                                        break;

                                    kmer_t mmer = sshash::util::string_to_uint_kmer<kmer_t>(
                                        this_query.data() + i_shift,
                                        m
                                    );

                                    // std::cerr << "\t\t" << i_shift << "\t" << sshash::util::uint_kmer_to_string<kmer_t>(mmer, m) << " vs. "
                                    //           << sshash::util::uint_kmer_to_string<kmer_t>(read_kmer, m) << "\n";
                                    if (read_kmer != mmer)
                                        break;

                                    // std::cerr << "\t\t\tnext\t" << contig_begin << "," << offset + w << "," << contig_end << "\n";
                                    // assert((res.kmer_id_in_contig + w + k <= contig_length)
                                    //         == (res.kmer_id_in_contig + w < res.contig_size));

                                    size_t i_in_contig = res.kmer_id_in_contig + w;

                                    if (i_in_contig < res.contig_size && i_shift + min_match_length <= this_query.size()) {
                                        // std::cerr << "\t\t\tfw\n";
                                        // forward
                                        assert(res.kmer_id_in_contig + w + k <= contig_length);
                                        DeBruijnGraph::node_index node = DBGSSHash::sshash_to_graph_index(res.kmer_id + w);
                                        // std::cerr << "\t" << graph.get_node_sequence(node) << "\n";
                                        if (i_shift >= nodes.size() || node != nodes[i_shift]) {
                                            std::string spelling = graph.get_node_sequence(node);
                                            assert(std::equal(this_query.begin() + i_shift, this_query.begin() + i_shift + m,
                                                              spelling.begin(), spelling.begin() + m));
                                            size_t this_end = std::min(i_shift + k, this_query.size());
                                            auto jt = std::mismatch(spelling.begin() + m, spelling.end(),
                                                                    this_query.begin() + i_shift + m, this_query.begin() + this_end).first;
                                            size_t full_match_size = jt - spelling.begin();
                                            assert(full_match_size <= k);
                                            if (full_match_size >= config_.min_seed_length) {
                                                auto [it, inserted] = cur_max_seeds[i_shift].try_emplace(node, Anchor());
                                                if (full_match_size > it->second.get_seed().size()) {
                                                    it.value() = Anchor(this_query,
                                                                        i_shift, i_shift + full_match_size,
                                                                        orientation,
                                                                        std::vector<Match::node_index>{ node },
                                                                        spelling.substr(full_match_size));
                                                    assert(it->second.get_path_spelling().size() == graph.get_k());
                                                    assert(it->second.get_spelling().size() + it->second.get_end_trim() == graph.get_k());
                                                    assert(it->second.is_spelling_valid(graph));
                                                }
                                            }
                                        }
                                    }

                                    if (sshash->get_mode() != DeBruijnGraph::BASIC && res.kmer_id_in_contig + w + m >= k && i_shift + m >= k) {
                                        // rev comp
                                        // std::cerr << "\t\t\trc\n";
                                        assert(sshash->get_mode() == DeBruijnGraph::CANONICAL);
                                        assert(res.kmer_id + w + m >= k);
                                        DeBruijnGraph::node_index node = DBGSSHash::sshash_to_graph_index(res.kmer_id + w + m - k);
                                        std::string spelling = graph.get_node_sequence(node);
                                        size_t i_back_shift = i_shift + m - k;
                                        assert(std::equal(this_query.begin() + i_shift, this_query.begin() + i_shift + m,
                                                          spelling.end() - m, spelling.end()));
                                        auto jt = std::mismatch(spelling.rbegin() + m, spelling.rend(),
                                                                std::make_reverse_iterator(this_query.begin() + i_shift),
                                                                std::make_reverse_iterator(this_query.begin() + i_back_shift)).first;
                                        size_t full_match_size = jt - spelling.rbegin();
                                        assert(full_match_size <= k);
                                        assert(full_match_size >= m);
                                        i_back_shift += k - full_match_size;
                                        if (full_match_size >= config_.min_seed_length) {
                                            assert(this_query.size() + k >= i_shift + m + full_match_size);
                                            size_t i_rc = this_query.size() - (i_back_shift + full_match_size);
                                            node = canonical ? canonical->reverse_complement(node)
                                                             : sshash->reverse_complement(node);
                                            auto [it, inserted] = max_seeds[!orientation][i_rc].try_emplace(node, Anchor());
                                            if (full_match_size > it->second.get_seed().size()) {
                                                ::reverse_complement(spelling.begin(), spelling.end());
                                                assert(graph.get_node_sequence(node) == spelling);
                                                assert(std::equal(query_.get_query(!orientation).begin() + i_rc,
                                                                  query_.get_query(!orientation).begin() + i_rc + full_match_size,
                                                                  spelling.begin(), spelling.begin() + full_match_size));
                                                it.value() = Anchor(query_.get_query(!orientation),
                                                                    i_rc, i_rc + full_match_size,
                                                                    !orientation,
                                                                    std::vector<Match::node_index>{ node },
                                                                    spelling.substr(full_match_size));
                                                assert(it->second.get_path_spelling().size() == graph.get_k());
                                                assert(it->second.get_spelling().size() + it->second.get_end_trim() == graph.get_k());
                                                assert(it->second.is_spelling_valid(graph));
                                            }
                                        }
                                    }

                                    if (w + 1 < window_size)
                                        read_kmer = bv_it.read_and_advance_by_char(kmer_t::bits_per_char * m);
                                }
                            }
                        }
                    }
                }, sshash->data());
            }
        }
    }

    std::vector<sdsl::bit_vector> coverage(2, sdsl::bit_vector(query_.get_query(false).size()));
    for (auto &cur_max_seeds : max_seeds) {
        for (auto &bucket : cur_max_seeds) {
            for (auto it = bucket.begin(); it != bucket.end(); ++it) {
                auto &cov = coverage[it->second.get_orientation()];
                std::fill(cov.begin() + it->second.get_clipping(),
                          cov.begin() + it->second.get_clipping() + it->second.get_seed().size(),
                          true);
                anchors.emplace_back(std::move(it.value()));
                assert(anchors.back().get_path_spelling().find(boss::BOSS::kSentinel) == std::string::npos);
                assert(anchors.back().get_path_spelling().size() == graph.get_k());
                assert(anchors.back().get_spelling().size() + anchors.back().get_end_trim() == graph.get_k());
            }
        }
    }

    size_t num_covered = std::max(sdsl::util::cnt_one_bits(coverage[false]),
                                  sdsl::util::cnt_one_bits(coverage[true]));

    if (static_cast<double>(num_covered) / query_.get_query(false).size() < config_.min_exact_match)
        anchors.clear();

    return anchors;
}

std::vector<Anchor> LabeledSeeder::get_anchors() const {
    auto anchors = ExactSeeder::get_anchors();

    common::logger->trace("Annotating anchors");

    for (const auto &anchor : anchors) {
        assert(anchor.get_path_spelling().find(boss::BOSS::kSentinel) == std::string::npos);
        anno_buffer_.queue_path(anchor.get_path());
    }

    anno_buffer_.fetch_queued_annotations();

    tsl::hopscotch_map<AnnotationBuffer::Column, sdsl::bit_vector> coverages;
    for (auto &anchor : anchors) {
        assert(anchor.get_path().size() == 1);
        if (const auto *labels_it = anno_buffer_.get_labels(anchor.get_path()[0])) {
            for (AnnotationBuffer::Column c : *labels_it) {
                auto [it, inserted] = coverages.try_emplace(c, sdsl::bit_vector());
                if (inserted)
                    it.value() = sdsl::bit_vector(query_.get_query(false).size());

                std::fill(it.value().begin() + anchor.get_clipping(),
                          it.value().begin() + anchor.get_clipping() + anchor.get_seed().size(),
                          true);
            }
        }
    }

    for (auto it = coverages.begin(); it != coverages.end(); ++it) {
        size_t num_covered = sdsl::util::cnt_one_bits(it->second);

        if (static_cast<double>(num_covered) / query_.get_query(false).size() < config_.min_exact_match) {
            // if the coverage is too low, clear out the vector
            it.value() = sdsl::bit_vector();
        } else {
            // otherwise, we don't need the coverage anymore, so just make it tiny
            it.value() = sdsl::bit_vector(1);
        }
    }

    std::vector<Anchor> labeled_anchors;
    labeled_anchors.reserve(anchors.size());
    for (auto &anchor : anchors) {
        auto [labels_it, coord_it] = anno_buffer_.get_labels_and_coords(anchor.get_path()[0]);
        if (labels_it) {
            for (size_t i = 0; i < labels_it->size(); ++i) {
                if (coverages[(*labels_it)[i]].size()) {
                    // only pick labels that have sufficient coverage
                    anchor.set_label_class(anno_buffer_.cache_column_set(labels_it->begin() + i,
                                                                         labels_it->begin() + i + 1));
                    if (coord_it) {
                        const auto &coords = (*coord_it)[i];
                        for (auto coord : coords) {
                            labeled_anchors.emplace_back(anchor).set_coord(coord);
                            assert(labeled_anchors.back().get_path_spelling().find(boss::BOSS::kSentinel) == std::string::npos);
                        }
                    } else {
                        labeled_anchors.emplace_back(anchor);
                        assert(labeled_anchors.back().get_path_spelling().find(boss::BOSS::kSentinel) == std::string::npos);
                    }
                }
            }
        }
        anchor = Anchor();
    }

    return labeled_anchors;
}

template <typename T>
class OffsetVector {
  public:
    OffsetVector() : offset_(0) {}

    template <typename... Args>
    OffsetVector(size_t offset, Args&&... args) : offset_(offset), null_(T()), data_(std::forward<Args>(args)...) {}

    size_t size() const { return offset_ + data_.size(); }
    bool empty() const { return offset_ == 0 && data_.empty(); }

    void resize(size_t new_size, const T &val = T()) {
        if (new_size >= offset_) {
            data_.resize(new_size - offset_, val);
            return;
        }

        new_size -= data_.size();
        data_.clear();
        assert(new_size <= offset_);
        offset_ -= new_size;
    }

    size_t offset() const { return offset_; }

    template <typename... Args>
    T& set(size_t i, Args&&... args) {
        get(i) = T(std::forward<Args>(args)...);
        return data_[i - offset_];
    }

    T& get(size_t i) {
        if (i < offset_) {
            data_.insert(data_.begin(), offset_ - i, null_);
            offset_ = i;
        }

        if (i - offset_ >= data_.size())
            data_.resize(i - offset_ + 1);

        return data_[i - offset_];
    }

    T& operator[](size_t i) {
        assert(i >= offset_);
        return data_[i - offset_];
    }

    const T& operator[](size_t i) const {
        return i >= offset_ ? data_[i - offset_] : null_;
    }

    using const_iterator = typename std::vector<T>::const_iterator;
    const_iterator begin() const { return data_.begin(); }
    const_iterator end() const { return data_.end(); }
    const_iterator cbegin() const { return data_.cbegin(); }
    const_iterator cend() const { return data_.cend(); }

  private:
    size_t offset_;
    const T null_;
    std::vector<T> data_;
};

template <typename Tuple>
using WaveFront = OffsetVector<VectorMap<DeBruijnGraph::node_index, Tuple>>; // S(cost)

template <typename Tuple>
using ScoreTable = std::vector<WaveFront<Tuple>>;

DBGAlignerConfig::score_t cost_to_score(size_t cost,
                                        size_t query_size,
                                        size_t match_size,
                                        DBGAlignerConfig::score_t match_score) {
    DBGAlignerConfig::score_t double_score = match_score * (query_size + match_size) - cost;
    assert(double_score % 2 == 0);
    return double_score / 2;
}

template <typename Tuple>
const Tuple* get_bucket(const ScoreTable<Tuple> &table,
                        size_t cost,
                        size_t query_dist,
                        DeBruijnGraph::node_index node) {
    if (cost >= table.size())
        return nullptr;

    if (query_dist >= table[cost].size() || query_dist < table[cost].offset())
        return nullptr;

    const auto &table_slot = table[cost][query_dist];
    auto it = table_slot.find(node);
    if (it == table_slot.end())
        return nullptr;

    return &it->second;
};

// dist, num_ops, last_node, last_op, last_char_of_this_node, num_matches
using SMap = std::tuple<size_t, size_t, DeBruijnGraph::node_index, Cigar::Operator, char, size_t, Anchor::label_class_t>;

// dist, num_ops
using EMap = std::tuple<size_t, size_t>;

// dist, num_ops, last_node
using FMap = std::tuple<size_t, size_t, DeBruijnGraph::node_index, char, Anchor::label_class_t>;

template <typename StrItr>
void align_impl(const std::function<size_t(DeBruijnGraph::node_index, size_t, Anchor::label_class_t)> &has_single_incoming,
                const std::function<void(DeBruijnGraph::node_index, StrItr, StrItr, size_t, Anchor::label_class_t, const TraverseCallback&, const Terminator&)> &traverse_back,
                const std::function<void(DeBruijnGraph::node_index, const EdgeCallback&, size_t, Anchor::label_class_t)> &call_incoming_kmers,
                const DBGAlignerConfig &config,
                const DeBruijnGraph::node_index start_node,
                size_t start_num_matches,
                const StrItr query_window_begin,
                const StrItr query_window_end,
                size_t max_dist,
                const std::function<void(std::string&&, std::vector<DeBruijnGraph::node_index>&&, Cigar&&, size_t, Anchor::label_class_t)> &callback,
                const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
                const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
                const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate,
                Anchor::label_class_t start_target) {
    using node_index = DeBruijnGraph::node_index;

    if (query_window_begin == query_window_end)
        return;

    size_t query_size = query_window_end - query_window_begin;

    const auto query_window_rbegin = std::make_reverse_iterator(query_window_end);
    const auto query_window_rend = std::make_reverse_iterator(query_window_begin);

    // derived from 2.4.1 in https://doi.org/10.1101/2022.01.12.476087 (v1)
    DBGAlignerConfig::score_t match_score = config.match_score("A");
    DBGAlignerConfig::score_t mismatch_score = config.score_sequences("A", "T");

    assert(config.gap_opening_penalty <= config.gap_extension_penalty);
    assert(config.gap_extension_penalty < 0);
    assert(mismatch_score < match_score);

    ssize_t mismatch_cost = (match_score - mismatch_score) * 2;
    ssize_t gap_ext = match_score - 2 * config.gap_extension_penalty;
    ssize_t gap_opn = match_score - 2 * config.gap_opening_penalty;

    assert(mismatch_cost > 0);
    assert(gap_ext > 0);
    assert(gap_opn >= gap_ext);

    assert(mismatch_cost % 2 == 0);
    assert(gap_ext % 2 == 0);
    assert(gap_opn % 2 == 0);

    common::logger->trace("x: {}\to: {}\te: {}", mismatch_cost, gap_opn, gap_ext);

    // S(cost, query_dist, node) = best_dist
    ScoreTable<SMap> S; // best cost
    ScoreTable<EMap> E; // best cost (last operation is insertion)
    ScoreTable<FMap> F; // best cost (last operation is deletion)

    auto set_value = [&](
            auto &table,
            size_t cost,
            size_t query_dist,
            node_index node,
            size_t dist,
            node_index last_node,
            size_t num_ops,
            Cigar::Operator last_op,
            char c,
            size_t num_matches,
            Anchor::label_class_t target) -> size_t {
        assert(target);
        assert(last_op == Cigar::MATCH
            || last_op == Cigar::MISMATCH
            || last_op == Cigar::INSERTION
            || last_op == Cigar::DELETION);
        bool inserted = false;
        if (cost >= table.size()) {
            inserted = true;
            table.resize(cost + 1);
        }

        if (query_dist >= table[cost].size()) {
            inserted = true;
            table[cost].resize(query_dist + 1);
        }

        auto &bucket = table[cost].get(query_dist)[node];
        using table_t = std::decay_t<decltype(bucket)>;
        static_assert(std::is_same_v<table_t, FMap>
                        || std::is_same_v<table_t, EMap>
                        || std::is_same_v<table_t, SMap>);

        inserted |= dist > std::get<0>(bucket);

        if constexpr(std::is_same_v<table_t, SMap>) {
            inserted |= std::get<3>(bucket) == Cigar::CLIPPED || std::make_pair(dist, num_matches) > std::make_pair(std::get<0>(bucket), std::get<5>(bucket));
        }

        if constexpr(std::is_same_v<table_t, FMap>)
            inserted |= std::get<2>(bucket) == DeBruijnGraph::npos;

        if constexpr(std::is_same_v<table_t, EMap>)
            inserted |= std::get<1>(bucket) == 0;

        if (inserted) {
            std::get<0>(bucket) = dist;
            std::get<1>(bucket) = num_ops;

            if constexpr(std::is_same_v<table_t, FMap>
                            || std::is_same_v<table_t, SMap>) {
                assert(target);
                assert(dist >= num_ops);
                std::get<2>(bucket) = last_node;
            }

            if constexpr(std::is_same_v<table_t, FMap>) {
                assert(c != '\0');
                std::get<3>(bucket) = c;
                std::get<4>(bucket) = target;
            }

            if constexpr(std::is_same_v<table_t, SMap>) {
                std::get<3>(bucket) = last_op;
                std::get<4>(bucket) = c;
                std::get<5>(bucket) = num_matches;
                std::get<6>(bucket) = target;
            }

            inserted = true;
        }

        if constexpr(std::is_same_v<table_t, SMap>) {
            assert(std::get<3>(bucket) == Cigar::MATCH
                || std::get<3>(bucket) == Cigar::MISMATCH
                || std::get<3>(bucket) == Cigar::INSERTION
                || std::get<3>(bucket) == Cigar::DELETION);
        }

        return inserted;
    };

    size_t num_explored_nodes = 0;
    size_t max_explored_dist = 0;
    {
        size_t query_dist = 0;
        DeBruijnGraph::node_index node = start_node;
        SMap data(0, 0, start_node, Cigar::MATCH, *(query_window_rbegin - 1), start_num_matches, start_target);
        auto &[best_dist, last_ext, last_node, last_op, c, num_matches, cur_target] = data;

        // std::cerr << "start ext" << std::endl;

        if (!has_single_incoming(node, best_dist, start_target)) {
            set_value(
                S,
                0,
                query_dist,
                node,
                best_dist,
                start_node,
                best_dist,
                Cigar::MATCH,
                best_dist > 0 ? *query_window_rbegin : '\0',
                num_matches,
                cur_target
            );
        }

        traverse_back(node, query_window_begin, query_window_end, best_dist, cur_target,
            [&](DeBruijnGraph::node_index prev, AlignmentGraph::label_class_t prev_target) {
                auto &[best_dist, last_ext, last_node, last_op, c, num_matches, cur_target] = data;
                ++num_explored_nodes;
                ++max_explored_dist;
                ++best_dist;
                ++query_dist;
                ++last_ext;
                ++num_matches;
                node = prev;
                cur_target = prev_target;

                if (!has_single_incoming(node, best_dist, cur_target)) {
                    set_value(
                        S,
                        0,
                        query_dist,
                        node,
                        best_dist,
                        start_node,
                        best_dist,
                        Cigar::MATCH,
                        best_dist > 0 ? *query_window_rbegin : '\0',
                        num_matches,
                        cur_target
                    );
                }
            },
            [&]() { return terminate_branch(0, data, query_dist, node) || std::get<0>(data) == max_dist; }
        );

        // std::cerr << "end ext" << std::endl;

        set_value(
            S,
            0,
            query_dist,
            node,
            best_dist,
            start_node,
            best_dist,
            Cigar::MATCH,
            best_dist > 0 ? *query_window_rbegin : '\0',
            num_matches,
            cur_target
        );

        // std::cerr << "init ext\n";
    }

    // size_t max_edit_cost = std::max({ gap_opn, gap_ext, mismatch_cost });

    // for (size_t cost = 0; cost < S.size() + max_edit_cost; ++cost) {
    for (size_t cost = 0; cost < S.size(); ++cost) {
        bool done = false;
        using HeapEntry = std::tuple<size_t, size_t, DeBruijnGraph::node_index>;
        std::priority_queue<HeapEntry> queue;
        for (size_t query_dist = S[cost].offset(); query_dist < S[cost].size(); ++query_dist) {
            for (const auto &[node, data] : S[cost][query_dist]) {
                size_t best_dist = std::get<0>(data);
                queue.emplace(best_dist, query_dist, node);
            }
        }

        size_t best_dist;
        size_t query_dist;
        DeBruijnGraph::node_index node;
        while (queue.size()) {
            std::tie(best_dist, query_dist, node) = queue.top();
            queue.pop();
            assert(query_dist <= query_size);

            auto it = S[cost][query_dist].find(node);
            if (it != S[cost][query_dist].end()) {
                // outdated entry
                if (std::get<0>(it->second) != best_dist)
                    continue;

                // node_index node = it->first;
                // size_t best_dist = std::get<0>(it->second);
                Cigar::Operator last_op = std::get<3>(it->second);
                size_t num_matches = std::get<5>(it->second);
                Anchor::label_class_t cur_target = std::get<6>(it->second);
                assert(last_op != Cigar::CLIPPED);

                // std::cerr << "Exp: " << cost << "\t" << node << "\t" << best_dist << "," << query_dist << "\n";

                if (terminate(cost, it->second, query_dist, node)) {
                    done = true;
                    continue;
                    // common::logger->info("Halted at cost {}", cost);
                    // break;
                }

                if (terminate_branch(cost, it->second, query_dist, node))
                    continue;

                // forward creation of insertions
                if (!done && query_dist < query_size) {
                    // extend a previous insertion
                    if (const auto *ins_ext_bucket = get_bucket(E, cost, query_dist, node)) {
                        auto [last_dist, last_num_ops] = *ins_ext_bucket;
                        assert(last_num_ops > 0);
                        assert(query_dist >= last_num_ops);
                        ssize_t next_ext_cost = cost + gap_ext;
                        if (!terminate_branch(next_ext_cost, SMap(last_dist, 0, node, Cigar::INSERTION, '\0', 0, Anchor::nannot), query_dist + 1, node)
                                && set_value(E, next_ext_cost, query_dist + 1, node, last_dist, node, last_num_ops + 1, Cigar::INSERTION, '\0', 0, Anchor::nannot)
                                && set_value(S, next_ext_cost, query_dist + 1, node, last_dist, node, 0, Cigar::INSERTION, '\0', 0, cur_target)) {
                            // it = S[cost][query_dist].begin() + it_dist;
                        }
                    }

                    // open an insertion
                    if (last_op != Cigar::DELETION) {
                        ssize_t next_opn_cost = cost + gap_opn;
                        if (!terminate_branch(next_opn_cost, SMap(best_dist, 0, node, Cigar::INSERTION, '\0', 0, Anchor::nannot), query_dist + 1, node)
                                && set_value(E, next_opn_cost, query_dist + 1, node, best_dist, node, 1, Cigar::INSERTION, '\0', 0, Anchor::nannot)
                                && set_value(S, next_opn_cost, query_dist + 1, node, best_dist, node, 0, Cigar::INSERTION, '\0', 0, cur_target)) {
                            // it = S[cost][query_dist].begin() + it_dist;
                        }
                    }
                }

                if (best_dist < max_dist) {
                    std::vector<std::tuple<node_index, char, AlignmentGraph::label_class_t>> prevs;
                    call_incoming_kmers(node, [&](auto prev, char c, AlignmentGraph::label_class_t prev_target) {
                        ++num_explored_nodes;
                        if (c != boss::BOSS::kSentinel)
                            prevs.emplace_back(prev, c, prev_target);
                    }, best_dist, cur_target);

                    if (prevs.empty())
                        continue;

                    max_explored_dist = std::max(max_explored_dist, best_dist + 1);

                    // deletion extension
                    if (!done) {
                        if (const auto *del_ext_bucket = get_bucket(F, cost, query_dist, node)) {
                            auto [last_dist, last_num_ops, last_node, del_char, last_target] = *del_ext_bucket;
                            assert(last_node != DeBruijnGraph::npos);
                            assert(last_num_ops > 0);
                            assert(last_dist >= last_num_ops);
                            if (last_dist < max_dist) {
                                ssize_t next_ext_cost = cost + gap_ext;
                                std::vector<std::tuple<DeBruijnGraph::node_index, char, AlignmentGraph::label_class_t>> last_prevs;
                                if (best_dist == last_dist) {
                                    last_prevs = prevs;
                                } else {
                                    call_incoming_kmers(node, [&](auto prev, char c, auto prev_target) {
                                        ++num_explored_nodes;
                                        if (c != boss::BOSS::kSentinel)
                                            last_prevs.emplace_back(prev, c, prev_target);
                                    }, last_dist, last_target);
                                }

                                for (const auto &[prev, c, prev_target] : last_prevs) {
                                    if (!terminate_branch(next_ext_cost, SMap(last_dist + 1, 0, node, Cigar::DELETION, c, 0, prev_target), query_dist, prev)
                                            && set_value(F, next_ext_cost, query_dist, prev, last_dist + 1, node, last_num_ops + 1, Cigar::DELETION, c, 0, prev_target)
                                            && set_value(S, next_ext_cost, query_dist, prev, last_dist + 1, node, 0, Cigar::DELETION, c, 0, prev_target)) {
                                        // it = S[cost][query_dist].begin() + it_dist;
                                    }
                                }
                            }
                        }

                        // deletion open
                        if (last_op != Cigar::INSERTION) {
                            ssize_t next_opn_cost = cost + gap_opn;
                            for (const auto &[prev, c, prev_target] : prevs) {
                                if (!terminate_branch(next_opn_cost, SMap(best_dist + 1, 0, node, Cigar::DELETION, c, 0, prev_target), query_dist, prev)
                                        && set_value(F, next_opn_cost, query_dist, prev, best_dist + 1, node, 1, Cigar::DELETION, c, 0, prev_target)
                                        && set_value(S, next_opn_cost, query_dist, prev, best_dist + 1, node, 0, Cigar::DELETION, c, 0, prev_target)) {
                                    // it = S[cost][query_dist].begin() + it_dist;
                                }
                            }
                        }
                    }

                    if (query_dist < query_size) {
                        // match
                        auto local_query_window_begin = query_window_begin;
                        auto local_query_window_end = query_window_end;
                        local_query_window_end -= query_dist;

                        std::reverse_iterator<StrItr> local_query_window_rbegin = std::make_reverse_iterator(local_query_window_end);
                        std::reverse_iterator<StrItr> local_query_window_rend = std::make_reverse_iterator(local_query_window_begin);

                        for (size_t i = 0; i < prevs.size(); ++i) {
                            DeBruijnGraph::node_index prev = std::get<0>(prevs[i]);
                            char c = std::get<1>(prevs[i]);
                            AlignmentGraph::label_class_t prev_target = std::get<2>(prevs[i]);

                            SMap data(best_dist + 1, 1, node, Cigar::MATCH, c, num_matches + 1, prev_target);
                            auto &[cur_best, num_ops, last_node, op, cur_c, cur_num_matches, cur_target] = data;
                            size_t cur_query_dist = query_dist + 1;
                            size_t cur_cost = cost;
                            if (c != *local_query_window_rbegin) {
                                cur_cost += mismatch_cost;
                                op = Cigar::MISMATCH;
                                cur_num_matches = 0;
                            } else if (last_op == Cigar::MATCH) {
                                continue;
                            }

                            if (!has_single_incoming(prev, cur_best, cur_target)) {
                                if (set_value(
                                    S,
                                    cur_cost,
                                    cur_query_dist,
                                    prev,
                                    cur_best,
                                    node,
                                    cur_query_dist - query_dist,
                                    op,
                                    c,
                                    cur_num_matches,
                                    cur_target
                                )) {
                                    if (cur_cost == cost) {
                                        queue.emplace(cur_best, cur_query_dist, prev);
                                    }
                                    // it = S[cost][query_dist].begin() + it_dist;
                                }
                            }

                            traverse_back(
                                prev,
                                local_query_window_rend.base(),
                                (local_query_window_rbegin + 1).base(),
                                cur_best, cur_target,
                                [&](DeBruijnGraph::node_index pprev, AlignmentGraph::label_class_t pprev_target) {
                                    auto &[cur_best, num_ops, last_node, op, cur_c, cur_num_matches, cur_target] = data;
                                    prev = pprev;
                                    ++num_explored_nodes;
                                    ++cur_best;
                                    max_explored_dist = std::max(max_explored_dist, cur_best);
                                    ++cur_query_dist;
                                    ++cur_num_matches;
                                    cur_target = pprev_target;

                                    if (!has_single_incoming(prev, cur_best, cur_target)) {
                                        if (set_value(
                                            S,
                                            cur_cost,
                                            cur_query_dist,
                                            prev,
                                            cur_best,
                                            node,
                                            cur_query_dist - query_dist,
                                            op,
                                            c,
                                            cur_num_matches,
                                            cur_target
                                        )) {
                                            if (cur_cost == cost) {
                                                queue.emplace(cur_best, cur_query_dist, prev);
                                            }
                                            // it = S[cost][query_dist].begin() + it_dist;
                                        }
                                    }
                                },
                                [&]() { return terminate_branch(cur_cost, data, cur_query_dist, prev) || std::get<0>(data) == max_dist; }
                            );

                            if (set_value(
                                S,
                                cur_cost,
                                cur_query_dist,
                                prev,
                                cur_best,
                                node,
                                cur_query_dist - query_dist,
                                op,
                                c,
                                cur_num_matches,
                                cur_target
                            )) {
                                if (cur_cost == cost) {
                                    queue.emplace(cur_best, cur_query_dist, prev);
                                }
                                // it = S[cost][query_dist].begin() + it_dist;
                            }
                        }
                    }
                }
            }
            // if (done)
            //     break;
        }

        if (done)
            break;
    }

    // backtrack
    for (size_t start_cost = 0; start_cost < S.size(); ++start_cost) {
        std::vector<std::tuple<size_t, ssize_t, size_t, DeBruijnGraph::node_index>> starts;
        for (size_t start_query_dist = S[start_cost].offset(); start_query_dist < S[start_cost].size(); ++start_query_dist) {
            assert(start_query_dist <= query_size);
            const auto &bucket = S[start_cost][start_query_dist];
            for (const auto &[node, data] : bucket) {
                const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches, cur_target] = data;
                assert(last_op != Cigar::CLIPPED);
                // std::cerr << "BTIter: " << start_cost << "\t" << node << "\t" << dist << "," << start_query_dist << std::endl;
                if (start_backtrack(start_cost, data, start_query_dist, node))
                    starts.emplace_back(dist, static_cast<ssize_t>(dist) - start_query_dist, start_query_dist, node);
            }
        }

        auto traverse_back_char = [&traverse_back](DeBruijnGraph::node_index node, char c, size_t dist, Anchor::label_class_t cur_target) {
            std::string c_str(1, c);
            StrItr begin;
            StrItr end;
            static_assert(std::is_same_v<StrItr,const char*>
                            || std::is_same_v<StrItr,std::reverse_iterator<const char*>>);
            if constexpr(std::is_same_v<StrItr,const char*>) {
                begin = c_str.c_str();
                end = begin + 1;
            } else if constexpr(std::is_same_v<StrItr,std::reverse_iterator<const char*>>) {
                begin = std::make_reverse_iterator(c_str.c_str() + 1);
                end = std::make_reverse_iterator(c_str.c_str());
            }

            DeBruijnGraph::node_index prev = DeBruijnGraph::npos;
            traverse_back(node, begin, end, dist, cur_target,
                          [&](DeBruijnGraph::node_index pprev, AlignmentGraph::label_class_t prev_target) {
                              prev = pprev;
                              cur_target = prev_target;
                          },
                          []() { return false; });
            return std::make_pair(prev, cur_target);
        };

        for (auto [dist, diag, query_dist, node] : starts) {
            size_t cost = start_cost;

            const auto *bt = get_bucket(S, cost, query_dist, node);
            assert(bt);

            std::vector<node_index> path;
            Cigar cigar;
            std::string spelling;
            Anchor::label_class_t final_target = std::get<6>(*bt);
            Anchor::label_class_t cur_target = final_target;

            do {
                assert(std::get<0>(*bt) == dist);

                if (std::get<3>(*bt) == Cigar::INSERTION) {
                    assert(std::get<1>(*bt) == 0);
                    assert(cost > 0);
                    const auto *et = get_bucket(E, cost, query_dist, node);
                    assert(et);

                    auto [ins_dist, num_ins] = *et;
                    assert(query_dist >= num_ins);
                    assert(num_ins > 0);
                    assert(ins_dist == dist);

                    size_t prev_cost = cost - gap_opn - (num_ins - 1) * gap_ext;
                    assert(prev_cost < cost);
                    size_t prev_query_dist = query_dist - num_ins;

                    const auto *check_bt = get_bucket(S, prev_cost, prev_query_dist, node);
                    assert(check_bt);
                    assert(std::get<0>(*check_bt) == dist);
                    assert(std::get<3>(*check_bt) != Cigar::CLIPPED);
                    assert(std::get<6>(*check_bt) == cur_target);

                    query_dist = prev_query_dist;
                    cost = prev_cost;
                    bt = check_bt;
                    cigar.append(Cigar::INSERTION, num_ins);
                } else if (std::get<3>(*bt) == Cigar::DELETION) {
                    assert(std::get<1>(*bt) == 0);
                    assert(cost > 0);
                    const auto *ft = get_bucket(F, cost, query_dist, node);
                    assert(ft);

                    auto [del_dist, num_del, prev_node, del_char, del_target] = *ft;
                    assert(prev_node != DeBruijnGraph::npos);
                    assert(del_char != '\0');
                    assert(del_dist == dist);
                    assert(dist >= num_del);
                    assert(del_target == cur_target);

                    cigar.append(Cigar::DELETION, num_del);

                    size_t prev_cost = cost;
                    size_t prev_dist = dist;
                    while (num_del > 1) {
                        path.emplace_back(node);
                        spelling += del_char;

                        prev_cost -= gap_ext;
                        --prev_dist;
                        node = prev_node;
                        --num_del;

                        const auto *prev_ft = get_bucket(F, prev_cost, query_dist, node);
                        assert(prev_ft);
                        assert(std::get<0>(*prev_ft) == prev_dist);
                        assert(std::get<1>(*prev_ft) == num_del);
                        ft = prev_ft;
                        prev_node = std::get<2>(*ft);
                        del_char = std::get<3>(*ft);
                        cur_target = std::get<4>(*ft);
                    }

                    path.emplace_back(node);
                    spelling += del_char;

                    prev_cost -= gap_opn;
                    --prev_dist;

                    assert(prev_cost < cost);

                    const auto *check_bt = get_bucket(S, prev_cost, query_dist, prev_node);
                    assert(check_bt);
                    assert(std::get<0>(*check_bt) == prev_dist);
                    assert(std::get<3>(*check_bt) != Cigar::CLIPPED);
                    assert(std::get<6>(*check_bt) == cur_target);

                    dist = prev_dist;
                    cost = prev_cost;
                    node = prev_node;
                    bt = check_bt;
                } else if (std::get<3>(*bt) == Cigar::MATCH || std::get<3>(*bt) == Cigar::MISMATCH) {
                    auto [cur_dist, num_match, last_node, last_op, mismatch_char, num_matches, mm_target] = *bt;

                    if (!num_match) {
                        assert(last_op == Cigar::MATCH);
                        bt = nullptr;
                        break;
                    }

                    assert(query_dist >= num_match);
                    assert(dist >= num_match);
                    assert(cur_dist == dist);
                    assert(cur_target == mm_target);

                    node_index traverse_node = last_node;

                    cigar.append(Cigar::MATCH, num_match - 1);
                    cigar.append(last_op);
                    size_t prev_dist = dist - num_match;
                    size_t prev_query = query_dist - num_match;

                    // reconstruct the backwards traversal from last_node to node
                    auto it = query_window_rbegin + prev_query;
                    std::ignore = query_window_rend;
                    assert(it != query_window_rend);
                    size_t cur_path_size = path.size();
                    size_t num_traversed = 0;
                    if (std::get<3>(*bt) == Cigar::MISMATCH) {
                        assert(mismatch_char != *it);
                        std::tie(traverse_node, cur_target) = traverse_back_char(traverse_node, mismatch_char, prev_dist, cur_target);
                        assert(traverse_node != DeBruijnGraph::npos);
                        path.emplace_back(traverse_node);
                        spelling += mismatch_char;
                        ++it;
                        --num_match;
                        ++num_traversed;
                    } else {
                        assert(mismatch_char == *it);
                    }

                    std::reverse_iterator<StrItr> it_end = it + num_match;
                    assert(it_end <= query_window_rend);
                    size_t old_path_size = path.size();
                    traverse_back(traverse_node, it_end.base(), it.base(),
                                  prev_dist + num_traversed,
                                  cur_target,
                                  [&](DeBruijnGraph::node_index prev, AlignmentGraph::label_class_t prev_target) {
                                      traverse_node = prev;
                                      cur_target = prev_target;
                                      path.emplace_back(traverse_node);
                                  },
                                  []() { return false; });
                    std::ignore = old_path_size;
                    assert(path.size() == old_path_size + num_match);
                    assert(traverse_node == node);
                    spelling += std::string(it, it_end);
                    num_traversed += num_match;
                    num_match = 0;

                    std::reverse(path.begin() + cur_path_size, path.end());
                    std::reverse(spelling.begin() + cur_path_size, spelling.end());

                    size_t prev_cost = cost - (last_op == Cigar::MATCH ? 0 : mismatch_cost);
                    assert(prev_cost <= cost);

                    const auto *check_bt = bt;
                    check_bt = nullptr;
                    if (prev_cost || prev_query || prev_dist) {
                        check_bt = get_bucket(S, prev_cost, prev_query, last_node);
                        assert(check_bt);
                        assert(std::get<0>(*check_bt) == prev_dist);
                        assert(std::get<3>(*check_bt) != Cigar::CLIPPED);
                        assert(std::get<6>(*check_bt) == cur_target);
                    } else {
                        assert(last_node == start_node);
                    }

                    dist = prev_dist;
                    cost = prev_cost;
                    query_dist = prev_query;
                    node = last_node;

                    if (!check_bt)
                        break;

                    bt = check_bt;
                } else {
                    assert(false && "Different operation found");
                    throw std::runtime_error("Inf. loop");
                }

            } while (bt);

            assert(dist == 0);
            assert(cost == 0);
            assert(query_dist == 0);
            // std::cerr << "Ext: " << start_cost << "\t" << cigar.to_string() << "\n";
            callback(std::move(spelling), std::move(path), std::move(cigar), start_cost, final_target);
        }
    }
    common::logger->trace("Explored {} nodes. q: {}\td: {} / {}",
                          num_explored_nodes, query_size, max_explored_dist, max_dist);
}

template <class T1, class T2>
bool overlap_with_diff(const T1 &tuple1, const T2 &tuple2, ssize_t diff) {
    auto a_begin = tuple1.begin();
    const auto a_end = tuple1.end();
    auto b_begin = tuple2.begin();
    const auto b_end = tuple2.end();

    assert(std::is_sorted(a_begin, a_end));
    assert(std::is_sorted(b_begin, b_end));

    while (a_begin != a_end && b_begin != b_end) {
        if (*a_begin + diff == *b_begin)
            return true;

        if (*a_begin + diff < *b_begin) {
            ++a_begin;
        } else {
            ++b_begin;
        }
    }

    return false;
}

Anchor::label_class_t subset_coords(AnnotationBuffer &anno_buffer,
                                    Anchor::label_class_t source,
                                    const AnnotationBuffer::CoordinateSet &source_coords,
                                    Anchor::label_class_t target,
                                    DeBruijnGraph::node_index target_node,
                                    ssize_t coord_offset) {
    assert(source);
    assert(target);
    assert(source != AnnotationBuffer::nannot);
    assert(target != AnnotationBuffer::nannot);

    const auto &source_columns = anno_buffer.get_cached_column_set(source);
    assert(source_columns.size() == source_coords.size());

    const auto &target_columns = anno_buffer.get_cached_column_set(target);
    const auto [target_columns_full, target_coords_full] = anno_buffer.get_labels_and_coords(target_node);
    assert(target_columns_full);
    assert(target_columns_full->size() >= target_columns.size());
    assert(target_coords_full);
    assert(target_coords_full->size() == target_columns_full->size());

#ifndef NDEBUG
    auto it = target_columns.begin();
#endif

    AnnotationBuffer::Columns intersection;
    utils::match_indexed_values(
        source_columns.begin(), source_columns.end(), source_coords.begin(),
        target_columns_full->begin(), target_columns_full->end(), target_coords_full->begin(),
        [&](AnnotationBuffer::Column c, const AnnotationBuffer::Tuple &source_coords, const AnnotationBuffer::Tuple &target_coords) {
            assert(source_coords.size());
            assert(target_coords.size());
#ifndef NDEBUG
            it = std::find(it, target_columns.end(), c);
            assert(it != target_columns.end());
#endif
            if (overlap_with_diff(source_coords, target_coords, coord_offset))
                intersection.emplace_back(c);
        }
    );

    assert(it + 1 == target_columns.end());

    if (intersection.size() == target_columns.size())
        return target;

    return anno_buffer.cache_column_set(std::move(intersection));
}

void align_bwd(const DeBruijnGraph &base_graph,
               const DBGAlignerConfig &config,
               const DeBruijnGraph::node_index start_node,
               size_t num_matches,
               std::string_view query_window,
               size_t max_dist,
               const std::function<void(std::string&&, std::vector<DeBruijnGraph::node_index>&&, Cigar&&, size_t, Anchor::label_class_t)> &callback,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate,
               AnnotationBuffer *anno_buffer = nullptr,
               Anchor::label_class_t target = Anchor::nannot,
               const AnnotationBuffer::CoordinateSet &coords = {}) {
    assert(!anno_buffer == (target == Anchor::nannot));
    assert(!anno_buffer || anno_buffer->node_has_labels(start_node, target));
    assert(!anno_buffer
            || (!anno_buffer->has_coordinates() && coords.empty())
            || anno_buffer->get_cached_column_set(target).size() == coords.size());

    using StrItr = std::string_view::iterator;

    align_impl<StrItr>(
        [&](DeBruijnGraph::node_index node, size_t dist, Anchor::label_class_t cur_target) {
            assert(cur_target);
            AlignmentGraph aln_graph(base_graph, anno_buffer, cur_target);
            if (coords.empty()) {
                return aln_graph.has_single_incoming(node);
            } else {
                assert(anno_buffer);
                size_t indegree = 0;
                size_t next_dist = dist + 1;
                aln_graph.call_incoming_kmers(
                    node,
                    [&](DeBruijnGraph::node_index prev, char, Anchor::label_class_t prev_target) {
                        if (subset_coords(*anno_buffer, target, coords, prev_target, prev,
                                          -static_cast<ssize_t>(next_dist)))
                            ++indegree;
                    }
                );
                return indegree == 1;
            }
        },
        [&](DeBruijnGraph::node_index node,
            StrItr begin, StrItr end,
            size_t dist, Anchor::label_class_t cur_target,
            const auto &callback, const auto &terminate) {

            assert(cur_target);
            if (begin != end) {
                AlignmentGraph aln_graph(base_graph, anno_buffer, cur_target);
                if (coords.empty()) {
                    aln_graph.traverse_back(
                        node,
                        std::string_view(&*begin, end - begin),
                        callback, terminate
                    );
                } else {
                    assert(anno_buffer);
                    bool stop_early = false;
                    size_t next_dist = dist + 1;
                    aln_graph.traverse_back(
                        node,
                        std::string_view(&*begin, end - begin),
                        [&](DeBruijnGraph::node_index prev, Anchor::label_class_t prev_target) {
                            if (subset_coords(*anno_buffer, target, coords, prev_target, prev,
                                              -static_cast<ssize_t>(next_dist++))) {
                                callback(prev, prev_target);
                            } else {
                                stop_early = true;
                            }
                        },
                        [&]() { return stop_early || terminate(); }
                    );
                }
            }
        },
        [&](DeBruijnGraph::node_index node,
            const auto &callback,
            size_t dist, Anchor::label_class_t cur_target) {

            assert(cur_target);

            AlignmentGraph aln_graph(base_graph, anno_buffer, cur_target);
            if (coords.empty()) {
                aln_graph.call_incoming_kmers(node, callback);
            } else {
                size_t next_dist = dist + 1;
                aln_graph.call_incoming_kmers(
                    node,
                    [&](DeBruijnGraph::node_index prev, char c, Anchor::label_class_t prev_target) {
                        if (Anchor::label_class_t prev_target_s = subset_coords(*anno_buffer, target, coords, prev_target, prev,
                                                                                -static_cast<ssize_t>(next_dist)))
                            callback(prev, c, prev_target_s);
                    }
                );
            }
        },
        config,
        start_node,
        num_matches,
        query_window.begin(), query_window.end(),
        max_dist,
        [&](auto&& spelling, auto&& path, auto&& cigar, size_t cost, Anchor::label_class_t final_target) {
            assert(final_target);
            assert(std::all_of(path.begin(), path.end(),
                               [&](auto node) {
                                   return AlignmentGraph(base_graph, anno_buffer, target).has_labels(node, final_target);
                                }));
            callback(std::move(spelling), std::move(path), std::move(cigar), cost, final_target);
        },
        start_backtrack,
        terminate_branch,
        terminate,
        target
    );
}

void align_fwd(const DeBruijnGraph &base_graph,
               const DBGAlignerConfig &config,
               const DeBruijnGraph::node_index start_node,
               size_t num_matches,
               std::string_view query_window,
               size_t max_dist,
               const std::function<void(std::string&&, std::vector<DeBruijnGraph::node_index>&&, Cigar&&, size_t, Anchor::label_class_t)> &callback,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &start_backtrack,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate_branch,
               const std::function<bool(size_t, const SMap&, size_t, DeBruijnGraph::node_index)> &terminate,
               std::string_view suffix = "",
               AnnotationBuffer *anno_buffer = nullptr,
               Anchor::label_class_t target = Anchor::nannot,
               const AnnotationBuffer::CoordinateSet &coords = {}) {
    assert(!anno_buffer == (target == Anchor::nannot));
    assert(!anno_buffer || anno_buffer->node_has_labels(start_node, target));
    assert(!anno_buffer
            || (!anno_buffer->has_coordinates() && coords.empty())
            || anno_buffer->get_cached_column_set(target).size() == coords.size());

    using StrItr = std::reverse_iterator<std::string_view::iterator>;

    align_impl<StrItr>(
        [&](DeBruijnGraph::node_index node, size_t dist, Anchor::label_class_t cur_target) {
            assert(cur_target);
            assert(dist > suffix.size() || node == start_node);
            if (dist < suffix.size())
                return true;

            AlignmentGraph aln_graph(base_graph, anno_buffer, cur_target);
            if (coords.empty()) {
                return aln_graph.has_single_outgoing(node);
            } else {
                assert(anno_buffer);
                size_t outdegree = 0;
                assert(dist >= suffix.size());
                size_t next_dist = dist - suffix.size() + 1;
                aln_graph.call_outgoing_kmers(
                    node,
                    [&](DeBruijnGraph::node_index next, char, Anchor::label_class_t next_target) {
                        if (subset_coords(*anno_buffer, target, coords, next_target, next, next_dist))
                            ++outdegree;
                    }
                );
                return outdegree == 1;
            }
        },
        [&](DeBruijnGraph::node_index node,
            StrItr rbegin, StrItr rend,
            size_t dist, Anchor::label_class_t cur_target,
            const auto &callback, const auto &terminate) {

            assert(cur_target);
            assert(dist > suffix.size() || node == start_node);

            auto begin = rend.base();
            auto end = rbegin.base();
            assert(begin <= end);

            for (size_t i = dist; !terminate() && i < suffix.size() && begin != end; ++i, ++begin) {
                if (*begin == suffix[i]) {
                    callback(node, cur_target);
                } else {
                    return;
                }
            }

            if (begin != end) {
                AlignmentGraph aln_graph(base_graph, anno_buffer, cur_target);
                if (coords.empty()) {
                    aln_graph.traverse(
                        node, std::string_view(&*begin, end - begin),
                        callback, terminate
                    );
                } else {
                    assert(anno_buffer);
                    bool stop_early = false;
                    assert(dist >= suffix.size());
                    size_t next_dist = dist - suffix.size() + 1;
                    aln_graph.traverse(
                        node,
                        std::string_view(&*begin, end - begin),
                        [&](DeBruijnGraph::node_index next, Anchor::label_class_t next_target) {
                            if (subset_coords(*anno_buffer, target, coords, next_target, next, next_dist++)) {
                                callback(next, next_target);
                            } else {
                                stop_early = true;
                            }
                        },
                        [&]() { return stop_early || terminate(); }
                    );
                }
            }
        },
        [&](DeBruijnGraph::node_index node, const auto &callback,
            size_t dist, Anchor::label_class_t cur_target) {

            assert(cur_target);
            if (dist < suffix.size()) {
                assert(node == start_node);
                callback(node, suffix[dist], cur_target);
            } else {
                AlignmentGraph aln_graph(base_graph, anno_buffer, cur_target);
                if (coords.empty()) {
                    aln_graph.call_outgoing_kmers(node, callback);
                } else {
                    assert(dist >= suffix.size());
                    size_t next_dist = dist - suffix.size() + 1;
                    aln_graph.call_outgoing_kmers(
                        node,
                        [&](DeBruijnGraph::node_index next, char c, Anchor::label_class_t next_target) {
                            if (Anchor::label_class_t next_target_s = subset_coords(*anno_buffer, target, coords, next_target, next, next_dist))
                                callback(next, c, next_target_s);
                        }
                    );
                }
            }
        },
        config,
        start_node,
        num_matches,
        query_window.rbegin(), query_window.rend(),
        max_dist,
        [&](auto&& spelling, auto&& path, auto&& cigar, size_t cost, Anchor::label_class_t final_target) {
            assert(final_target);
            path.resize(path.size() - std::min(path.size(), suffix.size()));
            assert(std::all_of(path.begin(), path.end(),
                               [&](auto node) {
                                   return AlignmentGraph(base_graph, anno_buffer, final_target).has_labels(node, target);
                                }));
            std::reverse(spelling.begin(), spelling.end());
            std::reverse(path.begin(), path.end());
            std::reverse(cigar.data().begin(), cigar.data().end());
            callback(std::move(spelling), std::move(path), std::move(cigar), cost, final_target);
        },
        start_backtrack,
        terminate_branch,
        terminate,
        target
    );
}

void Extender::extend(const Alignment &aln,
                      const std::function<void(Alignment&&)> &callback,
                      bool no_bwd, bool no_fwd) const {
    // std::cerr << "Base " << aln << "\n";
    DBGAlignerConfig::score_t match_score = config_.match_score("A");

    auto set_coords = [&](Alignment &aln) {
        assert(aln.get_label_classes().size());
        if (auto anno_buffer = make_aln_graph(aln.get_label_classes()[0]).get_anno_buffer()) {
            if (anno_buffer->has_coordinates()) {
                auto coords = anno_buffer->get_label_path_coords(aln.get_path(), aln.get_label_classes());
                assert(coords.size());
                aln.set_coords(std::move(coords));
            }
        }
    };

    std::vector<Alignment> fwd_exts;
    if (no_fwd || !aln.get_end_clipping()) {
        // std::cerr << "skip fwd\t" << aln << "\t" << aln.get_label_classes().back() << "\n";
        fwd_exts.emplace_back(aln);
    } else {
        // std::cerr << "fwd extending\t" << aln << "\t" << aln.get_label_classes().back() << "\n";
        std::string_view query_window = aln.get_query();
        query_window.remove_prefix(aln.get_clipping() + aln.get_seed().size());

        DBGAlignerConfig::score_t aln_score = score_match(aln, config_);
        auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
            return aln_score + cost_to_score(cost, query_dist, dist, match_score)
                    + (query_dist == query_window.size() ? config_.right_end_bonus : 0);
        };

        auto it = std::make_reverse_iterator(aln.get_cigar().data().end());
        auto end = std::make_reverse_iterator(aln.get_cigar().data().begin());
        std::ignore = end;
        assert(it != end);
        if (it->first == Cigar::CLIPPED) {
            ++it;
            assert(it != end);
        }
        size_t num_matches = it->first == Cigar::MATCH ? it->second : 0;
        size_t max_dist = 0;
        DBGAlignerConfig::score_t best_score = aln_score;
        size_t max_query_dist = 0;
        Anchor::label_class_t target = aln.get_label_classes().back();
        auto *anno_buffer = make_aln_graph(target).get_anno_buffer();
        align_fwd(
            query_.get_graph(),
            config_, aln.get_path().back(), num_matches,
            query_window, std::numeric_limits<size_t>::max(),
            [&](auto&& spelling, auto&& path, auto&& cigar, size_t, Anchor::label_class_t final_target) {
                size_t new_end_trim = aln.get_end_trim() - std::min(aln.get_end_trim(), spelling.size());

                assert(path.empty() == (spelling.size() <= aln.get_end_trim()));

                std::string new_spelling = spelling.size() > aln.get_end_trim()
                    ? std::string(aln.get_spelling()) + spelling
                    : aln.get_path_spelling();

                auto ext_path = aln.get_path();
                ext_path.insert(ext_path.end(), path.begin(), path.end());

                size_t query_dist = cigar.get_num_query();
                assert(query_dist <= query_window.size());
                cigar.append(Cigar::CLIPPED, query_window.size() - query_dist);

                Cigar ext_cigar = aln.get_cigar();
                ext_cigar.trim_end_clipping();
                ext_cigar.append(std::move(cigar));

                fwd_exts.emplace_back(
                    query_.get_graph(),
                    aln.get_query(),
                    aln.get_orientation(),
                    std::move(ext_path),
                    std::move(ext_cigar),
                    std::move(new_spelling),
                    new_end_trim,
                    final_target
                );

                set_coords(fwd_exts.back());

                // std::cerr << "\tend" << fwd_exts.back() << "\t" << fwd_exts.back().get_label_classes().back() << std::endl;

                // assert(fwd_exts.back().get_score() == best_score);
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // start backtrack
                size_t dist = std::get<0>(data);
                // size_t num_matches = std::get<5>(data);
                // if (!num_matches || dist < aln.get_end_trim() || (dist == 0 && query_dist == 0))
                // if (!num_matches || (dist == 0 && query_dist == 0))
                if (dist == 0 && query_dist == 0)
                    return false;

                // if (best_score == get_score(cost, dist, query_dist))
                    // std::cerr << "BT?: " << get_score(cost, dist, query_dist) << " vs. " << best_score << "\t" << query_dist << " vs. " << max_query_dist << std::endl;
                if (query_dist < max_query_dist)
                    return false;

                auto score = get_score(cost, dist, query_dist);
                assert(score <= best_score);
                return score == best_score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate branch
                //if (query_dist == query_window.size())
                if (query_dist > query_window.size())
                    return true;

                size_t dist = std::get<0>(data);
                auto score = get_score(cost, dist, query_dist);
                    // std::cerr << "TB?: " << score << " vs. " << best_score << "\t" << query_dist << " vs. " << max_query_dist << std::endl;
                return score <= 0 || config_.xdrop < best_score - score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate
                size_t dist = std::get<0>(data);
                max_dist = std::max(dist, max_dist);
                auto score = get_score(cost, dist, query_dist);
                if (score >= best_score) {
                    best_score = score;
                    max_query_dist = std::max(max_query_dist, query_dist);
                }
                // if (query_dist == query_window.size())
                //     std::cerr << "Exppp: " << cost << "\n";
                // return query_dist == query_window.size() && dist == max_dist;
                return query_dist == query_window.size();
            },
            aln.get_trim_spelling(),
            anno_buffer,
            target,
            anno_buffer ? anno_buffer->get_label_path_coords({ aln.get_path().back() },
                                                             { target })
                        : AnnotationBuffer::CoordinateSet()
        );

        if (fwd_exts.empty())
            fwd_exts.emplace_back(aln);
    }

    for (auto &fwd_ext : fwd_exts) {
        fwd_ext.trim_end();
        if (no_bwd || !fwd_ext.get_clipping()) {
            // std::cerr << "\tskipbwd\t" << fwd_ext << "\n";
            set_coords(fwd_ext);
            callback(std::move(fwd_ext));
            continue;
        }

        std::string_view query_window = fwd_ext.get_query();
        query_window.remove_suffix(fwd_ext.get_end_clipping() + fwd_ext.get_seed().size());

        DBGAlignerConfig::score_t fwd_ext_score = score_match(fwd_ext, config_);
        auto get_score = [&](size_t cost, size_t dist, size_t query_dist) {
            return fwd_ext_score + cost_to_score(cost, query_dist, dist, match_score)
                    + (query_dist == query_window.size() ? config_.left_end_bonus : 0);
        };

        DBGAlignerConfig::score_t best_score = fwd_ext_score;

        auto it = aln.get_cigar().data().begin();
        auto end = aln.get_cigar().data().end();
        std::ignore = end;
        assert(it != end);
        if (it->first == Cigar::CLIPPED) {
            ++it;
            assert(it != end);
        }
        size_t num_matches = it->first == Cigar::MATCH ? it->second : 0;

        bool found = false;
        size_t max_dist = 0;
        size_t max_query_dist = 0;
        // std::cerr << "bwd extending\t" << fwd_ext << "\t" << fwd_ext.get_label_classes()[0] << "\n";
        auto *anno_buffer = make_aln_graph(fwd_ext.get_label_classes()[0]).get_anno_buffer();
        align_bwd(
            query_.get_graph(),
            config_, fwd_ext.get_path()[0], num_matches,
            query_window, std::numeric_limits<size_t>::max(),
            [&](auto&& spelling, auto&& path, auto&& cigar, size_t, Anchor::label_class_t final_target) {
                path.insert(path.end(), fwd_ext.get_path().begin(), fwd_ext.get_path().end());
                size_t query_dist = cigar.get_num_query();
                assert(query_dist <= query_window.size());
                Cigar ext_cigar(Cigar::CLIPPED, query_window.size() - query_dist);
                ext_cigar.append(std::move(cigar));

                Cigar aln_cigar = fwd_ext.get_cigar();
                aln_cigar.trim_clipping();
                ext_cigar.append(std::move(aln_cigar));

                Alignment next_aln(query_.get_graph(),
                    fwd_ext.get_query(),
                    fwd_ext.get_orientation(),
                    std::move(path),
                    std::move(ext_cigar),
                    spelling + fwd_ext.get_path_spelling(),
                    fwd_ext.get_end_trim(),
                    final_target
                );
                next_aln.trim_end();

                // assert(next_aln.get_score() == best_score);
                // std::cerr << "found\t" << next_aln << "\n";
                set_coords(next_aln);
                callback(std::move(next_aln));
                found = true;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // start backtrack
                size_t dist = std::get<0>(data);
                // size_t num_matches = std::get<5>(data);
                // if ((dist == 0 && query_dist == 0) || !num_matches)
                //     return false;

                if (dist == 0 && query_dist == 0)
                    return false;

                // if (best_score == get_score(cost, dist, query_dist))
                //     std::cerr << "BT?: " << get_score(cost, dist, query_dist) << " vs. " << best_score << "\t" << query_dist << " vs. " << max_query_dist << std::endl;
                if (query_dist < max_query_dist)
                    return false;

                auto score = get_score(cost, dist, query_dist);
                assert(score <= best_score);
                return score == best_score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate branch
                //if (query_dist == query_window.size())
                if (query_dist > query_window.size())
                    return true;

                size_t dist = std::get<0>(data);
                auto score = get_score(cost, dist, query_dist);
                return score <= 0 || config_.xdrop < best_score - score;
            },
            [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                // terminate
                size_t dist = std::get<0>(data);
                max_dist = std::max(dist, max_dist);
                auto score = get_score(cost, dist, query_dist);
                if (score >= best_score) {
                    best_score = score;
                    max_query_dist = std::max(max_query_dist, query_dist);
                }
                // if (query_dist == query_window.size())
                //     std::cerr << "Exppp: " << cost << "\t" << score << " vs. " << best_score << "\t" << Cigar::opt_to_char(std::get<3>(data)) << "\n";
                // return false;
                // return query_dist == query_window.size() && dist == max_dist;
                return query_dist == query_window.size();
            },
            anno_buffer,
            fwd_ext.get_label_classes()[0],
            anno_buffer ? anno_buffer->get_label_path_coords({ fwd_ext.get_path()[0] },
                                                             { fwd_ext.get_label_classes()[0] })
                        : AnnotationBuffer::CoordinateSet()
        );

        if (!found) {
            set_coords(fwd_ext);
            callback(std::move(fwd_ext));
        }
    }
}

std::vector<Alignment> ExactSeeder::get_inexact_anchors(bool align) const {
    common::logger->trace("Seeding");
    std::vector<Anchor> anchors = get_anchors();

    const auto &graph = query_.get_graph();
    const auto *sshash = dynamic_cast<const DBGSSHash*>(&graph);

    auto full_it = std::partition(anchors.begin(), anchors.end(), [](const auto &a) { return a.get_end_trim(); });
    std::sort(full_it, anchors.end(), AnchorLess<Anchor>());
    if (anchors.end() - full_it > 1) {
        common::logger->trace("Merging {} anchors", anchors.size());
        assert(full_it->get_spelling().size() + full_it->get_end_trim() == graph.get_k());
        for (auto it = full_it; it + 1 != anchors.end(); ++it) {
            auto &a_i = *it;
            auto &a_j = *(it + 1);
            assert(!a_i.get_end_trim());
            assert(!a_j.get_end_trim());
            assert(a_j.get_spelling().size() + a_j.get_end_trim() == graph.get_k());
            if (a_i.get_orientation() != a_j.get_orientation() || a_i.get_label_class() != a_j.get_label_class())
                continue;

            if (a_i.get_seed().size() < config_.max_seed_length
                    && a_i.get_end_clipping() == a_j.get_end_clipping() + 1
                    && graph.has_single_outgoing(a_i.get_path().back())
                    && graph.has_single_incoming(a_j.get_path()[0])) {
                // merge
                a_i.append(a_j, &graph);
                std::swap(a_i, a_j);
                a_i = Anchor();
            }
        }

        full_it = std::remove_if(full_it, anchors.end(), [](const auto &a) { return a.get_path().empty(); });
        common::logger->trace("Merged {} anchors down to {}", anchors.size(), full_it - anchors.begin());
        anchors.erase(full_it, anchors.end());
    }

    using AnchorIt = std::vector<Anchor>::iterator;
    std::sort(anchors.begin(), anchors.end(), AnchorLess<Anchor>());

    // std::cerr << "Anchors\n";
    // for (const auto &a : anchors) {
    //     std::cerr << "\t" << a << "\t" << a.get_path_spelling() << "\t" << fmt::format("{}", fmt::join(a.get_path(),",")) << "\n";
    // }

    std::vector<Alignment> alignments;

    if (anchors.empty())
        return alignments;

    ssize_t k = graph.get_k();

    assert(!dynamic_cast<const LabeledSeeder*>(this) || std::all_of(anchors.begin(), anchors.end(),
                                          [](const auto &a) { return a.get_label_class() != Anchor::nannot; }));

    // sdsl::bit_vector selected(anchors.size(), false);

    bool first_chain = true;
    DBGAlignerConfig::score_t first_chain_score = DBGAlignerConfig::ninf;
    size_t num_chains = 0;

    std::string dummy(query_.get_query().size(), '$');

    tsl::hopscotch_set<Anchor::label_class_t> found_labels;

    bool ext_success = false;
    Anchor::score_t last_chain_score = 0;
    common::logger->trace("Chaining");
    chain_anchors<AnchorIt>(query_, config_, anchors.begin(), anchors.end(),
        [this,k,sshash,match_score=config_.match_score("A"),&dummy](
                const Anchor &a_j,
                ssize_t max_dist,
                AnchorIt rbegin,
                AnchorIt rend,
                typename ChainScores<AnchorIt>::iterator rchain_scores,
                const ScoreUpdater<AnchorIt> &update_score) {
            auto rchain_scores_end = rchain_scores + (rend - rbegin);

            auto begin = std::make_reverse_iterator(rend);
            auto end = std::make_reverse_iterator(rbegin);
            auto chain_scores = std::make_reverse_iterator(rchain_scores_end);

            std::string_view query_j = a_j.get_seed();
            // const DBGAlignerConfig::score_t &score_j = std::get<0>(*(chain_scores + (end - begin)));
            const DBGAlignerConfig::score_t &score_j = std::get<0>(*(rchain_scores_end));
            // const DBGAlignerConfig::score_t dist_j = std::get<2>(*(rchain_scores_end));

            auto try_chain = [&](auto it, auto chain_scores) {
                const Anchor &a_i = *it;
                assert(a_i.get_label_class() == a_j.get_label_class());
                assert(a_i.get_orientation() == a_j.get_orientation());

                auto get_sshash_distance = [&](DeBruijnGraph::node_index node, DeBruijnGraph::node_index target_node) -> size_t {
                    DeBruijnGraph::seq_index contig = sshash->get_sequence_ids(node)[0];
                    DeBruijnGraph::seq_index target_contig = sshash->get_sequence_ids(target_node)[0];
                    if (contig == target_contig) {
                        if (node < sshash->reverse_complement(node)) {
                            // canonical strand
                            if (target_node >= node) {
#ifndef NDEBUG
                                if (target_node != node) {
                                    std::vector<DeBruijnGraph::node_index> path(target_node - node + 1);
                                    std::iota(path.begin(), path.end(), node);
                                    assert(path.back() == target_node);
                                    assert(spell_path(*sshash, path) != "");
                                }
#endif
                                return target_node - node;
                            }
                        } else if (node >= target_node && target_node > sshash->reverse_complement(target_node) && node > sshash->reverse_complement(node)) {
                            // rc strand
#ifndef NDEBUG
                            if (target_node != node) {
                                std::vector<DeBruijnGraph::node_index> path(node - target_node + 1);
                                std::iota(path.begin(), path.end(), target_node);
                                assert(path.back() == node);
                                assert(spell_path(*sshash, path) != "");
                            }
#endif
                            return node - target_node;
                        }
                    }

                    return std::numeric_limits<size_t>::max();
                };

                auto graph = make_aln_graph(a_i.get_label_class());

                std::string_view query_i = a_i.get_seed();

                DBGAlignerConfig::score_t dist = query_j.end() - query_i.end();

                if (dist < 0 || query_i.begin() >= query_j.begin())
                    return;

                DBGAlignerConfig::score_t base_score = std::get<0>(*chain_scores);

                // std::cerr << "chh\t" << a_i << "\t" << a_i.get_path_spelling() << " -> " << a_j << "\t" << a_j.get_path_spelling() << "\t" << dist << "\n";
                if (dist == 0 && a_j.get_end_trim()) {
                    assert(a_j.get_path().size() == 1);
                    if (base_score <= score_j)
                        return;

                    if (a_i.get_coord() != Anchor::ncoord && a_j.get_coord() != Anchor::ncoord) {
                        size_t a_i_pos = a_i.get_coord() + a_i.get_seed().size();
                        size_t a_j_pos = a_j.get_coord() + a_j.get_seed().size();
                        if (a_i_pos == a_j_pos)
                            update_score(base_score, it, 0);

                        return;
                    }

                    std::string_view a_j_trim_spelling = a_j.get_trim_spelling();
                    // std::cerr << "extleft" << a_i << "\t" << a_i.get_path_spelling() << " -> " << a_j << "\t" << a_j.get_path_spelling() << "\t" << a_j_trim_spelling.size() << "\n";
                    if (a_i.get_end_trim()) {
                        std::string_view a_i_trim_spelling = a_i.get_trim_spelling();
                        if (!std::equal(a_i_trim_spelling.begin(), a_i_trim_spelling.end(),
                                        a_j_trim_spelling.begin(), a_j_trim_spelling.begin() + a_i_trim_spelling.size()))
                            return;
                        a_j_trim_spelling = a_j_trim_spelling.substr(a_i_trim_spelling.size());
                    }
                    DeBruijnGraph::node_index node = a_i.get_path().back();
                    size_t traverse = sshash
                        ? get_sshash_distance(node, a_j.get_path()[0])
                        : std::numeric_limits<size_t>::max();
                    if (traverse == std::numeric_limits<size_t>::max() || traverse != a_j_trim_spelling.size()) {
                        traverse = 0;
                        graph.traverse(a_i.get_path().back(), a_j_trim_spelling,
                                    [&](DeBruijnGraph::node_index next, Anchor::label_class_t next_target) {
                                        assert(next_target == a_i.get_label_class());
                                        ++traverse;
                                        node = next;
                                    });
                        assert(traverse != a_j_trim_spelling.size() || node == a_j.get_path()[0]);
                    } else {
                        node = a_j.get_path()[0];
                    }
                    if (traverse == a_j_trim_spelling.size()) {
                        // std::cerr << "\tworked!\t" << base_score << "\n";
                        update_score(base_score, it, 0);
                    }
                }

                auto added_seq_start = std::max(query_j.begin(), query_i.end());
                std::string_view added_seq(added_seq_start, query_j.end() - added_seq_start);

                DBGAlignerConfig::score_t score = base_score + config_.match_score(added_seq)
                                                + (config_.right_end_bonus * !a_j.get_end_clipping());

                // if (score < score_j || (score == score_j && dist >= dist_j))
                if (score <= score_j)
                    return;

                if (a_i.get_coord() != Anchor::ncoord && a_j.get_coord() != Anchor::ncoord) {
                    size_t a_i_pos = a_i.get_coord() + a_i.get_seed().size();
                    size_t a_j_pos = a_j.get_coord() + a_j.get_seed().size();
                    if (a_j_pos > a_i_pos) {
                        size_t traversed = a_j_pos - a_i_pos;
                        size_t gap = std::abs(static_cast<ssize_t>(traversed) - static_cast<ssize_t>(dist));
                        if (gap) {
                            score -= (0.01 * config_.min_seed_length * gap + log2l(static_cast<double>(gap)) / 2) * match_score;
                        }

                        if (score > score_j)
                            update_score(score, it, traversed);
                    }

                    return;
                }
                // std::cerr << "\tchhs\n";

                size_t i_seed_size = k - a_i.get_end_trim();
                auto last_node_i_query = query_i.end() - i_seed_size;

                if (last_node_i_query >= query_j.begin()) {
                    // there should be overlapping nodes
                    // std::cerr << "olN" << a_i << " -> " << a_j << "\n";
                    size_t j_start = last_node_i_query - query_j.begin();
                    if (a_i.get_path().size() >= j_start
                            && a_j.get_path().size() >= j_start
                            && std::equal(a_j.get_path().begin(), a_j.get_path().begin() + j_start,
                                        a_i.get_path().end() - j_start, a_i.get_path().end())) {
                        // std::cerr << "\t\tworked! " << score << "\n";
                        update_score(score, it, dist);
                    }
                    return;
                } else if (query_i.end() > query_j.begin()) {
                    // overlapping nucleotides
                    // std::cerr << "oln" << a_i << " -> " << a_j << "\n";
                    assert(query_j.data() > query_i.data());
                    if (!a_i.get_end_trim() && !a_j.get_end_trim() && query_i.end() - query_j.begin() == k - 1) {
                        // no need to check
                        update_score(score, it, dist);
                    }
                    return;
                    // TODO: this check is expensive, skip it
                    // std::string_view left(query_i.data(), query_j.data() - query_i.data());
                    // DeBruijnGraph::node_index node = a_j.get_path()[0];
                    // ssize_t num_matches = 0;
                    // for (auto ait = left.rbegin(); ait != left.rend() && node != DeBruijnGraph::npos; ++ait) {
                    //     node = graph.traverse_back(node, *ait);
                    //     ++num_matches;
                    // }
                    // if (num_matches == static_cast<ssize_t>(left.size()) && node == a_i.get_path()[0]) {
                    //     assert(dist == num_matches - static_cast<ssize_t>(query_i.size()) + static_cast<ssize_t>(query_j.size()));
                    //     update_score(score, it, dist);
                    // }
                    // return;
                } else if (a_i.get_end_trim() && a_i.get_clipping() + a_i.get_seed().size() + a_i.get_end_trim() >= a_j.get_clipping()) {
                    // overlap in end parts
                    std::string_view a_i_trim = a_i.get_trim_spelling();

                    size_t olap = std::min(a_i.get_query().size(),
                                           a_i.get_clipping() + a_i.get_seed().size() + a_i.get_end_trim()) - a_j.get_clipping();
                    auto a_j_spelling = a_j.get_path_spelling().substr(0, k);
                    assert(a_i_trim.size() >= olap);
                    for (size_t del = 0; del < olap; ++del) {
                        assert(a_i_trim.size() >= olap - del);
                        // std::cerr << "olap o:" << olap - del << "\td: " << del << "\t" << a_i << " " << a_i.get_path_spelling() << " -> " << a_j << " " << a_j.get_path_spelling() << std::endl;

                        // if the character before a_j also matches, then skip this
                        // since we should work with an earlier seed instead
                        assert(a_j.get_clipping());
                        if (a_i_trim.size() >= olap - del + 1 && *(a_i_trim.end() - olap + del - 1) == *(a_j.get_seed().data() - 1))
                            continue;

                        if (std::equal(a_i_trim.end() - olap + del, a_i_trim.end(), a_j_spelling.begin(), a_j_spelling.begin() + olap - del)) {
                            // std::cerr << "\tfound\n";
                            DeBruijnGraph::node_index node = a_i.get_path().back();
                            size_t traversed = sshash
                                ? get_sshash_distance(node, a_j.get_path()[0])
                                : std::numeric_limits<size_t>::max();

                            if (traversed == std::numeric_limits<size_t>::max() || traversed != a_j_spelling.size() - olap + del) {
                                traversed = 0;
                                graph.traverse(node, a_j_spelling.substr(olap - del),
                                    [&](DeBruijnGraph::node_index next, Anchor::label_class_t next_target) {
                                        assert(next_target == a_i.get_label_class());
                                        node = next;
                                        ++traversed;
                                    }
                                );
                                assert(traversed != a_j_spelling.size() - olap + del || node == a_j.get_path()[0]);
                            } else {
                                node = a_j.get_path()[0];
                            }

                            if (traversed == a_j_spelling.size() - olap + del) {
                                traversed += a_i.get_path().size() - 1;
                                assert(traversed >= del);
                                assert(static_cast<ssize_t>(traversed - del) == query_j.begin() - query_i.begin());
                                DBGAlignerConfig::score_t cur_score = score;
                                size_t nmismatch = query_j.begin() - query_i.end();
                                std::string_view query_trim_window(query_i.data() + query_i.size(), nmismatch);
                                auto cur_mismatch = config_.score_sequences(query_trim_window,
                                                                            del == 0 ? std::string_view(a_i_trim.begin(), nmismatch)
                                                                                        : std::string_view(dummy.c_str(), nmismatch));
                                cur_score += cur_mismatch;
                                if (del)
                                    cur_score += config_.gap_opening_penalty + (del - 1) * config_.gap_extension_penalty;
                                if (cur_score > score_j) {
                                    // std::cerr << "\t\t\tworked! " << score << "\n";
                                    update_score(cur_score, it, traversed - query_i.size() + query_j.size());
                                    break;
                                }
                            }
                        }
                    }
                    return;
                } else if (!a_i.get_end_trim() && query_i.end() != query_j.begin()) {
                    // detect insertions
                    // std::cerr << "ins\tq: " << query_dist << "\td: " << traversed << "\tg: " << gap << "\t" << a_i << " " << a_i.get_path_spelling() << "\t" << a_j << " " << a_j.get_path_spelling() << "\n";

                    // don't try inserting if the previous character matches
                    if (!graph.traverse_back(a_j.get_path()[0], *(query_j.data() - 1)).first) {
                        DeBruijnGraph::node_index node = a_i.get_path().back();
                        size_t traversed = sshash
                            ? get_sshash_distance(node, a_j.get_path()[0])
                            : std::numeric_limits<size_t>::max();
                        if (traversed == std::numeric_limits<size_t>::max() || traversed != query_i.size() - a_i.get_path().size() + 1) {
                            traversed = 0;
                            auto a_j_spelling = a_j.get_path_spelling().substr(0, k);
                            graph.traverse(node, a_j_spelling,
                                [&](DeBruijnGraph::node_index next, Anchor::label_class_t next_target) {
                                    assert(next_target == a_i.get_label_class());
                                    node = next;
                                    ++traversed;
                                }
                            );
                            assert(traversed != query_i.size() - a_i.get_path().size() + 1 || node == a_j.get_path()[0]);
                        } else {
                            node = a_j.get_path()[0];
                        }
                        if (traversed == query_i.size() - a_i.get_path().size() + 1) {
                            size_t query_dist = a_j.get_clipping() - (a_i.get_clipping() + a_i.get_path().size() - 1);
                            assert(query_dist >= traversed);
                            size_t gap = query_dist - traversed;
                            if (gap > 0)
                                score += config_.gap_opening_penalty + (gap - 1) * config_.gap_extension_penalty;

                            if (score > score_j) {
                                // std::cerr << "\t\tworked! " << score << "\n";
                                update_score(score, it, traversed - (query_i.size() - a_i.get_path().size() + 1) + query_j.size());
                            }
                        }
                    }
                }

                // TODO: if a_i.get_end_trim() != 0, then we need to deal with
                //       a combination of D from the trim, X to the trim, then I of the query
                if (!sshash)
                    return;

                size_t traversed = get_sshash_distance(a_i.get_path().back(), a_j.get_path()[0]);
                if (traversed > 0 && traversed != std::numeric_limits<size_t>::max()) {
                    traversed += a_j.get_path().size() - 1 + a_i.get_end_trim();
                    size_t gap = std::abs(static_cast<ssize_t>(traversed) - static_cast<ssize_t>(dist));
                    if (gap) {
                        score -= (0.01 * config_.min_seed_length * gap + log2l(static_cast<double>(gap)) / 2) * match_score;
                    }

                    if (score > score_j)
                        update_score(score, it, traversed);
                }
            };

            // for (auto it = rbegin; it != rend; ++it, ++rchain_scores) {
            //     try_chain(it, rchain_scores);
            // }

            // std::cerr << "round 1\n";
            for (auto rit = begin; rit != end; ++rit, ++chain_scores) {
                try_chain(rit.base() - 1, chain_scores);
            }

            // std::cerr << "round 2\n";

        },
        [&](const AnchorChain<AnchorIt> &chain, const std::vector<DBGAlignerConfig::score_t> &score_traceback) {
            assert(chain.size());

            ext_success = false;

            bool ret_val = (!chain[0].first->get_clipping() && !chain[0].first->get_end_clipping())
                            || chain.size() > 1
                            || std::any_of(chain.begin(), chain.end(),
                                           [this](const auto &a) { return a.first->get_seed().size() > config_.min_seed_length
                                                                            || a.first->get_path().size() > 1; });

            // ret_val &= std::all_of(chain.begin(), chain.end(),
            //                        [&](const auto &a) { return !selected[a.first - anchors.begin()]
            //                                                     && (a.first->get_label_class() == Anchor::nannot
            //                                                             || !found_labels.count(a.first->get_label_class())); });

            ret_val &= std::all_of(chain.begin(), chain.end(),
                                   [&](const auto &a) { return a.first->get_label_class() == Anchor::nannot
                                                                        || !found_labels.count(a.first->get_label_class()); });

            // common::logger->info("Try Chain\t{}\t{}", score_traceback.back(), ret_val);
            if (!ret_val)
                return false;

            last_chain_score = score_traceback.back();
            ret_val = score_traceback.back() >= first_chain_score * config_.rel_score_cutoff;

            if (ret_val) {
                common::logger->trace("Chain\t{}\t{}", score_traceback.back(), ret_val);
                // for (const auto &[it, dist] : chain) {
                //     std::cerr << "\t" << *it << "\t" << dist << "\t" << it->get_label_class() << "\t" << it->get_path_spelling() << "\n";
                // }
            }

            // if (ret_val) {
            //     ++num_chains;
            //     for (const auto &[it, dist] : chain) {
            //         // selected[it - anchors.begin()] = true;
            //         if (it->get_label_class() != Anchor::nannot)
            //             found_labels.emplace(it->get_label_class());
            //     }
            // }

            return ret_val;
        },
        [this,align,match_score=config_.match_score("A")](
                AnchorIt last,
                AnchorIt next,
                Alignment&& aln,
                size_t next_to_last_dist,
                DBGAlignerConfig::score_t score_up_to_now,
                const AlignmentCallback &callback) {
            // std::cerr << "try\t" << *next << " -> " << aln << std::endl;
            size_t traversal_dist = next_to_last_dist - last->get_seed().size() + next->get_seed().size();
            const auto &graph = query_.get_graph();
            if (next->get_seed().begin() + next->get_path().size() == last->get_seed().begin() && traversal_dist == next->get_path().size()) {
                // simple extension
                // std::cerr << "\tSimple extension" << std::endl;
#ifndef NDEBUG
                size_t k = graph.get_k();
                if (graph.traverse(next->get_path().back(), aln.get_path_spelling()[k-1]) != aln.get_path()[0]) {
                    std::cerr << "\t" << next->get_path().back() << " > " << aln.get_path_spelling()[k-1] << "\t" << graph.traverse(next->get_path()[0], aln.get_path_spelling()[k-1]) << " vs. " << aln.get_path()[0] << std::endl;
                    assert(false);
                }
#endif
                // common::logger->info("Simple extension {} -> {}", *next, aln);
                std::vector<DeBruijnGraph::node_index> path;
                path.resize(aln.get_path().size() + next->get_path().size());
                std::copy(next->get_path().begin(), next->get_path().end(), path.begin());
                std::copy(aln.get_path().begin(), aln.get_path().end(), path.begin() + next->get_path().size());

                Cigar cigar(Cigar::CLIPPED, next->get_clipping());

                Cigar aln_cigar = aln.get_cigar();
                aln_cigar.trim_clipping();

                cigar.append(Cigar::MATCH, next->get_path().size());
                cigar.append(std::move(aln_cigar));

                Alignment next_aln(graph,
                                aln.get_query(),
                                aln.get_orientation(),
                                std::move(path),
                                std::move(cigar),
                                std::string(next->get_seed().substr(0, next->get_path().size())) + aln.get_path_spelling(),
                                aln.get_end_trim(),
                                next->get_label_class());
                // std::cerr << "\tnew aln: " << next_aln << std::endl;
                callback(std::move(next_aln));
                return;
            }

            if (!align) {
                // connect anchors without actually aligning or traversing
                size_t query_dist = last->get_clipping() - next->get_clipping();
                std::vector<DeBruijnGraph::node_index> path;
                path.resize(aln.get_path().size() + traversal_dist);
                std::copy(aln.get_path().rbegin(), aln.get_path().rend(), path.rbegin());

                Cigar cigar(Cigar::CLIPPED, next->get_clipping());
                std::string spelling;

                if (next->get_seed().end() > last->get_seed().begin()) {
                    // overlap
                    // std::cerr << "foo\t" << next->get_seed().end() - last->get_seed().begin() << "\t" << *next << " " << next->get_path_spelling() << " -> " << *last << " " << last->get_path_spelling() << std::endl;
                    assert(traversal_dist == query_dist);
                    // TODO: don't bother traversing, trust that the connection worked.
                    //       keep intermediate nodes as DeBruijnGraph::npos
#ifndef NDEBUG
                    DeBruijnGraph::node_index node = *(path.rbegin() + aln.get_path().size() - 1);
                    // auto jt = path.rbegin() + aln.get_path().size();
                    // for (auto it = next->get_seed().rend() - query_dist; it != next->get_seed().rend(); ++it, ++jt) {
                    for (auto it = next->get_seed().rend() - query_dist; it != next->get_seed().rend(); ++it) {
                        // *jt = graph.traverse_back(*(jt - 1), *it);
                        // assert(graph.in_graph(*jt));
                        node = graph.traverse_back(node, *it);
                        assert(graph.in_graph(node));
                    }
                    assert(node == next->get_path()[0]);
#endif
                    // assert(path[0] == next->get_path()[0]);
                    cigar.append(Cigar::MATCH, query_dist);
                    assert(next->get_query().begin() + next->get_clipping() + query_dist <= next->get_query().end());
                    spelling += std::string(next->get_seed().data(), next->get_seed().data() + query_dist);
                } else {
                    cigar.append(Cigar::MATCH, next->get_seed().size());
                    spelling += next->get_seed();
                    if (traversal_dist == query_dist) {
                        // no gap
                        assert(query_dist > next->get_seed().size());
                        cigar.append(Cigar::MISMATCH, query_dist - next->get_seed().size());
                        // assert(query_dist - next->get_seed().size() <= next->get_end_trim());
                        // spelling += next->get_trim_spelling().substr(0, query_dist - next->get_seed().size());
                        spelling += std::string(query_dist - next->get_seed().size(), '$');
                    } else if (traversal_dist > query_dist) {
                        // deletion
                        // TODO: which first, mismatch or deletion?
                        cigar.append(Cigar::DELETION, traversal_dist - query_dist);
                        spelling += std::string(traversal_dist - query_dist, '$');
                        if (query_dist > next->get_seed().size()) {
                            cigar.append(Cigar::MISMATCH, query_dist - next->get_seed().size());
                            spelling += std::string(query_dist - next->get_seed().size(), '$');
                        }
                    } else {
                        // insertion
                        // TODO: which first, mismatch or insert?
                        cigar.append(Cigar::INSERTION, query_dist - traversal_dist);
                        if (traversal_dist > next->get_seed().size()) {
                            cigar.append(Cigar::MISMATCH, traversal_dist - next->get_seed().size());
                            spelling += std::string(traversal_dist - next->get_seed().size(), '$');
                        }
                    }

                    if (next->get_end_trim()) {
                        assert(next->get_path().size() == 1);
                        std::string_view next_trim_spelling = next->get_trim_spelling();
                        auto end = next_trim_spelling.begin() + std::min(next_trim_spelling.size(), spelling.size() - next->get_seed().size());
                        assert(next->get_seed().size() + (end - next_trim_spelling.begin()) <= spelling.size());
                        std::copy(next_trim_spelling.begin(), end,
                                  spelling.begin() + next->get_seed().size());
                    }
                }

                std::copy(next->get_path().begin(), next->get_path().end(), path.begin());

                Cigar aln_cigar = aln.get_cigar();
                aln_cigar.trim_clipping();
                cigar.append(std::move(aln_cigar));

                Alignment next_aln(graph,
                                aln.get_query(),
                                aln.get_orientation(),
                                std::move(path),
                                std::move(cigar),
                                spelling + aln.get_path_spelling(),
                                aln.get_end_trim(),
                                next->get_label_class());
                callback(std::move(next_aln));
                return;
            }

            // std::cerr << "\talignment\n";
            DeBruijnGraph::node_index target_node = next->get_path().back();
            assert(traversal_dist > next->get_path().size() - 1);
            traversal_dist -= next->get_path().size() - 1;

            std::string_view query_window(
                next->get_seed().data() + next->get_path().size() - 1,
                aln.get_seed().begin() - next->get_seed().begin() - (next->get_path().size() - 1)
            );

            bool aln_found = false;
            auto *anno_buffer = make_aln_graph(last->get_label_class()).get_anno_buffer();
            align_bwd(
                query_.get_graph(),
                config_,
                aln.get_path()[0],
                last->get_seed().size(),
                query_window,
                traversal_dist,
                [&](auto&& spelling, auto&& bt_path, auto&& bt_cigar, size_t cost, Anchor::label_class_t final_target) { // callback
                    assert(bt_path[0] == next->get_path().back());
                    bt_path.insert(bt_path.begin(), next->get_path().begin(), next->get_path().end() - 1);

                    Cigar cigar(Cigar::CLIPPED, next->get_clipping());
                    cigar.append(Cigar::MATCH, next->get_path().size() - 1);
                    cigar.append(std::move(bt_cigar));
                    // std::cerr << "\textension: " << *next << "\t" << cigar.to_string() << "\t" << aln << std::endl;

                    Cigar aln_cigar = aln.get_cigar();
                    aln_cigar.trim_clipping();
                    cigar.append(std::move(aln_cigar));
                    bt_path.insert(bt_path.end(), aln.get_path().begin(), aln.get_path().end());

                    Alignment next_aln(query_.get_graph(),
                                    aln.get_query(),
                                    aln.get_orientation(),
                                    std::move(bt_path),
                                    std::move(cigar),
                                    std::string(next->get_spelling().substr(0, next->get_path().size() - 1)) + spelling + aln.get_path_spelling(),
                                    aln.get_end_trim(),
                                    final_target);
                    // std::cerr << "\t\t" << next_aln << std::endl;
                    aln_found = true;
                    callback(std::move(next_aln));
                },
                [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                    // start_backtracking
                    const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches, cur_target] = data;
                    // common::logger->info("Checkbt: n: {}, d: {}, qd: {} / {}, s: {}, nm: {}",
                    //                      query_.get_graph().get_node_sequence(node), dist, query_dist, query_window.size(),
                    //                      get_score(cost, query_dist, dist), num_matches);

                    if (aln_found || dist == 0 || query_dist == 0)
                        return false;

                    bool ret_val = query_dist == query_window.size() && node == target_node;
                    if (ret_val)
                        common::logger->trace("Final cost: {}", cost);

                    return ret_val;
                },
                [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                    // terminate branch
                    const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches, cur_target] = data;
                    // common::logger->info("Check: n: {}, d: {}, qd: {} / {}, s: {}",
                    //                      query_.get_graph().get_node_sequence(node), dist, query_dist, query_window.size(),
                    //                      get_score(cost, query_dist, dist));

                    return dist > traversal_dist || query_dist > query_window.size();

                    if (dist == 0 && query_dist == 0)
                        return false;

                    if (query_dist == query_window.size() || dist == traversal_dist) {
                        return true;
                    }

                    return false;
                },
                [&](size_t cost, const SMap &data, size_t query_dist, DeBruijnGraph::node_index node) {
                    // terminate
                    std::ignore = cost;
                    std::ignore = data;
                    std::ignore = query_dist;
                    // const auto &[dist, last_ext, last_node, last_op, mismatch_char, num_matches] = data;
                    // common::logger->info("Checkt: n: {}, d: {}, qd: {} / {}, s: {}",
                    //                  graph.get_node_sequence(node), dist, query_dist, query_window.size(),
                    //                  get_score(cost, query_dist, dist));
                    return aln_found || (query_dist == query_window.size() && node == target_node);
                },
                anno_buffer,
                last->get_label_class(),
                anno_buffer ? anno_buffer->get_label_path_coords({ aln.get_path()[0] },
                                                                 { last->get_label_class() })
                            : AnnotationBuffer::CoordinateSet()
            );
            if (!aln_found) {
                std::cerr << "q: " << aln.get_clipping() - next->get_clipping() << "\td: " << next_to_last_dist - last->get_seed().size() + next->get_seed().size() << "\t" << *next << " " << next->get_path_spelling() << " -> " << aln << std::endl;
                throw std::runtime_error("Failed to connect anchors");
            }
        },
        [this,&num_chains,&alignments,&ext_success,&last_chain_score,&first_chain,&first_chain_score,&found_labels](Alignment&& aln) {
            aln.trim_end();

            assert(aln.get_label_classes().size());
            if (auto anno_buffer = make_aln_graph(aln.get_label_classes()[0]).get_anno_buffer()) {
                if (anno_buffer->has_coordinates()) {
                    auto coords = anno_buffer->get_label_path_coords(aln.get_path(), aln.get_label_classes());
                    assert(coords.size());
                    aln.set_coords(std::move(coords));
                }
            }

            if (!ext_success) {
                ++num_chains;
                ext_success = true;
                if (first_chain) {
                    first_chain_score = last_chain_score;
                    first_chain = false;
                }
                for (auto label : aln.get_label_classes()) {
                    if (label != Anchor::nannot)
                        found_labels.emplace(label);
                }
            }

            if (common::get_verbose()) {
                common::logger->trace("Aln: {}\t{}\t{}\t{}",
                                      aln.get_orientation() ? "-" : "+",
                                      score_match(aln, config_), aln.get_cigar().to_string(), aln.get_label_classes()[0]);
            }
            // std::cerr << "\tInit aln\t" << aln << "\t" << aln.get_label_classes()[0] << std::endl;
            alignments.emplace_back(std::move(aln));
        }
    );

    // size_t num_selected = sdsl::util::cnt_one_bits(selected);
    // common::logger->info("Selected {} / {} anchors from {} chains. Produced {} alignments",
    //                         num_selected,
    //                         selected.size(),
    //                         num_chains, alignments.size());

    return alignments;
}

AlignmentGraph ExactSeeder::make_aln_graph(Anchor::label_class_t) const {
    return AlignmentGraph(query_.get_graph());
}

AlignmentGraph LabeledSeeder::make_aln_graph(Anchor::label_class_t target) const {
    return AlignmentGraph(query_.get_graph(), &anno_buffer_, target);
}

AlignmentGraph Extender::make_aln_graph(Anchor::label_class_t) const {
    return AlignmentGraph(query_.get_graph());
}

AlignmentGraph LabeledExtender::make_aln_graph(Anchor::label_class_t target) const {
    return AlignmentGraph(query_.get_graph(), &anno_buffer_, target);
}


} // namespace mtg::graph::align
