#include "aligner_extender_methods.hpp"

#include "dbg_aligner.hpp"

#include "common/utils/simd_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "common/logger.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/rc_dbg.hpp"
#include "kmer/kmer_extractor.hpp"


namespace mtg {
namespace graph {
namespace align {

using score_t = Alignment::score_t;
const score_t ninf = Alignment::ninf;
using kmer::KmerExtractorBOSS;

void SeedFilteringExtender
::extend_seed_end(const Alignment &seed,
                  const std::function<void(Alignment&&)> &callback,
                  bool force_fixed_seed,
                  score_t min_path_score) {
    if (!force_fixed_seed) {
        for (auto&& extension : get_extensions(seed, min_path_score, false)) {
            assert(extension.get_score() >= min_path_score);
            callback(std::move(extension));
        }

        return;
    }

    auto [left, next] = seed.split_seed(graph_->get_k() - 1, config_);
    auto extensions = get_extensions(next, min_path_score, true, 0, DeBruijnGraph::npos, false);

    bool called = false;
    for (auto&& extension : extensions) {
        if (extension.get_end_clipping() < seed.get_end_clipping()) {
            assert(extension.get_nodes().front() == next.get_nodes().front());
            auto aln = left;
            aln.splice(std::move(extension));
            if (aln.size() && aln.get_score() >= min_path_score) {
                assert(aln.get_clipping() == seed.get_clipping());
                if (aln.get_offset() > seed.get_offset())
                    aln.trim_offset(aln.get_offset() - seed.get_offset());

                assert(aln.get_offset() == seed.get_offset());
                callback(std::move(aln));
                called = true;
            }
        }
    }

    if (!called && seed.get_score() >= min_path_score)
        callback(Alignment(seed));
}

void SeedFilteringExtender
::rc_extend_rc(const Alignment &seed,
               const std::function<void(Alignment&&)> &callback,
               bool force_fixed_seed,
               score_t min_path_score) {
    // TODO: the workflow for force_fixed_seed is hard to do when coordinates are involved
    // TODO: I believe this is still guaranteed to get the same result
    force_fixed_seed &= seed.label_coordinates.empty();

    bool called = false;

#if _PROTEIN_GRAPH
    logger->warn("Front extension not supported for amino acid graphs");
#else
    auto rc_graph = std::shared_ptr<const DeBruijnGraph>(
        std::shared_ptr<const DeBruijnGraph>{}, graph_);
    if (graph_->get_mode() != DeBruijnGraph::CANONICAL)
        rc_graph = std::make_shared<RCDBG>(rc_graph);

    const DeBruijnGraph *old_graph = graph_;
    auto rev = seed;
    rev.reverse_complement(*rc_graph, get_query());
    if (rev.size() && rev.get_nodes().back()) {
        set_graph(*rc_graph);
        assert(rev.get_end_clipping());
        extend_seed_end(rev, [&](Alignment&& aln) {
            aln.reverse_complement(*rc_graph, seed.get_full_query_view());
            assert(aln.get_offset() == seed.get_offset());
            assert(aln.get_end_clipping() == seed.get_end_clipping());
            if (aln.size()) {
                assert(aln.get_score() >= min_path_score);
                callback(std::move(aln));
                called = true;
            }
        }, force_fixed_seed, min_path_score);
        set_graph(*old_graph);
    }
#endif

    if (!called && seed.get_score() >= min_path_score)
        callback(Alignment(seed));
}

std::vector<Alignment> SeedFilteringExtender::connect_seeds(const Alignment &first,
                                                            const Alignment &second,
                                                            int64_t coord_dist,
                                                            score_t min_path_score) {
    auto [left, next] = first.split_seed(graph_->get_k() - 1, config_);
    DEBUG_LOG("Split seed: {}\n\t{}\n\t{}", first, left, next);
    coord_dist += second.get_sequence().size() + next.get_sequence().size()
            - first.get_sequence().size();
    assert(coord_dist > 0);

    auto extensions = get_extensions(next, min_path_score,
        true,                                // force_fixed_seed
        coord_dist,                          // target_length
        second.get_nodes().back(),           // target_node
        false,                               // trim_offset_after_extend
        second.get_end_clipping(),           // trim_query_suffix
        first.get_score() - next.get_score() // added_xdrop
    );

    std::vector<Alignment> alignments;

    for (auto&& extension : extensions) {
        if (extension.get_end_clipping() == second.get_end_clipping()
                && extension.get_nodes().back() == second.get_nodes().back()) {
            auto alignment = left;
            assert(extension.get_nodes().front() == next.get_nodes().front());
            alignment.splice(std::move(extension));
            assert(alignment.is_valid(*graph_, &config_));
            alignments.emplace_back(std::move(alignment));
        }
    }

    return alignments;
}

DefaultColumnExtender::DefaultColumnExtender(const DeBruijnGraph &graph,
                                             const DBGAlignerConfig &config,
                                             std::string_view query)
      : SeedFilteringExtender(&graph, config, query), query_(query) {
    // compute exact-match scores for all suffixes of the query
    partial_sums_.reserve(query_.size() + 1);
    partial_sums_.resize(query_.size(), 0);
    std::transform(query_.begin(), query_.end(), partial_sums_.begin(),
                   [&](char c) { return config_.score_matrix[c][c]; });

    std::partial_sum(partial_sums_.rbegin(), partial_sums_.rend(), partial_sums_.rbegin());
    assert(query_.empty() || config_.match_score(query_) == partial_sums_.front());
    assert(query_.empty() || config_.score_matrix[query_.back()][query_.back()] == partial_sums_.back());

    partial_sums_.push_back(0);

    profile_score_.resize(KmerExtractorBOSS::alphabet.size() + 1);
    profile_op_.resize(KmerExtractorBOSS::alphabet.size() + 1);

    for (size_t i = 0; i < profile_score_.size(); ++i) {
        // precompute profiles to store match/mismatch scores and Cigar::Operators
        // in contiguous arrays
        profile_score_[i] = AlignedVector<score_t>(query_.size() + kPadding);
        profile_op_[i] = AlignedVector<Cigar::Operator>(query_.size() + kPadding);
        char c = i != KmerExtractorBOSS::alphabet.size()
            ? KmerExtractorBOSS::decode(i)
            : '\0';
        const auto &row = config_.score_matrix[c];
        const auto &op_row = kCharToOp[c];

        // the first cell in a DP table row is one position before the first character,
        // so we need to shift the indices of profile_score_ and profile_op_
        std::transform(query_.begin(), query_.end(), profile_score_[i].begin() + 1,
                       [&row](char q) { return row[q]; });

        std::transform(query_.begin(), query_.end(), profile_op_[i].begin() + 1,
                       [&op_row](char q) { return op_row[q]; });
    }
}

DefaultColumnExtender::DefaultColumnExtender(const IDBGAligner &aligner,
                                             std::string_view query)
      : DefaultColumnExtender(aligner.get_graph(), aligner.get_config(), query) {}

bool SeedFilteringExtender::check_seed(const Alignment &seed) const {
    if (seed.empty())
        return false;

    assert(seed.get_nodes().back());
    assert(seed.get_cigar().size());

    node_index node = seed.get_nodes().back();
    if (dynamic_cast<const RCDBG*>(graph_))
        node += graph_->max_index();

    auto it = conv_checker_.find(node);

    if (it == conv_checker_.end())
        return true;

    size_t pos = seed.get_query_view().size() + seed.get_clipping() - 1;
    const auto &[start, vec] = it->second;

    return pos < start
            || pos - start >= vec.size()
            || vec[pos - start] < seed.get_score();
}

bool SeedFilteringExtender::set_seed(const Alignment &seed) {
    assert(seed.get_query_view().size() + seed.get_clipping() + seed.get_end_clipping()
            == query_size_);
    DEBUG_LOG("Seed: {}", seed);
    assert(seed.is_valid(*graph_, &config_));
    seed_ = &seed;
    clear_conv_checker();
    return seed_;
}

score_t SeedFilteringExtender::update_seed_filter(node_index node,
                                                  size_t query_start,
                                                  const score_t *s_begin,
                                                  const score_t *s_end) {
    if (node == DeBruijnGraph::npos)
        return *std::max_element(s_begin, s_end);

    if (dynamic_cast<const RCDBG*>(graph_))
        node += graph_->max_index();

    assert(s_end >= s_begin);
    assert(query_start + (s_end - s_begin) <= query_size_);

    size_t size = s_end - s_begin;

    auto it = conv_checker_.find(node);

    if (it == conv_checker_.end()) {
        conv_checker_.emplace(node, ScoreVec(query_start, { s_begin, s_end }));
        return *std::max_element(s_begin, s_end);
    }

    auto &[start, vec] = it.value();
    if (query_start + size <= start) {
        vec.insert(vec.begin(), start - query_start, ninf);
        std::copy(s_begin, s_end, vec.begin());
        start = query_start;
        return *std::max_element(s_begin, s_end);
    }

    if (query_start >= start + vec.size()) {
        vec.reserve(query_start + size - start);
        vec.insert(vec.end(), query_start - start - vec.size(), ninf);
        vec.insert(vec.end(), s_begin, s_end);
        return *std::max_element(s_begin, s_end);
    }

    // overlap
    if (query_start < start) {
        vec.insert(vec.begin(), start - query_start, ninf);
        start = query_start;
    }

    if (query_start + size > start + vec.size())
        vec.resize(query_start + size - start, ninf);

    score_t max_changed_value = ninf;
    score_t *v = vec.data() + query_start - start;
    for (size_t j = 0; j < size; ++j) {
        if (s_begin[j] > v[j] * config_.rel_score_cutoff) {
            v[j] = std::max(v[j], s_begin[j]);
            max_changed_value = std::max(max_changed_value, v[j]);
        }
    }

    return max_changed_value;
}

bool SeedFilteringExtender
::filter_nodes(node_index node, size_t query_start, size_t query_end) {
    assert(query_end >= query_start);
    assert(query_end <= query_size_);
    constexpr score_t mscore = -ninf;
    size_t size = query_end - query_start;

    auto it = conv_checker_.find(node);
    if (it == conv_checker_.end()) {
        conv_checker_.emplace(
            node, ScoreVec(query_start, AlignedVector<score_t>(size, mscore))
        );
        return true;
    }

    auto &[start, vec] = it.value();
    if (query_start + size <= start) {
        vec.insert(vec.begin(), start - query_start, ninf);
        std::fill(vec.begin(), vec.begin() + size, mscore);
        start = query_start;
        return true;
    }

    if (query_start >= start + vec.size()) {
        vec.reserve(query_start + size - start);
        vec.insert(vec.end(), query_start - start - vec.size(), ninf);
        vec.insert(vec.end(), size, mscore);
        return true;
    }

    // overlap
    if (query_start < start) {
        vec.insert(vec.begin(), start - query_start, ninf);
        start = query_start;
    }

    if (query_start + size > start + vec.size())
        vec.resize(query_start + size - start, ninf);

    bool converged = true;
    score_t *v = vec.data() + query_start - start;
    for (size_t j = 0; j < size; ++j) {
        if (mscore > v[j]) {
            converged = false;
            v[j] = mscore;
        }
    }

    return !converged;
}

void update_column(size_t prev_end,
                   const score_t *S_prev_v,
                   const score_t *F_prev_v,
                   AlignedVector<score_t> &S_v,
                   AlignedVector<score_t> &E_v,
                   AlignedVector<score_t> &F_v,
                   const score_t *profile_scores,
                   score_t xdrop_cutoff,
                   const DBGAlignerConfig &config_,
                   score_t init_score,
                   size_t offset) {
    static_assert(DefaultColumnExtender::kPadding == 5);
    constexpr size_t width = DefaultColumnExtender::kPadding - 1;
    const simde__m128i gap_open = simde_mm_set1_epi32(config_.gap_opening_penalty);
    const simde__m128i gap_extend = simde_mm_set1_epi32(config_.gap_extension_penalty);
    const simde__m128i xdrop_v = simde_mm_set1_epi32(xdrop_cutoff - 1);
    const simde__m128i ninf_v = simde_mm_set1_epi32(ninf);
    const simde__m128i score_v = simde_mm_set1_epi32(init_score);
    for (size_t j = 0; j < prev_end; j += width) {
        // ensure that nothing will access out of bounds
        assert(j + DefaultColumnExtender::kPadding <= S_v.capacity());

        // match = j ? S_prev_v[j - 1] + profile_scores[j] : ninf;
        simde__m128i match;
        if (j) {
            match = simde_mm_add_epi32(simde_mm_loadu_si128((simde__m128i*)&S_prev_v[j - 1]),
                                       simde_mm_loadu_si128((simde__m128i*)&profile_scores[j]));
            match = simde_mm_add_epi32(match, score_v);
        } else {
            // rotate elements to the right, then insert ninf in first cell
            match = simde_mm_shuffle_epi32(simde_mm_loadu_si128((simde__m128i*)&S_prev_v[j]), 0b10010000);
            match = simde_mm_add_epi32(match, simde_mm_loadu_si128((simde__m128i*)&profile_scores[j]));
            match = simde_mm_add_epi32(match, score_v);
            match = simde_mm_insert_epi32(match, ninf, 0);
        }

        // del_score = std::max(del_open, del_extend);
        simde__m128i del_score;
        if (offset > 1) {
            del_score = simde_mm_max_epi32(
                simde_mm_add_epi32(simde_mm_loadu_si128((simde__m128i*)&S_prev_v[j]), gap_open),
                simde_mm_add_epi32(simde_mm_loadu_si128((simde__m128i*)&F_prev_v[j]), gap_extend)
            );
            del_score = simde_mm_add_epi32(del_score, score_v);
        } else {
            del_score = ninf_v;
        }

        // F_v[j] = del_score
        simde_mm_store_si128((simde__m128i*)&F_v[j], del_score);

        // match = max(match, del_score)
        match = simde_mm_max_epi32(match, del_score);

        // E_v[j + 1] = S[j] + gap_open
        simde__m128i ins_open_next = simde_mm_add_epi32(match, gap_open);
        simde_mm_storeu_si128((simde__m128i*)&E_v[j + 1], ins_open_next);

        // This rolling update is hard to vectorize
        // E_v[j + 1] = max(E_v[j + 1], E_v[j] + gap_extend)
        E_v[j + 1] = std::max(E_v[j] + config_.gap_extension_penalty, E_v[j + 1]);
        E_v[j + 2] = std::max(E_v[j + 1] + config_.gap_extension_penalty, E_v[j + 2]);
        E_v[j + 3] = std::max(E_v[j + 2] + config_.gap_extension_penalty, E_v[j + 3]);
        E_v[j + 4] = std::max(E_v[j + 3] + config_.gap_extension_penalty, E_v[j + 4]);

        // S_v[j] = max(match, E_v[j])
        match = simde_mm_max_epi32(match, simde_mm_load_si128((simde__m128i*)&E_v[j]));

        // match >= xdrop_cutoff
        simde__m128i mask = simde_mm_cmpgt_epi32(match, xdrop_v);
        match = simde_mm_blendv_epi8(ninf_v, match, mask);

        simde_mm_store_si128((simde__m128i*)&S_v[j], match);
    }

    if (S_v.size() > std::max(size_t{1}, prev_end)) {
        size_t j = S_v.size() - 1;
        score_t match = std::max(S_prev_v[j - 1] + init_score + profile_scores[j], E_v[j]);
        if (match >= xdrop_cutoff)
            S_v[j] = match;
    }
}

// add insertions to the end of the array until the score drops too low
void extend_ins_end(AlignedVector<score_t> &S,
                    AlignedVector<score_t> &E,
                    AlignedVector<score_t> &F,
                    size_t max_size,
                    score_t xdrop_cutoff,
                    const DBGAlignerConfig &config_) {
    if (S.size() < max_size) {
        score_t ins_score = std::max(
            S.back() + config_.gap_opening_penalty,
            E.back() + config_.gap_extension_penalty
        );

        if (ins_score >= xdrop_cutoff) {
            S.push_back(ins_score);
            E.push_back(ins_score);
            F.push_back(ninf);

            while (E.back() + config_.gap_extension_penalty >= xdrop_cutoff
                    && E.size() < max_size) {
                E.push_back(E.back() + config_.gap_extension_penalty);
                S.push_back(E.back());
                F.push_back(ninf);
            }

            // allocate and initialize enough space to allow the SIMD code to access these
            // vectors in 16 byte blocks without reading out of bounds
            S.reserve(S.size() + DefaultColumnExtender::kPadding);
            E.reserve(E.size() + DefaultColumnExtender::kPadding);
            F.reserve(F.size() + DefaultColumnExtender::kPadding);

            std::fill(S.data() + S.size(), S.data() + S.capacity(), ninf);
            std::fill(E.data() + E.size(), E.data() + E.capacity(), ninf);
            std::fill(F.data() + F.size(), F.data() + F.capacity(), ninf);
        }
    }
}

void DefaultColumnExtender
::call_outgoing(node_index node,
                size_t /* max_prefetch_distance */,
                const std::function<void(node_index, char, score_t)> &callback,
                size_t table_i,
                bool force_fixed_seed) {
    assert(node == table[table_i].node);
    size_t next_offset = table[table_i].offset + 1;
    size_t seed_pos = next_offset - this->seed_->get_offset();
    bool in_seed = seed_pos < this->seed_->get_sequence().size();

    // Get the next node(s) from the graph. If the current node is
    // part of the seed and the arguments tell us to force using
    // the full seed, then pick the next node from the seed.
    if (in_seed && next_offset < graph_->get_k()) {
        assert(this->seed_->get_nodes().front());
        callback(this->seed_->get_nodes().front(),
                 this->seed_->get_sequence()[seed_pos],
                 0);
    } else if (in_seed && force_fixed_seed) {
        size_t node_i = next_offset - graph_->get_k() + 1;
        assert(node_i);
        assert(node_i < this->seed_->get_nodes().size());
        assert(node == this->seed_->get_nodes()[node_i - 1]);
        node_index next_node = this->seed_->get_nodes()[node_i];
        char next_c = this->seed_->get_sequence()[seed_pos];
        callback(next_node, next_c,
            // (next_node ? 0 : (!node ? config_.gap_extension_penalty : config_.gap_opening_penalty)) +
                (node_i - 1 < this->seed_->extra_scores.size() ? this->seed_->extra_scores[node_i - 1] : 0)

        );
        assert(!node || next_c == boss::BOSS::kSentinel || !next_node ||
                graph_->traverse(node, next_c) == next_node);
    } else {
        assert(node);
        graph_->call_outgoing_kmers(node, [&](node_index next, char c) {
            if (c != boss::BOSS::kSentinel)
                callback(next, c, 0);
        });
    }
}

template <typename... RestArgs>
DefaultColumnExtender::DPTColumn
DefaultColumnExtender::DPTColumn::create(size_t size, RestArgs&&... args) {
    DPTColumn column { {}, {}, {}, std::forward<RestArgs>(args)... };
    // allocate and initialize enough space to allow the SIMD code to access these
    // vectors in 16 byte blocks without reading out of bounds
    column.S.reserve(size + kPadding);
    column.E.reserve(size + kPadding);
    column.F.reserve(size + kPadding);

    // the size is set properly to allow for AlignedVector methods (size(), push_back())
    // to function properly
    column.S.resize(size, ninf);
    column.E.resize(size, ninf);
    column.F.resize(size, ninf);

    std::fill(column.S.data() + column.S.size(), column.S.data() + column.S.capacity(), ninf);
    std::fill(column.E.data() + column.E.size(), column.E.data() + column.E.capacity(), ninf);
    std::fill(column.F.data() + column.F.size(), column.F.data() + column.F.capacity(), ninf);

    return column;
}

std::vector<Alignment> DefaultColumnExtender::extend(score_t min_path_score,
                                                     bool force_fixed_seed,
                                                     size_t target_length,
                                                     node_index target_node,
                                                     bool trim_offset_after_extend,
                                                     size_t trim_query_suffix,
                                                     score_t added_xdrop) {
    assert(this->seed_);

    ++num_extensions_;

    min_path_score = std::max(0, min_path_score);

    explored_nodes_previous_ += table.size();
    table.clear();
    prev_starts.clear();

    // saturating addition
    score_t xdrop = std::min(config_.xdrop, std::numeric_limits<score_t>::max() - added_xdrop)
                        + added_xdrop;
    assert(xdrop > 0);

    xdrop_cutoffs_.assign(1, std::make_pair(0u, std::max(-xdrop, ninf + 1)));
    assert(xdrop_cutoffs_[0].second < 0);

    if (!config_.global_xdrop)
        scores_reached_.assign(1, 0);

    typedef typename decltype(xdrop_cutoffs_)::value_type xdrop_cutoff_v;
    table_size_bytes_ = sizeof(table) + sizeof(xdrop_cutoffs_) + sizeof(xdrop_cutoff_v)
        + sizeof(scores_reached_) + sizeof(score_t);

    size_t start = this->seed_->get_clipping();

    // the sequence to align (the suffix of the query starting from the seed)
    std::string_view window(this->seed_->get_query_view().data(),
                            query_.data() + query_.size()
                                - this->seed_->get_query_view().data() - trim_query_suffix);
    score_t partial_sum_offset = partial_sums_.at(start + window.size());
    assert(partial_sums_.at(start) - partial_sum_offset == config_.match_score(window));

    ssize_t seed_offset = static_cast<ssize_t>(this->seed_->get_offset()) - 1;

    // initialize the root of the tree
    table.emplace_back(DPTColumn::create(
            1, this->seed_->get_nodes().front(), static_cast<size_t>(-1),
            '\0', seed_offset, 0, 0, 0u, 0));

    {
        auto &[S, E, F, node, i_prev, c, offset, max_pos, trim,
               xdrop_cutoff_i, score] = table[0];
        S[0] = config_.left_end_bonus && !seed_->get_clipping() ? config_.left_end_bonus : 0;
        extend_ins_end(S, E, F, window.size() + 1 - trim,
                       xdrop_cutoffs_[xdrop_cutoff_i].second, config_);

        table_size_bytes_ = sizeof(decltype(table)::value_type) * table.capacity()
                                + S.capacity() * sizeof(decltype(S)::value_type)
                                + E.capacity() * sizeof(decltype(E)::value_type)
                                + F.capacity() * sizeof(decltype(F)::value_type);
    }

    // The nodes in the traversal (with corresponding score columns) are sorted by
    // 1) their score (higher is better), then by
    // 2) the absolute distance of their highest scoring index from the score
    //    matrix diagonal (lower is better), finally by
    // 3) Their index in the table vector (higher is better, for better cache locality)
    using TableIt = std::tuple<score_t, /* max changed value */
                               ssize_t, /* negative off_diag */
                               size_t, /* table idx */
                               score_t /* max score */>;
    min_cell_score_ = 0;
    score_t best_score = 0;

    std::priority_queue<TableIt> queue;

    // Initialize the node traversal heap with the root.
    queue.emplace(0, 0, 0, 0);

    std::vector<size_t> tips;

    while (queue.size()) {
        size_t i = std::get<2>(queue.top());
        queue.pop();
        bool has_single_outgoing = true;

        while (has_single_outgoing) {
            has_single_outgoing = false;
            std::vector<std::tuple<node_index, char, score_t>> outgoing;
            size_t next_offset = table[i].offset + 1;

            ssize_t begin = 0;
            ssize_t prev_end = window.size() + 1;
            size_t prev_xdrop_cutoff_i = table[i].xdrop_cutoff_i;
            score_t prev_xdrop_cutoff = xdrop_cutoffs_[prev_xdrop_cutoff_i].second;

            bool in_seed = next_offset - this->seed_->get_offset()
                                < this->seed_->get_sequence().size();

            {
                const auto &[S, E, F, node, i_prev, c, offset, max_pos, trim,
                             xdrop_cutoff_i, score] = table[i];

                double node_counter = config_.global_xdrop
                    ? table.size()
                    : next_offset - seed_offset;

                // if we are currently not on the optimal path, check early cutoff criteria
                if (S[max_pos - trim] < best_score) {
                    // if too many nodes have been explored, give up
                    if (node_counter / window.size() >= config_.max_nodes_per_seq_char) {
                        DEBUG_LOG("Position {}: Alignment node limit reached, stopping this branch",
                                  next_offset - seed_->get_offset());
                        if (config_.global_xdrop)
                            queue = std::priority_queue<TableIt>();

                        continue;
                    }

                    if (static_cast<double>(table_size_bytes_) / 1'000'000
                            > config_.max_ram_per_alignment) {
                        DEBUG_LOG("Position {}: Alignment RAM limit reached, stopping extension",
                                  next_offset - seed_->get_offset());
                        queue = std::priority_queue<TableIt>();
                        continue;
                    }
                }

                // determine maximal range within the xdrop score cutoff
                auto in_range = [prev_xdrop_cutoff](score_t s) {
                    return s >= prev_xdrop_cutoff;
                };

                begin = std::find_if(S.begin(), S.end(), in_range) - S.begin() + trim;
                prev_end = std::find_if(S.rbegin(), S.rend(), in_range).base() - S.begin() + trim;

                if (prev_end <= begin) {
                    DEBUG_LOG("Position {}: Bandwidth reached 0", next_offset - seed_->get_offset());
                    continue;
                }

                call_outgoing(node, window.size() + 1 - offset - S.size(),
                              [&](node_index next, char c, score_t s) {
                                  c = toupper(c);
                                  outgoing.emplace_back(next, c, s);
                              }, i, force_fixed_seed);

                if (outgoing.empty()) {
                    DEBUG_LOG("Positon {}: Tip reached", next_offset - seed_->get_offset());
                    tips.push_back(i);
                    continue;
                }
            }

            size_t end = std::min(static_cast<size_t>(prev_end), window.size()) + 1;

            for (const auto &[next, c, score] : outgoing) {
                bool forked = outgoing.size() > 1;
                bool forked_xdrop = !config_.global_xdrop && forked;
                size_t xdrop_cutoffs_sizediff = xdrop_cutoffs_.capacity();
                if (forked_xdrop) {
                    xdrop_cutoffs_.emplace_back(table.size(), prev_xdrop_cutoff);
                    xdrop_cutoffs_sizediff = xdrop_cutoffs_.capacity() - xdrop_cutoffs_sizediff;
                } else {
                    xdrop_cutoffs_sizediff = 0;
                }

                size_t table_sizediff = table.capacity();
                table.emplace_back(DPTColumn::create(end - begin, next, i, c,
                    static_cast<ssize_t>(next_offset),
                    begin, begin,
                    forked_xdrop ? xdrop_cutoffs_.size() - 1 : prev_xdrop_cutoff_i,
                    score
                ));

                const auto &[S_prev, E_prev, F_prev, node_prev, i_prev, c_prev,
                             offset_prev, max_pos_prev, trim_prev,
                             xdrop_cutoff_i_prev, score_prev] = table[i];

                auto &[S, E, F, node_cur, i_cur, c_stored, offset, max_pos, trim,
                       xdrop_cutoff_i, score_cur] = table.back();

                score_t &xdrop_cutoff = xdrop_cutoffs_[xdrop_cutoff_i].second;

                assert(i_cur == i);
                assert(node_cur == next);
                assert(c_stored == c);
                assert(offset == offset_prev + 1);
                assert(!node_cur || c == graph_->get_node_sequence(node_cur)[
                    std::min(static_cast<ssize_t>(graph_->get_k()) - 1, offset)]);

                // compute column scores
                update_column(prev_end - trim,
                              S_prev.data() + trim - trim_prev,
                              F_prev.data() + trim - trim_prev,
                              S, E, F,
                              profile_score_[KmerExtractorBOSS::encode(c)].data() + start + trim,
                              xdrop_cutoff, config_, score, offset);

                extend_ins_end(S, E, F, window.size() + 1 - trim, xdrop_cutoff, config_);

                assert(max_pos >= trim);
                assert(static_cast<size_t>(max_pos - trim) < S.size());

                // find the maximal scoring position which is closest to the diagonal
                // TODO: this can be done with SIMD, but it's not a bottleneck
                ssize_t cur_offset = begin;
                ssize_t diag_i = offset - seed_offset;
                bool has_extension = in_seed;
                const score_t *partial_sums = &partial_sums_[start + trim];
                score_t extension_cutoff
                    = best_score * config_.rel_score_cutoff + partial_sum_offset;
                score_t max_diff = ninf;

                size_t scores_reached_sizediff = 0;
                bool scores_reached_cutoff = true;
                if (!config_.global_xdrop) {
                    scores_reached_sizediff = scores_reached_.capacity();
                    scores_reached_.resize(S.size() + trim + 1, ninf);
                    scores_reached_sizediff = scores_reached_.capacity() - scores_reached_sizediff;
                }

                for (size_t j = 0; j < S.size(); ++j, ++cur_offset) {
                    if (S[j] != ninf)
                        min_cell_score_ = std::min(min_cell_score_, S[j]);

                    if (std::make_pair(S[j], std::abs(max_pos - diag_i))
                            > std::make_pair(S[max_pos - begin], std::abs(cur_offset - diag_i))) {
                        max_pos = j + begin;
                    }

                    if (!config_.global_xdrop) {
                        scores_reached_[trim + j] = std::max(scores_reached_[trim + j], S[j]);
                        scores_reached_cutoff = (S[j] >= scores_reached_[trim + j] * config_.rel_score_cutoff);
                    }

                    // check if this node can be extended to get a better alignment
                    assert(partial_sums[j] - partial_sum_offset
                            == config_.match_score(window.substr(j + trim)));
                    if (!has_extension && scores_reached_cutoff
                            && S[j] + partial_sums[j] >= extension_cutoff) {
                        has_extension = true;
                    }

                    if (static_cast<size_t>(trim - trim_prev) < S_prev.size()
                            && S[j] - S_prev[j + trim - trim_prev] > max_diff) {
                        max_diff = S[j] - S_prev[j + trim - trim_prev];
                    }
                }

                assert(max_pos >= trim);
                assert(static_cast<size_t>(max_pos - trim) < S.size());

                score_t max_val = S[max_pos - trim];

                if (static_cast<size_t>(offset - seed_offset) < target_length + 1) {
                    has_extension = true;
                } else if (target_length) {
                    DEBUG_LOG("Target distance reached");
                    has_extension = false;
                }

                if (!in_seed && max_val < xdrop_cutoff) {
                    DEBUG_LOG("Position {}: x-drop: {} < {}",
                              next_offset - seed_->get_offset(), max_val, xdrop_cutoff);
                    pop(table.size() - 1);
                    if (forked_xdrop)
                        xdrop_cutoffs_.pop_back();

                    continue;
                }

                if (!in_seed && !has_extension) {
                    DEBUG_LOG("Position {}: no extension possible from score {}. "
                              "Best score so far is {}",
                              next_offset - seed_->get_offset(), max_val,
                              best_score);
                    pop(table.size() - 1);
                    if (forked_xdrop)
                        xdrop_cutoffs_.pop_back();

                    continue;
                }

                table_sizediff = table.capacity() - table_sizediff;

                table_size_bytes_ += sizeof(decltype(table)::value_type) * table_sizediff
                    + S.capacity() * sizeof(decltype(S)::value_type)
                    + E.capacity() * sizeof(decltype(E)::value_type)
                    + F.capacity() * sizeof(decltype(F)::value_type)
                    + sizeof(decltype(scores_reached_)::value_type) * scores_reached_sizediff
                    + sizeof(xdrop_cutoff_v) * xdrop_cutoffs_sizediff;

                if (max_val - xdrop_cutoff > xdrop)
                    xdrop_cutoff = max_val - xdrop;

                // if the best score in this column is above the xdrop score
                // then check if the extension can continue
                best_score = std::max(best_score, max_val);

                // skip the first index since it corresponds to the position
                // before the query start
                size_t vec_offset = start + begin - static_cast<bool>(begin);
                score_t *s_begin = S.data() + !begin;
                score_t *s_end = S.data() + S.size();

                assert(s_begin <= s_end);
                assert(vec_offset + (s_end - s_begin) <= query_.size());

                // if this node has not been reached by a different
                // alignment with a better score, continue
                score_t converged_score = update_seed_filter(
                    next, vec_offset, s_begin, s_end
                );
                if (converged_score != ninf) {
                    // if this next node is the only next option, take it without pushing
                    // to the queue
                    if (outgoing.size() == 1) {
                        has_single_outgoing = true;
                        i = table.size() - 1;
                    } else {
                        queue.emplace(converged_score, -std::abs(max_pos - diag_i),
                                      table.size() - 1, max_val);
                    }
                } else {
                    DEBUG_LOG("Position {}: Dropped due to convergence",
                              next_offset - seed_->get_offset());
                }
            }
        }
    }

    if (config_.no_backtrack)
        return { *seed_ };

    std::sort(tips.begin(), tips.end());

    auto extensions = backtrack(min_path_score, window,
                                trim_query_suffix ? 0 : config_.right_end_bonus,
                                tips, target_node);
    if (trim_offset_after_extend) {
        for (auto &extension : extensions) {
            assert(extension.is_valid(*graph_, &config_));
            extension.trim_offset();
            assert(extension.is_valid(*graph_, &config_));
        }
    }

    return extensions;
}

Alignment DefaultColumnExtender::construct_alignment(Cigar cigar,
                                                     size_t clipping,
                                                     std::string_view window,
                                                     std::vector<node_index> final_path,
                                                     std::string match,
                                                     score_t score,
                                                     size_t offset,
                                                     const std::vector<score_t> &score_trace,
                                                     score_t extra_score) const {
    assert(final_path.size());
    assert(cigar.size());
    cigar.append(Cigar::CLIPPED, clipping);

    std::reverse(cigar.data().begin(), cigar.data().end());
    std::reverse(final_path.begin(), final_path.end());
    std::reverse(match.begin(), match.end());

    Alignment extension(window, std::move(final_path), std::move(match), score,
                        std::move(cigar), 0, this->seed_->get_orientation(), offset);
    extension.extend_query_begin(query_.data());
    extension.extend_query_end(query_.data() + query_.size());
    extension.extra_score = extra_score;
    if (extra_score) {
        auto score_it = score_trace.rend() - extension.get_nodes().size() + 1;
        assert(!*(score_it - 1));
        extension.extra_scores = std::vector<score_t>(score_it, score_trace.rend());
    }
    assert(extension.is_valid(*this->graph_, &config_));

    return extension;
}

std::vector<Alignment> DefaultColumnExtender::backtrack(score_t min_path_score,
                                                        std::string_view window,
                                                        score_t right_end_bonus,
                                                        const std::vector<size_t> &tips,
                                                        node_index target_node) {
    assert(std::is_sorted(tips.begin(), tips.end()));
    std::vector<Alignment> extensions;
    size_t seed_clipping = this->seed_->get_clipping();
    ssize_t seed_offset = static_cast<ssize_t>(this->seed_->get_offset() - 1);
    ssize_t k_minus_1 = graph_->get_k() - 1;
    ssize_t last_pos = window.size();
    ssize_t seed_dist = std::max(graph_->get_k(), this->seed_->get_sequence().size()) - 1;
    score_t min_start_score = target_node ? ninf : min_path_score;
    size_t min_trace_length = this->graph_->get_k() - this->seed_->get_offset();

    std::vector<std::tuple<score_t, ssize_t, ssize_t, ssize_t>> indices;
    indices.reserve(table.size());
    auto it = tips.begin();
    for (size_t i = 1; i < table.size(); ++i) {
        while (it != tips.end() && i > *it) {
            ++it;
        }

        auto check_and_add_pos = [&](ssize_t start_pos, bool is_tip) {
            const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim,
                         xdrop_cutoff_i, score] = table[i];
            const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p,
                         xdrop_cutoff_i_p, score_p] = table[j_prev];
            if (start_pos < trim_p + 1)
                return;

            size_t pos = start_pos - trim;
            size_t pos_p = start_pos - trim_p - 1;

            if (S[pos] == ninf || S_p[pos_p] == ninf)
                return;

            score_t end_bonus = start_pos == last_pos ? right_end_bonus : 0;

            KmerExtractorBOSS::TAlphabet s = KmerExtractorBOSS::encode(c);
            if (target_node) {
                if (is_tip || (node == target_node
                        && S[pos] == S_p[pos_p] + score + profile_score_[s][seed_clipping + start_pos])) {
                    indices.emplace_back(S[pos] + end_bonus, -std::abs(start_pos - offset + seed_offset),
                                         -static_cast<ssize_t>(i), start_pos);
                }
            } else if (S[pos] + end_bonus >= min_start_score) {
                bool is_match = S[pos] == S_p[pos_p] + score + profile_score_[s][seed_clipping + start_pos]
                    && profile_op_[s][seed_clipping + start_pos] == Cigar::MATCH;
                if (is_match || start_pos == last_pos || is_tip) {
                    indices.emplace_back(S[pos] + end_bonus, -std::abs(start_pos - offset + seed_offset),
                                         -static_cast<ssize_t>(i), start_pos);
                }
            }
        };

        if (table[i].offset < seed_dist)
            continue;

        bool is_tip = (it != tips.end() && i == *it);
        if (!target_node || is_tip)
            check_and_add_pos(table[i].max_pos, is_tip);

        if (table[i].S.size() + table[i].trim == window.size() + 1
                && (target_node || table[i].max_pos != last_pos)) {
            check_and_add_pos(last_pos, is_tip);
        }
    }

    size_t num_backtracks = 0;

    // find highest scoring which is closest to the diagonal
    // use heap sort to make this run in O(n + (num_alternative_paths) * log(n)) time
    std::make_heap(indices.begin(), indices.end());

    score_t best_score = std::numeric_limits<score_t>::min();

    for (auto it = indices.rbegin(); it != indices.rend(); ++it) {
        std::pop_heap(indices.begin(), it.base());
        const auto &[start_score, neg_off_diag, neg_j_start, start_pos] = *it;
        DEBUG_LOG("Trying backtracking from score {}", start_score);

        if (terminate_backtrack_start(extensions))
            break;

        size_t j = -neg_j_start;

        if (skip_backtrack_start(j)) {
            DEBUG_LOG("\tskipped");
            continue;
        }

        std::vector<DeBruijnGraph::node_index> path;
        std::vector<size_t> trace;
        std::vector<score_t> score_trace;
        Cigar ops;
        std::string seq;
        score_t score = start_score;

        if (score - min_cell_score_ < best_score)
            break;

        ++num_backtracks;

        size_t dummy_counter = 0;

        ssize_t pos = start_pos;
        ssize_t end_pos = pos;
        size_t align_offset = this->seed_->get_offset();

        score_t extra_score = 0;
        auto append_node = [&](node_index node, char c, ssize_t offset, Cigar::Operator op) {
            seq += c;
            ops.append(op);
            if (offset >= k_minus_1) {
                path.emplace_back(node);
                if (!node) {
                    ++dummy_counter;
                } else if (dummy_counter) {
                    ops.append(Cigar::NODE_INSERTION, dummy_counter);
                    // extra_score -= config_.gap_opening_penalty + (dummy_counter - 1) * config_.gap_extension_penalty;
                    dummy_counter = 0;
                }
            }
        };

        while (j) {
            assert(j != static_cast<size_t>(-1));
            const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim,
                         xdrop_cutoff_i, score_cur] = table[j];
            const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p,
                         xdrop_cutoff_i_p, score_cur_p] = table[j_prev];

            assert(pos >= trim);
            assert(*std::max_element(S.begin(), S.end()) == S[max_pos - trim]);
            assert(!node || c == graph_->get_node_sequence(node)[std::min(k_minus_1, offset)]);

            align_offset = std::min(offset, k_minus_1);

            if (pos == max_pos)
                prev_starts.emplace(j);

            KmerExtractorBOSS::TAlphabet s = KmerExtractorBOSS::encode(c);

            if (S[pos - trim] == ninf) {
                j = 0;

            } else if (pos && S[pos - trim] == E[pos - trim] && (ops.empty()
                    || ops.data().back().first != Cigar::DELETION)) {
                // prefer insertions over match/mismatch to prevent premature
                // clipping of the alignment beginning
                Cigar::Operator last_op = Cigar::INSERTION;
                while (last_op == Cigar::INSERTION) {
                    ops.append(last_op);

                    assert(E[pos - trim] == E[pos - trim - 1] + config_.gap_extension_penalty
                        || E[pos - trim] == S[pos - trim - 1] + config_.gap_opening_penalty);

                    last_op = E[pos - trim] == E[pos - trim - 1] + config_.gap_extension_penalty
                        ? Cigar::INSERTION
                        : Cigar::MATCH;

                    --pos;
                }
            } else if (pos && pos >= trim_p + 1
                    && S[pos - trim] == S_p[pos - trim_p - 1] + score_cur
                        + profile_score_[s][seed_clipping + pos]) {
                // match/mismatch
                trace.emplace_back(j);
                score_trace.emplace_back(score_cur);

                extra_score += score_cur;
                append_node(node, c, offset, profile_op_[s][seed_clipping + pos]);
                --pos;
                assert(j_prev != static_cast<size_t>(-1));
                j = j_prev;

            } else if (S[pos - trim] == F[pos - trim] && (ops.empty()
                    || ops.data().back().first != Cigar::INSERTION)) {
                // deletion
                Cigar::Operator last_op = Cigar::DELETION;
                while (last_op == Cigar::DELETION && j) {
                    const auto &[S, E, F, node, j_prev, c, offset, max_pos, trim,
                                 xdrop_cutoff_i, score_cur] = table[j];
                    const auto &[S_p, E_p, F_p, node_p, j_prev_p, c_p, offset_p, max_pos_p, trim_p,
                                 xdrop_cutoff_i_p, score_cur_p] = table[j_prev];

                    align_offset = std::min(offset, k_minus_1);

                    assert(pos >= trim_p);

                    assert(F[pos - trim] == F_p[pos - trim_p] + score_cur + config_.gap_extension_penalty
                        || F[pos - trim] == S_p[pos - trim_p] + score_cur + config_.gap_opening_penalty);

                    last_op = F[pos - trim] == F_p[pos - trim_p] + score_cur + config_.gap_extension_penalty
                        ? Cigar::DELETION
                        : Cigar::MATCH;

                    trace.emplace_back(j);
                    score_trace.emplace_back(score_cur);
                    extra_score += score_cur;
                    append_node(node, c, offset, Cigar::DELETION);

                    assert(j_prev != static_cast<size_t>(-1));
                    j = j_prev;
                }
            } else {
                DEBUG_LOG("\tBacktracking failed, trying next start point");
                break;
            }
        }

        if (trace.size() >= min_trace_length && path.size() && path.back()) {
            assert(!dummy_counter);
            score_t cur_cell_score = table[j].S[pos - table[j].trim];
            best_score = std::max(best_score, score - cur_cell_score);
            if (score - min_cell_score_ < best_score)
                break;

            assert(extra_score == std::accumulate(score_trace.begin(), score_trace.end(),
                                                    score_t(0)));

            if (score >= min_start_score
                    && (!pos || cur_cell_score == 0)
                    && (pos || cur_cell_score == table[0].S[0])
                    && (config_.allow_left_trim || !j)) {
                call_alignments(score, path, trace, score_trace, ops, pos, align_offset,
                                window.substr(pos, end_pos - pos), seq, extra_score,
                                [&](Alignment&& alignment) {
                    DEBUG_LOG("Extension: {}", alignment);
                    extensions.emplace_back(std::move(alignment));
                });
            }
        }
    }

    DEBUG_LOG("Backtracked from {}/{} indices", num_backtracks, indices.size());

    if (extensions.empty() && this->seed_->get_score() >= min_path_score)
        extensions.emplace_back(*this->seed_);

    return extensions;
}

} // namespace align
} // namespace graph
} // namespace mtg
