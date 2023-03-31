#ifndef __ALT_ALIGNER__
#define __ALT_ALIGNER__

class AltCigar {
  public:
    enum Operator : int8_t {
        CLIPPED,
        MISMATCH,
        MATCH,
        DELETION,
        INSERTION,
        NODE_INSERTION
    };

    template <typename... Args>
    AltCigar(Args&&... args) : data_(std::forward<Args>(args)...) {}

    const std::vector<Operator>& data() const { return data_; }
    std::vector<Operator>& data() { return data_; }

    std::string to_string() const {
        std::string ret_val;
        ret_val.reserve(data_.size());
        for (Operator op : data_) {
            ret_val.emplace_back(op_str_[op]);
        }

        return ret_val;
    }

    size_t size() const { return data_.size(); }
    bool empty() const { return data_.empty(); }

  private:
    std::vector<Operator> data_;
    static constexpr std::string op_str_ = "SX=DIG";
};

class AltAnchor {
  public:
    using node_index = DeBruijnGraph::node_index;
    using Column = annot::binmat::BinaryMatrix::Column;
    using Columns = size_t;
    using Tuple = SmallVector<int64_t>;
    using CoordinateSet = Vector<Tuple>;
    using score_t = DBGAlignerConfig::score_t;

    virtual ~AltAnchor() {}

  protected:
    AltAnchor(const AltAnchor&) = default;
    AltAnchor(AltAnchor&&) noexcept = default;

    AltAnchor(std::string_view query_view = {},
              std::vector<node_index>&& nodes = {},
              bool orientation = false,
              score_t score = 0,
              size_t clipping = 0,
              size_t end_clipping = 0,
              AnnotationBuffer *buffer = nullptr,
              std::vector<Columns>&& columns = {},
              CoordinateSet&& coords = {})
          : query_view_(query_view), orientation_(orientation),
            score_(score), clipping_(clipping), end_clipping_(end_clipping), buffer_(buffer),
            columns_(std::move(columns)), coords_(std::move(coords)) {}

    operator bool() { return !nodes_.empty(); }

    std::string_view get_query_view() const { return query_view_ }
    const std::vector<node_index> get_nodes() const { return nodes_; }
    size_t size() const { return nodes_.size(); }

    virtual std::string_view get_spelling() const = 0;

    const Vector<Column>& get_columns(size_t i = size() - 1) const {
        assert(columns_.size() == nodes_.size());
        return buffer_ ? buffer_->get_cached_column_set(columns_[i]) : dummy_;
    }

    const CoordinateSet& get_coordinates() const { return coords_; }

    void set_columns(Columns columns, size_t i = size() - 1) { columns_[i] = columns; }

    void set_coordinates(const CoordinateSet &coords) {
        assert(buffer_);
        assert(coords.size() == get_columns(size() - 1).size());
        coords_ = coords;
        ssize_t offset = nodes_.size() - 1;
        for (auto &tuple : coords_) {
            for (auto &c : tuple) {
                c -= offset;
            }
        }
    }

  protected:
    std::string_view query_view_;
    std::vector<node_index> nodes_;
    bool orientation_;
    score_t score_;
    size_t clipping_;
    size_t end_clipping_;
    AnnotationBuffer *buffer_;
    std::vector<Columns> columns_;
    CoordinateSet coords_;

    static Vector<Column> dummy_;
}

class ExactMatch : public AltAnchor {
  public:
    template <typename... Args>
    ExactMatch(Args&&... args) : AltAnchor(std::forward<Args>(args)...) {
        assert(query_view_.size() >= nodes_.size());
        nodes_.insert(nodes_.begin(), query_view_.size() - nodes_.size(), DeBruijnGraph::npos);
        if (columns_.empty()) {
            columns_.resize(nodes_.size(), 0);
        } else {
            assert(columns_.size() == nodes_.size());
        }
    }

    std::string_view get_spelling() const override final { return query_view_; }
};

class AltAlignment : public AltAnchor {
  public:
    AltAlignment(const ExactMatch &match)
          : AltAnchor(match), spelling_(match.get_spelling()),
            cigar_(match.get_spelling().size(), AltCigar::MATCH) {
        assert(spelling_.size() == nodes_.size());
        if (extra_scores_.empty()) {
            extra_scores_.resize(nodes_.size(), 0);
        } else {
            assert(extra_scores_.size() == nodes_.size());
        }
    }

    AltAlignment(ExactMatch&& match)
          : AltAnchor(std::move(match)), spelling_(get_query_view()),
            cigar_(get_query_view().size(), AltCigar::MATCH) {
        assert(spelling_.size() == nodes_.size());
        if (extra_scores_.empty()) {
            extra_scores_.resize(nodes_.size(), 0);
        } else {
            assert(extra_scores_.size() == nodes_.size());
        }
    }

    AltAlignment(const ExactMatch &match, Column alt_col, int64_t alt_coord)
          : AltAnchor(match.query_view_, std::vector<node_index>(match.get_nodes()),
                      match.orientation_, match.score_, match.clipping_, match.end_clipping_,
                      match.buffer_),
            spelling_(match.get_spelling()), cigar_(match.get_spelling().size(), AltCigar::MATCH) {
        assert(spelling_.size() == nodes_.size());
        if (extra_scores_.empty()) {
            extra_scores_.resize(nodes_.size(), 0);
        } else {
            assert(extra_scores_.size() == nodes_.size());
        }

        columns_.assign(nodes_.size(),
            alt_col != std::numeric_limits<Column>::max()
                ? buffer_->cache_column_set(1, alt_col)
                : Columns(0)
        );

        if (alt_coord != std::numeric_limits<int64_t>::max()) {
            coords_.resize(1);
            coords_[0].assign(1, alt_coord);
        }
    }

    template <typename... Args>
    AltAlignment(std::string&& spelling,
                 AltCigar&& cigar,
                 std::vector<score_t>&& extra_scores,
                 Args&&... args)
          : AltAnchor(std::forward<Args>(args)...),
            spelling_(std::move(spelling)),
            cigar_(std::move(cigar)),
            extra_scores_(std::move(extra_scores)) {
        assert(spelling_.size() >= nodes_.size());
        nodes_.insert(nodes_.begin(), spelling_.size() - nodes_.size(), DeBruijnGraph::npos);
        if (extra_scores_.empty()) {
            extra_scores_.resize(nodes_.size(), 0);
        } else {
            assert(extra_scores_.size() == nodes_.size());
        }

        if (columns_.empty()) {
            columns_.resize(nodes_.size(), 0);
        } else {
            assert(columns_.size() == nodes_.size());
        }
    }

    std::string_view get_spelling() const override final { return spelling_; }
    const AltCigar& get_cigar() const { return cigar_; }

  private:
    std::string spelling_;
    AltCigar cigar_;
    std::vector<score_t> extra_scores_;
}

class AltAlignmentResult {
  public:
    AltAlignmentResult(std::string_view header, std::string_view query)
          : header_(header), fwd_(query), rc_(query) {
        reverse_complement(rc_.begin(), rc_.end());
    }

    std::string_view get_header() const { return header_; }
    std::string_view get_query_view(bool orientation) const { return orientation ? rc_ : fwd_; }
    const std::vector<AltAlignment> get_alignments() const { return alignments_; }

    template <typename... Args>
    AltAlignment& emplace_back(Args&&... args) {
        alignments_.emplace_back(std::forward<Args>(args)...);
    }

    template <typename MoveIt>
    void extend(MoveIt begin, MoveIt end) {
        alignments_.insert(alignments_.end(), begin, end);
    }

  private:
    std::vector<AltAlignment> alignments_;
    std::string header_;
    std::string fwd_;
    std::string rc_;
}

class AltDBGAligner : public IDBGAligner {
  public:
    const DeBruijnGraph& get_graph() const override final { return graph_; }
    const DBGAlignerConfig& get_config() const override final { return config_; }

    // Main aligner
    void align_batch(const std::vector<Query> &,
                     const AlignmentCallback &,
                     size_t = 0) const {}

    std::unique_ptr<SeedFilteringExtender> make_extender(std::string_view) const {
        return {};
    }

    bool has_coordinates() const { return buffer_ && buffer->has_coordinates(); }

    // Construct a full alignment from a chain by aligning the query against
    // the graph in the regions of the query in between the chain seeds.
    void extend_chain(Chain&&,
                      SeedFilteringExtender &,
                      const std::function<void(Alignment&&)> &,
                      bool = true) const {}

    AltDBGAligner(const DeBruijnGraph &graph,
                  const DBGAlignerConfig &config,
                  AnnotationBuffer *buffer = nullptr)
          : graph_(graph), config_(config), buffer_(buffer) {}

    std::vector<AltAlignmentResult> align_batch(const std::vector<Query> &seq_batch) const {
        std::vector<AltAlignmentResult> results;
        results.reserve(seq_batch.size());
        for (const auto &[header, query] : seq_batch) {
            results.emplace_back(header, query);
        }

        // anchoring
        std::vector<std::vector<ExactMatch>> fwd_anchors(seq_batch.size());
        std::vector<std::vector<ExactMatch>> rc_anchors(seq_batch.size());
        ProgressBar progress_bar(seq_batch.size(), "Anchoring sequences",
                                 std::cerr, !common::get_verbose());
        std::vector<node_index> all_nodes;
        for (size_t i = 0; i < results.size(); ++i) {
            std::string_view fwd = results[i].get_query_view(false);
            std::string_view rc = results[i].get_query_view(true);

            if (config_.min_seed_length >= graph_.get_k()) {
                auto fwd_walk = map_to_nodes_sequentially(graph_, fwd);
                for (size_t j = 0; j + graph_.get_k() <= fwd.size(); ++j) {
                    if (fwd_walk[j]) {
                        if (buffer_)
                            all_nodes.emplace_back(fwd_walk[j]);

                        std::string_view fwd_view(fwd.begin(), graph_.get_k());
                        fwd_anchors[i].emplace_back(fwd_view,
                                                    std::vector<node_index>{ fwd_walk[j] },
                                                    false,
                                                    config_.match_score(fwd_view),
                                                    j,
                                                    fwd.size() - j - graph_.get_k(),
                                                    buffer_);
                    }
                }

                auto rc_walk = map_to_nodes_sequentially(graph_, bwd);
                for (size_t j = 0; j + graph_.get_k() <= rc.size(); ++j) {
                    if (rc_walk[j]) {
                        if (buffer_)
                            all_nodes.emplace_back(rc_walk[j]);

                        std::string_view rc_view(rc.begin(), graph_.get_k());
                        rc_anchors[i].emplace_back(rc_view,
                                                   std::vector<node_index>{ rc_walk[j] },
                                                   true,
                                                   config_.match_score(rc_view),
                                                   j,
                                                   rc.size() - j - graph_.get_k(),
                                                   buffer_);
                    }
                }

            } else if (const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph_)) {
                fwd_seeds[i].reserve(fwd.size() - config_.min_seed_length + 1);
                rc_seeds[i].reserve(rc.size() - config_.min_seed_length + 1);
                for (size_t j = 0; j + config_.min_seed_length <= fwd.size(); ++j) {
                    std::string_view fwd_view(fwd.begin() + j, config_.min_seed_length);
                    score_t score = config_.match_score(fwd_view);
                    dbg_succ->call_nodes_with_suffix_matching_longest_prefix(
                        fwd_view,
                        [&](node_index node, size_t) {
                            if (path_index && dbg_succ->get_node_sequence(node).find(boss::BOSS::kSentinel) != std::string::npos)
                                return;

                            if (buffer_)
                                all_nodes.push_back(node);

                            fwd_anchors[i].emplace_back(fwd_view,
                                                        std::vector<node_index>{ node },
                                                        false,
                                                        score,
                                                        j,
                                                        fwd.size() - j - graph_.get_k(),
                                                        buffer_);
                        },
                        config_.min_seed_length
                    );

                    std::string_view rc_view(rc.begin() + j, config_.min_seed_length);
                    score = config_.match_score(rc_view);
                    dbg_succ->call_nodes_with_suffix_matching_longest_prefix(
                        rc_view,
                        [&](node_index node, size_t) {
                            if (path_index && dbg_succ->get_node_sequence(node).find(boss::BOSS::kSentinel) != std::string::npos)
                                return;

                            if (buffer_)
                                all_nodes.push_back(node);

                            rc_anchors[i].emplace_back(rc_view,
                                                       std::vector<node_index>{ node },
                                                       true,
                                                       score,
                                                       j,
                                                       rc.size() - j - graph_.get_k(),
                                                       buffer_);
                        },
                        config_.min_seed_length
                    );
                }

            } else {
                throw std::runtime_error("Sub-k seeding only supported with DBGSuccinct");
            }

            ++progress_bar;
        }

        if (buffer_) {
            // fetching annotations
            logger->trace("Preparing nodes to fetch annotations");
            buffer_->queue_path(std::move(all_nodes));

            logger->trace("Fetching annotations");
            buffer_->fetch_queued_annotations();

            logger->trace("Annotating anchors");
            auto it = all_nodes.begin();
            for (size_t i = 0; i < results.size(); ++i) {
                for (auto &anchor : fwd_anchors[i]) {
                    if (size_t columns = buffer_->get_labels_id(anchor.get_nodes().back())) {
                        anchor.set_columns(columns);
                    } else {
                        assert(it != all_nodes.end());
                        *it = DeBruijnGraph::npos;
                        anchor = ExactMatch();
                    }
                    ++it;
                }

                for (auto &anchor : rc_anchors[i]) {
                    if (size_t columns = buffer_->get_labels_id(anchor.get_nodes().back())) {
                        anchor.set_columns(columns);
                    } else {
                        assert(it != all_nodes.end());
                        *it = DeBruijnGraph::npos;
                        anchor = ExactMatch();
                    }
                    ++it;
                }

                if (has_coordinates()) {
                    for (auto &anchor : fwd_anchors[i]) {
                        if (anchor.size())
                            anchor.set_coordinates(*buffer_->get_lables_and_coords(anchor.get_nodes().back()).second);
                    }

                    for (auto &anchor : rc_anchors[i]) {
                        if (anchor.size())
                            anchor.set_coordinates(*buffer_->get_lables_and_coords(anchor.get_nodes().back()).second);
                    }
                }

                fwd_anchors[i].erase(std::remove_if(fwd_anchors[i].begin(),
                                                    fwd_anchors[i].end(),
                                                    [](const auto &a) { return a.empty(); }),
                                     fwd_anchors[i].end());

                rc_anchors[i].erase(std::remove_if(rc_anchors[i].begin(),
                                                   rc_anchors[i].end(),
                                                   [](const auto &a) { return a.empty(); }),
                                    rc_anchors[i].end());
            }
            assert(it == all_nodes.end());
            all_nodes.erase(std::remove_if(all_nodes.begin(), all_nodes.end(),
                                           [](node_index n) { return n == DeBruijnGraph::npos; }),
                            all_nodes.end());
        }

        const auto *path_index = graph_.get_extension_threadsafe<IPathIndex>();
        for (size_t i = 0; i < results.size(); ++i) {
            std::string_view fwd = results[i].get_query_view(false);
            std::string_view rc = results[i].get_query_view(true);
            auto &fwd_anchor = fwd_anchors[i];
            auto &rc_anchor = rc_anchors[i];

            std::vector<AltAlignment> alignments;
            if (has_coordinates() || path_index) {

            } else {

            }

            tsl::hopscotch_map<node_index, std::vector<size_t>> node_to_anchor;
            logger->trace("Extending alignment ends");
            for (auto &alignment : alignments) {
                throw std::runtime_error("Implement end extension");
            }

            if (config_.allow_jump) {
                logger->trace("Chaining alignments");
                throw std::runtime_error("Implement alignment chaining");
            }

            results[i].extend(std::make_move_iterator(alignments.begin()),
                              std::make_move_iterator(alignments.end()));
        }

        return results;
    }

    AltAlignment connect_anchors(const AltAnchor &first,
                                 const AltAnchor &second,
                                 size_t dist) const {
        size_t n = dist + first.get_sequence().size();
        std::string_view query_window(
            first.get_query_view().begin(),
            second.get_query_view().end() - first.get_query_view().begin()
        );
        size_t m = query_window.size();

        // WFA
        using traceback_t = std::tuple<node_index, int64_t, size_t>;
        using WFS = tsl::hopscotch_map<node_index, tsl::hopscotch_map<int64_t, std::pair<bool, traceback_t>>>;
        std::vector<WFS> table;
        table.emplace_back();

        size_t gap_open = -2 * config_.gap_opening_penalty;
        size_t gap_extend = -2 * config_.gap_extension_penalty + config_.match_score("A");
        size_t mismatch = 2 * config_.match_score("A") + 2 * config_.score_sequences("A", "C");

        std::vector<std::pair<node_index, int64_t>> stack;
        stack.emplace_back(first.get_nodes().back(), 0);
        size_t cost = 0;

        while (cost < table.size()) {
            auto &wavefront = table[cost];
            //extend
            {
                std::vector<std::pair<node_index, int64_t>> stack_next;
                while (stack.size()) {
                    auto [n, k] = stack.back();
                    stack.pop_back();
                    node_index last = n;
                    if (!wavefront[n][k].first) {
                        size_t i = k;
                        graph_.traverse(
                            n, query_window.data() + i, query_window.data() + query_window.size(),
                            [&](node_index next) {
                                wavefront[next][k].first = true;
                                wavefront[next][k].second = std::make_tuple(last, k, cost);
                                last = next;
                            },
                            [&]() {
                                return graph_.has_multiple_outgoing(last)
                                    || graph_.indegree(last) > 1;
                            }
                        );
                    }
                    stack_next.emplace_back(last, k);
                    graph_.call_outgoing_kmers(last, [&](node_index next, char c) {
                        if (c == boss::BOSS::kSentinel)
                            return;

                        if (!wavefront[next][k + 1].first) {
                            wavefront[next][k + 1].second = std::make_tuple(last, k, cost);
                            stack.emplace_back(next, k + 1);
                        }
                    });
                }

                std::swap(stack, stack_next);
            }

            //check
            auto find_n = table[cost].find(second.get_nodes().back());
            if (find_n != table[cost].end()) {
                auto find_k = find_n->second.find(m - 1);
                if (find_k != find_n->second.end()) {
                    if (find_k->second.first)
                        break;
                }
            }

            // expand
            {
                std::vector<WFS> table_next;
                std::vector<std::tuple<size_t, node_index, int64_t>> stack_next;
                while (stack.size()) {
                    auto [n, k] = stack.back();
                    stack.pop_back();

                    if (cost + gap_open <= table.size())
                        table.resize(cost + gap_open + 1);

                    if (cost + gap_open <= table_next.size())
                        table_next.resize(cost + gap_open + 1);

                    if (cost + mismatch <= table.size())
                        table.resize(cost + mismatch + 1);

                    if (cost + mismatch <= table_next.size())
                        table_next.resize(cost + mismatch + 1);

                    // insert
                    {
                        size_t cost_next = cost + gap_open;
                        auto &w = table[cost_next];
                        auto &wp = table_next[cost_next];
                        size_t i = k + w[n][k].first;
                        if (i < m) {
                            if (!wp[n][k + 1].first) {
                                if (w[n][k + 1].first < w[n][k].first) {
                                    wp[n][k + 1].first = true;
                                    wp[n][k + 1].second = std::make_tuple(n, k, cost);
                                    stack_next.emplace_back(cost_next, n, k + 1);
                                }
                            } else {
                                wp[n][k + 1].first |= w[n][k].first;
                                if (w[n][k].first)
                                    wp[n][k + 1].second = std::make_tuple(n, k, cost);
                            }
                        }
                    }

                    // delete
                    {
                        size_t cost_next = cost + gap_open;
                        auto &w = table[cost_next];
                        auto &wp = table_next[cost_next];
                        size_t i = k + w[n][k].first;
                        if (!w[n][k].first) {
                            if (!wp[n][k - 1].first) {
                                if (!w[n][k - 1].first) {
                                    wp[n][k - 1].first = true;
                                    wp[n][k - 1].second = std::make_tuple(n, k, cost);
                                    stack_next.emplace_back(cost_next, n, k - 1);
                                }
                            } else {
                                wp[n][k - 1].first = true;
                                wp[n][k - 1].second = std::make_tuple(n, k, cost);
                            }
                        }
                    }

                    // mismatch
                    {
                        size_t cost_next = cost + mismatch;
                        auto &w = table[cost_next];
                        auto &wp = table_next[cost_next];
                        size_t i = k + w[n][k].first;
                        if (i < m && !w[n][k].first) {
                            if (!wp[n][k].first)
                                stack_next.emplace_back(cost_next, n, k);

                            wp[n][k].first = true;
                            wp[n][k].second = std::make_tuple(n, k, cost);
                        }
                    }
                }

                for (const auto &[cost_next, n, k] : stack_next) {
                    table[cost_next][n][k] = table_next[cost_next][n][k];
                }
            }

            ++cost;
            while (cost < table.size() && table[cost].empty()) {
                ++cost;
            }
        }

        if (cost >= table.size())
            throw std::runtime_error("No alignment found");

        // traceback
        node_index cur = second.get_nodes().back();
        int64_t k = m - 1;
        AltCigar cigar;
        std::string spelling;
        auto it = query_window.rbegin();
        while (cur) {
            auto [last_n, last_k, last_cost] = table[cost][cur][m - 1].second;
            if (last_cost == cost) {
                // match
                assert(last_k == k);
                spelling.emplace_back(*it);
                ++it;
            } else if (last_k == k + 1) {
                // delete
                assert(last_cost == cost - gap_open);
            } else if (last_k == k - 1) {
                // insert
                assert(last_cost == cost - gap_open);
            } else {
                // mismatch
                assert(last_cost == cost - mismatch);
            }
            std::tie(cur, k, cost) = std::tie(last_n, last_k, last_cost);
        }

        throw std::runtime_error("Backtracking not implemented");

        return {};
    }

  private:
    const DeBruijnGraph &graph_;
    DBGAlignerConfig config_;
    mutable AnnotationBuffer *buffer_;
};

// Myers
// std::vector<uint8_t> encoder(255, graph_.alphabet().size() + 1);
//         for (size_t i = 0; i < graph_.alphabet().size(); ++i) {
//             encoder[alphabet[c]] = i;
//         }

//         std::vector<mpz_class> B;
//         B.resize(graph_.alphabet().size(), mpz_class("0"));
//         for (size_t j = 0; j < m; ++j) {
//             mpz_class id("1");
//             id <<= m - j;
//             B[query_window[j]] |= id;
//         }

//         struct TableElem {
//             size_t i;
//             size_t last;
//             mpz_class VP;
//             mpz_class VN;
//             mpz_class D0;
//             mpz_class HN;
//             mpz_class HP;
//             score_t score;
//             node_index node;
//             uint8_t c_enc;
//         };

//         std::vector<TableElem> table;
//         table.emplace_back(TableElem{
//             .i = std::numeric_limits<size_t>::max(),
//             .last = std::numeric_limits<size_t>::max(),
//             .VP = (mpz_class("1") << m) - 1,
//             .VN = mpz_class("0"),
//             .score = static_cast<score_t>(m),
//             .node = first.get_nodes().back(),
//             .c_enc = encoder[first.get_sequence().back()]
//         });

//         std::priority_queue<std::pair<score_t, uint32_t>> queue;
//         queue.emplace(-table.back().score, 0);

//         std::vector<std::pair<score_t, size_t>> endpoints;

//         while (queue.size()) {
//             auto [last_nscore, last_table_i] = queue.top();
//             queue.pop();

//             size_t i = table[last_table_i].i + 1;
//             score_t score = table[last_table_i].score;
//             graph_.call_outgoing_kmers(table[last_table_i].node, [&](node_index next, char c) {
//                 if (c == boss::BOSS::kSentinel)
//                     return;

//                 if (i + 1 == dist && next != second.get_nodes().back())
//                     return;

//                 uint8_t c_enc = encoder[c];

//                 auto &next = table.emplace_back(TableElem{
//                     .i = i,
//                     .last = last_table_i,
//                     .VP = mpz_class(),
//                     .VN = mpz_class(),
//                     .score = score,
//                     .node = next,
//                     .c_enc = c_enc
//                 });

//                 const auto &last = table[last_table_i];

//                 mpz_class X = B[c_enc] | last.VN;
//                 mpz_class D0 = ((last.VP + (X & last.VP)) ^ last.VP) | X;
//                 mpz_class HN = last.VP & D0;
//                 mpz_class HP = last.VN | ~(last.VP | D0);
//                 X = HP << 1;
//                 next.VN = X & D0;
//                 next.VP = (HN << 1) | ~(X | D0);

//                 if (HP & 1) {
//                     ++next.score;
//                 } else if (HN & 1) {
//                     --next.score;
//                 }

//                 if (i + 1 == dist) {
//                     // reached the end
//                     endpoints.emplace_back(next.score, table.size() - 1);
//                 } else {
//                     queue.emplace(-next.score, table.size() - 1);
//                 }
//             });
//         }

#endif // __ALT_ALIGNER__
