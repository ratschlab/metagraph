#include "walker.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>
#include <folly/hash/FnvHash.h>

#include "resolve.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/annotation_key_mapper.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/unix_tools.hpp"
#include "common/vector.hpp"


namespace mtg {
namespace graph {
namespace traversal {

const char* to_string(Arm arm) {
    return arm == Arm::LEFT ? "left" : "right";
}

const char* to_string(EndReason reason) {
    switch (reason) {
        case EndReason::DEAD_END: return "dead_end";
        case EndReason::LABEL_LOST: return "label_lost";
        case EndReason::RECORD_END: return "record_end";
        case EndReason::BRANCH: return "branch";
        case EndReason::EDGE_REUSE: return "edge_reuse";
        case EndReason::EDGE_REUSE_RC: return "edge_reuse_rc";
        case EndReason::REACHED_SEED: return "rejoined_seed";
        case EndReason::MAX_EXTENSION: return "max_extension_bp";
        case EndReason::MAX_STEPS: return "max_steps";
        case EndReason::MAX_LIVE_PATHS: return "max_live_paths";
        case EndReason::MAX_PATHS: return "max_paths";
        case EndReason::MAX_SPLITS: return "max_splits_per_path";
        case EndReason::TIME_BUDGET: return "time_budget";
    }
    return "unknown";
}

const char* to_string(SeedDerivationError::Cause cause) {
    switch (cause) {
        case SeedDerivationError::NO_CARRIER: return "no_carrier";
        case SeedDerivationError::TOO_WIDE: return "too_wide";
        case SeedDerivationError::TIME_BUDGET: return "time_budget";
        case SeedDerivationError::AMBIGUOUS_HEADER: return "ambiguous_header";
    }
    return "unknown";
}

namespace {

constexpr char kSentinel = boss::BOSS::kSentinel;
// the k-mers whose rows a derivation decodes at once; the narrowest of the first window
// starts the intersection
constexpr size_t kDeriveWindow = 64;
// the rows the check of named seed labels decodes at once
constexpr size_t kSeedBatchRows = 4096;

// One live label on one path head
struct Entry {
    LabelId label = 0;
    uint32_t branches = 0;
    EndReason reason = EndReason::DEAD_END;   // set when the label ends
    // A HEADER label is followed in coordinate order: |coord| is the coordinate of the
    // head's k-mer in the label's record and |step| (+1 or -1) what one step of this arm
    // adds to it. A COLUMN label has step 0.
    int8_t step = 0;
    uint64_t coord = 0;
};
// sorted by label id
using State = std::vector<Entry>;

std::vector<LabelId> labels_of(const State &state) {
    std::vector<LabelId> out;
    out.reserve(state.size());
    for (const Entry &e : state) {
        out.push_back(e.label);
    }
    return out;
}

// Where a HEADER label's record holds the seed: seed k-mer j at coordinate origin + dir * j
// (dir -1: the record holds the seed's reverse complement). dir 0: nowhere.
struct Place {
    uint64_t origin = 0;
    int8_t dir = 0;
};

// A live path head
struct Item {
    size_t segment = 0;
    node_index node = npos;
    std::string kmer;               // spelling of |node| (uppercase, natural orientation)
    uint64_t ext_bp = 0;
    State state;
    uint32_t splits = 0;
    size_t path_id = 0;
};

struct Succ {
    node_index node;
    char ch;
};

// An edge ((k+1)-mer): 2-bit packed when it fits 64 bits, otherwise a 128-bit FNV-1a pair
// (two distinct edges colliding on it are treated as one: conservative, it blocks)
struct EdgeKey {
    uint64_t lo;
    uint64_t hi;
    bool operator==(const EdgeKey &o) const { return lo == o.lo && hi == o.hi; }
};
struct EdgeKeyHash {
    size_t operator()(const EdgeKey &k) const {
        return k.lo ^ (k.hi * 0x9E3779B97F4A7C15ULL) ^ (k.hi >> 29);
    }
};
// one use of an edge: the segment that took it and the orientation bit
struct EdgeUse {
    size_t segment;
    bool rc;
};
// the same use keyed the other way round, for probing "did THIS segment take the edge"
struct EdgeSegKey {
    EdgeKey edge;
    size_t segment;
    bool operator==(const EdgeSegKey &o) const { return edge == o.edge && segment == o.segment; }
};
struct EdgeSegKeyHash {
    size_t operator()(const EdgeSegKey &k) const {
        return EdgeKeyHash()(k.edge) ^ (k.segment * 0xC2B2AE3D27D4EB4FULL);
    }
};

// A successor under evaluation
struct Cand {
    const Succ *succ = nullptr;
    std::string kmer;
    State state;                    // the head's labels that may continue on the successor
    State barred;                   // its COLUMN labels there that the walk rule stops
    bool ruled = false;             // a COLUMN label is on the successor: the walk rule applies
    bool skipped = false;           // a hairpin, not followed by COLUMN labels
    bool blocked = false;
    EndReason block_reason = EndReason::DEAD_END;
    EdgeKey edge { 0, 0 };
    bool rc = false;
    bool followed = false;
    bool split_limited = false;     // refused by the path's split allowance
    bool admissible() const { return !skipped && !blocked; }
};

// a blocked successor's reason, the strongest of several winning
uint8_t block_rank(EndReason reason) {
    switch (reason) {
        case EndReason::REACHED_SEED: return 3;
        case EndReason::EDGE_REUSE_RC: return 2;
        case EndReason::EDGE_REUSE: return 1;
        default: return 0;
    }
}

struct LeafInfo {
    bool is_leaf = false;
    std::optional<EndReason> path_reason;
    std::vector<PathEnd> end_labels;
};

struct ArmState {
    Arm arm = Arm::RIGHT;
    ArmResult result;
    std::vector<Item> frontier;
    std::vector<Item> next;
    std::vector<std::string> walk_seq;      // per segment, in walking order
    std::vector<LeafInfo> leaves;           // per segment
    // The per-path edge-reuse record, kept both ways: every use of an edge (scanned when
    // the edge has few uses) and, per (edge, segment), the orientation bits 1 << rc of that
    // segment's uses (probed per ancestor when it has many: every live path through a
    // shared region records the same edge, and scanning those uses once per path would be
    // quadratic in the number of paths)
    tsl::hopscotch_map<EdgeKey, std::vector<EdgeUse>, EdgeKeyHash> used_edges;
    tsl::hopscotch_map<EdgeSegKey, uint8_t, EdgeSegKeyHash> used_by_segment;
    // the ancestors of |marked_segment|, as marks and as a list
    std::vector<uint32_t> visit_mark;
    uint32_t visit_epoch = 0;
    size_t marked_segment = SIZE_MAX;
    std::vector<size_t> ancestors;
    size_t finished_leaves = 0;
    size_t next_path_id = 1;
    bool stopped = false;
    // extension depth of the first head not expanded: its level is partial
    uint64_t boundary = std::numeric_limits<uint64_t>::max();
};


class Walker {
  public:
    Walker(const LabelOracle &oracle, const Seed &seed, const Strategy &strategy)
          : oracle_(oracle), graph_(oracle.graph()), key_mapper_(graph_), seed_(seed),
            strategy_(strategy), k_(graph_.get_k()), regime_(oracle.regime()),
            canonical_(dynamic_cast<const CanonicalDBG*>(&graph_)) {
        size_t n = 0;
        code_.fill(255);
        for (char c : graph_.alphabet()) {
            if (c != kSentinel)
                code_[static_cast<unsigned char>(c)] = n++;
        }
        bits_ = 1;
        while ((size_t(1) << bits_) < n) {
            ++bits_;
        }
        packable_ = (k_ + 1) * bits_ <= 64;
    }

    SeedResult run();

  private:
    const LabelOracle &oracle_;
    const DeBruijnGraph &graph_;
    const AnnotationKeyMapper key_mapper_;
    const Seed &seed_;
    const Strategy &strategy_;
    const size_t k_;
    const Regime regime_;
    const CanonicalDBG *canonical_;
    Timer timer_;

    std::array<uint8_t, 256> code_;
    size_t bits_;
    bool packable_;

    SeedResult result_;
    std::string seed_upper_;
    std::vector<node_index> nodes_;
    tsl::hopscotch_set<node_index> seed_nodes_;
    // per label of the dictionary: where a HEADER label's record holds the seed
    std::vector<Place> places_;
    std::unique_ptr<LabelQuery> query_;
    std::array<ArmState, 2> arms_;
    uint64_t steps_total_ = 0;
    bool seed_stopped_ = false;
    double cap_demand_ = 0;

    double elapsed_ms() const { return timer_.elapsed() * 1000; }
    bool time_exceeded() const {
        return strategy_.time_budget_ms <= 0 || elapsed_ms() >= strategy_.time_budget_ms;
    }

    // seed phase
    void validate_seed();
    void scan_named_labels(const std::vector<LabelRef> &refs,
                           const std::vector<node_index> &keys,
                           std::vector<std::vector<KmerInterval>> *runs,
                           std::vector<Place> *places) const;
    void derive_seed_labels(const std::vector<node_index> &keys, std::vector<LabelRef> *refs);
    void init_arm(ArmState &arm);

    // graph access
    std::vector<Succ> successors(Arm arm, node_index node, const std::string &kmer) const;
    void succ_kmer(Arm arm, const std::string &kmer, char c, std::string *out) const;
    void step_kmer(Arm arm, const std::string &kmer, char c, std::string *out) const;
    bool is_hairpin(const std::string &kmer_u, const std::string &kmer_v,
                    const std::string &step) const;
    std::pair<EdgeKey, bool> edge_key(const std::string &step) const;

    // bookkeeping
    size_t new_segment(ArmState &arm, std::vector<size_t> parents, uint64_t from_bp,
                       std::vector<LabelId> labels_start);
    void mark_ancestors(ArmState &arm, size_t segment);
    bool is_ancestor_or_self(const ArmState &arm, size_t anc, size_t segment) const;
    void record_edge(ArmState &arm, size_t segment, const Cand &c);
    void end_label(ArmState &arm, Item &item, Entry &e, EndReason reason, const char *text);
    void finish_path(ArmState &arm, Item &item, std::optional<EndReason> path_reason);
    void censor_item(ArmState &arm, Item &item, EndReason reason);

    // the walk
    // whether a HEADER label's record ends at the head's k-mer, in the walking direction
    bool at_record_end(const Entry &e) const {
        const std::pair<uint64_t, uint64_t> &record = query_->record_coords(e.label);
        return e.coord == (e.step > 0 ? record.second : record.first);
    }
    void run_level(ArmState &arm, uint64_t depth);
    void check_structure(ArmState &arm, const Item &item, Cand &c);
    std::optional<EndReason> cap_check(const ArmState &arm, size_t nf, size_t remaining);
    std::optional<EndReason> process_item(ArmState &arm, Item &item,
                                          const std::vector<Succ> &succs,
                                          const LabelQuery::Hits *const *hits,
                                          size_t remaining_in_level);
    void stop_arm(ArmState &arm, EndReason reason, std::vector<Item> *items, size_t from);
    void stop_frontier(ArmState &arm, EndReason reason);
    void trip(ArmState &arm, EndReason reason, std::vector<Item> &items, size_t index);
    void finalize(ArmState &arm);
};


/******************************** seed phase ********************************/

void Walker::validate_seed() {
    check_seed(oracle_, seed_);
    const std::string &seq = seed_.sequence;
    seed_upper_ = in_graph_case(seq);

    nodes_ = map_to_nodes_sequentially(graph_, seed_upper_);
    assert(nodes_.size() == seq.size() - k_ + 1);
    if (std::find(nodes_.begin(), nodes_.end(), npos) != nodes_.end()) {
        sdsl::bit_vector present(nodes_.size());
        for (size_t i = 0; i < nodes_.size(); ++i) {
            present[i] = nodes_[i] != npos;
        }
        throw SeedNotInGraph("Seed is not fully present in the graph; graph runs: "
                             + encode_presence_mask(present), nodes_.size(),
                             sdsl::util::cnt_one_bits(present));
    }
    result_.seed_id = seed_.seed_id;
    result_.length_bp = seq.size();
    result_.num_kmers = nodes_.size();

    // the seed's nodes, in both orientations in the canonical regimes: a step into any of
    // them re-enters the seed
    for (node_index n : nodes_) {
        seed_nodes_.insert(n);
    }
    if (regime_ != Regime::BASIC) {
        std::string rc = seed_upper_;
        ::reverse_complement(rc);
        for (node_index n : map_to_nodes_sequentially(graph_, rc)) {
            if (n != npos)
                seed_nodes_.insert(n);
        }
    }

    places_.clear();
    const std::vector<node_index> keys = map_to_nodes(graph_, seed_upper_);
    assert(keys.size() == nodes_.size());
    std::vector<LabelRef> refs;
    if (seed_.labels.empty()) {
        // every derived label carries the seed
        derive_seed_labels(keys, &refs);
        result_.label_dict = refs;
    } else {
        // A label that does not carry the seed is dropped, with the runs of k-mers it is
        // on (check_seed refused duplicates). A column carries the seed when it is on
        // every k-mer, a record when it holds the k-mers at consecutive coordinates.
        for (const auto &name : seed_.labels) {
            refs.push_back(oracle_.resolve_label(name));
        }
        std::vector<std::vector<KmerInterval>> runs;
        std::vector<Place> places;
        scan_named_labels(refs, keys, &runs, &places);
        const KmerInterval whole { 0, nodes_.size() };
        for (size_t l = 0; l < refs.size(); ++l) {
            const bool carries = refs[l].kind == LabelKind::HEADER
                    ? places[l].dir != 0
                    : runs[l].size() == 1 && runs[l][0] == whole;
            if (carries) {
                result_.label_dict.push_back(std::move(refs[l]));
                places_.push_back(places[l]);
            } else {
                result_.dropped_labels.push_back({ refs[l].name, std::move(runs[l]) });
            }
        }
        if (result_.label_dict.empty()) {
            throw SeedDerivationError(SeedDerivationError::NO_CARRIER,
                                      "No seed label supports every k-mer of the seed",
                                      nodes_.size(), nodes_.size());
        }
        result_.labels_supporting_total = result_.label_dict.size();
    }
    assert(places_.size() == result_.label_dict.size());

    std::vector<std::string> names;
    for (const LabelRef &ref : result_.label_dict) {
        names.push_back(ref.name);
    }
    result_.validated_seed_id = make_seed_id(seed_upper_, regime_ != Regime::BASIC, names);
    result_.seed_id_mismatch = !seed_.seed_id.empty()
                                && seed_.seed_id != result_.validated_seed_id;
}

// For every named label, the runs of seed k-mers it is on, and for a HEADER label where
// its record holds the seed (Place). Every coordinate of the record at the first k-mer is
// a candidate origin, read forwards and, in the canonical regimes, backwards; each next
// k-mer keeps the origins whose next coordinate it holds. Of the places left, the one with
// the lowest coordinate is taken (its origin forwards, its last k-mer's coordinate
// backwards; forwards first when equal). So a seed and its reverse complement, which are
// one seed in the canonical regimes, are placed on the same occurrence.
void Walker::scan_named_labels(const std::vector<LabelRef> &refs,
                               const std::vector<node_index> &keys,
                               std::vector<std::vector<KmerInterval>> *runs,
                               std::vector<Place> *places) const {
    LabelQuery query(oracle_, refs);
    runs->assign(refs.size(), {});
    places->assign(refs.size(), Place());
    // per header label, the candidate origins, ascending: forwards [0] and backwards [1]
    std::vector<std::array<std::vector<uint64_t>, 2>> origins(refs.size());
    for (size_t begin = 0; begin < keys.size(); begin += kSeedBatchRows) {
        const size_t end = std::min(keys.size(), begin + kSeedBatchRows);
        const std::vector<const LabelQuery::Hits*> hits = query.fetch_hits(
                std::vector<node_index>(keys.begin() + begin, keys.begin() + end));
        for (uint64_t j = begin; j < end; ++j) {
            const LabelQuery::Hits &hit = *hits[j - begin];
            for (LabelId l : hit.labels) {
                std::vector<KmerInterval> &run = (*runs)[l];
                if (!run.empty() && run.back().end == j) {
                    run.back().end = j + 1;
                } else {
                    run.push_back(KmerInterval{ j, j + 1 });
                }
            }
            auto it = hit.coords.begin();
            for (LabelId l = 0; l < refs.size(); ++l) {
                if (refs[l].kind != LabelKind::HEADER)
                    continue;
                // the label's coordinates at this k-mer
                while (it != hit.coords.end() && it->label < l) {
                    ++it;
                }
                const auto first = it;
                while (it != hit.coords.end() && it->label == l) {
                    ++it;
                }
                const auto last = it;
                std::vector<uint64_t> &forwards = origins[l][0];
                std::vector<uint64_t> &backwards = origins[l][1];
                if (!j) {
                    for (auto c = first; c != last; ++c) {
                        forwards.push_back(c->coord);
                        if (regime_ != Regime::BASIC)
                            backwards.push_back(c->coord);
                    }
                    continue;
                }
                auto holds = [&](uint64_t coord) {
                    return std::binary_search(first, last, LabelQuery::HeaderCoord{ l, coord });
                };
                forwards.erase(std::remove_if(forwards.begin(), forwards.end(),
                                              [&](uint64_t o) { return !holds(o + j); }),
                               forwards.end());
                backwards.erase(std::remove_if(backwards.begin(), backwards.end(),
                                               [&](uint64_t o) { return o < j || !holds(o - j); }),
                                backwards.end());
            }
        }
    }
    // a backward origin is the place's highest coordinate (every one left is >= span)
    const uint64_t span = keys.size() - 1;
    for (size_t l = 0; l < refs.size(); ++l) {
        const std::vector<uint64_t> &forwards = origins[l][0];
        const std::vector<uint64_t> &backwards = origins[l][1];
        if (!forwards.empty()
                && (backwards.empty() || forwards.front() <= backwards.front() - span)) {
            (*places)[l] = Place{ forwards.front(), 1 };
        } else if (!backwards.empty()) {
            (*places)[l] = Place{ backwards.front(), -1 };
        }
    }
}

// The permitted set of a seed without named labels. COLUMN labels: the columns on every
// seed k-mer. HEADER labels: the records holding the seed's k-mers at consecutive
// coordinates (k-mer j at origin + dir * j, see Place), each at its place with the lowest
// coordinate, as scan_named_labels places a named record.
void Walker::derive_seed_labels(const std::vector<node_index> &keys,
                                std::vector<LabelRef> *refs) {
    using Key = std::pair<Column, uint64_t>;
    // a candidate place of the seed in a column's coordinates
    struct Spot {
        Column column;
        uint64_t origin;
        int8_t dir;
        // ascending column and origin
        bool operator<(const Spot &o) const {
            return std::make_tuple(column, origin, dir)
                    < std::make_tuple(o.column, o.origin, o.dir);
        }
    };

    const LabelKind kind = strategy_.seed_label_kind.value_or(
            oracle_.coord_to_header() ? LabelKind::HEADER : LabelKind::COLUMN);
    if (kind == LabelKind::HEADER && !oracle_.coord_to_header()) {
        throw std::invalid_argument("Deriving sequence header labels from the seed requires a "
                                    "CoordToHeader index; set seed_label_kind to \"column\"");
    }
    // a seed of exactly k bases has no intersection to narrow the candidates: bound what
    // the first row may materialise
    const size_t max_candidates = std::max<size_t>(
            size_t(1) << 16, 64 * std::min<size_t>(strategy_.max_seed_labels, size_t(1) << 40));

    // the labels carrying every k-mer read so far, ascending
    std::vector<Key> live, next, cur;
    // HEADER: the places consistent with every k-mer read so far, ascending
    std::vector<Spot> spots, next_spots;
    std::vector<uint64_t> sorted;
    size_t read = 0;
    for (size_t begin = 0; begin < keys.size(); begin += kDeriveWindow) {
        const size_t end = std::min(keys.size(), begin + kDeriveWindow);
        std::vector<Row> distinct;
        for (size_t i = begin; i < end; ++i) {
            distinct.push_back(AnnotatedDBG::graph_to_anno_index(keys[i]));
        }
        std::sort(distinct.begin(), distinct.end());
        distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
        std::vector<annot::matrix::BinaryMatrix::SetBitPositions> plain;
        std::vector<annot::matrix::MultiIntMatrix::RowTuples> tuples;
        if (kind == LabelKind::HEADER) {
            tuples = oracle_.get_row_tuples(distinct);
        } else {
            plain = oracle_.get_rows(distinct);
        }
        auto row_of = [&](size_t i) {
            return static_cast<size_t>(
                    std::lower_bound(distinct.begin(), distinct.end(),
                                     AnnotatedDBG::graph_to_anno_index(keys[i]))
                    - distinct.begin());
        };
        // the size of a row: its columns, and for a tuple row its coordinates
        auto cost_of = [&](size_t i) {
            const size_t r = row_of(i);
            if (kind == LabelKind::COLUMN)
                return plain[r].size();
            size_t n = 0;
            for (const auto &entry : tuples[r]) {
                n += 1 + entry.second.size();
            }
            return n;
        };

        std::vector<size_t> order;
        for (size_t i = begin; i < end; ++i) {
            order.push_back(i);
        }
        if (!read) {
            // the narrowest k-mer of the first window starts the intersection
            std::vector<size_t> cost;
            for (size_t i : order) {
                cost.push_back(cost_of(i));
            }
            const size_t best = std::min_element(cost.begin(), cost.end()) - cost.begin();
            std::iter_swap(order.begin(), order.begin() + best);
            if (cost[best] > max_candidates) {
                throw SeedDerivationError(SeedDerivationError::TOO_WIDE,
                        "The permitted set cannot be derived from this seed: the narrowest of "
                        "its first " + std::to_string(end - begin) + " k-mers alone has "
                        + std::to_string(cost[best]) + " annotation entries, over the limit of "
                        + std::to_string(max_candidates) + ". Start the seed in a less "
                          "repetitive k-mer, or name the labels explicitly.",
                        static_cast<double>(max_candidates), static_cast<double>(cost[best]));
            }
        }

        for (size_t i : order) {
            const size_t r = row_of(i);
            if (kind == LabelKind::COLUMN) {
                cur.clear();
                for (Column c : plain[r]) {
                    // after the first k-mer, keep only the live columns
                    if (!read || std::binary_search(live.begin(), live.end(), Key(c, 0)))
                        cur.emplace_back(c, 0);
                }
                std::sort(cur.begin(), cur.end());
                cur.erase(std::unique(cur.begin(), cur.end()), cur.end());
                if (!read) {
                    live.swap(cur);
                } else {
                    next.clear();
                    std::set_intersection(live.begin(), live.end(), cur.begin(), cur.end(),
                                          std::back_inserter(next));
                    live.swap(next);
                }
            } else if (!read) {
                // every coordinate of the first k-mer read starts a place, in both
                // directions in the canonical regimes
                for (const auto &[c, coords] : tuples[r]) {
                    if (!oracle_.has_header_records(c))
                        continue;  // no records, so no header labels
                    for (uint64_t coord : coords) {
                        if (coord >= i)
                            spots.push_back(Spot{ c, coord - i, 1 });
                        if (regime_ != Regime::BASIC)
                            spots.push_back(Spot{ c, coord + i, -1 });
                    }
                }
                std::sort(spots.begin(), spots.end());
            } else {
                // keep the places whose coordinate for this k-mer its row holds, so that a
                // repetitive k-mer cannot grow the candidates
                next_spots.clear();
                for (const auto &[c, coords] : tuples[r]) {
                    auto spot = std::lower_bound(spots.begin(), spots.end(), c,
                                                 [](const Spot &s, Column column) {
                                                     return s.column < column;
                                                 });
                    if (spot == spots.end() || spot->column != c)
                        continue;
                    const uint64_t *coords_begin = coords.data();
                    const uint64_t *coords_end = coords_begin + coords.size();
                    if (!std::is_sorted(coords_begin, coords_end)) {
                        sorted.assign(coords_begin, coords_end);
                        std::sort(sorted.begin(), sorted.end());
                        coords_begin = sorted.data();
                        coords_end = coords_begin + sorted.size();
                    }
                    for ( ; spot != spots.end() && spot->column == c; ++spot) {
                        if (spot->dir < 0 && spot->origin < i)
                            continue;
                        const uint64_t coord = spot->dir > 0 ? spot->origin + i
                                                             : spot->origin - i;
                        if (std::binary_search(coords_begin, coords_end, coord))
                            next_spots.push_back(*spot);
                    }
                }
                std::sort(next_spots.begin(), next_spots.end());
                spots.swap(next_spots);
            }
            ++read;
            if (kind == LabelKind::COLUMN ? live.empty() : spots.empty()) {
                throw SeedDerivationError(SeedDerivationError::NO_CARRIER,
                                          "No label supports every k-mer of the seed (the "
                                          "permitted set derived from the seed is empty at "
                                          "k-mer " + std::to_string(i) + " of "
                                          + std::to_string(keys.size()) + ")",
                                          static_cast<double>(keys.size()),
                                          static_cast<double>(read));
            }
            if (read < keys.size() && strategy_.time_budget_ms > 0
                    && elapsed_ms() >= strategy_.time_budget_ms) {
                throw SeedDerivationError(SeedDerivationError::TIME_BUDGET,
                        "The time budget (" + std::to_string(strategy_.time_budget_ms)
                        + " ms) ran out deriving the permitted set from the seed, after "
                        + std::to_string(read) + " of its " + std::to_string(keys.size())
                        + " k-mers; raise bounds.time_budget_ms or name the labels explicitly",
                        strategy_.time_budget_ms, static_cast<double>(read));
            }
        }
    }

    std::vector<Place> live_places;
    if (kind == LabelKind::HEADER) {
        // A place counts when it lies in one record. In the order of their lowest
        // coordinate (the origin forwards, the last k-mer's coordinate backwards; every
        // backward origin left is >= span) a record's places are together: its first is
        // kept, forwards first when equal.
        const uint64_t span = keys.size() - 1;
        std::sort(spots.begin(), spots.end(), [span](const Spot &a, const Spot &b) {
            return std::make_tuple(a.column, a.dir > 0 ? a.origin : a.origin - span, -a.dir)
                    < std::make_tuple(b.column, b.dir > 0 ? b.origin : b.origin - span, -b.dir);
        });
        for (const Spot &spot : spots) {
            const uint64_t other = spot.dir > 0 ? spot.origin + span : spot.origin - span;
            const uint64_t seq_id = oracle_.map_coord(spot.column, spot.origin).first;
            if (oracle_.map_coord(spot.column, other).first != seq_id)
                continue;
            const Key key(spot.column, seq_id);
            if (!live.empty() && live.back() == key)
                continue;
            live.push_back(key);
            live_places.push_back(Place{ spot.origin, spot.dir });
        }
        if (live.empty()) {
            throw SeedDerivationError(SeedDerivationError::NO_CARRIER,
                                      "No sequence record holds the k-mers of the seed at "
                                      "consecutive coordinates (the permitted set derived "
                                      "from the seed is empty)",
                                      static_cast<double>(keys.size()),
                                      static_cast<double>(read));
        }
    }

    result_.labels_from_seed = true;
    result_.labels_supporting_total = live.size();
    auto name_of = [&](const Key &key) -> const std::string& {
        return kind == LabelKind::COLUMN ? oracle_.column_name(key.first)
                                         : oracle_.header_name(key.first, key.second);
    };
    if (kind == LabelKind::HEADER) {
        // A derived list must be resubmittable verbatim as an explicit one, and
        // resolve_label refuses a header that also names a column or another record. Only
        // the labels kept under max_seed_labels are returned, so only they are checked.
        const size_t kept = std::min(live.size(), strategy_.max_seed_labels);
        for (size_t i = 0; i < kept; ++i) {
            const Key &key = live[i];
            const std::string &name = name_of(key);
            const bool column_too = oracle_.find_column(name).has_value();
            if (column_too || oracle_.coord_to_header()->is_repeated_header(name)) {
                throw SeedDerivationError(SeedDerivationError::AMBIGUOUS_HEADER,
                        "The labels derived from the seed are not resubmittable: the sequence "
                        "header '" + name + "' "
                        + (column_too ? "is also the name of an annotation column"
                                      : "is also the header of another sequence record")
                        + ". Set seed_label_kind to \"column\" or name the labels explicitly.",
                        0, 0, name);
            }
        }
    }
    if (live.size() > strategy_.max_seed_labels) {
        result_.labels_dropped = live.size() - strategy_.max_seed_labels;
        live.resize(strategy_.max_seed_labels);
    }
    for (size_t i = 0; i < live.size(); ++i) {
        const Key &key = live[i];
        LabelRef ref;
        ref.kind = kind;
        ref.column = key.first;
        ref.seq_id = kind == LabelKind::HEADER ? key.second : 0;
        ref.name = name_of(key);
        refs->push_back(std::move(ref));
        places_.push_back(kind == LabelKind::HEADER ? live_places[i] : Place());
    }
}

void Walker::init_arm(ArmState &arm) {
    arm.result.arm = arm.arm;
    arm.result.requested = strategy_.direction == Strategy::BOTH
        || (strategy_.direction == Strategy::LEFT) == (arm.arm == Arm::LEFT);
    if (!arm.result.requested)
        return;

    Item root;
    root.node = arm.arm == Arm::RIGHT ? nodes_.back() : nodes_.front();
    root.kmer = arm.arm == Arm::RIGHT ? seed_upper_.substr(seed_upper_.size() - k_)
                                      : seed_upper_.substr(0, k_);
    const uint64_t span = nodes_.size() - 1;
    for (LabelId l = 0; l < result_.label_dict.size(); ++l) {
        Entry e;
        e.label = l;
        if (result_.label_dict[l].kind == LabelKind::HEADER) {
            // the right arm starts at the seed's last k-mer and follows the record in its
            // direction, the left arm at the first k-mer and against it
            const Place &place = places_[l];
            assert(place.dir);
            if (arm.arm == Arm::RIGHT) {
                e.coord = place.dir > 0 ? place.origin + span : place.origin - span;
                e.step = place.dir;
            } else {
                e.coord = place.origin;
                e.step = static_cast<int8_t>(-place.dir);
            }
        }
        root.state.push_back(e);
    }
    root.segment = new_segment(arm, {}, 0, labels_of(root.state));
    arm.frontier.push_back(std::move(root));
}


/******************************* graph access *******************************/

std::vector<Succ> Walker::successors(Arm arm, node_index node, const std::string &kmer) const {
    std::vector<Succ> out;
    auto cb = [&](node_index n, char c) {
        if (c != kSentinel && n != npos)
            out.push_back({ n, c });
    };
    try {
        if (arm == Arm::RIGHT) {
            if (canonical_) {
                canonical_->call_outgoing_kmers(node, kmer, cb);
            } else {
                graph_.call_outgoing_kmers(node, cb);
            }
        } else {
            if (canonical_) {
                canonical_->call_incoming_kmers(node, kmer, cb);
            } else {
                graph_.call_incoming_kmers(node, cb);
            }
        }
    } catch (const std::runtime_error &e) {
        // CanonicalDBG throws an empty runtime_error on an inconsistent primary graph
        throw std::runtime_error("Inconsistent primary graph at node " + std::to_string(node)
                                 + " (" + kmer + ") on arm " + to_string(arm)
                                 + (e.what()[0] ? std::string(": ") + e.what() : std::string()));
    }
    std::sort(out.begin(), out.end(), [](const Succ &a, const Succ &b) { return a.ch < b.ch; });
    return out;
}

// the k-mer of the successor reached by |c| (|out| must not alias |kmer|)
void Walker::succ_kmer(Arm arm, const std::string &kmer, char c, std::string *out) const {
    out->clear();
    if (arm == Arm::RIGHT) {
        out->append(kmer, 1, k_ - 1);
        out->push_back(c);
    } else {
        out->push_back(c);
        out->append(kmer, 0, k_ - 1);
    }
}

// the (k+1)-mer of a step in natural orientation
void Walker::step_kmer(Arm arm, const std::string &kmer, char c, std::string *out) const {
    out->clear();
    if (arm == Arm::RIGHT) {
        out->append(kmer);
        out->push_back(c);
    } else {
        out->push_back(c);
        out->append(kmer);
    }
}

// A self-reverse-complementary step: the (k+1)-mer is its own reverse complement, or for
// even k either k-mer is
bool Walker::is_hairpin(const std::string &kmer_u, const std::string &kmer_v,
                        const std::string &step) const {
    auto palindromic = [](const std::string &s) {
        std::string rc = s;
        ::reverse_complement(rc);
        return rc == s;
    };
    return palindromic(step) || (k_ % 2 == 0 && (palindromic(kmer_v) || palindromic(kmer_u)));
}

// the edge of a step, canonical in the canonical regimes, and whether it was flipped
std::pair<EdgeKey, bool> Walker::edge_key(const std::string &step) const {
    std::string rc;
    std::string_view s = step;
    bool flip = false;
    if (regime_ != Regime::BASIC) {
        rc = step;
        ::reverse_complement(rc);
        if (rc < step) {
            s = rc;
            flip = true;
        }
    }
    EdgeKey key { 0, 0 };
    if (packable_) {
        for (char c : s) {
            key.lo = (key.lo << bits_) | code_[static_cast<unsigned char>(c)];
        }
    } else {
        key.lo = folly::hash::fnva64_buf(s.data(), s.size());
        key.hi = folly::hash::fnva64_buf(s.data(), s.size(), 0x84222325cbf29ce4ULL);
    }
    return { key, flip };
}


/******************************** bookkeeping *******************************/

size_t Walker::new_segment(ArmState &arm, std::vector<size_t> parents, uint64_t from_bp,
                           std::vector<LabelId> labels_start) {
    Segment seg;
    seg.id = arm.result.segments.size();
    seg.parents = std::move(parents);
    seg.from_bp = from_bp;
    seg.labels_start = std::move(labels_start);
    arm.result.segments.push_back(std::move(seg));
    arm.walk_seq.emplace_back();
    arm.leaves.emplace_back();
    return arm.result.segments.size() - 1;
}

// Mark the ancestors of |segment| (its chain of parents). The parents of a segment never
// change, so the marks stay valid until another segment is marked.
void Walker::mark_ancestors(ArmState &arm, size_t segment) {
    if (arm.marked_segment == segment)
        return;
    const auto &segments = arm.result.segments;
    if (arm.visit_mark.size() < segments.size())
        arm.visit_mark.resize(segments.size(), 0);
    ++arm.visit_epoch;
    arm.ancestors.clear();
    for (size_t s = segment; !segments[s].parents.empty(); ) {
        s = segments[s].parents[0];
        arm.visit_mark[s] = arm.visit_epoch;
        arm.ancestors.push_back(s);
    }
    arm.marked_segment = segment;
}

// requires mark_ancestors(arm, segment)
bool Walker::is_ancestor_or_self(const ArmState &arm, size_t anc, size_t segment) const {
    assert(arm.marked_segment == segment);
    // segments are numbered in creation order
    return anc == segment || (anc < segment && arm.visit_mark[anc] == arm.visit_epoch);
}

void Walker::record_edge(ArmState &arm, size_t segment, const Cand &c) {
    arm.used_edges[c.edge].push_back(EdgeUse{ segment, c.rc });
    arm.used_by_segment[EdgeSegKey{ c.edge, segment }] |= static_cast<uint8_t>(1u << c.rc);
}

void Walker::end_label(ArmState &arm, Item &item, Entry &e, EndReason reason,
                       const char *text) {
    e.reason = reason;
    arm.result.segments[item.segment].label_ends.push_back(
            LabelEnd{ item.ext_bp, e.label, reason, text });
}

// |item| is a leaf: every label on it has ended
void Walker::finish_path(ArmState &arm, Item &item, std::optional<EndReason> path_reason) {
    arm.result.segments[item.segment].labels_end = labels_of(item.state);
    LeafInfo &leaf = arm.leaves[item.segment];
    leaf.is_leaf = true;
    leaf.path_reason = path_reason;
    for (const Entry &e : item.state) {
        leaf.end_labels.push_back(PathEnd{ e.label, e.reason, e.branches });
    }
    arm.finished_leaves++;
}

void Walker::censor_item(ArmState &arm, Item &item, EndReason reason) {
    for (Entry &e : item.state) {
        end_label(arm, item, e, reason, "");
    }
    finish_path(arm, item, reason);
}


/********************************** the walk ********************************/

void Walker::check_structure(ArmState &arm, const Item &item, Cand &c) {
    std::string step;
    step_kmer(arm.arm, item.kmer, c.succ->ch, &step);
    if (regime_ != Regime::BASIC && is_hairpin(item.kmer, c.kmer, step)) {
        c.skipped = true;
        return;
    }
    if (seed_nodes_.count(c.succ->node)) {
        c.blocked = true;
        c.block_reason = EndReason::REACHED_SEED;
        return;
    }
    std::tie(c.edge, c.rc) = edge_key(step);
    auto it = arm.used_edges.find(c.edge);
    if (it == arm.used_edges.end())
        return;
    // Was the edge taken by this path (the head's segment or an ancestor), and in which
    // orientation? Either scan the edge's uses or probe (edge, segment) for the path's own
    // segments, whichever are fewer: the same uses, so the same answer.
    mark_ancestors(arm, item.segment);
    const std::vector<EdgeUse> &uses = it->second;
    bool same = false, opposite = false;
    if (uses.size() <= arm.ancestors.size() + 1) {
        for (const EdgeUse &u : uses) {
            if (is_ancestor_or_self(arm, u.segment, item.segment))
                (u.rc == c.rc ? same : opposite) = true;
        }
    } else {
        auto probe = [&](size_t segment) {
            auto jt = arm.used_by_segment.find(EdgeSegKey{ c.edge, segment });
            if (jt == arm.used_by_segment.end())
                return;
            same |= (jt->second >> c.rc) & 1;
            opposite |= (jt->second >> !c.rc) & 1;
        };
        probe(item.segment);
        for (size_t a : arm.ancestors) {
            probe(a);
        }
    }
    if (opposite) {
        c.blocked = true;
        c.block_reason = EndReason::EDGE_REUSE_RC;
    } else if (same) {
        c.blocked = true;
        c.block_reason = EndReason::EDGE_REUSE;
    }
}

// The caps, decided before anything of a head is committed: |nf| successors would be
// followed, |remaining| heads of the level are still to come after this one
std::optional<EndReason> Walker::cap_check(const ArmState &arm, size_t nf, size_t remaining) {
    if (!nf)
        return std::nullopt;
    auto over = [&](uint64_t demand, EndReason reason) {
        cap_demand_ = static_cast<double>(demand);
        return reason;
    };
    if (steps_total_ + nf > strategy_.max_steps)
        return over(steps_total_ + nf, EndReason::MAX_STEPS);
    const size_t live_after = remaining + arm.next.size() + nf;
    if (nf >= 2 && arm.finished_leaves + live_after > strategy_.max_paths)
        return over(arm.finished_leaves + live_after, EndReason::MAX_PATHS);
    if (live_after > strategy_.max_live_paths)
        return over(live_after, EndReason::MAX_LIVE_PATHS);
    return std::nullopt;
}

std::optional<EndReason> Walker::process_item(ArmState &arm, Item &item,
                                              const std::vector<Succ> &succs,
                                              const LabelQuery::Hits *const *hits,
                                              size_t remaining_in_level) {
    const uint64_t at = item.ext_bp;
    if (succs.empty()) {
        for (Entry &e : item.state) {
            end_label(arm, item, e, EndReason::DEAD_END, "");
        }
        finish_path(arm, item, std::nullopt);
        return std::nullopt;
    }

    // ---- the head's labels on each successor, and the walk rule
    std::vector<Cand> cands(succs.size());
    for (size_t i = 0; i < succs.size(); ++i) {
        Cand &c = cands[i];
        c.succ = &succs[i];
        succ_kmer(arm.arm, item.kmer, succs[i].ch, &c.kmer);
        // all sorted by label
        const LabelQuery::Hits &hit = *hits[i];
        auto h = hit.labels.begin();
        auto hc = hit.coords.begin();
        for (const Entry &e : item.state) {
            if (!e.step) {
                // a column label: on the successor's row
                while (h != hit.labels.end() && *h < e.label) {
                    ++h;
                }
                if (h != hit.labels.end() && *h == e.label) {
                    c.state.push_back(e);
                    c.ruled = true;
                }
                continue;
            }
            // a header label: the successor holding the record's next coordinate
            if (at_record_end(e))
                continue;
            const LabelQuery::HeaderCoord next { e.label, e.step > 0 ? e.coord + 1
                                                                      : e.coord - 1 };
            // usually the entry after the one the previous label matched
            if (hc != hit.coords.end() && *hc < next && ++hc != hit.coords.end() && *hc < next)
                hc = std::lower_bound(hc, hit.coords.end(), next);
            if (hc != hit.coords.end() && *hc == next) {
                c.state.push_back(e);
                c.state.back().coord = next.coord;
            }
        }
        // The walk rule binds column labels only: where it forbids the successor, they
        // are stopped and the header labels, finite by coordinate order, go on
        if (c.ruled) {
            check_structure(arm, item, c);
            if (!c.admissible()) {
                auto columns = std::stable_partition(c.state.begin(), c.state.end(),
                                                     [](const Entry &e) { return e.step != 0; });
                c.barred.assign(columns, c.state.end());
                c.state.erase(columns, c.state.end());
            }
        }
    }

    // ---- per-lineage branching: a label on several admissible successors continues on
    // all of them within its branch allowance, else it ends on all of them (BRANCH)
    tsl::hopscotch_map<LabelId, uint32_t> continuations;
    for (const Cand &c : cands) {
        for (const Entry &e : c.state) {
            continuations[e.label]++;
        }
    }
    tsl::hopscotch_set<LabelId> excluded;
    for (const Entry &src : item.state) {
        auto it = continuations.find(src.label);
        if (it != continuations.end() && it->second >= 2
                && strategy_.max_label_branches != Strategy::kUnlimited
                && src.branches + 1 > strategy_.max_label_branches) {
            excluded.insert(src.label);
        }
    }
    bool ambiguous = false;
    std::vector<Cand*> followed;
    for (Cand &c : cands) {
        if (!excluded.empty()) {
            c.state.erase(std::remove_if(c.state.begin(), c.state.end(), [&](const Entry &e) {
                return excluded.count(e.label);
            }), c.state.end());
        }
        for (Entry &e : c.state) {
            if (continuations[e.label] >= 2) {
                e.branches++;
                ambiguous = true;
            }
        }
        if (!c.state.empty()) {
            c.followed = true;
            followed.push_back(&c);
        }
    }
    // The split limit is a safety cap, not a label rule: the walks past this fork are not
    // explored, so the arm is complete only up to here
    bool split_capped = false;
    if (followed.size() >= 2 && item.splits + 1 > strategy_.max_splits_per_path) {
        for (Cand *c : followed) {
            c->followed = false;
            c->split_limited = true;
        }
        followed.clear();
        split_capped = true;
        if (!arm.result.cap_trigger) {
            arm.result.cap_trigger = CapTrigger{ EndReason::MAX_SPLITS, at,
                                                 static_cast<double>(item.splits + 1) };
        }
        arm.boundary = std::min(arm.boundary, at);
        arm.result.status = ArmResult::TRUNCATED;
    }
    const size_t nf = followed.size();

    // ---- the caps, before anything is committed
    if (auto cap = cap_check(arm, nf, remaining_in_level))
        return cap;

    // ---- the labels that end here, in state order
    for (Entry &src : item.state) {
        const LabelId l = src.label;
        bool continues = false, split_limited = false, hairpin = false;
        uint8_t blocked = 0;
        auto has = [&](const State &state) {
            return std::binary_search(state.begin(), state.end(), src,
                    [](const Entry &a, const Entry &b) { return a.label < b.label; });
        };
        for (const Cand &c : cands) {
            if (has(c.state)) {
                continues |= c.followed;
                split_limited |= c.split_limited;
            } else if (has(c.barred)) {
                hairpin |= c.skipped;
                if (c.blocked)
                    blocked = std::max(blocked, block_rank(c.block_reason));
            }
        }
        if (continues)
            continue;
        if (excluded.count(l)) {
            end_label(arm, item, src, EndReason::BRANCH, "");
        } else if (split_limited) {
            end_label(arm, item, src, EndReason::MAX_SPLITS, "");
        } else if (blocked) {
            end_label(arm, item, src, blocked == 3 ? EndReason::REACHED_SEED
                                     : blocked == 2 ? EndReason::EDGE_REUSE_RC
                                                    : EndReason::EDGE_REUSE, "");
        } else if (hairpin) {
            end_label(arm, item, src, EndReason::DEAD_END, "hairpin");
        } else if (src.step && at_record_end(src)) {
            end_label(arm, item, src, EndReason::RECORD_END, "");
        } else {
            end_label(arm, item, src, EndReason::LABEL_LOST, "");
        }
    }

    if (!nf) {
        finish_path(arm, item, split_capped ? std::optional<EndReason>(EndReason::MAX_SPLITS)
                                            : std::nullopt);
        return std::nullopt;
    }

    // ---- steps
    arm.result.steps += nf;
    steps_total_ += nf;

    if (nf == 1) {
        Cand &c = *followed[0];
        arm.walk_seq[item.segment].push_back(c.succ->ch);
        arm.result.segments[item.segment].length_bp++;
        // only column labels look edges up, and theirs is an admissible step
        if (c.ruled && c.admissible())
            record_edge(arm, item.segment, c);
        item.node = c.succ->node;
        item.kmer = std::move(c.kmer);
        item.ext_bp = at + 1;
        item.state = std::move(c.state);
        arm.next.push_back(std::move(item));
        return std::nullopt;
    }

    // ---- split
    Split split { at, item.segment, {}, ambiguous };
    arm.result.segments[item.segment].labels_end = labels_of(item.state);
    for (Cand *c : followed) {
        size_t child = new_segment(arm, { item.segment }, at, labels_of(c->state));
        arm.result.segments[item.segment].children.push_back(child);
        split.children.push_back(child);
        arm.walk_seq[child].push_back(c->succ->ch);
        arm.result.segments[child].length_bp = 1;
        if (c->ruled && c->admissible())
            record_edge(arm, child, *c);
        Item ni;
        ni.segment = child;
        ni.node = c->succ->node;
        ni.kmer = std::move(c->kmer);
        ni.ext_bp = at + 1;
        ni.state = std::move(c->state);
        ni.splits = item.splits + 1;
        ni.path_id = arm.next_path_id++;
        arm.next.push_back(std::move(ni));
    }
    arm.result.splits.push_back(std::move(split));
    return std::nullopt;
}

// Stop the arm at a cap: the heads from items[from] on and every head of the next level are
// censored; a head that already reached the extension length ends with MAX_EXTENSION
void Walker::stop_arm(ArmState &arm, EndReason reason, std::vector<Item> *items, size_t from) {
    const size_t heads = arm.next.size() + (items ? items->size() - from : 0);
    if (!heads)
        return;
    const uint64_t radius = strategy_.max_extension_bp;
    size_t live = 0;
    for (const Item &item : arm.next) {
        live += item.ext_bp < radius;
    }
    if (items) {
        for (size_t j = from; j < items->size(); ++j) {
            live += (*items)[j].ext_bp < radius;
        }
    }
    const Item &first = items && from < items->size() ? (*items)[from] : arm.next.front();
    if (!arm.result.cap_trigger)
        arm.result.cap_trigger = CapTrigger{ reason, first.ext_bp, cap_demand_ };
    arm.boundary = std::min(arm.boundary, first.ext_bp);
    arm.result.frontier_live_paths = live;
    arm.result.status = ArmResult::TRUNCATED;
    auto censor = [&](Item &item) {
        censor_item(arm, item, item.ext_bp >= radius ? EndReason::MAX_EXTENSION : reason);
    };
    if (items) {
        for (size_t j = from; j < items->size(); ++j) {
            censor((*items)[j]);
        }
    }
    for (Item &item : arm.next) {
        censor(item);
    }
    arm.next.clear();
    arm.frontier.clear();
    arm.stopped = true;
}

void Walker::stop_frontier(ArmState &arm, EndReason reason) {
    if (arm.frontier.empty())
        return;
    std::vector<Item> rest;
    for (Item &item : arm.frontier) {
        if (item.ext_bp >= strategy_.max_extension_bp) {
            censor_item(arm, item, EndReason::MAX_EXTENSION);
        } else {
            rest.push_back(std::move(item));
        }
    }
    arm.frontier.swap(rest);
    if (!arm.frontier.empty())
        stop_arm(arm, reason, &arm.frontier, 0);
}

// A cap at items[index]: the steps and the time stop the whole seed, the others the arm
void Walker::trip(ArmState &arm, EndReason reason, std::vector<Item> &items, size_t index) {
    stop_arm(arm, reason, &items, index);
    if (reason == EndReason::MAX_STEPS || reason == EndReason::TIME_BUDGET) {
        seed_stopped_ = true;
        for (ArmState &other : arms_) {
            if (&other != &arm)
                stop_frontier(other, reason);
        }
    }
}

void Walker::run_level(ArmState &arm, uint64_t depth) {
    std::vector<Item> items;
    items.swap(arm.frontier);
    std::sort(items.begin(), items.end(),
              [](const Item &a, const Item &b) { return a.path_id < b.path_id; });

    // read the deadline before the level's fetch too, which can be long on a heavy
    // annotation; it ends the heads as the check before each head below would
    if (depth > 0 && time_exceeded()) {
        for (size_t i = 0; i < items.size(); ++i) {
            if (items[i].ext_bp < strategy_.max_extension_bp) {
                cap_demand_ = elapsed_ms();
                trip(arm, EndReason::TIME_BUDGET, items, i);
                return;
            }
            censor_item(arm, items[i], EndReason::MAX_EXTENSION);
        }
        return;
    }

    // the successors of every head below the extension length, and their labels in one fetch
    std::vector<std::vector<Succ>> succs(items.size());
    std::vector<node_index> keys;
    std::string kmer;
    for (size_t i = 0; i < items.size(); ++i) {
        if (items[i].ext_bp >= strategy_.max_extension_bp)
            continue;
        succs[i] = successors(arm.arm, items[i].node, items[i].kmer);
        for (const Succ &s : succs[i]) {
            // the spelling is only needed by native canonical graphs
            if (regime_ == Regime::CANONICAL)
                succ_kmer(arm.arm, items[i].kmer, s.ch, &kmer);
            keys.push_back(key_mapper_.key(s.node, kmer));
        }
    }
    const std::vector<const LabelQuery::Hits*> hits = query_->fetch_hits(keys);

    size_t offset = 0;
    for (size_t i = 0; i < items.size(); ++i) {
        Item &item = items[i];
        if (item.ext_bp >= strategy_.max_extension_bp) {
            censor_item(arm, item, EndReason::MAX_EXTENSION);
            continue;
        }
        // the deadline is read before every head; never at depth 0, so that a time budget
        // of 0 still delivers the result complete to the first level
        if (depth > 0 && time_exceeded()) {
            cap_demand_ = elapsed_ms();
            trip(arm, EndReason::TIME_BUDGET, items, i);
            return;
        }
        auto cap = process_item(arm, item, succs[i], hits.data() + offset, items.size() - i - 1);
        offset += succs[i].size();
        if (cap) {
            trip(arm, *cap, items, i);
            return;
        }
    }
    arm.frontier.swap(arm.next);
}

void Walker::finalize(ArmState &arm) {
    ArmResult &res = arm.result;
    // every walk up to the first unexpanded head's depth is present; an arm not requested
    // has no walk
    if (res.requested)
        res.complete_to_bp = std::min(arm.boundary, strategy_.max_extension_bp);
    for (size_t s = 0; s < res.segments.size(); ++s) {
        Segment &seg = res.segments[s];
        seg.sequence = std::move(arm.walk_seq[s]);
        if (arm.arm == Arm::LEFT)
            std::reverse(seg.sequence.begin(), seg.sequence.end());
        std::stable_sort(seg.label_ends.begin(), seg.label_ends.end(),
                         [](const LabelEnd &a, const LabelEnd &b) { return a.at_bp < b.at_bp; });
    }
    // one path per leaf, in segment order
    for (size_t s = 0; s < res.segments.size(); ++s) {
        LeafInfo &leaf = arm.leaves[s];
        if (!leaf.is_leaf)
            continue;
        PathResult path;
        path.id = res.paths.size();
        path.leaf = s;
        path.length_bp = res.segments[s].from_bp + res.segments[s].length_bp;
        path.path_reason = leaf.path_reason;
        path.end_labels = std::move(leaf.end_labels);
        res.paths.push_back(std::move(path));
    }
}

SeedResult Walker::run() {
    if (!strategy_.max_seed_labels)
        throw std::invalid_argument("max_seed_labels must be positive");
    if (!std::isfinite(strategy_.time_budget_ms) || strategy_.time_budget_ms < 0)
        throw std::invalid_argument("time_budget_ms must be finite and not negative");
    if (!strategy_.max_steps || !strategy_.max_live_paths || !strategy_.max_paths
            || !strategy_.max_extension_bp) {
        throw std::invalid_argument("The step, path and extension caps must be positive");
    }

    validate_seed();
    query_ = std::make_unique<LabelQuery>(oracle_, result_.label_dict);

    arms_[0].arm = Arm::LEFT;
    arms_[1].arm = Arm::RIGHT;
    for (ArmState &arm : arms_) {
        init_arm(arm);
    }

    for (uint64_t depth = 0; !seed_stopped_; ++depth) {
        if (std::all_of(arms_.begin(), arms_.end(),
                        [](const ArmState &arm) { return arm.frontier.empty(); }))
            break;
        if (depth > 0 && time_exceeded()) {
            // only heads below the extension length are censored by the deadline: when every
            // remaining head reached it, the walk is complete and nothing was stopped
            cap_demand_ = elapsed_ms();
            for (ArmState &arm : arms_) {
                stop_frontier(arm, EndReason::TIME_BUDGET);
            }
            break;
        }
        for (ArmState &arm : arms_) {
            if (!arm.frontier.empty() && !arm.stopped)
                run_level(arm, depth);
            if (seed_stopped_)
                break;
        }
    }

    for (ArmState &arm : arms_) {
        finalize(arm);
    }
    result_.arms[static_cast<size_t>(Arm::LEFT)] = std::move(arms_[0].result);
    result_.arms[static_cast<size_t>(Arm::RIGHT)] = std::move(arms_[1].result);
    result_.elapsed_seconds = elapsed_ms() / 1000;
    return std::move(result_);
}

} // namespace


void check_seed(const LabelOracle &oracle, const Seed &seed) {
    const std::string &seq = seed.sequence;
    if (seq.size() < oracle.get_k())
        throw std::invalid_argument("Seed shorter than k");
    const std::string upper = in_graph_case(seq);
    const std::string &alphabet = oracle.graph().alphabet();
    for (size_t i = 0; i < upper.size(); ++i) {
        if (upper[i] == kSentinel || alphabet.find(upper[i]) == std::string::npos) {
            throw std::invalid_argument("Invalid character '" + std::string(1, seq[i])
                                        + "' at position " + std::to_string(i) + " of the seed");
        }
    }
    tsl::hopscotch_set<std::string_view> named;
    for (const auto &name : seed.labels) {
        if (!named.insert(name).second)
            throw std::invalid_argument("Duplicate seed label '" + name + "'");
        oracle.resolve_label(name);
    }
}

SeedResult traverse_seed(const LabelOracle &oracle, const Seed &seed, const Strategy &strategy) {
    return Walker(oracle, seed, strategy).run();
}

} // namespace traversal
} // namespace graph
} // namespace mtg
