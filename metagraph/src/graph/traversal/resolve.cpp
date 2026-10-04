#include "resolve.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <numeric>

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "graph/annotated_dbg.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/logger.hpp"


namespace mtg {
namespace graph {
namespace traversal {

using mtg::common::logger;

uint64_t fnv1a64(std::string_view data, uint64_t hash) {
    for (unsigned char c : data) {
        hash ^= c;
        hash *= 0x100000001b3ULL;
    }
    return hash;
}

std::string hex64(uint64_t x) {
    static const char *digits = "0123456789abcdef";
    std::string s(16, '0');
    for (int i = 15; i >= 0; --i) {
        s[i] = digits[x & 15];
        x >>= 4;
    }
    return s;
}

std::string make_seed_id(const std::string &release_id,
                         std::string_view sequence,
                         bool canonical_orientation,
                         std::vector<std::string> labels) {
    // defined over the sequence after the build's case mapping (§4.3 / §6.1), so a
    // seed resubmitted in another case keeps its id
    std::string seq(sequence);
#if ! _DNA_CASE_SENSITIVE_GRAPH
    for (char &c : seq) {
        c = std::toupper(static_cast<unsigned char>(c));
    }
#endif
    if (canonical_orientation) {
        std::string rc(seq);
        ::reverse_complement(rc);
        if (rc < seq)
            seq = rc;
    }
    std::sort(labels.begin(), labels.end());
    uint64_t h = fnv1a64(release_id);
    h = fnv1a64("\t", h);
    h = fnv1a64(seq, h);
    h = fnv1a64("\t", h);
    for (const auto &l : labels) {
        h = fnv1a64(std::to_string(l.size()) + ":" + l, h);
    }
    return hex64(h);
}

std::string encode_runs(const std::vector<KmerInterval> &runs, uint64_t num_kmers) {
    std::string out;
    uint64_t pos = 0;
    for (const auto &run : runs) {
        if (run.begin > pos)
            out += "o" + std::to_string(run.begin - pos);
        out += "x" + std::to_string(run.size());
        pos = run.end;
    }
    if (num_kmers > pos)
        out += "o" + std::to_string(num_kmers - pos);
    return out;
}

static std::vector<KmerInterval> runs_of(const std::vector<bool> &mask) {
    std::vector<KmerInterval> runs;
    for (uint64_t i = 0; i < mask.size(); ++i) {
        if (!mask[i])
            continue;
        if (runs.empty() || runs.back().end != i) {
            runs.push_back({ i, i + 1 });
        } else {
            runs.back().end = i + 1;
        }
    }
    return runs;
}


SupportProfile resolve_support(LabelOracle &oracle,
                               std::string_view query,
                               const ResolveOptions &options) {
    SupportProfile profile;
    profile.k = oracle.get_k();
    profile.regime = oracle.regime();
    profile.support = options.support;

    if (options.labels.empty() == !options.discover) {
        throw std::invalid_argument("Specify either an explicit label list or discovery, not both");
    }
    if (options.support == Support::TRACE) {
        if (!oracle.has_coordinates())
            throw std::invalid_argument("Trace support requires an annotation with k-mer coordinates");
        if (oracle.regime() != Regime::BASIC) {
            throw std::invalid_argument("Trace support is only defined for BASIC (forward-strand) "
                                        "graphs: coordinates carry no strand in canonical indexes");
        }
    }
    if (query.size() < profile.k)
        throw std::invalid_argument("Query shorter than k");

    // a phase boundary: the caller can abandon the request here (ResolveOptions::stop)
    auto checkpoint = [&]() {
        if (options.stop && options.stop() && options.abandon)
            options.abandon();
    };
    profile.num_kmers = query.size() - profile.k + 1;
    checkpoint();
    std::vector<node_index> keys = oracle.keys_of_sequence(query);
    assert(keys.size() == profile.num_kmers);

    std::vector<bool> in_graph(profile.num_kmers);
    std::vector<node_index> present_keys;
    for (uint64_t i = 0; i < keys.size(); ++i) {
        in_graph[i] = keys[i] != npos;
        if (in_graph[i])
            present_keys.push_back(keys[i]);
    }
    profile.graph_runs = runs_of(in_graph);
    size_t num_present = present_keys.size();

    // ---- which labels to profile
    std::vector<LabelRef> refs;
    if (!options.labels.empty()) {
        for (const auto &name : options.labels) {
            refs.push_back(oracle.resolve_label(name));
        }
    } else {
        if (options.discover_kind == LabelKind::HEADER && !oracle.coord_to_header())
            throw std::invalid_argument("Header discovery requires a CoordToHeader index");
        checkpoint();

        std::vector<Row> rows;
        rows.reserve(present_keys.size());
        for (node_index key : present_keys) {
            rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
        }
        // (column, seq_id or 0) -> supported k-mers
        std::map<std::pair<Column, uint64_t>, uint64_t> counts;
        if (options.discover_kind == LabelKind::COLUMN) {
            for (const auto &row : oracle.get_rows(rows)) {
                for (Column c : row) {
                    counts[{ c, 0 }]++;
                }
            }
        } else {
            tsl::hopscotch_set<uint64_t> seen;
            for (const auto &row : oracle.get_row_tuples(rows)) {
                for (const auto &[c, coords] : row) {
                    seen.clear();
                    for (Coord coord : coords) {
                        uint64_t seq_id = oracle.map_coord(c, coord).first;
                        if (seen.insert(seq_id).second)
                            counts[{ c, seq_id }]++;
                    }
                }
            }
        }
        checkpoint();
        std::vector<std::pair<std::pair<Column, uint64_t>, uint64_t>> ranked(counts.begin(), counts.end());
        // more k-mers first; ties by column id, then seq_id
        std::stable_sort(ranked.begin(), ranked.end(),
                         [](const auto &a, const auto &b) { return a.second > b.second; });
        if (ranked.size() > options.discover_max_labels) {
            LabelTruncation trunc;
            trunc.total = ranked.size();
            trunc.kept = options.discover_max_labels;
            trunc.min_kept_kmers = ranked[trunc.kept - 1].second;
            trunc.max_dropped_kmers = ranked[trunc.kept].second;
            for (size_t i = trunc.kept; i < ranked.size(); ++i) {
                if (ranked[i].second == num_present)
                    trunc.dropped_full_length++;
            }
            ranked.resize(trunc.kept);
            profile.labels_truncated = trunc;
        }
        for (const auto &[id, count] : ranked) {
            LabelRef ref;
            ref.kind = options.discover_kind;
            ref.column = id.first;
            ref.seq_id = id.second;
            ref.name = ref.kind == LabelKind::COLUMN ? oracle.column_name(ref.column)
                                                     : oracle.header_name(ref.column, ref.seq_id);
            refs.push_back(ref);
        }
    }

    profile.labels.reserve(refs.size());
    for (const auto &ref : refs) {
        LabelProfile lp;
        lp.label = ref;
        profile.labels.push_back(lp);
    }
    if (refs.empty())
        return profile;

    // ---- support per k-mer
    const bool with_coords = options.support == Support::TRACE;
    LabelQuery query_labels(oracle, refs, with_coords);
    checkpoint();
    auto hits = query_labels.fetch(keys);
    checkpoint();

    // Presence (no coordinates) is scattered in one pass over the k-mers: iterating per
    // label and scanning each k-mer's hit list would cost O(labels x k-mers x hits),
    // which at discovery-scale label counts dominates everything else.
    if (!with_coords) {
        std::vector<std::vector<bool>> supported(refs.size(),
                                                 std::vector<bool>(profile.num_kmers, false));
        for (uint64_t i = 0; i < hits.size(); ++i) {
            for (const auto &h : hits[i]) {
                assert(h.label < refs.size());
                supported[h.label][i] = true;
            }
        }
        for (LabelId l = 0; l < refs.size(); ++l) {
            profile.labels[l].runs = runs_of(supported[l]);
            profile.labels[l].kmers_supported
                = std::count(supported[l].begin(), supported[l].end(), true);
        }
    }

    for (LabelId l = 0; l < refs.size() && with_coords; ++l) {
        LabelProfile &lp = profile.labels[l];
        std::vector<bool> supported(profile.num_kmers, false);
        {
            // trace-consistent runs: a chain of coordinates increasing by one per k-mer.
            // Column labels have coordinates in the column frame, header labels in the
            // sequence frame; either way consecutive k-mers must have consecutive coords.
            std::vector<Coord> live;
            for (uint64_t i = 0; i < hits.size(); ++i) {
                const SmallVector<Coord> *coords = nullptr;
                for (const auto &h : hits[i]) {
                    if (h.label == l) {
                        coords = &h.coords;
                        break;
                    }
                }
                if (!coords || coords->empty()) {
                    live.clear();
                    continue;
                }
                std::vector<Coord> next;
                for (Coord c : *coords) {
                    if (c && std::binary_search(live.begin(), live.end(), c - 1))
                        next.push_back(c);
                }
                if (next.empty()) {
                    if (!live.empty())
                        lp.trace_breaks.push_back(i);  // supported, but the trace jumped
                    // no coordinate continues the chain: start a new trace run here
                    next.assign(coords->begin(), coords->end());
                    lp.runs.push_back({ i, i + 1 });
                } else if (!lp.runs.empty() && lp.runs.back().end == i) {
                    lp.runs.back().end = i + 1;
                } else {
                    lp.runs.push_back({ i, i + 1 });
                }
                live.swap(next);
                supported[i] = true;
            }
        }
        lp.kmers_supported = std::count(supported.begin(), supported.end(), true);
    }

    // ---- seed candidates: identical maximal runs grouped
    std::map<KmerInterval, std::vector<LabelId>> groups;
    for (LabelId l = 0; l < profile.labels.size(); ++l) {
        for (const auto &run : profile.labels[l].runs) {
            if (run.size() >= options.min_block_kmers)
                groups[run].push_back(l);
        }
    }
    for (auto &[iv, labels] : groups) {
        std::sort(labels.begin(), labels.end());
        profile.candidates.push_back({ iv, std::move(labels) });
    }
    std::sort(profile.candidates.begin(), profile.candidates.end(),
              [](const SeedCandidate &a, const SeedCandidate &b) {
                  if (a.kmers.size() != b.kmers.size())
                      return a.kmers.size() > b.kmers.size();
                  if (a.labels.size() != b.labels.size())
                      return a.labels.size() > b.labels.size();
                  return a.kmers.begin < b.kmers.begin;
              });
    return profile;
}


/********************************* selection *********************************/

namespace {

// A Fenwick tree of counts over positions 0..n-1: add, the sum over positions < i, and
// the position of the k-th unit, each in O(log n).
struct Fenwick {
    std::vector<uint32_t> tree;
    size_t n;
    size_t top = 1;   // the largest power of two <= n (1 when n == 0)

    explicit Fenwick(size_t size) : tree(size + 1, 0), n(size) {
        while (top * 2 <= n) {
            top *= 2;
        }
    }
    void add(size_t i, uint32_t v) {
        for (++i; i <= n; i += i & (~i + 1)) {
            tree[i] += v;
        }
    }
    uint32_t prefix(size_t i) const {   // sum over positions < i
        uint32_t s = 0;
        for (; i > 0; i -= i & (~i + 1)) {
            s += tree[i];
        }
        return s;
    }
    size_t kth(uint32_t k) const {   // the 0-based position holding the k-th unit (1-based k)
        size_t pos = 0;
        for (size_t step = top; step; step >>= 1) {
            if (pos + step <= n && tree[pos + step] < k) {
                pos += step;
                k -= tree[pos];
            }
        }
        return pos;
    }
};

// does label l have a run containing [begin, end)?
bool covers(const LabelProfile &lp, const KmerInterval &iv) {
    for (const auto &run : lp.runs) {
        if (run.begin <= iv.begin && run.end >= iv.end)
            return true;
    }
    return false;
}

std::vector<LabelId> covering_labels(const SupportProfile &profile, const KmerInterval &iv) {
    std::vector<LabelId> out;
    for (LabelId l = 0; l < profile.labels.size(); ++l) {
        if (covers(profile.labels[l], iv))
            out.push_back(l);
    }
    return out;
}

struct Picked {
    KmerInterval kmers;
    std::vector<LabelId> labels;
    std::vector<std::string> labels_not_covering;
};

} // namespace

SeedSelection select_seeds(const SupportProfile &profile,
                           std::string_view query,
                           const SelectionPolicy &policy_in,
                           bool canonical_orientation) {
    SeedSelection selection;
    selection.policy = policy_in;
    SelectionPolicy &policy = selection.policy;
    if (policy.min_block_bp < profile.k)
        policy.min_block_bp = profile.k;
    const uint64_t min_kmers = policy.min_block_bp - profile.k + 1;
    selection.num_candidates = profile.candidates.size();

    std::vector<Picked> picked;

    switch (policy.policy) {
        case SelectionPolicy::LONGEST_FIRST: {
            for (const auto &cand : profile.candidates) {
                if (cand.kmers.size() < min_kmers)
                    continue;
                selection.num_eligible++;
                if (picked.size() < policy.max_seeds)
                    picked.push_back({ cand.kmers, cand.labels, {} });
            }
            break;
        }
        case SelectionPolicy::MAX_SUPPORT: {
            // Maximize the number of labels whose support run *contains* the chosen
            // interval [a, b): a run [s, e) qualifies iff s <= a and e >= b.
            //
            // For a fixed a, the count |{e >= b}| is non-increasing in b, so the best
            // count is reached at the shortest admissible b = a + min_kmers, and among
            // the intervals with that count the longest one ends at the smallest
            // qualifying end, i.e. b = min{e : e >= a + min_kmers}. So one binary
            // search per candidate a suffices. Both conditions only change at
            // endpoints, so a ranges over run starts and over (e - min_kmers).
            // This is O(R log R) per round for R runs, instead of scanning pairs.
            struct Run { uint64_t begin, end; };
            std::vector<Run> runs;
            for (const auto &lp : profile.labels) {
                for (const auto &run : lp.runs) {
                    if (run.size() >= min_kmers)
                        runs.push_back({ run.begin, run.end });
                }
            }
            const std::vector<Run> all_runs = [&runs]() {
                std::sort(runs.begin(), runs.end(),
                          [](const Run &a, const Run &b) { return a.begin < b.begin; });
                return runs;
            }();

            // intervals already claimed by a picked seed, as sorted disjoint ranges
            std::vector<KmerInterval> taken;

            while (picked.size() < policy.max_seeds) {
                // Rebuild the runs with the claimed intervals subtracted, so later rounds
                // see clipped runs and their endpoints. Reusing the original endpoints and
                // merely rejecting overlapping intervals loses the best remaining interval
                // (e.g. with A=[0,10), B=[0,5) and [0,5) taken, [5,10) is only reachable
                // through the clipped run's own endpoints).
                runs.clear();
                for (const Run &run : all_runs) {
                    uint64_t begin = run.begin;
                    for (const auto &t : taken) {
                        if (t.end <= begin || t.begin >= run.end)
                            continue;
                        if (t.begin > begin && t.begin - begin >= min_kmers)
                            runs.push_back({ begin, t.begin });
                        begin = std::max(begin, t.end);
                    }
                    if (begin < run.end && run.end - begin >= min_kmers)
                        runs.push_back({ begin, run.end });
                }
                std::sort(runs.begin(), runs.end(),
                          [](const Run &a, const Run &b) { return a.begin < b.begin; });
                if (runs.empty())
                    break;

                std::vector<uint64_t> candidate_a;
                for (const auto &run : runs) {
                    candidate_a.push_back(run.begin);
                    if (run.end >= min_kmers)
                        candidate_a.push_back(run.end - min_kmers);
                }
                std::sort(candidate_a.begin(), candidate_a.end());
                candidate_a.erase(std::unique(candidate_a.begin(), candidate_a.end()),
                                  candidate_a.end());

                std::optional<Picked> best;
                size_t best_count = 0;
                size_t next_run = 0;
                // the ends of the runs with begin <= a, as counts over the compressed run
                // ends: an insertion, "how many ends >= x" and "the smallest end >= x" each
                // cost O(log r) (keeping them in a sorted vector made the sweep O(r^2))
                std::vector<uint64_t> coords;
                coords.reserve(runs.size());
                for (const Run &run : runs) {
                    coords.push_back(run.end);
                }
                std::sort(coords.begin(), coords.end());
                coords.erase(std::unique(coords.begin(), coords.end()), coords.end());
                auto index_of = [&](uint64_t x) {
                    return static_cast<size_t>(std::lower_bound(coords.begin(), coords.end(), x)
                                               - coords.begin());
                };
                Fenwick ends(coords.size());
                uint32_t inserted = 0;
                for (uint64_t a : candidate_a) {
                    while (next_run < runs.size() && runs[next_run].begin <= a) {
                        ends.add(index_of(runs[next_run].end), 1);
                        inserted++;
                        next_run++;
                    }
                    if (a + min_kmers > profile.num_kmers)
                        break;
                    const size_t lo = index_of(a + min_kmers);
                    const uint32_t below = lo < coords.size() ? ends.prefix(lo) : inserted;
                    if (below == inserted)
                        continue;
                    KmerInterval iv { a, coords[ends.kth(below + 1)] };
                    size_t count = inserted - below;
                    bool better = !best
                        || count > best_count
                        || (count == best_count
                            && (iv.size() > best->kmers.size()
                                || (iv.size() == best->kmers.size() && iv.begin < best->kmers.begin)));
                    if (better) {
                        // only the count is needed to compare candidates; the label
                        // list is materialised once, for the winner, below (doing it on
                        // every improvement made a round O(r^2) in the number of runs)
                        best = Picked{ iv, {}, {} };
                        best_count = count;
                    }
                }
                if (!best)
                    break;
                best->labels = covering_labels(profile, best->kmers);
                assert(best->labels.size() == best_count);
                taken.push_back(best->kmers);
                std::sort(taken.begin(), taken.end());
                picked.push_back(*best);
            }
            selection.num_eligible = picked.size();
            break;
        }
        case SelectionPolicy::EXPLICIT: {
            for (const auto &ex : policy.explicit_seeds) {
                if (ex.kmers.end > profile.num_kmers || ex.kmers.begin >= ex.kmers.end)
                    throw std::invalid_argument("Explicit seed interval out of range");
                bool in_one_run = false;
                for (const auto &run : profile.graph_runs) {
                    if (run.begin <= ex.kmers.begin && run.end >= ex.kmers.end)
                        in_one_run = true;
                }
                if (!in_one_run)
                    throw std::invalid_argument("Explicit seed interval is not fully present in the graph");
                Picked p;
                p.kmers = ex.kmers;
                for (const auto &name : ex.labels) {
                    auto it = std::find_if(profile.labels.begin(), profile.labels.end(),
                                           [&](const LabelProfile &lp) { return lp.label.name == name; });
                    if (it == profile.labels.end())
                        throw std::invalid_argument("Explicit seed label was not profiled: '" + name + "'");
                    if (covers(*it, ex.kmers)) {
                        p.labels.push_back(it - profile.labels.begin());
                    } else {
                        p.labels_not_covering.push_back(name);
                    }
                }
                selection.num_eligible++;
                if (!p.labels.empty())
                    picked.push_back(std::move(p));
            }
            break;
        }
    }

    // merge overlapping picks onto their common interval
    if (policy.merge_overlapping && policy.policy != SelectionPolicy::EXPLICIT) {
        std::sort(picked.begin(), picked.end(),
                  [](const Picked &a, const Picked &b) { return a.kmers < b.kmers; });
        std::vector<Picked> merged;
        for (auto &p : picked) {
            if (!merged.empty() && merged.back().kmers.end > p.kmers.begin) {
                Picked &m = merged.back();
                KmerInterval common { std::max(m.kmers.begin, p.kmers.begin),
                                      std::min(m.kmers.end, p.kmers.end) };
                if (common.size() >= min_kmers) {
                    m.kmers = common;
                    m.labels = covering_labels(profile, common);
                    continue;
                }
            }
            merged.push_back(std::move(p));
        }
        picked.swap(merged);
    }

    // freeze
    for (auto &p : picked) {
        FrozenSeed seed;
        seed.kmers = p.kmers;
        seed.sequence = std::string(query.substr(p.kmers.begin, p.kmers.size() + profile.k - 1));
        seed.labels_not_covering = p.labels_not_covering;

        std::vector<std::pair<std::string, LabelId>> ordered;
        for (LabelId l : p.labels) {
            ordered.emplace_back(profile.labels[l].label.name, l);
        }
        switch (policy.label_order) {
            case SelectionPolicy::HASH:
                std::sort(ordered.begin(), ordered.end(), [&](const auto &a, const auto &b) {
                    std::string sa = a.first + "\t" + std::to_string(policy.sample_seed);
                    std::string sb = b.first + "\t" + std::to_string(policy.sample_seed);
                    auto ha = fnv1a64(sa), hb = fnv1a64(sb);
                    return std::tie(ha, a.first) < std::tie(hb, b.first);
                });
                break;
            case SelectionPolicy::COLUMN_ID:
                std::sort(ordered.begin(), ordered.end(), [&](const auto &a, const auto &b) {
                    const auto &ra = profile.labels[a.second].label, &rb = profile.labels[b.second].label;
                    return std::tie(ra.column, ra.seq_id) < std::tie(rb.column, rb.seq_id);
                });
                break;
            case SelectionPolicy::KMERS_SUPPORTED:
                std::sort(ordered.begin(), ordered.end(), [&](const auto &a, const auto &b) {
                    auto ka = profile.labels[a.second].kmers_supported;
                    auto kb = profile.labels[b.second].kmers_supported;
                    return std::tie(kb, a.first) < std::tie(ka, b.first);  // more first
                });
                break;
        }
        seed.population.supporting_total = ordered.size();
        seed.population.order = policy.label_order;
        seed.population.sample_seed = policy.sample_seed;
        uint64_t digest = 0xcbf29ce484222325ULL;
        for (size_t i = 0; i < ordered.size(); ++i) {
            if (i < policy.max_labels_per_seed) {
                seed.labels.push_back(ordered[i].first);
            } else {
                seed.population.dropped_count++;
                digest = fnv1a64(ordered[i].first + "\n", digest);
                if (ordered.size() - policy.max_labels_per_seed <= 1000)
                    seed.population.dropped.push_back(ordered[i].first);
            }
        }
        seed.population.included = seed.labels.size();
        if (seed.population.dropped_count)
            seed.population.dropped_digest = hex64(digest);
        std::sort(seed.labels.begin(), seed.labels.end());
        seed.seed_id = make_seed_id(policy.release_id, seed.sequence, canonical_orientation, seed.labels);
        selection.seeds.push_back(std::move(seed));
    }

    for (size_t i = 0; i < selection.seeds.size(); ++i) {
        for (size_t j = 0; j < selection.seeds.size(); ++j) {
            if (i == j)
                continue;
            const auto &a = selection.seeds[i].kmers, &b = selection.seeds[j].kmers;
            if (a.begin < b.end && b.begin < a.end)
                selection.seeds[i].overlaps_with.push_back(j);
        }
    }
    return selection;
}

} // namespace traversal
} // namespace graph
} // namespace mtg
