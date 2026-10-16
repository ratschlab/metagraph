#include "resolve.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "graph/annotated_dbg.hpp"
#include "common/seq_tools/reverse_complement.hpp"

namespace mtg::graph::traversal {

// Selection and identity logic extracted from gr/labeled-traversal 07dfa4875fcf.
namespace {
uint64_t fnv1a64(std::string_view data, uint64_t hash = 0xcbf29ce484222325ULL) {
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

} // namespace

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

SupportProfile resolve_support(const AnnotatedDBG &graph, std::string_view query,
                               const ResolveOptions &options) {
    if (options.labels.empty() == !options.discover)
        throw std::invalid_argument("Specify either labels or discover");
    if (!options.discover_max_labels || !options.min_block_kmers)
        throw std::invalid_argument("Resolve limits must be positive");
    SupportProfile profile;
    profile.k = graph.get_graph().get_k();
    if (query.size() < profile.k)
        throw std::invalid_argument("Query shorter than k");
    profile.num_kmers = query.size() - profile.k + 1;
    const auto &encoder = graph.get_annotator().get_label_encoder();
    std::unordered_map<uint64_t, LabelId> slots;
    auto add_column = [&](uint64_t column) -> LabelProfile& {
        auto [it, added] = slots.try_emplace(column, profile.labels.size());
        if (added)
            profile.labels.push_back({{column, encoder.decode(column)}, 0, {}});
        return profile.labels[it->second];
    };
    for (const auto &label : options.labels) {
        auto column = encoder.encode(label);
        if (slots.contains(column))
            throw std::invalid_argument("Duplicate resolve label: " + label);
        add_column(column);
    }
    using Node = AnnotatedDBG::node_index;
    std::vector<Node> keys;
    graph.get_graph().map_to_nodes(query, [&](Node node) { keys.push_back(node); });
    profile.graph_runs = runs_of(keys.size(), [&](size_t i) { return keys[i] != 0; });
    const auto &matrix = graph.get_annotator().get_matrix();
    // The prototype's column-membership path accumulates runs in query order.
    // Use the existing binary matrix interface, without LabelOracle's traversal,
    // coordinate, deadline and cache dependencies. Repeated rows share a decode
    // within each batch of at most 64 distinct rows.
    for (size_t begin = 0; begin < keys.size(); ) {
        std::vector<annot::matrix::BinaryMatrix::Row> rows;
        std::unordered_map<Node, size_t> positions;
        size_t end = begin;
        for (; end < keys.size(); ++end) {
            if (!keys[end] || positions.contains(keys[end]))
                continue;
            if (rows.size() == 64)
                break;
            positions.emplace(keys[end], rows.size());
            rows.push_back(AnnotatedDBG::graph_to_anno_index(keys[end]));
        }
        const auto decoded = matrix.get_rows(rows);
        for (size_t i = begin; i < end; ++i) {
            if (!keys[i])
                continue;
            for (uint64_t column : decoded[positions.at(keys[i])]) {
                if (!options.discover && !slots.contains(column))
                    continue;
                auto &label = add_column(column);
                ++label.kmers_supported;
                if (label.runs.empty() || label.runs.back().end != i)
                    label.runs.push_back({i, i + 1});
                else
                    label.runs.back().end = i + 1;
            }
        }
        begin = end;
    }
    if (options.discover) {
        std::sort(profile.labels.begin(), profile.labels.end(), [](const auto &a, const auto &b) {
            return a.kmers_supported != b.kmers_supported
                ? a.kmers_supported > b.kmers_supported : a.label.column < b.label.column;
        });
        if (profile.labels.size() > options.discover_max_labels) {
            const auto kept = options.discover_max_labels;
            uint64_t present = 0;
            for (const auto &run : profile.graph_runs) present += run.size();
            size_t full = 0;
            for (size_t i = kept; i < profile.labels.size(); ++i)
                full += profile.labels[i].kmers_supported == present;
            profile.labels_truncated = LabelTruncation{kept, profile.labels.size(),
                profile.labels[kept - 1].kmers_supported, profile.labels[kept].kmers_supported, full};
            profile.labels.resize(kept);
        }
    }
    // The prototype's identical maximal-run grouping and deterministic ordering.
    std::map<KmerInterval, std::vector<LabelId>> groups;
    for (LabelId l = 0; l < profile.labels.size(); ++l)
        for (const auto &run : profile.labels[l].runs)
            if (run.size() >= options.min_block_kmers)
                groups[run].push_back(l);
    for (auto &[iv, labels] : groups)
        profile.candidates.push_back({iv, std::move(labels)});
    std::sort(profile.candidates.begin(), profile.candidates.end(), [](const auto &a, const auto &b) {
        if (a.kmers.size() != b.kmers.size()) return a.kmers.size() > b.kmers.size();
        if (a.labels.size() != b.labels.size()) return a.labels.size() > b.labels.size();
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
                    return ra.column < rb.column;
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

} // namespace mtg::graph::traversal
