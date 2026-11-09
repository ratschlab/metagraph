#include "resolve.hpp"

#include <algorithm>
#include <cctype>
#include <optional>

#include <tsl/hopscotch_map.h>

#include "graph/annotated_dbg.hpp"
#include "annotation/coord_to_header.hpp"
#include "common/seq_tools/reverse_complement.hpp"


namespace mtg {
namespace graph {
namespace traversal {

using annot::matrix::BinaryMatrix;
using annot::matrix::MultiIntMatrix;

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

std::string make_seed_id(std::string_view sequence,
                         bool canonical_orientation,
                         std::vector<std::string> labels) {
    // defined over the sequence after the build's case mapping, so that a seed resubmitted
    // in another case keeps its id
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
    // a leading tab, as in the ids of the labeled-traversal prototype
    uint64_t h = fnv1a64("\t");
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


namespace {

// The rows a profile of explicit labels decodes at once
constexpr size_t kResolveBatchRows = 4096;

// The support of one label, built k-mer by k-mer in ascending order (each k-mer at most once)
struct LabelSupport {
    uint64_t kmers = 0;
    std::vector<KmerInterval> runs;

    void add(uint64_t i) {
        // present at the k-mer once, however many hits name it
        if (!runs.empty() && runs.back().end == i + 1)
            return;
        kmers++;
        if (!runs.empty() && runs.back().end == i) {
            runs.back().end = i + 1;
        } else {
            runs.push_back({ i, i + 1 });
        }
    }
};

} // namespace

SupportProfile resolve_support(const LabelOracle &oracle,
                               std::string_view query,
                               const ResolveOptions &options) {
    SupportProfile profile;
    profile.k = oracle.get_k();
    profile.regime = oracle.regime();

    if (options.labels.empty() == !options.discover)
        throw std::invalid_argument("Specify either an explicit label list or discovery, not both");
    if (query.size() < profile.k)
        throw std::invalid_argument("Query shorter than k");
    if (options.discover && !options.discover_max_labels)
        throw std::invalid_argument("Discovery needs discover_max_labels of at least 1");

    profile.num_kmers = query.size() - profile.k + 1;
    const std::vector<node_index> keys = oracle.keys_of_sequence(query);
    assert(keys.size() == profile.num_kmers);

    profile.graph_runs = runs_of(profile.num_kmers, [&](uint64_t i) { return keys[i] != npos; });

    std::vector<LabelRef> refs;
    std::vector<LabelSupport> support;
    if (!options.labels.empty()) {
        for (const auto &name : options.labels) {
            refs.push_back(oracle.resolve_label(name));
        }
        LabelQuery query_labels(oracle, refs);
        support.resize(refs.size());
        // decode the rows kResolveBatchRows at a time
        for (size_t begin = 0; begin < keys.size(); begin += kResolveBatchRows) {
            const size_t end = std::min(keys.size(), begin + kResolveBatchRows);
            const auto hits = query_labels.fetch(
                    std::vector<node_index>(keys.begin() + begin, keys.begin() + end));
            for (size_t i = 0; i < hits.size(); ++i) {
                for (LabelId l : hits[i]) {
                    support[l].add(begin + i);
                }
            }
        }
    } else {
        if (options.discover_kind == LabelKind::HEADER && !oracle.coord_to_header())
            throw std::invalid_argument("Header discovery requires a CoordToHeader index");
        // the ranking of /search with signatures, so that both report the same labels,
        // counts and runs
        const bool headers = options.discover_kind == LabelKind::HEADER;
        size_t total = 0;
        auto top = oracle.anno_graph().get_top_label_signatures(
                query, options.discover_max_labels, 0.0, 0.0, &total, headers);
        for (const auto &[name, count, presence] : top) {
            // a name /traverse would refuse as ambiguous (a column and a header, or the
            // header of several records) is left out, so that every label can be resubmitted
            LabelRef ref;
            try {
                ref = oracle.resolve_label(name);
            } catch (const std::invalid_argument &) {
                continue;
            }
            assert(ref.kind == options.discover_kind);
            LabelSupport found;
            for (uint64_t i = 0; i < presence.size(); ++i) {
                if (presence[i])
                    found.add(i);
            }
            refs.push_back(std::move(ref));
            support.push_back(std::move(found));
        }
        // kept: the labels returned, after the ambiguous ones are left out
        if (total > options.discover_max_labels)
            profile.labels_truncated = LabelTruncation { refs.size(), total };
    }

    for (LabelId l = 0; l < refs.size(); ++l) {
        LabelProfile lp;
        lp.label = refs[l];
        lp.kmers_supported = support[l].kmers;
        lp.runs = std::move(support[l].runs);
        profile.labels.push_back(std::move(lp));
    }
    return profile;
}


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
};

// Repeatedly pick the interval of at least |min_kmers| k-mers contained in the
// most support runs, outside the intervals picked before.
//
// For a fixed start a, the count |{runs [s, e) : s <= a, e >= b}| is non-increasing in b, so
// the best count is reached at the shortest admissible b = a + min_kmers, and the longest
// interval with that count ends at b = min{e : e >= a + min_kmers}. So one search per
// candidate a suffices, and a ranges over run starts and over (e - min_kmers): O(R log R) per
// round for R runs.
std::vector<Picked> pick_max_support(const SupportProfile &profile, uint64_t min_kmers,
                                     size_t max_seeds) {
    struct Run { uint64_t begin, end; };
    std::vector<Run> all_runs;
    for (const auto &lp : profile.labels) {
        for (const auto &run : lp.runs) {
            if (run.size() >= min_kmers)
                all_runs.push_back({ run.begin, run.end });
        }
    }
    std::sort(all_runs.begin(), all_runs.end(),
              [](const Run &a, const Run &b) { return a.begin < b.begin; });

    std::vector<Picked> picked;
    // intervals already claimed by a picked seed, as sorted disjoint ranges
    std::vector<KmerInterval> taken;
    std::vector<Run> runs;
    while (picked.size() < max_seeds) {
        // the runs with the claimed intervals subtracted, so that later rounds see the
        // clipped runs' own endpoints (with A=[0,10), B=[0,5) and [0,5) taken, [5,10) is
        // only reachable through them)
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
        candidate_a.erase(std::unique(candidate_a.begin(), candidate_a.end()), candidate_a.end());

        // the ends of the runs with begin <= a, as counts over the compressed run ends
        std::vector<uint64_t> ends_at;
        for (const Run &run : runs) {
            ends_at.push_back(run.end);
        }
        std::sort(ends_at.begin(), ends_at.end());
        ends_at.erase(std::unique(ends_at.begin(), ends_at.end()), ends_at.end());
        auto index_of = [&](uint64_t x) {
            return static_cast<size_t>(std::lower_bound(ends_at.begin(), ends_at.end(), x)
                                       - ends_at.begin());
        };
        Fenwick ends(ends_at.size());
        uint32_t inserted = 0;
        size_t next_run = 0;
        std::optional<KmerInterval> best;
        size_t best_count = 0;
        for (uint64_t a : candidate_a) {
            while (next_run < runs.size() && runs[next_run].begin <= a) {
                ends.add(index_of(runs[next_run].end), 1);
                inserted++;
                next_run++;
            }
            if (a + min_kmers > profile.num_kmers)
                break;
            const size_t lo = index_of(a + min_kmers);
            const uint32_t below = lo < ends_at.size() ? ends.prefix(lo) : inserted;
            if (below == inserted)
                continue;
            const KmerInterval iv { a, ends_at[ends.kth(below + 1)] };
            const size_t count = inserted - below;
            if (!best || count > best_count
                    || (count == best_count && (iv.size() > best->size()
                            || (iv.size() == best->size() && iv.begin < best->begin)))) {
                best = iv;
                best_count = count;
            }
        }
        if (!best)
            break;
        picked.push_back({ *best, covering_labels(profile, *best) });
        assert(picked.back().labels.size() == best_count);
        taken.push_back(*best);
        std::sort(taken.begin(), taken.end());
    }
    return picked;
}

} // namespace

std::vector<FrozenSeed> select_seeds(const SupportProfile &profile,
                                     std::string_view query,
                                     const SelectionPolicy &policy_in,
                                     bool canonical_orientation) {
    std::vector<FrozenSeed> seeds;
    SelectionPolicy policy = policy_in;
    if (policy.min_block_bp < profile.k)
        policy.min_block_bp = profile.k;
    const uint64_t min_kmers = policy.min_block_bp - profile.k + 1;

    // the picks never overlap: each is outside the earlier ones
    const std::vector<Picked> picked = pick_max_support(profile, min_kmers, policy.max_seeds);

    // freeze
    for (auto &p : picked) {
        FrozenSeed seed;
        seed.kmers = p.kmers;
        seed.sequence = std::string(query.substr(p.kmers.begin, p.kmers.size() + profile.k - 1));

        // the labels in the order of the hash of their name (then the name)
        std::vector<std::pair<uint64_t, std::string>> ordered;
        for (LabelId l : p.labels) {
            const std::string &name = profile.labels[l].label.name;
            ordered.emplace_back(fnv1a64(name + "\t0"), name);
        }
        std::sort(ordered.begin(), ordered.end());

        seed.population.supporting_total = ordered.size();
        for (size_t i = 0; i < ordered.size(); ++i) {
            if (i < policy.max_labels_per_seed) {
                seed.labels.push_back(ordered[i].second);
            } else {
                seed.population.dropped_count++;
                if (ordered.size() - policy.max_labels_per_seed <= 1000)
                    seed.population.dropped.push_back(ordered[i].second);
            }
        }
        std::sort(seed.labels.begin(), seed.labels.end());
        seed.seed_id = make_seed_id(seed.sequence, canonical_orientation, seed.labels);
        seeds.push_back(std::move(seed));
    }
    return seeds;
}

} // namespace traversal
} // namespace graph
} // namespace mtg
