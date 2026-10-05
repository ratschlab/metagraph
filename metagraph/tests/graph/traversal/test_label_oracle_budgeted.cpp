// The budget-aware reads of LabelQuery and LabelRecorder (stage 3 of
// DESIGN-traverse-graphlet.md §14.1): the same answers as the unbudgeted reads, all or
// nothing, every key admitted against its demand whether it is cached or not.
#include "gtest/gtest.h"

#include <chrono>
#include <future>
#include <memory>
#include <random>
#include <set>
#include <thread>

#include "tests/annotation/test_annotated_dbg_helpers.hpp"

#include "graph/traversal/label_oracle.hpp"
#include "graph/traversal/walker.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/coord_to_header.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"

#if USE_JEMALLOC
extern "C" int mallctl(const char *name, void *oldp, size_t *oldlenp, void *newp, size_t newlen);
#endif


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;
using annot::matrix::DecodeBudget;

const size_t kK = 11;

std::vector<std::string> sequences(size_t n, size_t len, uint32_t seed) {
    std::mt19937 gen(seed);
    std::vector<std::string> seqs;
    for (size_t i = 0; i < n; ++i) {
        std::string s(len, 'A');
        for (char &c : s) {
            c = "ACGT"[gen() % 4];
        }
        seqs.push_back(s);
    }
    return seqs;
}

// Labels sharing long stretches (rows of several labels, diffs that cancel), two records
// in column F for the header labels
struct Fixture {
    std::unique_ptr<AnnotatedDBG> anno;
    std::unique_ptr<annot::CoordToHeader> cth;
    std::vector<std::string> seqs;
    std::vector<std::string> labels;

    explicit Fixture(bool coordinates) {
        auto rnd = sequences(8, 90, 5);
        const std::string shared = rnd[7];
        for (size_t i = 0; i < 7; ++i) {
            seqs.push_back(rnd[i].substr(0, 30) + shared.substr(i, 60) + rnd[i].substr(30));
            labels.push_back(i < 2 ? "F" : "L" + std::to_string(i));
        }
        std::vector<uint64_t> starts;
        uint64_t f = 0;
        for (size_t i = 0; i < seqs.size(); ++i) {
            if (labels[i] == "F") {
                starts.push_back(f);
                f += seqs[i].size() - kK + 1;
            } else {
                starts.push_back(0);
            }
        }
        anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
                kK, seqs, labels, DeBruijnGraph::BASIC, coordinates,
                coordinates ? starts : std::vector<uint64_t>{});
        if (coordinates) {
            const auto &encoder = anno->get_annotator().get_label_encoder();
            std::vector<std::vector<std::string>> headers(encoder.size());
            std::vector<std::vector<uint64_t>> num_kmers(encoder.size());
            for (size_t i = 0; i < seqs.size(); ++i) {
                const size_t c = encoder.encode(labels[i]);
                headers[c].push_back(labels[i] == "F" ? "acc" + std::to_string(i) : labels[i] + "_h");
                num_kmers[c].push_back(seqs[i].size() - kK + 1);
            }
            cth = std::make_unique<annot::CoordToHeader>(std::move(headers), std::move(num_kmers));
        }
    }

    // keys of every k-mer of every sequence, in batches as a walk asks for them
    std::vector<std::vector<node_index>> batches(const LabelOracle &oracle) const {
        std::vector<node_index> all;
        for (const auto &s : seqs) {
            for (node_index k : oracle.keys_of_sequence(s)) all.push_back(k);
        }
        std::vector<std::vector<node_index>> out;
        std::mt19937 gen(3);
        for (size_t begin = 0; begin < all.size(); ) {
            const size_t n = 1 + gen() % 40;
            std::vector<node_index> b(all.begin() + begin, all.begin() + std::min(all.size(), begin + n));
            if (b.size() > 2) {
                b.push_back(b[1]);             // a duplicate
                b.push_back(npos);             // a key that is not in the graph
            }
            out.push_back(std::move(b));
            begin += n;
        }
        return out;
    }
};

// the label sets a query is checked with: column labels (rows), column labels with
// coordinates (tuples), header labels (tuples, coordinates mapped to records)
std::vector<std::pair<std::vector<std::string>, bool>> label_sets(bool coordinates) {
    std::vector<std::pair<std::vector<std::string>, bool>> sets {
        { { "F", "L3", "L5" }, false },
    };
    if (coordinates) {
        sets.push_back({ { "F", "L2", "L6" }, true });
        sets.push_back({ { "acc0", "acc1", "L4_h" }, true });
        sets.push_back({ { "acc1", "L2" }, false });
    }
    return sets;
}

std::vector<LabelRef> refs(const LabelOracle &oracle, const std::vector<std::string> &names) {
    std::vector<LabelRef> out;
    for (const auto &n : names) out.push_back(oracle.resolve_label(n));
    return out;
}

// Q1, Q3: the budget-aware fetch and warm return the unbudgeted hits, and a key's costs are
// the same whether a fetch decoded it or found it cached
TEST(LabelOracleBudgetedQuery, SameHitsAsTheUnbudgetedFetch) {
    for (bool coordinates : { false, true }) {
        Fixture fx(coordinates);
        for (const auto &[names, with_coords] : label_sets(coordinates)) {
            LabelOracle oracle(*fx.anno, fx.cth.get());
            ASSERT_TRUE(oracle.decode_charged());
            LabelQuery plain(oracle, refs(oracle, names), with_coords);
            LabelQuery budgeted(oracle, refs(oracle, names), with_coords);
            budgeted.set_max_cache_bytes(1 << 20);
            std::map<node_index, KeyCost> seen;
            size_t from_cache = 0;
            for (const auto &batch : fx.batches(oracle)) {
                const auto expected = plain.fetch(batch);
                std::vector<LabelQuery::NodeHits> hits;
                std::vector<KeyCost> costs;
                hits.reserve(batch.size());
                costs.reserve(batch.size());
                DecodeBudget budget;
                size_t refused = 0;
                ASSERT_TRUE(budgeted.fetch(batch.data(), batch.size(), budget, &hits, &costs, &refused));
                ASSERT_EQ(expected, hits) << names[0];
                uint64_t held = 0;
                for (size_t i = 0; i < batch.size(); ++i) {
                    held += LabelQuery::held_bytes(hits[i]);
                    if (batch[i] == npos)
                        continue;
                    EXPECT_GE(costs[i].demand, LabelQuery::held_bytes(hits[i]));
                    auto [it, inserted] = seen.emplace(batch[i], costs[i]);
                    if (!inserted) {
                        from_cache++;
                        EXPECT_EQ(it->second.demand, costs[i].demand);
                        EXPECT_EQ(it->second.dependency_units, costs[i].dependency_units);
                    }
                }
                // what the call holds is what it returned
                EXPECT_EQ(held, budget.held());
                EXPECT_LE(budgeted.cache_bytes(), uint64_t(1) << 20);
                // the lookahead's read answers like a fetch's
                DecodeBudget warm_budget;
                budgeted.warm(batch, warm_budget);
                EXPECT_EQ(0u, warm_budget.held());
            }
            EXPECT_GT(from_cache, 10u);
        }
    }
}

// The efficiency pass, the row-diff path cache (LabelOracle::path_cache): the reads of an
// oracle with the cache answer every fetch, budgeted and not, query and recorder, exactly as
// an oracle without it — the hits and lists, the costs (each key's whole path), what a
// budgeted call holds, and every counter (requests, cache hits, rows fetched, mappings) —
// and the paths stop at cached rows. With a shared bound (make_room), the path cache keeps
// within what the label cache leaves of it.
TEST(LabelOraclePathCache, SameAnswersAndCounters) {
    auto same_counters = [](const LabelOracle::Counters &a, const LabelOracle::Counters &b) {
        return a.keys_mapped == b.keys_mapped && a.rows_requested == b.rows_requested
            && a.cache_hits == b.cache_hits && a.rows_fetched == b.rows_fetched
            && a.tuple_rows_fetched == b.tuple_rows_fetched
            && a.direct_reads == b.direct_reads && a.coords_mapped == b.coords_mapped;
    };
    for (bool coordinates : { false, true }) {
        Fixture fx(coordinates);
        for (const auto &[names, with_coords] : label_sets(coordinates)) {
            for (bool shared : { false, true }) {
                LabelOracle without(*fx.anno, fx.cth.get());
                LabelOracle with(*fx.anno, fx.cth.get());
                with.set_path_cache_max(uint64_t(1) << 20);
                ASSERT_TRUE(with.path_cached());
                LabelQuery q0(without, refs(without, names), with_coords);
                LabelQuery q1(with, refs(with, names), with_coords);
                LabelQuery b0(without, refs(without, names), with_coords);
                LabelQuery b1(with, refs(with, names), with_coords);
                b0.set_max_cache_bytes(1 << 16);
                b1.set_max_cache_bytes(1 << 16);
                if (shared) {
                    // the walker's sharing under a memory budget: the path cache within what
                    // the budgeted query's cache leaves of 64 KiB
                    with.path_cache().set_bound(1 << 16, [&]() {
                        return b1.cache_bytes() < (1 << 16) ? (1 << 16) - b1.cache_bytes() : 0;
                    });
                }
                size_t calls = 0;
                fx.batches(without);     // keys_mapped alike
                for (const auto &batch : fx.batches(with)) {
                    ASSERT_EQ(q0.fetch(batch), q1.fetch(batch)) << names[0];
                    std::vector<LabelQuery::NodeHits> h0, h1;
                    std::vector<KeyCost> c0, c1;
                    h0.reserve(batch.size()); h1.reserve(batch.size());
                    c0.reserve(batch.size()); c1.reserve(batch.size());
                    DecodeBudget d0, d1;
                    size_t r0 = 0, r1 = 0;
                    ASSERT_TRUE(b0.fetch(batch.data(), batch.size(), d0, &h0, &c0, &r0));
                    ASSERT_TRUE(b1.fetch(batch.data(), batch.size(), d1, &h1, &c1, &r1));
                    ASSERT_EQ(h0, h1);
                    for (size_t i = 0; i < batch.size(); ++i) {
                        EXPECT_EQ(c0[i].demand, c1[i].demand) << names[0] << " call " << calls;
                        EXPECT_EQ(c0[i].dependency_units, c1[i].dependency_units);
                    }
                    EXPECT_EQ(d0.held(), d1.held());
                    if (shared) {
                        EXPECT_LE(with.path_cache().bytes() + b1.cache_bytes(),
                                  std::max<uint64_t>(1 << 16, b1.cache_bytes()));
                    }
                    calls++;
                }
                EXPECT_TRUE(same_counters(without.counters(), with.counters())) << names[0];
                EXPECT_GT(with.path_cache().hits(), 0u) << names[0];
            }
        }
        // annotate: the recorder's lists and names
        for (LabelKind kind : { LabelKind::COLUMN, LabelKind::HEADER }) {
            if (kind == LabelKind::HEADER && !coordinates)
                continue;
            LabelOracle without(*fx.anno, fx.cth.get());
            LabelOracle with(*fx.anno, fx.cth.get());
            with.set_path_cache_max(uint64_t(1) << 20);
            LabelRecorder r0(without, kind, 3), r1(with, kind, 3);
            fx.batches(without);     // keys_mapped alike
            for (const auto &batch : fx.batches(with)) {
                const auto l0 = r0.fetch(batch);
                const auto l1 = r1.fetch(batch);
                ASSERT_EQ(l0.size(), l1.size());
                for (size_t i = 0; i < l0.size(); ++i) {
                    EXPECT_EQ(l0[i].labels, l1[i].labels);
                    EXPECT_EQ(l0[i].total, l1[i].total);
                }
            }
            ASSERT_EQ(r0.labels().size(), r1.labels().size());
            for (size_t i = 0; i < r0.labels().size(); ++i) {
                EXPECT_EQ(r0.labels()[i].name, r1.labels()[i].name);
            }
            EXPECT_TRUE(same_counters(without.counters(), with.counters()));
            EXPECT_GT(with.path_cache().hits(), 0u);
        }
    }
}

// Q2, Q4: a key is admitted exactly when its demand fits what the call has left, cached or
// not; a refusal returns nothing and leaves the cache and the counters as they were
TEST(LabelOracleBudgetedQuery, AdmitsExactlyAndRefusesWhole) {
    Fixture fx(true);
    for (const auto &[names, with_coords] : label_sets(true)) {
        LabelOracle oracle(*fx.anno, fx.cth.get());
        LabelQuery probe(oracle, refs(oracle, names), with_coords);
        std::vector<node_index> keys = oracle.keys_of_sequence(fx.seqs[2]);
        keys.resize(12);
        std::vector<LabelQuery::NodeHits> hits;
        std::vector<KeyCost> costs;
        hits.reserve(keys.size());
        costs.reserve(keys.size());
        DecodeBudget unlimited;
        size_t refused = 0;
        ASSERT_TRUE(probe.fetch(keys.data(), keys.size(), unlimited, &hits, &costs, &refused));
        for (bool cached : { false, true }) {
            // The smallest budget that admits keys 0 .. i: key k needs its demand beside what
            // the keys before it returned. One byte less refuses at the first key that needs
            // all of it.
            uint64_t before = 0, need = 0;
            size_t first = 0;
            for (size_t i = 0; i < keys.size(); ++i) {
                if (before + costs[i].demand > need) {
                    need = before + costs[i].demand;
                    first = i;
                }
                LabelQuery query(oracle, refs(oracle, names), with_coords);
                if (cached) {
                    std::vector<LabelQuery::NodeHits> h;
                    std::vector<KeyCost> c;
                    h.reserve(keys.size());
                    c.reserve(keys.size());
                    DecodeBudget b;
                    ASSERT_TRUE(query.fetch(keys.data(), keys.size(), b, &h, &c, &refused));
                }
                for (uint64_t extra : { uint64_t(0), uint64_t(1) }) {
                    const uint64_t cache_bytes = query.cache_bytes();
                    const auto counters = oracle.counters();
                    DecodeBudget budget(need - extra);
                    budget.charge(0);
                    std::vector<LabelQuery::NodeHits> h { LabelQuery::NodeHits{} };
                    std::vector<KeyCost> c { KeyCost{} };
                    h.reserve(1 + i + 1);
                    c.reserve(1 + i + 1);
                    const bool ok = query.fetch(keys.data(), i + 1, budget, &h, &c, &refused);
                    const std::string what = names[0] + " key " + std::to_string(i)
                                           + (cached ? " cached" : " decoded");
                    if (!extra) {
                        EXPECT_TRUE(ok) << what;
                        if (ok) {
                            EXPECT_EQ(hits[i], h.back()) << what;
                        }
                    } else {
                        ASSERT_FALSE(ok) << what;
                        EXPECT_EQ(first, refused) << what;
                        EXPECT_EQ(1u, h.size()) << what;
                        EXPECT_EQ(0u, budget.held()) << what;
                        EXPECT_EQ(counters.rows_requested, oracle.counters().rows_requested);
                        EXPECT_EQ(counters.cache_hits, oracle.counters().cache_hits);
                        // the physical counters (timing) count the decoding a refused fetch
                        // did: they never go back (review of stage 3, F5: the spec now says so)
                        EXPECT_LE(counters.rows_fetched + counters.tuple_rows_fetched,
                                  oracle.counters().rows_fetched + oracle.counters().tuple_rows_fetched);
                        EXPECT_LE(counters.fetch_seconds, oracle.counters().fetch_seconds);
                        // the refusal says which key and why: a demand that does not fit
                        // what was left at its position, or a read that did not fit alone
                        const FetchRefusal &why = query.refusal();
                        EXPECT_EQ(first, why.position) << what;
                        EXPECT_EQ(need - extra - why.held, why.left) << what;
                        if (why.cause == FetchRefusal::DEMAND) {
                            EXPECT_EQ(costs[first].demand, why.demand) << what;
                            EXPECT_GT(why.demand, why.left) << what;
                        } else {
                            EXPECT_EQ(FetchRefusal::DECODE, why.cause) << what;
                            EXPECT_GT(why.need, why.left) << what;
                        }
                        EXPECT_EQ(cache_bytes, query.cache_bytes()) << what;
                        EXPECT_EQ(need - costs[first].demand, query.refused_held()) << what;
                    }
                }
                before += LabelQuery::held_bytes(hits[i]);
            }
        }
    }
}

// Q1, Q2 for the recorder: the same lists and dictionary as the unbudgeted reads, names
// priced inside the call and given only when the whole call is admitted
TEST(LabelOracleBudgetedRecorder, SameListsAndNames) {
    Fixture fx(true);
    for (LabelKind kind : { LabelKind::COLUMN, LabelKind::HEADER }) {
        for (size_t cap : { size_t(1), size_t(3), size_t(64) }) {
            LabelOracle oracle(*fx.anno, fx.cth.get());
            LabelRecorder plain(oracle, kind, cap);
            LabelRecorder budgeted(oracle, kind, cap);
            budgeted.set_max_cache_bytes(1 << 20);
            uint64_t priced = 0;
            auto name_bytes = [&](std::string_view name) {
                return 100 + name.size();
            };
            for (const auto &batch : fx.batches(oracle)) {
                const auto expected = plain.fetch(batch);
                std::vector<LabelRecorder::NodeLabels> lists;
                std::vector<KeyCost> costs;
                lists.reserve(batch.size());
                costs.reserve(batch.size());
                const size_t named = budgeted.labels().size();
                // a refused call (a budget below the first key's demand) names nothing
                {
                    DecodeBudget none(1);
                    std::vector<LabelRecorder::NodeLabels> l;
                    std::vector<KeyCost> c;
                    l.reserve(batch.size());
                    c.reserve(batch.size());
                    size_t refused = 0;
                    if (batch[0] != npos) {
                        EXPECT_FALSE(budgeted.fetch(batch.data(), batch.size(), none, &l, &c,
                                                    &refused, name_bytes));
                        EXPECT_TRUE(l.empty());
                        EXPECT_EQ(named, budgeted.labels().size());
                    }
                }
                DecodeBudget budget;
                size_t refused = 0;
                ASSERT_TRUE(budgeted.fetch(batch.data(), batch.size(), budget, &lists, &costs,
                                           &refused, name_bytes));
                ASSERT_EQ(expected.size(), lists.size());
                for (size_t i = 0; i < lists.size(); ++i) {
                    EXPECT_EQ(expected[i].labels, lists[i].labels);
                    EXPECT_EQ(expected[i].total, lists[i].total);
                }
                ASSERT_EQ(plain.labels().size(), budgeted.labels().size());
                uint64_t names = 0;
                for (size_t id = named; id < budgeted.labels().size(); ++id) {
                    EXPECT_EQ(plain.labels()[id].name, budgeted.labels()[id].name);
                    names += 100 + budgeted.labels()[id].name.size();
                }
                EXPECT_EQ(names, budgeted.last_names_bytes());
                // each new label is charged its provisional naming beside its name, which
                // bounds what the call's pending labels held (review of stage 3, F1)
                const uint64_t named_now = budgeted.labels().size() - named;
                EXPECT_EQ(named_now * LabelRecorder::kNamingBytes, budgeted.last_naming_bytes());
                EXPECT_LE(LabelRecorder::pending_bytes(named_now), budgeted.last_naming_bytes());
                priced += names;
                uint64_t held = names + budgeted.last_naming_bytes();
                for (const auto &l : lists) {
                    held += LabelRecorder::held_bytes(l);
                }
                EXPECT_EQ(held, budget.held());
                EXPECT_LE(budgeted.cache_bytes(), uint64_t(1) << 20);
                DecodeBudget warm_budget;
                budgeted.warm(batch, warm_budget);
                EXPECT_EQ(0u, warm_budget.held());
            }
            EXPECT_GT(priced, 0u);
        }
    }
}

// F2: the recorder says why a key was refused — its demand, or the dictionary labels it would
// name first — so that a stop by labels is not reported as a row that does not fit
TEST(LabelOracleBudgetedRecorder, RefusalNamesItsCause) {
    Fixture fx(true);
    for (LabelKind kind : { LabelKind::COLUMN, LabelKind::HEADER }) {
        LabelOracle oracle(*fx.anno, fx.cth.get());
        auto name_bytes = [](std::string_view name) { return 1000 + name.size(); };
        // the key of the shared stretch carrying the most labels
        node_index key = npos;
        size_t most = 0;
        for (node_index k : oracle.keys_of_sequence(fx.seqs[0])) {
            LabelRecorder probe(oracle, kind, 64);
            const auto lists = probe.fetch(std::vector<node_index>{ k });
            if (lists[0].labels.size() > most) {
                most = lists[0].labels.size();
                key = k;
            }
        }
        ASSERT_GE(most, 3u);
        LabelRecorder probe(oracle, kind, 64);
        std::vector<LabelRecorder::NodeLabels> l;
        std::vector<KeyCost> c;
        l.reserve(1);
        c.reserve(1);
        DecodeBudget unlimited;
        size_t refused = 0;
        ASSERT_TRUE(probe.fetch(&key, 1, unlimited, &l, &c, &refused, name_bytes));
        const uint64_t names = probe.last_names_bytes() + probe.last_naming_bytes();
        EXPECT_EQ(most * LabelRecorder::kNamingBytes, probe.last_naming_bytes());
        for (uint64_t budget_bytes : { c[0].demand + names - 1, c[0].demand }) {
            LabelRecorder recorder(oracle, kind, 64);
            std::vector<LabelRecorder::NodeLabels> out;
            std::vector<KeyCost> costs;
            out.reserve(1);
            costs.reserve(1);
            DecodeBudget budget(budget_bytes);
            ASSERT_FALSE(recorder.fetch(&key, 1, budget, &out, &costs, &refused, name_bytes));
            const FetchRefusal &why = recorder.refusal();
            EXPECT_EQ(FetchRefusal::NAMES, why.cause) << budget_bytes;
            EXPECT_EQ(c[0].demand, why.demand);
            EXPECT_EQ(most, why.labels);
            EXPECT_EQ(names, why.names_bytes);
            EXPECT_TRUE(recorder.labels().empty());
        }
        // one byte below the demand: the row itself
        LabelRecorder recorder(oracle, kind, 64);
        std::vector<LabelRecorder::NodeLabels> out;
        std::vector<KeyCost> costs;
        out.reserve(1);
        costs.reserve(1);
        DecodeBudget budget(c[0].demand - 1);
        ASSERT_FALSE(recorder.fetch(&key, 1, budget, &out, &costs, &refused, name_bytes));
        EXPECT_NE(FetchRefusal::NAMES, recorder.refusal().cause);
    }
}

// F1: what naming labels provisionally holds is bounded by the per-label charge, for any
// number of labels (the pending table and list, their growth transients included)
TEST(LabelOracleBudgetedRecorder, NamingIsChargedPerLabel) {
    for (uint64_t m = 1; m < 200000; m = m < 300 ? m + 1 : m * 17 / 16) {
        EXPECT_LE(LabelRecorder::pending_bytes(m), m * LabelRecorder::kNamingBytes) << m;
    }
    EXPECT_EQ(LabelRecorder::kNamingBytes, LabelRecorder::pending_bytes(1));
}

// F1: jemalloc's peak of a budget-aware fetch that names many labels stays within what the
// call charged (with names priced at the walker's model of a dictionary label without its
// delivery), up to the dictionary map's first bucket array, a constant per recorder: the
// pending labels were a map per label (~1.5 KB each) held uncharged
TEST(LabelOracleBudgetedRecorder, NamingWithinTheAllocatorsPeak) {
#if USE_JEMALLOC
    for (LabelKind kind : { LabelKind::COLUMN, LabelKind::HEADER }) {
        // 400 records sharing one stretch: its k-mers carry 400 labels (columns, or the
        // headers of 400 columns)
        const std::string shared = sequences(1, 40, 77)[0];
        std::vector<std::string> seqs, labels;
        for (size_t i = 0; i < 400; ++i) {
            seqs.push_back(sequences(1, 10, 100 + i)[0] + shared + sequences(1, 10, 900 + i)[0]);
            labels.push_back("C" + std::to_string(i));
        }
        auto anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
                kK, seqs, labels, DeBruijnGraph::BASIC, true, std::vector<uint64_t>(400, 0));
        const auto &encoder = anno->get_annotator().get_label_encoder();
        std::vector<std::vector<std::string>> headers(encoder.size());
        std::vector<std::vector<uint64_t>> num_kmers(encoder.size());
        for (size_t i = 0; i < seqs.size(); ++i) {
            const size_t c = encoder.encode(labels[i]);
            headers[c].push_back("h" + std::to_string(i));
            num_kmers[c].push_back(seqs[i].size() - kK + 1);
        }
        annot::CoordToHeader cth(std::move(headers), std::move(num_kmers));
        LabelOracle oracle(*anno, &cth);
        std::vector<node_index> keys = oracle.keys_of_sequence(shared);
        keys.resize(8);
        // a dictionary label as the walker's account models it, its delivery apart
        auto name_bytes = [](std::string_view name) {
            return 2 * sizeof(LabelRef) + 4 * sizeof(LabelArmSummary) + 96 + 2 * name.size();
        };
        LabelRecorder recorder(oracle, kind, 1000);
        recorder.set_max_cache_bytes(0);
        std::vector<LabelRecorder::NodeLabels> out;
        std::vector<KeyCost> costs;
        out.reserve(keys.size());
        costs.reserve(keys.size());
        DecodeBudget budget;
        size_t refused = 0;
        size_t sz = sizeof(uint64_t);
        if (mallctl("thread.peak.reset", nullptr, nullptr, nullptr, 0))
            GTEST_SKIP() << "jemalloc without thread.peak";
        ASSERT_TRUE(recorder.fetch(keys.data(), keys.size(), budget, &out, &costs, &refused,
                                   name_bytes));
        uint64_t peak = 0;
        mallctl("thread.peak.read", &peak, &sz, nullptr, 0);
        if (budget.peak() > 4096 && !peak)
            GTEST_SKIP() << "jemalloc is not the process allocator";
        ASSERT_EQ(400u, recorder.labels().size());
        // the dictionary map's first bucket array (tsl: 63 buckets), once per recorder
        const uint64_t first_buckets = 2048;
        EXPECT_LE(peak, budget.peak() + first_buckets) << to_string(kind);
        std::cerr << to_string(kind) << ": jemalloc peak " << peak << ", charged peak "
                  << budget.peak() << std::endl;
    }
#else
    GTEST_SKIP() << "needs jemalloc";
#endif
}

// Q5: which indexes have the budget-aware reads
TEST(LabelOracleBudgeted, DecodeChargedTruthTable) {
    auto seqs = sequences(3, 60, 9);
    std::vector<std::string> labels { "A", "B", "C" };
    auto column = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(kK, seqs, labels, DeBruijnGraph::BASIC);
    auto row_diff = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(kK, seqs, labels, DeBruijnGraph::BASIC);
    auto row_diff_coord = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(kK, seqs, labels, DeBruijnGraph::BASIC, true);
    auto row_disk = test::build_anno_graph<DBGSuccinct, annot::RowDiffDiskAnnotator>(kK, seqs, labels, DeBruijnGraph::BASIC);
    auto row_flat = test::build_anno_graph<DBGSuccinct, annot::RowFlatAnnotator>(kK, seqs, labels, DeBruijnGraph::BASIC);
    EXPECT_FALSE(LabelOracle(*column).decode_charged());
    EXPECT_FALSE(LabelOracle(*column).row_diff());
    EXPECT_TRUE(LabelOracle(*row_diff).decode_charged());
    EXPECT_TRUE(LabelOracle(*row_diff_coord).decode_charged());
    EXPECT_FALSE(LabelOracle(*row_disk).decode_charged());
    EXPECT_TRUE(LabelOracle(*row_disk).row_diff());
    EXPECT_FALSE(LabelOracle(*row_flat).decode_charged());
}

// Runs |body| on a thread of its own and fails if it has not returned within |seconds|: a
// warm that never advances would otherwise hang the whole suite. On a timeout the thread is
// left behind (it owns what it uses), and the process ends with the suite.
template <class Body>
void returns_within(Body body, int seconds, const std::string &what) {
    auto done = std::make_shared<std::promise<void>>();
    std::future<void> finished = done->get_future();
    std::thread([done, body]() mutable {
        body();
        done->set_value();
    }).detach();
    if (finished.wait_for(std::chrono::seconds(seconds)) != std::future_status::ready)
        ADD_FAILURE() << what << " did not return within " << seconds << " s";
}

// Review of stage 3, F1: a budgeted warm of a cache whose capacity is zero returned never
// (its runs were at most the capacity long and did not advance); it returns at once, warms
// nothing and leaves the budget as it was, for the query and the recorder alike
TEST(LabelOracleBudgetedQuery, ZeroCacheWarmReturns) {
    for (bool coordinates : { false, true }) {
        auto fx = std::make_shared<Fixture>(coordinates);
        returns_within([fx, coordinates]() {
            LabelOracle oracle(*fx->anno, fx->cth.get());
            for (const auto &[names, with_coords] : label_sets(coordinates)) {
                LabelQuery query(oracle, refs(oracle, names), with_coords,
                                 LabelOracle::Access::AUTO, 0);
                std::vector<node_index> keys = oracle.keys_of_sequence(fx->seqs[0]);
                DecodeBudget budget;
                query.warm(keys, budget);
                EXPECT_EQ(0u, budget.held());
                EXPECT_EQ(0u, query.cache_bytes());
                // and a fetch still answers
                std::vector<LabelQuery::NodeHits> hits;
                std::vector<KeyCost> costs;
                hits.reserve(keys.size());
                costs.reserve(keys.size());
                DecodeBudget fetch_budget;
                size_t refused = 0;
                EXPECT_TRUE(query.fetch(keys.data(), keys.size(), fetch_budget, &hits, &costs,
                                        &refused));
                EXPECT_EQ(keys.size(), hits.size());
            }
        }, 30, "LabelQuery::warm with a cache of capacity zero");
    }
}

TEST(LabelOracleBudgetedRecorder, ZeroCacheWarmReturns) {
    for (bool coordinates : { false, true }) {
        auto fx = std::make_shared<Fixture>(coordinates);
        returns_within([fx, coordinates]() {
            LabelOracle oracle(*fx->anno, fx->cth.get());
            std::vector<LabelKind> kinds { LabelKind::COLUMN };
            if (coordinates)
                kinds.push_back(LabelKind::HEADER);
            for (LabelKind kind : kinds) {
                LabelRecorder recorder(oracle, kind, 64, 0);
                std::vector<node_index> keys = oracle.keys_of_sequence(fx->seqs[0]);
                DecodeBudget budget;
                recorder.warm(keys, budget);
                EXPECT_EQ(0u, budget.held());
                EXPECT_EQ(0u, recorder.cache_bytes());
                std::vector<LabelRecorder::NodeLabels> lists;
                std::vector<KeyCost> costs;
                lists.reserve(keys.size());
                costs.reserve(keys.size());
                DecodeBudget fetch_budget;
                size_t refused = 0;
                EXPECT_TRUE(recorder.fetch(keys.data(), keys.size(), fetch_budget, &lists, &costs,
                                           &refused, [](std::string_view n) { return 100 + n.size(); }));
                EXPECT_EQ(keys.size(), lists.size());
            }
        }, 30, "LabelRecorder::warm with a cache of capacity zero");
    }
}

// Review of stage 3, F2: the ordinary and the budget-aware path evict through one helper,
// rows and costs together. With a cache of one row: budgeted fetch A, ordinary fetch (or
// warm) B, budgeted fetch (or warm) B — the recorder threw std::out_of_range for B, whose row
// was cached without its cost beside A's stale cost. Each answer equals the unbudgeted one.
TEST(LabelOracleBudgetedQuery, AlternatingFetchPathsKeepCosts) {
    Fixture fx(true);
    LabelOracle oracle(*fx.anno, fx.cth.get());
    const std::vector<node_index> keys = oracle.keys_of_sequence(fx.seqs[0]);
    ASSERT_GT(keys.size(), 2u);
    ASSERT_NE(keys[0], keys[1]);
    for (const auto &[names, with_coords] : label_sets(true)) {
        LabelQuery plain(oracle, refs(oracle, names), with_coords);
        for (int variant = 0; variant < 4; ++variant) {
            LabelQuery query(oracle, refs(oracle, names), with_coords,
                             LabelOracle::Access::AUTO, 1);
            std::vector<LabelQuery::NodeHits> hits;
            std::vector<KeyCost> costs;
            hits.reserve(8);
            costs.reserve(8);
            size_t refused = 0;
            DecodeBudget first;
            ASSERT_TRUE(query.fetch(&keys[0], 1, first, &hits, &costs, &refused));
            if (variant % 2) {
                query.warm(std::vector<node_index>{ keys[1] });
            } else {
                EXPECT_EQ(plain.fetch(std::vector<node_index>{ keys[1] }),
                          query.fetch(std::vector<node_index>{ keys[1] }));
            }
            if (variant >= 2) {
                DecodeBudget warm;
                EXPECT_NO_THROW(query.warm(std::vector<node_index>{ keys[1], keys[2] }, warm));
            }
            DecodeBudget next;
            EXPECT_NO_THROW(EXPECT_TRUE(query.fetch(&keys[1], 1, next, &hits, &costs, &refused)))
                << names[0] << " variant " << variant;
            ASSERT_EQ(2u, hits.size());
            EXPECT_EQ(plain.fetch(std::vector<node_index>{ keys[1] })[0], hits[1]) << names[0] << " variant " << variant;
        }
    }
}

TEST(LabelOracleBudgetedRecorder, AlternatingFetchPathsKeepCosts) {
    Fixture fx(true);
    LabelOracle oracle(*fx.anno, fx.cth.get());
    const std::vector<node_index> keys = oracle.keys_of_sequence(fx.seqs[0]);
    ASSERT_GT(keys.size(), 2u);
    ASSERT_NE(keys[0], keys[1]);
    auto name_bytes = [](std::string_view n) { return 100 + n.size(); };
    for (LabelKind kind : { LabelKind::COLUMN, LabelKind::HEADER }) {
        for (int variant = 0; variant < 4; ++variant) {
            LabelRecorder plain(oracle, kind, 64);
            LabelRecorder recorder(oracle, kind, 64, 1);
            std::vector<LabelRecorder::NodeLabels> lists;
            std::vector<KeyCost> costs;
            lists.reserve(8);
            costs.reserve(8);
            size_t refused = 0;
            DecodeBudget first;
            ASSERT_TRUE(recorder.fetch(&keys[0], 1, first, &lists, &costs, &refused, name_bytes));
            if (variant % 2) {
                recorder.warm(std::vector<node_index>{ keys[1] });
            } else {
                recorder.fetch(std::vector<node_index>{ keys[1] });
            }
            if (variant >= 2) {
                DecodeBudget warm;
                EXPECT_NO_THROW(recorder.warm(std::vector<node_index>{ keys[1], keys[2] }, warm));
            }
            DecodeBudget next;
            EXPECT_NO_THROW(EXPECT_TRUE(recorder.fetch(&keys[1], 1, next, &lists, &costs,
                                                       &refused, name_bytes)))
                << static_cast<int>(kind) << " variant " << variant;
            ASSERT_EQ(2u, lists.size());
            // the same labels as an unbudgeted read of the same keys in the same order
            const auto expected = plain.fetch(std::vector<node_index>{ keys[0], keys[1] });
            EXPECT_EQ(expected[1].total, lists[1].total);
            ASSERT_EQ(expected[1].labels.size(), lists[1].labels.size());
            for (size_t i = 0; i < lists[1].labels.size(); ++i) {
                EXPECT_EQ(plain.labels()[expected[1].labels[i]].name,
                          recorder.labels()[lists[1].labels[i]].name);
            }
        }
    }
}


/************* pass 5: the chunked deadlines (LabelOracle::pacer, ReadPacing) *************/

// the pacer: a read the deadline cannot fall into is one piece (no deadline; predicted at the
// slowest rate seen to take less than 1/64 of the time left; its rest, at its own previous
// chunk's rate, less than 1/4); a read that may reach the deadline starts with at most
// first_rows (8) rows whatever the rate known, then chunks sized to the target (or the time
// left) at the previous chunk's rate, growing at most 4 times per chunk; off at target 0
TEST(LabelOraclePacing, PacerSizesChunksByTime) {
    const double inf = std::numeric_limits<double>::infinity();
    DecodePacer p;
    EXPECT_EQ(1000u, p.next(1000, 5, 0, 0));       // off: one piece
    p.target_ms = 50;
    EXPECT_EQ(1000u, p.next(1000, inf, 0, 0));     // no deadline: one piece
    EXPECT_EQ(8u, p.next(1000, 1e9, 0, 0));        // no rate known yet: a first chunk measures
    EXPECT_EQ(5u, p.next(5, 1e9, 0, 0));
    p.record(64, 6.4);                             // 0.1 ms per row
    EXPECT_NEAR(0.1, p.ms_per_row, 1e-12);
    // far from the deadline: 100,000 rows predicted at 10 s, x 64 = 640 s, less than the time left
    EXPECT_EQ(100000u, p.next(100000, 1e9, 0, 0));
    EXPECT_EQ(100000u, p.next(100000, 640000.1, 0, 0));
    // near it: a first chunk of at most 8 rows, though 500 would take the target at 0.1 ms
    EXPECT_EQ(8u, p.next(100000, 639999.9, 0, 0));
    // then sized at the previous chunk's own rate, at most 4 x its rows
    EXPECT_EQ(32u, p.next(100000, 1000, 8, 0.8));  // 0.1 ms/row: 500 rows, at most 32
    EXPECT_EQ(25u, p.next(100000, 1000, 8, 16));   // 2 ms/row (slower than any before): 25
    EXPECT_EQ(5u, p.next(100000, 5, 8, 8));        // 5 ms left at 1 ms/row
    // the rest predicted at its own rate (0.1 ms/row: 200 ms) x 4 below the time left: one piece
    EXPECT_EQ(2000u, p.next(2000, 800.1, 64, 6.4));
    EXPECT_EQ(256u, p.next(2000, 799.9, 64, 6.4)); // else 50 ms at 0.1 ms/row, at most 4 x 64
    EXPECT_EQ(1u, p.next(100000, -5, 0, 0));       // past the deadline: one row, the stop decides
    EXPECT_EQ(1u, p.next(100000, -5, 8, 0.8));
    p.record(10, 10);                              // a slower piece raises the rate at once
    EXPECT_NEAR(1.0, p.ms_per_row, 1e-12);
    p.record(100, 1);                              // a faster one does not lower it
    EXPECT_NEAR(1.0, p.ms_per_row, 1e-12);
    EXPECT_NEAR(10.0, p.max_read_ms, 1e-12);       // the longest piece
    // the tests' mode: infinite factors split every read into the smallest chunks
    DecodePacer t;
    t.target_ms = 1e-9;
    t.first_rows = 1;
    t.far_factor = t.rest_factor = inf;
    EXPECT_EQ(1u, t.next(1000, inf, 0, 0));
    EXPECT_EQ(1u, t.next(999, inf, 1, 0.5));
    EXPECT_EQ(4u, t.next(999, inf, 1, 0));         // a chunk too fast to measure: grown 4 x
}

namespace {

// a paced read that never stops, in chunks of one row once a rate is known
ReadPacing never_stop() {
    ReadPacing p;
    p.ms_left = []() { return std::numeric_limits<double>::infinity(); };
    p.stop = []() { return false; };
    return p;
}

void expect_same_counters(const LabelOracle::Counters &a, const LabelOracle::Counters &b,
                          const std::string &what) {
    // (keys_mapped is not the reads': the test maps its keys with the first oracle)
    EXPECT_EQ(a.rows_requested, b.rows_requested) << what;
    EXPECT_EQ(a.cache_hits, b.cache_hits) << what;
    EXPECT_EQ(a.rows_fetched, b.rows_fetched) << what;
    EXPECT_EQ(a.tuple_rows_fetched, b.tuple_rows_fetched) << what;
    EXPECT_EQ(a.direct_reads, b.direct_reads) << what;
    EXPECT_EQ(a.coords_mapped, b.coords_mapped) << what;
}

// an oracle whose reads are paced in chunks of one row (a rate far above the target), every
// read split, deadline or none
void one_row_chunks(LabelOracle &oracle) {
    oracle.pacer().target_ms = 1e-9;
    oracle.pacer().first_rows = 1;
    oracle.pacer().far_factor = std::numeric_limits<double>::infinity();
    oracle.pacer().rest_factor = std::numeric_limits<double>::infinity();
}

} // namespace

// A paced read answers exactly as one read: the same hits, labels and ids, and the same
// counters, cache and per-call bytes, on every path (direct, rows, tuples, budget-aware), with
// a small cache forcing evictions, in chunks of one row
TEST(LabelOraclePacing, PacedReadsAnswerAsOneRead) {
    for (bool coordinates : { false, true }) {
        Fixture fx(coordinates);
        // the direct path needs an annotation with single-cell reads
        auto column = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                kK, fx.seqs, fx.labels);
        for (const auto &[names, with_coords] : label_sets(coordinates)) {
            for (int variant = 0; variant < 3; ++variant) {
                const AnnotatedDBG &index = variant == 2 ? *column : *fx.anno;
                if (variant == 2 && (coordinates || with_coords))
                    continue;
                LabelOracle a(index, variant == 2 ? nullptr : fx.cth.get());
                LabelOracle b(index, variant == 2 ? nullptr : fx.cth.get());
                one_row_chunks(b);
                const auto access = variant == 2 ? LabelOracle::Access::DIRECT
                                  : variant == 1 && !with_coords ? LabelOracle::Access::ROWS
                                                                 : LabelOracle::Access::AUTO;
                bool has_headers = false;
                for (const auto &n : names) {
                    has_headers |= n.rfind("acc", 0) == 0 || n.find("_h") != std::string::npos;
                }
                if (variant == 2 && has_headers)
                    continue;
                LabelQuery qa(a, refs(a, names), with_coords, access, 10);
                LabelQuery qb(b, refs(b, names), with_coords, access, 10);
                qa.set_max_cache_bytes(1 << 12);
                qb.set_max_cache_bytes(1 << 12);
                const std::string what = names[0] + " variant " + std::to_string(variant)
                                       + " path " + qa.access_path();
                const auto batches = fx.batches(a);
                for (size_t i = 0; i < batches.size(); ++i) {
                    ReadPacing pacing = never_stop();
                    ASSERT_EQ(qa.fetch(batches[i]), qb.fetch(batches[i], &pacing)) << what;
                    EXPECT_FALSE(pacing.interrupted);
                    EXPECT_EQ(qa.last_call_bytes(), qb.last_call_bytes()) << what;
                    EXPECT_EQ(qa.cache_bytes(), qb.cache_bytes()) << what;
                    if (i + 1 < batches.size()) {
                        qa.warm(batches[i + 1]);
                        qb.warm(batches[i + 1], &pacing);
                        EXPECT_EQ(qa.cache_bytes(), qb.cache_bytes()) << what;
                    }
                }
                expect_same_counters(a.counters(), b.counters(), what);
                EXPECT_GT(b.pacer().max_read_ms, 0);
            }
        }
        // the budget-aware path (row-diff): the same hits, costs, held bytes and counters
        for (const auto &[names, with_coords] : label_sets(coordinates)) {
            LabelOracle a(*fx.anno, fx.cth.get());
            LabelOracle b(*fx.anno, fx.cth.get());
            one_row_chunks(b);
            LabelQuery qa(a, refs(a, names), with_coords, LabelOracle::Access::AUTO, 10);
            LabelQuery qb(b, refs(b, names), with_coords, LabelOracle::Access::AUTO, 10);
            qa.set_max_cache_bytes(1 << 12);
            qb.set_max_cache_bytes(1 << 12);
            for (const auto &batch : fx.batches(a)) {
                std::vector<LabelQuery::NodeHits> ha, hb;
                std::vector<KeyCost> ca, cb;
                ha.reserve(batch.size());
                hb.reserve(batch.size());
                ca.reserve(batch.size());
                cb.reserve(batch.size());
                DecodeBudget ba, bb;
                size_t ra = 0, rb = 0;
                ReadPacing pacing = never_stop();
                ASSERT_TRUE(qa.fetch(batch.data(), batch.size(), ba, &ha, &ca, &ra));
                ASSERT_TRUE(qb.fetch(batch.data(), batch.size(), bb, &hb, &cb, &rb, &pacing));
                ASSERT_EQ(ha, hb) << names[0];
                ASSERT_EQ(ca.size(), cb.size());
                for (size_t i = 0; i < ca.size(); ++i) {
                    EXPECT_EQ(ca[i].demand, cb[i].demand);
                    EXPECT_EQ(ca[i].dependency_units, cb[i].dependency_units);
                }
                EXPECT_EQ(ba.held(), bb.held());
                EXPECT_EQ(qa.cache_bytes(), qb.cache_bytes());
                DecodeBudget wa, wb;
                qa.warm(batch, wa);
                qb.warm(batch, wb, &pacing);
                EXPECT_EQ(qa.cache_bytes(), qb.cache_bytes());
            }
            expect_same_counters(a.counters(), b.counters(), names[0] + " budgeted");
        }
    }
}

// LabelRecorder: the same lists, the same dictionary in the same order, the same counters,
// unbudgeted and budget-aware, in chunks of one row
TEST(LabelOraclePacing, PacedRecordingAnswersAsOneRead) {
    Fixture fx(true);
    auto name_bytes = [](std::string_view n) { return 100 + n.size(); };
    for (LabelKind kind : { LabelKind::COLUMN, LabelKind::HEADER }) {
        for (bool budgeted : { false, true }) {
            LabelOracle a(*fx.anno, fx.cth.get());
            LabelOracle b(*fx.anno, fx.cth.get());
            one_row_chunks(b);
            LabelRecorder ra(a, kind, 3, 10);
            LabelRecorder rb(b, kind, 3, 10);
            const auto batches = fx.batches(a);
            for (size_t i = 0; i < batches.size(); ++i) {
                const auto &batch = batches[i];
                ReadPacing pacing = never_stop();
                if (!budgeted) {
                    const auto la = ra.fetch(batch);
                    const auto lb = rb.fetch(batch, &pacing);
                    ASSERT_EQ(la.size(), lb.size());
                    for (size_t j = 0; j < la.size(); ++j) {
                        EXPECT_EQ(la[j].labels, lb[j].labels);
                        EXPECT_EQ(la[j].total, lb[j].total);
                    }
                    if (i + 1 < batches.size()) {
                        ra.warm(batches[i + 1]);
                        rb.warm(batches[i + 1], &pacing);
                    }
                } else {
                    std::vector<LabelRecorder::NodeLabels> la, lb;
                    std::vector<KeyCost> ca, cb;
                    la.reserve(batch.size());
                    lb.reserve(batch.size());
                    ca.reserve(batch.size());
                    cb.reserve(batch.size());
                    DecodeBudget ba, bb;
                    size_t xa = 0, xb = 0;
                    ASSERT_TRUE(ra.fetch(batch.data(), batch.size(), ba, &la, &ca, &xa, name_bytes));
                    ASSERT_TRUE(rb.fetch(batch.data(), batch.size(), bb, &lb, &cb, &xb, name_bytes,
                                         &pacing));
                    for (size_t j = 0; j < la.size(); ++j) {
                        EXPECT_EQ(la[j].labels, lb[j].labels);
                        EXPECT_EQ(la[j].total, lb[j].total);
                        EXPECT_EQ(ca[j].demand, cb[j].demand);
                    }
                    EXPECT_EQ(ba.held(), bb.held());
                    DecodeBudget wa, wb;
                    ra.warm(batch, wa);
                    rb.warm(batch, wb, &pacing);
                }
                EXPECT_EQ(ra.cache_bytes(), rb.cache_bytes());
            }
            ASSERT_EQ(ra.labels().size(), rb.labels().size());
            for (size_t id = 0; id < ra.labels().size(); ++id) {
                EXPECT_EQ(ra.labels()[id].name, rb.labels()[id].name);
            }
            expect_same_counters(a.counters(), b.counters(), budgeted ? "budgeted" : "plain");
        }
    }
}

// An interrupted read returns nothing and changes no counter, no result and (budget-aware)
// neither the cache nor the budget; it states the work its finished chunks decoded, and the
// same read run again answers as one read
TEST(LabelOraclePacing, InterruptedReadsChangeNothing) {
    Fixture fx(true);
    const auto names = label_sets(true)[1].first;    // F, L2, L6 with coordinates
    LabelOracle a(*fx.anno, fx.cth.get());
    LabelOracle b(*fx.anno, fx.cth.get());
    one_row_chunks(b);
    LabelQuery qa(a, refs(a, names), true);
    LabelQuery qb(b, refs(b, names), true);
    std::vector<node_index> keys = a.keys_of_sequence(fx.seqs[0]);
    ASSERT_GT(keys.size(), 10u);
    // stopped before the third chunk
    auto stop_at = [](int n) {
        ReadPacing p = never_stop();
        auto calls = std::make_shared<int>(0);
        p.stop = [calls, n]() { return ++*calls >= n; };
        return p;
    };
    const LabelOracle::Counters before = b.counters();
    ReadPacing p = stop_at(3);
    EXPECT_TRUE(qb.fetch(keys, &p).empty());
    EXPECT_TRUE(p.interrupted);
    EXPECT_GE(p.units, 2 * 8u);                     // two rows decoded, each at least its key
    EXPECT_EQ(before.rows_requested, b.counters().rows_requested);
    EXPECT_EQ(before.cache_hits, b.counters().cache_hits);
    // the same call again, unstopped: what one read answers (the decoded chunks were cached,
    // so this call's own hits differ, as a lookahead's warming makes them differ)
    ReadPacing go = never_stop();
    EXPECT_EQ(qa.fetch(keys), qb.fetch(keys, &go));
    EXPECT_EQ(a.counters().rows_requested, b.counters().rows_requested);

    // budget-aware: refused whole with cause INTERRUPTED, budget and cache as on entry
    LabelQuery qc(b, refs(b, names), true);
    std::vector<LabelQuery::NodeHits> hits;
    std::vector<KeyCost> costs;
    hits.reserve(keys.size());
    costs.reserve(keys.size());
    DecodeBudget budget;
    ASSERT_TRUE(budget.charge(123));
    const uint64_t cache_before = qc.cache_bytes();
    const LabelOracle::Counters c_before = b.counters();
    size_t refused = 0;
    ReadPacing q = stop_at(3);
    EXPECT_FALSE(qc.fetch(keys.data(), keys.size(), budget, &hits, &costs, &refused, &q));
    EXPECT_TRUE(q.interrupted);
    EXPECT_EQ(FetchRefusal::INTERRUPTED, qc.refusal().cause);
    EXPECT_GT(q.units, 0u);
    EXPECT_TRUE(hits.empty());
    EXPECT_TRUE(costs.empty());
    EXPECT_EQ(123u, budget.held());
    EXPECT_EQ(cache_before, qc.cache_bytes());
    EXPECT_EQ(c_before.rows_requested, b.counters().rows_requested);
    EXPECT_EQ(c_before.cache_hits, b.counters().cache_hits);

    // a recorder names nothing when interrupted, by either path
    for (bool budgeted : { false, true }) {
        LabelRecorder rec(b, LabelKind::HEADER, 64);
        ReadPacing r = stop_at(3);
        if (!budgeted) {
            EXPECT_TRUE(rec.fetch(keys, &r).empty());
        } else {
            std::vector<LabelRecorder::NodeLabels> lists;
            std::vector<KeyCost> lc;
            lists.reserve(keys.size());
            lc.reserve(keys.size());
            DecodeBudget lb;
            EXPECT_FALSE(rec.fetch(keys.data(), keys.size(), lb, &lists, &lc, &refused,
                                   [](std::string_view) { return uint64_t(1); }, &r));
            EXPECT_EQ(FetchRefusal::INTERRUPTED, rec.refusal().cause);
            EXPECT_EQ(0u, lb.held());
        }
        EXPECT_TRUE(r.interrupted);
        EXPECT_GE(r.units, 16u);
        EXPECT_TRUE(rec.labels().empty());
    }
    // a warming stopped ends silently
    LabelQuery qw(b, refs(b, names), true);
    ReadPacing w = stop_at(2);
    qw.warm(keys, &w);
    EXPECT_TRUE(w.interrupted);
}

// Review of pass 5, finding 5 (the reviewer's pacing_probe.cpp): the budget-aware lookahead
// (warm) cut its runs from the globally sorted missing keys, which put rows of as many paths
// as rows into each run. Its runs are taken in the walk's order now, each piece sorted: on 32
// independent paths of 64 walked rows each, in runs of 8 (a cache of 8 keys), the decoder's
// charges are those of decoding the walk-ordered runs (28,736 in the review's measurement)
// and a third or less of those of the sorted runs (161,679). Counted, not timed
TEST(LabelOracleBudgeted, LookaheadRunsFollowTheWalk) {
    std::mt19937 gen(1282);
    std::vector<std::string> seqs, names;
    for (size_t i = 0; i < 32; ++i) {
        std::string s(240, 'A');
        for (char &c : s) {
            c = "ACGT"[gen() % 4];
        }
        seqs.push_back(s);
        names.push_back("p" + std::to_string(i));
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(21, seqs, names);
    LabelOracle oracle(*anno);
    ASSERT_TRUE(oracle.decode_charged());
    // the walk: 64 consecutive nodes of each path, path after path
    std::vector<node_index> walking;
    for (const std::string &seq : seqs) {
        const auto keys = oracle.keys_of_sequence(seq);
        walking.insert(walking.end(), keys.begin() + 70, keys.begin() + 134);
    }
    std::vector<node_index> sorted = walking;
    std::sort(sorted.begin(), sorted.end());
    ASSERT_EQ(walking.size(), std::set<node_index>(walking.begin(), walking.end()).size());
    // the decoder's charges for runs of 8 cut from |order|, each run sorted
    auto charges_of = [&](const std::vector<node_index> &order) {
        uint64_t charges = 0;
        for (size_t begin = 0; begin < order.size(); begin += 8) {
            std::vector<Row> chunk;
            for (size_t i = begin; i < std::min(begin + 8, order.size()); ++i) {
                chunk.push_back(AnnotatedDBG::graph_to_anno_index(order[i]));
            }
            std::sort(chunk.begin(), chunk.end());
            DecodeBudget budget;
            std::vector<annot::matrix::BinaryMatrix::SetBitPositions> out;
            std::vector<annot::matrix::RowCost> costs;
            std::vector<uint64_t> held;
            EXPECT_EQ(annot::matrix::DecodeStatus::OK,
                      oracle.get_rows(chunk, budget, &out, &costs, &held));
            charges += budget.charges();
        }
        return charges;
    };
    const uint64_t walk_order = charges_of(walking), sorted_order = charges_of(sorted);
    EXPECT_GT(sorted_order, 3 * walk_order);
    std::vector<LabelRef> labels;
    for (const std::string &name : names) {
        labels.push_back(oracle.resolve_label(name));
    }
    // the warms' own charges: the decoder's for the walk-ordered runs, plus each run's
    // containers and the hits built (a few per run), far below the sorted runs'
    for (bool recorder : { false, true }) {
        DecodeBudget budget;
        if (recorder) {
            LabelRecorder rec(oracle, LabelKind::COLUMN, 64, 8);
            rec.warm(walking, budget, nullptr);
        } else {
            LabelQuery query(oracle, labels, false, LabelOracle::Access::ROWS, 8);
            query.warm(walking, budget, nullptr);
        }
        std::cerr << (recorder ? "recorder" : "query") << " warm: " << budget.charges()
                  << " decoder charges; runs cut from the walk " << walk_order
                  << ", from the sorted keys " << sorted_order << std::endl;
        EXPECT_GE(budget.charges(), walk_order) << recorder;
        EXPECT_LT(budget.charges(), walk_order + 64 * (walking.size() / 8)) << recorder;
        EXPECT_LT(3 * budget.charges(), sorted_order) << recorder;
    }
}

} // namespace
