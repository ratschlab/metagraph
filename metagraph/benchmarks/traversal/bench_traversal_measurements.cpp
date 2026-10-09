// Measurements of the labeled traversal and of the pattern search, kept out of the unit tests:
// they print numbers and check no contract the unit tests do not. The depth at a memory stop
// with and without record coordinates and what completing takes (M1), the coordinate share of a
// wide fixture's delivery reserve (M2) -- the tables of DESIGN-traverse-graphlet.md §26.5 --,
// the speed of the coordinate mappers, and the dummy-fraction sampler on a large graph. Each
// runs once and prints its table to stderr; a consistency check that fails ends it with an
// error.
//
//     cd metagraph/build && make benchmarks
//     ./benchmarks --benchmark_filter=Traversal
//
// The mini index (scripts/traversal/build_mini_refseq.sh) is read from $MINI_REFSEQ (default:
// mini_refseq in the working directory), the fixtures of scripts/traversal/
// make_column_coord_fixtures.sh and make_wide_coord_fixture.sh from the working directory.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <benchmark/benchmark.h>

#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "cli/server_checks.hpp"
#include "cli/traverse.hpp"
#include "cli/traverse_attempts.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/unix_tools.hpp"
#include "graph/alignment/pattern_search.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "graph/traversal/walker.hpp"


namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;

// a failed consistency check ends the measurement with an error
#define MEASURE_CHECK(state, condition)                                                    \
    do {                                                                                 \
        if (!(condition)) {                                                              \
            (state).SkipWithError("check failed: " #condition);                         \
            return;                                                                      \
        }                                                                                \
    } while (false)

// the mini index, as the MiniRefSeq unit tests load it, and its blaNDM-1 query
struct Mini {
    std::unique_ptr<AnnotatedDBG> anno;
    std::unique_ptr<LabelOracle> oracle;
    std::string query;
};

// loaded once; nullptr, the benchmark ended with an error, when it is not there
const Mini* mini_index(benchmark::State &state) {
    static std::unique_ptr<Mini> mini;
    static bool tried = false;
    if (!tried) {
        tried = true;
        const char *env = std::getenv("MINI_REFSEQ");
        const std::string dir = env ? env : "mini_refseq";
        const std::string base = dir + "/annotation.relaxed.relabeled";
        auto graph = std::make_shared<DBGSuccinct>(2);
        auto annotation = std::make_unique<annot::RowDiffBRWTCoordAnnotator>();
        auto cth = std::make_unique<annot::CoordToHeader>();
        std::ifstream fasta(dir + "/sanity/ndm1.fa");
        if (graph->load(dir + "/graph_k31.dbg")
                && annotation->load(base + annot::RowDiffBRWTCoordAnnotator::kExtension)
                && cth->load(base) && fasta.good()) {
            using Matrix = annot::RowDiffBRWTCoordAnnotator::binary_matrix_type;
            const_cast<Matrix&>(annotation->get_matrix()).set_graph(graph.get());
            auto m = std::make_unique<Mini>();
            m->anno = std::make_unique<AnnotatedDBG>(std::move(graph), std::move(annotation),
                                                     false, std::move(cth));
            m->oracle = std::make_unique<LabelOracle>(*m->anno);
            std::string line;
            while (std::getline(fasta, line)) {
                if (!line.empty() && line[0] != '>')
                    m->query += line;
            }
            if (m->anno->check_compatibility() && !m->query.empty())
                mini = std::move(m);
        }
    }
    if (!mini)
        state.SkipWithError("the mini index is not in $MINI_REFSEQ or ./mini_refseq "
                            "(scripts/traversal/build_mini_refseq.sh)");
    return mini.get();
}


// ---------------------------------------------------------------- the coordinate mappers

// a column of |num_seqs| sequences of random lengths, as test_label_oracle.cpp builds it
struct CoordColumn {
    std::vector<uint64_t> lengths;
    std::vector<uint64_t> starts;
    std::unique_ptr<annot::CoordToHeader> cth;
};

CoordColumn coord_column(size_t num_seqs, uint64_t min_len, uint64_t spread, uint64_t seed) {
    std::mt19937_64 rng(seed);
    CoordColumn col;
    std::vector<std::string> names(num_seqs);
    for (size_t i = 0; i < num_seqs; ++i) {
        names[i] = "s" + std::to_string(i);
        col.lengths.push_back(min_len + (spread ? rng() % spread : 0));
        col.starts.push_back(i ? col.starts.back() + col.lengths[i - 1] : 0);
    }
    std::vector<uint64_t> lengths = col.lengths;
    col.cth = std::make_unique<annot::CoordToHeader>(
            std::vector<std::vector<std::string>>{ names },
            std::vector<std::vector<uint64_t>>{ lengths });
    return col;
}

// the coordinates of the patterns of a tuple row: |kind| 0 one per sequence in consecutive
// sequences, 1 runs of 7 in every third sequence, 2 one in scattered sequences
std::vector<uint64_t> coord_pattern(const CoordColumn &col, int kind, std::mt19937_64 &rng) {
    std::vector<uint64_t> coords;
    for (size_t i = 0; i < col.lengths.size(); ++i) {
        if (kind == 0 && i % 10 != 9) {
            coords.push_back(col.starts[i] + rng() % col.lengths[i]);
        } else if (kind == 1 && i % 3 == 0) {
            for (int j = 0; j < 7; ++j) {
                coords.push_back(col.starts[i] + rng() % col.lengths[i]);
            }
        } else if (kind == 2 && rng() % 20 == 0) {
            coords.push_back(col.starts[i] + rng() % col.lengths[i]);
        }
    }
    std::sort(coords.begin(), coords.end());
    return coords;
}


// The run mapper (LabelOracle::CoordRuns), on a column of 2,000 sequences of 100,000
// to 300,000 k-mers: map_single_coord per coordinate, the sequence ranges used always, and
// the run mapper, for the three patterns
void coord_runs(benchmark::State &state) {
    std::mt19937_64 rng(1);
    const CoordColumn col = coord_column(2000, 100'000, 200'000, 1);
    for (int kind : { 0, 1, 2 }) {
        const std::vector<uint64_t> coords = coord_pattern(col, kind, rng);
        const int reps = 2000;
        uint64_t sums[3] = { 0, 0, 0 };
        double ns[3];
        for (int method = 0; method < 3; ++method) {
            const auto t0 = std::chrono::steady_clock::now();
            for (int rep = 0; rep < reps; ++rep) {
                if (method == 0) {
                    for (uint64_t c : coords) {
                        const auto [seq, local] = col.cth->map_single_coord(0, c);
                        sums[0] += seq + local;
                    }
                } else if (method == 1) {
                    annot::CoordToHeader::SequenceRange range { 0, 1, 0 };
                    for (uint64_t c : coords) {
                        if (c < range.first || c > range.last)
                            range = col.cth->sequence_range(0, c);
                        sums[1] += range.seq_id + (c - range.first);
                    }
                } else {
                    LabelOracle::CoordRuns runs(*col.cth, 0, coords.size());
                    for (uint64_t c : coords) {
                        const auto [seq, local] = runs.map(c);
                        sums[2] += seq + local;
                    }
                }
            }
            ns[method] = std::chrono::duration<double, std::nano>(
                    std::chrono::steady_clock::now() - t0).count() / (double(coords.size()) * reps);
        }
        MEASURE_CHECK(state, sums[0] == sums[1]);
        MEASURE_CHECK(state, sums[0] == sums[2]);
        std::cerr << "pattern " << kind << " (" << coords.size() << " coordinates): map_single_coord "
                  << ns[0] << " ns, sequence ranges " << ns[1] << " ns, run mapper " << ns[2]
                  << " ns a coordinate" << std::endl;
    }
}

// The coordinate mapping by runs (LabelOracle::map_coords) against
// one map_single_coord per coordinate, on the tuple rows of the blaNDM query (same results,
// checked; the time printed)
void coord_mapping(benchmark::State &state, const Mini &mini) {
    std::vector<Row> rows;
    for (node_index key : mini.oracle->keys_of_sequence(mini.query)) {
        if (key != npos)
            rows.push_back(AnnotatedDBG::graph_to_anno_index(key));
    }
    const auto tuples = mini.oracle->get_row_tuples(rows);
    uint64_t coords = 0;
    for (const auto &row : tuples) {
        for (const auto &entry : row) {
            coords += entry.second.size();
        }
    }
    MEASURE_CHECK(state, coords > 0u);
    const int reps = 200;
    uint64_t sum_a = 0, sum_b = 0;
    Timer timer;
    for (int rep = 0; rep < reps; ++rep) {
        for (const auto &row : tuples) {
            for (const auto &[c, cs] : row) {
                for (Coord coord : cs) {
                    const auto [seq, local] = mini.oracle->map_coord(c, coord);
                    sum_a += seq * 31 + local;
                }
            }
        }
    }
    const double per_coord_a = timer.elapsed() * 1e9 / (double(coords) * reps);
    timer.reset();
    for (int rep = 0; rep < reps; ++rep) {
        for (const auto &row : tuples) {
            for (const auto &[c, cs] : row) {
                mini.oracle->map_coords(c, cs.data(), cs.size(), [&](Coord, uint64_t seq, Coord local) {
                    sum_b += seq * 31 + local;
                });
            }
        }
    }
    const double per_coord_b = timer.elapsed() * 1e9 / (double(coords) * reps);
    MEASURE_CHECK(state, sum_a == sum_b);
    std::cerr << rows.size() << " rows, " << coords << " coordinates: map_single_coord "
              << per_coord_a << " ns a coordinate, map_coords " << per_coord_b
              << " ns (x" << per_coord_a / per_coord_b << ")" << std::endl;
}


// ---------------------------------------------------------------- M1 and M2

namespace {

// a seed's result as the server writes it (compact JSON; a graphlet with its body), the
// coordinates' share of that text, and the seconds building it took
struct Written {
    uint64_t text = 0, coordinate_text = 0;
    double seconds = 0;
};

Written write_result(const SeedResult &r, const Seed &seed, const Strategy &st,
                     const std::string &detail, const LabelOracle &oracle) {
    Timer timer;
    Json::Value rj = mtg::cli::seed_result_to_json(r, st, detail, false);
    if (detail == "graphlet") {
        mtg::cli::GraphletContext ctx;
        ctx.k = oracle.get_k();
        ctx.regime = to_string(oracle.regime());
        ctx.alphabet = oracle.graph().alphabet();
        ctx.identity.meta_fp = "*";
        size_t lines = 0;
        std::string text = mtg::cli::graphlet_text(r, seed, st, ctx, rj, &lines);
        rj["graphlet_bytes"] = Json::UInt64(text.size());
        rj["graphlet_lines"] = Json::UInt64(lines);
        rj["graphlet"] = std::move(text);
    }
    Written w;
    w.coordinate_text = mtg::cli::coordinates_text_bytes(rj);
    w.text = mtg::cli::json_text(rj, true).size();
    w.seconds = timer.elapsed();
    return w;
}

// the delivery costs the server sets for |st| (process_traverse_request)
void price(Strategy *st, const LabelChangeCost &cost, const std::string &detail) {
    st->delivery = mtg::cli::delivery_costs(detail, st->sequences, mtg::cli::mgt_float_width(*st, cost),
                                            st->coordinates ? mtg::cli::CoordinatesOutput::BLOCK
                                                            : mtg::cli::CoordinatesOutput::NONE);
}

double median(std::vector<double> v) {
    if (v.empty())
        return 0;
    std::sort(v.begin(), v.end());
    return v.size() % 2 ? v[v.size() / 2] : (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2;
}

} // namespace

namespace {

uint64_t requested_depth(const SeedResult &r) {
    uint64_t d = 0;
    for (const ArmResult &arm : r.arms) {
        if (arm.requested)
            d += arm.complete_to_bp;
    }
    return d;
}

constexpr uint64_t kMiB = uint64_t(1) << 20;

// The smallest whole-MiB budget (bounds.max_memory_mb) whose admitted account holds an
// unbudgeted walk whose account peaked at |peak|: a budget's caches' allotments are part of its
// account (memory_allotments), so the budget is above the peak
uint64_t completion_budget(uint64_t peak) {
    uint64_t mib = std::max<uint64_t>(1, (peak + kMiB - 1) / kMiB);
    while (peak + memory_allotments(mib * kMiB) > mib * kMiB) {
        ++mib;
    }
    return mib * kMiB;
}

// One cell of the depth gate (whether coordinates may be asked under a memory budget): a trace walk
// unbudgeted, and under memory budgets of 50% and 75% of its own opt-out peak (exact bytes),
// with and without coordinates (cap 16): the depth at the stop. That metric is relative to the
// opt-out walk and drops what fails at depth 0 without coordinates, so a cell also states what
// it takes to COMPLETE: the smallest whole-MiB budget that holds the walk with and without
// coordinates (completion_budget, checked by walking at it and a MiB below), and what the walk
// with coordinates does at the budget that completes it without them (complete, stopped, failed
// at depth 0). |project|: coordinates priced as on refseq33m's taxid columns, whose runs carry
// thousands of chains — every run's and every seed label's list at the cap of 16, whatever the
// fixture's chains (the walker adds 2 * sizeof(Coord) for each occurrence it really holds: at
// most 256 B a list over the projection)
struct GateCell {
    bool ok = true;
    uint64_t peak_off = 0, peak_on = 0;      // unbudgeted account peaks, exact bytes
    uint64_t need_off = 0, need_on = 0;      // completion budgets, whole MiB in bytes
    uint64_t depth_free = 0;
    uint64_t d_off[2] {}, d_on[2] {};        // depth at 50% and 75% of peak_off
    bool f_off[2] {}, f_on[2] {}, s_off[2] {}, s_on[2] {};
    int at_need_off = 0;                     // with coordinates at need_off: 0 complete, 1 stop, 2 failed
    uint64_t at_need_off_depth = 0;
    std::vector<uint64_t> chains;            // unbudgeted, with coordinates: each run's true count
    size_t lower = 0, switch_runs = 0, switch_lower = 0;
};

GateCell gate_cell(const AnnotatedDBG &anno, const Seed &seed, const Strategy &base,
                   const LabelChangeCost &cost, const std::string &detail, bool project) {
    GateCell c;
    auto walk = [&](bool coordinates, uint64_t budget, bool *failed) {
        Strategy st = base;
        st.coordinates = coordinates;
        st.max_memory_bytes = budget;
        price(&st, cost, detail);
        if (coordinates && project) {
            const uint64_t occurrence = 2 * sizeof(Coord) + st.delivery.occurrence;
            st.delivery.coordinate_run += 16 * occurrence;
            st.delivery.coordinate_seed += 16 * occurrence;
            st.delivery.occurrence = 0;
        }
        LabelOracle oracle(anno);
        *failed = false;
        try {
            return traverse_seed(oracle, seed, st, cost);
        } catch (const SeedBudgetError &) {
            *failed = true;
            return SeedResult();
        }
    };
    bool failed = false;
    const SeedResult free_off = walk(false, 0, &failed);
    if (failed) {
        c.ok = false;
        return c;
    }
    const SeedResult free_on = walk(true, 0, &failed);
    if (failed) {
        c.ok = false;
        return c;
    }
    c.peak_off = free_off.account.memory_peak;
    c.peak_on = free_on.account.memory_peak;
    c.depth_free = requested_depth(free_off);
    for (const ArmResult &arm : free_on.arms) {
        for (size_t i = 0; i < arm.run_coordinates.size(); ++i) {
            c.chains.push_back(arm.run_coordinates[i].total);
            c.lower += arm.run_coordinates[i].lower_bound;
            if (arm.runs[i].entered_by_switch) {
                c.switch_runs++;
                c.switch_lower += arm.run_coordinates[i].lower_bound;
            }
        }
    }
    for (int i = 0; i < 2; ++i) {
        const double fraction = i ? 0.75 : 0.5;
        const uint64_t budget = static_cast<uint64_t>(fraction * c.peak_off);
        const SeedResult off = walk(false, budget, &c.f_off[i]);
        const SeedResult on = walk(true, budget, &c.f_on[i]);
        c.d_off[i] = c.f_off[i] ? 0 : requested_depth(off);
        c.d_on[i] = c.f_on[i] ? 0 : requested_depth(on);
        c.s_off[i] = !c.f_off[i] && off.resource_stop.has_value();
        c.s_on[i] = !c.f_on[i] && on.resource_stop.has_value();
    }
    auto completes = [&](bool coordinates, uint64_t budget) {
        bool f = false;
        const SeedResult r = walk(coordinates, budget, &f);
        return !f && !r.resource_stop.has_value();
    };
    auto need = [&](bool coordinates, uint64_t peak) {
        uint64_t b = completion_budget(peak);
        while (!completes(coordinates, b)) {
            b += kMiB;
        }
        while (b > kMiB && completes(coordinates, b - kMiB)) {
            b -= kMiB;
        }
        return b;
    };
    c.need_off = need(false, c.peak_off);
    c.need_on = need(true, c.peak_on);
    const SeedResult at = walk(true, c.need_off, &failed);
    c.at_need_off = failed ? 2 : at.resource_stop.has_value() ? 1 : 0;
    c.at_need_off_depth = failed ? 0 : requested_depth(at);
    return c;
}

// the rows of the gate's tables: per detail and budget fraction (and strategy, or regime)
struct GateRow {
    std::vector<double> ratio;               // depth with / without, where without > 0
    size_t cells = 0, shallower = 0, deeper = 0, failed_on = 0, failed_off = 0;
    size_t stopped_off = 0, stopped_on = 0;
};

// per detail (and strategy, or regime): what completing takes
struct NeedRow {
    std::vector<double> peak, need;          // with / without: account peaks, completion budgets
    size_t cells = 0, complete = 0, stopped = 0, failed = 0;   // with coordinates at need_off
    size_t runs = 0, runs16 = 0;
    uint64_t max_chains = 0;
};

void add_gate(const GateCell &c, const std::string &detail, const std::string &suffix,
              std::map<std::string, GateRow> *rows, std::map<std::string, NeedRow> *needs) {
    for (int i = 0; i < 2; ++i) {
        GateRow &row = (*rows)[detail + (i ? " 75%" : " 50%") + suffix];
        row.cells++;
        row.failed_off += c.f_off[i];
        row.failed_on += c.f_on[i];
        row.stopped_off += c.s_off[i];
        row.stopped_on += c.s_on[i];
        row.shallower += c.d_on[i] < c.d_off[i];
        row.deeper += c.d_on[i] > c.d_off[i];
        if (c.d_off[i])
            row.ratio.push_back(double(c.d_on[i]) / double(c.d_off[i]));
    }
    NeedRow &n = (*needs)[detail + suffix];
    n.cells++;
    n.peak.push_back(double(c.peak_on) / double(c.peak_off));
    n.need.push_back(double(c.need_on) / double(c.need_off));
    n.complete += c.at_need_off == 0;
    n.stopped += c.at_need_off == 1;
    n.failed += c.at_need_off == 2;
    for (uint64_t chains : c.chains) {
        n.runs++;
        n.runs16 += chains >= 16;
        n.max_chains = std::max(n.max_chains, chains);
    }
}

void print_gate(const std::map<std::string, GateRow> &rows,
                const std::map<std::string, NeedRow> &needs) {
    std::cerr << "depth at the stop: detail fraction [group] | cells | stopped off/on | failed at "
                 "depth 0 off/on | shallower | median depth with/without | mean | worst" << std::endl;
    for (const auto &[key, row] : rows) {
        double mean = 0, worst = 1;
        for (double x : row.ratio) {
            mean += x / row.ratio.size();
            worst = std::min(worst, x);
        }
        std::cerr << "  " << key << " | " << row.cells << " | " << row.stopped_off << "/"
                  << row.stopped_on << " | " << row.failed_off << "/" << row.failed_on << " | "
                  << row.shallower << " | " << median(row.ratio) << " | " << mean << " | "
                  << worst << std::endl;
    }
    std::cerr << "completing: detail [group] | cells | runs (>= 16 chains) max chains | account "
                 "peak with/without median, max | completion budget (MiB) with/without median, "
                 "max | with coordinates at the budget completing without: complete / stopped / "
                 "failed at depth 0" << std::endl;
    for (const auto &[key, n] : needs) {
        std::cerr << "  " << key << " | " << n.cells << " | " << n.runs << " (" << n.runs16
                  << ") " << n.max_chains << " | " << median(n.peak) << ", "
                  << (n.peak.empty() ? 0 : *std::max_element(n.peak.begin(), n.peak.end()))
                  << " | " << median(n.need) << ", "
                  << (n.need.empty() ? 0 : *std::max_element(n.need.begin(), n.need.end()))
                  << " | " << n.complete << " / " << n.stopped << " / " << n.failed << std::endl;
    }
}

} // namespace

// M1, the depth gate on mini_refseq: every trace cell — 7 seeds (blaNDM
// both ways, three 200-bp windows of its carriers, two repeat windows) x 8 strategies (branch
// limit 0 and 2, constant switch costs 0.5 and 1 at limits 0 and 2 within a loss budget of 2,
// column labels at limit 2, with and without a switch), out to 3,000 bp, in detail full and
// graphlet — through gate_cell: the depth at the stop under 50% and 75% of its own opt-out peak
// with against without coordinates, and what completing takes with and without; unbudgeted, how
// often a run is a lower bound. mini_refseq's runs carry at most 12 chains
// (taxid columns of a few genomes each), so the column cells are also walked under the
// refseq33m projection (gate_cell's |project|: every list at the cap of 16, as refseq33m's
// taxid columns give). Budgets are bytes, exact, through the C++ API (the request's knob is
// whole MiB). The tables are in DESIGN §26.5
void coordinates_depth_at_the_stop(benchmark::State &state, const Mini &mini) {
    const std::string repeat = "CAAAGTTAGCGATGAGGCAGCCTTTTGTCTTATTCAAAGGCCTTACATTTCAAAAACTCTGCTTACC"
                               "AGGCGCATTTCGCCCAGGGGATCACCATAATAAAATGCTGAGGCCTGGCCTTTGCGTAGTGCACGCAT"
                               "CACCTCAATACCTTT";
    std::string ndm1_rc = mini.query;
    reverse_complement(ndm1_rc.begin(), ndm1_rc.end());
    const std::vector<std::pair<std::string, std::string>> seeds {
        { "ndm1", mini.query }, { "ndm1_rc", ndm1_rc },
        { "win200_00", "ACCCGACCAAGGTCACCCGCACCGCGCTGCAGAACGCCGCGTCGATCGCGGGCCTGATGATCACCACCGAAGC"
                       "GATGGTGGCCGAGGCCCCGAAGAAGGACGAGCCGGCGATGCCGGCCGGCGGCGGCATGGGCGGCATGGGCGG"
                       "CATGGATTTCTAAGCCCCGCGATCCATCAAGCAAGACCACAAAGCCCGGCCTCGT" },
        { "win200_03", "GCGATCCTTCCAACTCGTCGCAAAGCCCAGCTTCGCATAAAACGCCTCTGTCACATCGAAATCGCGCGATGG"
                       "CAGATTGGGGGTGACGTGGTCAGCCATGGCTCAGCGCAGCTTGTCGGCCATGCGGGCCGTATGAGTGATTGC"
                       "GGCGCGGCTATCGGGGGCGGAATGGCTCATCACGATCATGCTGGCCTTGGGGAACG" },
        { "win200_05", "CGCCCCGTGCGGTTACGTCGAATGTCGCGGGCGCTTTGACATCGCGCGCAGCTGGCCAGATCGCCATGGTCG"
                       "GTTTGTTCGTCGATGCGGATGATGCTGTCATCGCCGACGCACTGGTGGCAGCCAAGCTGAACGCGCTGCAGC"
                       "TGCACGGTTCGGAATCGCCCGAACGCGTGGCCCAGTTGCGCGCGCGGTTTGGCAAG" },
        { "rep0", repeat },
        { "rep3", "AGCGGTAAATCGTGGAGTGATCGACATTCACTCCGCGTTCAGCCAGCATCTCCTGCAGCTCACGGTAACTGATG"
                  "CCGTATTTGCAGTACCAGCGTACGGCCCACAGAATGATGTCACGCTGAAAATGCCGGCCTTTGAATGGGTTCATGT" },
    };
    struct Kind {
        std::string name;
        bool column;
        double cost;       // 0: forbid
        size_t limit;
    };
    const std::vector<Kind> kinds {
        { "b0", false, 0, 0 }, { "b2", false, 0, 2 },
        { "sw0.5b0", false, 0.5, 0 }, { "sw0.5b2", false, 0.5, 2 },
        { "sw1b0", false, 1, 0 }, { "sw1b2", false, 1, 2 },
        { "col_b2", true, 0, 2 }, { "col_sw0.5b2", true, 0.5, 2 },
    };
    std::map<std::string, GateRow> rows;     // by detail and fraction, and by strategy
    std::map<std::string, NeedRow> needs;
    std::map<std::string, GateRow> projected;
    std::map<std::string, NeedRow> projected_needs;
    size_t runs = 0, lower = 0, switch_runs = 0, switch_lower = 0, strictly_cells = 0;
    std::map<std::string, std::pair<size_t, size_t>> lower_by_kind;
    for (const char *detail : { "full", "graphlet" }) {
        for (const auto &[seed_name, sequence] : seeds) {
            for (const Kind &kind : kinds) {
                Strategy base;
                base.support = Support::TRACE;
                base.merge_reconverge = false;
                base.max_extension_bp = 3000;
                base.max_label_branches = kind.limit;
                if (kind.column)
                    base.seed_label_kind = LabelKind::COLUMN;
                LabelChangeCost cost = LabelChangeCost::forbid();
                if (kind.cost > 0) {
                    cost = LabelChangeCost::constant(kind.cost);
                    base.loss_budget = 2;
                }
                Seed seed;
                seed.sequence = sequence;
                const GateCell c = gate_cell(*mini.anno, seed, base, cost, detail, false);
                MEASURE_CHECK(state, c.ok);
                if (std::string(detail) == "full") {
                    // lower bounds, unbudgeted
                    runs += c.chains.size();
                    lower += c.lower;
                    switch_runs += c.switch_runs;
                    switch_lower += c.switch_lower;
                    strictly_cells += c.lower > 0;
                    lower_by_kind[kind.name].first += c.chains.size();
                    lower_by_kind[kind.name].second += c.lower;
                }
                add_gate(c, detail, "", &rows, &needs);
                add_gate(c, detail, " " + kind.name, &rows, &needs);
                for (int i = 0; i < 2; ++i) {
                    MEASURE_CHECK(state, c.d_on[i] <= c.d_off[i]);
                }
                MEASURE_CHECK(state, c.need_on >= c.need_off);
                if (kind.column) {
                    const GateCell p = gate_cell(*mini.anno, seed, base, cost, detail, true);
                    MEASURE_CHECK(state, p.ok);
                    add_gate(p, detail, " column, refseq33m projection", &projected,
                             &projected_needs);
                    add_gate(p, detail, " " + kind.name + ", refseq33m projection", &projected,
                             &projected_needs);
                    for (int i = 0; i < 2; ++i) {
                        MEASURE_CHECK(state, p.d_on[i] <= c.d_on[i]);
                    }
                }
            }
        }
    }
    std::cerr << "M1, mini_refseq:" << std::endl;
    print_gate(rows, needs);
    std::cerr << "M1, the column cells under the refseq33m projection (every list at 16):" << std::endl;
    print_gate(projected, projected_needs);
    std::cerr << lower << " of " << runs << " runs lower bounds (" << switch_lower
              << " of " << switch_runs << " switch-entered), in " << strictly_cells
              << " of 56 cells" << std::endl;
    for (const auto &[kind, n] : lower_by_kind) {
        std::cerr << "  " << kind << ": " << n.second << " of " << n.first << std::endl;
    }
}

// M1 by regime (column fixtures with at least 16 chains a run, and header-heavy seeds): the
// gate's cells on the fixtures of make_column_coord_fixtures.sh —
// coord_lockstep (200 columns of one shared sequence, 16 records each: one path, 200 runs of 16
// chains an arm) and coord_divcol (100 columns of their own sequences around a shared 600-bp
// core, 16 records each: a split into a branch a column at each end of the core; seeds in the
// core and in a flank) — and on make_wide_coord_fixture.sh's wide_coord (5,000 records of one
// sequence: column labels one run of 5,000 chains an arm, header labels 5,000 runs of one chain
// each). Column labels, branch limits 0 and 2 (header labels at 0), details full and graphlet,
// to 3,000 bp. Run from the build directory after building
// the fixtures there; the tables are in DESIGN §26.5
void coordinates_depth_by_regime(benchmark::State &state) {
    struct Fixture {
        std::string dir, annotation, fasta;
        std::vector<std::pair<std::string, size_t>> seeds;   // name, offset of 100 bp
        std::vector<LabelKind> kinds;
    };
    const std::vector<Fixture> fixtures {
        { "coord_lockstep", "cols", "fa/c000.fa", { { "1400", 1400 } }, { LabelKind::COLUMN } },
        { "coord_divcol", "cols", "fa/c000.fa", { { "core 1400", 1400 }, { "flank 200", 200 } },
          { LabelKind::COLUMN } },
        { "wide_coord", "wide", "wide.fa", { { "1400", 1400 } },
          { LabelKind::COLUMN, LabelKind::HEADER } },
    };
    std::map<std::string, GateRow> rows;
    std::map<std::string, NeedRow> needs;
    size_t built = 0;
    for (const Fixture &fx : fixtures) {
        if (!std::filesystem::exists(fx.dir + "/graph_k31.dbg")) {
            std::cerr << "skipped " << fx.dir << ": build it with ../scripts/traversal/make_"
                      << (fx.dir == "wide_coord" ? "wide_coord_fixture.sh" : "column_coord_fixtures.sh")
                      << std::endl;
            continue;
        }
        built++;
        auto graph = std::make_shared<DBGSuccinct>(2);
        MEASURE_CHECK(state, graph->load(fx.dir + "/graph_k31.dbg"));
        auto annotation = std::make_unique<annot::RowDiffBRWTCoordAnnotator>();
        MEASURE_CHECK(state, annotation->load(fx.dir + "/" + fx.annotation
                                              + annot::RowDiffBRWTCoordAnnotator::kExtension));
        using Matrix = annot::RowDiffBRWTCoordAnnotator::binary_matrix_type;
        const_cast<Matrix&>(annotation->get_matrix()).set_graph(graph.get());
        std::unique_ptr<annot::CoordToHeader> cth;
        if (std::filesystem::exists(fx.dir + "/" + fx.annotation + ".seqs")) {
            cth = std::make_unique<annot::CoordToHeader>();
            MEASURE_CHECK(state, cth->load(fx.dir + "/" + fx.annotation));
        }
        AnnotatedDBG anno(std::move(graph), std::move(annotation), false, std::move(cth));
        std::string sequence;
        {
            std::ifstream in(fx.dir + "/" + fx.fasta);
            std::string line;
            std::getline(in, line);
            std::getline(in, sequence);
        }
        MEASURE_CHECK(state, sequence.size() >= 1500u);
        for (const auto &[seed_name, offset] : fx.seeds) {
            Seed seed;
            seed.sequence = sequence.substr(offset, 100);
            for (LabelKind kind : fx.kinds) {
                const bool column = kind == LabelKind::COLUMN;
                for (size_t limit : { size_t(0), size_t(2) }) {
                    if (!column && limit)
                        continue;
                    Strategy base;
                    base.support = Support::TRACE;
                    base.merge_reconverge = false;
                    base.max_extension_bp = 3000;
                    base.max_seed_labels = 5000;
                    base.max_label_branches = limit;
                    base.seed_label_kind = kind;
                    for (const char *detail : { "full", "graphlet" }) {
                        const GateCell c = gate_cell(anno, seed, base, LabelChangeCost::forbid(),
                                                     detail, false);
                        const std::string name = fx.dir + " " + seed_name + " "
                                               + (column ? "column" : "header") + " b"
                                               + std::to_string(limit) + " " + detail;
                        MEASURE_CHECK(state, c.ok);
                        uint64_t max_chains = 0;
                        for (uint64_t x : c.chains) {
                            max_chains = std::max(max_chains, x);
                        }
                        std::cerr << name << ": runs " << c.chains.size() << ", max chains "
                                  << max_chains << " | depth free " << c.depth_free
                                  << " | 50%: " << c.d_off[0] << (c.f_off[0] ? "F" : c.s_off[0] ? "s" : "")
                                  << " / " << c.d_on[0] << (c.f_on[0] ? "F" : c.s_on[0] ? "s" : "")
                                  << " | 75%: " << c.d_off[1] << (c.f_off[1] ? "F" : c.s_off[1] ? "s" : "")
                                  << " / " << c.d_on[1] << (c.f_on[1] ? "F" : c.s_on[1] ? "s" : "")
                                  << " | peak " << c.peak_on << " / " << c.peak_off
                                  << " | completes at " << (c.need_on >> 20) << " / "
                                  << (c.need_off >> 20) << " MiB | with coordinates at "
                                  << (c.need_off >> 20) << " MiB: "
                                  << (c.at_need_off == 0 ? "complete" : c.at_need_off == 1 ? "stopped" : "failed at depth 0")
                                  << " (" << c.at_need_off_depth << ")" << std::endl;
                        const std::string regime = column
                                ? (max_chains >= 16 ? " column, >= 16 chains a run" : " column, < 16 chains")
                                : " header, many runs";
                        add_gate(c, detail, regime, &rows, &needs);
                        add_gate(c, detail, " " + fx.dir + (column ? " column" : " header"),
                                 &rows, &needs);
                        for (int i = 0; i < 2; ++i) {
                            MEASURE_CHECK(state, c.d_on[i] <= c.d_off[i]);
                        }
                    }
                }
            }
        }
    }
    if (!built) {
        state.SkipWithError("no fixture built (scripts/traversal/make_column_coord_fixtures.sh, "
                            "make_wide_coord_fixture.sh)");
        return;
    }
    print_gate(rows, needs);
}

// M2: the wide coordinate fixture (scripts/traversal/
// make_wide_coord_fixture.sh: one 3,000-bp sequence in 5,000 records under one column, sparse
// row-diff anchors) walked under trace from a 100-bp seed in its middle: column labels (one run
// an arm with 5,000 chains) and header labels (5,000 runs an arm with one chain each), caps 16
// and "unlimited", details full and graphlet. Per cell: the text and its coordinates' share,
// the account and its coordinate share, the time to build the text, the delivery reserve's
// estimate of the text — the coordinate share at kCoordinateAccountPerTextByte must be at least
// the coordinates' real text — and the walk-until it gives an attempt of one seed
// at 30 s against the same request without coordinates. Run from
// the build directory after building the fixture there
void coordinate_reserve_on_the_wide_fixture(benchmark::State &state) {
    const std::string dir = "wide_coord";
    if (!std::filesystem::exists(dir + "/graph_k31.dbg")) {
        state.SkipWithError("build it with scripts/traversal/make_wide_coord_fixture.sh");
        return;
    }
    auto graph = std::make_shared<DBGSuccinct>(2);
    MEASURE_CHECK(state, graph->load(dir + "/graph_k31.dbg"));
    auto annotation = std::make_unique<annot::RowDiffBRWTCoordAnnotator>();
    MEASURE_CHECK(state, annotation->load(dir + "/wide" + annot::RowDiffBRWTCoordAnnotator::kExtension));
    using Matrix = annot::RowDiffBRWTCoordAnnotator::binary_matrix_type;
    const_cast<Matrix&>(annotation->get_matrix()).set_graph(graph.get());
    auto cth = std::make_unique<annot::CoordToHeader>();
    MEASURE_CHECK(state, cth->load(dir + "/wide"));
    AnnotatedDBG anno(std::move(graph), std::move(annotation), false, std::move(cth));
    std::string sequence;
    {
        std::ifstream in(dir + "/wide.fa");
        std::string line;
        std::getline(in, line);
        std::getline(in, sequence);
    }
    MEASURE_CHECK(state, sequence.size() >= 1500u);
    Seed seed;
    seed.sequence = sequence.substr(1400, 100);
    std::cerr << "M2: kind cap detail | walk s | text B (coordinates B) | account B (coordinates "
                 "B) | build s | coordinate estimate / text | whole estimate / text (configured "
                 "ratio) | pre-split estimate / text | reserve ms with / without | walk_until ms "
                 "with / without" << std::endl;
    for (bool column : { true, false }) {
        for (size_t cap : { size_t(16), Strategy::kUnlimited }) {
            for (const char *detail : { "full", "graphlet" }) {
                Strategy st;
                st.support = Support::TRACE;
                st.merge_reconverge = false;
                st.max_extension_bp = 3000;
                st.max_seed_labels = 5000;
                st.seed_label_kind = column ? LabelKind::COLUMN : LabelKind::HEADER;
                st.max_coordinate_occurrences = cap;
                const LabelChangeCost cost = LabelChangeCost::forbid();
                Written w[2];
                uint64_t account[2] = { 0, 0 }, coordinates[2] = { 0, 0 };
                double walk_s[2] = { 0, 0 }, reserve[2] = { 0, 0 }, until[2] = { 0, 0 };
                for (bool on : { false, true }) {
                    Strategy s = st;
                    s.coordinates = on;
                    price(&s, cost, detail);
                    LabelOracle oracle(anno);
                    Timer timer;
                    const SeedResult r = traverse_seed(oracle, seed, s, cost);
                    walk_s[on] = timer.elapsed();
                    w[on] = write_result(r, seed, s, detail, oracle);
                    account[on] = r.account.memory_final;
                    coordinates[on] = r.account.coordinates;
                    // the reserve an attempt of this one seed keeps for it (configured rates and
                    // ratios: what a server that has measured nothing yet uses)
                    mtg::cli::AttemptSettings settings;
                    auto attempt = std::make_shared<mtg::cli::Attempt>(
                            0, std::chrono::system_clock::time_point(), settings, "0", true);
                    mtg::cli::AttemptIds ids;
                    ids.attempt_id = "m2";
                    attempt->set_ids(ids);
                    attempt->set_delivery_detail(detail);
                    attempt->set_bound(1, 30'000);
                    attempt->progress(account[on], coordinates[on]);
                    reserve[on] = attempt->reserve_ms();
                    until[on] = attempt->walk_until_ms();
                }
                const double ratio = std::string(detail) == "graphlet"
                                   ? mtg::cli::AttemptSettings().account_per_text_byte_graphlet
                                   : mtg::cli::AttemptSettings().account_per_text_byte_json;
                const uint64_t k = mtg::cli::kCoordinateAccountPerTextByte;
                const double coordinate_estimate = std::ceil(double(coordinates[1]) / k);
                const double whole = std::ceil(double(account[1] - coordinates[1]) / ratio)
                                   + coordinate_estimate;
                // what the reserve estimated before the split, at the ratio the rest measures
                const double rest_ratio = double(account[1] - coordinates[1])
                                        / double(w[1].text - w[1].coordinate_text);
                const double before = double(account[1]) / rest_ratio;
                std::cerr << "  " << (column ? "column" : "header") << " "
                          << (cap == Strategy::kUnlimited ? std::string("unlimited") : std::to_string(cap))
                          << " " << detail << " | " << walk_s[1] << " (" << walk_s[0] << ") | "
                          << w[1].text << " (" << w[1].coordinate_text << ") | " << account[1]
                          << " (" << coordinates[1] << ") | " << w[1].seconds << " ("
                          << w[0].seconds << ") | " << coordinate_estimate / w[1].coordinate_text
                          << " | " << whole / w[1].text << " | " << before / w[1].text << " | "
                          << reserve[1] << " / " << reserve[0] << " | " << until[1] << " / "
                          << until[0] << std::endl;
                // the coordinate share's estimate bounds its real text
                MEASURE_CHECK(state, coordinate_estimate >= double(w[1].coordinate_text));
                // the rest is the opt-out walk's, exactly
                MEASURE_CHECK(state, account[0] == account[1] - coordinates[1]);
                MEASURE_CHECK(state, w[0].text == w[1].text - w[1].coordinate_text);
            }
        }
    }
}

// ---------------------------------------------------------------- the dummy fraction

// graph::pattern::sample_real_fraction on a random graph of 300 records of 100,000 bases at
// k = 31 (about 30 million k-mers, few sources: f near 1): its time, and the exact fraction
// (the entries with W != $ whose node holds no $) inside its interval
void real_fraction_on_a_large_graph(benchmark::State &state) {
    const size_t k = 31;
    std::mt19937 rng(3);
    std::vector<std::string> records(300, std::string(100'000, 'A'));
    for (std::string &record : records) {
        for (char &c : record) {
            c = "ACGT"[rng() % 4];
        }
    }
    boss::BOSSConstructor constructor(k - 1, false, 0, "", 0, 1, 0,
                                      kmer::ContainerType::VECTOR,
                                      std::filesystem::temp_directory_path().string());
    constructor.add_sequences(std::move(records));
    DBGSuccinct graph(new boss::BOSS(&constructor), DeBruijnGraph::BASIC);
    const auto start = std::chrono::steady_clock::now();
    const pattern::RealFraction f = pattern::sample_real_fraction(graph);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    const boss::BOSS &boss = graph.get_boss();
    sdsl::bit_vector source(boss.get_W().size(), false);
    boss.mark_source_dummy_edges(&source, 1);
    uint64_t entries = 0, dummies = 0;
    for (uint64_t e = 1; e < source.size(); ++e) {
        if (boss.get_W(e) % boss.alph_size) {
            ++entries;
            dummies += source[e];
        }
    }
    const double exact = entries ? double(entries - dummies) / double(entries) : 1;
    std::cerr << "random graph: " << graph.max_index() << " edges, sampled f " << f.value
              << " [" << f.lower << ", " << f.upper << "] in " << seconds * 1000
              << " ms; exact f " << exact << std::endl;
    MEASURE_CHECK(state, f.lower <= exact && exact <= f.upper);
}


// ---------------------------------------------------------------- the benchmarks

void BM_TraversalCoordRuns(benchmark::State &state) {
    for (auto _ : state) {
        coord_runs(state);
    }
}

void BM_TraversalCoordMappingOnTheMini(benchmark::State &state) {
    for (auto _ : state) {
        if (const Mini *mini = mini_index(state))
            coord_mapping(state, *mini);
    }
}

void BM_TraversalCoordinatesDepthAtTheStop(benchmark::State &state) {
    for (auto _ : state) {
        if (const Mini *mini = mini_index(state))
            coordinates_depth_at_the_stop(state, *mini);
    }
}

void BM_TraversalCoordinatesDepthByRegime(benchmark::State &state) {
    for (auto _ : state) {
        coordinates_depth_by_regime(state);
    }
}

void BM_TraversalCoordinateReserveOnTheWideFixture(benchmark::State &state) {
    for (auto _ : state) {
        coordinate_reserve_on_the_wide_fixture(state);
    }
}

void BM_PatternRealFractionOnALargeGraph(benchmark::State &state) {
    for (auto _ : state) {
        real_fraction_on_a_large_graph(state);
    }
}

BENCHMARK(BM_TraversalCoordRuns)->Iterations(1)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_TraversalCoordMappingOnTheMini)->Iterations(1)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_TraversalCoordinatesDepthAtTheStop)->Iterations(1)->Unit(benchmark::kSecond);
BENCHMARK(BM_TraversalCoordinatesDepthByRegime)->Iterations(1)->Unit(benchmark::kSecond);
BENCHMARK(BM_TraversalCoordinateReserveOnTheWideFixture)->Iterations(1)->Unit(benchmark::kSecond);
BENCHMARK(BM_PatternRealFractionOnALargeGraph)->Iterations(1)->Unit(benchmark::kSecond);

} // namespace
