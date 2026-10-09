#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <json/json.h>
#include <sdust.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"

#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/pattern.hpp"
#include "cli/pattern_retrieval_impl.hpp"
#include "cli/pattern_support.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/traversal/label_oracle.hpp"


// The supported-path search (long_search "supported_paths", SPEC §20): the label and trace
// trackers (cli::PathTracker) driving the engine's extension (PatternSearch::count with a
// SupportTracker and a PathSink), the sink's supported paths with their labels from the frames
// (cli::SupportedPathSink), the release rule (cli::supported_release), and two engine rules
// beside them (SupportTracker::prepare; the low-complexity note beside threshold stops).
// Tiny indexes built from explicit records in labelled columns with their coordinates and
// record mapping. The expectations come from oracles that never ask the engine, the trackers
// or support_step:
//  - a WALK-TREE oracle over the k-mers of the records (as deposited on BASIC graphs; with their
//    reverse complements on CANONICAL and PRIMARY ones): from every k-mer that is a prefix of an
//    instance of the oriented pattern (this file's own IUPAC table and standard genetic code),
//    every extension by a base whose next k-mer exists, the support of each prefix computed from
//    the records — the columns carrying every k-mer of it (label level), and the records holding
//    it whole (record level) — and a branch followed only while it has support: the walks
//    completed, the supported ones, the branches pruned, the branches entered, the rows read;
//  - a RECORD-SCAN oracle: where a column's records hold a path whole (1-based starts), and,
//    without the record mapping, where a column's coordinates (its records' k-mers numbered one
//    after another) hold its k-mers consecutively;
//  - the paths of long_search "paths" with output.labels "all" (the labelled retrieval of every
//    walk, its labels verified by a join of the k-mers' coordinate lists): an independent
//    implementation the supported paths must agree with, path by path.

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::pattern;
using namespace mtg::cli;
using graph::traversal::LabelOracle;
using graph::traversal::Support;
using Clock = Deadline::Clock;
typedef DeBruijnGraph::node_index node_index;


// ------------------------------------------------------------------ the test's own tables

const std::map<char, std::string> kCodes = {
    { 'A', "A" }, { 'C', "C" }, { 'G', "G" }, { 'T', "T" }, { 'R', "AG" }, { 'Y', "CT" },
    { 'S', "CG" }, { 'W', "AT" }, { 'K', "GT" }, { 'M', "AC" }, { 'B', "CGT" },
    { 'D', "AGT" }, { 'H', "ACT" }, { 'V', "ACG" }, { 'N', "ACGT" },
};

char complement_code(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'R': return 'Y';
        case 'Y': return 'R';
        case 'K': return 'M';
        case 'M': return 'K';
        case 'B': return 'V';
        case 'V': return 'B';
        case 'D': return 'H';
        case 'H': return 'D';
        default: return c;     // S, W, N
    }
}

std::string rc(const std::string &s) {
    std::string out(s.rbegin(), s.rend());
    for (char &c : out) {
        c = complement_code(c);
    }
    return out;
}

bool admits(char code, char base) {
    return kCodes.at(code).find(base) != std::string::npos;
}

// the standard genetic code (NCBI table 1), its 64 codons in TCAG order
const std::string kStandardCode
        = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

char amino_acid(const std::string &codon) {
    auto rank = [](char b) { return b == 'T' ? 0 : b == 'C' ? 1 : b == 'A' ? 2 : 3; };
    return kStandardCode[16 * rank(codon[0]) + 4 * rank(codon[1]) + rank(codon[2])];
}

bool residue_admits(char residue, char aa) {
    switch (residue) {
        case 'X': return aa != '*';
        case 'B': return aa == 'D' || aa == 'N';
        case 'Z': return aa == 'E' || aa == 'Q';
        case 'J': return aa == 'I' || aa == 'L';
        default: return aa == residue;
    }
}

std::vector<std::string> codons_of(char residue) {
    std::vector<std::string> out;
    const std::string bases = "ACGT";
    for (char a : bases) {
        for (char b : bases) {
            for (char c : bases) {
                const std::string codon { a, b, c };
                if (residue_admits(residue, amino_acid(codon)))
                    out.push_back(codon);
            }
        }
    }
    return out;
}

std::string translate(const std::string &dna) {
    std::string out;
    for (size_t i = 0; i + 3 <= dna.size(); i += 3) {
        out.push_back(amino_acid(dna.substr(i, 3)));
    }
    return out;
}

/**
 * One oriented pattern as a predicate on strings: s is a prefix of an instance. An IUPAC
 * pattern position by position; a peptide P codon by codon, read forward (FORWARD), or as the
 * reverse complement of its instances (REVERSE: rc(s) is a suffix of an instance of P).
 */
struct Oriented {
    Orientation orientation = Orientation::FORWARD;
    bool peptide = false;
    std::string q;
    size_t length = 0;

    bool prefix_ok(const std::string &s) const {
        if (s.size() > length)
            return false;
        if (!peptide) {
            for (size_t i = 0; i < s.size(); ++i) {
                if (!admits(q[i], s[i]))
                    return false;
            }
            return true;
        }
        const size_t m = s.size();
        std::string known(length, '?');
        for (size_t i = 0; i < m; ++i) {
            if (orientation == Orientation::REVERSE) {
                known[length - 1 - i] = complement_code(s[i]);
            } else {
                known[i] = s[i];
            }
        }
        for (size_t c = 0; c < q.size(); ++c) {
            bool any = false;
            for (const std::string &codon : codons_of(q[c])) {
                bool ok = true;
                for (size_t j = 0; j < 3 && ok; ++j) {
                    ok = known[3 * c + j] == '?' || known[3 * c + j] == codon[j];
                }
                if (ok) {
                    any = true;
                    break;
                }
            }
            if (!any)
                return false;
        }
        return true;
    }
};

std::vector<Oriented> orientations(const std::string &text, bool peptide, Strands strands) {
    std::vector<Oriented> out;
    const size_t length = peptide ? 3 * text.size() : text.size();
    if (!peptide && rc(text) == text)
        return { Oriented { Orientation::PALINDROMIC, false, text, length } };
    if (strands != Strands::REVERSE)
        out.push_back(Oriented { Orientation::FORWARD, peptide, text, length });
    if (strands != Strands::FORWARD) {
        out.push_back(Oriented { Orientation::REVERSE, peptide, peptide ? text : rc(text),
                                 length });
    }
    return out;
}

uint8_t orientation_rank(Orientation o) {
    return o == Orientation::FORWARD ? 0 : o == Orientation::REVERSE ? 1 : 2;
}


// ------------------------------------------------------------------ the index

struct Record {
    std::string column;
    std::string header;
    std::string seq;
};

struct Index {
    size_t k = 0;
    DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC;
    bool coordinates = false;
    std::vector<Record> records;
    std::unique_ptr<AnnotatedDBG> anno;
    std::unique_ptr<annot::CoordToHeader> cth;
};

// each record's k-mers numbered from its column's running count (as `annotate --coordinates`
// numbers them), and the record mapping built from the same records
template <class Annotation>
Index build(size_t k, const std::vector<Record> &records, bool coordinates,
            DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
    Index idx;
    idx.k = k;
    idx.mode = mode;
    idx.coordinates = coordinates;
    idx.records = records;
    std::vector<std::string> seqs, labels;
    std::vector<uint64_t> starts;
    std::map<std::string, uint64_t> next;
    for (const Record &r : records) {
        EXPECT_GE(r.seq.size(), k);
        seqs.push_back(r.seq);
        labels.push_back(r.column);
        starts.push_back(next[r.column]);
        next[r.column] += r.seq.size() - k + 1;
    }
    idx.anno = test::build_anno_graph<DBGSuccinct, Annotation>(
            k, seqs, labels, mode, coordinates, coordinates ? starts : std::vector<uint64_t>{});
    if (coordinates) {
        const auto &encoder = idx.anno->get_annotator().get_label_encoder();
        std::vector<std::vector<std::string>> headers(encoder.size());
        std::vector<std::vector<uint64_t>> num_kmers(encoder.size());
        for (const Record &r : records) {
            const size_t c = encoder.encode(r.column);
            headers[c].push_back(r.header);
            num_kmers[c].push_back(r.seq.size() - k + 1);
        }
        idx.cth = std::make_unique<annot::CoordToHeader>(std::move(headers),
                                                         std::move(num_kmers));
    }
    return idx;
}

Index build_rd(size_t k, const std::vector<Record> &records,
               DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
    return build<annot::RowDiffColumnAnnotator>(k, records, mode == DeBruijnGraph::BASIC, mode);
}


// ------------------------------------------------------------------ the oracles

// per label of a path: its support, whether its occurrences are listed, and which ("seq:start"
// with the record mapping, "g:coordinate" without)
using LabelView = std::tuple<std::string, bool, std::set<std::string>>;
using PathLabels = std::map<std::string, LabelView>;

struct ExpectedPath {
    Orientation orientation = Orientation::FORWARD;
    std::string sequence;
    node_index anchor = 0;
    PathLabels labels;
};

struct Expected {
    // the supported walks, in answer order (anchor node, orientation, sequence)
    std::vector<ExpectedPath> supported;
    std::map<Orientation, uint64_t> supported_by;
    uint64_t walks = 0;
    uint64_t branches_pruned = 0;
    uint64_t candidates = 0;
    bool pruned_early = false;
    uint64_t anchors = 0;
    // the rows (k-mers; the canonical one of a pair on CANONICAL and PRIMARY graphs) of the
    // anchors and of the branches entered after them
    std::set<std::string> anchor_rows;
    std::set<std::string> pushed_rows;
};

struct Oracle {
    const Index &idx;
    const size_t k;
    // a k-mer and its reverse complement are one (CANONICAL, PRIMARY)
    const bool both;
    std::set<std::string> kmers;
    std::map<std::string, std::set<std::string>> carried;
    std::map<std::string, std::vector<const Record*>> columns;
    // per column, the k-mer at each column coordinate
    std::map<std::string, std::vector<std::string>> coordinates;

    explicit Oracle(const Index &idx)
          : idx(idx), k(idx.k), both(idx.mode != DeBruijnGraph::BASIC) {
        for (const Record &r : idx.records) {
            columns[r.column].push_back(&r);
            for (size_t i = 0; i + k <= r.seq.size(); ++i) {
                const std::string x = r.seq.substr(i, k);
                kmers.insert(x);
                carried[r.column].insert(x);
                coordinates[r.column].push_back(x);
                if (both) {
                    kmers.insert(rc(x));
                    carried[r.column].insert(rc(x));
                }
            }
        }
    }

    std::string row(const std::string &x) const {
        return both ? std::min(x, rc(x)) : x;
    }

    std::set<std::string> carriers(const std::string &s) const {
        std::set<std::string> out;
        for (const auto &[column, set] : carried) {
            bool all = true;
            for (size_t i = 0; i + k <= s.size() && all; ++i) {
                all = set.count(s.substr(i, k));
            }
            if (all)
                out.insert(column);
        }
        return out;
    }

    // (seq_id, 1-based start) of |s| in the records of |column|, as deposited
    std::set<std::string> holders(const std::string &column, const std::string &s) const {
        std::set<std::string> out;
        const auto &recs = columns.at(column);
        for (uint64_t id = 0; id < recs.size(); ++id) {
            const std::string &seq = recs[id]->seq;
            for (size_t i = 0; i + s.size() <= seq.size(); ++i) {
                if (seq.compare(i, s.size(), s) == 0)
                    out.insert(std::to_string(id) + ":" + std::to_string(i + 1));
            }
        }
        return out;
    }

    // the column coordinates c with c + i holding the i-th k-mer of |s| (records concatenated)
    std::set<std::string> chains(const std::string &column, const std::string &s) const {
        std::set<std::string> out;
        const auto &at = coordinates.at(column);
        const size_t n = s.size() - k + 1;
        for (size_t c = 0; c + n <= at.size(); ++c) {
            bool all = true;
            for (size_t i = 0; i < n && all; ++i) {
                all = at[c + i] == s.substr(i, k);
            }
            if (all)
                out.insert("g:" + std::to_string(c));
        }
        return out;
    }

    bool supported(const std::string &s, bool records) const {
        for (const std::string &c : carriers(s)) {
            if (!records || holders(c, s).size())
                return true;
        }
        return false;
    }

    /**
     * The supported-path search of |text| as the tracker should run it: |records| the record
     * level (else the label level), |occurrences| the labels' occurrences listed (record level:
     * their records; |global|: the chains of the column coordinates).
     */
    Expected expect(const std::string &text, bool peptide, Strands strands, bool records,
                    bool occurrences, bool global) const {
        Expected e;
        const DeBruijnGraph &graph = idx.anno->get_graph();
        struct Found {
            node_index anchor;
            uint8_t rank;
            ExpectedPath path;
        };
        std::vector<Found> found;
        for (const Oriented &o : orientations(text, peptide, strands)) {
            const size_t L = o.length;
            std::function<void(const std::string&, node_index)> dfs
                    = [&](const std::string &s, node_index anchor) {
                for (char b : std::string("ACGT")) {
                    const std::string t = s + b;
                    if (!o.prefix_ok(t) || !kmers.count(t.substr(t.size() - k)))
                        continue;
                    ++e.candidates;
                    e.pushed_rows.insert(row(t.substr(t.size() - k)));
                    const bool alive = supported(t, records);
                    if (t.size() == L) {
                        ++e.walks;
                        if (!alive) {
                            ++e.branches_pruned;
                            continue;
                        }
                        ExpectedPath p;
                        p.orientation = o.orientation;
                        p.sequence = t;
                        p.anchor = anchor;
                        for (const std::string &c : carriers(t)) {
                            std::set<std::string> held = holders(c, t);
                            const bool verified = records && held.size();
                            if (!occurrences) {
                                p.labels[c] = LabelView(verified ? "record_verified"
                                                                 : "label_intersection",
                                                        false, {});
                            } else if (records) {
                                p.labels[c] = LabelView(verified ? "record_verified"
                                                                 : "label_intersection",
                                                        true, held);
                            } else if (global) {
                                p.labels[c] = LabelView("label_intersection", true,
                                                        chains(c, t));
                            } else {
                                p.labels[c] = LabelView("label_intersection", false, {});
                            }
                        }
                        e.supported_by[o.orientation]++;
                        found.push_back(Found { anchor, orientation_rank(o.orientation), p });
                    } else if (!alive) {
                        ++e.branches_pruned;
                        e.pruned_early = true;
                    } else {
                        dfs(t, anchor);
                    }
                }
            };
            for (const std::string &x : kmers) {
                if (L <= k || !o.prefix_ok(x))
                    continue;
                ++e.anchors;
                e.anchor_rows.insert(row(x));
                node_index anchor = 0;
                graph.map_to_nodes_sequentially(x, [&](node_index v) { anchor = v; });
                if (!supported(x, records)) {
                    ++e.branches_pruned;
                    e.pruned_early = true;
                    continue;
                }
                dfs(x, anchor);
            }
        }
        std::stable_sort(found.begin(), found.end(), [](const Found &a, const Found &b) {
            return std::tie(a.anchor, a.rank, a.path.sequence)
                    < std::tie(b.anchor, b.rank, b.path.sequence);
        });
        for (Found &f : found) {
            e.supported.push_back(std::move(f.path));
        }
        return e;
    }
};

PathLabels labels_of(const Json::Value &labels) {
    PathLabels out;
    for (const Json::Value &l : labels) {
        std::set<std::string> occ;
        const bool listed = l.isMember("occurrence_list") && l["occurrence_list"].isArray();
        if (listed) {
            for (const Json::Value &o : l["occurrence_list"]) {
                if (o.isMember("seq_id")) {
                    const std::string nt = o["nt_coords"].asString();
                    occ.insert(std::to_string(o["seq_id"].asUInt64()) + ":"
                               + nt.substr(0, nt.find('-')));
                } else {
                    EXPECT_EQ(0u, o["offset"].asUInt64());
                    occ.insert("g:" + std::to_string(o["kmer_coord"].asUInt64()));
                }
            }
        }
        out[l["column"].asString()] = LabelView(l["support"].asString(), listed, occ);
    }
    return out;
}


// ------------------------------------------------------------------ running the search

struct Ask {
    Mode mode = Mode::PARTIAL;
    Support level = Support::TRACE;
    bool labels = true;
    bool require_verified = false;
    bool occurrences = true;
    // the record mapping given to the oracle (the index's .seqs)
    bool records = true;
    Strands strands = Strands::BOTH;
    uint64_t max_paths = 1'000'000;
    uint64_t max_anchors = 1'000'000;
    bool stop_at_threshold = false;
    uint64_t max_memory = uint64_t(1) << 30;
    uint64_t max_work = std::numeric_limits<uint64_t>::max() / 4;
    uint64_t max_steps = 1'000'000'000;
    uint64_t max_labels = 1'000'000;
    uint64_t max_occurrences = 1'000'000;
    // a virtual clock and the time budget it is read against (none: no deadline)
    std::function<Clock::time_point()> clock;
    double time_budget_ms = std::numeric_limits<double>::infinity();
    std::function<bool(uint64_t)> deny;
};

struct Outcome {
    Result result;
    std::vector<Context> paths;
    std::optional<Extraction> x;
    uint64_t named = 0;
    std::optional<LabelsAnswer> labels;
    PathSupportWork work;
    std::string stop;
    Json::Value rows_refused;
    uint64_t units = 0;
    uint64_t peak = 0;
    uint64_t held_after = 0;
    uint64_t accepted = 0;
    bool output_cut = false;
    // the request's annotation units before each annotation read began
    std::vector<uint64_t> units_before_reads;
};

Pattern pattern_of(const std::string &text, bool peptide) {
    return Pattern::parse(peptide ? PatternKind::PROTEIN : PatternKind::IUPAC, text);
}

Outcome run_support(const Index &idx, const std::string &text, bool peptide, const Ask &ask) {
    Outcome out;
    LabelOracle oracle(*idx.anno, ask.records ? idx.cth.get() : nullptr);
    retrieval::Account account(ask.max_memory);
    const Clock::time_point start = Clock::now();
    Budget budget(ask.max_steps,
                  Deadline(start, ask.time_budget_ms, 0,
                           ask.clock ? ask.clock
                                     : std::function<Clock::time_point()>(&Clock::now)));
    uint64_t units = 0;
    RetrievalLimits limits;
    limits.max_annotation_work = ask.max_work;
    limits.max_memory_bytes = ask.max_memory;
    limits.occurrences = ask.occurrences;
    limits.max_labels = ask.max_labels;
    limits.max_occurrences_per_label = ask.max_occurrences;
    const DeBruijnGraph &graph = idx.anno->get_graph();
    oracle.test_read_hook = [&](size_t) { out.units_before_reads.push_back(units); };
    PathSupportEnv env { oracle, budget, account, units, limits,
                         PatternSearch::support(graph).mode, oracle.decode_charged(), nullptr,
                         ask.deny };
    PathSupportOptions options;
    options.level = ask.level;
    options.labels = ask.labels;
    options.require_verified = ask.require_verified;
    PathTracker tracker(env, options);
    SupportedPathSink sink(tracker);
    PatternSearch search(graph);
    Request rq;
    rq.mode = ask.mode;
    rq.strands = ask.strands;
    rq.extend_paths = true;
    rq.max_paths = ask.max_paths;
    rq.max_anchors = ask.max_anchors;
    rq.stop_at_threshold = ask.stop_at_threshold;
    rq.min_information_bits = 0;
    rq.support = &tracker;
    rq.sink = &sink;
    const Pattern p = pattern_of(text, peptide);
    tracker.begin_pattern(p.length());
    sink.begin_pattern(ask.mode, ask.max_paths);
    out.result = search.count(p, rq, budget);
    out.x = supported_release(out.result, rq, sink, &out.named);
    if (ask.labels && out.x)
        out.labels = sink.labels_answer(*out.x, out.named, p.length(), Json::Value("g"));
    out.paths = sink.paths();
    out.work = tracker.work();
    out.accepted = sink.accepted();
    out.output_cut = sink.output_cut();
    if (tracker.stop_reason()) {
        out.stop = tracker.stop_reason();
    } else if (sink.stop_reason()) {
        out.stop = sink.stop_reason();
    }
    out.rows_refused = tracker.rows_refused();
    sink.end_pattern(out.x ? out.x->returned : 0);
    tracker.end_pattern();
    out.units = units;
    out.peak = account.peak();
    out.held_after = account.held();
    return out;
}

// the supported paths of |out| against the oracle's: the counts, the paths in answer order (on
// BASIC graphs; as a set elsewhere) and each path's labels; |out| completed without a stop
void check_complete(const Index &idx, const Expected &e, const Outcome &out, const Ask &ask) {
    ASSERT_FALSE(out.result.refusal);
    ASSERT_TRUE(out.result.anchors);
    const AnchorCounts &a = *out.result.anchors;
    ASSERT_FALSE(out.result.stop) << out.stop;
    EXPECT_EQ(e.anchors ? Extension::COMPLETED : Extension::NO_ANCHORS, a.extension);
    EXPECT_EQ(Relation::EXACT, a.supported.relation);
    EXPECT_EQ(e.supported.size(), a.supported.value);
    for (const auto &[o, c] : a.supported_by_orientation) {
        EXPECT_EQ(Relation::EXACT, c.relation);
        EXPECT_EQ(e.supported_by.count(o) ? e.supported_by.at(o) : 0, c.value);
    }
    EXPECT_EQ(e.walks, a.paths.value);
    EXPECT_EQ(e.pruned_early ? Relation::AT_LEAST : Relation::EXACT, a.paths.relation);
    EXPECT_EQ(e.branches_pruned, a.branches_pruned);
    EXPECT_EQ(e.pruned_early, a.pruned_before_completion);
    EXPECT_EQ(e.candidates, a.candidates_examined);
    EXPECT_EQ(e.supported.size(), out.accepted);

    if (ask.mode == Mode::COUNT) {
        EXPECT_TRUE(out.paths.empty());
        EXPECT_FALSE(out.x);
        return;
    }
    ASSERT_EQ(e.supported.size(), out.paths.size());
    ASSERT_TRUE(out.x);
    EXPECT_TRUE(out.x->complete);
    std::vector<std::pair<uint8_t, std::string>> got, want;
    for (size_t i = 0; i < out.paths.size(); ++i) {
        got.emplace_back(orientation_rank(out.paths[i].orientation), out.paths[i].sequence);
        want.emplace_back(orientation_rank(e.supported[i].orientation), e.supported[i].sequence);
    }
    if (idx.mode == DeBruijnGraph::BASIC) {
        EXPECT_EQ(want, got);
    } else {
        using Bag = std::multiset<std::pair<uint8_t, std::string>>;
        EXPECT_EQ(Bag(want.begin(), want.end()), Bag(got.begin(), got.end()));
    }
    for (const Context &c : out.paths) {
        EXPECT_EQ(c.sequence.substr(0, idx.k), idx.anno->get_graph().get_node_sequence(c.node));
        ASSERT_EQ(c.sequence.size() - idx.k + 1, c.path.size());
        for (size_t j = 0; j < c.path.size(); ++j) {
            EXPECT_EQ(c.sequence.substr(j, idx.k),
                      idx.anno->get_graph().get_node_sequence(c.path[j]));
        }
    }
    if (!ask.labels)
        return;
    ASSERT_TRUE(out.labels);
    const LabelsAnswer &l = *out.labels;
    EXPECT_FALSE(l.withheld);
    EXPECT_FALSE(l.stop);
    EXPECT_TRUE(l.complete);
    ASSERT_EQ(out.paths.size(), l.result_fields.size());
    std::map<std::pair<uint8_t, std::string>, PathLabels> expected_labels;
    for (const ExpectedPath &p : e.supported) {
        PathLabels pl = p.labels;
        if (ask.require_verified) {
            for (auto it = pl.begin(); it != pl.end(); ) {
                it = std::get<0>(it->second) == "record_verified" ? std::next(it) : pl.erase(it);
            }
        }
        expected_labels[{ orientation_rank(p.orientation), p.sequence }] = pl;
    }
    for (size_t i = 0; i < out.paths.size(); ++i) {
        const Json::Value &f = l.result_fields[i];
        EXPECT_EQ("complete", f["labels_status"].asString());
        const auto key = std::make_pair(orientation_rank(out.paths[i].orientation),
                                        out.paths[i].sequence);
        EXPECT_EQ(expected_labels[key], labels_of(f["labels"])) << out.paths[i].sequence;
    }
}


// ------------------------------------------------------------------ hand-made cases

TEST(PatternSupport, TwoRecordsOneLabel) {
    // k = 3, A holds ACG and CGT in two records: the walk ACGT carries A on both k-mers, and no
    // record holds it whole
    const Index idx = build_rd(3, { { "A", "a0", "ACG" }, { "A", "a1", "CGT" } });
    const Oracle oracle(idx);
    Ask label;
    label.level = Support::KMER;
    Outcome l = run_support(idx, "ACGT", false, label);
    check_complete(idx, oracle.expect("ACGT", false, Strands::BOTH, false, false, false), l,
                   label);
    EXPECT_EQ(1u, l.result.anchors->supported.value);
    ASSERT_EQ(1u, l.paths.size());
    EXPECT_EQ("label_intersection", l.labels->result_fields[0]["support"].asString());

    Ask record;
    Outcome r = run_support(idx, "ACGT", false, record);
    check_complete(idx, oracle.expect("ACGT", false, Strands::BOTH, true, true, false), r,
                   record);
    EXPECT_EQ(Relation::EXACT, r.result.anchors->paths.relation);
    EXPECT_EQ(1u, r.result.anchors->paths.value);
    EXPECT_EQ(0u, r.result.anchors->supported.value);
    EXPECT_EQ(1u, r.result.anchors->branches_pruned);
    EXPECT_FALSE(r.result.anchors->pruned_before_completion);
    EXPECT_TRUE(r.paths.empty());
    EXPECT_TRUE(r.x->complete);
}

TEST(PatternSupport, OneRowAOtherRowB) {
    // ACG only in A, CGT only in B: the walk ACGT exists, its intersection is empty
    const Index idx = build_rd(3, { { "A", "a0", "ACGA" }, { "B", "b0", "TCGT" } });
    for (Support level : { Support::KMER, Support::TRACE }) {
        Ask ask;
        ask.level = level;
        Outcome out = run_support(idx, "ACGT", false, ask);
        check_complete(idx, Oracle(idx).expect("ACGT", false, Strands::BOTH,
                                               level == Support::TRACE, true, false), out, ask);
        EXPECT_EQ(1u, out.result.anchors->paths.value);
        EXPECT_EQ(0u, out.result.anchors->supported.value);
        EXPECT_TRUE(out.paths.empty());
    }
}

TEST(PatternSupport, RecombinantBubble) {
    // the mosaics of a bubble die at their last k-mer: complete walks, counted
    const Index idx = build_rd(4, { { "A", "a0", "ACGTTGCA" }, { "B", "b0", "TCGTTGCC" } });
    for (Support level : { Support::KMER, Support::TRACE }) {
        Ask ask;
        ask.level = level;
        Outcome out = run_support(idx, "WCGTTGCM", false, ask);
        check_complete(idx, Oracle(idx).expect("WCGTTGCM", false, Strands::BOTH,
                                               level == Support::TRACE, true, false), out, ask);
        EXPECT_EQ(Relation::EXACT, out.result.anchors->paths.relation);
        EXPECT_EQ(4u, out.result.anchors->paths.value);
        EXPECT_EQ(2u, out.result.anchors->supported.value);
        EXPECT_EQ(2u, out.result.anchors->branches_pruned);
        ASSERT_EQ(2u, out.paths.size());
        EXPECT_EQ("ACGTTGCA", out.paths[0].sequence);
        EXPECT_EQ("TCGTTGCC", out.paths[1].sequence);
    }
}

TEST(PatternSupport, MosaicDiesBeforeL) {
    const Index idx = build_rd(4, { { "A", "a0", "ACGTTGCAA" }, { "B", "b0", "TCGTTGCCA" } });
    Ask ask;
    Outcome out = run_support(idx, "WCGTTGCMA", false, ask);
    check_complete(idx, Oracle(idx).expect("WCGTTGCMA", false, Strands::BOTH, true, true, false),
                   out, ask);
    EXPECT_EQ(Relation::AT_LEAST, out.result.anchors->paths.relation);
    EXPECT_EQ(2u, out.result.anchors->paths.value);
    EXPECT_EQ(Relation::EXACT, out.result.anchors->supported.relation);
    EXPECT_EQ(2u, out.result.anchors->supported.value);
    EXPECT_TRUE(out.result.anchors->pruned_before_completion);
}

TEST(PatternSupport, CrossRecordChain) {
    // B's records b0 and b1 make the column coordinates 2, 3, 4 hold GTA, TAC, ACC: a chain
    // across two records, which is no occurrence; A holds the walk GTACC whole
    const std::vector<Record> records {
        { "A", "a0", "TTGTACCTT" },
        { "B", "b0", "ACGTA" },
        { "B", "b1", "TACC" },
    };
    const Index idx = build_rd(3, records);
    const Oracle oracle(idx);
    Ask record;
    Outcome r = run_support(idx, "GTACC", false, record);
    check_complete(idx, oracle.expect("GTACC", false, Strands::BOTH, true, true, false), r,
                   record);
    ASSERT_EQ(1u, r.paths.size());
    const PathLabels got = labels_of(r.labels->result_fields[0]["labels"]);
    ASSERT_EQ(2u, got.size());
    EXPECT_EQ("record_verified", std::get<0>(got.at("A")));
    EXPECT_EQ(std::set<std::string>({ "0:3" }), std::get<2>(got.at("A")));
    EXPECT_EQ("label_intersection", std::get<0>(got.at("B")));
    EXPECT_TRUE(std::get<2>(got.at("B")).empty());
    EXPECT_EQ("mixed", r.labels->result_fields[0]["support"].asString());

    // without the record mapping: the label level, B's chain listed by its column coordinate
    Ask global;
    global.level = Support::KMER;
    global.records = false;
    Outcome g = run_support(idx, "GTACC", false, global);
    check_complete(idx, oracle.expect("GTACC", false, Strands::BOTH, false, true, true), g,
                   global);
    ASSERT_EQ(1u, g.paths.size());
    const PathLabels gl = labels_of(g.labels->result_fields[0]["labels"]);
    EXPECT_EQ(std::set<std::string>({ "g:2" }), std::get<2>(gl.at("B")));
    EXPECT_NE(std::find(g.labels->notes.begin(), g.labels->notes.end(),
                        "record_bounds_unknown"), g.labels->notes.end());
    // the record level needs the record mapping
    EXPECT_THROW(run_support(idx, "GTACC", false, [] {
        Ask a;
        a.records = false;
        return a;
    }()), std::invalid_argument);
}

TEST(PatternSupport, ALabelSupportsAWalkInOneOrientationOnly) {
    // the walk ATCGG (k = 3: ATC, TCG, CGG) exists through B. A holds ATC and TCG as spelled and
    // CGG only on its other strand (its record CCG): on a BASIC graph A carries the walk's first
    // k-mers on + and its last on -, and supports it in neither orientation
    const std::vector<Record> records {
        { "A", "a0", "ATCG" },
        { "A", "a1", "CCG" },
        { "B", "b0", "ATCGG" },
    };
    const Index idx = build_rd(3, records);
    for (Support level : { Support::KMER, Support::TRACE }) {
        Ask ask;
        ask.level = level;
        ask.strands = Strands::FORWARD;
        Outcome out = run_support(idx, "ATCGG", false, ask);
        check_complete(idx, Oracle(idx).expect("ATCGG", false, Strands::FORWARD,
                                               level == Support::TRACE, true, false), out, ask);
        ASSERT_EQ(1u, out.paths.size());
        const PathLabels got = labels_of(out.labels->result_fields[0]["labels"]);
        EXPECT_EQ(1u, got.size());
        EXPECT_TRUE(got.count("B"));
        EXPECT_FALSE(got.count("A"));
    }
    // on a CANONICAL graph a k-mer and its reverse complement share one row: A carries every
    // k-mer of the walk there (no strand is known), at the label level
    const Index canonical = build<annot::ColumnCompressed<>>(3, records, false,
                                                             DeBruijnGraph::CANONICAL);
    Ask ask;
    ask.level = Support::KMER;
    ask.records = false;
    ask.strands = Strands::FORWARD;
    Outcome out = run_support(canonical, "ATCGG", false, ask);
    check_complete(canonical, Oracle(canonical).expect("ATCGG", false, Strands::FORWARD, false,
                                                       true, false), out, ask);
    ASSERT_EQ(1u, out.paths.size());
    EXPECT_EQ(2u, labels_of(out.labels->result_fields[0]["labels"]).size());
}

TEST(PatternSupport, HomopolymerChainsAreOneRunPerStep) {
    // a record of 3,000 A: its A^k k-mer has 2,990 consecutive coordinates, one run; each step
    // of the 140 merges a run with a run
    const size_t k = 11;
    const Index idx = build_rd(k, { { "H", "h0", std::string(3000, 'A') } });
    Ask ask;
    ask.strands = Strands::FORWARD;
    Outcome out = run_support(idx, std::string(150, 'N'), false, ask);
    check_complete(idx, Oracle(idx).expect(std::string(150, 'N'), false, Strands::FORWARD,
                                           true, true, false), out, ask);
    ASSERT_EQ(1u, out.paths.size());
    const PathLabels got = labels_of(out.labels->result_fields[0]["labels"]);
    EXPECT_EQ(3000u - 150 + 1, std::get<2>(got.at("H")).size());
    // the rows: the anchor's whole row and its coordinate row, read once (the 139 steps after
    // the anchor find it in the cache)
    EXPECT_EQ(2u, out.work.rows);
    EXPECT_EQ(139u, out.work.cache_hits);
    // a merge per step compares a few labels and runs, not the 2,990 coordinates; the reads
    // (each coordinate once) and turning them into runs dominate
    const uint64_t merges = out.units - 2 * (3000 - k + 1);
    EXPECT_LT(merges, 140u * 8 + 2'000);
    // and each of the 139 steps is charged its merge: two labels and two runs compared
    EXPECT_GE(merges, 139u * 4);

    // a deadline in the middle of the search, at every reading of the clock: a stated time
    // stop, of the extension once discovery is done
    uint64_t in_extension = 0;
    for (uint64_t late_after = 0; late_after < 20; ++late_after) {
        SCOPED_TRACE(late_after);
        auto readings = std::make_shared<uint64_t>(0);
        const Clock::time_point start = Clock::now();
        Ask timed = ask;
        timed.time_budget_ms = 1000;
        timed.clock = [=]() {
            return start + std::chrono::seconds(++*readings > late_after ? 100 : 0);
        };
        Outcome t = run_support(idx, std::string(150, 'N'), false, timed);
        if (!t.result.stop) {
            // the search completed; the clock may still have stopped the labels' output
            ASSERT_TRUE(t.labels);
            if (t.labels->stop) {
                EXPECT_EQ(std::make_pair(std::string("output"), std::string("time")),
                          *t.labels->stop);
                EXPECT_TRUE(t.labels->time_limited);
            } else {
                check_complete(idx, Oracle(idx).expect(std::string(150, 'N'), false,
                                                       Strands::FORWARD, true, true, false),
                               t, ask);
            }
            continue;
        }
        EXPECT_EQ(StopReason::TIME, t.result.stop->reason);
        EXPECT_TRUE(t.result.time_limited);
        ASSERT_TRUE(t.x);
        EXPECT_EQ(StopReason::TIME, t.x->cut);
        EXPECT_FALSE(t.x->complete);
        if (t.result.stop->phase == StopPhase::EXTENSION) {
            ++in_extension;
            EXPECT_EQ(Relation::AT_LEAST, t.result.anchors->supported.relation);
            EXPECT_TRUE(t.paths.empty());
        }
    }
    EXPECT_LT(1u, in_extension);
}

TEST(PatternSupport, AMergeTheClockStoppedIsNeverUsed) {
    // N^6000 on a record of 8,000 A: 5,990 steps of a few rounds each, so that the merges'
    // clock (every 4,096 rounds across the calls) reads in the middle of a step. At every
    // reading made late: a stated time stop, never a frame the stop left incomplete used (the
    // next step of the homopolymer finds its row in the cache and steps that frame)
    const size_t k = 11;
    const Index idx = build_rd(k, { { "H", "h0", std::string(8000, 'A') } });
    const std::string text(6000, 'N');
    Ask ask;
    ask.strands = Strands::FORWARD;
    ask.labels = false;
    auto count = std::make_shared<uint64_t>(0);
    {
        Ask a = ask;
        a.time_budget_ms = 1000;
        const Clock::time_point start = Clock::now();
        a.clock = [=]() { ++*count; return start; };
        const Outcome full = run_support(idx, text, false, a);
        ASSERT_FALSE(full.result.stop);
        ASSERT_EQ(1u, full.paths.size());
    }
    uint64_t in_extension = 0;
    for (uint64_t n = 0; n < *count; ++n) {
        SCOPED_TRACE(n);
        auto readings = std::make_shared<uint64_t>(0);
        Ask a = ask;
        a.time_budget_ms = 1000;
        const Clock::time_point start = Clock::now();
        a.clock = [=]() {
            return start + std::chrono::seconds(++*readings > n ? 100 : 0);
        };
        const Outcome out = run_support(idx, text, false, a);
        if (!out.result.stop) {
            // a search the clock did not stop is the complete one (only the output of its
            // labels may have been stopped)
            EXPECT_FALSE(out.result.time_limited);
            EXPECT_EQ(Relation::EXACT, out.result.anchors->supported.relation);
            EXPECT_EQ(1u, out.result.anchors->supported.value);
            EXPECT_EQ(1u, out.paths.size());
            continue;
        }
        EXPECT_EQ(StopReason::TIME, out.result.stop->reason);
        if (out.result.stop->phase == StopPhase::EXTENSION) {
            ++in_extension;
            EXPECT_TRUE(out.paths.empty());
            EXPECT_EQ(0u, out.result.anchors->supported.value);
        }
    }
    EXPECT_LT(10u, in_extension);
}


// ------------------------------------------------------------------ random indexes

std::string random_bases(std::mt19937 &rng, size_t n) {
    std::string s;
    for (size_t i = 0; i < n; ++i) {
        s.push_back("ACGT"[rng() % 4]);
    }
    return s;
}

// 2-6 columns of 1-3 records each, half of the records with pieces of earlier ones planted
// (either strand): repeats, bubbles and mosaics
std::vector<Record> random_records(std::mt19937 &rng, size_t k) {
    std::vector<Record> records;
    const size_t columns = 2 + rng() % 5;
    for (size_t c = 0; c < columns; ++c) {
        const size_t n = 1 + rng() % 3;
        for (size_t r = 0; r < n; ++r) {
            std::string seq = random_bases(rng, k + 2 + rng() % 26);
            if (records.size() && rng() % 2) {
                const std::string &from = records[rng() % records.size()].seq;
                const size_t len = std::min(from.size(), k + rng() % (k + 4));
                std::string piece = from.substr(rng() % (from.size() - len + 1), len);
                if (rng() % 3 == 0)
                    piece = rc(piece);
                const size_t at = rng() % (seq.size() + 1);
                seq.insert(at, piece);
                // a mismatch inside the copy now and then: a bubble
                if (rng() % 2)
                    seq[at + rng() % piece.size()] = "ACGT"[rng() % 4];
            }
            records.push_back(Record { "c" + std::to_string(c),
                                       "c" + std::to_string(c) + "_r" + std::to_string(r),
                                       seq });
        }
    }
    return records;
}

// a pattern of L in (k, 3k] bases cut from a record (either strand), some positions made
// ambiguous; or a peptide translated from a record's bases, now and then with an X
std::pair<std::string, bool> random_pattern(std::mt19937 &rng, const std::vector<Record> &records,
                                            size_t k) {
    for (int attempt = 0; attempt < 100; ++attempt) {
        const std::string &seq = records[rng() % records.size()].seq;
        if (rng() % 4 == 0) {
            const size_t m = (k + 3) / 3 + rng() % 3;
            if (3 * m > seq.size())
                continue;
            std::string peptide = translate(seq.substr(rng() % (seq.size() - 3 * m + 1), 3 * m));
            if (rng() % 2)
                peptide[rng() % peptide.size()] = 'X';
            return { peptide, true };
        }
        const size_t L = k + 1 + rng() % (2 * k);
        if (L > seq.size())
            continue;
        std::string p = seq.substr(rng() % (seq.size() - L + 1), L);
        if (rng() % 3 == 0)
            p = rc(p);
        for (char &c : p) {
            if (rng() % 7 == 0) {
                static const std::string kAmbiguous = "RYSWKMBDHVN";
                for (int t = 0; t < 10; ++t) {
                    const char code = kAmbiguous[rng() % kAmbiguous.size()];
                    if (admits(code, c)) {
                        c = code;
                        break;
                    }
                }
            }
        }
        return { p, false };
    }
    return { std::string(k + 1, 'N'), false };
}

TEST(PatternSupport, RandomIndexesAgainstTheOracles) {
    std::mt19937 rng(20261009);
    uint64_t cases = 0, supported = 0, pruned_early = 0, peptides = 0, mixed = 0;
    for (int t = 0; t < 120; ++t) {
        const size_t k = 3 + rng() % 5;
        const std::vector<Record> records = random_records(rng, k);
        const Index idx = build_rd(k, records);
        const Oracle oracle(idx);
        for (int j = 0; j < 3; ++j) {
            const auto [text, peptide] = random_pattern(rng, records, k);
            for (Support level : { Support::TRACE, Support::KMER }) {
                for (Strands strands : { Strands::BOTH, Strands::FORWARD }) {
                    SCOPED_TRACE(std::to_string(t) + " " + text + " level "
                                 + std::to_string(level == Support::TRACE)
                                 + " strands " + to_string(strands));
                    Ask ask;
                    ask.level = level;
                    ask.strands = strands;
                    const Expected e = oracle.expect(text, peptide, strands,
                                                     level == Support::TRACE,
                                                     level == Support::TRACE, false);
                    const Outcome out = run_support(idx, text, peptide, ask);
                    check_complete(idx, e, out, ask);
                    // the rows: the anchors' whole rows, then each row the walk entered once (the
                    // allotment holds every row of these indexes); at the record level the
                    // anchors' coordinate rows too
                    std::set<std::string> walk_rows = e.pushed_rows;
                    if (level == Support::TRACE) {
                        if (e.anchors)
                            walk_rows.insert(e.anchor_rows.begin(), e.anchor_rows.end());
                        EXPECT_EQ(e.anchor_rows.size() + walk_rows.size(), out.work.rows);
                    } else {
                        walk_rows.insert(e.anchor_rows.begin(), e.anchor_rows.end());
                        EXPECT_EQ(walk_rows.size(), out.work.rows);
                    }
                    EXPECT_EQ(e.anchor_rows.size(), out.work.anchor_rows);
                    EXPECT_EQ(0u, out.work.evictions);
                    EXPECT_LE(out.peak, ask.max_memory);
                    // what stays in the account: the returned paths' descriptors and their
                    // labels built for the answer (their summary), nothing of the search
                    if (out.paths.empty()) {
                        EXPECT_EQ(0u, out.held_after);
                    }
                    ++cases;
                    supported += e.supported.size();
                    pruned_early += e.pruned_early;
                    peptides += peptide;
                    for (const ExpectedPath &p : e.supported) {
                        bool v = false, u = false;
                        for (const auto &[c, view] : p.labels) {
                            (std::get<0>(view) == "record_verified" ? v : u) = true;
                        }
                        mixed += v && u;
                    }
                }
            }
        }
    }
    // the cases exercise what they are for
    EXPECT_LT(1000u, cases);
    EXPECT_LT(300u, supported);
    EXPECT_LT(50u, pruned_early);
    EXPECT_LT(100u, peptides);
    EXPECT_LT(10u, mixed);
}

/**
 * long_search "paths" with output.labels "all" (every walk, its labels read after the
 * extension and verified by a join of the coordinate lists), partial with every cap raised:
 * the entry of |text|.
 */
Json::Value increment4(const Index &idx, const std::string &text, bool peptide, Strands strands,
                       bool occurrences, bool records) {
    PatternLimits l;
    l.min_information_bits = 0;
    l.max_anchors = l.max_paths = 1'000'000;
    l.max_labels = l.max_occurrences_per_label = 1'000'000;
    l.max_labels_per_anchor = 1'000;
    const std::string body = std::string("{\"patterns\": [{\"")
            + (peptide ? "protein" : "iupac") + "\": \"" + text + "\"}], "
            "\"long_search\": \"paths\", \"mode\": \"partial\", \"strands\": \""
            + to_string(strands) + "\", \"max_paths\": 1000000, \"max_anchors\": 1000000, "
            "\"max_labels\": 1000000, \"max_occurrences_per_label\": 1000000, "
            "\"output\": {\"labels\": \"all\", \"occurrences\": "
            + (occurrences ? "true" : "false") + "}}";
    RetrievalHooks hooks;
    hooks.coord_to_header = records ? idx.cth.get() : nullptr;
    return process_pattern_request(parse_pattern_body(body), *idx.anno, l, "rel", nullptr,
                                   nullptr, nullptr, &hooks)["patterns"][0];
}

TEST(PatternSupport, TheAccountKeepsTheReturnedPathsOnly) {
    // after a pattern the account holds what the answer holds: without labels the returned
    // paths' descriptors (as many bytes as path_descriptor_bytes says), nothing of the search
    // (its frames, rows, the row cache's allotment, the dictionary); with labels, more
    std::mt19937 rng(13);
    uint64_t returned = 0;
    for (int t = 0; t < 40; ++t) {
        const size_t k = 3 + rng() % 4;
        const std::vector<Record> records = random_records(rng, k);
        const Index idx = build_rd(k, records);
        const auto [text, peptide] = random_pattern(rng, records, k);
        const size_t L = peptide ? 3 * text.size() : text.size();
        for (Support level : { Support::TRACE, Support::KMER }) {
            SCOPED_TRACE(text);
            Ask ask;
            ask.level = level;
            ask.labels = false;
            const Outcome out = run_support(idx, text, peptide, ask);
            ASSERT_TRUE(out.x);
            EXPECT_EQ(out.x->returned * path_descriptor_bytes(k, L), out.held_after);
            EXPECT_LE(out.peak, ask.max_memory);
            // the search held the row cache's allotment and its frames at once
            if (out.work.frames_peak_bytes) {
                EXPECT_GE(out.peak, PathTracker::cache_allotment(ask.max_memory)
                                        + out.work.frames_peak_bytes);
            }
            returned += out.x->returned;
            Ask labelled = ask;
            labelled.labels = true;
            const Outcome l = run_support(idx, text, peptide, labelled);
            EXPECT_LE(out.held_after, l.held_after);
        }
    }
    EXPECT_LT(20u, returned);
}

TEST(PatternSupport, DifferentialAgainstTheLabelsOfAllPaths) {
    // record level: the supported paths are exactly the paths with a record_verified label,
    // each with the same labels, supports and occurrences; label level: the paths with a label
    // at all, with the same labels
    std::mt19937 rng(77);
    uint64_t compared = 0, dropped = 0;
    for (int t = 0; t < 60; ++t) {
        const size_t k = 3 + rng() % 5;
        const std::vector<Record> records = random_records(rng, k);
        const Index idx = build_rd(k, records);
        for (int j = 0; j < 3; ++j) {
            const auto [text, peptide] = random_pattern(rng, records, k);
            for (Support level : { Support::TRACE, Support::KMER }) {
                SCOPED_TRACE(std::to_string(t) + " " + text);
                const bool trace = level == Support::TRACE;
                const Json::Value all = increment4(idx, text, peptide, Strands::BOTH, trace, true);
                ASSERT_FALSE(all.isMember("error")) << all;
                std::vector<std::pair<std::string, PathLabels>> want;
                for (const Json::Value &r : all["results"]) {
                    const PathLabels pl = labels_of(r["labels"]);
                    bool keep = !pl.empty();
                    if (trace) {
                        keep = false;
                        for (const auto &[c, view] : pl) {
                            keep |= std::get<0>(view) == "record_verified";
                        }
                    }
                    if (keep) {
                        want.emplace_back(r["sequence"].asString(), pl);
                    } else {
                        ++dropped;
                    }
                }
                Ask ask;
                ask.level = level;
                ask.occurrences = trace;
                const Outcome out = run_support(idx, text, peptide, ask);
                ASSERT_FALSE(out.result.stop);
                ASSERT_EQ(want.size(), out.paths.size());
                for (size_t i = 0; i < want.size(); ++i) {
                    EXPECT_EQ(want[i].first, out.paths[i].sequence);
                    EXPECT_EQ(want[i].second,
                              labels_of(out.labels->result_fields[i]["labels"]));
                }
                compared += want.size();
            }
        }
    }
    EXPECT_LT(200u, compared);
    EXPECT_LT(20u, dropped);
}

TEST(PatternSupport, GlobalPlacementAgainstTheLabelsOfAllPaths) {
    // coordinates without the record mapping: the label level, the labels' chains listed by
    // their column coordinates as long_search "paths" lists them
    std::mt19937 rng(5);
    uint64_t compared = 0;
    for (int t = 0; t < 30; ++t) {
        const size_t k = 3 + rng() % 4;
        const std::vector<Record> records = random_records(rng, k);
        const Index idx = build_rd(k, records);
        const Oracle oracle(idx);
        for (int j = 0; j < 3; ++j) {
            const auto [text, peptide] = random_pattern(rng, records, k);
            SCOPED_TRACE(std::to_string(t) + " " + text);
            Ask ask;
            ask.level = Support::KMER;
            ask.records = false;
            const Outcome out = run_support(idx, text, peptide, ask);
            check_complete(idx, oracle.expect(text, peptide, Strands::BOTH, false, true, true),
                           out, ask);
            const Json::Value all = increment4(idx, text, peptide, Strands::BOTH, true, false);
            size_t i = 0;
            for (const Json::Value &r : all["results"]) {
                const PathLabels pl = labels_of(r["labels"]);
                if (pl.empty())
                    continue;
                ASSERT_LT(i, out.paths.size());
                EXPECT_EQ(r["sequence"].asString(), out.paths[i].sequence);
                EXPECT_EQ(pl, labels_of(out.labels->result_fields[i]["labels"]));
                ++i;
                ++compared;
            }
            EXPECT_EQ(i, out.paths.size());
        }
    }
    EXPECT_LT(50u, compared);
}

TEST(PatternSupport, CanonicalAndPrimaryGraphsAtTheLabelLevel) {
    std::mt19937 rng(31);
    uint64_t supported = 0;
    for (auto mode : { DeBruijnGraph::CANONICAL, DeBruijnGraph::PRIMARY }) {
        for (int t = 0; t < 25; ++t) {
            const size_t k = 3 + rng() % 4;
            const std::vector<Record> records = random_records(rng, k);
            const Index idx = build<annot::ColumnCompressed<>>(k, records, false, mode);
            const Oracle oracle(idx);
            for (int j = 0; j < 3; ++j) {
                const auto [text, peptide] = random_pattern(rng, records, k);
                SCOPED_TRACE(std::to_string(t) + " " + text);
                Ask ask;
                ask.level = Support::KMER;
                ask.records = false;
                const Outcome out = run_support(idx, text, peptide, ask);
                const Expected e = oracle.expect(text, peptide, Strands::BOTH, false, true,
                                                 false);
                check_complete(idx, e, out, ask);
                supported += e.supported.size();
                // unbudgeted reads, stated; no coordinates read
                EXPECT_NE(std::find(out.labels->notes.begin(), out.labels->notes.end(),
                                    "annotation_unbudgeted"), out.labels->notes.end());
                EXPECT_NE(std::find(out.labels->notes.begin(), out.labels->notes.end(),
                                    "label_intersection_only"), out.labels->notes.end());
            }
            // no record level where coordinates carry no strand
            Ask trace;
            EXPECT_THROW(run_support(idx, "ACGTACGT", false, trace), std::invalid_argument);
        }
    }
    EXPECT_LT(50u, supported);
}

TEST(PatternSupport, UnbudgetedCoordinatesAndNoCoordinates) {
    std::mt19937 rng(8);
    for (int t = 0; t < 20; ++t) {
        const size_t k = 3 + rng() % 4;
        const std::vector<Record> records = random_records(rng, k);
        const Index coords = build<annot::ColumnCompressed<>>(k, records, true);
        const Index plain = build<annot::ColumnCompressed<>>(k, records, false);
        const auto [text, peptide] = random_pattern(rng, records, k);
        SCOPED_TRACE(text);
        Ask trace;
        check_complete(coords, Oracle(coords).expect(text, peptide, Strands::BOTH, true, true,
                                                     false),
                       run_support(coords, text, peptide, trace), trace);
        Ask label;
        label.level = Support::KMER;
        check_complete(plain, Oracle(plain).expect(text, peptide, Strands::BOTH, false, true,
                                                   false),
                       run_support(plain, text, peptide, label), label);
        EXPECT_THROW(run_support(plain, text, peptide, trace), std::invalid_argument);
    }
}


// ------------------------------------------------------------------ release, modes, caps

// a random index and pattern with at least |min| supported paths at the record level
std::tuple<Index, std::string, bool, Expected> rich_case(uint32_t seed, size_t min) {
    std::mt19937 rng(seed);
    for (int t = 0; ; ++t) {
        const size_t k = 3 + rng() % 3;
        std::vector<Record> records = random_records(rng, k);
        Index idx = build_rd(k, records);
        for (int j = 0; j < 10; ++j) {
            auto [text, peptide] = random_pattern(rng, records, k);
            Expected e = Oracle(idx).expect(text, peptide, Strands::BOTH, true, true, false);
            if (e.supported.size() >= min && e.pruned_early)
                return { std::move(idx), text, peptide, std::move(e) };
        }
        EXPECT_LT(t, 1000);
    }
}

TEST(PatternSupport, ModesAndThresholds) {
    auto [idx, text, peptide, e] = rich_case(3, 4);
    const uint64_t S = e.supported.size();
    // count: the counts, nothing kept
    {
        Ask ask;
        ask.mode = Mode::COUNT;
        Outcome out = run_support(idx, text, peptide, ask);
        check_complete(idx, e, out, ask);
    }
    // all_or_count: all of them, or none above max_paths (count_above_threshold)
    for (uint64_t max_paths : { S, S - 1 }) {
        Ask ask;
        ask.mode = Mode::ALL_OR_COUNT;
        ask.max_paths = max_paths;
        Outcome out = run_support(idx, text, peptide, ask);
        ASSERT_TRUE(out.x);
        EXPECT_EQ(S, out.result.anchors->supported.value);
        if (max_paths == S) {
            EXPECT_TRUE(out.x->complete);
            EXPECT_EQ(S, out.paths.size());
        } else {
            EXPECT_EQ(Withheld::COUNT_ABOVE_THRESHOLD, out.x->withheld);
            EXPECT_TRUE(out.paths.empty());
            EXPECT_EQ(0u, out.held_after);
        }
    }
    // partial: the first max_paths in answer order, the cut stated
    for (uint64_t max_paths = 0; max_paths <= S; ++max_paths) {
        Ask ask;
        ask.max_paths = max_paths;
        Outcome out = run_support(idx, text, peptide, ask);
        ASSERT_EQ(max_paths, out.paths.size());
        for (size_t i = 0; i < out.paths.size(); ++i) {
            EXPECT_EQ(e.supported[i].sequence, out.paths[i].sequence);
        }
        EXPECT_EQ(max_paths == S, out.x->complete);
        if (max_paths < S) {
            EXPECT_EQ(StopReason::MAX_PATHS, out.x->cut);
        }
        ASSERT_TRUE(out.labels);
        EXPECT_EQ(max_paths, out.labels->result_fields.size());
    }
    // stop_at_threshold: the search ends once more than max_paths supported paths are complete
    for (Mode mode : { Mode::PARTIAL, Mode::ALL_OR_COUNT }) {
        Ask ask;
        ask.mode = mode;
        ask.max_paths = S - 2;
        ask.stop_at_threshold = true;
        Outcome out = run_support(idx, text, peptide, ask);
        ASSERT_TRUE(out.result.stop);
        EXPECT_EQ(StopReason::MAX_PATHS, out.result.stop->reason);
        EXPECT_EQ(Relation::AT_LEAST, out.result.anchors->supported.relation);
        EXPECT_EQ(S - 1, out.result.anchors->supported.value);
        if (mode == Mode::PARTIAL) {
            EXPECT_EQ(StopReason::MAX_PATHS, out.x->cut);
            EXPECT_EQ(S - 2, out.paths.size());
        } else {
            EXPECT_EQ(Withheld::THRESHOLD_CROSSED, out.x->withheld);
        }
    }
    // the extension not admitted
    {
        Ask ask;
        ask.max_anchors = 0;
        Outcome out = run_support(idx, text, peptide, ask);
        EXPECT_EQ(Extension::NOT_ADMITTED, out.result.anchors->extension);
        EXPECT_EQ(Withheld::ANCHORS_ABOVE_THRESHOLD, out.x->withheld);
        EXPECT_EQ(0u, out.work.rows);
    }
}

TEST(PatternSupport, RequireVerifiedListsTheVerifiedLabels) {
    std::mt19937 rng(17);
    uint64_t excluded = 0;
    for (int t = 0; t < 80; ++t) {
        const size_t k = 3 + rng() % 4;
        const std::vector<Record> records = random_records(rng, k);
        const Index idx = build_rd(k, records);
        const auto [text, peptide] = random_pattern(rng, records, k);
        Ask ask;
        ask.require_verified = true;
        const Expected e = Oracle(idx).expect(text, peptide, Strands::BOTH, true, true, false);
        const Outcome out = run_support(idx, text, peptide, ask);
        check_complete(idx, e, out, ask);
        for (size_t i = 0; i < out.paths.size(); ++i) {
            uint64_t unverified = 0;
            for (const auto &[c, view] : e.supported[i].labels) {
                unverified += std::get<0>(view) != "record_verified";
            }
            EXPECT_EQ(e.supported[i].labels.size(),
                      out.labels->result_fields[i]["labels_total"].asUInt64());
            EXPECT_EQ(unverified,
                      out.labels->result_fields[i]["labels_excluded_unverified"].asUInt64());
            excluded += unverified;
        }
    }
    EXPECT_LT(5u, excluded);
    // only at the record level
    const Index idx = build_rd(3, { { "A", "a0", "ACGT" } });
    Ask label;
    label.level = Support::KMER;
    label.require_verified = true;
    EXPECT_THROW(run_support(idx, "ACGT", false, label), std::invalid_argument);
}


// ------------------------------------------------------------------ budgets and stops

// a stopped search against the complete one: the supported paths completed before the stop are
// a prefix of the complete search's, in answer order (partial lists them), at least as many as
// the counts claim, and the labels of each as the complete search gives them
void check_prefix(const Outcome &full, const Outcome &stopped) {
    ASSERT_TRUE(stopped.result.anchors);
    const AnchorCounts &a = *stopped.result.anchors;
    if (!stopped.result.stop) {
        // the search completed: its counts are the complete ones; the list is all of them, or
        // (the allowance could not hold them) a prefix
        if (!stopped.output_cut) {
            EXPECT_EQ(full.paths.size(), stopped.paths.size());
        }
        EXPECT_EQ(full.result.anchors->supported.value, a.supported.value);
        EXPECT_EQ(Relation::EXACT, a.supported.relation);
    } else {
        EXPECT_NE(Relation::EXACT, a.supported.relation);
    }
    ASSERT_LE(stopped.paths.size(), full.paths.size());
    EXPECT_LE(a.supported.value, full.result.anchors->supported.value);
    for (size_t i = 0; i < stopped.paths.size(); ++i) {
        EXPECT_EQ(full.paths[i].sequence, stopped.paths[i].sequence);
        EXPECT_EQ(full.paths[i].orientation, stopped.paths[i].orientation);
        if (stopped.labels && i < stopped.labels->result_fields.size()
                && stopped.labels->result_fields[i]["labels"].isArray()) {
            EXPECT_EQ(labels_of(full.labels->result_fields[i]["labels"]),
                      labels_of(stopped.labels->result_fields[i]["labels"]));
        }
    }
}

TEST(PatternSupport, WorkBudgetSweep) {
    auto [idx, text, peptide, e] = rich_case(11, 3);
    Ask ask;
    const Outcome full = run_support(idx, text, peptide, ask);
    check_complete(idx, e, full, ask);
    uint64_t last = 0, stops = 0;
    for (uint64_t work = 0; work <= full.units + 1; ++work) {
        SCOPED_TRACE(work);
        Ask a = ask;
        a.max_work = work;
        const Outcome out = run_support(idx, text, peptide, a);
        check_prefix(full, out);
        if (out.result.stop) {
            ++stops;
            EXPECT_EQ(StopPhase::EXTENSION, out.result.stop->phase);
            EXPECT_EQ(StopReason::EXTERNAL, out.result.stop->reason);
            EXPECT_EQ("max_annotation_work", out.stop);
            EXPECT_EQ(StopReason::EXTERNAL, out.x->cut);
            // the units reached the budget (a read or a merge is gated before it, charged
            // after it)
            EXPECT_GE(out.units, work);
            EXPECT_LE(out.units, full.units);
        }
        // no read begins once the work reached the budget
        for (uint64_t before : out.units_before_reads) {
            EXPECT_LT(before, work);
        }
        if (out.result.stop) {
        } else {
            // completed: the last read or merge may pass the budget, never another after it
            EXPECT_EQ(full.units, out.units);
        }
        // more work never completes fewer supported paths
        EXPECT_LE(last, out.result.anchors->supported.value);
        last = out.result.anchors->supported.value;
    }
    EXPECT_LT(10u, stops);
    // all_or_count: a stop withholds everything, stated as the tracker's
    Ask all = ask;
    all.mode = Mode::ALL_OR_COUNT;
    all.max_work = full.units / 2;
    const Outcome out = run_support(idx, text, peptide, all);
    EXPECT_EQ(Withheld::EXTERNAL, out.x->withheld);
}

TEST(PatternSupport, StepsAndTimeSweeps) {
    auto [idx, text, peptide, e] = rich_case(23, 3);
    Ask ask;
    const Outcome full = run_support(idx, text, peptide, ask);
    const uint64_t steps = full.result.work.steps;
    for (uint64_t max_steps = 1; max_steps <= steps + 1; ++max_steps) {
        Ask a = ask;
        a.max_steps = max_steps;
        const Outcome out = run_support(idx, text, peptide, a);
        check_prefix(full, out);
        if (out.result.stop) {
            EXPECT_EQ(StopReason::MAX_STEPS, out.result.stop->reason);
        }
    }
    // a virtual clock late from its n-th reading on, for every n up to the readings of the
    // complete search
    auto count = std::make_shared<uint64_t>(0);
    {
        Ask a = ask;
        a.time_budget_ms = 1000;
        const Clock::time_point start = Clock::now();
        a.clock = [=]() { ++*count; return start; };
        run_support(idx, text, peptide, a);
    }
    uint64_t time_stops = 0;
    for (uint64_t n = 0; n <= *count; ++n) {
        SCOPED_TRACE(n);
        auto readings = std::make_shared<uint64_t>(0);
        Ask a = ask;
        a.time_budget_ms = 1000;
        const Clock::time_point start = Clock::now();
        a.clock = [=]() {
            return start + std::chrono::seconds(++*readings > n ? 100 : 0);
        };
        const Outcome out = run_support(idx, text, peptide, a);
        if (!out.result.stop || out.result.stop->phase != StopPhase::EXTENSION)
            continue;
        ++time_stops;
        check_prefix(full, out);
        EXPECT_EQ(StopReason::TIME, out.result.stop->reason);
        EXPECT_TRUE(out.result.time_limited);
        EXPECT_EQ(StopReason::TIME, out.x->cut);
    }
    EXPECT_LT(0u, time_stops);
}

TEST(PatternSupport, MemorySweep) {
    auto [idx, text, peptide, e] = rich_case(41, 8);
    Ask ask;
    const Outcome full = run_support(idx, text, peptide, ask);
    std::map<std::string, uint64_t> outcomes;
    for (uint64_t mem = 64; mem <= full.peak * 4 + 4096; mem += std::max<uint64_t>(16, mem / 40)) {
        SCOPED_TRACE(mem);
        Ask a = ask;
        a.max_memory = mem;
        const Outcome out = run_support(idx, text, peptide, a);
        // never past the account
        EXPECT_LE(out.peak, mem);
        // the kept paths within half of the account, the other half for the search
        EXPECT_LE(out.paths.size() * path_descriptor_bytes(idx.k, full.paths.size()
                                                                  ? full.paths[0].sequence.size()
                                                                  : 0),
                  mem / 2);
        check_prefix(full, out);
        if (out.result.stop) {
            EXPECT_EQ("max_memory", out.stop);
            ++outcomes["search"];
        } else if (out.output_cut) {
            EXPECT_LT(out.paths.size(), full.paths.size());
            ASSERT_TRUE(out.labels);
            ASSERT_TRUE(out.labels->stop);
            EXPECT_EQ(std::make_pair(std::string("output"), std::string("max_memory")),
                      *out.labels->stop);
            EXPECT_EQ("max_memory", out.labels->cut.value_or(""));
            ++outcomes["list"];
        } else if (out.labels && out.labels->stop) {
            EXPECT_EQ("output", out.labels->stop->first);
            ++outcomes["labels"];
        } else {
            EXPECT_EQ(full.paths.size(), out.paths.size());
            ++outcomes["complete"];
        }
        // a refused row is stated
        if (out.result.stop && out.rows_refused.size()) {
            EXPECT_EQ("extension", out.rows_refused[0]["phase"].asString());
            EXPECT_EQ("max_memory", out.rows_refused[0]["reason"].asString());
        }
    }
    EXPECT_TRUE(outcomes.count("search"));
    EXPECT_TRUE(outcomes.count("complete"));
    // all_or_count: a list the allowance cannot hold is withheld whole (output_budget)
    for (uint64_t mem = 64; mem <= full.peak * 4; mem += std::max<uint64_t>(16, mem / 40)) {
        Ask a = ask;
        a.mode = Mode::ALL_OR_COUNT;
        a.max_memory = mem;
        const Outcome out = run_support(idx, text, peptide, a);
        EXPECT_LE(out.peak, mem);
        if (!out.result.stop && out.output_cut) {
            EXPECT_TRUE(out.paths.empty());
            ASSERT_TRUE(out.labels);
            EXPECT_EQ("output_budget", out.labels->withheld.value_or(""));
        }
    }
}

TEST(PatternSupport, RefusedRowsStopTheSearch) {
    // every budget-aware read refused from its n-th charge on: the search stops at the first
    // refused row (max_memory, stated), never takes it as empty
    auto [idx, text, peptide, e] = rich_case(53, 2);
    Ask ask;
    const Outcome full = run_support(idx, text, peptide, ask);
    uint64_t refused = 0;
    for (uint64_t n = 0; n < 400; n += 3) {
        Ask a = ask;
        a.deny = [n](uint64_t ordinal) { return ordinal >= n; };
        const Outcome out = run_support(idx, text, peptide, a);
        check_prefix(full, out);
        if (!out.result.stop)
            continue;
        ++refused;
        EXPECT_EQ("max_memory", out.stop);
        ASSERT_EQ(1u, out.rows_refused.size());
        EXPECT_EQ("extension", out.rows_refused[0]["phase"].asString());
    }
    EXPECT_LT(5u, refused);
}

TEST(PatternSupport, TheRowCacheChangesNoAnswer) {
    // a small account: a small allotment, rows evicted and read again (charged again); the
    // supported paths and their labels are those of a large account
    std::mt19937 rng(99);
    uint64_t evicted = 0;
    for (int t = 0; t < 40; ++t) {
        const size_t k = 3 + rng() % 3;
        const std::vector<Record> records = random_records(rng, k);
        const Index idx = build_rd(k, records);
        const auto [text, peptide] = random_pattern(rng, records, k);
        Ask large;
        const Outcome full = run_support(idx, text, peptide, large);
        for (uint64_t mem : { 2'000, 4'000, 8'000, 16'000 }) {
            Ask small = large;
            small.max_memory = mem;
            const Outcome out = run_support(idx, text, peptide, small);
            EXPECT_LE(out.peak, mem);
            if (out.result.stop || out.output_cut || (out.labels && out.labels->stop))
                continue;
            check_prefix(full, out);
            EXPECT_EQ(full.paths.size(), out.paths.size());
            EXPECT_GE(out.work.rows, full.work.rows);
            EXPECT_GE(out.units, full.units);
            if (!out.work.evictions && !out.work.uncached_rows) {
                EXPECT_EQ(full.work.rows, out.work.rows);
                EXPECT_EQ(full.units, out.units);
            } else if (out.work.rows > full.work.rows) {
                // rows read again after an eviction, charged again
                ++evicted;
                EXPECT_LT(full.units, out.units);
            }
        }
    }
    EXPECT_LT(5u, evicted);
}

TEST(PatternSupport, Deterministic) {
    auto [idx, text, peptide, e] = rich_case(61, 3);
    Ask ask;
    const Outcome a = run_support(idx, text, peptide, ask);
    const Outcome b = run_support(idx, text, peptide, ask);
    ASSERT_EQ(a.paths.size(), b.paths.size());
    for (size_t i = 0; i < a.paths.size(); ++i) {
        EXPECT_EQ(a.paths[i].sequence, b.paths[i].sequence);
        EXPECT_EQ(a.paths[i].path, b.paths[i].path);
        EXPECT_EQ(a.labels->result_fields[i], b.labels->result_fields[i]);
    }
    EXPECT_EQ(a.units, b.units);
    EXPECT_EQ(a.peak, b.peak);
    EXPECT_EQ(a.labels->fields, b.labels->fields);
    EXPECT_EQ(a.labels->labels_count, b.labels->labels_count);
    EXPECT_EQ(a.labels->occurrences_count, b.labels->occurrences_count);
}

TEST(PatternSupport, LabelsAnswerCountsAndPartialCuts) {
    auto [idx, text, peptide, e] = rich_case(71, 3);
    Ask ask;
    const Outcome full = run_support(idx, text, peptide, ask);
    const LabelsAnswer &l = *full.labels;
    // counts.labels over the returned paths, by support; by_label in label order
    std::map<std::string, std::pair<uint64_t, uint64_t>> per_label;
    std::map<std::string, std::set<std::string>> unions;
    for (const ExpectedPath &p : e.supported) {
        for (const auto &[c, view] : p.labels) {
            per_label[c].first++;
            per_label[c].second += std::get<0>(view) == "record_verified";
            const auto &occ = std::get<2>(view);
            unions[c].insert(occ.begin(), occ.end());
        }
    }
    EXPECT_EQ(per_label.size(), l.labels_count["value"].asUInt64());
    EXPECT_EQ("exact", l.labels_count["relation"].asString());
    uint64_t verified = 0, occurrences = 0;
    for (const auto &[c, v] : per_label) {
        verified += v.second > 0;
        occurrences += unions[c].size();
    }
    EXPECT_EQ(verified, l.labels_count["by_support"]["record_verified"]["value"].asUInt64());
    EXPECT_EQ(occurrences, l.occurrences_count["value"].asUInt64());
    ASSERT_EQ(per_label.size(), l.fields["by_label"].size());
    uint64_t previous = std::numeric_limits<uint64_t>::max();
    for (const Json::Value &b : l.fields["by_label"]) {
        const std::string c = b["column"].asString();
        EXPECT_EQ(per_label[c].first, b["paths"]["value"].asUInt64());
        EXPECT_EQ(per_label[c].second, b["paths_record_verified"]["value"].asUInt64());
        EXPECT_EQ(unions[c].size(), b["occurrences"]["value"].asUInt64());
        EXPECT_LE(per_label[c].first, previous);
        previous = per_label[c].first;
    }
    // partial's caps: one label listed, one occurrence per label
    Ask capped = ask;
    capped.max_labels = 1;
    capped.max_occurrences = 1;
    const Outcome c = run_support(idx, text, peptide, capped);
    if (per_label.size() > 1) {
        EXPECT_EQ("max_labels", c.labels->fields["labels_cut"]["reason"].asString());
    }
    for (const Json::Value &f : c.labels->result_fields) {
        EXPECT_LE(f["labels"].size(), 1u);
        for (const Json::Value &label : f["labels"]) {
            EXPECT_LE(label["occurrence_list"].size(), 1u);
        }
    }
    EXPECT_FALSE(c.labels->complete);
}


// ------------------------------------------------------------------ the engine's rules

// a tracker recording what the engine asks: every walk supported
class Recording : public SupportTracker {
  public:
    std::vector<std::vector<node_index>> prepared;
    std::vector<node_index> opened;
    bool stop_in_prepare = false;
    uint64_t depth = 0;

    bool prepare(const std::vector<Context> &anchors) override {
        EXPECT_TRUE(opened.empty());
        std::vector<node_index> nodes;
        for (const Context &a : anchors) {
            nodes.push_back(a.node);
        }
        prepared.push_back(nodes);
        return !stop_in_prepare;
    }
    Verdict open(const SearchState &anchor) override {
        EXPECT_EQ(0u, depth);
        opened.push_back(anchor.node);
        ++depth;
        return Verdict::ALIVE;
    }
    Verdict push(const SearchState &) override {
        ++depth;
        return Verdict::ALIVE;
    }
    void pop() override { --depth; }
    Verdict complete(const PathView &) override { return Verdict::ALIVE; }
    const char* stop_reason() const override { return stop_in_prepare ? "test" : nullptr; }
};

TEST(PatternSupport, PrepareSeesEveryAnchorBeforeTheFirstOpens) {
    std::mt19937 rng(101);
    uint64_t with_anchors = 0;
    for (int t = 0; t < 40; ++t) {
        const size_t k = 3 + rng() % 4;
        const std::vector<Record> records = random_records(rng, k);
        const Index idx = build_rd(k, records);
        const auto [text, peptide] = random_pattern(rng, records, k);
        PatternSearch search(idx.anno->get_graph());
        Request rq;
        rq.extend_paths = true;
        rq.min_information_bits = 0;
        rq.max_anchors = 1'000'000;
        Recording tracker;
        rq.support = &tracker;
        Budget budget(1'000'000'000, Deadline(Clock::now(), 1e9, 0));
        const Result r = search.count(pattern_of(text, peptide), rq, budget);
        ASSERT_TRUE(r.anchors);
        if (r.anchors->extension == Extension::NO_ANCHORS) {
            EXPECT_TRUE(tracker.prepared.empty());
            continue;
        }
        ++with_anchors;
        // once, with exactly the anchors opened, in their order
        ASSERT_EQ(1u, tracker.prepared.size());
        EXPECT_EQ(tracker.opened, tracker.prepared[0]);
        EXPECT_EQ(r.anchors->total.value, tracker.opened.size());
        EXPECT_EQ(0u, tracker.depth);

        // not admitted: never asked
        Recording none;
        Request small = rq;
        small.max_anchors = 0;
        small.support = &none;
        Budget b2(1'000'000'000, Deadline(Clock::now(), 1e9, 0));
        const Result n = search.count(pattern_of(text, peptide), small, b2);
        EXPECT_EQ(Extension::NOT_ADMITTED, n.anchors->extension);
        EXPECT_TRUE(none.prepared.empty());

        // stopped in prepare: no anchor opened, the extension stopped before its first
        Recording stopper;
        stopper.stop_in_prepare = true;
        Request stopping = rq;
        stopping.support = &stopper;
        Budget b3(1'000'000'000, Deadline(Clock::now(), 1e9, 0));
        const Result s = search.count(pattern_of(text, peptide), stopping, b3);
        ASSERT_TRUE(s.stop);
        EXPECT_EQ(StopPhase::EXTENSION, s.stop->phase);
        EXPECT_EQ(StopReason::EXTERNAL, s.stop->reason);
        EXPECT_TRUE(stopper.opened.empty());
        EXPECT_EQ(Extension::STOPPED, s.anchors->extension);
        EXPECT_EQ(Relation::AT_LEAST, s.anchors->paths.relation);
        EXPECT_EQ(0u, s.anchors->paths.value);
        EXPECT_EQ(Relation::AT_LEAST, s.anchors->supported.relation);
        EXPECT_EQ(0u, s.anchors->supported.value);
        EXPECT_EQ(0u, s.work.extension_anchors);
    }
    EXPECT_LT(20u, with_anchors);
}

// the oracle: sdust over the whole text, as the seeder calls it (T = 20, W = 64)
bool sdust_flags(const std::string &text) {
    int n = 0;
    uint64_t *intervals = sdust(0, reinterpret_cast<const uint8_t*>(text.data()),
                                static_cast<int>(text.size()), 20, 64, &n);
    std::free(intervals);
    return n > 0;
}

bool noted_low_complexity(const Result &r) {
    return std::find(r.notes.begin(), r.notes.end(), std::string(kNoteLowComplexity))
            != r.notes.end();
}

// a path sink that stops the extension at the first path it is handed, as a selecting sink
// does at its threshold (|threshold|: stop_at_threshold) or at a budget of its own
class StoppingSink : public PathSink {
  public:
    explicit StoppingSink(bool threshold) : threshold_(threshold) {}
    bool accept(const PathView &, const SupportTracker *) override {
        ++accepted;
        return false;
    }
    const char* stop_reason() const override {
        return accepted ? (threshold_ ? "max_paths" : "max_memory") : nullptr;
    }
    bool stopped_at_threshold() const override { return accepted && threshold_; }
    uint64_t accepted = 0;

  private:
    bool threshold_;
};

TEST(PatternSupport, LowComplexityNoteBesideThresholdStopsOnly) {
    // the note is diagnosed beside a threshold stop (stop_at_threshold: max_contexts,
    // max_anchors, max_paths, a path sink's threshold), never beside a budget stop (max_steps,
    // time, a tracker's or a sink's own)
    const std::string low = "GCCGCCGCCGCCGCCGCCGC";
    const std::string plain = "ACGTTGCATGACCTAGGTCA";
    ASSERT_TRUE(sdust_flags(low));
    ASSERT_FALSE(sdust_flags(plain));
    std::mt19937 rng(3);
    std::vector<Record> records;
    for (int i = 0; i < 4; ++i) {
        records.push_back(Record { "c", "r" + std::to_string(i),
                                   random_bases(rng, 20) + low + low + random_bases(rng, 20)
                                   + plain + random_bases(rng, 5) + plain });
    }
    // k = 15: the 20-base patterns are long (anchors, paths); k = 31: they have contexts
    const Index idx = build_rd(15, records);
    const Index wide = build_rd(31, records);
    PatternSearch search(idx.anno->get_graph());
    PatternSearch wide_search(wide.anno->get_graph());
    auto count = [&](const std::string &text, Request rq, uint64_t max_steps = 1'000'000'000,
                     const PatternSearch *engine = nullptr) {
        rq.min_information_bits = 0;
        Budget budget(max_steps, Deadline(Clock::now(), 1e9, 0));
        return (engine ? *engine : search).count(Pattern::parse(PatternKind::DNA, text), rq,
                                                 budget);
    };
    for (const std::string &text : { low, plain }) {
        SCOPED_TRACE(text);
        const bool flagged = sdust_flags(text);
        // completed
        Request full;
        full.max_contexts = full.max_anchors = full.max_paths = 1'000'000;
        full.extend_paths = true;
        EXPECT_EQ(flagged, noted_low_complexity(count(text, full)));
        // max_anchors (L > k), with stop_at_threshold
        Request anchors = full;
        anchors.max_anchors = 0;
        anchors.stop_at_threshold = true;
        const Result a = count(text, anchors);
        ASSERT_TRUE(a.stop);
        EXPECT_EQ(StopReason::MAX_ANCHORS, a.stop->reason);
        EXPECT_EQ(flagged, noted_low_complexity(a));
        // max_paths in the extension
        Request paths = full;
        paths.max_paths = 0;
        paths.stop_at_threshold = true;
        const Result p = count(text, paths);
        ASSERT_TRUE(p.stop);
        EXPECT_EQ(StopReason::MAX_PATHS, p.stop->reason);
        EXPECT_EQ(flagged, noted_low_complexity(p));
        // max_contexts (L <= k)
        Request contexts = full;
        contexts.max_contexts = 0;
        contexts.stop_at_threshold = true;
        const Result c = count(text, contexts, 1'000'000'000, &wide_search);
        ASSERT_TRUE(c.contexts);
        ASSERT_TRUE(c.stop);
        EXPECT_EQ(StopReason::MAX_CONTEXTS, c.stop->reason);
        EXPECT_EQ(flagged, noted_low_complexity(c));
        EXPECT_EQ(flagged, noted_low_complexity(count(text, full, 1'000'000'000,
                                                      &wide_search)));
        // a budget stop: max_steps
        const Result s = count(text, full, 1);
        ASSERT_TRUE(s.stop);
        EXPECT_EQ(StopReason::MAX_STEPS, s.stop->reason);
        EXPECT_FALSE(noted_low_complexity(s));
        // a tracker's own stop
        Recording stopper;
        stopper.stop_in_prepare = true;
        Request tracked = full;
        tracked.support = &stopper;
        const Result t = count(text, tracked);
        ASSERT_TRUE(t.stop);
        EXPECT_EQ(StopReason::EXTERNAL, t.stop->reason);
        EXPECT_FALSE(noted_low_complexity(t));
        // a path sink's stop: its threshold, or a budget of its own
        for (const bool threshold : { true, false }) {
            SCOPED_TRACE(threshold ? "sink threshold" : "sink budget");
            StoppingSink sink(threshold);
            Request sunk = full;
            sunk.sink = &sink;
            const Result r = count(text, sunk);
            ASSERT_EQ(1u, sink.accepted);
            ASSERT_TRUE(r.stop);
            EXPECT_EQ(StopReason::EXTERNAL, r.stop->reason);
            EXPECT_EQ(flagged && threshold, noted_low_complexity(r));
        }
    }
}

} // namespace
