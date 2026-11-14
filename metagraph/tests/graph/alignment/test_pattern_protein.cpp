/**
 * The pattern search on protein graphs (docs/DESIGN-pattern-protein.md, E16): the residue
 * language (Pattern::parse on the Protein build), its bits, one strand, and the engine's counts
 * and lists on small protein graphs against a brute-force oracle over the graph's k-mers.
 *
 * The oracle reads a pattern through this file's own table (oracle_sets: the letters of each
 * position, classes and negations included), never through Pattern::positions(); the parser's
 * sets are checked against it position by position (ParseMatchesTheOracleTable). It reads the
 * graph's k-mers through call_kmers, and holds a stored ambiguity letter (B, J, Z, X) to be
 * matched by nothing, at no mismatch level.
 *
 * Alphabet: the Protein build only (_PROTEIN_GRAPH). The DNA builds' suites are the other
 * test_pattern_*.cpp files.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "../../test_helpers.hpp"
#include "../all/test_dbg_helpers.hpp"
#include "pattern_test_support.hpp"

#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace {

#if _PROTEIN_GRAPH

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::pattern;
using mtg::test::build_graph;
using mtg::test::unbounded_deadline;
using mtg::test::without_mask;

typedef DeBruijnGraph::node_index node_index;

constexpr uint64_t kManySteps = 1'000'000'000;

// the residues read as themselves and the ambiguity letters, written out here (the oracle's)
const std::string kLiterals = "ACDEFGHIKLMNPQRSTVWYUO";
const std::string kAmbiguous = "BJZX";

std::string letter_set(char c) {
    switch (c) {
        case 'X': return kLiterals;
        case 'B': return "DN";
        case 'Z': return "EQ";
        case 'J': return "IL";
        default: return kLiterals.find(c) != std::string::npos ? std::string(1, c) : "";
    }
}

// the oracle's reading of a pattern text: per position the residues it admits (sorted)
std::vector<std::string> oracle_sets(const std::string &text) {
    std::vector<std::string> sets;
    for (size_t i = 0; i < text.size(); ++i) {
        const char c = std::toupper(static_cast<unsigned char>(text[i]));
        std::string set;
        if (c == '[') {
            const bool negated = text[i + 1] == '^';
            size_t j = i + 1 + negated;
            std::string inside;
            for (; text[j] != ']'; ++j) {
                inside += letter_set(std::toupper(static_cast<unsigned char>(text[j])));
            }
            for (char r : kLiterals) {
                if ((inside.find(r) != std::string::npos) != negated)
                    set.push_back(r);
            }
            i = j;
        } else {
            set = letter_set(c);
        }
        std::sort(set.begin(), set.end());
        set.erase(std::unique(set.begin(), set.end()), set.end());
        sets.push_back(set);
    }
    return sets;
}

std::string sets_of(const Pattern &pattern, size_t i) {
    std::string set;
    for (char c = 'A'; c <= 'Z'; ++c) {
        if (pattern.positions()[i] & (BaseSet(1) << (c - 'A')))
            set.push_back(c);
    }
    return set;
}

Budget unbounded_budget(uint64_t max_steps = kManySteps) {
    return Budget(max_steps, unbounded_deadline());
}

Request request_of(Scope scope, uint8_t max_mismatches = 0) {
    Request r;
    r.scope = scope;
    r.min_information_bits = 0;
    r.max_mismatches = max_mismatches;
    r.max_contexts = 1'000'000;
    r.max_anchors = 1'000'000;
    r.max_paths = 1'000'000;
    return r;
}

// one context as the oracle keys it: (k-mer, offset, substitutions)
typedef std::tuple<std::string, uint32_t, uint32_t> Key;

/**
 * The contexts of a pattern of L <= k in the graph's k-mers: every (k-mer, offset) whose
 * instance holds at most |m| substitutions, a position matching when the k-mer's residue is in
 * its set and substituted when it is another residue read as itself; a stored ambiguity letter
 * at a pattern position rules the offset out.
 */
std::set<Key> oracle_contexts(const DeBruijnGraph &graph, const std::vector<std::string> &sets,
                              Scope scope, uint32_t m) {
    const size_t k = graph.get_k();
    const size_t L = sets.size();
    std::set<Key> found;
    graph.call_kmers([&](node_index, const std::string &kmer) {
        const size_t first = scope == Scope::SUFFIX ? k - L : 0;
        for (size_t p = first; p + L <= k; ++p) {
            uint32_t subs = 0;
            bool possible = true;
            for (size_t i = 0; i < L && possible; ++i) {
                const char c = kmer[p + i];
                if (kLiterals.find(c) == std::string::npos) {
                    possible = false;
                } else if (sets[i].find(c) == std::string::npos) {
                    ++subs;
                }
            }
            if (possible && subs <= m)
                found.emplace(kmer, static_cast<uint32_t>(p), subs);
        }
    });
    return found;
}

// the walks of L > k residues spelling the pattern: from every k-mer instantiating its first k
// positions, along outgoing edges whose residue is in the next position's set
uint64_t oracle_paths(const DeBruijnGraph &graph, const std::vector<std::string> &sets,
                      uint64_t *anchors) {
    const size_t k = graph.get_k();
    uint64_t paths = 0;
    *anchors = 0;
    std::function<void(node_index, size_t)> walk = [&](node_index node, size_t position) {
        if (position == sets.size()) {
            ++paths;
            return;
        }
        graph.call_outgoing_kmers(node, [&](node_index next, char c) {
            if (sets[position].find(c) != std::string::npos)
                walk(next, position + 1);
        });
    };
    graph.call_kmers([&](node_index node, const std::string &kmer) {
        for (size_t i = 0; i < k; ++i) {
            if (sets[i].find(kmer[i]) == std::string::npos)
                return;
        }
        ++*anchors;
        walk(node, k);
    });
    return paths;
}

// random proteins over the 22 residues, a few stored ambiguity letters among them, with the
// motifs planted (each a few times, some with a residue changed)
std::vector<std::string> random_records(std::mt19937 &rng, size_t n, size_t length,
                                        const std::vector<std::string> &motifs) {
    std::uniform_int_distribution<size_t> pick(0, kLiterals.size() - 1);
    std::uniform_int_distribution<size_t> percent(0, 99);
    std::vector<std::string> records;
    for (size_t r = 0; r < n; ++r) {
        std::string s;
        for (size_t i = 0; i < length; ++i) {
            s.push_back(percent(rng) < 4 ? kAmbiguous[pick(rng) % 4] : kLiterals[pick(rng)]);
        }
        for (const std::string &motif : motifs) {
            std::string planted = motif;
            if (percent(rng) < 50)
                planted[pick(rng) % planted.size()] = kLiterals[pick(rng)];
            const size_t at = pick(rng) * (length - planted.size()) / kLiterals.size();
            s.replace(at, planted.size(), planted);
        }
        records.push_back(s);
    }
    return records;
}


// ---------------------------------------------------------------- the language

TEST(PatternProtein, ParseMatchesTheOracleTable) {
    for (const std::string text : { "MKVLAAG", "mkvlaag", "UOX", "BZJ", "[STA]", "[^P]",
                                    "A[ST]X[^PG]B", "[x]", "[BZ]", "[^bzj]K", "[UO]W",
                                    "[A][C]", "[^ACDEFGHIKLMNPQRSTVWYU]" }) {
        SCOPED_TRACE(text);
        const Pattern pattern = Pattern::parse(PatternKind::PROTEIN, text);
        const std::vector<std::string> sets = oracle_sets(text);
        ASSERT_EQ(sets.size(), pattern.length());
        for (size_t i = 0; i < sets.size(); ++i) {
            EXPECT_EQ(sets[i], sets_of(pattern, i)) << i;
        }
        EXPECT_FALSE(pattern.codon_automaton());
        EXPECT_FALSE(pattern.is_palindromic());
        EXPECT_TRUE(pattern.has_instances());
        EXPECT_EQ(PatternKind::PROTEIN, pattern.kind());
    }
    // the text in upper case, classes kept; position_text names each position as written
    const Pattern p = Pattern::parse(PatternKind::PROTEIN, "a[st]x[^pg]b");
    EXPECT_EQ("A[ST]X[^PG]B", p.text());
    ASSERT_EQ(5u, p.length());
    EXPECT_EQ("A", p.position_text(0));
    EXPECT_EQ("[ST]", p.position_text(1));
    EXPECT_EQ("X", p.position_text(2));
    EXPECT_EQ("[^PG]", p.position_text(3));
    EXPECT_EQ("B", p.position_text(4));
    // the stored ambiguity letters are in no set; X is the 22 residues
    EXPECT_EQ(kAllResidues, Pattern::parse(PatternKind::PROTEIN, "X").positions()[0]);
    for (const char *text : { "X", "[^P]", "[BJZX]" }) {
        for (char c : kAmbiguous) {
            EXPECT_FALSE(Pattern::parse(PatternKind::PROTEIN, text).positions()[0]
                            & residue_bit(c)) << text << " " << c;
        }
    }
    // a one-residue class is exact
    EXPECT_TRUE(Pattern::parse(PatternKind::PROTEIN, "[S]TA").is_exact());
    EXPECT_FALSE(Pattern::parse(PatternKind::PROTEIN, "[ST]TA").is_exact());
    // the genetic code is not read
    EXPECT_EQ(Pattern::parse(PatternKind::PROTEIN, "MW", GeneticCode::get(11)).positions(),
              Pattern::parse(PatternKind::PROTEIN, "MW").positions());
}

TEST(PatternProtein, ParseRefusals) {
    // (text, the offending character's 0-based place)
    const std::vector<std::pair<std::string, size_t>> refused = {
        { "MK*", 2 }, { "MK1", 2 }, { "M-K", 1 }, { "M K", 1 }, { "M.K", 1 }, { "x(2,4)", 1 },
        { "{P}A", 0 }, { "<MK", 0 }, { "MK>", 2 }, { "[]", 0 }, { "[^]", 0 }, { "A[^X]", 1 },
        { "[ST", 0 }, { "A[S[T]]", 3 }, { "AS]", 2 }, { "A^S", 1 }, { "[S^T]", 2 },
        { "[S*]", 2 }, { "[S1]", 2 },
    };
    for (const auto &[text, at] : refused) {
        SCOPED_TRACE(text);
        try {
            Pattern::parse(PatternKind::PROTEIN, text);
            ADD_FAILURE() << "parsed";
        } catch (const PatternError &e) {
            EXPECT_EQ("bad_alphabet", e.code());
            EXPECT_NE(std::string::npos,
                      std::string(e.what()).find("at position " + std::to_string(at) + " "))
                << e.what();
        }
    }
    try {
        Pattern::parse(PatternKind::PROTEIN, "");
        ADD_FAILURE() << "parsed an empty text";
    } catch (const PatternError &e) {
        EXPECT_EQ("bad_alphabet", e.code());
    }
    // the nucleotide kinds: a protein graph holds no nucleotides
    for (PatternKind kind : { PatternKind::DNA, PatternKind::IUPAC }) {
        try {
            Pattern::parse(kind, "ACGT");
            ADD_FAILURE() << "parsed " << to_string(kind);
        } catch (const PatternError &e) {
            EXPECT_EQ("bad_alphabet", e.code());
            EXPECT_NE(std::string::npos, std::string(e.what()).find("protein graph"));
        }
    }
}

TEST(PatternProtein, InformationBits) {
    const double residue = std::log2(22.0);
    EXPECT_DOUBLE_EQ(residue, Pattern::parse(PatternKind::PROTEIN, "W").information_bits());
    EXPECT_DOUBLE_EQ(0, Pattern::parse(PatternKind::PROTEIN, "X").information_bits());
    EXPECT_DOUBLE_EQ(std::log2(11.0), Pattern::parse(PatternKind::PROTEIN, "B").information_bits());
    EXPECT_DOUBLE_EQ(std::log2(22.0 / 3),
                     Pattern::parse(PatternKind::PROTEIN, "[STA]").information_bits());
    EXPECT_DOUBLE_EQ(std::log2(22.0 / 21),
                     Pattern::parse(PatternKind::PROTEIN, "[^P]").information_bits());
    const Pattern p = Pattern::parse(PatternKind::PROTEIN, "MKV[ST]X");
    EXPECT_NEAR(3 * residue + std::log2(11.0), p.information_bits(), 1e-9);
    EXPECT_NEAR(2 * residue, p.information_bits(1, 3), 1e-9);
    // within m substitutions: the ball over the 22 residues (21 substitutes for a residue, 20
    // for a pair, none for X)
    const Pattern w = Pattern::parse(PatternKind::PROTEIN, "WWX");
    EXPECT_NEAR(3 * residue - std::log2(22.0 * (1 + 2 * 21)), w.mismatch_information_bits(1),
                1e-9);
    EXPECT_NEAR(3 * residue - std::log2(22.0 * (1 + 2 * 21 + 21 * 21)),
                w.mismatch_information_bits(2), 1e-9);
    EXPECT_DOUBLE_EQ(w.information_bits(), w.mismatch_information_bits(0));
}

TEST(PatternProtein, OneStrand) {
    const Pattern p = Pattern::parse(PatternKind::PROTEIN, "MKVK");
    EXPECT_FALSE(p.is_palindromic());
    EXPECT_THROW(p.reverse_complement(), std::logic_error);
    // a context's substitutions are residues, read forward
    const auto subs = substitutions(Pattern::parse(PatternKind::PROTEIN, "M[ST]VK"),
                                    Orientation::FORWARD, "MAVR");
    ASSERT_EQ(2u, subs.size());
    EXPECT_EQ(1u, subs[0].position);
    EXPECT_EQ('[', subs[0].expected);
    EXPECT_EQ('A', subs[0].found);
    EXPECT_EQ(3u, subs[1].position);
    EXPECT_EQ('K', subs[1].expected);
    EXPECT_EQ('R', subs[1].found);
    EXPECT_THROW(substitutions(p, Orientation::REVERSE, "MKVK"), std::invalid_argument);
    // a stored ambiguity letter is no residue of an instance
    EXPECT_THROW(substitutions(p, Orientation::FORWARD, "MKXK"), std::invalid_argument);
    EXPECT_THROW(residue_substitutions(p, Orientation::FORWARD, "MKVK"),
                 std::invalid_argument);
}


// ---------------------------------------------------------------- the engine

TEST(PatternProtein, SupportOfAProteinGraph) {
    auto graph = build_graph<DBGSuccinct>(5, { "MKVLAAGIVGLLLAQSEWHT" });
    const GraphSupport support = PatternSearch::support(*graph);
    EXPECT_TRUE(support.supported) << support.reason;
    EXPECT_EQ("$ABCDEFGHIJKLMNOPQRSTUVWYZX", support.alphabet);
    EXPECT_EQ(GraphMode::BASIC, support.mode);
    EXPECT_TRUE(support.strand_stated);
    EXPECT_TRUE(support.mask_present);
    EXPECT_EQ(5u, support.k);
}

// count() and enumerate() against the oracle: every scope, 0 to 3 substitutions, with and
// without the dummy-edge mask
TEST(PatternProtein, ContextsMatchTheOracle) {
    std::mt19937 rng(20261107);
    const std::vector<std::string> motifs = { "DNPQE", "STAMK", "WHEDN", "GGSTA" };
    const std::vector<std::string> records = random_records(rng, 12, 60, motifs);
    const size_t k = 7;
    const std::vector<std::string> patterns = {
        "DNPQE", "dnpqe", "STAMK", "[ST]TAM", "S[^P]AMK", "XNPQE", "BNPQZ", "W[HK]EDN", "J",
        "GG", "[GS]G[ST]", "X", "XX", "WHEDNPQ", "[^W]HEDN", "U", "O[UO]", "MKV", "AMKXX",
    };
    for (bool masked : { true, false }) {
        auto graph = build_graph<DBGSuccinct>(k, records);
        // the oracle reads the masked graph's k-mers: the true ones
        auto reference = build_graph<DBGSuccinct>(k, records);
        if (!masked)
            graph = without_mask(graph);
        PatternSearch search(*graph);
        for (const std::string &text : patterns) {
            const Pattern pattern = Pattern::parse(PatternKind::PROTEIN, text);
            const std::vector<std::string> sets = oracle_sets(text);
            for (Scope scope : { Scope::ANY_OFFSET, Scope::SUFFIX }) {
                for (uint8_t m = 0; m <= 3; ++m) {
                    if (m > pattern.length())
                        continue;
                    SCOPED_TRACE(text + (scope == Scope::SUFFIX ? " suffix" : " any_offset")
                                 + " m=" + std::to_string(m) + (masked ? "" : " unmasked"));
                    const std::set<Key> expected = oracle_contexts(*reference, sets, scope, m);
                    Request request = request_of(scope, m);
                    // without the mask: every candidate checked, the counts exact
                    request.max_checked_entries = masked ? 0 : 1'000'000;
                    Budget budget = unbounded_budget();
                    std::set<Key> released;
                    const Result result = search.enumerate(pattern, request, budget,
                            [&](const Context &c) {
                                EXPECT_EQ(Orientation::FORWARD, c.orientation);
                                released.emplace(graph->get_node_sequence(c.node), c.offset,
                                                 c.mismatches);
                            });
                    ASSERT_FALSE(result.refusal) << result.refusal->message;
                    ASSERT_TRUE(result.contexts);
                    EXPECT_EQ(Relation::EXACT, result.contexts->total.relation);
                    EXPECT_EQ(expected.size(), result.contexts->total.value);
                    EXPECT_EQ(expected, released);
                    ASSERT_EQ(std::vector<Orientation> { Orientation::FORWARD },
                              result.searched);
                    EXPECT_FALSE(result.palindromic);
                    ASSERT_EQ(1u, result.contexts->by_orientation.size());
                    EXPECT_EQ(expected.size(), result.contexts->by_orientation.at(
                            Orientation::FORWARD).value);
                    if (m) {
                        ASSERT_EQ(m + 1u, result.contexts->by_mismatches.size());
                        for (uint32_t t = 0; t <= m; ++t) {
                            const uint64_t at_level = std::count_if(
                                    expected.begin(), expected.end(),
                                    [&](const Key &key) { return std::get<2>(key) == t; });
                            EXPECT_EQ(at_level, result.contexts->by_mismatches[t].value) << t;
                        }
                    }
                }
            }
        }
    }
}

// without the mask and without the check of a few candidates: bounds around the true count
TEST(PatternProtein, UnmaskedBoundsHoldTheCount) {
    std::mt19937 rng(7);
    const std::vector<std::string> records = random_records(rng, 6, 40, { "KDNPQEL" });
    const size_t k = 6;
    auto reference = build_graph<DBGSuccinct>(k, records);
    auto graph = without_mask(build_graph<DBGSuccinct>(k, records));
    PatternSearch search(*graph);
    for (const std::string text : { "DNPQ", "X", "[DE]", "KDN" }) {
        SCOPED_TRACE(text);
        const Pattern pattern = Pattern::parse(PatternKind::PROTEIN, text);
        const uint64_t truth = oracle_contexts(*reference, oracle_sets(text),
                                               Scope::ANY_OFFSET, 0).size();
        Request request = request_of(Scope::ANY_OFFSET);
        request.mode = Mode::COUNT;
        Budget budget = unbounded_budget();
        const Result result = search.count(pattern, request, budget);
        ASSERT_TRUE(result.contexts);
        const Count &c = result.contexts->total;
        if (c.relation == Relation::EXACT) {
            EXPECT_EQ(truth, c.value);
        } else {
            ASSERT_EQ(Relation::BOUNDS, c.relation);
            EXPECT_LE(c.lower, truth);
            EXPECT_GE(c.upper, truth);
        }
    }
}

// patterns longer than k: the anchors (P[0, k)) and, extended, the walks spelling the pattern
TEST(PatternProtein, LongPatternsAnchorsAndPaths) {
    const std::vector<std::string> records = {
        "MKVLAAGIVGLLLAQSEWHTRDNPQELKYSTAMKV",
        "MKVLAAGPVGRWHEDNQSTAEEKLIQRSTW",
        "QQMKVLAAGIVGLLRAQSEWHT",
        "GGSTARDNPQEIKYSTAMUOPWQ",
    };
    const size_t k = 6;
    auto graph = build_graph<DBGSuccinct>(k, records);
    PatternSearch search(*graph);
    for (const std::string text : { "MKVLAAGIVGLL", "MKVLAAG[IP]VG", "MKVLAAGXVG", "AGIVGLL[LR]AQSEW",
                                    "RDNPQE[LI]KYSTAM", "DNPQEJKY", "MKVLAAGIVGLLW" }) {
        SCOPED_TRACE(text);
        const Pattern pattern = Pattern::parse(PatternKind::PROTEIN, text);
        ASSERT_GT(pattern.length(), k);
        uint64_t anchors = 0;
        const uint64_t paths = oracle_paths(*graph, oracle_sets(text), &anchors);
        Request request = request_of(Scope::ANY_OFFSET);
        request.extend_paths = true;
        Budget budget = unbounded_budget();
        std::vector<Context> released;
        const Result result = search.enumerate(pattern, request, budget,
                                               [&](const Context &c) { released.push_back(c); });
        ASSERT_FALSE(result.refusal) << result.refusal->message;
        EXPECT_EQ(Scope::LONG, result.scope);
        ASSERT_TRUE(result.anchors);
        EXPECT_EQ(Relation::EXACT, result.anchors->total.relation);
        EXPECT_EQ(anchors, result.anchors->total.value);
        EXPECT_EQ(Relation::EXACT, result.anchors->paths.relation);
        EXPECT_EQ(paths, result.anchors->paths.value);
        EXPECT_EQ(paths, released.size());
        for (const Context &c : released) {
            EXPECT_EQ(Orientation::FORWARD, c.orientation);
            ASSERT_EQ(pattern.length(), c.sequence.size());
            const std::vector<std::string> sets = oracle_sets(text);
            for (size_t i = 0; i < sets.size(); ++i) {
                EXPECT_NE(std::string::npos, sets[i].find(c.sequence[i])) << c.sequence;
            }
            EXPECT_EQ(pattern.length() - k + 1, c.path.size());
        }
    }
}

// the request's strands: both and forward search the one strand; reverse is refused
TEST(PatternProtein, StrandsOnAProteinGraph) {
    auto graph = build_graph<DBGSuccinct>(5, { "MKVLAAGIVGLLLAQSEWHT" });
    PatternSearch search(*graph);
    const Pattern pattern = Pattern::parse(PatternKind::PROTEIN, "AAGIV");
    for (Strands strands : { Strands::BOTH, Strands::FORWARD }) {
        Request request = request_of(Scope::ANY_OFFSET);
        request.strands = strands;
        Budget budget = unbounded_budget();
        const Result result = search.count(pattern, request, budget);
        ASSERT_TRUE(result.contexts);
        EXPECT_EQ(1u, result.contexts->total.value);
        EXPECT_EQ(std::vector<Orientation> { Orientation::FORWARD }, result.searched);
    }
    Request request = request_of(Scope::ANY_OFFSET);
    request.strands = Strands::REVERSE;
    Budget budget = unbounded_budget();
    EXPECT_THROW(search.count(pattern, request, budget), std::invalid_argument);
}

// a leading run of X branches 22 ways per position (no shortcut skips it: a stored ambiguity
// letter is no residue of X): max_steps stops its discovery, the count a lower bound of the truth
TEST(PatternProtein, MaxStepsStopsAnXRun) {
    std::mt19937 rng(15);
    const std::vector<std::string> records = random_records(rng, 20, 80, { "MKWWWWWW" });
    const size_t k = 15;
    auto graph = build_graph<DBGSuccinct>(k, records);
    PatternSearch search(*graph);
    const std::string text = "XXXXXXXXXWWWWWW";
    const Pattern pattern = Pattern::parse(PatternKind::PROTEIN, text);
    const uint64_t truth = oracle_contexts(*graph, oracle_sets(text), Scope::ANY_OFFSET, 0).size();
    {
        Request request = request_of(Scope::ANY_OFFSET);
        request.mode = Mode::COUNT;
        Budget budget = unbounded_budget(20);
        const Result result = search.count(pattern, request, budget);
        ASSERT_TRUE(result.stop);
        EXPECT_EQ(StopPhase::DISCOVERY, result.stop->phase);
        EXPECT_EQ(StopReason::MAX_STEPS, result.stop->reason);
        ASSERT_TRUE(result.contexts);
        EXPECT_EQ(Relation::AT_LEAST, result.contexts->total.relation);
        EXPECT_LE(result.contexts->total.value, truth);
        EXPECT_LE(result.work.steps, 20u);
    }
    {
        // and with the budget it needs, the oracle's count
        Request request = request_of(Scope::ANY_OFFSET);
        request.mode = Mode::COUNT;
        Budget budget = unbounded_budget();
        const Result result = search.count(pattern, request, budget);
        EXPECT_FALSE(result.stop);
        ASSERT_TRUE(result.contexts);
        EXPECT_EQ(Relation::EXACT, result.contexts->total.relation);
        EXPECT_EQ(truth, result.contexts->total.value);
        EXPECT_GT(result.work.steps, 20u);
    }
}

// the information floor gates residues as it gates bases; no low-complexity note on residues
TEST(PatternProtein, FloorAndNotes) {
    auto graph = build_graph<DBGSuccinct>(7, { "MKVLAAGIVGLLLAQSEWHTQQQQQQQQQQQQQQQQ" });
    PatternSearch search(*graph);
    Request request = request_of(Scope::ANY_OFFSET);
    request.min_information_bits = 24;
    {
        // five residues: 22.3 bits, below the floor
        Budget budget = unbounded_budget();
        const Result result = search.count(Pattern::parse(PatternKind::PROTEIN, "AAGIV"),
                                           request, budget);
        ASSERT_TRUE(result.refusal);
        EXPECT_EQ("information_below_floor", result.refusal->code);
    }
    {
        // six: 26.8 bits; a run of Q is exact but carries no low-complexity note here
        Budget budget = unbounded_budget();
        const Result result = search.count(Pattern::parse(PatternKind::PROTEIN, "QQQQQQ"),
                                           request, budget);
        ASSERT_FALSE(result.refusal);
        EXPECT_EQ(oracle_contexts(*graph, oracle_sets("QQQQQQ"), Scope::ANY_OFFSET, 0).size(),
                  result.contexts->total.value);
        EXPECT_LT(0u, result.contexts->total.value);
        EXPECT_TRUE(std::find(result.notes.begin(), result.notes.end(), kNoteLowComplexity)
                        == result.notes.end());
    }
}


// ---------------------------------------------------------------- more than four children

// the walks of L > k residues spelling the pattern (oracle_paths), as their sequences
std::multiset<std::string> oracle_walks(const DeBruijnGraph &graph,
                                        const std::vector<std::string> &sets) {
    const size_t k = graph.get_k();
    std::multiset<std::string> walks;
    std::function<void(node_index, const std::string&)> walk
            = [&](node_index node, const std::string &spelled) {
        if (spelled.size() == sets.size()) {
            walks.insert(spelled);
            return;
        }
        graph.call_outgoing_kmers(node, [&](node_index next, char c) {
            if (sets[spelled.size()].find(c) != std::string::npos)
                walk(next, spelled + c);
        });
    };
    graph.call_kmers([&](node_index node, const std::string &kmer) {
        for (size_t i = 0; i < k; ++i) {
            if (sets[i].find(kmer[i]) == std::string::npos)
                return;
        }
        walk(node, kmer);
    });
    return walks;
}

// the k-mer |prefix| followed by each of |residues| (one record each)
std::vector<std::string> fan_out(const std::string &prefix, const std::string &residues) {
    std::vector<std::string> records;
    for (char r : residues) {
        records.push_back(prefix + r);
    }
    return records;
}

// a node with more than four allowed children (the review of 2026-11-14, [P1]: the extension
// kept four per level, a Release build dropping the others and their walks from an exact count):
// every allowed child is extended, at one level and at two in a row, up to the 22 residues a
// position admits; the stored ambiguity letters among a node's children are admitted by none
TEST(PatternProtein, ExtensionKeepsEveryChild) {
    const size_t k = 6;
    // CDEFGK: the review's five children; MNPQRS: all 26 symbols of the alphabet after it;
    // WYTVLI: six children, each with six more
    std::vector<std::string> records = fan_out("CDEFGK", "ACDFH");
    for (const std::string &r : fan_out("MNPQRS", kLiterals + kAmbiguous)) {
        records.push_back(r);
    }
    for (char a : std::string("ACDFHK")) {
        for (const std::string &r : fan_out(std::string("WYTVLI") + a, "EGMNPQ")) {
            records.push_back(r);
        }
    }
    auto graph = build_graph<DBGSuccinct>(k, records);
    PatternSearch search(*graph);
    // (pattern, its walks)
    const std::vector<std::pair<std::string, uint64_t>> cases = {
        { "CDEFGKX", 5 }, { "CDEFGK[^A]", 4 }, { "CDEFGK[CDFHW]", 4 }, { "CDEFGK[FH]", 2 },
        { "MNPQRSX", 22 }, { "MNPQRS[^W]", 21 }, { "MNPQRS[UOY]", 3 }, { "MNPQRSB", 2 },
        { "WYTVLIXX", 36 }, { "WYTVLIX[^E]", 30 }, { "XYTVLIXX", 36 }, { "WYTVLI[HK]X", 12 },
    };
    for (const auto &[text, expected] : cases) {
        SCOPED_TRACE(text);
        const Pattern pattern = Pattern::parse(PatternKind::PROTEIN, text);
        ASSERT_GT(pattern.length(), k);
        const std::multiset<std::string> walks = oracle_walks(*graph, oracle_sets(text));
        ASSERT_EQ(expected, walks.size());
        Request request = request_of(Scope::ANY_OFFSET);
        request.extend_paths = true;
        {
            Budget budget = unbounded_budget();
            std::multiset<std::string> released;
            const Result result = search.enumerate(pattern, request, budget,
                    [&](const Context &c) {
                        EXPECT_EQ(pattern.length() - k + 1, c.path.size());
                        released.insert(c.sequence);
                    });
            ASSERT_FALSE(result.refusal) << result.refusal->message;
            EXPECT_FALSE(result.stop);
            ASSERT_TRUE(result.anchors);
            EXPECT_EQ(Relation::EXACT, result.anchors->paths.relation);
            EXPECT_EQ(expected, result.anchors->paths.value);
            EXPECT_EQ(walks, released);
        }
        {
            // count(): the same count
            request.mode = Mode::COUNT;
            Budget budget = unbounded_budget();
            const Result result = search.count(pattern, request, budget);
            ASSERT_FALSE(result.refusal) << result.refusal->message;
            ASSERT_TRUE(result.anchors);
            EXPECT_EQ(Relation::EXACT, result.anchors->paths.relation);
            EXPECT_EQ(expected, result.anchors->paths.value);
        }
    }
}

// the segment search's arms on a protein graph: the anchor's node has five (then 22) allowed
// children or parents, and only the alphabetically last one's walk goes on, so that it alone
// sets the extent's end (right arm) or start (left arm)
TEST(PatternProtein, SegmentArmsKeepEveryChild) {
    const size_t k = 6;
    struct Case {
        std::string what;
        std::vector<std::string> records;
        std::string text;
        WindowChoice choice;
        // the information floor (window 0 below it for the left arm: window 2's root then
        // extends left from position 2 to 0)
        double floor;
        uint32_t window;
        uint32_t from;
        uint32_t to;
    };
    std::vector<Case> cases;
    {
        // right: CDEFGK's children A C D F H, only H followed by W
        std::vector<std::string> records = fan_out("CDEFGK", "ACDF");
        records.push_back("CDEFGKHW");
        cases.push_back({ "right, five", records, "CDEFGKXW", WindowChoice::FIRST, 0, 0, 0, 8 });
        // all 26 symbols after MNPQRS, only Y followed by W
        records = fan_out("MNPQRS", kLiterals + kAmbiguous);
        records.push_back("MNPQRSYW");
        cases.push_back({ "right, 22", records, "MNPQRSXW", WindowChoice::FIRST, 0, 0, 0, 8 });
    }
    for (const auto &[what, parents] : { std::make_pair(std::string("left, five"),
                                                        std::string("ACDFH")),
                                         std::make_pair(std::string("left, 22"),
                                                        kLiterals + kAmbiguous) }) {
        // left: CDEFGK's parents, only the last allowed one (H; Y) preceded by W
        std::vector<std::string> records;
        for (char r : parents) {
            records.push_back(std::string(1, r) + "CDEFGK");
        }
        const char last = what == "left, five" ? 'H' : 'Y';
        records.push_back(std::string("W") + last + "CDEFGK");
        cases.push_back({ what, records, "WXCDEFGK", WindowChoice::ENDS, 24, 2, 0, 8 });
    }
    for (const Case &c : cases) {
        SCOPED_TRACE(c.what + " " + c.text);
        auto graph = build_graph<DBGSuccinct>(k, c.records);
        const Pattern pattern = Pattern::parse(PatternKind::PROTEIN, c.text);
        Request request = request_of(Scope::ANY_OFFSET);
        request.min_information_bits = c.floor;
        request.segments.enabled = true;
        request.segments.choice = c.choice;
        GraphSegmentStore store(static_cast<uint32_t>(pattern.length()),
                                std::numeric_limits<uint64_t>::max() / 2);
        request.segment_sink = &store;
        Budget budget = unbounded_budget();
        PatternSearch search(*graph);
        const Result result = search.count(pattern, request, budget);
        ASSERT_FALSE(result.refusal) << result.refusal->message;
        EXPECT_FALSE(result.stop);
        ASSERT_TRUE(result.segments);
        ASSERT_EQ(1u, result.segments->windows.count(Orientation::FORWARD));
        bool extended = false;
        for (const WindowCount &w : result.segments->windows.at(Orientation::FORWARD)) {
            if (w.offset == c.window) {
                EXPECT_EQ(WindowState::EXTENDED, w.state);
                EXPECT_EQ(1u, w.roots);
                extended = true;
            } else {
                EXPECT_EQ(WindowState::BELOW_FLOOR, w.state) << w.offset;
            }
        }
        EXPECT_TRUE(extended);
        EXPECT_EQ(1u, result.segments->reported);
        std::vector<std::tuple<uint32_t, uint32_t, uint64_t, uint32_t>> segments;
        store.for_each([&](const GraphSegment &s) {
            EXPECT_EQ(Orientation::FORWARD, s.orientation);
            segments.emplace_back(s.from, s.to, s.roots, s.window);
        });
        EXPECT_EQ((std::vector<std::tuple<uint32_t, uint32_t, uint64_t, uint32_t>> {
                      { c.from, c.to, 1, c.window } }), segments);
        // the memoised states released: the store holds its one segment's entry
        EXPECT_EQ(GraphSegmentStore::kGraphSegmentBytes * store.size(), store.held());
    }
}

#endif // _PROTEIN_GRAPH

} // namespace
