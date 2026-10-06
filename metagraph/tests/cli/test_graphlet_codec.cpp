#include "gtest/gtest.h"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <thread>
#include <unistd.h>

#include <json/json.h>

#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "cli/server_checks.hpp"
#include "cli/traverse.hpp"
#include "cli/traverse_attempts.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/seq_tools/reverse_complement.hpp"


/*
 * The MGT v1 codec against the shared golden vectors (DESIGN-traverse-graphlet.md §2.7,
 * the freeze gate): api/python/tests/data/traverse/codec_vectors.tsv, read here and by
 * the Python library's tests, byte-exactly on both sides. The file's header is the
 * contract: TAB-separated <tag> <input JSON> <expected JSON> [<note>]; for an encode tag
 * encode(input) == expected AND decode(expected) == input, for a decode tag a null
 * expectation means decoding fails and otherwise encode(decode(input)) == expected (a
 * valid non-canonical token dumps canonical). The COUNTS line says how many of each tag
 * a reader must have seen.
 */
namespace {

using namespace mtg::cli;
using namespace mtg::cli::mgt;

// run from the build directory, like every test reading ../tests/data
const char kVectors[] = "../api/python/tests/data/traverse/codec_vectors.tsv";

Json::Value parse_cell(const std::string &cell) {
    Json::CharReaderBuilder builder;
    builder["strictRoot"] = false;      // a cell may be a bare string or null
    builder["allowComments"] = false;
    builder["failIfExtra"] = true;
    builder["rejectDupKeys"] = true;
    std::unique_ptr<Json::CharReader> reader(builder.newCharReader());
    Json::Value v;
    std::string errs;
    if (!reader->parse(cell.data(), cell.data() + cell.size(), &v, &errs))
        throw std::runtime_error("bad JSON cell " + cell + ": " + errs);
    return v;
}

double from_bits(const std::string &hex) {
    const uint64_t bits = std::stoull(hex, nullptr, 16);
    double x;
    std::memcpy(&x, &bits, sizeof(x));
    return x;
}

uint64_t bits_of(double x) {
    uint64_t bits;
    std::memcpy(&bits, &x, sizeof(bits));
    return bits;
}

std::vector<uint64_t> ids_of(const Json::Value &v) {
    std::vector<uint64_t> out;
    for (const Json::Value &x : v) out.push_back(x.asUInt64());
    return out;
}

KValue kvalue_of(const Json::Value &v) {
    KValue k;
    const std::string type = v[0].asString();
    if (type == "i") {
        k.type = KValue::INTEGER;
        size_t used = 0;
        k.integer = std::stoull(v[1].asString(), &used, 10);
        EXPECT_EQ(v[1].asString().size(), used);
    } else if (type == "f") {
        k.type = KValue::FLOAT;
        k.real = from_bits(v[1].asString());
    } else if (type == "s") {
        k.type = KValue::STRING;
        k.string = v[1].asString();
    } else {
        EXPECT_EQ("u", type);
        k.type = KValue::UNLIMITED;
    }
    return k;
}

bool same_kvalue(const KValue &a, const KValue &b) {
    if (a.type != b.type)
        return false;
    switch (a.type) {
        case KValue::INTEGER: return a.integer == b.integer;
        // the -0 vectors decode "0" to +0
        case KValue::FLOAT: return a.real == b.real && (a.real != 0 || bits_of(b.real) == 0);
        case KValue::STRING: return a.string == b.string;
        case KValue::UNLIMITED: return true;
    }
    return false;
}

// encode(decode(token)) for a decode tag, or throws FormatError
std::string redump(const std::string &tag, const Json::Value &input) {
    if (tag == "float_decode") return encode_float(decode_float(input.asString()));
    if (tag == "ranges_decode") return encode_ranges(decode_ranges(input.asString()));
    if (tag == "setexpr_decode") {
        std::vector<uint64_t> base = ids_of(input[1]);
        const std::vector<uint64_t> *b = input[1].isNull() ? nullptr : &base;
        return encode_setexpr(decode_setexpr(input[0].asString(), b), b);
    }
    if (tag == "pct_decode") return pct_escape(pct_unescape(input.asString()));
    if (tag == "front_decode") {
        const std::string previous = input[0].asString();
        return front_code(previous, front_decode(previous, input[1].asString()));
    }
    if (tag == "kvalue_decode") return encode_kvalue(decode_kvalue(input.asString()));
    throw std::runtime_error("unknown tag " + tag);
}

} // namespace


TEST(GraphletCodec, GoldenVectors) {
    std::ifstream in(kVectors, std::ios::binary);
    ASSERT_TRUE(in.good()) << kVectors << " not found (run from the build directory): the "
                              "golden vectors are the freeze gate of MGT v1";
    std::map<std::string, size_t> seen, expected_counts;
    size_t failures = 0;
    auto report = [&](size_t line_no, const std::string &line, const std::string &what) {
        if (++failures <= 20)
            ADD_FAILURE() << "line " << line_no << ": " << what << "\n  " << line;
    };
    std::string line;
    size_t line_no = 0;
    while (std::getline(in, line)) {
        line_no++;
        if (line.rfind("#   float=", 0) == 0) {
            // the COUNTS line: tag=n ... total=n
            std::istringstream ss(line.substr(4));
            std::string item;
            while (ss >> item) {
                const size_t eq = item.find('=');
                expected_counts[item.substr(0, eq)] = std::stoull(item.substr(eq + 1));
            }
            continue;
        }
        if (line.empty() || line[0] == '#')
            continue;
        std::vector<std::string> cells;
        size_t pos = 0;
        while (true) {
            const size_t tab = line.find('\t', pos);
            cells.push_back(line.substr(pos, tab == std::string::npos ? std::string::npos : tab - pos));
            if (tab == std::string::npos)
                break;
            pos = tab + 1;
        }
        ASSERT_TRUE(cells.size() == 3 || cells.size() == 4) << "line " << line_no;
        const std::string &tag = cells[0];
        const Json::Value input = parse_cell(cells[1]);
        const Json::Value expected = parse_cell(cells[2]);
        seen[tag]++;
        try {
            if (tag.size() > 7 && tag.compare(tag.size() - 7, 7, "_decode") == 0) {
                if (expected.isNull()) {
                    try {
                        const std::string got = redump(tag, input);
                        report(line_no, line, "decoded to '" + got + "', must fail");
                    } catch (const FormatError &) {}
                    continue;
                }
                const std::string got = redump(tag, input);
                if (got != expected.asString())
                    report(line_no, line, "dumps as '" + got + "'");
                continue;
            }
            const std::string want = expected.asString();
            if (tag == "float") {
                const double x = from_bits(input.asString());
                const std::string got = encode_float(x);
                const double back = decode_float(want);
                if (got != want)
                    report(line_no, line, "encoded '" + got + "'");
                // the same double, except that -0 is written "0" and decodes to +0
                else if (x == 0 ? bits_of(back) != 0 : back != x)
                    report(line_no, line, "decodes to other bits");
            } else if (tag == "ranges") {
                const std::vector<uint64_t> ids = ids_of(input);
                const std::string got = encode_ranges(ids);
                if (got != want || decode_ranges(want) != ids)
                    report(line_no, line, "encoded '" + got + "'");
            } else if (tag == "setexpr") {
                const std::vector<uint64_t> ids = ids_of(input[0]), base = ids_of(input[1]);
                const std::vector<uint64_t> *b = input[1].isNull() ? nullptr : &base;
                const std::string got = encode_setexpr(ids, b);
                if (got != want || decode_setexpr(want, b) != ids)
                    report(line_no, line, "encoded '" + got + "'");
            } else if (tag == "pct") {
                const std::string got = pct_escape(input.asString());
                if (got != want || pct_unescape(want) != input.asString())
                    report(line_no, line, "encoded '" + got + "'");
            } else if (tag == "front") {
                const std::string previous = input[0].asString(), name = input[1].asString();
                const std::string got = front_code(previous, name);
                if (got != want || front_decode(previous, want) != name)
                    report(line_no, line, "encoded '" + got + "'");
            } else if (tag == "kvalue") {
                const KValue k = kvalue_of(input);
                const std::string got = encode_kvalue(k);
                if (got != want || !same_kvalue(k, decode_kvalue(want)))
                    report(line_no, line, "encoded '" + got + "'");
            } else {
                report(line_no, line, "unknown tag");
            }
        } catch (const std::exception &e) {
            report(line_no, line, std::string("threw: ") + e.what());
        }
    }
    EXPECT_EQ(0u, failures);
    // every vector was seen: the file's own COUNTS line
    ASSERT_FALSE(expected_counts.empty()) << "no COUNTS line";
    size_t total = 0;
    for (const auto &[tag, n] : seen) {
        EXPECT_EQ(expected_counts[tag], n) << tag;
        total += n;
    }
    EXPECT_EQ(expected_counts["total"], total);
    EXPECT_EQ(expected_counts.size(), seen.size() + 1);
}

// FIPS 180-4 examples: the index fingerprint of §3.1 is a sha256 the deployment can
// recompute with any tool
TEST(GraphletCodec, Sha256) {
    EXPECT_EQ("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", sha256_hex(""));
    EXPECT_EQ("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", sha256_hex("abc"));
    EXPECT_EQ("248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
              sha256_hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"));
    EXPECT_EQ("cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0",
              sha256_hex(std::string(1000000, 'a')));
}


// Label names come from FASTA headers and file names, which need not be UTF-8. No output
// carries such a name verbatim, and a REPLACED one (U+FFFD) can be another label's name in
// the index -- a continuation that resubmitted it went on under that label (GPT review,
// finding 1) -- so a seed with such a name is refused (Graphlet.NamesThatAreNotUtf8...);
// valid_utf8 is the test applied there and in the writers' backstops.
TEST(GraphletCodec, NamesThatAreNotUtf8AreRecognised) {
    const std::vector<std::pair<std::string, bool>> cases = {
        { "abc", true },
        { "caf\xC3\xA9", true },
        { "\xF0\x9F\x98\x80" "z", true },
        { "bad\xEF\xBF\xBD", true },          // U+FFFD itself is a valid name
        { "", true },
        { "\xC0" "abc", false },                 // C0: never a lead byte
        { "\xE9" "ab", false },                  // Latin-1: a lead without its tail
        { "\xED\xA0\x80", false },              // a surrogate
        { "\xE0\x80\x80", false },              // overlong
        { "\xF4\x90\x80\x80", false },          // beyond U+10FFFF
        { "\xF0\x9F\x98", false },              // truncated
        { "\xE2\x82" "x", false },
        { "bad\xFF", false },
        { "\x80\x80", false },
    };
    for (const auto &[raw, valid] : cases) {
        EXPECT_EQ(valid, valid_utf8(raw)) << raw;
    }
}

// RANGES at the top of the id space: nothing follows UINT64_MAX. MAX + 1 wraps to 0, so
// the decoder accepted "18446744073709551615,1" (the descending {MAX, 1}; the Python
// reader rejects it) and the encoder wrote {MAX, 0} as the descending run "MAX-0" (GPT
// review, finding 8). The shared golden vectors have no case at MAX (codec_vectors.tsv is
// generated by scripts/traversal/make_codec_vectors.py; not edited here).
TEST(GraphletCodec, RangesAtTheTopOfTheIdSpace) {
    const uint64_t MAX = std::numeric_limits<uint64_t>::max();
    const std::string max = "18446744073709551615", below = "18446744073709551614";
    EXPECT_THROW(decode_ranges(max + ",1"), FormatError);
    EXPECT_THROW(decode_ranges(max + ",0"), FormatError);
    EXPECT_THROW(decode_ranges(max + "," + max), FormatError);
    EXPECT_THROW(decode_ranges(below + "-" + max + ",0"), FormatError);
    EXPECT_THROW(decode_ranges(below + "," + max), FormatError);      // a run of two
    EXPECT_EQ(std::vector<uint64_t>({ MAX }), decode_ranges(max));
    EXPECT_EQ(std::vector<uint64_t>({ MAX - 1, MAX }), decode_ranges(below + "-" + max));
    EXPECT_EQ(std::vector<uint64_t>({ 0, MAX }), decode_ranges("0," + max));
    EXPECT_THROW(encode_ranges({ MAX, 0 }), std::logic_error);
    EXPECT_THROW(encode_ranges({ MAX - 1, MAX, 0 }), std::logic_error);
    EXPECT_THROW(encode_ranges({ MAX, MAX }), std::logic_error);
    EXPECT_EQ(below + "-" + max, encode_ranges({ MAX - 1, MAX }));
    EXPECT_EQ("0," + max, encode_ranges({ 0, MAX }));
    EXPECT_EQ(max, encode_ranges({ MAX }));
    // the delta parts of a SETEXPR are RANGES too
    const std::vector<uint64_t> base = { 0 };
    EXPECT_THROW(decode_setexpr("!+" + max + ",1", &base), FormatError);
    EXPECT_THROW(encode_setexpr({ MAX, 0 }, nullptr), std::logic_error);
}

// A RANGES token is checked against the id space before a run is expanded: a 20-byte
// token must not allocate gigabytes before it is rejected.
TEST(GraphletCodec, RangesAreBoundedBeforeExpansion) {
    EXPECT_EQ(std::vector<uint64_t>({ 0, 1, 2, 7 }), decode_ranges("0-2,7", 8));
    EXPECT_THROW(decode_ranges("0-2,8", 8), FormatError);
    EXPECT_THROW(decode_ranges("0-3000000000", 3), FormatError);
    EXPECT_THROW(decode_ranges("0-18446744073709551615", 3), FormatError);
    const std::vector<uint64_t> base = { 0, 1 };
    EXPECT_THROW(decode_setexpr("!+2-4000000000", &base, 3), FormatError);
    EXPECT_EQ(std::vector<uint64_t>({ 0, 2 }), decode_setexpr("!1+2", &base, 3));
}

/*
 * Whole documents (DESIGN-traverse-graphlet.md §2.2-§2.5, the C++ side of T37/T40). An
 * independent reader in this file parses the MGT text of `detail: graphlet` and rebuilds
 * from it, with the summary beside it, the `detail: full` result of the SAME request: every
 * `*` the writer emitted must decode by its rule to the walker's value, every derived
 * record (paths, splits, label ends, reconvergences, end labels, continuations, the label
 * summary) must come out as the server's own JSON, and the counts in A and Z must hold.
 * The fixtures cover merges and keeps, both label modes, switching with silent
 * switch-source ends, quorum and branch limits, caps, cut label lists, a beam, followed
 * hairpins, sequences: false and a batch with a duplicate and a failed derivation.
 */

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::graph::traversal;

std::string compact(const Json::Value &v) {
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    return Json::writeString(b, v);
}

Json::Value parse_json(const std::string &text) {
    Json::Value v;
    std::istringstream in(text);
    Json::CharReaderBuilder b;
    std::string errs;
    if (!Json::parseFromStream(b, in, &v, &errs))
        throw std::runtime_error(errs);
    return v;
}

std::vector<std::string> fields(const std::string &line, size_t max_fields = SIZE_MAX) {
    std::vector<std::string> out;
    size_t pos = 0;
    while (out.size() + 1 < max_fields) {
        const size_t space = line.find(' ', pos);
        if (space == std::string::npos)
            break;
        out.push_back(line.substr(pos, space - pos));
        pos = space + 1;
    }
    out.push_back(line.substr(pos));
    return out;
}

std::vector<std::string> split_on(const std::string &s, char sep) {
    std::vector<std::string> out;
    size_t pos = 0;
    while (true) {
        const size_t next = s.find(sep, pos);
        out.push_back(s.substr(pos, next == std::string::npos ? std::string::npos : next - pos));
        if (next == std::string::npos)
            return out;
        pos = next + 1;
    }
}

uint64_t u(const std::string &s) {
    size_t used = 0;
    const uint64_t x = std::stoull(s, &used);
    if (used != s.size() || (s.size() > 1 && s[0] == '0'))
        throw std::runtime_error("not a canonical integer: " + s);
    return x;
}

const char kCodes[] = "DLBRUVJTXSPNOMWY";

EndReason reason_of(char code) {
    const char *p = std::strchr(kCodes, code);
    if (!p || !*p)
        throw std::runtime_error(std::string("unknown reason code ") + code);
    return static_cast<EndReason>(p - kCodes);
}

std::string qualifier_text(char code, char q) {
    if (!q) return "";
    if (code == 'R' && q == 'm') return "minority";
    if (code == 'R' && q == 'b') return "below_min_labels";
    if (code == 'R' && q == 's') return "split_limit";
    if (code == 'D' && q == 'h') return "hairpin";
    if (code == 'L' && q == 's') return "superseded";
    if (code == 'L' && q == 'x') return "switch_sources";
    throw std::runtime_error(std::string("unknown qualifier ") + code + q);
}

Json::Value ids_json(const std::vector<uint64_t> &ids) {
    Json::Value a(Json::arrayValue);
    for (uint64_t x : ids) a.append(Json::UInt64(x));
    return a;
}

std::vector<uint64_t> unite(const std::vector<std::vector<uint64_t>> &sets) {
    std::vector<uint64_t> out, next;
    for (const auto &s : sets) {
        next.clear();
        std::set_union(out.begin(), out.end(), s.begin(), s.end(), std::back_inserter(next));
        out.swap(next);
    }
    return out;
}

Json::Value kvalue_json(const KValue &k) {
    switch (k.type) {
        case KValue::INTEGER: return Json::UInt64(k.integer);
        case KValue::FLOAT: return k.real;
        case KValue::STRING: return k.string;
        case KValue::UNLIMITED: return "unlimited";
    }
    return Json::Value();
}

struct RunRec {
    size_t segment;
    uint64_t label, from, to, route;
    std::string end;               // code (+ qualifier) | Lw | m
    bool by_switch = false;
    uint64_t from_label = 0;
    double cost = 0;
    std::optional<uint64_t> prev;
    std::optional<uint64_t> structural;
    uint64_t branches = 0;
    double loss = 0;
    std::optional<double> needed;
};

struct SegRec {
    std::vector<uint64_t> parents;
    uint64_t from = 0, length = 0;
    std::string entry, entry_total, end, partition, split, first_base, bases;
    std::vector<std::vector<std::string>> presence;   // P fields
    std::vector<std::vector<std::string>> events;     // E fields
    std::optional<std::vector<std::string>> leaf;     // T fields
    std::optional<std::vector<std::string>> cont;     // C fields
    // resolved
    std::vector<uint64_t> entry_set, end_set;
    std::vector<std::vector<uint64_t>> via;
    std::vector<std::vector<uint64_t>> sets;          // P sets
    std::vector<uint64_t> totals;                     // P totals
    uint64_t total = 0;
};

struct ArmRec {
    std::vector<std::string> a;
    std::vector<std::vector<std::string>> bins, events;
    std::vector<SegRec> segs;
    std::vector<RunRec> runs;
    Json::Value limitations = Json::Value(Json::arrayValue);
};

// The reader: the body (and its summary) back to today's `detail: full` result.
struct Reader {
    // H
    bool annotate = false, sequences_absent = false;
    uint64_t cap = 0, nsl = 0, seed_len = 0;
    std::string seq;
    std::vector<Json::Value> dict;
    Json::Value dropped = Json::Value(Json::arrayValue);
    Json::Value seed_limitations = Json::Value(Json::arrayValue);
    Json::Value outcome;
    Json::Value resource_stop;     // null: no Q record
    std::map<char, ArmRec> arms;
    std::vector<char> order;
    size_t lines = 0, z = 0;

    explicit Reader(const std::string &text) {
        EXPECT_FALSE(text.empty());
        EXPECT_EQ('\n', text.back()) << "a document ends with a line break";
        std::vector<std::string> ls = split_on(text.substr(0, text.size() - 1), '\n');
        lines = ls.size();
        std::string previous_name;
        ArmRec *arm = nullptr;
        for (const std::string &line : ls) {
            const char type = line[0];
            EXPECT_EQ(' ', line[1]) << line;
            switch (type) {
                case 'H': {
                    auto f = fields(line);
                    EXPECT_EQ(16u, f.size()) << line;
                    EXPECT_EQ("mgt", f[1]);
                    EXPECT_EQ("1", f[2]);
                    EXPECT_EQ("walk", f[12]);
                    annotate = f[6] == "a";
                    cap = u(f[9]);
                    break;
                }
                case 'S': {
                    auto f = fields(line);
                    seed_len = u(f[2]);
                    nsl = u(f[4]);
                    seq = f[5];
                    EXPECT_EQ(seed_len, seq.size());
                    break;
                }
                case 'X': {
                    auto f = fields(line, 4);
                    Json::Value d;
                    d["label"] = pct_unescape(f[3]);
                    d["reason"] = f[1];
                    d["runs"] = Json::Value(Json::arrayValue);
                    if (f[2] != ".") {
                        for (const std::string &r : split_on(f[2], ',')) {
                            auto ab = split_on(r, '-');
                            Json::Value iv(Json::arrayValue);
                            iv.append(Json::UInt64(u(ab[0])));
                            iv.append(Json::UInt64(u(ab[1])));
                            d["runs"].append(iv);
                        }
                    }
                    d["runs_kind"] = "presence";
                    dropped.append(d);
                    break;
                }
                case 'L': {
                    auto f = fields(line, 5);
                    Json::Value l;
                    const std::string name = front_decode(previous_name, f[4]);
                    // canonical: the writer's prefix is the longest
                    EXPECT_EQ(f[4], front_code(previous_name, name));
                    l["name"] = name;
                    l["kind"] = f[1] == "h" ? "header" : "column";
                    l["column"] = Json::UInt64(u(f[2]));
                    if (f[1] == "h")
                        l["seq_id"] = Json::UInt64(u(f[3]));
                    else
                        EXPECT_EQ("*", f[3]);
                    dict.push_back(l);
                    previous_name = name;
                    break;
                }
                case 'O': {
                    auto f = fields(line);
                    auto pick = [](const std::string &c, std::map<std::string, std::string> m) {
                        EXPECT_TRUE(m.count(c)) << c;
                        return m[c];
                    };
                    outcome["walks"] = pick(f[1], { { "c", "complete" }, { "p", "partial" }, { "f", "failed" } });
                    outcome["branch_diagnostics"] = pick(f[2], { { "c", "complete" }, { "x", "cut" } });
                    outcome["label_evidence"] = pick(f[3], { { "c", "complete" }, { "l", "lower_bound" }, { "q", "qualified" } });
                    outcome["delivery"] = pick(f[4], { { "i", "inline" }, { "s", "spooled" }, { "p", "paged" } });
                    break;
                }
                case 'Q': {
                    // a resource stop (DESIGN §14): typed amounts (* = not stated), the
                    // actions, and the message as the free-text last field
                    auto f = fields(line, 10);
                    EXPECT_EQ(10u, f.size()) << line;
                    EXPECT_TRUE(resource_stop.isNull()) << "one Q record at most";
                    Json::Value q;
                    q["scope"] = f[1];
                    q["resource"] = f[2];
                    q["phase"] = f[3];
                    const char *amounts[] = { "requested", "effective", "used", "remaining" };
                    for (size_t i = 0; i < 4; ++i) {
                        q[amounts[i]] = f[4 + i] == "*" ? Json::Value()
                                                        : kvalue_json(decode_kvalue(f[4 + i]));
                    }
                    q["actions"] = Json::Value(Json::arrayValue);
                    if (f[8] != ".") {
                        for (const std::string &a : split_on(f[8], ',')) q["actions"].append(a);
                    }
                    q["message"] = pct_unescape(f[9]);
                    resource_stop = q;
                    break;
                }
                case 'K': {
                    auto f = fields(line, 9);
                    Json::Value l;
                    l["kind"] = f[2];
                    l["knob"] = f[3];
                    l["limit"] = kvalue_json(decode_kvalue(f[4]));
                    l["observed"] = kvalue_json(decode_kvalue(f[5]));
                    if (f[6] != "*")
                        l["complete_to_bp"] = Json::UInt64(u(f[6]));
                    if (f[7] != ".") {
                        for (const std::string &item : split_on(f[7], ',')) {
                            const size_t eq = item.find('=');
                            l[item.substr(0, eq)] = kvalue_json(decode_kvalue(item.substr(eq + 1)));
                        }
                    }
                    l["effect"] = pct_unescape(f[8]);
                    if (f[1] == "*") {
                        seed_limitations.append(l);
                    } else {
                        arms[f[1][0]].limitations.append(l);
                    }
                    break;
                }
                case 'A': {
                    auto f = fields(line);
                    EXPECT_EQ(20u, f.size()) << line;
                    arm = &arms[f[1][0]];
                    arm->a = f;
                    order.push_back(f[1][0]);
                    break;
                }
                case 'B': arm->bins.push_back(fields(line)); break;
                case 'V': arm->events.push_back(fields(line)); break;
                case 'G': {
                    auto f = fields(line);
                    EXPECT_EQ(11u, f.size()) << line;
                    SegRec s;
                    if (f[1] != "*") {
                        for (const std::string &p : split_on(f[1], ',')) s.parents.push_back(u(p));
                    }
                    s.from = u(f[2]);
                    s.length = u(f[3]);
                    s.entry = f[4]; s.entry_total = f[5]; s.end = f[6]; s.partition = f[7];
                    s.split = f[8]; s.first_base = f[9]; s.bases = f[10];
                    sequences_absent = f[10] == "*";
                    arm->segs.push_back(s);
                    break;
                }
                case 'P': arm->segs.back().presence.push_back(fields(line)); break;
                case 'E': arm->segs.back().events.push_back(fields(line)); break;
                case 'T': arm->segs.back().leaf = fields(line); break;
                case 'C': arm->segs.back().cont = fields(line); break;
                case 'R': {
                    auto f = fields(line);
                    RunRec r;
                    r.segment = u(f[1]); r.label = u(f[2]); r.from = u(f[3]); r.to = u(f[4]);
                    r.end = f[5]; r.route = u(f[6]);
                    if (f[7] != "*") {
                        auto lc = split_on(f[7], ':');
                        r.by_switch = true;
                        r.from_label = u(lc[0]);
                        r.cost = decode_float(lc[1]);
                    }
                    if (f[8] != "*") r.prev = u(f[8]);
                    if (f[9] != "*") r.structural = u(f[9]);
                    r.branches = u(f[10]);
                    r.loss = decode_float(f[11]);
                    if (f.size() > 12) r.needed = decode_float(f[12]);
                    EXPECT_EQ(r.end[0] == 'B' ? 13u : 12u, f.size()) << line;
                    arm->runs.push_back(r);
                    break;
                }
                case 'Z': {
                    z = u(fields(line)[1]);
                    break;
                }
                default:
                    ADD_FAILURE() << "unknown record " << line;
            }
        }
        EXPECT_EQ(lines, z) << "Z counts the lines of the document";
        for (auto &[side, a] : arms) {
            resolve(side, a);
        }
    }

    std::string bases_of(const SegRec &s) const { return s.bases == "." ? "" : s.bases; }

    // the §2.3 rules, applied where the writer wrote '*'
    void resolve(char side, ArmRec &arm) {
        (void)side;
        for (size_t i = 0; i < arm.segs.size(); ++i) {
            SegRec &s = arm.segs[i];
            const bool root = s.parents.empty();
            if (s.partition == "*") {
                if (s.parents.size() > 1 && annotate)
                    s.via.resize(s.parents.size());
            } else {
                for (const std::string &p : split_on(s.partition, '|')) s.via.push_back(decode_ranges(p));
                EXPECT_EQ(s.parents.size(), s.via.size());
            }
            std::vector<uint64_t> base;
            if (s.parents.size() == 1) {
                base = arm.segs[s.parents[0]].end_set;
            } else if (s.parents.size() > 1) {
                std::vector<std::vector<uint64_t>> ends;
                for (uint64_t p : s.parents) ends.push_back(arm.segs[p].end_set);
                base = unite(ends);
            }
            if (s.entry == "*") {
                EXPECT_FALSE(annotate) << "G.entry is never * in annotate mode";
                if (root) {
                    for (uint64_t l = 0; l < nsl; ++l) s.entry_set.push_back(l);
                } else {
                    EXPECT_GT(s.parents.size(), 1u);
                    s.entry_set = unite(s.via);
                }
            } else {
                s.entry_set = decode_setexpr(s.entry, root ? nullptr : &base);
                // canonical: the shortest form, and never explicit where the rule holds
                EXPECT_EQ(s.entry, encode_setexpr(s.entry_set, root ? nullptr : &base));
                if (!annotate && root) {
                    std::vector<uint64_t> seeds;
                    for (uint64_t l = 0; l < nsl; ++l) seeds.push_back(l);
                    EXPECT_NE(seeds, s.entry_set) << "G.entry explicit where * holds";
                } else if (!annotate && s.parents.size() > 1) {
                    EXPECT_NE(unite(s.via), s.entry_set) << "G.entry explicit where * holds";
                }
            }
            if (s.partition != "*") {
                EXPECT_FALSE(annotate || s.parents.size() < 2) << "G.partition explicit where * holds";
            }
            s.total = s.entry_total == "*" ? s.entry_set.size() : u(s.entry_total);
            if (s.entry_total != "*") {
                EXPECT_NE(s.total, s.entry_set.size());
            }
            std::vector<uint64_t> previous = s.entry_set;
            for (const auto &p : s.presence) {
                std::vector<uint64_t> set = decode_setexpr(p[4], &previous);
                EXPECT_EQ(p[4], encode_setexpr(set, &previous));
                s.totals.push_back(u(p[3]));
                s.sets.push_back(set);
                previous = std::move(set);
            }
            std::vector<uint64_t> rule;
            if (annotate) {
                rule = s.sets.empty() ? s.entry_set : s.sets.back();
            } else {
                // chronological: ends before switch-ins at one position
                std::vector<std::tuple<uint64_t, int, uint64_t>> changes;
                for (const RunRec &r : arm.runs) {
                    if (r.segment == i && r.to < s.from + s.length)
                        changes.emplace_back(r.to, 0, r.label);
                }
                for (const auto &e : s.events) {
                    if (e[2] == "s")
                        changes.emplace_back(u(e[1]), 1, u(e[4]));
                }
                std::sort(changes.begin(), changes.end());
                std::set<uint64_t> cur(s.entry_set.begin(), s.entry_set.end());
                for (const auto &[at, kind, label] : changes) {
                    if (kind == 0) cur.erase(label); else cur.insert(label);
                }
                rule.assign(cur.begin(), cur.end());
            }
            if (s.end == "*") {
                s.end_set = rule;
            } else {
                s.end_set = decode_setexpr(s.end, &s.entry_set);
                EXPECT_NE(rule, s.end_set) << "G.end written explicitly where the rule holds";
                EXPECT_EQ(s.end, encode_setexpr(s.end_set, &s.entry_set));
            }
        }
    }

    // today's per-seed `detail: full` result, from the body and the summary
    Json::Value to_json(const Json::Value &summary) const {
        Json::Value j;
        Json::Value seed = summary["seed"];
        seed.removeMember("num_labels");
        seed.removeMember("num_seed_labels");
        Json::Value labels(Json::arrayValue);
        for (uint64_t i = 0; i < nsl; ++i) labels.append(dict[i]["name"]);
        seed["labels"] = labels;
        seed["dropped_labels"] = dropped;
        j["seed"] = seed;
        j["limitations"] = seed_limitations;
        j["outcome"] = outcome;
        if (!resource_stop.isNull())
            j["resource_stop"] = resource_stop;
        j["label_mode"] = annotate ? "annotate" : "constrain";
        j["label_dict"] = Json::Value(Json::arrayValue);
        for (const auto &l : dict) j["label_dict"].append(l);
        j["annotation"] = summary["annotation"];
        if (summary.isMember("duplicate")) j["duplicate"] = summary["duplicate"];
        // record coordinates are no part of the body (MGT v1 is frozen, §18): the envelope
        // carries them, the same block in every detail
        for (const char *key : { "coordinates", "coordinates_reason" }) {
            if (summary.isMember(key))
                j[key] = summary[key];
        }
        std::vector<std::array<Json::Value, 2>> ls(dict.size());
        for (const auto &[side, arm] : arms) {
            const std::string name = side == 'l' ? "left" : "right";
            j["arms"][name] = arm_json(side, arm, &ls, name);
        }
        Json::Value summary_json(Json::arrayValue);
        for (size_t l = 0; l < dict.size(); ++l) {
            Json::Value lj;
            lj["label"] = Json::UInt64(l);
            for (const auto &[side, arm] : arms) {
                lj[side == 'l' ? "left" : "right"] = ls[l][side == 'l' ? 0 : 1];
            }
            summary_json.append(lj);
        }
        j["label_summary"] = summary_json;
        return j;
    }

    Json::Value arm_json(char side, const ArmRec &arm, std::vector<std::array<Json::Value, 2>> *ls,
                         const std::string &name) const {
        const auto &a = arm.a;
        const auto &segs = arm.segs;
        Json::Value j;
        j["status"] = a[2] == "c" ? "complete" : a[2] == "t" ? "truncated" : "pruned";
        j["complete_to_bp"] = Json::UInt64(u(a[3]));
        j["completeness_scope"] = a[4] == "u" ? "united_history" : "per_path";
        j["frontier_remaining"]["live_paths"] = Json::UInt64(u(a[5]));
        j["frontier_remaining"]["live_labels"] = Json::UInt64(u(a[6]));
        j["frontier_remaining"]["exact"] = a[7] == "1";
        j["labels_per_node"]["cap"] = Json::UInt64(cap);
        j["labels_per_node"]["max_seen"] = Json::UInt64(u(a[8]));
        j["labels_per_node"]["nodes_truncated"] = Json::UInt64(u(a[9]));
        Json::Value counters(Json::objectValue);
        for (const std::string &item : split_on(a[10], ',')) {
            const size_t eq = item.find('=');
            counters[item.substr(0, eq)] = Json::UInt64(u(item.substr(eq + 1)));
        }
        j["counters"] = counters;
        if (a[11] != "*") {
            auto c = split_on(a[11], ',');
            j["cap_trigger"]["reason"] = to_string(reason_of(c[0][0]));
            j["cap_trigger"]["at_bp"] = Json::UInt64(u(c[1]));
            j["cap_trigger"]["segment"] = Json::UInt64(u(c[2]));
            j["cap_trigger"]["live_paths"] = Json::UInt64(u(c[3]));
            j["cap_trigger"]["live_labels"] = Json::UInt64(u(c[4]));
            j["cap_trigger"]["exact"] = c[5] == "1";
            decode_float(c[6]);
        }
        j["evidence"]["complete"] = a[13] == "*";
        j["evidence"]["complete_to_bp"] = a[13] == "*" ? Json::Value() : Json::Value(Json::UInt64(u(a[13])));
        j["limitations"] = arm.limitations;
        // the counts that validate the body
        size_t leaves = 0, splits = 0, merges = 0;
        uint64_t bases = 0;
        for (const SegRec &s : segs) {
            leaves += bool(s.leaf);
            splits += s.split != "*";
            merges += s.parents.size() > 1;
            bases += s.length;
            if (!sequences_absent) {
                EXPECT_EQ(s.length, bases_of(s).size());
            }
        }
        EXPECT_EQ(u(a[14]), segs.size());
        EXPECT_EQ(u(a[15]), arm.runs.size());
        EXPECT_EQ(u(a[16]), leaves);
        EXPECT_EQ(u(a[17]), splits);
        EXPECT_EQ(u(a[18]), merges);
        EXPECT_EQ(u(a[19]), bases);
        // ---- segments
        auto walk_of = [&](size_t s) { return bases_of(segs[s]); };
        Json::Value sj_all(Json::arrayValue);
        std::vector<std::vector<size_t>> children(segs.size());
        for (size_t i = 0; i < segs.size(); ++i) {
            const SegRec &s = segs[i];
            if (s.parents.size() == 1)
                children[s.parents[0]].push_back(i);
            Json::Value sj;
            sj["id"] = Json::UInt64(i);
            sj["parents"] = ids_json(s.parents);
            if (s.parents.size() > 1) {
                sj["labels_via_parent"] = Json::Value(Json::arrayValue);
                for (const auto &v : s.via) sj["labels_via_parent"].append(ids_json(v));
            }
            sj["from_bp"] = Json::UInt64(s.from);
            sj["length_bp"] = Json::UInt64(s.length);
            if (!sequences_absent) {
                std::string natural = walk_of(i);
                if (side == 'l') std::reverse(natural.begin(), natural.end());
                sj["sequence"] = natural;
            }
            sj["labels"] = ids_json(s.entry_set);
            sj["labels_total"] = Json::UInt64(s.total);
            sj["labels_truncated"] = s.total > s.entry_set.size();
            sj["labels_at_end"] = ids_json(s.end_set);
            Json::Value evs(Json::arrayValue);
            for (const auto &e : s.events) {
                Json::Value ej;
                ej["at_bp"] = Json::UInt64(u(e[1]));
                if (e[2] == "s") {
                    ej["type"] = "switch";
                    ej["from"] = Json::UInt64(u(e[3]));
                    ej["to"] = Json::UInt64(u(e[4]));
                    ej["cost"] = decode_float(e[5]);
                } else if (e[2] == "b") {
                    ej["type"] = "blocked";
                    ej["char"] = e[3];
                    ej["reason"] = to_string(reason_of(e[4][0]));
                    ej["labels"] = ids_json(decode_ranges(e[6]));
                    ej["labels_total"] = Json::UInt64(u(e[5]));
                    ej["truncated"] = u(e[5]) > decode_ranges(e[6]).size();
                } else if (e[2] == "h") {
                    ej["type"] = "hairpin";
                    ej["char"] = e[3];
                    ej["labels"] = ids_json(decode_ranges(e[5]));
                    ej["labels_total"] = Json::UInt64(u(e[4]));
                    ej["truncated"] = u(e[4]) > decode_ranges(e[5]).size();
                    ej["followed"] = e[6] == "f";
                } else if (e[2] == "v") {
                    ej["type"] = "revisit";
                    ej["segments"] = ids_json({ u(e[3]) });
                    ej["length_bp"] = Json::UInt64(e[4] == "=" ? 0 : u(e[4]));
                    ej["same_distance"] = e[4] == "=";
                } else {
                    ADD_FAILURE() << "unexpected event " << e[2];
                }
                evs.append(ej);
            }
            for (const RunRec &r : arm.runs) {
                if (r.segment != i || r.end == "Lw" || r.end == "m")
                    continue;
                Json::Value ej;
                ej["at_bp"] = Json::UInt64(r.to);
                ej["type"] = "label_end";
                ej["label"] = Json::UInt64(r.label);
                const std::string text = qualifier_text(r.end[0], r.end.size() > 1 ? r.end[1] : 0);
                ej["reason"] = text.empty() ? to_string(reason_of(r.end[0])) : text.c_str();
                if (r.end[0] == 'B') ej["needed_budget"] = *r.needed;
                ej["structural_successors"] = Json::UInt64(*r.structural);
                evs.append(ej);
            }
            if (s.parents.size() > 1) {
                Json::Value ej;
                ej["at_bp"] = Json::UInt64(s.from);
                ej["type"] = "reconverge";
                ej["segments"] = ids_json(s.parents);
                evs.append(ej);
            }
            sj["events"] = evs;
            if (annotate) {
                Json::Value sets(Json::arrayValue);
                for (size_t k = 0; k < s.presence.size(); ++k) {
                    Json::Value rj;
                    rj["from_bp"] = Json::UInt64(u(s.presence[k][1]));
                    rj["to_bp"] = Json::UInt64(u(s.presence[k][2]));
                    rj["labels"] = ids_json(s.sets[k]);
                    rj["labels_total"] = Json::UInt64(s.totals[k]);
                    rj["truncated"] = s.totals[k] > s.sets[k].size();
                    sets.append(rj);
                }
                sj["label_sets"] = sets;
            }
            sj_all.append(sj);
        }
        j["segments"] = sj_all;
        // ---- splits: single-parent children grouped by parent, by (at_bp, first child)
        std::vector<std::pair<std::pair<uint64_t, size_t>, size_t>> split_order;
        for (size_t i = 0; i < segs.size(); ++i) {
            EXPECT_EQ(segs[i].split != "*", !children[i].empty()) << "G.split on a split";
            if (!children[i].empty())
                split_order.push_back({ { segs[i].from + segs[i].length, children[i][0] }, i });
        }
        std::sort(split_order.begin(), split_order.end());
        Json::Value splits_json(Json::arrayValue);
        for (const auto &[key, i] : split_order) {
            const SegRec &p = segs[i];
            Json::Value sj;
            sj["at_bp"] = Json::UInt64(key.first);
            sj["prefix_bp"] = Json::UInt64(key.first);
            sj["segment"] = Json::UInt64(i);
            sj["children"] = Json::Value(Json::arrayValue);
            for (size_t c : children[i]) sj["children"].append(Json::UInt64(c));
            sj["kind"] = p.split == "1" ? "ambiguous" : "divergence";
            sj["labels_before"] = Json::UInt64(annotate ? (p.totals.empty() ? p.total : p.totals.back())
                                                        : p.end_set.size());
            sj["branches"] = Json::Value(Json::arrayValue);
            for (size_t c : children[i]) {
                const SegRec &ch = segs[c];
                Json::Value bj;
                bj["segment"] = Json::UInt64(c);
                bj["char"] = sequences_absent ? ch.first_base : walk_of(c).substr(0, 1);
                bj["labels_distinct"] = Json::UInt64(ch.total);
                std::vector<uint64_t> cut(ch.entry_set.begin(),
                                          ch.entry_set.begin() + std::min<size_t>(cap, ch.entry_set.size()));
                bj["labels"] = ids_json(cut);
                bj["labels_truncated"] = ch.total > cut.size();
                sj["branches"].append(bj);
            }
            splits_json.append(sj);
        }
        j["splits"] = splits_json;
        // ---- paths: leaf ordinal in segment order, first-parent chain
        Json::Value paths(Json::arrayValue);
        for (size_t i = 0; i < segs.size(); ++i) {
            const SegRec &s = segs[i];
            if (!s.leaf)
                continue;
            Json::Value pj;
            pj["id"] = Json::UInt64(paths.size());
            std::vector<uint64_t> chain;
            for (size_t c = i; ; c = segs[c].parents[0]) {
                chain.insert(chain.begin(), c);
                if (segs[c].parents.empty()) break;
            }
            pj["segments"] = ids_json(chain);
            const uint64_t length = s.from + s.length;
            pj["length_bp"] = Json::UInt64(length);
            std::map<uint64_t, std::tuple<double, uint64_t, uint64_t>> extras;
            if ((*s.leaf)[2] != ".") {
                for (const std::string &x : split_on((*s.leaf)[2], ',')) {
                    auto f = split_on(x, ':');
                    extras[u(f[0])] = { decode_float(f[1]), u(f[2]), u(f[3]) };
                }
            }
            Json::Value reasons(Json::objectValue), ends(Json::arrayValue);
            std::vector<std::pair<uint64_t, size_t>> alive;
            for (size_t r = 0; r < arm.runs.size(); ++r) {
                if (arm.runs[r].segment == i && arm.runs[r].to == length)
                    alive.emplace_back(arm.runs[r].label, r);
            }
            std::sort(alive.begin(), alive.end());
            for (const auto &[label, r] : alive) {
                const std::string reason = to_string(reason_of(arm.runs[r].end[0]));
                reasons[reason] = reasons[reason].asUInt() + 1;
                Json::Value ej;
                auto it = extras.find(label);
                ej["label"] = Json::UInt64(label);
                ej["loss"] = it == extras.end() ? 0.0 : std::get<0>(it->second);
                ej["branches"] = Json::UInt64(it == extras.end() ? 0 : std::get<1>(it->second));
                ej["run"] = Json::UInt64(r);
                ej["route_bp"] = Json::UInt64(it == extras.end() ? 0 : std::get<2>(it->second));
                ends.append(ej);
            }
            pj["end_reasons"] = reasons;
            if ((*s.leaf)[1] != "*") pj["path_reason"] = to_string(reason_of((*s.leaf)[1][0]));
            pj["n_labels"] = Json::UInt64(alive.size());
            pj["end_labels"] = ends;
            if (s.cont) {
                const auto &c = *s.cont;
                const uint64_t n = u(c[1]);
                std::string sequence;
                if (c.size() > 5) {
                    sequence = c[5];
                    EXPECT_TRUE(sequences_absent && n > 0) << "C.sequence written where it derives";
                } else {
                    // right: the last n of seed + natural flank; left: the first n of
                    // natural flank + seed
                    std::string walk;
                    for (uint64_t x : chain) walk += walk_of(x);
                    if (side == 'r') {
                        const std::string all = seq + walk;
                        sequence = all.substr(all.size() - n);
                    } else {
                        std::string natural(walk.rbegin(), walk.rend());
                        sequence = (natural + seq).substr(0, n);
                    }
                }
                EXPECT_EQ(n, sequence.size());
                pj["continuation"]["sequence"] = sequence;
                pj["continuation"]["labels"] = ids_json(decode_setexpr(c[4], &s.end_set));
                pj["continuation"]["loss_used"] = decode_float(c[2]);
                pj["continuation"]["branches_used"] = Json::UInt64(u(c[3]));
            }
            paths.append(pj);
        }
        j["paths"] = paths;
        // ---- runs
        Json::Value runs(Json::arrayValue);
        std::vector<double> needed;
        for (const RunRec &r : arm.runs) {
            Json::Value rj;
            rj["label"] = Json::UInt64(r.label);
            rj["from_bp"] = Json::UInt64(r.from);
            rj["to_bp"] = Json::UInt64(r.to);
            rj["entered_by"] = r.by_switch ? "switch" : "seed";
            rj["route_bp"] = Json::UInt64(r.route);
            if (r.by_switch) { rj["from"] = Json::UInt64(r.from_label); rj["cost"] = r.cost; }
            if (r.end != "m") rj["end_reason"] = to_string(reason_of(r.end[0]));
            rj["segment"] = Json::UInt64(r.segment);
            rj["branches"] = Json::UInt64(r.branches);
            rj["loss"] = r.loss;
            if (r.needed) needed.push_back(*r.needed);
            runs.append(rj);
        }
        j["runs"] = runs;
        std::sort(needed.begin(), needed.end());
        j["needed_budgets"] = Json::Value(Json::arrayValue);
        for (double x : needed) j["needed_budgets"].append(x);
        // ---- growth and branch events
        Json::Value growth(Json::arrayValue);
        for (const auto &b : arm.bins) {
            Json::Value g;
            const char *keys[] = { "from_bp", "max_live_paths", "distinct_live_labels", "live_pairs" };
            for (size_t k = 0; k < 4; ++k) g[keys[k]] = Json::UInt64(u(b[1 + k]));
            g["exact"] = b[5] == "1";
            const char *more[] = { "steps", "divergences", "ambiguous_branches", "splits",
                                   "reconvergences", "bubbles", "tips", "blocked_repeat" };
            for (size_t k = 0; k < 8; ++k) g[more[k]] = Json::UInt64(u(b[6 + k]));
            g["label_ends"] = Json::Value(Json::objectValue);
            if (b[14] != ".") {
                for (const std::string &x : split_on(b[14], ',')) {
                    g["label_ends"][to_string(reason_of(x[0]))] = Json::UInt64(u(x.substr(2)));
                }
            }
            growth.append(g);
        }
        j["growth"] = growth;
        Json::Value bes(Json::arrayValue);
        for (const auto &v : arm.events) {
            Json::Value bj;
            bj["at_bp"] = Json::UInt64(u(v[1]));
            bj["segment"] = Json::UInt64(u(v[2]));
            bj["successors"] = v[3] == "." ? "" : v[3];
            bj["labels_per_successor"] = Json::Value(Json::arrayValue);
            if (v[4] != ".") {
                for (const std::string &x : split_on(v[4], ',')) bj["labels_per_successor"].append(Json::UInt64(u(x)));
            }
            bj["ambiguous"] = ids_json(decode_ranges(v[5]));
            bj["dropped"] = ids_json(decode_ranges(v[6]));
            bj["refused"] = Json::Value(Json::arrayValue);
            if (v[7] != ".") {
                for (const std::string &x : split_on(v[7], ';')) {
                    auto f = split_on(x, ':');
                    Json::Value rj;
                    rj["char"] = f[0];
                    rj["cause"] = f[1];
                    rj["labels"] = ids_json(decode_ranges(f[2]));
                    bj["refused"].append(rj);
                }
            }
            bes.append(bj);
        }
        j["branch_events"] = bes;
        j["branch_events_truncated"] = Json::UInt64(u(a[12]) - arm.events.size());
        // ---- label_summary, the normative pseudocode of §5.1
        const size_t ai = side == 'l' ? 0 : 1;
        std::vector<uint64_t> direct(dict.size()), reach(dict.size()), reentries(dict.size());
        std::vector<std::vector<uint64_t>> label_runs(dict.size());
        if (!annotate) {
            std::vector<uint64_t> root(arm.runs.size());
            for (size_t r = 0; r < arm.runs.size(); ++r) {
                const RunRec &x = arm.runs[r];
                root[r] = x.by_switch && x.prev ? root[*x.prev] : x.label;
            }
            for (size_t r = 0; r < arm.runs.size(); ++r) {
                const RunRec &x = arm.runs[r];
                label_runs[x.label].push_back(r);
                if (!x.by_switch && x.from == 0) direct[x.label] = std::max(direct[x.label], x.to);
                reach[root[r]] = std::max(reach[root[r]], x.to);
            }
            for (const SegRec &s : segs) {
                for (const auto &e : s.events) {
                    if (e[2] == "s") reentries[u(e[4])]++;
                }
            }
        } else {
            for (const SegRec &s : segs) {
                for (size_t k = 0; k < s.sets.size(); ++k) {
                    for (uint64_t l : s.sets[k]) reach[l] = std::max(reach[l], u(s.presence[k][2]));
                }
            }
            std::vector<std::vector<uint64_t>> alive_end(segs.size());
            for (size_t i = 0; i < segs.size(); ++i) {
                const SegRec &s = segs[i];
                std::vector<uint64_t> alive;
                if (s.parents.empty()) {
                    alive = s.entry_set;
                } else {
                    std::vector<std::vector<uint64_t>> ends;
                    for (uint64_t p : s.parents) ends.push_back(alive_end[p]);
                    alive = unite(ends);
                }
                for (size_t k = 0; k < s.sets.size() && !alive.empty(); ++k) {
                    std::vector<uint64_t> still;
                    std::set_intersection(alive.begin(), alive.end(), s.sets[k].begin(), s.sets[k].end(),
                                          std::back_inserter(still));
                    for (uint64_t l : alive) {
                        if (!std::binary_search(still.begin(), still.end(), l))
                            direct[l] = std::max(direct[l], u(s.presence[k][1]));
                    }
                    alive.swap(still);
                }
                for (uint64_t l : alive) direct[l] = std::max(direct[l], s.from + s.length);
                alive_end[i] = alive;
            }
        }
        for (size_t l = 0; l < dict.size(); ++l) {
            Json::Value sj;
            sj["direct_bp"] = Json::UInt64(direct[l]);
            sj["reach_bp"] = Json::UInt64(reach[l]);
            sj["reentries"] = Json::UInt64(reentries[l]);
            sj["runs"] = ids_json(label_runs[l]);
            (*ls)[l][ai] = sj;
        }
        (void)name;
        return j;
    }
};

// events within one segment are unordered among equal at_bp (§2.5), needed_budgets is a
// histogram, timing is not part of the result
void normalise(Json::Value *result) {
    result->removeMember("timing");
    for (const char *side : { "left", "right" }) {
        if (!(*result)["arms"].isMember(side))
            continue;
        Json::Value &arm = (*result)["arms"][side];
        for (Json::Value &seg : arm["segments"]) {
            std::vector<std::pair<uint64_t, std::string>> evs;
            for (const Json::Value &e : seg["events"]) evs.emplace_back(e["at_bp"].asUInt64(), compact(e));
            std::sort(evs.begin(), evs.end());
            Json::Value sorted(Json::arrayValue);
            for (const auto &[at, text] : evs) sorted.append(parse_json(text));
            seg["events"] = sorted;
        }
        std::vector<double> nb;
        for (const Json::Value &x : arm["needed_budgets"]) nb.push_back(x.asDouble());
        std::sort(nb.begin(), nb.end());
        arm["needed_budgets"] = Json::Value(Json::arrayValue);
        for (double x : nb) arm["needed_budgets"].append(x);
    }
}

// first difference between two JSON values, as a path
std::string first_difference(const Json::Value &a, const Json::Value &b, const std::string &path = "") {
    if (a.isObject() && b.isObject()) {
        std::set<std::string> keys;
        for (const auto &k : a.getMemberNames()) keys.insert(k);
        for (const auto &k : b.getMemberNames()) keys.insert(k);
        for (const auto &k : keys) {
            if (!a.isMember(k) || !b.isMember(k))
                return path + "." + k + (a.isMember(k) ? " (only in the graphlet's)" : " (only in full)");
            std::string d = first_difference(a[k], b[k], path + "." + k);
            if (!d.empty()) return d;
        }
        return "";
    }
    if (a.isArray() && b.isArray()) {
        if (a.size() != b.size())
            return path + ": " + std::to_string(a.size()) + " vs " + std::to_string(b.size()) + " items";
        for (Json::ArrayIndex i = 0; i < a.size(); ++i) {
            std::string d = first_difference(a[i], b[i], path + "[" + std::to_string(i) + "]");
            if (!d.empty()) return d;
        }
        return "";
    }
    return compact(a) == compact(b) ? "" : path + ": " + compact(a) + " vs " + compact(b);
}

std::string random_seq(size_t len, uint32_t seed) {
    std::mt19937 gen(seed);
    std::string s(len, 'A');
    for (char &c : s) c = "ACGT"[gen() % 4];
    return s;
}

struct Coverage {
    std::map<std::string, size_t> seen;     // record kinds and field forms met
    void note(const std::string &text) {
        for (const std::string &line : split_on(text.substr(0, text.size() - 1), '\n')) {
            seen[line.substr(0, 1)]++;
            auto f = fields(line, line[0] == 'L' || line[0] == 'X' ? 4 : line[0] == 'K' ? 9
                                  : line[0] == 'Q' ? 10 : SIZE_MAX);
            if (line[0] == 'R') seen["R:" + f[5]]++;
            if (line[0] == 'E') seen["E:" + f[2]]++;
            if (line[0] == 'G') {
                seen[std::string("G.end:") + (f[6] == "*" ? "*" : "explicit")]++;
                seen[std::string("G.entry:") + (f[4] == "*" ? "*" : f[4][0] == '!' ? "delta" : "explicit")]++;
                if (f[7] != "*") seen["G.partition:explicit"]++;
                if (f[8] == "1") seen["G.split:ambiguous"]++;
                if (f[9] != "*") seen["G.first_base"]++;
            }
            if (line[0] == 'C') seen[f.size() > 5 ? "C:sequence" : "C:derived"]++;
        }
    }
};

// One request, three ways: graphlet twice (byte-equal) and full; the graphlet's reader
// must rebuild the full result, and the body must be under half the compact full JSON (T40)
void check_request(const AnnotatedDBG &anno, Json::Value request, const std::string &what,
                   Coverage *coverage) {
    request["strategy"]["output"]["detail"] = "graphlet";
    request["strategy"]["output"]["timing"] = false;
    const Json::Value g = process_traverse_request(request, anno, "");
    EXPECT_EQ(compact(g), compact(process_traverse_request(request, anno, "")))
        << what << ": not deterministic";
    EXPECT_EQ(1, g["capabilities"]["graphlet_format"].asInt());
    EXPECT_EQ("graphlet", g["strategy"]["output"]["detail"].asString());
    request["strategy"]["output"]["detail"] = "full";
    const Json::Value full = process_traverse_request(request, anno, "");
    ASSERT_EQ(full["results"].size(), g["results"].size()) << what;
    for (Json::ArrayIndex i = 0; i < g["results"].size(); ++i) {
        const Json::Value &summary = g["results"][i];
        const std::string at = what + " seed " + std::to_string(i);
        if (summary.isMember("error")) {
            // a failed derivation keeps today's shape, no graphlet
            EXPECT_FALSE(summary.isMember("graphlet")) << at;
            EXPECT_EQ(compact(full["results"][i]), compact(summary)) << at;
            continue;
        }
        ASSERT_TRUE(summary["graphlet"].isString()) << at;
        const std::string text = summary["graphlet"].asString();
        coverage->note(text);
        EXPECT_EQ(text.size(), summary["graphlet_bytes"].asUInt64()) << at;
        Reader reader(text);
        EXPECT_EQ(reader.lines, summary["graphlet_lines"].asUInt64()) << at;
        Json::Value rebuilt = reader.to_json(summary);
        Json::Value expected = full["results"][i];
        normalise(&rebuilt);
        normalise(&expected);
        const std::string diff = first_difference(rebuilt, expected);
        EXPECT_TRUE(diff.empty()) << at << ": " << diff << "\n" << text;
        // the summary's counts and the body agree
        for (const auto &[side, arm] : reader.arms) {
            const Json::Value &counts = summary["arms"][side == 'l' ? "left" : "right"]["counts"];
            EXPECT_EQ(counts["segments"].asUInt64(), arm.segs.size()) << at;
            EXPECT_EQ(counts["runs"].asUInt64(), arm.runs.size()) << at;
        }
        EXPECT_LT(2 * text.size(), compact(full["results"][i]).size()) << at << ": the body is not under half";
    }
}

Json::Value request_for(const std::vector<std::pair<std::string, std::vector<std::string>>> &seeds,
                        const std::string &strategy) {
    Json::Value r;
    for (const auto &[sequence, labels] : seeds) {
        Json::Value s;
        s["sequence"] = sequence;
        for (const auto &l : labels) s["labels"].append(l);
        r["seeds"].append(s);
    }
    r["strategy"] = parse_json(strategy);
    return r;
}

} // namespace


TEST(Graphlet, ReaderRebuildsTheFullResult) {
    Coverage coverage;
    // dense random graphs (k = 3): splits, merges of several parents, revisits, blocked
    // successors, switches and silent switch-source ends, hairpins in the stranded modes
    std::vector<DeBruijnGraph::Mode> modes { DeBruijnGraph::BASIC };
#if ! _PROTEIN_GRAPH
    modes.push_back(DeBruijnGraph::CANONICAL);
    modes.push_back(DeBruijnGraph::PRIMARY);
#endif
    const std::vector<std::string> strategies {
        R"({"bounds": {"max_extension_bp": 8}})",
        R"({"exhaustive": true, "bounds": {"max_extension_bp": 7}})",
        R"({"direction": "left", "bounds": {"max_extension_bp": 8},
            "branching": {"min_successor_labels": 2, "max_label_branches": 1}})",
        R"({"bounds": {"max_extension_bp": 8, "max_steps": 12}})",
        R"({"bounds": {"max_extension_bp": 8}, "output": {"sequences": false, "max_branch_events": 1}})",
        R"({"bounds": {"max_extension_bp": 8}, "output": {"sequences": false, "continuation_bp": 0}})",
        R"({"bounds": {"max_extension_bp": 8}, "branching": {"hairpins": "follow", "on_reconverge": "keep"}})",
        R"({"bounds": {"max_extension_bp": 8}, "labels": {"max_labels_per_node": 1}})",
        R"({"bounds": {"max_extension_bp": 8}, "output": {"continuation_bp": 3}})",
        // the budgets of §14: a work stop (Q, resource_limit ends), a memory budget that
        // admits everything (memory_bound_soft only)
        R"({"bounds": {"max_extension_bp": 8, "max_work_units": 120}})",
        R"({"bounds": {"max_extension_bp": 8, "max_memory_mb": 1}, "output": {"sequences": false}})",
    };
    const std::vector<std::string> switching {
        R"({"bounds": {"max_extension_bp": 9}, "labels": {"extra": ["E", "F"], "loss_budget": 2,
            "change_cost": {"model": "constant", "value": 1}}, "branching": {"max_label_branches": 1}})",
        R"({"bounds": {"max_extension_bp": 9, "max_live_paths": 40}, "labels": {"extra": ["E", "F"],
            "loss_budget": 1.5, "switch_on": "any", "change_cost": {"model": "constant", "value": 0.5}},
            "branching": {"on_reconverge": "keep", "max_label_branches": "unlimited"}})",
        // one switch fits the budget, a second does not: loss_budget ends after a switch
        R"({"bounds": {"max_extension_bp": 9, "max_live_paths": 60}, "labels": {"extra": ["E", "F"],
            "loss_budget": 1, "change_cost": {"model": "constant", "value": 1}},
            "branching": {"max_label_branches": "unlimited"}})",
        // a switch above the budget: loss_budget ends with their needed budget
        R"({"bounds": {"max_extension_bp": 9}, "labels": {"extra": ["E", "F"], "loss_budget": 1,
            "change_cost": {"model": "table", "default": "forbid",
                            "entries": [["C", "E", 0.5], ["D", "F", 0.5], ["C", "F", 2.5],
                                        ["D", "E", 3]]}}})",
    };
    const std::vector<std::string> annotate {
        R"({"exhaustive": true, "labels": {"mode": "annotate"}, "bounds": {"max_extension_bp": 5}})",
        R"({"labels": {"mode": "annotate", "max_labels_per_node": 1}, "bounds": {"max_extension_bp": 5}})",
        R"({"labels": {"mode": "annotate"}, "frontier": {"on_overflow": "beam", "order": "most_supported_first"},
            "bounds": {"max_extension_bp": 8, "max_live_paths": 2}})",
        R"({"labels": {"mode": "annotate"}, "branching": {"on_reconverge": "merge"},
            "bounds": {"max_extension_bp": 8}, "output": {"sequences": false}})",
        R"({"labels": {"mode": "annotate"}, "bounds": {"max_extension_bp": 8, "max_work_units": 600}})",
    };
    for (auto mode : modes) {
        for (uint32_t s = 1; s <= 4; ++s) {
            std::vector<std::string> seqs, labels;
            for (uint32_t i = 0; i < 6; ++i) {
                seqs.push_back("AAA" + random_seq(14, s * 11 + i));
                labels.push_back(std::string(1, "CDECDF"[i]));
            }
            auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(3, seqs, labels, mode);
            const std::string where = "mode " + std::to_string(mode) + " graph " + std::to_string(s);
            for (const std::string &st : strategies)
                check_request(*anno, request_for({ { "AAA", { "C", "D", "E" } } }, st), where + " " + st, &coverage);
            for (const std::string &st : switching)
                check_request(*anno, request_for({ { "AAA", { "C", "D" } } }, st), where + " " + st, &coverage);
            for (const std::string &st : annotate)
                check_request(*anno, request_for({ { "AAA", {} } }, st), where + " " + st, &coverage);
            // a batch: a derived set, a duplicate of it, an explicit list
            check_request(*anno, request_for({ { "AAA", {} }, { "AAA", {} }, { "AAA", { "F" } } },
                                             R"({"bounds": {"max_extension_bp": 6}})"),
                          where + " batch", &coverage);
        }
    }
    // bubbles of equal length closing into Y, so that walks merge: C follows all three
    // branches (merge closures of its runs, a three-parent merge), A, B, D one each
    for (uint32_t s = 0; s < 2; ++s) {
        std::string X, P, Q, R, Y;
        for (uint32_t t = 200 + 10 * s; ; t += 5) {
            X = random_seq(30, t); P = random_seq(20, t + 1); Q = random_seq(20, t + 2);
            R = random_seq(20, t + 3); Y = random_seq(30, t + 4);
            if (P[0] != Q[0] && Q[0] != R[0] && P[0] != R[0])
                break;
        }
        auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                11, { X + P + Y, X + Q + Y, X + P + Y, X + Q + Y, X + R + Y, X + R + Y + X.substr(0, 5) },
                { "A", "B", "C", "C", "C", "D" }, modes[s % modes.size()]);
        for (const std::string &st : {
                std::string(R"({"branching": {"max_label_branches": "unlimited"}})"),
                std::string(R"({"branching": {"max_label_branches": 1, "min_successor_labels": 2}})"),
                std::string(R"({"exhaustive": true})"),
                std::string(R"({"labels": {"mode": "annotate"}, "branching": {"on_reconverge": "merge"}})") }) {
            Json::Value request = request_for({ { X, { "A", "B", "C", "D" } } }, st);
            if (st.find("annotate") != std::string::npos)
                request["seeds"][0].removeMember("labels");
            check_request(*anno, request, "bubble " + std::to_string(s) + " " + st, &coverage);
        }
    }
    // a seed present in the graph that no label carries in full: a failed derivation
    // between two good seeds; and names that need front coding with spaces and escapes
    {
        const std::string X = random_seq(30, 101), P = random_seq(30, 102), Y = random_seq(30, 103),
                          Z = random_seq(30, 104);
        auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                11, { X + P, P + Y, X + P + Z }, { "acc 1%", "acc 2\r", "éfoo" }, DeBruijnGraph::BASIC);
        // X·P·Y is in the graph but only acc 1% and acc 2 carry it, each in part: the
        // middle seed's derivation fails; the last seed names a label that is dropped
        check_request(*anno, request_for({ { X.substr(10) + P, {} },
                                           { X.substr(25) + P + Y.substr(0, 5), {} },
                                           { P + Y.substr(0, 10), { "acc 2\r", "acc 1%" } } },
                                         R"({"bounds": {"max_extension_bp": 40}})"),
                      "names", &coverage);
        check_request(*anno, request_for({ { P, {} } }, R"({"labels": {"mode": "annotate"},
                                                           "bounds": {"max_extension_bp": 40}})"),
                      "names annotate", &coverage);
    }
    // record coordinates (opt-in, §18): the envelope carries the block (or null with its
    // reason) in every detail, and a cut list adds its K record to the body — the reader
    // rebuilds the full result from the two, the block copied from the summary
    {
        const std::string S = random_seq(30, 131), P = random_seq(20, 132), Q = random_seq(20, 133);
        auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                11, { S + P + S + Q, S + P }, { "C", "D" }, DeBruijnGraph::BASIC, true, { 0, 0 });
        for (const std::string cap : { "16", "1", "\"unlimited\"" }) {
            check_request(*anno, request_for({ { S, { "C", "D" } } },
                                             R"({"support": "trace", "branching": {"on_reconverge": "keep",
                                                 "max_label_branches": 2}, "output": {"coordinates": true,
                                                 "max_coordinate_occurrences": )" + cap + "}}"),
                          "coordinates, cap " + cap, &coverage);
        }
        check_request(*anno, request_for({ { S, {} } }, R"({"output": {"coordinates": true}})"),
                      "coordinates under support kmer", &coverage);
    }
    // not vacuous: every record kind and every end token the walker can produce today
    for (const char *kind : { "H", "S", "L", "O", "Q", "K", "A", "B", "V", "G", "P", "E", "T", "C", "R", "Z",
                              "R:Y",
                              "E:s", "E:b", "E:h", "E:v", "R:Lw", "R:m", "R:D", "R:L", "R:X", "R:R",
                              "R:Rm", "R:B", "R:Ls", "G.end:*", "G.entry:*",
                              "G.entry:delta", "G.entry:explicit", "G.partition:explicit",
                              "G.split:ambiguous", "G.first_base", "C:derived", "C:sequence" }) {
        EXPECT_GT(coverage.seen[kind], 0u) << kind << " never written";
    }
}


/*
 * The index identity (DESIGN-traverse-graphlet.md §3.1). index_meta_fp is a negative
 * check only: two annotations of the same records with the same column names and counts
 * but swapped memberships get the SAME meta fingerprint. index_fp, the digest of the
 * bundle's manifest, tells them apart because their files differ, and a manifest that
 * does not describe the loaded files is refused.
 */


// Pass 5, W6: the server writes each seed's result as text once it is built and assembles the
// response from the texts; the decompressed bytes are those of the whole tree written at once —
// every detail, failed seeds beside walked ones, with and without an attempt's usage (its
// clock-dependent values are the same in both writings of one response)
TEST(GraphletServer, AssembledResponseIsByteIdentical) {
    std::vector<std::string> seqs, labels;
    std::mt19937 gen(11);
    for (size_t i = 0; i < 6; ++i) {
        std::string s(40, 'A');
        for (char &c : s) c = "ACGT"[gen() % 4];
        seqs.push_back(s);
        labels.push_back("L" + std::to_string(i % 3));
    }
    // P+Q and Q+R in two labels: a path P..R through Q that no label carries whole
    const std::string P = "GATTACAGGCATTAC", Q = "CCGTTAGCAT", R = "TTGACCAGTAGGCTA";
    seqs.push_back(P + Q);
    labels.push_back("X");
    seqs.push_back(Q + R);
    labels.push_back("Y");
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(5, seqs, labels);
    size_t checked = 0;
    for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
        for (bool managed : { false, true }) {
            Json::Value r;
            r["seeds"][0]["sequence"] = seqs[0].substr(0, 12);
            r["seeds"][1]["sequence"] = seqs[1].substr(3, 14);
            r["seeds"][1]["labels"].append("L1");
            // no carrier of every k-mer: a failed derivation beside the walked seeds
            r["seeds"][2]["sequence"] = P.substr(10) + Q + R.substr(0, 5);
            r["strategy"] = parse_json(R"({"bounds": {"max_extension_bp": 6}})");
            r["strategy"]["output"]["detail"] = detail;
            r["strategy"]["output"]["timing"] = false;
            if (managed)
                r["attempt_id"] = "asm-1";
            // the tree, and the texts of one run: equal when written whole
            const Json::Value tree = process_traverse_request(r, *anno, "");
            ResultTexts texts;
            Json::Value envelope = process_traverse_request(r, *anno, "", {}, nullptr, nullptr,
                                                            &texts);
            ASSERT_TRUE(texts.active);
            ASSERT_EQ(3u, texts.texts.size());
            EXPECT_TRUE(parse_json(texts.texts[2]).isMember("error"));
            EXPECT_FALSE(envelope.isMember("results"));
            int checks = 0;
            const std::string assembled = assemble_traverse_response(envelope, texts.texts,
                                                                     [&]() { checks++; });
            // the same response: the texts are the results' compact texts
            Json::Value whole = envelope;
            for (const std::string &t : texts.texts) {
                whole["results"].append(parse_json(t));
            }
            EXPECT_EQ(json_text(whole, true), assembled) << detail;
            if (!managed) {
                // no clock in it: the response of the tree, byte for byte
                EXPECT_EQ(json_text(tree, true), assembled) << detail;
            } else {
                EXPECT_TRUE(envelope.isMember("usage"));
            }
            checked++;
        }
    }
    EXPECT_EQ(8u, checked);
    // members on one side of "results" only
    Json::Value only_before;
    only_before["a"] = 1;
    EXPECT_EQ("{\"a\":1,\"results\":[{}]}", assemble_traverse_response(only_before, { "{}" }));
    Json::Value only_after;
    only_after["z"] = "x";
    EXPECT_EQ("{\"results\":[1,2],\"z\":\"x\"}", assemble_traverse_response(only_after, { "1", "2" }));
    EXPECT_EQ("{\"results\":[]}", assemble_traverse_response(Json::Value(Json::objectValue), {}));
    Json::Value bad;
    bad["results"] = 1;
    EXPECT_THROW(assemble_traverse_response(bad, {}), std::logic_error);
}

TEST(Graphlet, IndexIdentity) {
    const std::string S1 = random_seq(60, 301), S2 = random_seq(60, 302);
    auto ab = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, { S1, S2 }, { "A", "B" }, DeBruijnGraph::BASIC);
    auto swapped = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, { S2, S1 }, { "A", "B" }, DeBruijnGraph::BASIC);
    const std::string meta = index_meta_fingerprint(LabelOracle(*ab));
    EXPECT_EQ(16u, meta.size());
    // not identity: the swap keeps k, rows, column names and their order
    EXPECT_EQ(meta, index_meta_fingerprint(LabelOracle(*swapped)));
    // ... but a different column name changes it (a mismatch proves different indexes)
    auto renamed = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, { S1, S2 }, { "A", "C" }, DeBruijnGraph::BASIC);
    EXPECT_NE(meta, index_meta_fingerprint(LabelOracle(*renamed)));

    // the responses state it; without a manifest index_fp is null (joins unverifiable)
    Json::Value request = request_for({ { S1.substr(0, 20), { "A" } } }, R"({"output": {"detail": "graphlet"}})");
    Json::Value out = process_traverse_request(request, *ab, "");
    EXPECT_TRUE(out["capabilities"]["index_fp"].isNull());
    EXPECT_TRUE(out["capabilities"]["index_ns"].isNull());
    EXPECT_EQ(meta, out["capabilities"]["index_meta_fp"].asString());
    auto header = [](const Json::Value &out) {
        const std::string h = split_on(out["results"][0]["graphlet"].asString(), '\n')[0];
        return h.substr(h.find(" walk ") + 1);
    };
    EXPECT_EQ("walk * * " + meta, header(out));
    IndexIdentity named { "refseq.v2", std::string(64, 'a'), meta };
    out = process_traverse_request(request, *ab, "", {}, &named);
    EXPECT_EQ("refseq.v2", out["capabilities"]["index_ns"].asString());
    EXPECT_EQ(std::string(64, 'a'), out["capabilities"]["index_fp"].asString());
    EXPECT_EQ("walk refseq.v2 " + std::string(64, 'a') + " " + meta, header(out));
    EXPECT_TRUE(valid_index_name("mini_refseq-1.0"));
    EXPECT_FALSE(valid_index_name("a b"));
    EXPECT_FALSE(valid_index_name(""));
    EXPECT_FALSE(valid_index_name("*"));

    // manifests: two bundles that differ in one file's content get different fingerprints
    const std::filesystem::path dir = "temp_graphlet_manifest_" + std::to_string(getpid());
    std::filesystem::create_directories(dir);
    auto write = [&](const std::string &name, const std::string &content) {
        std::ofstream(dir / name, std::ios::binary) << content;
        return (dir / name).string();
    };
    const std::string graph = write("g.dbg", "graph bytes"), anno = write("a.annodbg", "annotation A");
    auto manifest = [&](const std::string &name, const std::string &anno_content,
                        const std::string &extra = "") {
        Json::Value m;
        Json::Value g, a;
        g["path"] = "g.dbg"; g["size"] = Json::UInt64(11); g["sha256"] = sha256_hex("graph bytes");
        a["path"] = "a.annodbg"; a["size"] = Json::UInt64(anno_content.size());
        a["sha256"] = sha256_hex(anno_content);
        m["files"].append(g);
        m["files"].append(a);
        if (!extra.empty()) m["index_fp"] = extra;
        return write(name, compact(m));
    };
    const std::string fp = index_manifest_fingerprint(manifest("m1.json", "annotation A"), { graph, anno });
    // the canonical file list, in byte order of path
    EXPECT_EQ(sha256_hex("a.annodbg\t12\t" + sha256_hex("annotation A") + "\n"
                         "g.dbg\t11\t" + sha256_hex("graph bytes") + "\n"), fp);
    // (another content of the same size: what sizes cannot tell apart, the digests do)
    EXPECT_NE(fp, index_manifest_fingerprint(manifest("m2.json", "annotation B"), { graph, anno }));
    // a stated index_fp must be the computed one
    EXPECT_EQ(fp, index_manifest_fingerprint(manifest("m3.json", "annotation A", fp), { graph, anno }));
    EXPECT_THROW(index_manifest_fingerprint(manifest("m4.json", "annotation A", std::string(64, '0')),
                                            { graph, anno }),
                 std::runtime_error);
    // the loaded annotation is not the one the manifest lists (another size)
    EXPECT_THROW(index_manifest_fingerprint(manifest("m5.json", "annotation AB"), { graph, anno }),
                 std::runtime_error);
    // malformed: no files, a bad digest, a path listed twice
    EXPECT_THROW(index_manifest_fingerprint(write("m6.json", "{\"files\": []}"), {}), std::runtime_error);
    EXPECT_THROW(index_manifest_fingerprint(
            write("m7.json", R"({"files": [{"path": "x", "size": 1, "sha256": "ABC"}]})"), {}),
            std::runtime_error);
    const std::string digest = sha256_hex("x");
    EXPECT_THROW(index_manifest_fingerprint(
            write("m8.json", R"({"files": [{"path": "x", "size": 1, "sha256": ")" + digest
                                + R"("}, {"path": "x", "size": 1, "sha256": ")" + digest + R"("}]})"), {}),
            std::runtime_error);
    EXPECT_THROW(index_manifest_fingerprint((dir / "missing.json").string(), {}), std::runtime_error);

    // Review of pass 5: a manifest describes one graph with one annotation — one that also
    // lists another annotation (written for a directory holding both) or another graph would
    // lend one fingerprint to several indexes, and is refused
    auto slurp = [](const std::filesystem::path &path) {
        std::ifstream in(path, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(in), {});
    };
    auto listing = [&](const std::string &name, const std::vector<std::string> &extra) {
        Json::Value m = parse_json(slurp(dir / "m1.json"));
        for (const std::string &path : extra) {
            Json::Value e;
            e["path"] = path;
            e["size"] = Json::UInt64(7);
            e["sha256"] = sha256_hex(path);
            m["files"].append(e);
        }
        return write(name, compact(m));
    };
    // const char *: a std::string reference would bind to a temporary per element
    // (GCC's -Wrange-loop-construct)
    for (const char *other : { "a2.annodbg", "other.dbg", "g2.orhashdbg" }) {
        try {
            index_manifest_fingerprint(listing("m9.json", { other }), { graph, anno });
            ADD_FAILURE() << "a manifest listing " << other << " was accepted";
        } catch (const std::runtime_error &e) {
            EXPECT_NE(std::string::npos, std::string(e.what()).find("which this index does not "
                                                                    "load")) << e.what();
        }
    }
    // other files (metadata the tool's --extra adds) are part of the identity, not refused
    EXPECT_NO_THROW(index_manifest_fingerprint(listing("m10.json", { "README" }), { graph, anno }));
    // ... and the sidecars the loader reads beside the graph and the annotation (the loader
    // dependency inventory) are checked as the two are: a manifest whose sidecars are another
    // build's, or that does not cover one, does not lend its identity
    EXPECT_EQ((std::vector<std::string> { graph, anno }), index_bundle_files(graph, anno));
    const std::string mask = write("g.edgemask", "mask"),
                      bloom = write("g.bloom", "bloom filter");
    write("a.annodbg.coords", "coords");   // not loaded: never in the inventory
    write("g.dbg.anchors", "anchors");      // a leftover beside a non-row_diff annotation
    const std::string coord_anno = write("x.row_diff_brwt_coord.annodbg", "coordinate annotation"),
                      seqs = write("x.seqs", "headers");
    const std::string rd_anno = write("r.row_diff.annodbg", "row diff annotation");
    const std::string anchors = (dir / "g.dbg.anchors").string(),
                      fork_succ = (dir / "g.dbg.rd_succ").string();
    EXPECT_EQ((std::vector<std::string> { graph, mask, bloom, anno }),
              index_bundle_files(graph, anno));
    EXPECT_EQ((std::vector<std::string> { graph, mask, bloom, coord_anno, seqs }),
              index_bundle_files(graph, coord_anno));
    // --no-coord-mapping: the headers are not loaded
    EXPECT_EQ((std::vector<std::string> { graph, mask, bloom, coord_anno }),
              index_bundle_files(graph, coord_anno, false));
    // row_diff (column): the anchors and fork successors beside the graph, required (listed
    // whether they exist or not: the loader fails without them)
    EXPECT_EQ((std::vector<std::string> { graph, mask, bloom, rd_anno, anchors, fork_succ }),
              index_bundle_files(graph, rd_anno));
    // the Bloom filter is read only with the mask
    std::filesystem::remove(dir / "g.edgemask");
    EXPECT_EQ((std::vector<std::string> { graph, anno }), index_bundle_files(graph, anno));
    write("g.edgemask", "mask");
    // a manifest that covers neither the mask nor the Bloom filter: refused, naming the file
    try {
        index_manifest_fingerprint((dir / "m1.json").string(), index_bundle_files(graph, anno));
        ADD_FAILURE() << "a manifest without the mask was accepted";
    } catch (const std::runtime_error &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("does not cover the file "
                                                                + mask)) << e.what();
    }
    Json::Value full = parse_json(slurp(dir / "m1.json"));
    for (const auto &[name, content] : std::vector<std::pair<std::string, std::string>> {
             { "g.edgemask", "mask" } }) {
        Json::Value e;
        e["path"] = name;
        e["size"] = Json::UInt64(content.size());
        e["sha256"] = sha256_hex(content);
        full["files"].append(e);
    }
    const std::string no_bloom = write("m11.json", compact(full));
    try {
        index_manifest_fingerprint(no_bloom, index_bundle_files(graph, anno));
        ADD_FAILURE() << "a manifest without the Bloom filter was accepted";
    } catch (const std::runtime_error &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("does not cover the file "
                                                                + bloom)) << e.what();
    }
    Json::Value bloom_entry;
    bloom_entry["path"] = "g.bloom";
    bloom_entry["size"] = Json::UInt64(12);
    bloom_entry["sha256"] = sha256_hex("bloom filter");
    full["files"].append(bloom_entry);
    const std::string complete = write("m12.json", compact(full));
    EXPECT_NO_THROW(index_manifest_fingerprint(complete, index_bundle_files(graph, anno)));
    write("g.bloom", "a Bloom filter of another build");
    try {
        index_manifest_fingerprint(complete, index_bundle_files(graph, anno));
        ADD_FAILURE() << "another build's Bloom filter was accepted";
    } catch (const std::runtime_error &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("is not listed with that size"))
                << e.what();
    }
    write("g.bloom", "bloom filter");

    // Review of the pass-5 fixes: a manifest written for a directory of bundles that share
    // base names (a/x.seqs, b/x.seqs) matched every one of them, and two indexes stated one
    // index_fp; the base names of a manifest's entries must be distinct
    Json::Value shared = parse_json(slurp(dir / "m12.json"));
    for (const char *path : { "a/x.seqs", "b/x.seqs" }) {
        Json::Value e;
        e["path"] = path;
        e["size"] = Json::UInt64(7);
        e["sha256"] = sha256_hex(path);
        shared["files"].append(e);
    }
    try {
        index_manifest_fingerprint(write("m13.json", compact(shared)),
                                   index_bundle_files(graph, anno));
        ADD_FAILURE() << "a manifest with two entries named x.seqs was accepted";
    } catch (const std::runtime_error &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("share the base name x.seqs"))
                << e.what();
    }
    // ... also when one of the two is the bundle's own file (sub/g.dbg beside g.dbg)
    try {
        index_manifest_fingerprint(listing("m14.json", { "sub/g.dbg" }), { graph, anno });
        ADD_FAILURE() << "a manifest with two entries named g.dbg was accepted";
    } catch (const std::runtime_error &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("share the base name g.dbg"))
                << e.what();
    }

    // ... and a manifest lists exactly the optional inventory files the pair loads: one that
    // lists a .seqs, a mask or a Bloom filter the server does not load (missing beside the
    // listed spelling, the mask unread, --no-coord-mapping) would state the identity of an
    // index that does load it
    Json::Value with_seqs = parse_json(slurp(dir / "m12.json"));
    with_seqs["files"][1]["path"] = "x.row_diff_brwt_coord.annodbg";
    with_seqs["files"][1]["size"] = Json::UInt64(21);
    with_seqs["files"][1]["sha256"] = sha256_hex("coordinate annotation");
    Json::Value seqs_entry;
    seqs_entry["path"] = "x.seqs";
    seqs_entry["size"] = Json::UInt64(7);
    seqs_entry["sha256"] = sha256_hex("headers");
    with_seqs["files"].append(seqs_entry);
    const std::string coord_manifest = write("m15.json", compact(with_seqs));
    EXPECT_NO_THROW(index_manifest_fingerprint(coord_manifest, index_bundle_files(graph, coord_anno),
                                               nullptr,
                                               index_unloaded_optional_files(graph, coord_anno)));
    EXPECT_TRUE(index_unloaded_optional_files(graph, coord_anno).empty());
    EXPECT_EQ((std::vector<std::string> { seqs }),
              index_unloaded_optional_files(graph, coord_anno, false));
    try {
        index_manifest_fingerprint(coord_manifest, index_bundle_files(graph, coord_anno, false),
                                   nullptr, index_unloaded_optional_files(graph, coord_anno, false));
        ADD_FAILURE() << "a manifest listing the unloaded .seqs was accepted";
    } catch (const std::runtime_error &e) {
        EXPECT_NE(std::string::npos, std::string(e.what()).find("which the server does not load"))
                << e.what();
    }
    // the mask unread: neither it nor the Bloom filter is loaded, and listing either is refused
    std::filesystem::remove(dir / "g.edgemask");
    EXPECT_EQ((std::vector<std::string> { mask, bloom }),
              index_unloaded_optional_files(graph, anno));
    EXPECT_THROW(index_manifest_fingerprint(complete, index_bundle_files(graph, anno), nullptr,
                                            index_unloaded_optional_files(graph, anno)),
                 std::runtime_error);
    std::filesystem::remove_all(dir);
}


// GPT review, finding 1: a label name that is not UTF-8 is never replaced. The replaced
// name "bad\uFFFD" is another label's VALID name here, and a continuation that named it
// went on under that label along other bases. The seed is refused per seed instead, in
// both modes, in the shape of a failed derivation (cause unrepresentable_label_name, the
// knob that avoids it, the column named in the effect, never the name's bytes); the other
// seeds of the request are traversed, and detail full and graphlet agree.
TEST(Graphlet, NamesThatAreNotUtf8AreRefusedPerSeed) {
    const std::string S = random_seq(40, 411), T = random_seq(40, 412), U = random_seq(40, 413);
    const std::string invalid = "bad\xFF", replaced = "bad\xEF\xBF\xBD";
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, { S + T, S + U, U }, { invalid, replaced, "ok" }, DeBruijnGraph::BASIC);
    auto column_of = [&](const std::string &name) {
        return anno->get_annotator().get_label_encoder().encode(name);
    };
    const std::string seed = S.substr(0, 20);
    struct Case {
        std::string what;
        std::vector<std::pair<std::string, std::vector<std::string>>> seeds;
        std::string strategy;
        std::string knob;            // of the refused seed 0
        Json::Value limit;
    };
    const std::vector<Case> cases = {
        { "constrain, derived", { { seed, {} }, { U.substr(0, 20), {} } },
          R"({"direction": "right", "bounds": {"max_extension_bp": 30}})",
          "seeds[].labels", Json::Value(Json::UInt64(2)) },
        { "constrain, named", { { seed, { invalid } }, { seed, { replaced } } },
          R"({"direction": "right", "bounds": {"max_extension_bp": 30}})",
          "seeds[].labels", Json::Value(Json::UInt64(1)) },
        { "annotate", { { seed, {} }, { U.substr(0, 20), {} } },
          R"({"direction": "right", "labels": {"mode": "annotate"}, "bounds": {"max_extension_bp": 30}})",
          "labels.mode", Json::Value("annotate") },
    };
    for (const Case &c : cases) {
        Json::Value request = request_for(c.seeds, c.strategy);
        std::vector<Json::Value> outs;
        for (const char *detail : { "full", "graphlet" }) {
            request["strategy"]["output"]["detail"] = detail;
            request["strategy"]["output"]["timing"] = false;
            outs.push_back(process_traverse_request(request, *anno, ""));
        }
        for (const Json::Value &out : outs) {
            ASSERT_EQ(2u, out["results"].size()) << c.what;
            const Json::Value &refused = out["results"][0];
            EXPECT_EQ("failed", refused["outcome"]["walks"].asString()) << c.what;
            EXPECT_FALSE(refused.isMember("arms")) << c.what;
            EXPECT_FALSE(refused.isMember("graphlet")) << c.what;
            ASSERT_EQ(1u, refused["limitations"].size()) << c.what;
            const Json::Value &d = refused["limitations"][0];
            EXPECT_EQ("derivation", d["kind"].asString()) << c.what;
            EXPECT_EQ("unrepresentable_label_name", d["cause"].asString()) << c.what;
            EXPECT_EQ(c.knob, d["knob"].asString()) << c.what;
            EXPECT_EQ(c.limit, d["limit"]) << c.what;
            EXPECT_EQ(1u, d["observed"].asUInt64()) << c.what;
            const std::string column = "column " + std::to_string(column_of(invalid));
            EXPECT_NE(std::string::npos, d["effect"].asString().find(column)) << c.what;
            EXPECT_NE(std::string::npos, refused["error"].asString().find(column)) << c.what;
            // no name at all, so neither its bytes nor a replacement
            EXPECT_EQ(std::string::npos, compact(refused).find("bad")) << c.what;
            EXPECT_TRUE(valid_utf8(compact(out))) << c.what;
            // the other seed is traversed, its names valid
            const Json::Value &other = out["results"][1];
            EXPECT_FALSE(other.isMember("error")) << c.what << ": " << other["error"];
            for (const Json::Value &l : other["label_dict"]) {
                EXPECT_TRUE(valid_utf8(l["name"].asString())) << c.what;
            }
        }
        // detail full and graphlet carry the same refusal
        EXPECT_EQ(compact(outs[0]["results"][0]), compact(outs[1]["results"][0])) << c.what;
    }
    // the replaced name resolves to its own label and walks its own bases: S then U
    Json::Value request = request_for({ { seed, { replaced } } },
                                      R"({"direction": "right", "bounds": {"max_extension_bp": 30}})");
    request["strategy"]["output"]["detail"] = "full";
    const Json::Value out = process_traverse_request(request, *anno, "");
    const Json::Value &arm = out["results"][0]["arms"]["right"];
    ASSERT_EQ(1u, arm["paths"].size());
    std::string spelled;
    for (const Json::Value &sid : arm["paths"][0]["segments"]) {
        spelled += arm["segments"][sid.asUInt()]["sequence"].asString();
    }
    EXPECT_EQ((S + U).substr(20, 30), spelled);
}


/*
 * Stage 2 of DESIGN-traverse-graphlet.md §14.1: every head is planned, admitted against
 * the request's budgets and then committed; a head that is not admitted is censored like
 * a head beyond a cap (resource_limit), and the result must stay one the writer accepts,
 * the reader rebuilds, and whose walks up to complete_to_bp are exactly the unbudgeted
 * walk's.
 */
namespace {

GraphletContext context_of(const LabelOracle &oracle) {
    GraphletContext ctx;
    ctx.k = oracle.get_k();
    ctx.regime = to_string(oracle.regime());
    ctx.alphabet = oracle.graph().alphabet();
    ctx.identity.meta_fp = index_meta_fingerprint(oracle);
    return ctx;
}

// One result serialised as the server does (full; the graphlet summary and its body) and
// read back: the writer must accept it (its checks reject inconsistent runs, splits,
// events and paths) and the reader must rebuild the full result. Returns the body.
std::string check_serialised(const SeedResult &r, const Seed &seed, const Strategy &st,
                             const LabelOracle &oracle, const std::string &what) {
    const Json::Value full = seed_result_to_json(r, st, "full", false);
    const Json::Value summary = seed_result_to_json(r, st, "graphlet", false);
    std::string text;
    try {
        text = graphlet_text(r, seed, st, context_of(oracle), summary);
    } catch (const std::exception &e) {
        ADD_FAILURE() << what << ": the writer refused the result: " << e.what();
        return "";
    }
    Reader reader(text);
    Json::Value rebuilt = reader.to_json(summary);
    Json::Value expected = full;
    normalise(&rebuilt);
    normalise(&expected);
    const std::string diff = first_difference(rebuilt, expected);
    EXPECT_TRUE(diff.empty()) << what << ": " << diff << "\n" << text;
    return text;
}

std::string walked(const ArmResult &arm, const Segment &s) {
    std::string w = s.sequence;
    if (arm.arm == Arm::LEFT)
        std::reverse(w.begin(), w.end());
    return w;
}

// Every walk the arm holds of at most |c| bases, along ANY route of the segment DAG (a
// merge's routes included): "P<its first c bases>" for a walk reaching c, "E<bases>|<end>"
// for a leaf before c with its end. Two arms that agree on this agree on every walk of at
// most c bases and on why each shorter one ended.
std::set<std::string> walks_upto(const ArmResult &arm, uint64_t c) {
    std::set<std::string> out;
    if (arm.segments.empty())
        return out;
    std::vector<int64_t> leaf_of(arm.segments.size(), -1);
    for (size_t i = 0; i < arm.paths.size(); ++i) {
        leaf_of[arm.paths[i].leaf] = i;
    }
    std::function<void(size_t, const std::string&)> visit = [&](size_t s, const std::string &pre) {
        const Segment &seg = arm.segments[s];
        const std::string w = pre + walked(arm, seg);
        if (seg.from_bp + seg.length_bp >= c) {
            out.insert("P" + w.substr(0, c));
            return;
        }
        if (seg.children.empty()) {
            std::string end = "E" + w + "|";
            if (leaf_of[s] >= 0) {
                const PathResult &p = arm.paths[leaf_of[s]];
                end += p.path_reason ? to_string(*p.path_reason) : "semantic";
                for (uint32_t x : p.end_reasons) end += "," + std::to_string(x);
            }
            out.insert(end);
            return;
        }
        for (size_t child : seg.children) {
            visit(child, w);
        }
    };
    visit(0, "");
    return out;
}

// the label ends the growth bins count are the runs that ended (§7.2)
void check_bins(const ArmResult &arm, const std::string &what) {
    std::array<uint32_t, kNumEndReasons> bins {}, runs {};
    for (const GrowthBin &b : arm.growth) {
        for (size_t r = 0; r < kNumEndReasons; ++r) bins[r] += b.label_ends[r];
    }
    for (const LabelRun &run : arm.runs) {
        if (run.ended) runs[static_cast<size_t>(run.end_reason)]++;
    }
    EXPECT_EQ(runs, bins) << what;
}

struct BudgetCase {
    std::string name;
    std::shared_ptr<graph::AnnotatedDBG> anno;
    Seed seed;
    Strategy st;
    LabelChangeCost cost = LabelChangeCost::forbid();
};

Seed seed_of(const std::string &sequence, std::vector<std::string> labels = {}) {
    Seed s;
    s.sequence = sequence;
    s.labels = std::move(labels);
    return s;
}

Strategy strategy_of(const std::string &json, LabelChangeCost *cost = nullptr) {
    Json::Value r;
    r["seeds"][0]["sequence"] = "ACGT";
    r["strategy"] = parse_json(json);
    TraverseRequest req = parse_traverse_request(r);
    if (cost) {
        *cost = req.cost.model == CostSpec::CONSTANT ? LabelChangeCost::constant(req.cost.value)
                                                     : LabelChangeCost::forbid();
    }
    return req.strategy;
}

// the cases of the denial sweep: dense random graphs (splits, merges of several parents,
// revisits, blocked successors, hairpins in the stranded modes) under constrain, switching
// and annotate strategies; a switch chain (a silent switch-source end at an ambiguous split
// whose switches sit on the children, the review's case); three bubbles closing into one
// node (a three-parent merge, merge closures)
std::vector<BudgetCase> budget_cases() {
    std::vector<BudgetCase> cases;
    std::vector<DeBruijnGraph::Mode> modes { DeBruijnGraph::BASIC };
#if ! _PROTEIN_GRAPH
    modes.push_back(DeBruijnGraph::CANONICAL);
    modes.push_back(DeBruijnGraph::PRIMARY);
#endif
    for (auto mode : modes) {
        for (uint32_t s = 1; s <= 2; ++s) {
            std::vector<std::string> seqs, labels;
            for (uint32_t i = 0; i < 6; ++i) {
                seqs.push_back("AAA" + random_seq(14, s * 11 + i));
                labels.push_back(std::string(1, "CDECDF"[i]));
            }
            std::shared_ptr<graph::AnnotatedDBG> anno
                = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(3, seqs, labels, mode);
            const std::string where = "mode " + std::to_string(mode) + " graph " + std::to_string(s);
            LabelChangeCost constant;
            cases.push_back({ where + " merge", anno, seed_of("AAA", { "C", "D", "E" }),
                              strategy_of(R"({"bounds": {"max_extension_bp": 7}})") });
            cases.push_back({ where + " exhaustive", anno, seed_of("AAA", { "C", "D", "E" }),
                              strategy_of(R"({"exhaustive": true, "bounds": {"max_extension_bp": 6}})") });
            Strategy sw = strategy_of(R"({"bounds": {"max_extension_bp": 8}, "labels": {"extra": ["E", "F"],
                "loss_budget": 2, "change_cost": {"model": "constant", "value": 1}},
                "branching": {"max_label_branches": 1}})", &constant);
            cases.push_back({ where + " switching", anno, seed_of("AAA", { "C", "D" }), sw, constant });
            cases.push_back({ where + " annotate", anno, seed_of("AAA"),
                              strategy_of(R"({"exhaustive": true, "labels": {"mode": "annotate"},
                                              "bounds": {"max_extension_bp": 5}})") });
            cases.push_back({ where + " annotate merge", anno, seed_of("AAA"),
                              strategy_of(R"({"labels": {"mode": "annotate"}, "branching": {"on_reconverge": "merge"},
                                              "bounds": {"max_extension_bp": 6}})") });
        }
    }
    {
        // A -> B -> {C, D}: A carries S T1, B the last 20 bp of T1 and T2, C and D the
        // last 20 bp of T2 and then T3 / T4
        std::string S, T1, T2, T3, T4;
        for (uint32_t t = 500; ; t += 5) {
            S = random_seq(30, t); T1 = random_seq(30, t + 1); T2 = random_seq(30, t + 2);
            T3 = random_seq(30, t + 3); T4 = random_seq(30, t + 4);
            if (T3[0] != T4[0])
                break;
        }
        std::shared_ptr<graph::AnnotatedDBG> anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                11, { S + T1, T1.substr(10) + T2, T2.substr(10) + T3, T2.substr(10) + T4 },
                { "A", "B", "C", "D" }, DeBruijnGraph::BASIC);
        LabelChangeCost constant;
        Strategy st = strategy_of(R"({"direction": "right", "labels": {"extra": ["B", "C", "D"],
            "change_cost": {"model": "constant", "value": 1}, "loss_budget": 3},
            "branching": {"max_label_branches": 1}, "bounds": {"max_extension_bp": 80}})", &constant);
        cases.push_back({ "switch chain", anno, seed_of(S, { "A" }), st, constant });
    }
    {
        std::string X, P, Q, R, Y;
        for (uint32_t t = 200; ; t += 5) {
            X = random_seq(30, t); P = random_seq(20, t + 1); Q = random_seq(20, t + 2);
            R = random_seq(20, t + 3); Y = random_seq(30, t + 4);
            if (P[0] != Q[0] && Q[0] != R[0] && P[0] != R[0])
                break;
        }
        std::shared_ptr<graph::AnnotatedDBG> anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                11, { X + P + Y, X + Q + Y, X + P + Y, X + Q + Y, X + R + Y, X + R + Y },
                { "A", "B", "C", "C", "C", "D" }, DeBruijnGraph::BASIC);
        cases.push_back({ "bubbles merge", anno, seed_of(X, { "A", "B", "C", "D" }),
                          strategy_of(R"({"branching": {"max_label_branches": "unlimited"},
                                          "bounds": {"max_extension_bp": 60}})") });
        cases.push_back({ "bubbles keep", anno, seed_of(X, { "A", "B", "C", "D" }),
                          strategy_of(R"({"exhaustive": true, "bounds": {"max_extension_bp": 60}})") });
    }
    return cases;
}

SeedResult run_case(const BudgetCase &c, const Strategy &st, const WalkerHooks *hooks = nullptr) {
    LabelOracle oracle(*c.anno);
    return traverse_seed(oracle, c.seed, st, c.cost, "", hooks);
}

bool ended_by_time(const ArmResult &arm) {
    if (arm.cap_trigger && arm.cap_trigger->reason == EndReason::TIME_BUDGET)
        return true;
    return std::any_of(arm.paths.begin(), arm.paths.end(), [](const PathResult &p) {
        return p.path_reason == EndReason::TIME_BUDGET;
    });
}

} // namespace


// GPT review (also requested): T37 compares two RUNS of one request (detail full, detail
// graphlet), which a time-budget stop cannot be held to -- where the clock stops the walk
// differs between runs -- so time_budget stops are left out there (and in the real
// suite's time_tight cells). Where the comparison is sound, ONE captured time-stopped
// SeedResult serialised both ways, it holds: the reader rebuilds detail full from the
// body of the same result (check_serialised), the stop stated as a walk_domain naming
// bounds.time_budget_ms.
TEST(Graphlet, TimeStopSerialisesBothWaysFromOneResult) {
    std::vector<BudgetCase> cases = budget_cases();
    {
        // wide enough that a level takes measurable time: stops fall mid-walk too
        std::vector<std::string> seqs, labels;
        for (uint32_t i = 0; i < 24; ++i) {
            seqs.push_back("ACGTACG" + random_seq(400, 900 + i));
            labels.push_back("L" + std::to_string(i % 6));
        }
        std::shared_ptr<graph::AnnotatedDBG> anno
                = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                        5, seqs, labels, DeBruijnGraph::BASIC);
        cases.push_back({ "wide annotate", anno, seed_of("ACGTACG"),
                          strategy_of(R"({"exhaustive": true, "labels": {"mode": "annotate"},
                                          "bounds": {"max_extension_bp": 60}})") });
        cases.push_back({ "wide constrain", anno, seed_of("ACGTACG"),
                          strategy_of(R"({"bounds": {"max_extension_bp": 60}})") });
    }
    size_t stopped = 0, mid_walk = 0;
    for (const BudgetCase &c : cases) {
        LabelOracle oracle(*c.anno);
        for (double budget : { 1e-6, 0.05, 0.5, 2.0 }) {
            Strategy st = c.st;
            st.time_budget_ms = budget;
            SeedResult r;
            try {
                r = run_case(c, st);
            } catch (const SeedDerivationError &e) {
                // the clock ran out while a permitted set was derived: a failed seed, no
                // SeedResult to serialise
                EXPECT_EQ(SeedDerivationError::TIME_BUDGET, e.cause()) << c.name;
                continue;
            }
            bool time = false;
            uint64_t depth = 0;
            for (const ArmResult &a : r.arms) {
                if (a.requested && ended_by_time(a)) {
                    time = true;
                    depth = std::max(depth, a.complete_to_bp);
                }
            }
            if (!time)
                continue;
            stopped++;
            mid_walk += depth > 0;
            const std::string what = c.name + " time_budget_ms " + std::to_string(budget);
            check_serialised(r, c.seed, st, oracle, what);
            const Json::Value j = seed_result_to_json(r, st, "full", false);
            bool stated = false;
            for (const char *side : { "left", "right" }) {
                for (const Json::Value &l : j["arms"][side]["limitations"]) {
                    stated |= l["kind"].asString() == "walk_domain"
                           && l["knob"].asString() == "bounds.time_budget_ms";
                }
            }
            EXPECT_TRUE(stated) << what;
            EXPECT_EQ("partial", j["outcome"]["walks"].asString()) << what;
        }
    }
    EXPECT_GT(stopped, 10u);
    // not a guarantee of the clock, so only reported: how many stops fell beyond depth 0
    std::cerr << "time stops: " << stopped << ", " << mid_walk << " beyond depth 0" << std::endl;
}


// Allocation denial at every admission (§14 freeze-gate fixtures: a denial around a
// switch, a split and a merge leaves runs, events, splits and paths consistent)
TEST(Graphlet, DenialLeavesAConsistentPrefix) {
    size_t denials = 0, mid_level = 0;
    for (const BudgetCase &c : budget_cases()) {
        std::vector<Admission> admissions;
        WalkerHooks record;
        record.deny = [&](const Admission &a) { admissions.push_back(a); return false; };
        const SeedResult base = run_case(c, c.st, &record);
        EXPECT_FALSE(base.resource_stop) << c.name;
        const size_t n = admissions.size();
        ASSERT_GT(n, 0u) << c.name;
        // every admission of the small cases, evenly spaced ones of the dense graphs
        const size_t most = c.name.rfind("mode ", 0) == 0 ? 8 : 40;
        std::vector<size_t> ordinals;
        for (size_t i = 0; i < std::min(n, most); ++i) {
            ordinals.push_back(n <= most ? i : i * (n - 1) / (most - 1));
        }
        LabelOracle oracle(*c.anno);
        for (size_t ordinal : ordinals) {
            const Admission &denied = admissions[ordinal];
            const std::string what = c.name + " denied #" + std::to_string(ordinal) + " ("
                                   + to_string(denied.arm) + " at " + std::to_string(denied.at_bp) + ")";
            WalkerHooks deny;
            deny.deny = [&](const Admission &a) { return a.ordinal == ordinal; };
            const SeedResult r = run_case(c, c.st, &deny);
            denials++;
            ASSERT_TRUE(r.resource_stop) << what;
            EXPECT_EQ(ResourceStop::MEMORY, r.resource_stop->resource) << what;
            EXPECT_EQ(denied.arm, r.resource_stop->arm) << what;
            EXPECT_EQ(denied.at_bp, r.resource_stop->at_bp) << what;
            const ArmResult &arm = r.arms[static_cast<size_t>(denied.arm)];
            EXPECT_EQ(ArmResult::TRUNCATED, arm.status) << what;
            EXPECT_EQ(denied.at_bp, arm.complete_to_bp) << what;
            ASSERT_TRUE(arm.cap_trigger) << what;
            EXPECT_EQ(EndReason::RESOURCE_LIMIT, arm.cap_trigger->reason) << what;
            EXPECT_EQ(denied.at_bp, arm.cap_trigger->at_bp) << what;
            mid_level += ordinal > 0 && admissions[ordinal - 1].arm == denied.arm
                       && admissions[ordinal - 1].at_bp == denied.at_bp;
            for (size_t a = 0; a < 2; ++a) {
                const ArmResult &ra = r.arms[a];
                if (!ra.requested)
                    continue;
                const std::string at = what + " arm " + to_string(ra.arm);
                check_bins(ra, at);
                // every walk up to the boundary is the unbudgeted walk's, and no other
                EXPECT_EQ(walks_upto(base.arms[a], ra.complete_to_bp),
                          walks_upto(ra, ra.complete_to_bp)) << at;
                // the other arm stops with the locus: at the same depth, or one level on
                // when it ran that level before the denied arm
                if (&ra != &arm && ra.status != ArmResult::COMPLETE) {
                    EXPECT_GE(ra.complete_to_bp, denied.at_bp) << at;
                    EXPECT_LE(ra.complete_to_bp, denied.at_bp + 1) << at;
                }
                // beyond the boundary only the censored heads and their siblings' ends
                for (const PathResult &p : ra.paths) {
                    if (p.length_bp >= ra.complete_to_bp && p.path_reason) {
                        EXPECT_TRUE(*p.path_reason == EndReason::RESOURCE_LIMIT
                                    || *p.path_reason == EndReason::MAX_EXTENSION
                                    || (c.st.label_mode == LabelMode::ANNOTATE
                                        && !is_resource_stop(*p.path_reason)))
                            << at << " path " << p.id << " " << to_string(*p.path_reason);
                    }
                }
            }
            const std::string text = check_serialised(r, c.seed, c.st, oracle, what);
            EXPECT_NE(std::string::npos, text.find("\nQ locus memory traversal ")) << what;
            const Json::Value j = seed_result_to_json(r, c.st, "summary", false);
            EXPECT_EQ("partial", j["outcome"]["walks"].asString()) << what;
            bool stated = false;
            for (const Json::Value &l : j["arms"][to_string(denied.arm)]["limitations"]) {
                stated |= l["kind"].asString() == "walk_domain"
                       && l["knob"].asString() == "bounds.max_memory_mb";
            }
            EXPECT_TRUE(stated) << what;
        }
    }
    EXPECT_GT(denials, 200u);
    EXPECT_GT(mid_level, 10u) << "denials between two heads of one level";
}

namespace {

// The cache allotments a memory budget B reserves of itself (Walker::init_budgets): what
// is left of B for the account an unbudgeted run records.
uint64_t usable(uint64_t budget) {
    return budget - std::min<uint64_t>(budget / 4, uint64_t(64) << 20)
                  - std::min<uint64_t>(budget / 16, uint64_t(64) << 20);
}

// the smallest budget whose usable part is at least |need|
uint64_t budget_for(uint64_t need) {
    uint64_t lo = need, hi = 4 * need + 64;
    while (lo < hi) {
        const uint64_t mid = lo + (hi - lo) / 2;
        if (usable(mid) >= need) hi = mid; else lo = mid + 1;
    }
    return lo;
}

// what an admission needs: the account before it and after it
uint64_t need_of(const Admission &a) {
    return std::max(a.total, a.total - a.reserve + a.cost);
}

} // namespace

// A memory budget between two heads of one level (§14 freeze gate, "mid-level traversal"):
// every head before the refused one is admitted, the refused one and the rest of its
// level are censored, and the stop is stated (resource_stop, Q, walk_domain on
// bounds.max_memory_mb, memory_bound_soft, 'Y' in T, R, A and B).
TEST(Graphlet, MemoryBudgetStopsMidLevel) {
    size_t checked = 0;
    for (const BudgetCase &c : budget_cases()) {
        Strategy st = c.st;
        st.delivery = delivery_costs("full", true);
        std::vector<Admission> admissions;
        WalkerHooks record;
        record.deny = [&](const Admission &a) { admissions.push_back(a); return false; };
        run_case(c, st, &record);
        // the first admission inside a level (not its arm's first there) that needs more
        // than everything before it
        uint64_t before = 0;
        std::optional<size_t> pick;
        for (size_t i = 0; i < admissions.size(); ++i) {
            const Admission &a = admissions[i];
            const bool inside = i > 0 && admissions[i - 1].arm == a.arm
                              && admissions[i - 1].at_bp == a.at_bp && a.at_bp > 0;
            if (inside && need_of(a) > before && usable(budget_for(before)) < need_of(a)) {
                pick = i;
                break;
            }
            before = std::max(before, need_of(a));
        }
        if (!pick)
            continue;
        const Admission &refused = admissions[*pick];
        st.max_memory_bytes = budget_for(before);
        ASSERT_LT(usable(st.max_memory_bytes), need_of(refused)) << c.name;
        const std::string what = c.name + " budget " + std::to_string(st.max_memory_bytes)
                               + " refusing #" + std::to_string(*pick);
        const SeedResult r = run_case(c, st);
        checked++;
        ASSERT_TRUE(r.resource_stop) << what;
        EXPECT_EQ(ResourceStop::MEMORY, r.resource_stop->resource) << what;
        EXPECT_EQ(refused.arm, r.resource_stop->arm) << what;
        EXPECT_EQ(refused.at_bp, r.resource_stop->at_bp) << what;
        EXPECT_LE(r.resource_stop->used, r.resource_stop->limit) << what;
        EXPECT_LE(r.account.memory_peak, st.max_memory_bytes) << what;
        EXPECT_EQ(st.max_memory_bytes, r.account.memory_limit) << what;
        const ArmResult &arm = r.arms[static_cast<size_t>(refused.arm)];
        EXPECT_EQ(refused.at_bp, arm.complete_to_bp) << what;
        ASSERT_TRUE(arm.cap_trigger) << what;
        EXPECT_EQ(EndReason::RESOURCE_LIMIT, arm.cap_trigger->reason) << what;
        LabelOracle oracle(*c.anno);
        const std::string text = check_serialised(r, c.seed, st, oracle, what);
        EXPECT_NE(std::string::npos, text.find("\nQ locus memory traversal ")) << what;
        EXPECT_NE(std::string::npos, text.find("\nT Y")) << what;
        EXPECT_NE(std::string::npos, text.find(" Y,")) << what << ": the A record's cap trigger";
        if (c.st.label_mode == LabelMode::CONSTRAIN) {
            // label ends (annotate mode has none: its leaves carry the path reason only)
            EXPECT_NE(std::string::npos, text.find("Y:")) << what << ": a B record's label ends";
            EXPECT_NE(std::string::npos, text.find(" Y 0 ")) << what << ": an R record ending Y";
        }
        const Json::Value j = seed_result_to_json(r, st, "full", false);
        EXPECT_EQ("partial", j["outcome"]["walks"].asString()) << what;
        EXPECT_EQ("memory", j["resource_stop"]["resource"].asString()) << what;
        EXPECT_EQ("locus", j["resource_stop"]["scope"].asString()) << what;
        EXPECT_EQ("traversal", j["resource_stop"]["phase"].asString()) << what;
        bool soft = false;
        for (const Json::Value &l : j["limitations"]) soft |= l["kind"].asString() == "memory_bound_soft";
        EXPECT_TRUE(soft) << what;
    }
    EXPECT_GT(checked, 10u);
}

// A budget that admits exactly what levels 0..d need completes them on both arms, one
// byte less does not (the admission is per head, deterministic, and includes what ending
// and delivering every admitted head costs)
TEST(Graphlet, MemoryBudgetAdmitsALevel) {
    size_t checked = 0;
    for (const BudgetCase &c : budget_cases()) {
        Strategy st = c.st;
        st.delivery = delivery_costs("graphlet", true);
        std::vector<Admission> admissions;
        WalkerHooks record;
        record.deny = [&](const Admission &a) { admissions.push_back(a); return false; };
        const SeedResult base = run_case(c, st, &record);
        uint64_t deepest = 0;
        for (const Admission &a : admissions) deepest = std::max(deepest, a.at_bp);
        for (uint64_t d = 1; d + 1 <= deepest; d += 2) {
            uint64_t need = 0;
            for (const Admission &a : admissions) {
                if (a.at_bp <= d) need = std::max(need, need_of(a));
            }
            const std::string what = c.name + " through level " + std::to_string(d);
            st.max_memory_bytes = budget_for(need);
            const SeedResult ok = run_case(c, st);
            for (size_t a = 0; a < 2; ++a) {
                if (ok.arms[a].requested) {
                    EXPECT_GE(ok.arms[a].complete_to_bp,
                              std::min(d + 1, base.arms[a].complete_to_bp)) << what;
                }
            }
            st.max_memory_bytes = budget_for(need) - 1;
            if (usable(st.max_memory_bytes) >= need)
                continue;
            const SeedResult cut = run_case(c, st);
            ASSERT_TRUE(cut.resource_stop) << what << " minus one";
            EXPECT_LE(cut.resource_stop->at_bp, d) << what << " minus one";
            const ArmResult &arm = cut.arms[static_cast<size_t>(cut.resource_stop->arm)];
            EXPECT_LE(arm.complete_to_bp, d) << what << " minus one";
            checked++;
        }
    }
    EXPECT_GT(checked, 20u);
}

// The work budget: charged work units, checked at least every kWorkCheckInterval units —
// also inside a derivation, so one head cannot run past it unchecked. 300 seed labels
// switching into 300 others on the step after the seed price 90,000 pairs in ONE
// derivation; a budget of 1,000 units beyond what the seed's validation spends is refused
// inside it, and the seed's root is censored (complete_to_bp 0), the overshoot bounded by
// the interval plus one target's pairs.
TEST(Graphlet, WorkBudgetIsCheckedInsideADerivation) {
    const std::string X = random_seq(30, 901), P = random_seq(30, 902);
    std::vector<std::string> seqs, labels, seed_labels, extra;
    for (size_t i = 0; i < 300; ++i) {
        seqs.push_back(X);
        labels.push_back("A" + std::to_string(i));
        seed_labels.push_back(labels.back());
        seqs.push_back(X.substr(X.size() - 10) + P);
        labels.push_back("B" + std::to_string(i));
        extra.push_back(labels.back());
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, seqs, labels, DeBruijnGraph::BASIC);
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.extra = extra;
    st.loss_budget = 1;
    st.max_switch_sources = Strategy::kUnlimited;
    st.max_label_branches = Strategy::kUnlimited;
    st.max_extension_bp = 20;
    const LabelChangeCost cost = LabelChangeCost::table({}, 1.0);
    LabelOracle oracle(*anno);
    const Seed seed = seed_of(X, seed_labels);
    const SeedResult free = traverse_seed(oracle, seed, st, cost);
    ASSERT_FALSE(free.resource_stop);
    EXPECT_GE(free.arms[1].pair_evaluations, 90'000u);
    // the seed's validation is charged too: 20 k-mers, each 8 units and its 300 hits
    EXPECT_EQ(20u * (8 + 300), free.account.work_seed);
    EXPECT_EQ(free.account.work_seed + free.arms[1].work_units, free.account.work_used);

    // 1,000 units beyond what the seed phase spends
    const uint64_t budget = free.account.work_seed + 1000;
    st.max_work_units = budget;
    const SeedResult r = traverse_seed(oracle, seed, st, cost);
    ASSERT_TRUE(r.resource_stop);
    EXPECT_EQ(ResourceStop::WORK, r.resource_stop->resource);
    EXPECT_EQ(0u, r.resource_stop->at_bp);
    EXPECT_EQ(0u, r.arms[1].complete_to_bp);
    EXPECT_EQ(0u, r.arms[1].steps);
    EXPECT_GT(r.resource_stop->used, static_cast<double>(budget));
    // checked when the interval passed, within one target's pricing of it
    EXPECT_LE(r.resource_stop->used, static_cast<double>(budget) + kWorkCheckInterval + 2 * 300 + 64);
    EXPECT_LT(r.arms[1].pair_evaluations, free.arms[1].pair_evaluations);
    EXPECT_EQ(r.account.work_seed + r.arms[1].work_units, r.account.work_used);
    const std::string text = check_serialised(r, seed, st, oracle, "work inside a derivation");
    const std::string b = std::to_string(budget);
    EXPECT_NE(std::string::npos, text.find("\nQ locus work traversal i:" + b + " i:" + b + " i:"));
    const Json::Value j = seed_result_to_json(r, st, "full", false);
    EXPECT_EQ("bounds.max_work_units", j["arms"]["right"]["limitations"][0]["knob"].asString());
    EXPECT_EQ(r.arms[1].work_units, j["arms"]["right"]["counters"]["work_units"].asUInt64());
    EXPECT_EQ("partial", j["outcome"]["walks"].asString());
}

// Budgets keep the determinism contract (§6.8): the same result for batch_kmers 1 and 64,
// for seeds in another order and for one seed per request — memory and work alike (the
// caches get fixed allotments and work is charged on consumption, never on fetches).
TEST(Graphlet, BudgetsAreDeterministic) {
    size_t stopped = 0;
    std::vector<DeBruijnGraph::Mode> modes { DeBruijnGraph::BASIC };
#if ! _PROTEIN_GRAPH
    modes.push_back(DeBruijnGraph::PRIMARY);
#endif
    for (auto mode : modes) {
        std::vector<std::string> seqs, labels;
        for (uint32_t i = 0; i < 6; ++i) {
            seqs.push_back("AAA" + random_seq(40, 70 + i));
            labels.push_back(std::string(1, "CDECDF"[i]));
        }
        auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(3, seqs, labels, mode);
        for (const char *budget : { R"("max_work_units": 150)", R"("max_work_units": 900)",
                                    R"("max_memory_mb": 1)" }) {
            for (const char *detail : { "full", "graphlet" }) {
                auto request = [&](const std::vector<std::string> &seeds, size_t batch) {
                    Json::Value r;
                    for (const std::string &s : seeds) {
                        Json::Value seed;
                        seed["sequence"] = s;
                        r["seeds"].append(seed);
                    }
                    r["strategy"] = parse_json(std::string(R"({"labels": {"mode": "annotate"},
                        "bounds": {"max_extension_bp": 12, )") + budget + R"(},
                        "output": {"timing": false, "detail": ")" + detail + R"("},
                        "annotation": {"batch_kmers": )" + std::to_string(batch) + "}}");
                    return process_traverse_request(r, *anno, "");
                };
                const std::vector<std::string> seeds { "AAA", seqs[1].substr(5, 8), seqs[2].substr(9, 7) };
                const Json::Value all = request(seeds, 64);
                const Json::Value one = request(seeds, 1);
                const Json::Value reversed = request({ seeds[2], seeds[1], seeds[0] }, 64);
                // the H record states the seed's position in its request: not part of the
                // result (§6.8 compares results across orders and splits)
                auto result = [](const Json::Value &out, size_t i) {
                    Json::Value r = out["results"][Json::ArrayIndex(i)];
                    if (r.isMember("graphlet")) {
                        std::string text = r["graphlet"].asString();
                        const size_t walk = text.find(" walk ");
                        const size_t index = text.rfind(' ', walk - 1);
                        text.replace(index + 1, walk - index - 1, "#");
                        r["graphlet"] = text;
                    }
                    return compact(r);
                };
                for (size_t i = 0; i < seeds.size(); ++i) {
                    const std::string what = std::string(budget) + " " + detail + " seed " + std::to_string(i);
                    stopped += all["results"][Json::ArrayIndex(i)].isMember("resource_stop");
                    EXPECT_EQ(result(all, i), result(one, i)) << what << " batch_kmers";
                    EXPECT_EQ(result(all, i), result(reversed, 2 - i)) << what << " order";
                    EXPECT_EQ(result(all, i), result(request({ seeds[i] }, 64), 0)) << what << " alone";
                }
            }
        }
    }
    EXPECT_GT(stopped, 4u);
}

// A comb-shaped trie: a 3,000 bp backbone with a one-base dead-end tip at every node
// (§14 freeze gate). Paths are lazy, so finalisation and the MGT body are linear in the
// trie and the walk completes under a budget; detail full spells every leaf's chain
// (quadratic: ~4.5 million entries), which the same budget refuses early — a valid
// prefix, its chains exact, not an explosion at the end.
TEST(Graphlet, CombTrieStaysLinear) {
    const size_t k = 15, n = 3000;
    std::string B;
    std::vector<std::string> tips;
    for (uint32_t t = 41; ; ++t) {
        B = random_seq(n + k, t);
        // no (k-1)-mer twice, so the backbone is one unbranched path apart from the tips
        std::set<std::string> seen;
        bool clean = true;
        for (size_t i = 0; clean && i + k - 1 <= B.size(); ++i) {
            clean = seen.insert(B.substr(i, k - 1)).second;
        }
        if (!clean)
            continue;
        tips.clear();
        for (size_t i = k; i + 1 < B.size() && clean; ++i) {
            // the tip leaves the backbone after B[0, i) with another base; its last
            // (k-1)-mer must be new, or the tip would continue somewhere
            std::string tip;
            for (char ch : std::string("ACGT")) {
                if (ch == B[i])
                    continue;
                std::string cand = B.substr(i - k + 1, k - 1) + ch;
                if (seen.insert(cand.substr(1)).second) {
                    tip = cand;
                    break;
                }
            }
            clean = !tip.empty();
            tips.push_back(tip);
        }
        if (clean)
            break;
    }
    std::vector<std::string> seqs { B };
    seqs.insert(seqs.end(), tips.begin(), tips.end());
    std::vector<std::string> labels(seqs.size(), "A");
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(k, seqs, labels,
                                                                               DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    const Seed seed = seed_of(B.substr(0, k), { "A" });
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.exhaustive = true;
    st.max_label_branches = Strategy::kUnlimited;
    st.max_splits_per_path = Strategy::kUnlimited;
    st.merge_reconverge = false;
    st.max_extension_bp = n;
    st.max_branch_events = Strategy::kUnlimited;
    st.continuation_bp = 0;
    st.max_memory_bytes = uint64_t(64) << 20;

    st.delivery = delivery_costs("graphlet", true);
    const SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    ASSERT_FALSE(r.resource_stop) << "the graphlet of the comb fits 64 MiB";
    const ArmResult &arm = r.arms[1];
    EXPECT_EQ(ArmResult::COMPLETE, arm.status);
    EXPECT_GT(arm.segments.size(), 2 * n - 100);
    EXPECT_GT(arm.paths.size(), n - 50);
    // the account is linear: a bounded charge per segment (beyond the 20 MiB of cache
    // allotments of a 64 MiB budget), whatever the chains sum to
    EXPECT_LE(r.account.memory_peak, st.max_memory_bytes);
    EXPECT_LE(r.account.memory_final, r.account.memory_peak);
    EXPECT_LT((r.account.memory_final - (uint64_t(20) << 20)) / arm.segments.size(), 8192u);
    uint64_t chains = 0;
    for (const PathResult &p : arm.paths) chains += path_segments(arm, p).size();
    EXPECT_GT(chains, uint64_t(n) * n / 3) << "a comb: the chains are quadratic";
    const std::string text = graphlet_text(r, seed, st, context_of(oracle),
                                           seed_result_to_json(r, st, "graphlet", false));
    EXPECT_LT(text.size(), 120 * arm.segments.size()) << "the body is linear";

    // detail full: the chains are output, charged per expansion, and refused early
    st.delivery = delivery_costs("full", true);
    const SeedResult full = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    ASSERT_TRUE(full.resource_stop);
    EXPECT_EQ(ResourceStop::MEMORY, full.resource_stop->resource);
    EXPECT_LT(full.arms[1].complete_to_bp, n / 2);
    EXPECT_GT(full.arms[1].complete_to_bp, 100u);
    EXPECT_LE(full.account.memory_peak, st.max_memory_bytes);
    // the delivered prefix's chains are exact (the reader rebuilds paths[].segments)
    check_serialised(full, seed, st, oracle, "comb full");
}

namespace {

// malloc's footprint of |n| requested bytes, at most: a 16-byte quantum and an 8-byte
// header, 32 bytes at least (glibc; macOS's tiny zone rounds to 16 without a header)
uint64_t heap_bytes(uint64_t n) {
    return std::max<uint64_t>(32, (n + 8 + 15) / 16 * 16);
}

// What jsoncpp allocates for the tree |v|, from this platform's types: per object or array
// its map, per member a node (rb-tree links and colour, the pair of CZString and Value) and
// its key's copy, per element a node, per string its length-prefixed buffer
uint64_t json_tree_bytes(const Json::Value &v) {
    static const uint64_t node
        = heap_bytes(sizeof(Json::Value::ObjectValues::value_type) + 4 * sizeof(void*));
    static const uint64_t map = heap_bytes(sizeof(Json::Value::ObjectValues));
    uint64_t bytes = 0;
    if (v.isObject()) {
        bytes += map;
        for (const std::string &name : v.getMemberNames()) {
            bytes += node + heap_bytes(name.size() + 1) + json_tree_bytes(v[name]);
        }
    } else if (v.isArray()) {
        bytes += map;
        for (const Json::Value &x : v) {
            bytes += node + json_tree_bytes(x);
        }
    } else if (v.isString()) {
        bytes += heap_bytes(sizeof(unsigned) + v.asString().size() + 1);
    }
    return bytes;
}

// The peak of delivering one seed's JSON |j| as the response holds it (a graphlet's with
// its body): the tree with Json::writeString's text three times (its stream's buffer,
// grown by doubling to less than twice the text, and the copy it returns); later the text,
// the deflate output (at most 0.03% and a few bytes larger) and zlib's 256 KiB state; or a
// phase before either (|before|: the graphlet writer's text and its JSON copy). The text is
// the larger of the server's compact form and the CLI's indented one, nested in a response.
uint64_t delivery_peak(const Json::Value &j, uint64_t before = 0) {
    Json::StreamWriterBuilder compact, indented;
    compact["indentation"] = "";
    indented["indentation"] = "  ";
    const std::string c = Json::writeString(compact, j), i = Json::writeString(indented, j);
    // a result is two levels deep in a response (the root object, the results array)
    const uint64_t lines = std::count(i.begin(), i.end(), '\n') + 1;
    const uint64_t text = std::max<uint64_t>(c.size(), i.size() + 4 * lines);
    return std::max({ json_tree_bytes(j) + 3 * text,
                      2 * text + text / 4096 + 64 + (uint64_t(1) << 18), before });
}

// The assumptions delivery_costs() derives its bounds from, checked on what the
// serialisers wrote: keys of at most 25 characters, at most ten levels in a response (the
// result itself is at level 2), numbers of at most 24 characters and fixed strings of at most
// 28; names, sequences, effects (640) and messages (1 KiB) are priced by their length
void check_widths(const Json::Value &v, size_t level, const std::string &key,
                  const std::string &what) {
    static const std::set<std::string> kPriced {
        "name", "label", "seed_id", "sequence", "graphlet", "labels", "error"
    };
    EXPECT_LE(level, 10u) << what << " " << key;
    if (v.isObject()) {
        for (const std::string &name : v.getMemberNames()) {
            EXPECT_LE(name.size(), 25u) << what << " key " << name;
            check_widths(v[name], level + 1, name, what);
        }
    } else if (v.isArray()) {
        for (const Json::Value &x : v) {
            check_widths(x, level + 1, key, what);
        }
    } else if (v.isString()) {
        const size_t n = v.asString().size();
        if (key == "effect") {
            EXPECT_LE(n, 640u) << what;
        } else if (key == "knob") {
            // a knob's name, the longest output.max_coordinate_occurrences (33): still a short
            // string's buffer (at most 35 characters in 48 B) and a member's text of 80 B
            EXPECT_LE(n, 35u) << what << " " << v.asString();
        } else if (key == "message") {
            EXPECT_LE(n, 1024u) << what;
        } else if (!kPriced.count(key)) {
            EXPECT_LE(n, 28u) << what << " " << key << ": " << v.asString();
        }
    } else if (!v.isNull()) {
        Json::StreamWriterBuilder compact;
        compact["indentation"] = "";
        EXPECT_LE(Json::writeString(compact, v).size(), 24u) << what << " " << key;
    }
}

// The cases of names whose delivered size is far from their length (GPT review of stage 2,
// finding 1): control characters (six bytes each in JSON), '%' (three in MGT), quotes and
// backslashes, tabs and DEL, and two-, three- and four-byte UTF-8 (six or twelve bytes per
// character in JSON) — as seed labels (written twice by detail full), extra and recorded
// labels, a dropped label and the seed_id
std::vector<BudgetCase> adversarial_name_cases() {
    // long enough that the names, not the fixed part, dominate what is delivered
    const size_t n = 24'000;
    auto repeat = [](const std::string &unit, size_t bytes) {
        std::string s;
        while (s.size() + unit.size() <= bytes) s += unit;
        return s;
    };
    const std::vector<std::string> names {
        std::string(n, '\x01'),
        std::string(n, '%'),
        repeat("\"\\", n),
        repeat("\t\x7f\x1f", n),
        repeat("\xc3\xa9", n),
        repeat("\xe2\x82\xac", n),
        repeat("\xf0\x9f\x98\x80%\x02", n),
    };
    std::vector<std::string> seqs, labels;
    for (uint32_t i = 0; i < names.size(); ++i) {
        seqs.push_back("AAA" + random_seq(14, 7700 + i));
        labels.push_back(names[i]);
    }
    std::shared_ptr<graph::AnnotatedDBG> anno
        = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(3, seqs, labels,
                                                                         DeBruijnGraph::BASIC);
    std::vector<BudgetCase> cases;
    // every label carries AAA; the last one's record is cut so that it does not carry the
    // seed AAAx and is dropped
    for (size_t first = 0; first + 1 < names.size(); ++first) {
        Seed seed = seed_of(seqs[first].substr(0, 6), { names[first], names[names.size() - 1] });
        seed.seed_id = std::string(64, '\x03') + names[first].substr(0, 40);
        LabelChangeCost constant;
        Strategy st = strategy_of(R"({"bounds": {"max_extension_bp": 8},
            "labels": {"loss_budget": 1, "change_cost": {"model": "constant", "value": 1}}})",
            &constant);
        st.extra = { names[(first + 1) % (names.size() - 1)] };
        cases.push_back({ "names " + std::to_string(first), anno, seed, st, constant });
    }
    cases.push_back({ "names annotate", anno, seed_of("AAA"),
                      strategy_of(R"({"labels": {"mode": "annotate"}, "bounds": {"max_extension_bp": 5}})") });
    return cases;
}

// Floats canonical MGT writes wide (review of the stage-2 recheck, P1: 1e-300 is written
// positionally, 302 characters, where the model assumed 24): on a graph whose k-mers
// alternate between labels A and B, so that every step switches, switches of 1e-300 (losses
// summed to some 300 digits), switches of 0.1 under a loss budget of 1e300 (sums of 17
// digits), and a time budget of 1e-300 (its time stop states it in K, Q and A)
std::vector<BudgetCase> wide_float_cases() {
    const std::string source = random_seq(1230, 9900);
    std::vector<std::string> kmers, labels;
    for (size_t i = 0; i + 31 <= source.size(); ++i) {
        kmers.push_back(source.substr(i, 31));
        labels.push_back(i % 2 ? "B" : "A");
    }
    std::shared_ptr<graph::AnnotatedDBG> anno
        = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(31, kmers, labels,
                                                                         DeBruijnGraph::BASIC);
    const Seed seed = seed_of(source.substr(0, 31), { "A" });
    const Strategy st = strategy_of(R"({"direction": "right", "labels": {"extra": ["B"],
        "loss_budget": 1, "change_cost": {"model": "constant", "value": 1}},
        "bounds": {"max_extension_bp": 1200}})");
    std::vector<BudgetCase> cases;
    cases.push_back({ "wide tiny switches", anno, seed, st, LabelChangeCost::table({}, 1e-300) });
    Strategy sums = st;
    sums.loss_budget = 1e300;
    cases.push_back({ "wide sums", anno, seed, sums, LabelChangeCost::constant(0.1) });
    Strategy late = st;
    late.time_budget_ms = 1e-300;
    late.max_work_units = 1'000'000'000;      // a time stop is stated (Q) under a budget
    cases.push_back({ "wide time budget", anno, seed, late, LabelChangeCost::constant(1) });
    return cases;
}

} // namespace

TEST(GraphletCodec, JsonEscapedSizeIsTheWriters) {
    Json::StreamWriterBuilder compact;
    compact["indentation"] = "";
    auto written = [&](const std::string &s) {
        return Json::writeString(compact, Json::Value(s)).size() - 2;
    };
    std::vector<std::string> cases;
    for (int c = 1; c < 128; ++c) {
        cases.push_back(std::string(3, static_cast<char>(c)));
    }
    cases.push_back("\xc3\xa9\xe2\x82\xac\xf0\x9f\x98\x80plain");
    std::mt19937 rng(17);
    for (size_t i = 0; i < 200; ++i) {
        std::string s;
        for (size_t j = 0; j < 40; ++j) {
            switch (rng() % 5) {
                case 0: s.push_back(static_cast<char>(rng() % 0x80)); break;
                case 1: s += "\xc3\xa9"; break;
                case 2: s += "\xe2\x82\xac"; break;
                case 3: s += "\xf0\x9f\x98\x80"; break;
                default: s.push_back(static_cast<char>(1 + rng() % 0x1f)); break;
            }
        }
        cases.push_back(s);
    }
    for (const std::string &s : cases) {
        if (s.find('\0') != std::string::npos)
            continue;
        EXPECT_EQ(written(s), json_escaped_size(s)) << s;
    }
    // bytes that are not UTF-8 (refused per seed before any output; priced at least as
    // the writer would write them)
    for (const std::string &s : { std::string("\xff\xfe"), std::string("\x80"), std::string("a\xc3"),
                                  std::string("\xe2\x82"), std::string("\xf8\x88\x80\x80\x80") }) {
        EXPECT_GE(json_escaped_size(s), written(s)) << s;
    }
}

namespace {

// what the delivery model charges for |r| in |detail| (DeliveryCosts per object), from
// the result's own objects: a lower bound of the walker's charge, which prices lists at
// their bounds
uint64_t modelled_delivery(const SeedResult &r, const DeliveryCosts &d) {
    uint64_t m = d.fixed;
    auto name = [&](const std::string &n, DeliveryCosts::Name use) {
        return d.name ? d.name(n, use) : 0;
    };
    m += name(r.seed_id, DeliveryCosts::Name::SEED_ID);
    for (size_t i = 0; i < r.label_dict.size(); ++i) {
        m += d.label + name(r.label_dict[i].name, i < r.num_seed_labels
                                                   ? DeliveryCosts::Name::SEED_LABEL
                                                   : DeliveryCosts::Name::LABEL);
    }
    for (const DroppedLabel &dl : r.dropped_labels) {
        m += d.dropped + dl.runs.size() * d.dropped_run + name(dl.name, DeliveryCosts::Name::DROPPED_LABEL);
    }
    for (const ArmResult &arm : r.arms) {
        for (const Segment &s : arm.segments) {
            m += d.segment + s.length_bp * d.base;
            size_t labels = s.labels_start.size() + s.labels_end.size();
            if (s.parents.size() > 1)
                m += s.parents.size() * d.merge_parent;
            for (const auto &via : s.labels_via_parent) labels += via.size();
            for (const LabelSetRun &run : s.label_sets) {
                labels += run.labels.size();
                m += d.presence_run;
            }
            m += labels * d.segment_label;
            for (const Event &ev : s.events) m += d.event + ev.labels.size() * d.event_label;
        }
        m += arm.runs.size() * d.run;
        for (const PathResult &p : arm.paths) {
            m += d.leaf + p.end_labels.size() * d.leaf_label
               + path_segments(arm, p).size() * d.chain_entry;
            if (p.continuation) {
                m += p.continuation->labels.size() * d.leaf_label
                   + p.continuation->sequence.size() * d.continuation_base;
            }
        }
        for (const Split &sp : arm.splits) {
            m += d.split;
            for (const SplitBranch &b : sp.branches) m += d.split_branch + b.labels.size() * d.segment_label;
        }
        for (const BranchEvent &be : arm.branch_events) {
            m += d.branch_event + (be.chars.size() + be.ambiguous.size() + be.dropped.size()) * d.branch_event_entry;
            for (const auto &rf : be.refused) m += d.refusal + rf.labels.size() * d.branch_event_entry;
        }
        m += arm.growth.size() * d.bin;
    }
    // record coordinates (opt-in): their entries and the occurrences listed
    if (r.coordinates_recorded) {
        for (const SeedCoordinates &sc : r.seed_coordinates) {
            m += d.coordinate_seed + sc.starts.size() * d.occurrence;
        }
        for (const ArmResult &arm : r.arms) {
            for (const RunCoordinates &rc : arm.run_coordinates) {
                m += d.coordinate_run + rc.ends.size() * d.occurrence;
            }
        }
    }
    // a set derived from part of the seed states its derivation (D3)
    if (r.derivation_partial)
        m += d.extra_limitation;
    return m;
}

} // namespace

// The delivery model bounds what the serialisers hold (§14: every committed prefix stays
// deliverable within what was reserved for it), as a demonstrated upper bound and not an
// average (the owner's answer to the stage-2 review): per detail, the peak of the JSON tree
// with writeString's three texts, the compressor's, or the graphlet writer's (its text and
// the JSON copy), computed from the real serialisers' output and this platform's jsoncpp
// types — with names whose delivered size is far from their length (GPT review of stage 2,
// finding 1). The walker charges at least the model (memory_final), the serialised output
// keeps to the widths the model assumes, and the model stays within a small factor of the
// output, or the budget would refuse far too early.
TEST(Graphlet, DeliveryCostsBoundTheOutput) {
    size_t checked = 0, names = 0, wide = 0;
    std::vector<BudgetCase> cases = budget_cases();
    for (BudgetCase &c : adversarial_name_cases()) {
        cases.push_back(std::move(c));
    }
    for (BudgetCase &c : wide_float_cases()) {
        cases.push_back(std::move(c));
    }
    for (const BudgetCase &c : cases) {
        LabelOracle oracle(*c.anno);
        // the widest the result's floats can be written, as process_traverse_request prices it
        const uint64_t width = mgt_float_width(c.st, c.cost);
        for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
            const std::string what = c.name + " " + detail;
            Strategy st = c.st;
            st.delivery = delivery_costs(detail, st.sequences, width);
            const SeedResult r = run_case(c, st);
            const DeliveryCosts &d = st.delivery;
            const uint64_t model = modelled_delivery(r, d);
            Json::Value j = seed_result_to_json(r, st, detail, false);
            uint64_t before = 0;
            if (std::string(detail) == "graphlet") {
                const std::string text = graphlet_text(r, c.seed, st, context_of(oracle), j);
                // the writer's text (less than twice the body by doubling) and its JSON
                // copy, with the writer's index (per segment, run and label end) and the
                // JSON summary
                uint64_t index = 0;
                for (const ArmResult &arm : r.arms) {
                    index += 40 * arm.segments.size() + 8 * arm.runs.size();
                    for (const Segment &s : arm.segments) {
                        for (const Event &ev : s.events) index += 96 * (ev.type == EventType::LABEL_END);
                    }
                }
                before = 3 * text.size() + index + json_tree_bytes(j);
                // the records whose floats the model widens keep to the widths it assumes
                const uint64_t extra = width - kMgtFloatWidth;
                for (const std::string &line : split_on(text.substr(0, text.size() - 1), '\n')) {
                    if (line[0] == 'R') {
                        EXPECT_LE(line.size(), 172 + 3 * extra) << what << ": " << line.substr(0, 80);
                    }
                    if (line.rfind("E ", 0) == 0 && line.find(" s ") != std::string::npos) {
                        EXPECT_LE(line.size(), 64 + extra) << what << ": " << line.substr(0, 80);
                    }
                    if (line[0] == 'A' || line[0] == 'Q' || line[0] == 'K') {
                        for (const std::string &field : split_on(line, ' ')) {
                            for (const std::string &part : split_on(field, ',')) {
                                const size_t colon = part.find(':');
                                const std::string number = colon == std::string::npos
                                    ? part : part.substr(colon + 1);
                                if (number.size() > 24 && number.find_first_not_of("0123456789.") == std::string::npos) {
                                    EXPECT_LE(number.size(), width) << what << ": " << line.substr(0, 80);
                                }
                            }
                        }
                    }
                }
                wide += width > kMgtFloatWidth;
                j["graphlet"] = text;
            }
            check_widths(j, 2, "", what);
            size_t limitations = j["limitations"].size();
            for (const std::string &side : j["arms"].getMemberNames()) {
                EXPECT_LE(j["arms"][side]["limitations"].size(), 8u) << what;
                limitations += j["arms"][side]["limitations"].size();
            }
            EXPECT_LE(limitations, 5u + 2 * 8) << what;
            const uint64_t actual = delivery_peak(j, before);
            EXPECT_LE(actual, model) << what;
            EXPECT_LE(model, r.account.memory_final) << what << ": the walker charges less than the model";
            // a bound fixed per object cannot know before the walk which fields will be
            // short, '*' or delta-coded (MGT), so the model is looser than any one output;
            // floats priced wider than 24 characters are priced at the widest the request
            // allows, which most of its floats need not reach (not a useful model there,
            // only a bound)
            const uint64_t factor = std::string(detail) == "graphlet" ? 10 : 6;
            if (width == kMgtFloatWidth) {
                EXPECT_LE(model, factor * actual) << what << ": the model is not useful";
            }
            checked++;
            names += c.name.rfind("names", 0) == 0;
        }
    }
    EXPECT_GT(checked, 40u);
    EXPECT_GE(names, 28u);
    EXPECT_GE(wide, 3u);
}

// The deadline is checked between the heads of a level too (§14: at least every W work
// units), not only between levels: a level that outlives the budget stops at its next
// head, censored with time_budget like any cap, a valid prefix. Stated as resource_stop
// only under a request budget (without one, a time stop keeps its earlier form).
TEST(Graphlet, DeadlineIsCheckedWithinALevel) {
    size_t checked = 0;
    for (const BudgetCase &c : budget_cases()) {
        std::vector<Admission> admissions;
        WalkerHooks record;
        record.deny = [&](const Admission &a) { admissions.push_back(a); return false; };
        const SeedResult base = run_case(c, c.st, &record);
        // a head followed by another of its level, at depth > 0 (depth 0 is never cut
        // by the deadline: a zero budget means one level, then the boundary check)
        std::optional<size_t> pick;
        for (size_t i = 0; i + 1 < admissions.size(); ++i) {
            if (admissions[i].at_bp > 0 && admissions[i + 1].arm == admissions[i].arm
                    && admissions[i + 1].at_bp == admissions[i].at_bp) {
                pick = i;
                break;
            }
        }
        if (!pick || checked >= 6)
            continue;
        const Admission &slow = admissions[*pick];
        Strategy st = c.st;
        st.time_budget_ms = 300;
        WalkerHooks sleep;
        sleep.deny = [&](const Admission &a) {
            if (a.ordinal == *pick)
                std::this_thread::sleep_for(std::chrono::milliseconds(400));
            return false;
        };
        const SeedResult r = run_case(c, st, &sleep);
        const std::string what = c.name + " slow #" + std::to_string(*pick);
        ASSERT_TRUE(r.resource_stop) << what;
        EXPECT_EQ(ResourceStop::TIME, r.resource_stop->resource) << what;
        EXPECT_EQ(slow.arm, r.resource_stop->arm) << what;
        EXPECT_EQ(slow.at_bp, r.resource_stop->at_bp) << what;
        const ArmResult &arm = r.arms[static_cast<size_t>(slow.arm)];
        ASSERT_TRUE(arm.cap_trigger) << what;
        EXPECT_EQ(EndReason::TIME_BUDGET, arm.cap_trigger->reason) << what;
        EXPECT_EQ(slow.at_bp, arm.complete_to_bp) << what;
        for (size_t a = 0; a < 2; ++a) {
            if (r.arms[a].requested) {
                EXPECT_EQ(walks_upto(base.arms[a], r.arms[a].complete_to_bp),
                          walks_upto(r.arms[a], r.arms[a].complete_to_bp)) << what;
            }
        }
        LabelOracle oracle(*c.anno);
        const std::string text = check_serialised(r, c.seed, st, oracle, what);
        // no request budget: the time stop keeps its form (walk_domain, no Q)
        EXPECT_EQ(std::string::npos, text.find("\nQ ")) << what;
        EXPECT_FALSE(seed_result_to_json(r, st, "full", false).isMember("resource_stop")) << what;
        st.max_work_units = 1'000'000'000;
        const Json::Value budgeted = seed_result_to_json(r, st, "full", false);
        EXPECT_EQ("time", budgeted["resource_stop"]["resource"].asString()) << what;
        EXPECT_EQ(300.0, budgeted["resource_stop"]["requested"].asDouble()) << what;
        checked++;
    }
    EXPECT_GE(checked, 4u);
}


/*
 * The review of stage 2 (DESIGN-traverse-graphlet.md §14): one regression test per
 * confirmed finding. Results without a budget stay as they were; each test states the
 * budgeted (or refused) case the finding was about.
 */
namespace {

std::string review_rc(std::string s) { ::reverse_complement(s); return s; }

// test_walker.cpp's clean blocks (k = 11): no repeated k-mer or (k-1)-mer in either
// orientation and no RC-palindromic (k-1)- or (k+1)-mer, so that a fixture's graph is
// exactly its records'
std::vector<std::string> review_blocks(const std::vector<size_t> &lengths, uint32_t seed) {
    const size_t k = 11;
    auto clean = [&](const std::string &s) {
        for (size_t len : { k - 1, k }) {
            std::set<std::string> seen;
            for (size_t i = 0; i + len <= s.size(); ++i) {
                const std::string km = s.substr(i, len);
                if (seen.count(km) || seen.count(review_rc(km)))
                    return false;
                seen.insert(km);
            }
        }
        for (size_t len : { k - 1, k + 1 }) {
            for (size_t i = 0; i + len <= s.size(); ++i) {
                if (review_rc(s.substr(i, len)) == s.substr(i, len))
                    return false;
            }
        }
        return true;
    };
    size_t total = 0;
    for (size_t l : lengths) total += l;
    std::string master;
    for (uint32_t s = seed; ; ++s) {
        master = random_seq(total, s);
        if (clean(master))
            break;
    }
    std::vector<std::string> blocks;
    size_t pos = 0;
    for (size_t l : lengths) {
        blocks.push_back(master.substr(pos, l));
        pos += l;
    }
    return blocks;
}

bool states(const Json::Value &limitations, const std::string &kind) {
    for (const Json::Value &l : limitations) {
        if (l["kind"].asString() == kind)
            return true;
    }
    return false;
}

const Json::Value* limitation_of(const Json::Value &limitations, const std::string &kind) {
    for (const Json::Value &l : limitations) {
        if (l["kind"].asString() == kind)
            return &l;
    }
    return nullptr;
}

// every label id an annotate result records (its runs, labels_start/end, the labels of
// the successors not taken, split branches and continuations)
std::set<LabelId> recorded_labels(const SeedResult &r) {
    std::set<LabelId> out;
    auto add = [&](const std::vector<LabelId> &labels) { out.insert(labels.begin(), labels.end()); };
    for (const ArmResult &arm : r.arms) {
        for (const Segment &seg : arm.segments) {
            add(seg.labels_start);
            add(seg.labels_end);
            for (const LabelSetRun &run : seg.label_sets) add(run.labels);
            for (const Event &ev : seg.events) {
                if (ev.type == EventType::BLOCKED || ev.type == EventType::HAIRPIN)
                    add(ev.labels);
            }
        }
        for (const Split &split : arm.splits) {
            for (const SplitBranch &b : split.branches) add(b.labels);
        }
        for (const PathResult &p : arm.paths) {
            if (p.continuation)
                add(p.continuation->labels);
        }
    }
    return out;
}

} // namespace

// Finding 1: the re-minimisation rounds were counted while the head was PLANNED, so a head
// the budget then refused still made its arm state greedy_losses (and label_evidence
// lower_bound) for a step that was never taken. Two fixtures, both re-minimising twice at
// the right root (A continues on both successors and is excluded by the branch limit 0;
// the derivation is re-run and the switches it then prices cascade):
//   - the review's repro (X+P+Y: A, X+Q+Z: A, X+P+Y: B; constant cost 1, loss budget 1):
//     every lineage is excluded and the root ends there (no successor followed);
//   - a third label C on P that no switch reaches (a table cost with only B -> A), which
//     goes on along P: the re-minimising head is followed.
// Refused (by a hook or by the memory budget) the arm counts no round and states no greedy
// loss. A cap at that head counts them as it always did (a result without a budget is
// unchanged), and a committed head counts them.
TEST(Stage2Review, RefusedHeadStatesNoGreedyLosses) {
    std::vector<std::string> b;
    for (uint32_t seed = 4; ; ++seed) {
        b = review_blocks({ 30, 25, 30, 25, 30, 30 }, seed);
        if (b[1][0] != b[3][0])
            break;
    }
    const std::string &X = b[0], &P = b[1], &Y = b[2], &Q = b[3], &Z = b[4];
    struct Case {
        std::string name;
        std::vector<std::string> records, labels, seed_labels;
        LabelChangeCost cost;
        bool followed;           // the re-minimising root follows a successor
    };
    const std::vector<Case> cases {
        { "the review's repro", { X + P + Y, X + Q + Z, X + P + Y }, { "A", "A", "B" },
          { "A", "B" }, LabelChangeCost::constant(1), false },
        // request indices: A = 0, B = 1, C = 2; every other switch is forbidden
        { "a followed head", { X + P + Y, X + Q + Z, X + P + Y, X + P + Y },
          { "A", "A", "B", "C" }, { "A", "B", "C" },
          LabelChangeCost::table({ { { 1, 0 }, 1.0 } }, kInfiniteLoss), true },
    };
    for (const Case &k : cases) {
        auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                11, k.records, k.labels, DeBruijnGraph::BASIC);
        Strategy st;
        st.direction = Strategy::RIGHT;
        st.loss_budget = 1;
        st.delivery = delivery_costs("full", true);
        const Seed seed = seed_of(X, k.seed_labels);
        LabelOracle oracle(*anno);

        std::vector<Admission> admissions;
        WalkerHooks record;
        record.deny = [&](const Admission &a) { admissions.push_back(a); return false; };
        const SeedResult base = traverse_seed(oracle, seed, st, k.cost, "", &record);
        ASSERT_EQ(2u, base.arms[1].max_reminimisation_rounds) << k.name << ": the fixture re-minimises twice";
        ASSERT_FALSE(admissions.empty()) << k.name;
        ASSERT_EQ(0u, admissions[0].at_bp) << k.name;
        ASSERT_EQ(k.followed, admissions[0].followed > 0) << k.name;
        EXPECT_EQ(k.followed, base.arms[1].steps > 0) << k.name;
        {
            const Json::Value j = seed_result_to_json(base, st, "full", false);
            EXPECT_TRUE(states(j["arms"]["right"]["limitations"], "greedy_losses")) << k.name;
            EXPECT_EQ("lower_bound", j["outcome"]["label_evidence"].asString()) << k.name;
        }

        auto refused_at_the_root = [&](const SeedResult &r, const Strategy &s, const std::string &what) {
            const ArmResult &arm = r.arms[1];
            ASSERT_TRUE(r.resource_stop) << what;
            EXPECT_EQ(0u, arm.steps) << what;
            EXPECT_EQ(0u, arm.complete_to_bp) << what;
            EXPECT_TRUE(arm.branch_events.empty()) << what;
            EXPECT_EQ(0u, arm.reminimisation_rounds) << what;
            EXPECT_EQ(0u, arm.max_reminimisation_rounds) << what;
            const Json::Value j = seed_result_to_json(r, s, "full", false);
            EXPECT_FALSE(states(j["arms"]["right"]["limitations"], "greedy_losses")) << what;
            EXPECT_EQ("complete", j["outcome"]["label_evidence"].asString()) << what;
            EXPECT_EQ("partial", j["outcome"]["walks"].asString()) << what;
            check_serialised(r, seed, s, oracle, what);
        };
        WalkerHooks deny;
        deny.deny = [](const Admission &a) { return a.ordinal == 0; };
        refused_at_the_root(traverse_seed(oracle, seed, st, k.cost, "", &deny), st,
                            k.name + ": denied by the hook");

        // the memory budget: one byte short of what admitting the root's step needs, and
        // enough for the depth-0 state
        Strategy tight = st;
        tight.max_memory_bytes = budget_for(need_of(admissions[0])) - 1;
        ASSERT_LT(usable(tight.max_memory_bytes), need_of(admissions[0])) << k.name;
        ASSERT_GE(usable(tight.max_memory_bytes), admissions[0].total) << k.name;
        refused_at_the_root(traverse_seed(oracle, seed, tight, k.cost), tight,
                            k.name + ": refused by the budget");

        if (!k.followed)
            continue;
        // a cap at the same head (a step cap below the successors it follows): counted, as
        // before stage 2, and stated
        Strategy capped = st;
        capped.max_steps = admissions[0].followed - 1;
        const SeedResult c = traverse_seed(oracle, seed, capped, k.cost);
        EXPECT_FALSE(c.resource_stop) << k.name;
        EXPECT_EQ(0u, c.arms[1].steps) << k.name;
        ASSERT_TRUE(c.arms[1].cap_trigger) << k.name;
        EXPECT_EQ(EndReason::MAX_STEPS, c.arms[1].cap_trigger->reason) << k.name;
        EXPECT_EQ(2u, c.arms[1].reminimisation_rounds) << k.name;
        EXPECT_EQ(2u, c.arms[1].max_reminimisation_rounds) << k.name;
        EXPECT_TRUE(states(seed_result_to_json(c, capped, "full", false)["arms"]["right"]["limitations"],
                           "greedy_losses")) << k.name;
        // the root committed and a later head refused: the root's rounds are counted
        if (admissions.size() > 1) {
            WalkerHooks later;
            later.deny = [](const Admission &a) { return a.ordinal == 1; };
            const SeedResult r = traverse_seed(oracle, seed, st, k.cost, "", &later);
            ASSERT_TRUE(r.resource_stop) << k.name;
            EXPECT_GT(r.arms[1].steps, 0u) << k.name;
            EXPECT_EQ(2u, r.arms[1].max_reminimisation_rounds) << k.name;
            check_serialised(r, seed, st, oracle, k.name + ": denied after the root");
        }
    }
}

// Finding 2: annotate mode names a label when a level's fetch returns it, before the
// heads of the level are admitted, so a walk stopped after the fetch held dictionary
// labels (L records, label_summary) that no segment, split or event recorded. Fixture: D
// is met only at depth 21; every refused admission, and every work budget, leaves a
// dictionary of exactly the recorded labels, in the unbudgeted walk's order. A cap keeps
// its dictionary as it was before stage 2 (a result without a budget is unchanged).
TEST(Stage2Review, StoppedAnnotateWalkNamesOnlyRecordedLabels) {
    const auto b = review_blocks({ 30, 40, 40 }, 777);
    const std::string &X = b[0], &P = b[1], &R = b[2];
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, { X + P, P.substr(10) + R }, { "A", "D" }, DeBruijnGraph::BASIC);
    const Seed seed = seed_of(X);
    const Strategy st = strategy_of(R"({"exhaustive": true, "direction": "right",
        "labels": {"mode": "annotate"}, "bounds": {"max_extension_bp": 200}})");
    LabelOracle oracle(*anno);
    std::vector<Admission> admissions;
    WalkerHooks record;
    record.deny = [&](const Admission &a) { admissions.push_back(a); return false; };
    const SeedResult base = traverse_seed(oracle, seed, st, LabelChangeCost::forbid(), "", &record);
    ASSERT_EQ(2u, base.label_dict.size());
    ASSERT_EQ(2u, recorded_labels(base).size());

    size_t stops = 0, cut = 0;
    auto check = [&](const SeedResult &r, const Strategy &s, const std::string &what) {
        ASSERT_TRUE(r.resource_stop) << what;
        stops++;
        const std::set<LabelId> recorded = recorded_labels(r);
        EXPECT_EQ(r.label_dict.size(), recorded.size()) << what << ": a dictionary label recorded nowhere";
        for (LabelId l : recorded) {
            ASSERT_LT(l, r.label_dict.size()) << what;
        }
        ASSERT_EQ(r.label_dict.size(), r.label_summary.size()) << what;
        // the walk's first-seen order: a prefix of the unbudgeted dictionary
        ASSERT_LE(r.label_dict.size(), base.label_dict.size()) << what;
        for (size_t l = 0; l < r.label_dict.size(); ++l) {
            EXPECT_EQ(base.label_dict[l].name, r.label_dict[l].name) << what;
        }
        cut += r.label_dict.size() < base.label_dict.size();
        if (stops % 8 == 1)
            check_serialised(r, seed, s, oracle, what);
    };
    for (const Admission &a : admissions) {
        WalkerHooks deny;
        deny.deny = [&](const Admission &x) { return x.ordinal == a.ordinal; };
        check(traverse_seed(oracle, seed, st, LabelChangeCost::forbid(), "", &deny), st,
              "denied #" + std::to_string(a.ordinal) + " at " + std::to_string(a.at_bp));
    }
    for (uint64_t w = 1; w <= base.account.work_used; w += 3) {
        Strategy s = st;
        s.max_work_units = w;
        const SeedResult r = traverse_seed(oracle, seed, s, LabelChangeCost::forbid());
        if (r.resource_stop)
            check(r, s, "max_work_units " + std::to_string(w));
    }
    EXPECT_GT(stops, 20u);
    EXPECT_GT(cut, 0u) << "a stop that met D before recording it";

    // the cap at the head that meets D: as before stage 2
    Strategy capped = st;
    capped.max_steps = 20;
    const SeedResult c = traverse_seed(oracle, seed, capped, LabelChangeCost::forbid());
    EXPECT_FALSE(c.resource_stop);
    EXPECT_EQ(base.label_dict.size(), c.label_dict.size());
}

// Finding 3: the depth-0 state (the dictionary and both roots, each reserved with what
// ending and delivering it costs) was never admitted: a seed carried by many labels was
// delivered whole under any budget (44 MiB under 1 MiB in the review), and the excess was
// reported as memory_bound_soft. Now a budget that does not hold it fails the seed, before
// anything per label is delivered; one that holds it admits it; and the account never
// exceeds the budget, at a level's start (its bins) included.
TEST(Stage2Review, DepthZeroStateIsAdmitted) {
    const std::string X = random_seq(30, 7101), P = random_seq(40, 7102);
    std::vector<std::string> seqs, labels;
    for (size_t i = 0; i < 200; ++i) {
        seqs.push_back(X + P);
        labels.push_back("L" + std::to_string(i));
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, seqs, labels, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    const Seed seed = seed_of(X);       // derived: all 200 labels
    Strategy st;
    st.delivery = delivery_costs("full", true);
    st.max_extension_bp = 20;

    // what the depth-0 state needs: the account at the first admission of a free walk
    std::vector<Admission> admissions;
    WalkerHooks record;
    record.deny = [&](const Admission &a) { admissions.push_back(a); return false; };
    traverse_seed(oracle, seed, st, LabelChangeCost::forbid(), "", &record);
    ASSERT_FALSE(admissions.empty());
    const uint64_t depth0 = admissions[0].total;
    ASSERT_GT(usable(uint64_t(1) << 20), 0u);
    ASSERT_GT(depth0, usable(uint64_t(1) << 20)) << "200 labels in detail full do not fit 1 MiB";

    st.max_memory_bytes = uint64_t(1) << 20;
    try {
        traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        FAIL() << "a depth-0 state over the budget was walked";
    } catch (const SeedBudgetError &e) {
        EXPECT_EQ(ResourceStop::MEMORY, e.stop().resource);
        EXPECT_EQ(0u, e.stop().at_bp);
        EXPECT_EQ(static_cast<double>(st.max_memory_bytes), e.stop().limit);
        // used: what the account held before any label (the seed, the fixed allotments)
        EXPECT_LE(e.stop().used, e.stop().limit);
        EXPECT_GT(e.stop().demand, e.stop().limit);
        EXPECT_EQ(200u, e.labels());
        EXPECT_TRUE(e.labels_from_seed());
        EXPECT_NE(std::string::npos, std::string(e.what()).find("bounds.max_memory_mb")) << e.what();
    }
    // exactly enough for it: admitted, and the account stays within the budget; one byte
    // less: refused
    for (uint64_t budget : { budget_for(depth0), budget_for(depth0) - 1 }) {
        st.max_memory_bytes = budget;
        const std::string what = "budget " + std::to_string(budget);
        if (usable(budget) >= depth0) {
            const SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
            EXPECT_LE(r.account.memory_peak, budget) << what;
            EXPECT_LE(r.account.memory_final, budget) << what;
            check_serialised(r, seed, st, oracle, what);
        } else {
            EXPECT_THROW(traverse_seed(oracle, seed, st, LabelChangeCost::forbid()), SeedBudgetError) << what;
        }
    }

    // the admitted account never exceeds a budget (constrain mode: nothing soft is charged
    // to it), at every level's start included, for budgets around every admission's need
    size_t runs = 0, failed = 0;
    for (const BudgetCase &c : budget_cases()) {
        if (c.st.label_mode != LabelMode::CONSTRAIN)
            continue;
        Strategy s = c.st;
        s.delivery = delivery_costs("full", true);
        std::vector<Admission> adm;
        WalkerHooks rec;
        rec.deny = [&](const Admission &a) { adm.push_back(a); return false; };
        run_case(c, s, &rec);
        std::vector<uint64_t> needs { adm[0].total };
        for (size_t i = 0; i < adm.size(); i += std::max<size_t>(1, adm.size() / 12)) {
            needs.push_back(need_of(adm[i]));
        }
        for (uint64_t need : needs) {
            for (uint64_t budget : { budget_for(need), budget_for(need) - 1 }) {
                s.max_memory_bytes = budget;
                try {
                    const SeedResult r = run_case(c, s);
                    EXPECT_LE(r.account.memory_peak, budget) << c.name << " budget " << budget;
                    EXPECT_LE(r.account.memory_final, budget) << c.name << " budget " << budget;
                    runs++;
                } catch (const SeedBudgetError &e) {
                    EXPECT_LT(usable(budget), adm[0].total) << c.name << ": " << e.what();
                    failed++;
                }
            }
        }
    }
    EXPECT_GT(runs, 50u);
    EXPECT_GT(failed, 0u);
}

// Finding 3, as a client sees it, and finding 7: the failed seed has the shape of a
// failed derivation (no arms, an error, outcome.walks failed) with its resource_stop and
// a seed-level walk_domain on the budget; memory_bound_soft does not report the depth-0
// state's own excess. The budget is per seed: the other seeds of the request are walked
// and each result equals the seed's alone (no request-level ledger).
TEST(Stage2Review, FailedSeedIsStatedPerSeed) {
    const std::string X = random_seq(30, 7101), P = random_seq(40, 7102), W = random_seq(60, 7103);
    std::vector<std::string> seqs, labels;
    for (size_t i = 0; i < 200; ++i) {
        seqs.push_back(X + P);
        labels.push_back("L" + std::to_string(i));
    }
    seqs.push_back(W);
    labels.push_back("solo");
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, seqs, labels, DeBruijnGraph::BASIC);
    auto request = [&](const std::vector<std::string> &which, const std::string &detail) {
        Json::Value r;
        for (const std::string &name : which) {
            Json::Value s;
            s["seed_id"] = name;
            if (name == "wide") {
                s["sequence"] = X;                     // derived: 200 labels
            } else {
                s["sequence"] = W.substr(name == "solo1" ? 0 : 20, 25);
                s["labels"].append("solo");
            }
            r["seeds"].append(s);
        }
        r["strategy"] = parse_json(R"({"bounds": {"max_extension_bp": 20, "max_memory_mb": 1},
            "output": {"timing": false}})");
        r["strategy"]["output"]["detail"] = detail;
        return process_traverse_request(r, *anno, "");
    };
    const Json::Value out = request({ "solo1", "wide", "solo2" }, "full");
    ASSERT_EQ(3u, out["results"].size());
    const Json::Value &f = out["results"][1];
    EXPECT_FALSE(f.isMember("arms"));
    EXPECT_TRUE(f.isMember("error"));
    EXPECT_EQ("wide", f["seed"]["seed_id"].asString());
    EXPECT_TRUE(f["seed"]["labels_from_seed"].asBool());
    EXPECT_EQ("failed", f["outcome"]["walks"].asString());
    EXPECT_EQ("complete", f["outcome"]["branch_diagnostics"].asString());
    EXPECT_EQ("complete", f["outcome"]["label_evidence"].asString());
    const Json::Value &q = f["resource_stop"];
    EXPECT_EQ("memory", q["resource"].asString());
    EXPECT_EQ("locus", q["scope"].asString());
    EXPECT_EQ("traversal", q["phase"].asString());
    EXPECT_EQ(1u, q["requested"].asUInt64());
    EXPECT_LE(q["used"].asUInt64(), 1u);
    std::set<std::string> actions;
    for (const Json::Value &a : q["actions"]) actions.insert(a.asString());
    EXPECT_TRUE(actions.count("raise_memory_budget"));
    EXPECT_TRUE(actions.count("use_graphlet"));
    EXPECT_TRUE(actions.count("lower_max_seed_labels"));
    EXPECT_FALSE(actions.count("continue_from_leaves")) << "a failed seed has no leaves";
    const Json::Value *wd = limitation_of(f["limitations"], "walk_domain");
    ASSERT_TRUE(wd);
    EXPECT_EQ("bounds.max_memory_mb", (*wd)["knob"].asString());
    EXPECT_EQ(1u, (*wd)["limit"].asUInt64());
    EXPECT_GT((*wd)["observed"].asUInt64(), 1u) << "the MiB the depth-0 state needs";
    const Json::Value *soft = limitation_of(f["limitations"], "memory_bound_soft");
    ASSERT_TRUE(soft) << "every response under a memory budget states it";
    EXPECT_EQ(0u, (*soft)["observed"].asUInt64()) << "not the depth-0 state's own excess";
    for (Json::ArrayIndex i : { 0u, 2u }) {
        EXPECT_TRUE(out["results"][i].isMember("arms")) << i;
        EXPECT_NE("failed", out["results"][i]["outcome"]["walks"].asString()) << i;
    }
    // per seed: each result is the seed's alone
    const char *names[] = { "solo1", "wide", "solo2" };
    for (Json::ArrayIndex i = 0; i < 3; ++i) {
        EXPECT_EQ(compact(out["results"][i]), compact(request({ names[i] }, "full")["results"][0]))
            << names[i];
    }
    // detail graphlet: no body for the failed seed, the others are delivered
    const Json::Value g = request({ "solo1", "wide", "solo2" }, "graphlet");
    EXPECT_FALSE(g["results"][1].isMember("graphlet"));
    EXPECT_EQ("failed", g["results"][1]["outcome"]["walks"].asString());
    EXPECT_TRUE(g["results"][0].isMember("graphlet"));
    EXPECT_TRUE(g["results"][2].isMember("graphlet"));
}

// Finding 3 in annotate mode: no permitted set, so the depth-0 state is the labels recorded
// at the roots, which labels.max_labels_per_node caps — the lever the failed seed names
// (not the derived set's cap, nor a named list, which annotate mode does not have).
TEST(Stage2Review, FailedAnnotateSeedNamesItsLever) {
    const std::string X = random_seq(30, 7301), P = random_seq(40, 7302);
    std::vector<std::string> seqs, labels;
    for (size_t i = 0; i < 2000; ++i) {
        seqs.push_back(X + P);
        labels.push_back("L" + std::to_string(i));
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, seqs, labels, DeBruijnGraph::BASIC);
    auto request = [&](uint64_t cap) {
        Json::Value r;
        r["seeds"][0]["sequence"] = X;
        r["strategy"] = parse_json(R"({"labels": {"mode": "annotate"}, "bounds":
            {"max_extension_bp": 20, "max_memory_mb": 1}, "output": {"detail": "full",
            "timing": false}})");
        r["strategy"]["labels"]["max_labels_per_node"] = Json::UInt64(cap);
        return process_traverse_request(r, *anno, "")["results"][0];
    };
    const Json::Value f = request(100'000);
    ASSERT_EQ("failed", f["outcome"]["walks"].asString()) << f.toStyledString();
    EXPECT_FALSE(f.isMember("arms"));
    std::set<std::string> actions;
    for (const Json::Value &a : f["resource_stop"]["actions"]) actions.insert(a.asString());
    EXPECT_TRUE(actions.count("lower_max_labels_per_node"));
    EXPECT_FALSE(actions.count("lower_max_seed_labels"));
    EXPECT_FALSE(actions.count("name_fewer_labels"));
    const Json::Value *wd = limitation_of(f["limitations"], "walk_domain");
    ASSERT_TRUE(wd);
    EXPECT_EQ("bounds.max_memory_mb", (*wd)["knob"].asString());
    // the lever works: a cap of 4 labels per node fits (the default 64 does not since the
    // delivery model prices every object at its demonstrated upper bound: a recorded label
    // costs some 8 KiB in detail full, its root reservations as much again)
    const Json::Value ok = request(4);
    EXPECT_TRUE(ok.isMember("arms")) << ok.toStyledString();
    EXPECT_NE("failed", ok["outcome"]["walks"].asString());
}

// Finding 4: plan_cost repeated commit_entries' scan of the runs a split had already
// taken (std::find per entry, O(|σ|²) per split) on every head, with no budget set. Both
// are O(1) per entry now (stamps by run id). Fixture: 300 labels on every branch of
// three splits, so that every run continues on both children of each split: the second
// child clones each one, and every run belongs to exactly one path (the plan's count of
// the clones is checked against the commit in debug builds).
TEST(Stage2Review, SplitClonesEveryRunContinuedTwice) {
    const auto b = review_blocks({ 30, 20, 20, 20, 20, 20, 20 }, 31);
    const std::string &X = b[0];
    // X -> {P, Q}, P -> {U, V}, Q -> {U', V'}: the branches start with distinct bases
    std::string P = b[1], Q = b[2], U = b[3], V = b[4], U2 = b[5], V2 = b[6];
    auto distinct_first = [](std::string &x, const std::string &y) {
        if (x[0] == y[0])
            x[0] = x[0] == 'A' ? 'C' : 'A';
    };
    distinct_first(Q, P);
    distinct_first(V, U);
    distinct_first(V2, U2);
    const size_t n = 300;
    std::vector<std::string> seqs, labels;
    for (size_t i = 0; i < n; ++i) {
        for (const std::string &s : { X + P + U, X + P + V, X + Q + U2, X + Q + V2 }) {
            seqs.push_back(s);
            labels.push_back("L" + std::to_string(i));
        }
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, seqs, labels, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    const Seed seed = seed_of(X);
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.exhaustive = true;
    st.max_label_branches = Strategy::kUnlimited;
    st.max_splits_per_path = Strategy::kUnlimited;
    st.merge_reconverge = false;
    st.max_extension_bp = 60;
    const SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    ASSERT_EQ(n, r.num_seed_labels);
    const ArmResult &arm = r.arms[1];
    ASSERT_EQ(3u, arm.splits.size());
    ASSERT_EQ(4u, arm.paths.size());
    // the roots' runs, and one clone per run at each split
    EXPECT_EQ(n + 3 * n, arm.runs.size());
    // every run is on exactly one path: no two leaves end the same run
    std::vector<size_t> owners(arm.runs.size(), 0);
    for (const PathResult &p : arm.paths) {
        EXPECT_EQ(n, p.end_labels.size());
        for (const LabelEnd &e : p.end_labels) owners[e.run]++;
    }
    for (uint32_t run = 0; run < arm.runs.size(); ++run) {
        // the first child of a split continues a run, the second its clone: no run ends
        // at a split, and each one at exactly one leaf
        EXPECT_EQ(1u, owners[run]) << "run " << run;
    }
    check_serialised(r, seed, st, oracle, "three wide splits");
}

// Finding 5 (stated, not charged): neither the work budget nor the deadline bounds the
// delivery of detail tree/full, whose paths spell every leaf's chain. Work stays the
// walk's own — the same walk, and the same stop, in every detail, which is what lets a
// graphlet rebuild the full response — and the size of the output (with the time to
// serialise it) is bounded by the memory budget, which charges the chains per expansion
// in the requested detail (CombTrieStaysLinear).
TEST(Stage2Review, WorkIsTheWalksInEveryDetail) {
    size_t stopped = 0;
    for (const BudgetCase &c : budget_cases()) {
        const SeedResult free = run_case(c, c.st);
        Strategy s = c.st;
        s.max_work_units = std::max<uint64_t>(1, free.account.work_used / 2);
        std::string reference;
        for (const char *detail : { "summary", "full", "graphlet" }) {
            s.delivery = delivery_costs(detail, s.sequences);
            const SeedResult r = run_case(c, s);
            stopped += r.resource_stop.has_value();
            const std::string j = compact(seed_result_to_json(r, s, "full", false));
            if (reference.empty()) {
                reference = j;
            } else {
                EXPECT_EQ(reference, j) << c.name << " " << detail;
            }
        }
    }
    EXPECT_GT(stopped, 10u);
}

// Finding 6: the seed phase (one annotation row per seed k-mer: the validation of an
// explicit set, or the derivation of one) was not charged as work, so a tiny work budget
// still read a whole long seed before its first head. It is charged now (8 per k-mer, 1
// per entry and coordinate) and checked like the walk, once every kWorkCheckInterval
// units: a seed phase within one interval is followed by a stop at the first head (a
// valid result complete to 0 bp); a longer one is cut within about the interval and the
// seed fails, with its resource_stop.
TEST(Stage2Review, SeedPhaseIsChargedAsWork) {
    const std::string L = random_seq(40'000, 7201);
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, { L }, { "A" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.max_extension_bp = 10;

    // unbudgeted: the seed phase is accounted beside the arms
    const Seed shortseed = seed_of(L.substr(100, 30), { "A" });
    const SeedResult free = traverse_seed(oracle, shortseed, st, LabelChangeCost::forbid());
    EXPECT_EQ(20u * (8 + 1), free.account.work_seed);
    EXPECT_EQ(free.account.work_seed + free.arms[0].work_units + free.arms[1].work_units,
              free.account.work_used);

    // within one interval: charged, and the first head stops the walk (a valid result)
    st.max_work_units = 50;
    const SeedResult within = traverse_seed(oracle, shortseed, st, LabelChangeCost::forbid());
    ASSERT_TRUE(within.resource_stop);
    EXPECT_EQ(ResourceStop::WORK, within.resource_stop->resource);
    EXPECT_EQ(0u, within.arms[1].complete_to_bp);
    EXPECT_GE(within.resource_stop->used, 180.0);
    check_serialised(within, shortseed, st, oracle, "seed phase within an interval");

    // a long seed, explicit and derived: cut within about two intervals
    for (bool derived : { false, true }) {
        const Seed longseed = seed_of(L.substr(0, 35'000), derived ? std::vector<std::string>()
                                                                   : std::vector<std::string> { "A" });
        const std::string what = derived ? "derived" : "explicit";
        try {
            traverse_seed(oracle, longseed, st, LabelChangeCost::forbid());
            FAIL() << what << ": a long seed phase ran past the work budget";
        } catch (const SeedBudgetError &e) {
            EXPECT_EQ(ResourceStop::WORK, e.stop().resource) << what;
            EXPECT_GT(e.stop().used, 50.0) << what;
            // compared every interval and cut at the first comparison that finds it W past the
            // budget (a seed phase less than W past it is let finish: review of stage 3, F9),
            // within one more chunk or row of it
            EXPECT_GE(e.stop().used - 50.0, static_cast<double>(kWorkCheckInterval)) << what;
            EXPECT_LE(e.stop().used, 2.0 * kWorkCheckInterval + 50 + 64) << what;
            EXPECT_EQ(e.stop().used, static_cast<double>(e.account().work_seed)) << what;
            EXPECT_EQ(derived, e.labels_from_seed()) << what;
        }
        // as a client sees it: a failed seed naming the work budget
        Json::Value r;
        r["seeds"][0]["sequence"] = longseed.sequence;
        for (const std::string &l : longseed.labels) r["seeds"][0]["labels"].append(l);
        r["strategy"] = parse_json(R"({"direction": "right", "bounds": {"max_extension_bp": 10,
            "max_work_units": 50}, "output": {"timing": false}})");
        const Json::Value res = process_traverse_request(r, *anno, "")["results"][0];
        EXPECT_EQ("failed", res["outcome"]["walks"].asString()) << what;
        EXPECT_EQ("work", res["resource_stop"]["resource"].asString()) << what;
        std::set<std::string> actions;
        for (const Json::Value &a : res["resource_stop"]["actions"]) actions.insert(a.asString());
        EXPECT_TRUE(actions.count("raise_work_budget")) << what;
        EXPECT_TRUE(actions.count("shorten_seed")) << what;
        const Json::Value *wd = limitation_of(res["limitations"], "walk_domain");
        ASSERT_TRUE(wd) << what;
        EXPECT_EQ("bounds.max_work_units", (*wd)["knob"].asString()) << what;
        EXPECT_FALSE(states(res["limitations"], "memory_bound_soft")) << what;
    }
}


/*
 * The external re-review of stage 2 (GPT, of 4d332729 and the fixes, confirmed on
 * 278a53dd): one regression per finding of the C++ side, each failing before its fix.
 */
namespace {

// |n| sequences |record|, one column "F", headers L0 ... (or |first| then L0 ...): the
// reviewer's indexes in-process (build -k 3 --mode basic, coordinates, a CoordToHeader)
struct HeaderIndexCase {
    std::shared_ptr<graph::AnnotatedDBG> anno;
    std::unique_ptr<annot::CoordToHeader> cth;
};

HeaderIndexCase header_index(size_t k, const std::vector<std::pair<std::string, std::string>> &records) {
    std::vector<std::string> seqs, labels, headers;
    std::vector<uint64_t> offsets, num_kmers;
    uint64_t offset = 0;
    for (const auto &[name, seq] : records) {
        seqs.push_back(seq);
        labels.push_back("F");
        headers.push_back(name);
        offsets.push_back(offset);
        num_kmers.push_back(seq.size() - k + 1);
        offset += seq.size() - k + 1;
    }
    HeaderIndexCase c;
    c.anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            k, seqs, labels, DeBruijnGraph::BASIC, true, offsets);
    std::vector<std::vector<std::string>> h { headers };
    std::vector<std::vector<uint64_t>> n { num_kmers };
    c.cth = std::make_unique<annot::CoordToHeader>(std::move(h), std::move(n));
    return c;
}

} // namespace

// Finding 2: the work budget exceeded its advertised maximum overrun. The level's rows were
// charged and checked in fixed batches of keys, whatever their width, and the root's row
// before any check: 25,000 headers on AAACAAAGAAAT under a budget of 1 used 125,044 units.
// Now the walk compares the budget after every charge — a fetch call's rows are charged
// when it returns them, and near the budget a call reads one key — so a stop overruns by
// what was charged since the previous comparison at most: here one row, as wide as the
// index makes it, stated with its number (ResourceAccount::largest_charge) rather than
// promised to be below W.
TEST(Stage2ReviewRecheck, WorkStopsWithinOneRowOfTheBudget) {
    const size_t n = 3000;
    std::vector<std::pair<std::string, std::string>> records;
    for (size_t i = 0; i < n; ++i) {
        records.emplace_back("L" + std::to_string(i), "AAACAAAGAAAT");
    }
    const HeaderIndexCase idx = header_index(3, records);
    LabelOracle oracle(*idx.anno, idx.cth.get());

    // annotate: the root's row (n labels) is charged before the first level reads anything
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.label_mode = LabelMode::ANNOTATE;
    st.seed_label_kind = LabelKind::HEADER;
    st.max_labels_per_node = 1;
    st.max_label_branches = Strategy::kUnlimited;
    st.max_splits_per_path = Strategy::kUnlimited;
    st.min_live_labels = 0;
    st.max_extension_bp = 10;
    st.max_work_units = 1;
    const SeedResult a = traverse_seed(oracle, seed_of("AAA"), st, LabelChangeCost::forbid());
    ASSERT_TRUE(a.resource_stop);
    EXPECT_EQ(ResourceStop::WORK, a.resource_stop->resource);
    EXPECT_EQ(0u, a.arms[1].complete_to_bp);
    // the root's row, the one charge before the first comparison (annotate reads no row
    // in its seed phase)
    EXPECT_EQ(8u + n, a.account.largest_charge);
    EXPECT_LE(a.resource_stop->used - 1, static_cast<double>(a.account.largest_charge));
    EXPECT_LE(a.resource_stop->used, 1.0 + kWorkCheckInterval);

    // constrain, a derived set of n labels: a budget that admits the seed phase and the
    // first enumeration stops at the first row the head consumes, not after the level's
    st = Strategy();
    st.direction = Strategy::RIGHT;
    st.max_seed_labels = n;
    st.max_extension_bp = 10;
    const Seed seed = seed_of("AAA");
    const SeedResult free = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    ASSERT_FALSE(free.resource_stop);
    ASSERT_EQ(n, free.label_dict.size());
    st.max_work_units = free.account.work_seed + 4 + 1;
    const SeedResult c = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    ASSERT_TRUE(c.resource_stop);
    EXPECT_EQ(ResourceStop::WORK, c.resource_stop->resource);
    EXPECT_EQ(0u, c.arms[1].complete_to_bp);
    const double over = c.resource_stop->used - static_cast<double>(st.max_work_units);
    EXPECT_GT(over, 0.0);
    EXPECT_LE(over, static_cast<double>(8 + n)) << "more than one row past the budget";
    EXPECT_LE(over, static_cast<double>(c.account.largest_charge));
    // the stop states the bound, with the most the seed charged between two comparisons
    const Json::Value j = seed_result_to_json(c, st, "full", false);
    const std::string message = j["resource_stop"]["message"].asString();
    EXPECT_NE(std::string::npos, message.find("after every charge")) << message;
    EXPECT_NE(std::string::npos, message.find("between two comparisons: "
                                              + std::to_string(c.account.largest_charge) + " units"))
        << message;
    check_serialised(c, seed, st, oracle, "work stop at a row");
}

// Finding 3: an interrupted fetch hid its memory excess. The level's fetch was checked
// against the work budget between its chunks, and a stop there skipped the observation of
// what the fetch held: 100,000 headers on the next node, max_memory_mb 1, max_work_units 100
// reported memory_bound_soft observed 0 while the cache alone held 1.6 MB beyond its
// allotment. Now what a fetch holds is observed before any check after it can throw.
TEST(Stage2ReviewRecheck, InterruptedFetchStatesWhatItHeld) {
    const size_t n = 20'000;
    std::vector<std::pair<std::string, std::string>> records { { "root", "AAAC" } };
    for (size_t i = 0; i < n; ++i) {
        records.emplace_back("L" + std::to_string(i), "AAC");
    }
    const HeaderIndexCase idx = header_index(3, records);
    LabelOracle oracle(*idx.anno, idx.cth.get());
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.label_mode = LabelMode::ANNOTATE;
    st.seed_label_kind = LabelKind::HEADER;
    st.max_labels_per_node = 100'000;
    st.max_label_branches = Strategy::kUnlimited;
    st.max_splits_per_path = Strategy::kUnlimited;
    st.min_live_labels = 0;
    st.max_extension_bp = 10;
    st.max_memory_bytes = uint64_t(1) << 20;
    st.max_work_units = 100;
    st.delivery = delivery_costs("full", true);
    const SeedResult r = traverse_seed(oracle, seed_of("AAA"), st, LabelChangeCost::forbid());
    ASSERT_TRUE(r.resource_stop);
    EXPECT_EQ(0u, r.arms[1].complete_to_bp);
    // the cache's kept keys alone (n rows' worth beyond a 256 KiB allotment) are held
    EXPECT_GT(r.account.soft_overshoot, n * 16 - (uint64_t(1) << 18));
    const Json::Value j = seed_result_to_json(r, st, "full", false);
    const Json::Value *soft = limitation_of(j["limitations"], "memory_bound_soft");
    ASSERT_TRUE(soft);
    EXPECT_GE((*soft)["observed"].asUInt64(), 1u);

    // a stop by work alone, at the first head after the fetch, states the excess too
    st.max_memory_bytes = uint64_t(1) << 30;
    st.max_work_units = 30;
    const SeedResult w = traverse_seed(oracle, seed_of("AAA"), st, LabelChangeCost::forbid());
    ASSERT_TRUE(w.resource_stop);
    EXPECT_EQ(ResourceStop::WORK, w.resource_stop->resource);
}

// Finding 8: a seed failed in its derivation (no carrier, the time budget) omitted
// memory_bound_soft under max_memory_mb, which budget failures and UTF-8 refusals state
TEST(Stage2ReviewRecheck, FailedDerivationsStateTheSoftBound) {
    const std::string X = random_seq(30, 8101), Y = random_seq(30, 8102), Z = random_seq(30, 8103);
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, { X + Y, Y + Z }, { "A", "B" }, DeBruijnGraph::BASIC);
    for (const char *bounds : { R"("max_extension_bp": 5, "max_memory_mb": 1)",
                                R"("max_extension_bp": 5, "max_memory_mb": 1, "time_budget_ms": 0.000000001)",
                                R"("max_extension_bp": 5)" }) {
        Json::Value r;
        r["seeds"][0]["sequence"] = X + Y + Z;
        r["strategy"] = parse_json(std::string(R"({"bounds": {)") + bounds
                                   + R"(}, "output": {"detail": "graphlet", "timing": false}})");
        const Json::Value res = process_traverse_request(r, *anno, "")["results"][0];
        const bool budget = std::string(bounds).find("max_memory_mb") != std::string::npos;
        const bool time = std::string(bounds).find("time_budget_ms") != std::string::npos;
        // a time budget spent after the first k-mer delivers the set derived so far (D3): X's
        // k-mers carry A alone, which the whole seed (no carrier) would exclude — walks partial,
        // label evidence qualified (overstated); otherwise no carrier fails the seed
        ASSERT_EQ(time ? "partial" : "failed", res["outcome"]["walks"].asString()) << bounds;
        const Json::Value *d = limitation_of(res["limitations"], "derivation");
        ASSERT_TRUE(d) << bounds;
        EXPECT_EQ(time ? "time_budget" : "no_carrier", (*d)["cause"].asString()) << bounds;
        const Json::Value *soft = limitation_of(res["limitations"], "memory_bound_soft");
        EXPECT_EQ(budget, soft != nullptr) << bounds;
        if (soft) {
            EXPECT_EQ(1u, (*soft)["limit"].asUInt64()) << bounds;
            EXPECT_TRUE((*soft)["observed"].isUInt64()) << bounds;
        }
        // memory_bound_soft is in no outcome class
        EXPECT_EQ(time ? "qualified" : "complete", res["outcome"]["label_evidence"].asString())
            << bounds;
    }
}

// N1: a traversal complete to its radius reported a time stop. With time_budget_ms 0 the
// deadline has passed at every depth after 0; when every remaining head has reached the
// radius, nothing is censored (the heads end as max_extension_bp), so no stop is stated
TEST(Stage2ReviewRecheck, CompleteWalkIsNoTimeStop) {
    const std::string seq = random_seq(60, 8201);
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            15, { seq }, { "F" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    const Seed seed = seed_of(seq.substr(0, 40), { "F" });
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.merge_reconverge = false;
    st.max_extension_bp = 1;
    st.time_budget_ms = 0;
    st.max_work_units = 1'000'000;
    st.delivery = delivery_costs("graphlet", true);
    const SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    EXPECT_FALSE(r.resource_stop);
    EXPECT_EQ(ArmResult::COMPLETE, r.arms[1].status);
    EXPECT_EQ(1u, r.arms[1].complete_to_bp);
    const Json::Value j = seed_result_to_json(r, st, "graphlet", false);
    EXPECT_FALSE(j.isMember("resource_stop"));
    EXPECT_EQ("complete", j["outcome"]["walks"].asString());
    const std::string text = check_serialised(r, seed, st, oracle, "radius 1, time 0");
    EXPECT_EQ(std::string::npos, text.find("\nQ "));
    // below the radius the deadline does censor a head, and says so
    st.max_extension_bp = 5;
    const SeedResult cut = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    ASSERT_TRUE(cut.resource_stop);
    EXPECT_EQ(ResourceStop::TIME, cut.resource_stop->resource);
    EXPECT_EQ(1u, cut.arms[1].complete_to_bp);
    EXPECT_EQ(ArmResult::TRUNCATED, cut.arms[1].status);
}

// Part A: a refusal injected by the denial hook stays a memory stop (limit "unlimited"
// without a budget, so that it is representable in Q and K), but no budget refused the
// head, and the statement says so
TEST(Stage2ReviewRecheck, InjectedRefusalIsNotABudget) {
    const std::string seq = random_seq(60, 8301);
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            15, { seq }, { "F" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    const Seed seed = seed_of(seq.substr(0, 20), { "F" });
    WalkerHooks hooks;
    hooks.deny = [](const Admission &a) { return a.at_bp == 2; };
    for (uint64_t budget : { uint64_t(0), uint64_t(64) << 20 }) {
        Strategy st;
        st.direction = Strategy::RIGHT;
        st.max_extension_bp = 10;
        st.max_memory_bytes = budget;
        const SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid(), "", &hooks);
        ASSERT_TRUE(r.resource_stop);
        EXPECT_TRUE(r.resource_stop->injected);
        const Json::Value j = seed_result_to_json(r, st, "full", false);
        const Json::Value &q = j["resource_stop"];
        EXPECT_EQ("memory", q["resource"].asString());
        EXPECT_EQ(budget ? Json::Value(Json::UInt64(64)) : Json::Value("unlimited"), q["effective"]);
        const std::string message = q["message"].asString();
        EXPECT_NE(std::string::npos, message.find("injected")) << message;
        EXPECT_EQ(std::string::npos, message.find("the memory budget (")) << message;
        // no budget lever helps against it
        ASSERT_EQ(1u, q["actions"].size());
        EXPECT_EQ("continue_from_leaves", q["actions"][0].asString());
        const Json::Value *wd = limitation_of(j["arms"]["right"]["limitations"], "walk_domain");
        ASSERT_TRUE(wd);
        EXPECT_NE(std::string::npos, (*wd)["effect"].asString().find("injected"));
        check_serialised(r, seed, st, oracle, "injected refusal");
    }
    // a refusal by the budget itself still names the budget
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.max_extension_bp = 10;
    st.max_memory_bytes = uint64_t(1) << 20;
    st.delivery = delivery_costs("full", true);
    const SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    if (r.resource_stop) {
        EXPECT_FALSE(r.resource_stop->injected);
        const std::string message = seed_result_to_json(r, st, "full", false)["resource_stop"]["message"].asString();
        EXPECT_NE(std::string::npos, message.find("the memory budget (")) << message;
    }
}

// Finding 1, as the reviewer ran it (detail full, a header of 180,000 control characters,
// a 4-base record, max_memory_mb 2; detail graphlet, a header of 200,000 '%' in annotate
// mode): the delivered result exceeded what the budget had reserved for it, since a name
// byte was charged at four bytes but JSON writes a control character as six. Now names are
// charged as escaped: the seed either fails at depth 0 (nothing per label delivered) or its
// delivery peak stays within the budget.
TEST(Stage2ReviewRecheck, EscapedNamesStayWithinTheBudget) {
    for (const auto &[name, mode, detail, mb] : {
             std::make_tuple(std::string(180'000, '\x01'), "constrain", "full", 2),
             std::make_tuple(std::string(200'000, '%'), "annotate", "graphlet", 2),
             std::make_tuple(std::string(70'000, '%'), "annotate", "graphlet", 1),
             std::make_tuple(std::string(60'000, '\x01'), "constrain", "full", 2) }) {
        auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                3, { "AAAC" }, { name }, DeBruijnGraph::BASIC);
        Json::Value r;
        r["seeds"][0]["sequence"] = "AAA";
        r["strategy"] = parse_json(std::string(R"({"direction": "right", "labels": {"mode": ")") + mode
                + R"(", "seed_label_kind": "column", "max_labels_per_node": 1}, "bounds":
                {"max_extension_bp": 1}, "output": {"detail": ")" + detail + R"(", "timing": false}})");
        r["strategy"]["bounds"]["max_memory_mb"] = mb;
        const Json::Value res = process_traverse_request(r, *anno, "")["results"][0];
        const std::string what = std::string(mode) + " " + detail + " " + std::to_string(name.size());
        const uint64_t budget = uint64_t(mb) << 20;
        EXPECT_LE(compact(res).size(), budget) << what;
        if (res["outcome"]["walks"].asString() == "failed") {
            EXPECT_EQ(std::string::npos, compact(res).find(name.substr(0, 100))) << what;
            EXPECT_EQ("memory", res["resource_stop"]["resource"].asString()) << what;
            continue;
        }
        EXPECT_LE(delivery_peak(res), budget) << what;
    }
}


/*
 * The adversarial review of the fixes to the re-review (round 3): what the fixes left open.
 * One regression per finding, each failing before its fix.
 */
namespace {

// two columns "a" and "b" holding one record each, both under the header |name|: a derived
// header label that cannot be resubmitted (ambiguous_header)
HeaderIndexCase same_header_in_two_columns(const std::string &name, const std::string &seq) {
    HeaderIndexCase c;
    c.anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            3, { seq, seq }, { "a", "b" }, DeBruijnGraph::BASIC, true, { 0, 0 });
    std::vector<std::vector<std::string>> h { { name }, { name } };
    std::vector<std::vector<uint64_t>> n { { seq.size() - 2 }, { seq.size() - 2 } };
    c.cth = std::make_unique<annot::CoordToHeader>(std::move(h), std::move(n));
    return c;
}

// the work bound a stop states: used exceeds the budget by at most the most the seed charged
// between two comparisons, and the message gives that number
void expect_stated_work_bound(const SeedResult &r, const Strategy &st, const std::string &what) {
    ASSERT_TRUE(r.resource_stop) << what;
    ASSERT_EQ(ResourceStop::WORK, r.resource_stop->resource) << what;
    const double over = r.resource_stop->used - static_cast<double>(st.max_work_units);
    EXPECT_GT(over, 0.0) << what;
    EXPECT_LE(over, static_cast<double>(r.account.largest_charge)) << what;
    const std::string message = seed_result_to_json(r, st, "full", false)["resource_stop"]["message"].asString();
    EXPECT_NE(std::string::npos, message.find("between two comparisons: "
                                              + std::to_string(r.account.largest_charge) + " units"))
        << what << ": " << message;
}

} // namespace

// F1: a failed seed's result echoed an index-supplied header (ambiguous_header: in the error
// and as observed) and the request's seed_id, neither charged: one failed result of 2.16 MB
// under 1 MiB with memory_bound_soft observed 0. Now, under a memory budget, a long index
// name is echoed as a bounded prefix with its length and where it is, and what the echoed
// seed_id holds beyond the budget is stated as memory_bound_soft's observed excess.
TEST(Stage2ReviewRound3, FailedResultsEchoWithinWhatTheyState) {
    // an index-supplied header (the CLI's view of the same case: integration
    // test_stage2_failed_results_echo_within_what_they_state)
    const std::string name(180'000, '\x01');
    const HeaderIndexCase idx = same_header_in_two_columns(name, "AAAC");
    LabelOracle oracle(*idx.anno, idx.cth.get());
    for (uint64_t budget : { uint64_t(1) << 20, uint64_t(0) }) {
        Strategy st;
        st.direction = Strategy::RIGHT;
        st.seed_label_kind = LabelKind::HEADER;
        st.max_extension_bp = 1;
        st.max_memory_bytes = budget;
        try {
            traverse_seed(oracle, seed_of("AAA"), st, LabelChangeCost::forbid());
            FAIL() << "the derived header is ambiguous";
        } catch (const SeedDerivationError &e) {
            EXPECT_EQ(SeedDerivationError::AMBIGUOUS_HEADER, e.cause());
            if (!budget) {
                // without a budget the header is echoed whole, as before
                EXPECT_EQ(name, e.subject());
                continue;
            }
            EXPECT_LE(e.subject().size(), 512u);
            EXPECT_NE(std::string::npos, e.subject().find("180000 bytes; column "));
            EXPECT_LE(std::string(e.what()).size(), 1024u);
            EXPECT_EQ(std::string::npos, std::string(e.what()).find(name.substr(0, 1000)));
        }
    }

    // the request's seed_id: a seed failed at depth 0 because of it still echoes it, and
    // the excess is stated
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            3, { "AAACGTTGCA" }, { "r1" }, DeBruijnGraph::BASIC);
    std::string seed_id;
    for (size_t i = 0; i < 100'000; ++i) {
        seed_id += "\xF0\x9F\x98\x80";     // U+1F600: twelve bytes in JSON (a surrogate pair)
    }
    Json::Value r;
    r["seeds"][0]["sequence"] = "AAACGT";
    r["seeds"][0]["seed_id"] = seed_id;
    r["strategy"] = parse_json(R"({"direction": "both", "labels": {"mode": "constrain",
        "seed_label_kind": "column"}, "bounds": {"max_extension_bp": 4, "max_memory_mb": 1},
        "output": {"detail": "graphlet", "timing": false}})");
    const Json::Value res = process_traverse_request(r, *anno, "")["results"][0];
    ASSERT_EQ("failed", res["outcome"]["walks"].asString());
    const Json::Value *soft = limitation_of(res["limitations"], "memory_bound_soft");
    ASSERT_TRUE(soft);
    const uint64_t budget = uint64_t(1) << 20;
    ASSERT_GT(compact(res).size(), budget);
    EXPECT_GE((*soft)["observed"].asUInt64() << 20, compact(res).size() - budget)
        << "the failed result exceeds the budget by more than memory_bound_soft states";
}

// F2: the explicit-label validation charged each row (which can fail the seed) before it
// observed what the fetched rows held: a row of 400,000 coordinates under 1 MiB and a work
// budget of 100 reported memory_bound_soft observed 0. Now a fetch call's rows and the
// query's cache are observed before the first charge.
TEST(Stage2ReviewRound3, ValidationStatesWhatItHeldBeforeItFails) {
    std::string seq = "CGAT";
    for (size_t i = 0; i < 200'000; ++i) {
        seq += "GAT";
    }
    const HeaderIndexCase idx = header_index(3, { { "r1", seq } });
    LabelOracle oracle(*idx.anno, idx.cth.get());
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.support = Support::TRACE;
    st.merge_reconverge = false;
    st.max_extension_bp = 5;
    st.max_memory_bytes = uint64_t(1) << 20;
    st.max_work_units = 100;
    st.delivery = delivery_costs("summary", true);
    try {
        traverse_seed(oracle, seed_of("GATG", { "r1" }), st, LabelChangeCost::forbid());
        FAIL() << "the first row alone is past the work budget";
    } catch (const SeedBudgetError &e) {
        EXPECT_EQ(ResourceStop::WORK, e.stop().resource);
        // GAT's row holds 200,000 coordinates (1.6 MB) under a 1 MiB budget
        EXPECT_GT(e.account().soft_overshoot, 200'000u * sizeof(uint64_t) - (uint64_t(1) << 20));
    }
}

// F3: the fix to finding 2 charged a row only when a head consumed it, so the rows a cut
// level had fetched were decoded but never charged, and `used` understated the decoding.
// Every row a fetch returns is charged again: with every row n labels wide, the work charged
// is at least rows_requested x (8 + n) when the level's second head is refused after the
// fetch returned its row.
TEST(Stage2ReviewRound3, EveryFetchedRowIsCharged) {
    const size_t n = 1000;
    const std::string S = random_seq(40, 9101);
    std::string X = random_seq(30, 9102), Y = random_seq(30, 9103);
    X[0] = 'A';
    Y[0] = 'C';
    std::vector<std::string> seqs, labels;
    for (size_t c = 0; c < n; ++c) {
        for (const std::string *tail : { &X, &Y }) {
            seqs.push_back(S + *tail);
            labels.push_back("c" + std::to_string(c));
        }
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            15, seqs, labels, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    Strategy st;
    st.direction = Strategy::RIGHT;
    st.label_mode = LabelMode::ANNOTATE;
    st.seed_label_kind = LabelKind::COLUMN;
    st.max_labels_per_node = 1;
    st.max_label_branches = Strategy::kUnlimited;
    st.max_splits_per_path = Strategy::kUnlimited;
    st.min_live_labels = 0;
    st.merge_reconverge = false;
    st.max_extension_bp = 30;
    st.max_work_units = 1'000'000'000;
    // the seed's last k-mer is S[5, 20); the head at depth 20 is S's last k-mer, which
    // splits into X and Y: the level at depth 21 has two heads, the second one refused
    const uint64_t D = 21;
    size_t at_d = 0;
    WalkerHooks hooks;
    hooks.deny = [&](const Admission &a) { return a.at_bp == D && at_d++ == 1; };
    const SeedResult r = traverse_seed(oracle, seed_of(S.substr(0, 20)), st,
                                       LabelChangeCost::forbid(), "", &hooks);
    ASSERT_TRUE(r.resource_stop);
    EXPECT_TRUE(r.resource_stop->injected);
    EXPECT_EQ(D, r.resource_stop->at_bp);
    EXPECT_EQ(2u, at_d);
    const uint64_t rows = r.annotation_counters.rows_requested;
    ASSERT_GT(rows, D);
    EXPECT_GE(r.arms[1].work_units, rows * (8 + n))
        << "a row the fetch returned was not charged";
}

// F4, F5, F7: "one indivisible charge, the widest row" was broken by the two roots' rows of
// direction both (charged before the first comparison: 50,016 under a budget of 1, stated
// 25,008), by a row's trace coordinates (charged as one sum after the row: stated 9, overrun
// 99,924) and could be by a label-state scan (which carried no number). Now the coordinates
// are part of their row's charge, and every stop states the most its seed charged between
// two comparisons, which bounds its overrun whatever the charge was.
TEST(Stage2ReviewRound3, WorkStopsStateTheirLargestCharge) {
    {
        // F4: both arms' roots, as one charge with the (empty) annotate seed phase
        const size_t n = 3000;
        std::vector<std::pair<std::string, std::string>> records;
        for (size_t i = 0; i < n; ++i) {
            records.emplace_back("L" + std::to_string(i), "AAACAAAGAAAT");
        }
        const HeaderIndexCase idx = header_index(3, records);
        LabelOracle oracle(*idx.anno, idx.cth.get());
        Strategy st;
        st.label_mode = LabelMode::ANNOTATE;
        st.seed_label_kind = LabelKind::HEADER;
        st.max_labels_per_node = 1;
        st.max_label_branches = Strategy::kUnlimited;
        st.max_splits_per_path = Strategy::kUnlimited;
        st.min_live_labels = 0;
        st.max_extension_bp = 10;
        st.max_work_units = 1;
        for (const char *seed : { "AAA", "AAAC", "AAACAAAG" }) {
            const SeedResult r = traverse_seed(oracle, seed_of(seed), st, LabelChangeCost::forbid());
            expect_stated_work_bound(r, st, std::string("both roots of ") + seed);
            EXPECT_EQ(2 * (8u + n), r.account.largest_charge) << seed;
        }
    }
    {
        // F5: a row with one label and 20,000 coordinates under support: trace
        std::string seq = "CGAT";
        for (size_t i = 0; i < 20'000; ++i) {
            seq += "GAT";
        }
        const HeaderIndexCase idx = header_index(3, { { "r1", seq } });
        LabelOracle oracle(*idx.anno, idx.cth.get());
        Strategy st;
        st.direction = Strategy::RIGHT;
        st.support = Support::TRACE;
        st.merge_reconverge = false;
        st.seed_label_kind = LabelKind::HEADER;
        st.max_extension_bp = 5;
        st.max_work_units = 100;
        const SeedResult r = traverse_seed(oracle, seed_of("CGA"), st, LabelChangeCost::forbid());
        expect_stated_work_bound(r, st, "trace coordinates");
        EXPECT_GE(r.account.largest_charge, 20'000u);
    }
    {
        // F7: a sweep of budgets over a constrain walk of 1,000 derived labels (rows of
        // 1,000 entries, label-state scans of 2,000): whatever charge a stop comes after,
        // its overrun is within the stated number
        const size_t n = 1000;
        std::vector<std::pair<std::string, std::string>> records;
        for (size_t i = 0; i < n; ++i) {
            records.emplace_back("L" + std::to_string(i), "AAACAAAGAAAT");
        }
        const HeaderIndexCase idx = header_index(3, records);
        LabelOracle oracle(*idx.anno, idx.cth.get());
        Strategy st;
        st.max_seed_labels = n;
        st.max_extension_bp = 10;
        st.max_label_branches = Strategy::kUnlimited;
        const Seed seed = seed_of("AAACA");
        const SeedResult free = traverse_seed(oracle, seed, st, LabelChangeCost::constant(1));
        ASSERT_FALSE(free.resource_stop);
        size_t stops = 0;
        const uint64_t walk = free.account.work_used - free.account.work_seed;
        for (uint64_t extra = 0; extra < walk; extra += walk / 12 + 1) {
            st.max_work_units = free.account.work_seed + extra + 1;
            const SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::constant(1));
            if (!r.resource_stop)
                continue;
            ++stops;
            expect_stated_work_bound(r, st, "budget " + std::to_string(st.max_work_units));
        }
        EXPECT_GT(stops, 3u);
    }
}

// F6: the depth-0 state built its label dictionary (names in each LabelRef and again in the
// query) before its admission; a seed the admission failed reported memory_bound_soft 0 with
// twelve names of 1 MB held under 1 MiB. Now what the dictionary holds is observed first.
TEST(Stage2ReviewRound3, FailedAdmissionStatesTheDictionaryItHeld) {
    std::vector<std::pair<std::string, std::string>> records;
    for (size_t i = 0; i < 6; ++i) {
        records.emplace_back("N" + std::to_string(i) + std::string(200'000, 'x'), "AAAC");
    }
    const HeaderIndexCase idx = header_index(3, records);
    LabelOracle oracle(*idx.anno, idx.cth.get());
    for (LabelMode mode : { LabelMode::CONSTRAIN, LabelMode::ANNOTATE }) {
        Strategy st;
        st.direction = Strategy::RIGHT;
        st.label_mode = mode;
        st.seed_label_kind = LabelKind::HEADER;
        st.max_labels_per_node = 100;
        st.max_extension_bp = 1;
        st.max_memory_bytes = uint64_t(1) << 20;
        st.delivery = delivery_costs("graphlet", true);
        try {
            traverse_seed(oracle, seed_of("AAA"), st, LabelChangeCost::forbid());
            FAIL() << to_string(mode) << ": six names of 200 KB fit 1 MiB";
        } catch (const SeedBudgetError &e) {
            EXPECT_EQ(ResourceStop::MEMORY, e.stop().resource) << to_string(mode);
            // the names alone, held twice (the LabelRef and the query's or recorder's copy)
            EXPECT_GE(e.account().soft_overshoot, 2 * 6 * 200'000u - (uint64_t(1) << 20))
                << to_string(mode);
        }
    }
}


/*************** stage 3: the budget-aware annotation reads (§14.1) ***************/

namespace {

// Row-diff budget cases: the dense random graphs of budget_cases() and the switch chain,
// annotated by the row-diff pipeline (RowDiff<ColumnMajor>, whose reads are budget-aware),
// under the same strategies; and a trace case on a coordinate row-diff annotation
std::vector<BudgetCase> rowdiff_budget_cases() {
    std::vector<BudgetCase> cases;
    std::vector<DeBruijnGraph::Mode> modes { DeBruijnGraph::BASIC };
#if ! _PROTEIN_GRAPH
    modes.push_back(DeBruijnGraph::PRIMARY);
#endif
    for (auto mode : modes) {
        std::vector<std::string> seqs, labels;
        for (uint32_t i = 0; i < 6; ++i) {
            seqs.push_back("AAA" + random_seq(14, 11 + i));
            labels.push_back(std::string(1, "CDECDF"[i]));
        }
        std::shared_ptr<graph::AnnotatedDBG> anno
            = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(3, seqs, labels, mode);
        const std::string where = "row-diff mode " + std::to_string(mode);
        LabelChangeCost constant;
        cases.push_back({ where + " merge", anno, seed_of("AAA", { "C", "D", "E" }),
                          strategy_of(R"({"bounds": {"max_extension_bp": 7}})") });
        cases.push_back({ where + " exhaustive", anno, seed_of("AAA", { "C", "D", "E" }),
                          strategy_of(R"({"exhaustive": true, "bounds": {"max_extension_bp": 6}})") });
        Strategy sw = strategy_of(R"({"bounds": {"max_extension_bp": 8}, "labels": {"extra": ["E", "F"],
            "loss_budget": 2, "change_cost": {"model": "constant", "value": 1}},
            "branching": {"max_label_branches": 1}})", &constant);
        cases.push_back({ where + " switching", anno, seed_of("AAA", { "C", "D" }), sw, constant });
        cases.push_back({ where + " annotate", anno, seed_of("AAA"),
                          strategy_of(R"({"exhaustive": true, "labels": {"mode": "annotate"},
                                          "bounds": {"max_extension_bp": 5}})") });
        cases.push_back({ where + " annotate merge", anno, seed_of("AAA"),
                          strategy_of(R"({"labels": {"mode": "annotate"}, "branching": {"on_reconverge": "merge"},
                                          "bounds": {"max_extension_bp": 6}})") });
        cases.push_back({ where + " derived", anno, seed_of("AAA"),
                          strategy_of(R"({"labels": {"seed_label_kind": "column"},
                                          "bounds": {"max_extension_bp": 6}})") });
    }
    {
        std::string S, T1, T2, T3, T4;
        for (uint32_t t = 500; ; t += 5) {
            S = random_seq(30, t); T1 = random_seq(30, t + 1); T2 = random_seq(30, t + 2);
            T3 = random_seq(30, t + 3); T4 = random_seq(30, t + 4);
            if (T3[0] != T4[0])
                break;
        }
        std::shared_ptr<graph::AnnotatedDBG> anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
                11, { S + T1, T1.substr(10) + T2, T2.substr(10) + T3, T2.substr(10) + T4 },
                { "A", "B", "C", "D" }, DeBruijnGraph::BASIC);
        LabelChangeCost constant;
        Strategy st = strategy_of(R"({"direction": "right", "labels": {"extra": ["B", "C", "D"],
            "change_cost": {"model": "constant", "value": 1}, "loss_budget": 3},
            "branching": {"max_label_branches": 1}, "bounds": {"max_extension_bp": 80}})", &constant);
        cases.push_back({ "row-diff switch chain", anno, seed_of(S, { "A" }), st, constant });
        // the same, by coordinates (TupleRowDiff), traced
        std::shared_ptr<graph::AnnotatedDBG> coords = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
                11, { S + T1, T1.substr(10) + T2, T2.substr(10) + T3, T2.substr(10) + T4 },
                { "A", "B", "C", "D" }, DeBruijnGraph::BASIC, true);
        cases.push_back({ "row-diff coordinates trace", coords, seed_of(S, { "A", "B" }),
                          strategy_of(R"({"direction": "both", "support": "trace",
                                          "branching": {"on_reconverge": "keep"},
                                          "bounds": {"max_extension_bp": 60}})") });
        cases.push_back({ "row-diff coordinates annotate", coords, seed_of(S),
                          strategy_of(R"({"labels": {"mode": "annotate", "seed_label_kind": "column"},
                                          "bounds": {"max_extension_bp": 40}})") });
    }
    return cases;
}

// §14 freeze gate, "a row-diff row whose dependencies are dense but whose result is tiny":
// along the path P (label A) the node u starts 300 records B0 .. B299, so the node v
// before it carries A alone while its row-diff successor carries 301 labels — v's row is
// {A}, read through a diff of 300 entries and the dense rows after it
struct DenseCase {
    std::shared_ptr<graph::AnnotatedDBG> anno;
    std::string P;
    size_t v = 30;          // v = P[v, v + k), u the k-mer after it
};
DenseCase dense_case() {
    DenseCase c;
    c.P = random_seq(80, 7000);
    std::vector<std::string> seqs { c.P }, labels { "A" };
    for (uint32_t i = 0; i < 300; ++i) {
        seqs.push_back(c.P.substr(c.v + 1, 11) + random_seq(15, 7100 + i));
        labels.push_back("B" + std::to_string(i));
    }
    c.anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
            11, seqs, labels, DeBruijnGraph::BASIC);
    return c;
}

// the result as the response states it, wall-clock values apart
std::string result_text(const SeedResult &r, const Strategy &st) {
    Json::StreamWriterBuilder w;
    w["indentation"] = "";
    return Json::writeString(w, seed_result_to_json(r, st, "full", false));
}

} // namespace

// W1: a refusal injected at any charge of a budget-aware read. A read of several keys that
// does not fit is retried in smaller runs, so a refused charge of a multi-key run changes
// nothing (the result is the unrefused one); a key that does not fit alone is the stop: in
// the seed phase or at an annotate root it fails the seed (phase annotation_decode), in a
// level's fetch it censors the level from its first head — the result a valid prefix,
// consistent and serialisable, its stop stated (Q ... annotation_decode). In the lookahead
// a refusal changes nothing.
TEST(GraphletStage3Decode, DecodeDenialLeavesAConsistentPrefix) {
    size_t seed_failures = 0, level_stops = 0, warm_denials = 0, retried = 0;
    for (const BudgetCase &c : rowdiff_budget_cases()) {
        Strategy st = c.st;
        st.max_memory_bytes = uint64_t(1) << 30;    // never trips: the reads are budget-aware
        std::vector<DecodeCharge> charges;
        WalkerHooks record;
        record.deny_decode = [&](const DecodeCharge &d) { charges.push_back(d); return false; };
        const SeedResult base = run_case(c, st, &record);
        ASSERT_FALSE(base.resource_stop) << c.name;
        ASSERT_TRUE(base.account.decode_charged) << c.name;
        ASSERT_GT(charges.size(), 0u) << c.name;
        const std::string base_text = result_text(base, st);
        const size_t most = 60;
        LabelOracle oracle(*c.anno);
        for (size_t s = 0; s < std::min(charges.size(), most); ++s) {
            const size_t ordinal = charges.size() <= most ? s : s * (charges.size() - 1) / (most - 1);
            const DecodeCharge &denied = charges[ordinal];
            const std::string what = c.name + " decode charge #" + std::to_string(ordinal) + " ("
                                   + std::to_string(denied.where) + " " + to_string(denied.arm)
                                   + " at " + std::to_string(denied.at_bp) + ")";
            WalkerHooks deny;
            deny.deny_decode = [&](const DecodeCharge &d) { return d.ordinal == ordinal; };
            if (denied.where == DecodeCharge::SEED || denied.where == DecodeCharge::ROOT) {
                try {
                    const SeedResult r = run_case(c, st, &deny);
                    // the run was retried in smaller runs
                    EXPECT_EQ(base_text, result_text(r, st)) << what;
                    retried++;
                } catch (const SeedBudgetError &e) {
                    EXPECT_STREQ("annotation_decode", e.stop().phase) << what;
                    EXPECT_TRUE(e.stop().injected) << what;
                    EXPECT_EQ(ResourceStop::MEMORY, e.stop().resource) << what;
                    seed_failures++;
                }
                continue;
            }
            const SeedResult r = run_case(c, st, &deny);
            if (denied.where == DecodeCharge::WARM) {
                // the lookahead gives up silently: nothing depends on what it cached
                EXPECT_EQ(base_text, result_text(r, st)) << what;
                warm_denials++;
                continue;
            }
            if (!r.resource_stop) {
                EXPECT_EQ(base_text, result_text(r, st)) << what;
                retried++;
                continue;
            }
            level_stops++;
            EXPECT_STREQ("annotation_decode", r.resource_stop->phase) << what;
            EXPECT_TRUE(r.resource_stop->injected) << what;
            EXPECT_EQ(denied.arm, r.resource_stop->arm) << what;
            EXPECT_EQ(denied.at_bp, r.resource_stop->at_bp) << what;
            const ArmResult &arm = r.arms[static_cast<size_t>(denied.arm)];
            EXPECT_EQ(ArmResult::TRUNCATED, arm.status) << what;
            EXPECT_EQ(denied.at_bp, arm.complete_to_bp) << what;
            for (size_t a = 0; a < 2; ++a) {
                const ArmResult &ra = r.arms[a];
                if (!ra.requested)
                    continue;
                const std::string at = what + " arm " + to_string(ra.arm);
                check_bins(ra, at);
                EXPECT_EQ(walks_upto(base.arms[a], ra.complete_to_bp),
                          walks_upto(ra, ra.complete_to_bp)) << at;
                if (&ra != &arm && ra.status != ArmResult::COMPLETE) {
                    EXPECT_GE(ra.complete_to_bp, denied.at_bp) << at;
                    EXPECT_LE(ra.complete_to_bp, denied.at_bp + 1) << at;
                }
            }
            const std::string text = check_serialised(r, c.seed, st, oracle, what);
            EXPECT_NE(std::string::npos, text.find("\nQ locus memory annotation_decode ")) << what;
            const Json::Value j = seed_result_to_json(r, st, "summary", false);
            EXPECT_EQ("partial", j["outcome"]["walks"].asString()) << what;
            EXPECT_EQ("annotation_decode", j["resource_stop"]["phase"].asString()) << what;
            // an injected refusal names no budget's levers
            Json::Value only(Json::arrayValue);
            only.append("continue_from_leaves");
            EXPECT_EQ(only, j["resource_stop"]["actions"]) << what;
        }
    }
    EXPECT_GT(seed_failures, 20u);
    EXPECT_GT(level_stops, 100u);
    EXPECT_GT(warm_denials, 5u);
    EXPECT_GT(retried, 5u);
    std::cerr << "decode denials: " << seed_failures << " seed failures, " << level_stops
              << " level stops, " << warm_denials << " in the lookahead, " << retried
              << " retried in smaller runs" << std::endl;
}

// W3, W4: memory and work budgets swept across the stops, with annotation.batch_kmers from 1
// to 1000: the same result, the same stop (phase, used, remaining, depth) and the same
// failures, whatever the lookahead decoded
TEST(GraphletStage3Decode, StopsDoNotDependOnBatchKmers) {
    size_t compared = 0, decode_stops = 0, decode_failures = 0;
    for (const BudgetCase &c : rowdiff_budget_cases()) {
        Strategy big = c.st;
        big.max_memory_bytes = uint64_t(1) << 30;
        const SeedResult full = run_case(c, big);
        const uint64_t peak = full.account.memory_peak;
        const uint64_t work = full.account.work_used;
        std::vector<std::pair<uint64_t, uint64_t>> budgets;    // memory, work
        for (uint64_t m = 8192; m < 8 * peak; m = m * 5 / 4 + 1) budgets.emplace_back(m, 0);
        for (uint64_t w = 1; w < 2 * work; w = w * 3 / 2 + 1) budgets.emplace_back(0, w);
        for (const auto &[memory, units] : budgets) {
            std::string reference;
            for (size_t batch : { 1, 3, 64, 1000 }) {
                Strategy st = c.st;
                st.max_memory_bytes = memory;
                st.max_work_units = units;
                st.batch_kmers = batch;
                std::string text;
                try {
                    const SeedResult r = run_case(c, st);
                    text = result_text(r, st);
                    if (batch == 64 && r.resource_stop
                            && std::string(r.resource_stop->phase) == "annotation_decode")
                        decode_stops++;
                } catch (const SeedBudgetError &e) {
                    text = std::string("failed: ") + e.what() + " " + e.stop().phase + " "
                         + std::to_string(e.stop().used) + " " + std::to_string(e.stop().demand)
                         + " " + std::to_string(e.account().work_used)
                         + " " + std::to_string(e.account().soft_overshoot);
                    if (batch == 64 && std::string(e.stop().phase) == "annotation_decode")
                        decode_failures++;
                } catch (const SeedDerivationError &e) {
                    text = std::string("derivation: ") + e.what();
                }
                if (reference.empty()) {
                    reference = text;
                } else {
                    EXPECT_EQ(reference, text) << c.name << " memory " << memory << " work "
                                               << units << " batch_kmers " << batch;
                }
                compared++;
            }
        }
    }
    EXPECT_GT(compared, 1000u);
    EXPECT_GT(decode_stops, 5u);
    EXPECT_GT(decode_failures, 5u);
}

// W5 (§14 freeze gate), W6, W7: a row whose dependencies are dense but whose result is tiny
// stops the walk at the level that reads it (partial, phase annotation_decode, complete to
// that level) under a budget the walk otherwise holds; read in the seed phase — to validate
// the seed, to derive its labels, or as an annotate root — it fails the seed. The work of a
// row includes its dependency rows.
TEST(GraphletStage3Decode, DenseDependenciesTinyResult) {
    const DenseCase dc = dense_case();
    LabelOracle oracle(*dc.anno);
    ASSERT_TRUE(oracle.decode_charged());
    // the walk: right from P[0, 20) under label A, past v
    const Seed walk = seed_of(dc.P.substr(0, 20), { "A" });
    const Strategy base_st = strategy_of(R"({"direction": "right", "bounds": {"max_extension_bp": 40}})");
    const SeedResult unbudgeted = traverse_seed(oracle, walk, base_st, LabelChangeCost::forbid());
    EXPECT_EQ(ArmResult::COMPLETE, unbudgeted.arms[static_cast<size_t>(Arm::RIGHT)].status);
    // v is first read by the level at depth v + 1 - 20 (the successor keys of its heads)
    const uint64_t v_level = dc.v + 11 - 20;
    size_t stopped_at_v = 0;
    uint64_t smallest_complete = 0;
    for (uint64_t memory = 16384; memory < (uint64_t(4) << 20); memory = memory * 9 / 8 + 1) {
        Strategy st = base_st;
        st.max_memory_bytes = memory;
        SeedResult r;
        try {
            r = traverse_seed(oracle, walk, st, LabelChangeCost::forbid());
        } catch (const SeedBudgetError &) {
            continue;
        }
        if (!r.resource_stop) {
            if (!smallest_complete)
                smallest_complete = memory;
            continue;
        }
        if (std::string(r.resource_stop->phase) != "annotation_decode")
            continue;
        const ArmResult &arm = r.arms[static_cast<size_t>(Arm::RIGHT)];
        EXPECT_EQ(arm.complete_to_bp, r.resource_stop->at_bp);
        EXPECT_EQ(ArmResult::TRUNCATED, arm.status);
        stopped_at_v += r.resource_stop->at_bp == v_level || r.resource_stop->at_bp == v_level + 1;
        // the dense rows are read whole only under a budget that holds them: what a stop
        // leaves is a valid prefix, stated
        const Json::Value j = seed_result_to_json(r, st, "summary", false);
        EXPECT_EQ("partial", j["outcome"]["walks"].asString());
        check_serialised(r, walk, st, oracle, "dense " + std::to_string(memory));
    }
    EXPECT_GT(stopped_at_v, 0u) << "no budget stopped the walk at the dense row";
    ASSERT_GT(smallest_complete, 0u);
    // the work of the walk counts v's dependency rows (300 entries and more)
    {
        Strategy st = base_st;
        st.max_memory_bytes = smallest_complete;
        const SeedResult r = traverse_seed(oracle, walk, st, LabelChangeCost::forbid());
        EXPECT_GT(r.arms[static_cast<size_t>(Arm::RIGHT)].work_units,
                  unbudgeted.arms[static_cast<size_t>(Arm::RIGHT)].work_units + 300);
        EXPECT_TRUE(r.account.decode_charged);
    }
    // in the seed phase the same read fails the seed, with the seed's levers
    struct SeedPhase { std::string what; Seed seed; Strategy st; };
    std::vector<SeedPhase> phases {
        { "validation", seed_of(dc.P.substr(dc.v - 5, 20), { "A" }), base_st },
        { "derivation", seed_of(dc.P.substr(dc.v - 5, 20)),
          strategy_of(R"({"direction": "right", "labels": {"seed_label_kind": "column"},
                          "bounds": {"max_extension_bp": 40}})") },
        { "annotate root", seed_of(dc.P.substr(dc.v - 9, 20)),
          strategy_of(R"({"direction": "right", "labels": {"mode": "annotate"},
                          "bounds": {"max_extension_bp": 10}})") },
    };
    for (const SeedPhase &p : phases) {
        size_t failed = 0;
        for (uint64_t memory = 8192; memory < (uint64_t(4) << 20); memory = memory * 9 / 8 + 1) {
            Strategy st = p.st;
            st.max_memory_bytes = memory;
            try {
                traverse_seed(oracle, p.seed, st, LabelChangeCost::forbid());
            } catch (const SeedBudgetError &e) {
                if (std::string(e.stop().phase) != "annotation_decode")
                    continue;
                failed++;
                EXPECT_FALSE(e.stop().injected) << p.what;
                // what admitting the refused row needed, beside what was held: more than the
                // budget (review of stage 3, F7: it was the budget plus one byte)
                EXPECT_GT(e.stop().demand, memory) << p.what;
                EXPECT_EQ(ResourceStop::READ_ROW, e.stop().cause) << p.what;
                EXPECT_TRUE(e.account().decode_charged) << p.what;
            }
        }
        EXPECT_GT(failed, 0u) << p.what;
    }
}

// The statements a decode stop adds fit the widths the delivery bound prices (kMessage,
// kEffect), at the widest numbers, and a work stop's on every annotation format
TEST(GraphletStage3Decode, StatementsFitTheirWidths) {
    const uint64_t big = std::numeric_limits<uint64_t>::max() / 2;
    for (bool annotate : { false, true }) {
        // 0, 1: a level's row (its demand known), injected; 2 .. 4: work stops; 5: a row whose
        // read alone was refused; 6: the labels a row names (phase traversal); 7: the level's
        // own lists (phase traversal)
        for (int variant = 0; variant < 8; ++variant) {
            SeedResult r;
            Strategy st;
            st.label_mode = annotate ? LabelMode::ANNOTATE : LabelMode::CONSTRAIN;
            st.max_memory_bytes = big & ~((uint64_t(1) << 20) - 1);
            st.max_work_units = big;
            ResourceStop q;
            const bool memory = variant < 2 || variant > 4;
            q.resource = memory ? ResourceStop::MEMORY : ResourceStop::WORK;
            q.phase = variant < 2 || variant == 5 ? "annotation_decode" : "traversal";
            q.cause = variant < 2 || variant == 5 ? ResourceStop::READ_ROW
                    : variant == 6 ? ResourceStop::LABEL_NAMES
                    : variant == 7 ? ResourceStop::LEVEL_LISTS : ResourceStop::HEAD;
            q.row_demand = variant == 5 ? 0 : big;
            q.left = q.held = q.labels = q.label_bytes = big;
            q.lower_bound = variant == 5 || variant == 7;
            q.injected = variant == 1;
            q.arm = Arm::RIGHT;
            q.at_bp = big;
            q.limit = static_cast<double>(st.max_memory_bytes);
            q.used = 1;
            q.demand = q.limit + 1;
            if (q.resource == ResourceStop::WORK) {
                q.limit = static_cast<double>(big);
                q.used = q.demand = static_cast<double>(big) * 1.5;
            }
            r.resource_stop = q;
            r.account.decode_charged = variant != 3;
            r.account.row_diff_uncounted = variant == 3;
            r.account.largest_charge = big;
            r.account.soft_overshoot = big;
            for (ArmResult &arm : r.arms) {
                arm.status = ArmResult::TRUNCATED;
                arm.complete_to_bp = big;
                arm.cap_trigger = CapTrigger{ EndReason::RESOURCE_LIMIT, big, 0, 1, 1, true,
                                              static_cast<double>(big) };
            }
            const Json::Value j = seed_result_to_json(r, st, "summary", false);
            check_widths(j, 2, "", "variant " + std::to_string(variant));
            EXPECT_EQ(q.phase, j["resource_stop"]["phase"].asString());
        }
    }
}


/********* stage 3, review fixes: true causes, the window, the seed phase *********/

namespace {

// Annotate mode against a row that names many labels with long names: along P (label A) the
// node u = P[v + 1, v + 12) starts |n| records, each its own column with a name of |name_bytes|
// bytes, so reading u names n labels whose dictionary entries and delivery cost megabytes
// while u's row itself (n entries) costs kilobytes
struct NamesCase {
    std::shared_ptr<graph::AnnotatedDBG> anno;
    std::string P;
    size_t v = 30;
};
NamesCase names_case(size_t n = 300, size_t name_bytes = 4000) {
    NamesCase c;
    c.P = random_seq(80, 8100);
    std::vector<std::string> seqs { c.P }, labels { "A" };
    for (uint32_t i = 0; i < n; ++i) {
        seqs.push_back(c.P.substr(c.v + 1, 11) + random_seq(15, 8200 + i));
        labels.push_back("B" + std::to_string(i) + std::string(name_bytes, 'x'));
    }
    c.anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
            11, seqs, labels, DeBruijnGraph::BASIC);
    return c;
}

Json::Value names_request(const NamesCase &c, const std::string &sequence, const std::string &direction,
                          uint64_t memory_mb, size_t max_labels_per_node) {
    Json::Value r;
    Json::Value s;
    s["seed_id"] = "s";
    s["sequence"] = sequence;
    r["seeds"].append(s);
    r["strategy"] = parse_json(R"({"labels": {"mode": "annotate"}, "output": {"timing": false,
                                   "detail": "full"}, "bounds": {"max_extension_bp": 30}})");
    r["strategy"]["direction"] = direction;
    r["strategy"]["bounds"]["max_memory_mb"] = Json::UInt64(memory_mb);
    r["strategy"]["labels"]["max_labels_per_node"] = Json::UInt64(max_labels_per_node);
    return process_traverse_request(r, *c.anno, "");
}

std::set<std::string> actions_of(const Json::Value &q) {
    std::set<std::string> out;
    for (const Json::Value &a : q["actions"]) out.insert(a.asString());
    return out;
}

} // namespace

// F2: in annotate mode a level's key whose row fits but whose new dictionary labels do not is
// a stop by those labels — phase traversal, with the levers that name fewer or cheaper labels —
// not a row that "needs more than the walk had left"; with fewer labels per node the same
// budget walks on
TEST(GraphletStage3Review, LabelsThatDoNotFitAreNotARowStop) {
    const NamesCase c = names_case();
    const Json::Value out = names_request(c, c.P.substr(0, 20), "right", 2, 1000);
    const Json::Value &res = out["results"][0];
    ASSERT_TRUE(res.isMember("resource_stop")) << res.toStyledString().substr(0, 2000);
    const Json::Value &q = res["resource_stop"];
    EXPECT_EQ("traversal", q["phase"].asString());
    EXPECT_EQ("memory", q["resource"].asString());
    const std::string message = q["message"].asString();
    EXPECT_NE(std::string::npos, message.find("new dictionary label(s)")) << message;
    EXPECT_EQ(std::string::npos, message.find("needs more")) << message;
    const auto actions = actions_of(q);
    for (const char *a : { "raise_memory_budget", "use_graphlet", "lower_max_labels_per_node",
                           "label_constrained_query", "continue_from_leaves" }) {
        EXPECT_TRUE(actions.count(a)) << a;
    }
    EXPECT_FALSE(actions.count("more_selective_seed"));
    const uint64_t stopped_at = res["arms"]["right"]["complete_to_bp"].asUInt64();
    EXPECT_GT(stopped_at, 0u);
    // the lever works: with one label per node the same budget walks past u (where its 300
    // branches, one new label each, stop it again)
    const Json::Value fewer = names_request(c, c.P.substr(0, 20), "right", 2, 1)["results"][0];
    EXPECT_GT(fewer["arms"]["right"]["complete_to_bp"].asUInt64(), stopped_at);
}

// F7: an annotate root whose row fits but whose labels (with their delivery) do not fails the
// seed as a depth-0 state — phase traversal, the depth-0 levers, its need stated as a lower
// bound — not as the row's decoding; use_graphlet and lower_max_labels_per_node turn it into
// a walk
TEST(GraphletStage3Review, RootLabelsThatDoNotFitFailAsDepthZero) {
    const NamesCase c = names_case();
    // the right root is u
    const std::string seed = c.P.substr(c.v + 1 - 9, 20);
    const Json::Value res = names_request(c, seed, "right", 2, 1000)["results"][0];
    ASSERT_EQ("failed", res["outcome"]["walks"].asString()) << res.toStyledString().substr(0, 2000);
    const Json::Value &q = res["resource_stop"];
    EXPECT_EQ("traversal", q["phase"].asString());
    const std::string error = res["error"].asString();
    EXPECT_NE(std::string::npos, error.find("does not hold the seed's depth-0 state")) << error;
    EXPECT_NE(std::string::npos, error.find("need at least")) << error;
    const auto actions = actions_of(q);
    for (const char *a : { "raise_memory_budget", "use_graphlet", "lower_max_labels_per_node" }) {
        EXPECT_TRUE(actions.count(a)) << a;
    }
    // observed: at least the labels' bytes (300 names of 4 KB, priced twice and delivered)
    bool walk_domain = false;
    for (const Json::Value &l : res["limitations"]) {
        if (l["kind"].asString() != "walk_domain")
            continue;
        walk_domain = true;
        EXPECT_GE(l["observed"].asUInt64(), 3u);
        EXPECT_NE(std::string::npos, l["effect"].asString().find("at least"));
    }
    EXPECT_TRUE(walk_domain);
    const Json::Value fewer = names_request(c, seed, "right", 2, 1)["results"][0];
    EXPECT_NE("failed", fewer["outcome"]["walks"].asString());
}

// F7, byte by byte: a root's budget-aware read is reported as a decode failure only while the
// account had room for it; when the depth-0 state built before it (the other arm's root, its
// labels and reservation) already reached the budget, the seed fails as a depth-0 state whose
// need is stated as a lower bound
TEST(GraphletStage3Review, RootReadIsBlamedOnlyWithRoomLeft) {
    const NamesCase c = names_case(120, 200);
    LabelOracle oracle(*c.anno);
    const Seed seed = seed_of(c.P.substr(c.v + 1, 11));     // both roots are u
    size_t decode = 0, unread = 0, names = 0;
    for (uint64_t memory = 16384; memory < (uint64_t(8) << 20); memory = memory * 21 / 20 + 1) {
        Strategy st = strategy_of(R"({"direction": "both", "labels": {"mode": "annotate"},
                                      "bounds": {"max_extension_bp": 3}})");
        st.max_labels_per_node = 1000;
        st.max_memory_bytes = memory;
        st.delivery = delivery_costs("full", true);
        try {
            traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        } catch (const SeedBudgetError &e) {
            const ResourceStop &q = e.stop();
            if (std::string(q.phase) == "annotation_decode") {
                EXPECT_EQ(ResourceStop::ROOT, q.where) << memory;
                EXPECT_LT(q.used, static_cast<double>(memory)) << memory;
                EXPECT_GT(q.demand, static_cast<double>(memory)) << memory;
                decode++;
            } else if (q.lower_bound) {
                EXPECT_GT(q.demand, static_cast<double>(memory)) << memory;
                const std::string what = e.what();
                unread += what.find("was not read") != std::string::npos;
                names += what.find("further label(s)") != std::string::npos;
            }
        }
    }
    EXPECT_GT(unread + names, 0u);
    std::cerr << "root failures: " << decode << " decode, " << unread << " root not read, "
              << names << " labels that do not fit" << std::endl;
}

// F6: a memory stop whose used exceeds the budget states that excess as memory_bound_soft
// (it was observed as 0: the level's lists were never observed before the read refused), and a
// level whose own lists leave nothing for its read is stopped by them (phase traversal)
TEST(GraphletStage3Review, ExcessAtAStopIsObserved) {
    size_t over = 0, stops = 0;
    for (const BudgetCase &c : rowdiff_budget_cases()) {
        Strategy big = c.st;
        big.max_memory_bytes = uint64_t(1) << 30;
        const uint64_t peak = run_case(c, big).account.memory_peak;
        for (uint64_t memory = 8192; memory < 4 * peak; memory = memory * 51 / 50 + 1) {
            Strategy st = c.st;
            st.max_memory_bytes = memory;
            SeedResult r;
            try {
                r = run_case(c, st);
            } catch (const SeedBudgetError &) {
                continue;
            }
            if (!r.resource_stop || r.resource_stop->resource != ResourceStop::MEMORY)
                continue;
            stops++;
            const ResourceStop &q = *r.resource_stop;
            if (q.cause == ResourceStop::LEVEL_LISTS) {
                EXPECT_STREQ("traversal", q.phase) << c.name << " " << memory;
            }
            if (q.used > q.limit) {
                over++;
                EXPECT_GE(static_cast<double>(r.account.soft_overshoot), q.used - q.limit)
                    << c.name << " memory " << memory << " cause " << q.cause;
            }
        }
    }
    EXPECT_GT(stops, 50u);
    std::cerr << "memory stops: " << stops << ", " << over << " with used above the budget"
              << std::endl;
}

// F9: a seed phase that has run past the work budget by less than W when it is compared is
// let finish, and the walk stops at its first head with a valid result complete to 0 bp (it
// failed the seed); one that has run past it by W or more fails the seed. The overrun of the
// result stays within the stretch it states.
TEST(GraphletStage3Review, SeedPhasePastTheBudgetByLessThanAnIntervalFinishes) {
    // a long seed whose validation charges between one and two intervals
    std::string S = random_seq(6000, 8800);
    auto anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
            11, { S, S.substr(1000, 3000) }, { "A", "B" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    Strategy st = strategy_of(R"({"direction": "right", "bounds": {"max_extension_bp": 10}})");
    st.max_work_units = uint64_t(1) << 40;
    size_t len = 0;
    uint64_t seed_work = 0;
    for (size_t l = 100; l <= S.size(); l += 10) {
        const SeedResult r = traverse_seed(oracle, seed_of(S.substr(0, l), { "A" }), st,
                                           LabelChangeCost::forbid());
        if (r.account.work_seed > kWorkCheckInterval + 4096
                && r.account.work_seed + 4096 < 2 * kWorkCheckInterval) {
            len = l;
            seed_work = r.account.work_seed;
            break;
        }
    }
    ASSERT_GT(len, 0u) << "no seed length charges between one and two intervals";
    const Seed seed = seed_of(S.substr(0, len), { "A" });
    // compared at about W, over the budget by less than W: let finish
    st.max_work_units = kWorkCheckInterval - 2048;
    const SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
    ASSERT_TRUE(r.resource_stop);
    EXPECT_EQ(ResourceStop::WORK, r.resource_stop->resource);
    EXPECT_EQ(0u, r.arms[static_cast<size_t>(Arm::RIGHT)].complete_to_bp);
    EXPECT_GE(r.account.work_used, seed_work);
    EXPECT_GE(r.account.largest_charge, r.account.work_used - st.max_work_units);
    const Json::Value j = seed_result_to_json(r, st, "summary", false);
    EXPECT_EQ("partial", j["outcome"]["walks"].asString());
    // compared at about W, over the budget by W or more at the next comparison: failed
    st.max_work_units = 1;
    try {
        const SeedResult r1 = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
        // a seed phase that ended before its second comparison is let finish too
        EXPECT_LT(r1.account.work_seed, 2 * kWorkCheckInterval);
        EXPECT_EQ(0u, r1.arms[static_cast<size_t>(Arm::RIGHT)].complete_to_bp);
    } catch (const SeedBudgetError &e) {
        EXPECT_GE(e.stop().used - 1, static_cast<double>(kWorkCheckInterval));
    }
    // a longer seed phase is cut once it has run W past the budget
    const Seed longer = seed_of(S, { "A" });
    try {
        traverse_seed(oracle, longer, st, LabelChangeCost::forbid());
        FAIL() << "a seed phase far past the budget finished";
    } catch (const SeedBudgetError &e) {
        EXPECT_EQ(ResourceStop::WORK, e.stop().resource);
        EXPECT_GE(e.stop().used - 1, static_cast<double>(kWorkCheckInterval));
        EXPECT_LT(e.stop().used, 3.0 * kWorkCheckInterval);
    }
}

// F3: the derivation holds its window of up to 64 rows at once (to choose the cheapest). Its
// rows are read in runs and admitted one by one beside the rows before them, and a window that
// does not fit says so — the refused row's demand, what the seed phase had left and what the
// window's earlier rows held — instead of claiming that a row read alone needs more than the
// seed phase had, when each row alone fits easily
TEST(GraphletStage3Review, DerivationWindowStatesWhatItHolds) {
    // every k-mer of S[0, 40) occurs in 60 copies of it: rows of 60 coordinates, whose
    // row-diff paths share the anchors (empty diffs between consecutive k-mers)
    const std::string S = random_seq(60, 9100);
    std::string R;
    for (uint32_t i = 0; i < 60; ++i) {
        R += S.substr(0, 40) + random_seq(12, 9200 + i);
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
            11, { S, R }, { "A", "B" }, DeBruijnGraph::BASIC, true);
    LabelOracle oracle(*anno);
    ASSERT_TRUE(oracle.decode_charged());
    const Seed seed = seed_of(S.substr(0, 40));
    const Strategy base = strategy_of(R"({"direction": "right", "support": "trace",
        "branching": {"on_reconverge": "keep"}, "labels": {"seed_label_kind": "column"},
        "bounds": {"max_extension_bp": 10}})");
    size_t failures = 0, beside_rows = 0;
    uint64_t passes_from = 0;
    for (uint64_t memory = 4096; memory < (uint64_t(4) << 20); memory = memory * 21 / 20 + 1) {
        Strategy st = base;
        st.max_memory_bytes = memory;
        try {
            traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
            if (!passes_from)
                passes_from = memory;
        } catch (const SeedBudgetError &e) {
            const ResourceStop &q = e.stop();
            if (q.where != ResourceStop::DERIVATION)
                continue;
            failures++;
            EXPECT_STREQ("annotation_decode", q.phase) << memory;
            EXPECT_GT(q.demand, static_cast<double>(memory)) << memory;
            EXPECT_LT(q.index, 30u) << memory;
            const std::string what = e.what();
            EXPECT_NE(std::string::npos, what.find("row(s) read before it")) << what;
            EXPECT_EQ(std::string::npos, what.find("needs more than what the seed phase had left"))
                << what;
            if (!q.lower_bound) {
                EXPECT_GT(q.row_demand, q.left) << memory;
                // the window's rows read before the refused one were held beside it: the row
                // alone would have fitted the seed phase's budget
                if (q.row_demand < memory / 2 && q.held > memory / 2)
                    beside_rows++;
            }
        }
    }
    EXPECT_GT(failures, 0u);
    EXPECT_GT(beside_rows, 0u) << "no window was refused for the rows it held before a row";
    EXPECT_GT(passes_from, 0u);
}


/********* the stage-2 recheck (P1, P2) and the stage-3 review (F3), on stage 3 *********/

namespace {

// The reviewer's seed-phase index (review of the stage-2 recheck, P2): k = 3, one record
// "CGAT" + "GAT" x |n|, so that GAT, ATG and TGA each occur about |n| times (rows of |n|
// coordinates), labelled r1, with coordinates
std::shared_ptr<graph::AnnotatedDBG> repeat_index(size_t n, bool rowdiff) {
    std::string record = "CGAT";
    for (size_t i = 0; i < n; ++i) {
        record += "GAT";
    }
    if (rowdiff) {
        return test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
                3, { record }, { "r1" }, DeBruijnGraph::BASIC, true);
    }
    return test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            3, { record }, { "r1" }, DeBruijnGraph::BASIC, true);
}

size_t occurrences(size_t n, const std::string &kmer) {
    std::string record = "CGAT";
    for (size_t i = 0; i < n; ++i) {
        record += "GAT";
    }
    size_t count = 0;
    for (size_t i = 0; i + kmer.size() <= record.size(); ++i) {
        count += record.compare(i, kmer.size(), kmer) == 0;
    }
    return count;
}

const char kRepeatStrategy[] = R"({"direction": "right", "support": "trace",
    "branching": {"on_reconverge": "keep"},
    "labels": {"seed_label_kind": "column", "max_seed_labels": 10000},
    "bounds": {"max_extension_bp": 5}})";

} // namespace

// P2: every row a seed-phase read returned is charged before a comparison can fail the seed.
// The validation reads in calls grown from one key, so a failure has charged exactly the rows
// read (annotation.rows_requested); the derivation holds its window of up to 64 rows, all of
// them charged at once (the reviewer's 400,019 units, reported as 200,009)
TEST(GraphletStage2Recheck, SeedFetchChargesEveryReturnedRow) {
    const size_t n = 70'000;
    const uint64_t gat = 8 + 1 + occurrences(n, "GAT"), atg = 8 + 1 + occurrences(n, "ATG");
    for (bool rowdiff : { false, true }) {
        auto anno = repeat_index(n, rowdiff);
        Strategy st = strategy_of(kRepeatStrategy);
        st.max_work_units = 100;
        // explicit labels: the validation
        {
            LabelOracle oracle(*anno);
            try {
                traverse_seed(oracle, seed_of("GATG", { "r1" }), st, LabelChangeCost::forbid());
                ADD_FAILURE() << "rowdiff " << rowdiff << ": the validation was not failed";
            } catch (const SeedBudgetError &e) {
                EXPECT_EQ(ResourceStop::WORK, e.stop().resource);
                const uint64_t read = oracle.counters().rows_requested;
                ASSERT_GE(read, 1u);
                const uint64_t rows = read == 1 ? gat : gat + atg;
                if (rowdiff) {
                    // each row with its row-diff dependency rows
                    EXPECT_GE(e.stop().used, static_cast<double>(rows)) << read;
                } else {
                    EXPECT_EQ(static_cast<double>(rows), e.stop().used) << read;
                }
                EXPECT_EQ(e.stop().used, static_cast<double>(e.account().work_seed));
            }
        }
        // derived labels: the derivation's window (both k-mers) is one charge
        {
            LabelOracle oracle(*anno);
            try {
                traverse_seed(oracle, seed_of("GATG"), st, LabelChangeCost::forbid());
                ADD_FAILURE() << "rowdiff " << rowdiff << ": the derivation was not failed";
            } catch (const SeedBudgetError &e) {
                EXPECT_EQ(ResourceStop::WORK, e.stop().resource);
                if (rowdiff) {
                    EXPECT_GE(e.stop().used, static_cast<double>(gat + atg));
                } else {
                    EXPECT_EQ(static_cast<double>(gat + atg), e.stop().used);
                }
                // what the seed phase charged between two comparisons: the window
                EXPECT_EQ(e.stop().used, static_cast<double>(e.account().largest_charge));
            }
        }
    }
}

// P2: a seed the work budget failed states how far its seed phase ran past the budget, with
// the most it charged between two comparisons, as a walk's work stop does; and the wording is
// the threshold's (review of stage 3, answer 3), never an exact overshoot
TEST(GraphletStage2Recheck, FailedWorkStopStatesLargestCharge) {
    auto anno = repeat_index(100'000, false);
    Json::Value r;
    Json::Value seed;
    seed["sequence"] = "GATG";
    r["seeds"].append(seed);
    r["strategy"] = parse_json(kRepeatStrategy);
    r["strategy"]["bounds"]["max_work_units"] = 100;
    r["strategy"]["output"]["timing"] = false;
    const Json::Value out = process_traverse_request(r, *anno, "");
    const Json::Value &res = out["results"][0];
    ASSERT_EQ("failed", res["outcome"]["walks"].asString()) << res;
    const Json::Value &q = res["resource_stop"];
    ASSERT_EQ("work", q["resource"].asString());
    const std::string message = q["message"].asString();
    EXPECT_NE(std::string::npos, message.find("fails at a comparison finding it at least "
                                              + std::to_string(kWorkCheckInterval)
                                              + " units over budget")) << message;
    EXPECT_NE(std::string::npos, message.find("the most this seed charged between two "
                                              "comparisons: " + q["used"].asString() + " units"))
        << message;
    EXPECT_EQ(std::string::npos, message.find("past the budget")) << message;
    check_widths(res, 2, "", "failed work seed");
}

// P2: the trace validation's coordinate copies — every label's live set and both arms'
// boundary coordinates, beside the hits — are observed before the depth-0 admission can fail
// the seed: the roots case (seed GAT, both arms, 1 MiB) reported 3 MiB where the copies held
// about 9.6 MB
TEST(GraphletStage2Recheck, ValidationCoordinatesAreObserved) {
    const size_t n = 100'000;
    auto anno = repeat_index(n, false);
    LabelOracle oracle(*anno);
    Strategy st = strategy_of(kRepeatStrategy);
    st.direction = Strategy::BOTH;
    st.max_memory_bytes = uint64_t(1) << 20;
    try {
        traverse_seed(oracle, seed_of("GAT", { "r1" }), st, LabelChangeCost::forbid());
        FAIL() << "the depth-0 state of two roots of " << n << " coordinates fit 1 MiB";
    } catch (const SeedBudgetError &e) {
        // the hits, the live set and both boundary sets: four copies of the coordinates
        const uint64_t copies = 4 * occurrences(n, "GAT") * sizeof(uint64_t);
        EXPECT_GE(e.account().soft_overshoot, copies - st.max_memory_bytes);
    }
}

// The stage-2 recheck's design answer: an extra label is accepted when a CHAIN of switches
// from a seed label enters it within the loss budget (the walk enforces the cumulative loss);
// A -> B = 1, B -> C = 1 under a budget of 2 rejected C. An extra label no chain reaches is
// still refused, by name.
TEST(GraphletStage2Recheck, ExtraLabelsReachableByAChain) {
    // A carries S T1, B the last 20 bp of T1 and T2, C the last 20 bp of T2 and T3
    std::string S, T1, T2, T3;
    for (uint32_t t = 700; ; t += 4) {
        S = random_seq(30, t); T1 = random_seq(30, t + 1); T2 = random_seq(30, t + 2);
        T3 = random_seq(30, t + 3);
        if (T1[0] != T2[0])
            break;
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, { S + T1, T1.substr(10) + T2, T2.substr(10) + T3 }, { "A", "B", "C" },
            DeBruijnGraph::BASIC);
    auto request = [&](double budget, const std::string &entries, const std::string &fallback) {
        Json::Value r;
        Json::Value seed;
        seed["sequence"] = S;
        seed["labels"].append("A");
        r["seeds"].append(seed);
        r["strategy"] = parse_json(R"({"direction": "right", "labels": {"extra": ["B", "C"],
            "change_cost": {"model": "table", "entries": )" + entries + R"(, "default": )"
            + fallback + R"(}}, "bounds": {"max_extension_bp": 80}, "output": {"timing": false}})");
        r["strategy"]["labels"]["loss_budget"] = budget;
        return r;
    };
    const std::string chain = R"([["A", "B", 1], ["B", "C", 1]])";
    // C is two switches away: accepted at 2, and the walk enters it at its cumulative loss
    {
        const Json::Value out = process_traverse_request(request(2, chain, "\"forbid\""), *anno, "");
        LabelOracle oracle(*anno);
        Strategy st = parse_traverse_request(request(2, chain, "\"forbid\"")).strategy;
        st.delivery = delivery_costs("full", true);
        const SeedResult r = traverse_seed(oracle, seed_of(S, { "A" }), st,
                                           LabelChangeCost::table({ { { 0, 1 }, 1.0 }, { { 1, 2 }, 1.0 } },
                                                                  kInfiniteLoss));
        bool entered_c = false;
        for (const LabelRun &run : r.arms[static_cast<size_t>(Arm::RIGHT)].runs) {
            if (run.label == 2) {
                entered_c = true;
                EXPECT_TRUE(run.entered_by_switch);
                EXPECT_LE(run.loss, 2.0);
                EXPECT_EQ(1u, run.from_label);
            }
        }
        EXPECT_TRUE(entered_c) << out;
        EXPECT_EQ(1u, out["results"].size());
    }
    // within 1.5 no chain reaches C: refused by name
    try {
        process_traverse_request(request(1.5, chain, "\"forbid\""), *anno, "");
        ADD_FAILURE() << "an unreachable extra label was accepted";
    } catch (const std::exception &e) {
        const std::string what = e.what();
        EXPECT_NE(std::string::npos, what.find("'C'")) << what;
        EXPECT_EQ(std::string::npos, what.find("'B'")) << what;
        EXPECT_NE(std::string::npos, what.find("no chain of switches")) << what;
    }
    // a finite default reaches C through B although the explicit A -> C entry is above the
    // budget (an explicit entry overrides the default of its own pair only)
    process_traverse_request(request(2, R"([["A", "C", 10]])", "1"), *anno, "");
    EXPECT_THROW(process_traverse_request(request(1.5, R"([["A", "C", 10]])", "1"), *anno, ""),
                 std::exception);
}

// switch_reach against a plain fixpoint over every pair, on random tables with and without a
// finite default, and a pool of many labels under a dense default (lazy, not quadratic)
TEST(GraphletStage2Recheck, SwitchReachIsTheCheapestChain) {
    std::mt19937 gen(4242);
    for (size_t round = 0; round < 300; ++round) {
        const size_t n = 2 + gen() % 12, sources = 1 + gen() % std::min<size_t>(3, n - 1);
        std::map<std::pair<LabelId, LabelId>, double> entries;
        for (size_t e = gen() % (n * 2); e > 0; --e) {
            const LabelId a = gen() % n, b = gen() % n;
            entries[{ a, b }] = gen() % 3 ? (gen() % 7) * 0.5 : kInfiniteLoss;
        }
        const double fallback = gen() % 2 ? kInfiniteLoss : (gen() % 6) * 0.5;
        const double budget = (gen() % 8) * 0.5;
        const LabelChangeCost cost = LabelChangeCost::table(entries, fallback);
        std::vector<double> naive(n, kInfiniteLoss);
        for (size_t s = 0; s < sources; ++s) naive[s] = 0;
        for (bool changed = true; changed; ) {
            changed = false;
            for (LabelId u = 0; u < n; ++u) {
                for (LabelId v = 0; v < n; ++v) {
                    const double loss = naive[u] + cost.cost(u, v);
                    if (naive[u] != kInfiniteLoss && loss <= budget && loss < naive[v]) {
                        naive[v] = loss;
                        changed = true;
                    }
                }
            }
        }
        EXPECT_EQ(naive, switch_reach(cost, n, sources, budget)) << "round " << round;
    }
    // constant and forbid
    EXPECT_EQ(std::vector<double>({ 0, 1, 1 }), switch_reach(LabelChangeCost::constant(1), 3, 1, 1));
    EXPECT_EQ(std::vector<double>({ 0, kInfiniteLoss }),
              switch_reach(LabelChangeCost::constant(2), 2, 1, 1));
    EXPECT_EQ(std::vector<double>({ 0, kInfiniteLoss }), switch_reach(LabelChangeCost::forbid(), 2, 1, 5));
    // 20,000 labels under a finite default with a few overrides: every one reached
    const size_t many = 20'000;
    std::map<std::pair<LabelId, LabelId>, double> few { { { 0, 5 }, 9.0 }, { { 0, 6 }, 9.0 } };
    const auto start = std::chrono::steady_clock::now();
    const std::vector<double> reach = switch_reach(LabelChangeCost::table(few, 1.0), many, 1, 2);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    EXPECT_LT(seconds, 2.0);
    EXPECT_EQ(2.0, reach[5]);
    EXPECT_EQ(1.0, reach[7]);
    EXPECT_EQ(many, static_cast<size_t>(std::count_if(reach.begin(), reach.end(),
                                                      [](double x) { return x <= 2; })));
}

// P1: the float width the delivery model prices is at least the width of every float a result
// under the request can hold: costs, sums of costs up to the loss budget plus one cost, and
// the time budgets with an elapsed time above them
TEST(GraphletStage2Recheck, FloatWidthBoundsEveryFloat) {
    std::mt19937_64 gen(77);
    for (double c : { 1e-300, 3e-17, 0.1, 0.30000000000000004, 1.0, 7.25, 1e20, 1e300 }) {
        for (double budget : { 0.0, 1.0, 1e10, 1e300 }) {
            Strategy st;
            st.loss_budget = budget;
            const uint64_t width = mgt_float_width(st, LabelChangeCost::constant(c));
            EXPECT_GE(width, kMgtFloatWidth);
            EXPECT_GE(width, encode_float(c).size()) << c;
            double loss = 0;
            for (size_t i = 0; i < 1000 && loss + c <= budget + c; ++i) {
                loss += c;                                    // a lineage's losses, and the
                EXPECT_GE(width, encode_float(loss).size());  // budget a loss-budget end needs
            }
            // values between the smallest cost and the budget plus a cost, 17 digits each
            for (size_t i = 0; i < 200; ++i) {
                const double lo = std::log10(c), hi = std::log10(budget + c);
                const double x = std::pow(10.0, lo + (hi - lo) * (gen() % 100000) / 100000.0)
                               * (1 + 1e-16 * (gen() % 7));
                if (x >= c && x <= budget + c) {
                    EXPECT_GE(width, encode_float(x).size()) << x;
                }
            }
        }
    }
    for (double t : { 1e-300, 1e-5, 0.5, 30000.0, 1e17, 1e300 }) {
        Strategy st;
        st.time_budget_ms = t;
        const uint64_t width = mgt_float_width(st, LabelChangeCost::forbid());
        EXPECT_GE(width, encode_float(t).size()) << t;
        // an elapsed time at least the budget, measured in ns ticks
        for (double elapsed : { t * 1.0000000000000002, t * 1.2345678901234567 }) {
            if (elapsed < 1e17) {
                EXPECT_GE(width, encode_float(elapsed).size()) << elapsed;
            }
        }
    }
    // a usual request: 24, the model's widths unchanged
    Strategy usual;
    usual.loss_budget = 2;
    EXPECT_EQ(kMgtFloatWidth, mgt_float_width(usual, LabelChangeCost::constant(1)));
    EXPECT_EQ(kMgtFloatWidth, mgt_float_width(usual, LabelChangeCost::forbid(), 1e9, 3e5));
}

// Review of stage 3, F3: a level with no key to read (the radius, where heads only end) is
// not admitted as a read: a completed radius-0 walk whose depth-0 state fills the budget
// exactly reported a memory stop and a truncated arm. Every budget from the first that holds
// the depth-0 state completes (radius 0), or stops only below the radius (radius 1).
TEST(GraphletStage3Review, RadiusOnlyLevelIsNotStopped) {
    auto anno = test::build_anno_graph<DBGSuccinct, annot::RowDiffColumnAnnotator>(
            3, { "AAACAAAGAAAT" }, { "A" }, DeBruijnGraph::BASIC);
    LabelOracle oracle(*anno);
    ASSERT_TRUE(oracle.decode_charged());
    const Seed seed = seed_of("AAA", { "A" });
    for (uint64_t radius : { 0, 1 }) {
        Strategy st;
        st.direction = Strategy::RIGHT;
        st.max_extension_bp = radius;
        st.delivery = delivery_costs("graphlet", true);
        auto admitted = [&](uint64_t memory) {
            st.max_memory_bytes = memory;
            try {
                traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
                return true;
            } catch (const SeedBudgetError &) {
                return false;
            }
        };
        uint64_t lo = 1, hi = uint64_t(1) << 22;
        while (lo < hi) {
            const uint64_t mid = lo + (hi - lo) / 2;
            if (admitted(mid)) hi = mid; else lo = mid + 1;
        }
        for (uint64_t memory = lo; memory < lo + 256; ++memory) {
            st.max_memory_bytes = memory;
            const SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
            const ArmResult &arm = r.arms[static_cast<size_t>(Arm::RIGHT)];
            if (radius == 0) {
                EXPECT_FALSE(r.resource_stop) << memory;
                EXPECT_EQ(ArmResult::COMPLETE, arm.status) << memory;
            } else if (r.resource_stop) {
                EXPECT_LT(r.resource_stop->at_bp, radius) << memory;
            } else {
                EXPECT_EQ(ArmResult::COMPLETE, arm.status) << memory;
            }
        }
    }
}

// Review of stage 3, answer 2: every "at least" — a depth-0 failure's message and effect, a
// refused read's effect at the seed or at a level, a level's lists — says that raising the
// budget to that value may still fail; and every statement, the failed seeds' with their
// caveats and largest charges included, keeps to its width (effect 640, message 1,024)
TEST(GraphletStage2Recheck, LowerBoundsAreNoPromiseAndFitTheirWidths) {
    size_t lower = 0, failed = 0, checked = 0;
    auto inspect = [&](const Json::Value &res, const std::string &what) {
        check_widths(res, 2, "", what);
        checked++;
        failed += res["outcome"]["walks"].asString() == "failed";
        std::vector<std::string> texts { res["error"].asString(),
                                         res["resource_stop"]["message"].asString() };
        for (const Json::Value &l : res["limitations"]) texts.push_back(l["effect"].asString());
        for (const std::string &side : res["arms"].getMemberNames()) {
            for (const Json::Value &l : res["arms"][side]["limitations"]) {
                texts.push_back(l["effect"].asString());
            }
        }
        for (const std::string &t : texts) {
            if (t.find("at least") == std::string::npos || t.find("units over budget") != std::string::npos)
                continue;
            lower++;
            EXPECT_NE(std::string::npos, t.find("may still fail")) << what << ": " << t;
        }
    };
    // annotate roots whose labels do not fit (depth-0 lower bounds), over budgets
    const NamesCase c = names_case();
    for (uint64_t mb : { 1, 2, 3, 5 }) {
        // const char *: a std::string reference would bind to a temporary per element
        // (GCC's -Wrange-loop-construct)
        for (const char *direction : { "right", "both" }) {
            const std::string seed = c.P.substr(c.v + 1 - 9, 20);
            inspect(names_request(c, seed, direction, mb, 1000)["results"][0],
                    std::string("names ") + direction + " " + std::to_string(mb));
        }
    }
    // level reads and level lists on row-diff budget cases (arm-level lower bounds)
    for (const BudgetCase &bc : rowdiff_budget_cases()) {
        for (uint64_t memory = 8192; memory < (uint64_t(1) << 22); memory = memory * 3 / 2) {
            Strategy st = bc.st;
            st.max_memory_bytes = memory;
            st.delivery = delivery_costs("full", st.sequences);
            try {
                const SeedResult r = run_case(bc, st);
                inspect(seed_result_to_json(r, st, "full", false), bc.name);
            } catch (const SeedBudgetError &) {
                // failed seeds are inspected through the request below
            }
        }
    }
    // failed seeds through the request: work in validation and derivation, memory in the
    // trace roots
    auto anno = repeat_index(70'000, false);
    for (bool named : { false, true }) {
        for (const char *budget : { "max_work_units", "max_memory_mb" }) {
            Json::Value r;
            Json::Value seed;
            seed["sequence"] = named ? "GAT" : "GATG";
            if (named)
                seed["labels"].append("r1");
            r["seeds"].append(seed);
            r["strategy"] = parse_json(kRepeatStrategy);
            r["strategy"]["direction"] = "both";
            r["strategy"]["bounds"][budget] = 1;
            r["strategy"]["output"]["timing"] = false;
            inspect(process_traverse_request(r, *anno, "")["results"][0],
                    std::string(named ? "validation " : "derivation ") + budget);
        }
    }
    EXPECT_GT(lower, 5u);
    EXPECT_GE(failed, 4u);
    EXPECT_GT(checked, 20u);
}


/********* stage 4, backend half: an attempt stopped from outside the walk *********/

namespace {

// A control whose poll answers |stop| from its |at|-th call on (0-based), counting every call:
// a cancel, the attempt's bound or a gone client injected at any checkpoint the walk polls
struct InjectedStop {
    uint64_t at;
    ExternalStop stop;
    uint64_t calls = 0;
    AttemptMeter meter;
    AttemptControl control;
    explicit InjectedStop(uint64_t at = std::numeric_limits<uint64_t>::max(),
                          ExternalStop stop = ExternalStop::CANCELLED)
          : at(at), stop(stop) {
        control.poll = [this]() {
            return calls++ >= this->at ? this->stop : ExternalStop::NONE;
        };
        control.elapsed_ms = [this]() { return 1000.5 + static_cast<double>(calls); };
        control.bound_ms = 70'000;
        control.meter = &meter;
    }
    InjectedStop(const InjectedStop&) = delete;
    InjectedStop& operator=(const InjectedStop&) = delete;
};

SeedResult run_stopped(const BudgetCase &c, const Strategy &st, InjectedStop *inject) {
    LabelOracle oracle(*c.anno);
    return traverse_seed(oracle, c.seed, st, c.cost, "", nullptr, &inject->control);
}

const char* stop_token(ExternalStop s) {
    return s == ExternalStop::CANCELLED ? "cancelled" : "attempt_deadline";
}

ResourceStop::Resource stop_resource(ExternalStop s) {
    return s == ExternalStop::CANCELLED ? ResourceStop::CANCELLED : ResourceStop::ATTEMPT_DEADLINE;
}

// A server-side attempt (traverse_attempts.hpp) whose |at|-th poll (1-based, the polls
// between seeds included) finds it cancelled, as POST /traverse/cancel would: the client check
// runs at every poll and asks for the cancel itself
struct CancelAtPoll {
    AttemptSettings settings;
    std::shared_ptr<Attempt> attempt;
    uint64_t calls = 0;
    uint64_t at;
    explicit CancelAtPoll(uint64_t at, const std::string &id = "t-1") : at(at) {
        settings.client_check_ms = 0;
        settings.poll_stride = 1;
        settings.allowance_ms = 1000;
        attempt = std::make_shared<Attempt>(0, std::chrono::system_clock::time_point(), settings,
                                            "0123456789abcdef", true, [this]() {
            if (++calls == this->at)
                attempt->cancel();
            return false;
        });
        AttemptIds ids;
        ids.attempt_id = id;
        attempt->set_ids(ids);
    }
};

Json::Value request_of(const std::vector<Seed> &seeds, const std::string &strategy,
                       const std::string &detail = "graphlet") {
    Json::Value r;
    for (const Seed &s : seeds) {
        Json::Value j;
        j["sequence"] = s.sequence;
        for (const std::string &l : s.labels) {
            j["labels"].append(l);
        }
        r["seeds"].append(j);
    }
    r["strategy"] = parse_json(strategy);
    r["strategy"]["output"]["detail"] = detail;
    r["strategy"]["output"]["timing"] = false;
    return r;
}

std::string compact_json(const Json::Value &v) {
    Json::StreamWriterBuilder b;
    b["indentation"] = "";
    return Json::writeString(b, v);
}

} // namespace

// A cancel, or the attempt's bound, at every checkpoint the walk polls leaves a consistent
// prefix (as an allocation denial does, §14): every walk up to each arm's complete_to_bp is
// the unstopped walk's, the stop is stated (resource_stop with scope attempt, Q, a walk_domain
// naming attempt_id), the result serialises both ways, and a later stop never walks less. A
// stop in the seed phase fails the seed (no result exists to deliver).
TEST(GraphletAttempt, StopAtEveryPollLeavesAConsistentPrefix) {
    size_t stops = 0, seed_phase = 0, mid_walk = 0;
    for (const BudgetCase &c : budget_cases()) {
        InjectedStop count;
        const SeedResult base = run_stopped(c, c.st, &count);
        EXPECT_FALSE(base.resource_stop) << c.name;
        const uint64_t polls = count.calls;
        ASSERT_GT(polls, 0u) << c.name;
        const size_t most = c.name.rfind("mode ", 0) == 0 ? 6 : 16;
        std::vector<uint64_t> ats;
        for (size_t i = 0; i < std::min<uint64_t>(polls, most); ++i) {
            ats.push_back(polls <= most ? i : i * (polls - 1) / (most - 1));
        }
        LabelOracle oracle(*c.anno);
        for (ExternalStop which : { ExternalStop::CANCELLED, ExternalStop::ATTEMPT_DEADLINE }) {
            std::array<uint64_t, 2> depth { 0, 0 };
            for (size_t i = 0; i < ats.size(); ++i) {
                // the bound's stop takes the same path as a cancel's: every other poll
                if (which == ExternalStop::ATTEMPT_DEADLINE && i % 2)
                    continue;
                const uint64_t at = ats[i];
                const std::string what = c.name + " " + stop_token(which) + " at poll "
                                       + std::to_string(at) + " of " + std::to_string(polls);
                InjectedStop inject(at, which);
                SeedResult r;
                try {
                    r = run_stopped(c, c.st, &inject);
                } catch (const SeedBudgetError &e) {
                    seed_phase++;
                    EXPECT_EQ(stop_resource(which), e.stop().resource) << what;
                    EXPECT_EQ(70'000, e.stop().limit) << what;
                    EXPECT_TRUE(inject.meter.walked) << what;
                    continue;
                }
                stops++;
                ASSERT_TRUE(r.resource_stop) << what;
                const ResourceStop &q = *r.resource_stop;
                EXPECT_EQ(stop_resource(which), q.resource) << what;
                EXPECT_EQ(70'000, q.limit) << what;
                EXPECT_EQ(std::ceil(q.used), q.used) << what << ": whole milliseconds";
                EXPECT_GE(q.used, 1001) << what;
                EXPECT_TRUE(inject.meter.walked) << what;
                EXPECT_EQ(r.account.work_used, inject.meter.work_units) << what;
                bool deeper = false;
                for (size_t a = 0; a < 2; ++a) {
                    const ArmResult &ra = r.arms[a];
                    if (!ra.requested)
                        continue;
                    const std::string arm = what + " arm " + to_string(ra.arm);
                    check_bins(ra, arm);
                    EXPECT_EQ(walks_upto(base.arms[a], ra.complete_to_bp),
                              walks_upto(ra, ra.complete_to_bp)) << arm;
                    EXPECT_GE(ra.complete_to_bp, depth[a]) << arm << ": a later stop walks no less";
                    depth[a] = ra.complete_to_bp;
                    deeper |= ra.complete_to_bp > 0;
                    for (const PathResult &p : ra.paths) {
                        if (p.length_bp >= ra.complete_to_bp && p.path_reason) {
                            EXPECT_TRUE(*p.path_reason == EndReason::RESOURCE_LIMIT
                                        || !is_resource_stop(*p.path_reason)
                                        || *p.path_reason == EndReason::BEAM_PRUNED)
                                << arm << " path " << p.id << " " << to_string(*p.path_reason);
                        }
                    }
                }
                mid_walk += deeper;
                const std::string text = check_serialised(r, c.seed, c.st, oracle, what);
                EXPECT_NE(std::string::npos,
                          text.find(std::string("\nQ attempt ") + stop_token(which) + " traversal "))
                    << what;
                const Json::Value j = seed_result_to_json(r, c.st, "summary", false);
                check_widths(j, 2, "", what);
                EXPECT_EQ("partial", j["outcome"]["walks"].asString()) << what;
                EXPECT_EQ("attempt", j["resource_stop"]["scope"].asString()) << what;
                EXPECT_EQ(stop_token(which), j["resource_stop"]["resource"].asString()) << what;
                EXPECT_EQ("continue_from_leaves", j["resource_stop"]["actions"][0].asString()) << what;
                EXPECT_EQ(1u, j["resource_stop"]["actions"].size()) << what;
                bool stated = false;
                for (const char *side : { "left", "right" }) {
                    for (const Json::Value &l : j["arms"][side]["limitations"]) {
                        if (l["knob"].asString() != "attempt_id")
                            continue;
                        stated = true;
                        EXPECT_EQ("walk_domain", l["kind"].asString()) << what;
                        EXPECT_EQ(70'000u, l["limit"].asUInt64()) << what;
                        EXPECT_EQ(std::string::npos, l["effect"].asString().find("raise the knob"))
                            << what << ": no budget ran out";
                    }
                }
                EXPECT_TRUE(stated) << what;
            }
        }
    }
    EXPECT_GT(stops, 200u);
    EXPECT_GT(seed_phase, 10u);
    EXPECT_GT(mid_walk, 100u);
}

// No stop: the walk is the one without a control, byte for byte, and the meter is its account
TEST(GraphletAttempt, NoStopIsByteIdentical) {
    for (const BudgetCase &c : budget_cases()) {
        LabelOracle oracle(*c.anno);
        for (bool budgeted : { false, true }) {
            Strategy st = c.st;
            if (budgeted) {
                st.max_memory_bytes = uint64_t(64) << 20;
                st.max_work_units = 1'000'000'000;
            }
            const SeedResult plain = run_case(c, st);
            InjectedStop never;
            const SeedResult controlled = run_stopped(c, st, &never);
            const std::string what = c.name + (budgeted ? " budgeted" : "");
            EXPECT_GT(never.calls, 0u) << what;
            EXPECT_EQ(compact_json(seed_result_to_json(plain, st, "full", false)),
                      compact_json(seed_result_to_json(controlled, st, "full", false))) << what;
            const Json::Value summary = seed_result_to_json(plain, st, "graphlet", false);
            EXPECT_EQ(graphlet_text(plain, c.seed, st, context_of(oracle), summary),
                      graphlet_text(controlled, c.seed, st, context_of(oracle), summary)) << what;
            EXPECT_TRUE(never.meter.walked) << what;
            EXPECT_EQ(plain.account.work_used, never.meter.work_units) << what;
            EXPECT_EQ(plain.account.work_seed, never.meter.work_seed) << what;
            EXPECT_GE(never.meter.memory_peak, plain.account.memory_peak) << what;
            EXPECT_EQ(plain.account.memory_final, never.meter.memory_final) << what;
            EXPECT_EQ(plain.account.soft_overshoot, never.meter.soft_excess) << what;
        }
    }
}

// A gone client abandons the walk wherever it is: nothing is finalised, and what it consumed
// until then still reaches the caller
TEST(GraphletAttempt, GoneClientAbandonsTheWalk) {
    for (const BudgetCase &c : budget_cases()) {
        InjectedStop count;
        const SeedResult base = run_stopped(c, c.st, &count);
        for (uint64_t at : { uint64_t(0), count.calls / 2, count.calls - 1 }) {
            InjectedStop gone(at, ExternalStop::CLIENT_GONE);
            EXPECT_THROW(run_stopped(c, c.st, &gone), AttemptAborted) << c.name << " at " << at;
            EXPECT_TRUE(gone.meter.walked) << c.name;
            EXPECT_LE(gone.meter.work_units, base.account.work_used) << c.name;
        }
    }
}

// A stop in the seed phase (a derivation's window, a validation's fetch call) fails the seed:
// no result exists yet, not even one complete to 0 bp
TEST(GraphletAttempt, SeedPhaseStopFailsTheSeed) {
    BudgetCase c = budget_cases()[0];
    for (bool derived : { true, false }) {
        c.seed = derived ? seed_of("AAA") : seed_of("AAA", { "C", "D", "E" });
        InjectedStop inject(0, ExternalStop::CANCELLED);
        try {
            run_stopped(c, c.st, &inject);
            ADD_FAILURE() << "derived " << derived << ": not stopped in its seed phase";
        } catch (const SeedBudgetError &e) {
            EXPECT_EQ(ResourceStop::CANCELLED, e.stop().resource);
            EXPECT_NE(std::string::npos,
                      std::string(e.what()).find(derived ? "derive" : "validated")) << e.what();
            EXPECT_EQ(e.stop().used, e.stop().demand);
        }
        // through the request: failed per seed with the attempt's stop, and its lever is a new
        // attempt (the polls between seeds come first: the walk's first poll is the second)
        CancelAtPoll cancel(2);
        const Json::Value out = process_traverse_request(
                request_of({ c.seed }, R"({"bounds": {"max_extension_bp": 7}})"), *c.anno, "", {},
                nullptr, cancel.attempt.get());
        const Json::Value &r = out["results"][0];
        const std::string what = std::string("derived ") + (derived ? "true" : "false") + ": "
                               + compact_json(r);
        check_widths(r, 2, "", what);
        EXPECT_EQ("failed", r["outcome"]["walks"].asString()) << what;
        EXPECT_FALSE(r.isMember("graphlet")) << what;
        EXPECT_EQ("attempt", r["resource_stop"]["scope"].asString()) << what;
        EXPECT_EQ("cancelled", r["resource_stop"]["resource"].asString()) << what;
        EXPECT_EQ("traversal", r["resource_stop"]["phase"].asString()) << what;
        EXPECT_EQ("retry_attempt", r["resource_stop"]["actions"][0].asString()) << what;
        EXPECT_EQ("walk_domain", r["limitations"][0]["kind"].asString()) << what;
        EXPECT_EQ("attempt_id", r["limitations"][0]["knob"].asString()) << what;
        EXPECT_EQ(cancel.attempt->bound_ms(), r["limitations"][0]["limit"].asDouble()) << what;
        const Json::Value &u = out["usage"];
        EXPECT_EQ("cancelled", u["reason"].asString()) << what;
        EXPECT_EQ("failed", u["per_seed"][0]["outcome"].asString()) << what;
        EXPECT_EQ("cancelled", u["per_seed"][0]["stopped_by"].asString()) << what;
        EXPECT_EQ(1u, u["seeds"]["finished"].asUInt64()) << what;
    }
}

// A stop between seeds leaves the rest unstarted: each failed with the stop that left it so
// (resource_stop phase not_started, a walk_domain naming attempt_id), the seeds walked before
// it delivered whole, and the usage stating which seeds ran
TEST(GraphletAttempt, UnstartedSeedsAreFailedWithTheStop) {
    const BudgetCase c = budget_cases()[0];
    const std::string strategy = R"({"bounds": {"max_extension_bp": 7}})";
    for (bool budgeted : { false, true }) {
        Json::Value one = request_of({ c.seed }, strategy);
        if (budgeted)
            one["strategy"]["bounds"]["max_memory_mb"] = 64;
        CancelAtPoll single(std::numeric_limits<uint64_t>::max());
        process_traverse_request(one, *c.anno, "", {}, nullptr, single.attempt.get());
        // the polls of one seed, the poll before it included (a budget polls between its fetch
        // calls too): the next poll is the one before the second seed
        const uint64_t polls = single.calls;
        ASSERT_GT(polls, 1u);
        Json::Value request = one;
        request["seeds"].append(one["seeds"][0]);
        request["seeds"].append(one["seeds"][0]);
        CancelAtPoll cancel(polls + 1);
        const Json::Value out = process_traverse_request(request, *c.anno, "", {}, nullptr,
                                                         cancel.attempt.get());
        ASSERT_EQ(3u, out["results"].size());
        const Json::Value &first = out["results"][0];
        EXPECT_FALSE(first.isMember("resource_stop"));
        EXPECT_TRUE(first.isMember("graphlet"));
        for (Json::ArrayIndex i = 1; i < 3; ++i) {
            const Json::Value &r = out["results"][i];
            const std::string what = "seed " + std::to_string(i) + ": " + compact_json(r);
            check_widths(r, 2, "", what);
            EXPECT_EQ(0u, r["error"].asString().rfind("not started: the attempt was cancelled", 0))
                << what;
            EXPECT_EQ("failed", r["outcome"]["walks"].asString()) << what;
            EXPECT_EQ("not_started", r["resource_stop"]["phase"].asString()) << what;
            EXPECT_EQ("cancelled", r["resource_stop"]["resource"].asString()) << what;
            EXPECT_EQ("attempt", r["resource_stop"]["scope"].asString()) << what;
            EXPECT_EQ("retry_attempt", r["resource_stop"]["actions"][0].asString()) << what;
            EXPECT_EQ("attempt_id", r["limitations"][0]["knob"].asString()) << what;
            // every response under a memory budget states it (§7.0)
            EXPECT_EQ(budgeted ? 2u : 1u, r["limitations"].size()) << what;
            if (budgeted) {
                EXPECT_EQ("memory_bound_soft", r["limitations"][1]["kind"].asString()) << what;
            }
        }
        const Json::Value &u = out["usage"];
        EXPECT_EQ("cancelled", u["reason"].asString());
        EXPECT_EQ("t-1", u["attempt_id"].asString());
        EXPECT_EQ(3u, u["seeds"]["requested"].asUInt64());
        EXPECT_EQ(1u, u["seeds"]["started"].asUInt64());
        EXPECT_EQ(1u, u["seeds"]["finished"].asUInt64());
        EXPECT_EQ(first["outcome"]["walks"].asString(), u["per_seed"][0]["outcome"].asString());
        EXPECT_EQ("not_started", u["per_seed"][1]["outcome"].asString());
        EXPECT_EQ(0u, u["per_seed"][2]["work_units"].asUInt64());
        EXPECT_EQ(u["work_units"], u["per_seed"][0]["work_units"]);
        EXPECT_EQ(budgeted, !u["memory"]["soft_excess_bytes"].isNull());
    }
}

// The usage block: per seed the walk's own account (the work of its seed phase and both arms,
// the modelled peak), the request's sums; only a request with attempt_id has one, and with
// it the response is otherwise byte for byte the one without
TEST(GraphletAttempt, UsageStatesWhatTheSeedsConsumed) {
    const BudgetCase c = budget_cases()[2];
    const BudgetCase d = budget_cases()[0];
    const std::string strategy = R"({"bounds": {"max_extension_bp": 7}})";
    for (bool budgeted : { false, true }) {
        Json::Value request = request_of({ d.seed, seed_of("AAA") }, strategy, "full");
        if (budgeted)
            request["strategy"]["bounds"]["max_memory_mb"] = 64;
        const Json::Value plain = process_traverse_request(request, *d.anno, "");
        EXPECT_FALSE(plain.isMember("usage"));
        request["attempt_id"] = "cli-1";
        request["budget_id"] = "budget-7";
        const Json::Value with = process_traverse_request(request, *d.anno, "");
        ASSERT_TRUE(with.isMember("usage"));
        Json::Value without = with;
        without.removeMember("usage");
        EXPECT_EQ(compact_json(plain), compact_json(without));
        const Json::Value &u = with["usage"];
        EXPECT_EQ("cli-1", u["attempt_id"].asString());
        EXPECT_EQ("budget-7", u["budget_id"].asString());
        EXPECT_FALSE(u.isMember("locus_id"));
        EXPECT_EQ("completed", u["reason"].asString());
        // the CLI states the bound but enforces none
        EXPECT_FALSE(u["bound"]["enforced"].asBool());
        EXPECT_EQ(2u, u["bound"]["seeds"].asUInt64());
        EXPECT_EQ(2u, u["seeds"]["finished"].asUInt64());
        uint64_t total = 0, held = 0, peak = 0;
        for (Json::ArrayIndex i = 0; i < 2; ++i) {
            const Json::Value &s = u["per_seed"][i];
            Strategy st = parse_traverse_request(request).strategy;
            st.delivery = delivery_costs("full", st.sequences);
            LabelOracle oracle(*d.anno);
            const Seed seed = i ? seed_of("AAA") : d.seed;
            const SeedResult r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
            EXPECT_EQ(r.account.work_used, s["work_units"].asUInt64()) << i;
            EXPECT_EQ(r.account.work_seed, s["work_seed"].asUInt64()) << i;
            EXPECT_EQ(r.account.memory_final, s["final_bytes"].asUInt64()) << i;
            EXPECT_EQ(with["results"][i]["outcome"]["walks"].asString(), s["outcome"].asString());
            EXPECT_EQ(budgeted, !s["soft_excess_bytes"].isNull()) << i;
            total += s["work_units"].asUInt64();
            peak = std::max(peak, held + s["peak_admitted_bytes"].asUInt64());
            held += s["final_bytes"].asUInt64();
        }
        EXPECT_EQ(total, u["work_units"].asUInt64());
        EXPECT_EQ(peak, u["memory"]["peak_admitted_bytes"].asUInt64());
        EXPECT_GT(peak, 0u);
    }
    (void)c;
}

// The statements of a stop from outside the walk fit the widths the delivery model prices
// (effects 640, messages 1 KiB), with the largest numbers they can carry
TEST(GraphletAttempt, StatementsFitTheirWidths) {
    const uint64_t big = std::numeric_limits<uint64_t>::max() / 2;
    for (ExternalStop which : { ExternalStop::CANCELLED, ExternalStop::ATTEMPT_DEADLINE }) {
        for (bool trigger : { true, false }) {
            SeedResult r;
            Strategy st;
            st.max_memory_bytes = big & ~((uint64_t(1) << 20) - 1);
            ResourceStop q;
            q.resource = stop_resource(which);
            q.arm = Arm::RIGHT;
            q.at_bp = big;
            q.limit = static_cast<double>(big);
            q.used = q.demand = static_cast<double>(big);
            r.resource_stop = q;
            r.account.soft_overshoot = big;
            for (ArmResult &arm : r.arms) {
                arm.status = ArmResult::TRUNCATED;
                arm.complete_to_bp = big;
                arm.cap_trigger = CapTrigger{ trigger ? EndReason::RESOURCE_LIMIT
                                                      : EndReason::BEAM_PRUNED,
                                              big, 0, 1, 1, true, static_cast<double>(big) };
                PathResult p;
                p.path_reason = EndReason::RESOURCE_LIMIT;
                arm.paths.push_back(p);
            }
            const Json::Value j = seed_result_to_json(r, st, "summary", false);
            const std::string what = std::string(stop_token(which)) + (trigger ? " trigger" : " later");
            check_widths(j, 2, "", what);
            bool stated = false;
            for (const Json::Value &l : j["arms"]["right"]["limitations"]) {
                stated |= l["knob"].asString() == "attempt_id";
            }
            EXPECT_TRUE(stated) << what;
        }
    }
}


/********* the review of the stage-3 fixes (P2, P3) and of the stage-4 backend (F1) *********/

namespace {

// the roots case of the stage-2 recheck (P2): seed GAT on the repeat index, both arms, trace,
// detail full; its depth-0 state needs about 7 MiB of its own beside the caches' allotments
Json::Value roots_request(uint64_t mb, const std::string &attempt_id = "") {
    Json::Value r;
    Json::Value seed;
    seed["sequence"] = "GAT";
    seed["labels"].append("r1");
    r["seeds"].append(seed);
    r["strategy"] = parse_json(kRepeatStrategy);
    r["strategy"]["direction"] = "both";
    r["strategy"]["bounds"]["max_memory_mb"] = Json::UInt64(mb);
    r["strategy"]["output"]["detail"] = "full";
    r["strategy"]["output"]["timing"] = false;
    if (!attempt_id.empty())
        r["attempt_id"] = attempt_id;
    return r;
}

const Json::Value* memory_walk_domain(const Json::Value &limitations) {
    for (const Json::Value &l : limitations) {
        if (l["kind"].asString() == "walk_domain" && l["knob"].asString() == "bounds.max_memory_mb")
            return &l;
    }
    return nullptr;
}

} // namespace

// P2: a memory stop states the knob value that holds what it needed, not the need at its own
// budget, which includes the caches' allotments of that budget (5/16 of it, at most 128 MiB):
// raised to the need, the seed failed again with a larger one (1 MiB: "need 7 MiB"; 7: 9; 8
// and 9: 10; 10 held it). The value is the smallest whole MiB that holds the need beside its
// own allotments: raised to it the state is held, one MiB less still fails
TEST(GraphletStage3Fixes, MemoryStopsStateTheBudgetThatHoldsThem) {
    constexpr uint64_t kMiB = uint64_t(1) << 20;
    std::mt19937_64 rng(5);
    for (int i = 0; i < 20'000; ++i) {
        const uint64_t budget = rng() % (uint64_t(1) << (10 + rng() % 30)) + 1;
        const uint64_t allotted = memory_allotments(budget);
        const uint64_t need = allotted + rng() % (uint64_t(1) << (rng() % 38));
        const uint64_t b = memory_budget_holding(need, allotted);
        ASSERT_EQ(0u, b % kMiB) << need << " " << allotted;
        if (!allotted) {
            // a budget below 4 bytes allots nothing: the need rounded up
            EXPECT_EQ(std::max<uint64_t>(1, (need + kMiB - 1) / kMiB) * kMiB, b) << need;
            continue;
        }
        EXPECT_LE(need - allotted + memory_allotments(b), b) << need << " " << allotted;
        if (b > kMiB) {
            EXPECT_GT(need - allotted + memory_allotments(b - kMiB), b - kMiB)
                << need << " " << allotted;
        }
    }
    // nothing allotted (the seed phase): the need rounded up, as before
    EXPECT_EQ(3 * kMiB, memory_budget_holding(5 * kMiB / 2, 0));
    EXPECT_EQ(kMiB, memory_budget_holding(1, 0));

    auto anno = repeat_index(100'000, false);
    size_t failures = 0;
    uint64_t holds_from = 0;
    for (uint64_t mb = 1; mb <= 12; ++mb) {
        const Json::Value out = process_traverse_request(roots_request(mb), *anno, "");
        const Json::Value &res = out["results"][0];
        const std::string what = std::to_string(mb) + " MiB";
        if (res["outcome"]["walks"].asString() != "failed") {
            if (!holds_from)
                holds_from = mb;
            continue;
        }
        EXPECT_FALSE(holds_from) << what << ": failed above a budget that held it";
        failures++;
        check_widths(res, 2, "", what);
        const Json::Value *wd = memory_walk_domain(res["limitations"]);
        ASSERT_TRUE(wd) << what;
        const uint64_t knob = (*wd)["observed"].asUInt64();
        EXPECT_GT(knob, mb) << what;
        EXPECT_NE(std::string::npos, res["error"].asString().find(
                "the smallest budget that holds them is " + std::to_string(knob) + " MiB"))
            << what << ": " << res["error"].asString();
        EXPECT_NE(std::string::npos, (*wd)["effect"].asString().find("caches' allotments"))
            << what;
        // raised to the stated value, the depth-0 state is held; one MiB less is not
        const Json::Value raised = process_traverse_request(roots_request(knob), *anno, "");
        EXPECT_NE("failed", raised["results"][0]["outcome"]["walks"].asString()) << what;
        const Json::Value below = process_traverse_request(roots_request(knob - 1), *anno, "");
        EXPECT_EQ("failed", below["results"][0]["outcome"]["walks"].asString()) << what;
    }
    EXPECT_GE(failures, 5u);
    EXPECT_GT(holds_from, 0u);

    // A head's stop states its knob value too (the cap trigger's demand and the walk_domain's
    // observed): raised to it, the walk admits the head that was refused
    size_t heads = 0;
    for (const BudgetCase &c : budget_cases()) {
        for (uint64_t memory = 8192; memory < (uint64_t(2) << 20); memory = memory * 3 / 2) {
            Strategy st = c.st;
            st.max_memory_bytes = memory;
            st.delivery = delivery_costs("full", st.sequences);
            size_t admitted = 0;
            WalkerHooks count;
            count.deny = [&](const Admission&) { admitted++; return false; };
            SeedResult r;
            try {
                r = run_case(c, st, &count);
            } catch (const SeedBudgetError &) {
                continue;
            }
            if (!r.resource_stop || r.resource_stop->resource != ResourceStop::MEMORY
                    || r.resource_stop->cause != ResourceStop::HEAD
                    || r.resource_stop->lower_bound || r.resource_stop->injected)
                continue;
            const ResourceStop &q = *r.resource_stop;
            ASSERT_GT(q.allotted, 0u);
            const uint64_t knob = memory_budget_holding(
                    static_cast<uint64_t>(std::ceil(q.demand)), q.allotted);
            const ArmResult &arm = r.arms[static_cast<size_t>(q.arm)];
            ASSERT_TRUE(arm.cap_trigger) << c.name;
            EXPECT_EQ(static_cast<double>(knob >> 20), arm.cap_trigger->demand) << c.name;
            Strategy raised = st;
            raised.max_memory_bytes = knob;
            size_t admitted_raised = 0;
            WalkerHooks count_raised;
            count_raised.deny = [&](const Admission&) { admitted_raised++; return false; };
            run_case(c, raised, &count_raised);
            EXPECT_GT(admitted_raised, admitted) << c.name << " " << memory;
            heads++;
        }
    }
    EXPECT_GT(heads, 10u);
}

// P3: a derivation found too_wide read its window whole before it was found so: the rows are
// the seed's work (the usage a ledger reconciles), added without a comparison so that
// too_wide stays the cause stated (400,019 units were reported as 0)
TEST(GraphletStage3Fixes, TooWideWindowIsTheSeedsWork) {
    const size_t n = 70'000;
    const uint64_t gat = 8 + 1 + occurrences(n, "GAT"), atg = 8 + 1 + occurrences(n, "ATG");
    for (bool rowdiff : { false, true }) {
        auto anno = repeat_index(n, rowdiff);
        for (uint64_t work : { uint64_t(0), uint64_t(1'000'000'000) }) {
            Strategy st = strategy_of(kRepeatStrategy);
            st.max_seed_labels = 10;     // at most 65,536 entries: both rows are wider
            st.max_work_units = work;
            AttemptMeter meter;
            AttemptControl control;
            control.poll = []() { return ExternalStop::NONE; };
            control.elapsed_ms = []() { return 0.0; };
            control.meter = &meter;
            LabelOracle oracle(*anno);
            const std::string what = std::string(rowdiff ? "row-diff" : "column")
                + (work ? ", work budget" : "");
            try {
                traverse_seed(oracle, seed_of("GATG"), st, LabelChangeCost::forbid(), "", nullptr,
                              &control);
                ADD_FAILURE() << what << ": the derivation was not too wide";
            } catch (const SeedDerivationError &e) {
                EXPECT_EQ(SeedDerivationError::TOO_WIDE, e.cause()) << what;
                if (rowdiff && work) {
                    // each row with its row-diff dependency rows (the budget-aware reads)
                    EXPECT_GE(meter.work_seed, gat + atg) << what;
                } else {
                    EXPECT_EQ(gat + atg, meter.work_seed) << what;
                }
                EXPECT_EQ(meter.work_seed, meter.work_units) << what;
            }
        }
    }
}

// F1 (review of the stage-4 backend): the usage of a failed or never started seed is what its
// result states and holds — its memory_bound_soft (the echo of seed_id included), the failed
// result as priced until the response is written, the demand a memory budget refused — and an
// admitted peak never above the budget; held_bound_bytes bounds the request as observed
TEST(GraphletAttempt, UsageStatesWhatFailedResultsHold) {
    constexpr uint64_t kMiB = uint64_t(1) << 20;
    const BudgetCase c = budget_cases()[0];
    const std::string strategy = R"({"bounds": {"max_extension_bp": 7, "max_memory_mb": 1}})";
    const std::string big_a(2 * kMiB, 'z'), big_b(2 * kMiB, 'q');
    std::vector<Seed> seeds { c.seed, c.seed, c.seed, c.seed };
    seeds[0].seed_id = big_a;   // failed at depth 0: its echo alone exceeds the budget
    seeds[3].seed_id = big_b;   // never started: its failed result's echo exceeds it
    // the polls of the first two seeds, the poll before each included: the next is the one
    // before the third seed
    Json::Value two = request_of({ seeds[0], seeds[1] }, strategy, "full");
    for (Json::ArrayIndex i = 0; i < 2; ++i) {
        two["seeds"][i]["seed_id"] = seeds[i].seed_id;
    }
    CancelAtPoll count(std::numeric_limits<uint64_t>::max());
    process_traverse_request(two, *c.anno, "", {}, nullptr, count.attempt.get());
    Json::Value request = request_of(seeds, strategy, "full");
    for (Json::ArrayIndex i = 0; i < 4; ++i) {
        if (!seeds[i].seed_id.empty())
            request["seeds"][i]["seed_id"] = seeds[i].seed_id;
    }
    CancelAtPoll cancel(count.calls + 1);
    const Json::Value out = process_traverse_request(request, *c.anno, "", {}, nullptr,
                                                     cancel.attempt.get());
    ASSERT_EQ(4u, out["results"].size());
    const Json::Value &u = out["usage"];
    EXPECT_EQ("failed", out["results"][0]["outcome"]["walks"].asString());
    EXPECT_NE("failed", out["results"][1]["outcome"]["walks"].asString());
    EXPECT_EQ("not_started", out["results"][2]["resource_stop"]["phase"].asString());
    EXPECT_EQ("not_started", out["results"][3]["resource_stop"]["phase"].asString());
    EXPECT_EQ(4u, u["seeds"]["requested"].asUInt64());
    EXPECT_EQ(2u, u["seeds"]["started"].asUInt64());
    EXPECT_EQ(2u, u["seeds"]["finished"].asUInt64());
    EXPECT_EQ(0u, u["seeds"]["abandoned"].asUInt64());

    Strategy st = parse_traverse_request(request).strategy;
    st.delivery = delivery_costs("full", st.sequences);
    auto priced = [&](const std::string &seed_id) {
        return st.delivery.fixed + seed_id.size()
            + (st.delivery.name ? st.delivery.name(seed_id, DeliveryCosts::Name::SEED_ID) : 0);
    };
    uint64_t held = 0, peak = 0, bound = 0, soft = 0;
    for (Json::ArrayIndex i = 0; i < 4; ++i) {
        const Json::Value &s = u["per_seed"][i];
        const Json::Value &res = out["results"][i];
        const std::string what = "seed " + std::to_string(i) + ": " + compact_json(s);
        // what the result's memory_bound_soft states, in bytes
        uint64_t stated = 0;
        for (const Json::Value &l : res["limitations"]) {
            if (l["kind"].asString() == "memory_bound_soft")
                stated = l["observed"].asUInt64();
        }
        EXPECT_EQ(stated, (s["soft_excess_bytes"].asUInt64() + kMiB - 1) / kMiB) << what;
        EXPECT_LE(s["peak_admitted_bytes"].asUInt64(), kMiB) << what;
        if (res["outcome"]["walks"].asString() == "failed") {
            EXPECT_EQ(priced(seeds[i].seed_id), s["final_bytes"].asUInt64()) << what;
        }
        const bool walked = i < 2;
        if (walked) {
            peak = std::max(peak, held + s["peak_admitted_bytes"].asUInt64());
            bound = std::max(bound, held + kMiB + s["soft_excess_bytes"].asUInt64());
        }
        soft = std::max(soft, s["soft_excess_bytes"].asUInt64());
        held += s["final_bytes"].asUInt64();
        peak = std::max(peak, held);
        bound = std::max(bound, held);
    }
    // the depth-0 failure: the refused demand, above the budget; the echo is the excess
    EXPECT_GT(u["per_seed"][0]["refused_bytes"].asUInt64(), kMiB);
    EXPECT_GT(u["per_seed"][0]["soft_excess_bytes"].asUInt64(), kMiB);
    EXPECT_TRUE(u["per_seed"][2]["refused_bytes"].isNull());
    EXPECT_GT(u["per_seed"][3]["soft_excess_bytes"].asUInt64(), kMiB);
    EXPECT_EQ(peak, u["memory"]["peak_admitted_bytes"].asUInt64());
    EXPECT_EQ(soft, u["memory"]["soft_excess_bytes"].asUInt64());
    EXPECT_EQ(bound, u["memory"]["held_bound_bytes"].asUInt64());

    // without a memory budget the account is no bound of what was held: no bound is stated
    Json::Value plain = request_of({ c.seed }, R"({"bounds": {"max_extension_bp": 7}})", "full");
    plain["attempt_id"] = "plain-1";
    const Json::Value pu = process_traverse_request(plain, *c.anno, "")["usage"];
    EXPECT_TRUE(pu["memory"]["held_bound_bytes"].isNull());
    EXPECT_TRUE(pu["memory"]["soft_excess_bytes"].isNull());
    EXPECT_TRUE(pu["per_seed"][0]["refused_bytes"].isNull());
}


// ---------------------------------------------------------------- record coordinates (§18)

namespace {

// The opt-out response from an opt-in one: without the two echo fields, each result's block
// (or null and reason), the `coordinates` limitation and, in a graphlet, the one K record of
// kind coordinates (the Z count, graphlet_lines and graphlet_bytes adjusted). |k_records|
// counts the K records removed
Json::Value strip_coordinates(Json::Value out, size_t *k_records = nullptr) {
    Json::Value &o = out["strategy"]["output"];
    o.removeMember("coordinates");
    o.removeMember("max_coordinate_occurrences");
    for (Json::Value &r : out["results"]) {
        r.removeMember("coordinates");
        r.removeMember("coordinates_reason");
        Json::Value kept(Json::arrayValue);
        for (const Json::Value &l : r["limitations"]) {
            if (l["kind"].asString() != "coordinates")
                kept.append(l);
        }
        r["limitations"] = kept;
        if (!r.isMember("graphlet"))
            continue;
        std::vector<std::string> lines = split_on(r["graphlet"].asString(), '\n');
        std::string text;
        size_t removed = 0;
        for (const std::string &line : lines) {
            if (line.rfind("K * coordinates ", 0) == 0) {
                removed++;
                continue;
            }
            if (line.empty())
                continue;
            text += line.rfind("Z ", 0) == 0 ? "Z " + std::to_string(u(line.substr(2)) - removed)
                                             : line;
            text += '\n';
        }
        r["graphlet"] = text;
        r["graphlet_bytes"] = Json::UInt64(text.size());
        r["graphlet_lines"] = Json::UInt64(r["graphlet_lines"].asUInt64() - removed);
        if (k_records)
            *k_records += removed;
    }
    return out;
}

// a column-label coordinate index: S twice in C (S·P·S·Q), once in D (S·P); a stretch no label
// carries whole (M·N in E, N·O in G)
struct CoordIndex {
    std::string S, P, Q, M, N, O;
    std::shared_ptr<graph::AnnotatedDBG> anno;
    CoordIndex() {
        for (uint32_t t = 0; ; ++t) {
            S = random_seq(30, 2101 + t); P = random_seq(25, 2201 + t); Q = random_seq(25, 2301 + t);
            if (P[0] != Q[0])
                break;
        }
        M = random_seq(25, 2401); N = random_seq(25, 2402); O = random_seq(25, 2403);
        anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                11, { S + P + S + Q, S + P, M + N, N + O }, { "C", "D", "E", "G" },
                DeBruijnGraph::BASIC, true, { 0, 0, 0, 0 });
    }
    // no label carries it whole: a failed derivation
    std::string orphan() const { return M.substr(10) + N + O.substr(0, 10); }
};

const char *const kTrace = R"({"support": "trace", "branching": {"on_reconverge": "keep",
                                  "max_label_branches": 2}, "bounds": {"max_extension_bp": 60}})";

Json::Value with_coordinates(Json::Value request, Json::Value cap = Json::Value()) {
    request["strategy"]["output"]["coordinates"] = true;
    if (!cap.isNull())
        request["strategy"]["output"]["max_coordinate_occurrences"] = cap;
    return request;
}

} // namespace

// The cap is refused without coordinates: true (it would change nothing: §7.0, decision C-N6),
// accepted with it as an integer of at least 1 or "unlimited" (default 16), and inert under
// support kmer (C-N7)
TEST(GraphletCoordinates, CapWithoutCoordinatesIsRefused) {
    auto parse = [](const std::string &output, const std::string &support = "trace") {
        return parse_traverse_request(parse_json(
                R"({"seeds": [{"sequence": "ACGTACGTACGTACGT"}], "strategy": {"support": ")"
                + support + R"(", "branching": {"on_reconverge": "keep"}, "output": )" + output + "}}"));
    };
    auto refused = [&](const std::string &output, const std::string &says) {
        try {
            parse(output);
            ADD_FAILURE() << output << " was accepted";
        } catch (const InvalidRequest &e) {
            EXPECT_NE(std::string::npos, std::string(e.what()).find(says)) << output << ": " << e.what();
            EXPECT_NE(std::string::npos,
                      std::string(e.what()).find("strategy.output.max_coordinate_occurrences")
                      + std::string(e.what()).find("strategy.output.coordinates") + 1)
                << e.what();
        }
    };
    refused(R"({"max_coordinate_occurrences": 4})", "only with output.coordinates true");
    refused(R"({"coordinates": false, "max_coordinate_occurrences": "unlimited"})",
            "only with output.coordinates true");
    // 0 and a number above the maximum state the same whole range, "unlimited" included: [0, ...]
    // misstated it, and a client following it was refused again (review of W1, finding 4)
    const std::string range = "out of range [1, " + std::to_string(Strategy::kUnlimited - 1)
                            + "] (or \"unlimited\")";
    refused(R"({"coordinates": true, "max_coordinate_occurrences": 0})", range);
    refused(R"({"coordinates": true, "max_coordinate_occurrences": 18446744073709551615})", range);
    refused(R"({"coordinates": true, "max_coordinate_occurrences": "many"})",
            "expected an integer or \"unlimited\"");
    refused(R"({"coordinates": true, "max_coordinate_occurrences": -1})", "non-negative");
    refused(R"({"coordinates": "yes"})", "expected a boolean");
    EXPECT_FALSE(parse(R"({})").strategy.coordinates);
    EXPECT_FALSE(parse(R"({"coordinates": false})").strategy.coordinates);
    const TraverseRequest def = parse(R"({"coordinates": true})");
    EXPECT_TRUE(def.strategy.coordinates);
    EXPECT_EQ(16u, def.strategy.max_coordinate_occurrences);
    EXPECT_EQ(Strategy::kUnlimited,
              parse(R"({"coordinates": true, "max_coordinate_occurrences": "unlimited"})")
                  .strategy.max_coordinate_occurrences);
    EXPECT_EQ(1u, parse(R"({"coordinates": true, "max_coordinate_occurrences": 1})")
                      .strategy.max_coordinate_occurrences);
    // inert under support kmer, accepted: a request that switches its support stays valid
    EXPECT_EQ(3u, parse(R"({"coordinates": true, "max_coordinate_occurrences": 3})", "kmer")
                      .strategy.max_coordinate_occurrences);
}

// Without coordinates nothing is added (decision C1): no echo field, no block, no reason. With
// them the echo carries both fields (the cap's default included), and is resubmittable: the
// echoed strategy gives the same response
TEST(GraphletCoordinates, EchoOnlyWhenRequestedAndResubmittable) {
    const CoordIndex ix;
    const Json::Value plain = request_of({ seed_of(ix.S, { "C", "D" }) }, kTrace, "full");
    const Json::Value off = process_traverse_request(plain, *ix.anno, "");
    EXPECT_FALSE(off["strategy"]["output"].isMember("coordinates"));
    EXPECT_FALSE(off["strategy"]["output"].isMember("max_coordinate_occurrences"));
    EXPECT_FALSE(off["results"][0].isMember("coordinates"));
    EXPECT_FALSE(off["results"][0].isMember("coordinates_reason"));
    const Json::Value on = process_traverse_request(with_coordinates(plain), *ix.anno, "");
    EXPECT_TRUE(on["strategy"]["output"]["coordinates"].asBool());
    EXPECT_EQ(16u, on["strategy"]["output"]["max_coordinate_occurrences"].asUInt64());
    EXPECT_TRUE(on["results"][0]["coordinates"].isObject());
    Json::Value again = plain;
    again["strategy"] = on["strategy"];
    again["strategy"].removeMember("clamped");
    EXPECT_EQ(compact_json(on), compact_json(process_traverse_request(again, *ix.anno, "")));
    const Json::Value unlimited = process_traverse_request(
            with_coordinates(plain, "unlimited"), *ix.anno, "");
    EXPECT_EQ("unlimited", unlimited["strategy"]["output"]["max_coordinate_occurrences"].asString());
    EXPECT_EQ("unlimited", unlimited["results"][0]["coordinates"]["max_occurrences"].asString());
}

// coordinates: null with the reason, in the order of §18.1: the index's ("index has no
// coordinates"), the support's ("support kmer", annotate mode included), then "no traversal"
// for a seed without a walk — a failed derivation, a seed a budget does not hold, a refused
// non-UTF-8 dictionary, a seed never started — and "partial derivation" for a trace seed whose
// permitted set was derived from part of it (D3)
TEST(GraphletCoordinates, NullReasons) {
    const CoordIndex ix;
    auto reason = [](const Json::Value &result) {
        EXPECT_TRUE(result.isMember("coordinates")) << compact_json(result);
        EXPECT_TRUE(result["coordinates"].isNull()) << compact_json(result);
        return result["coordinates_reason"].asString();
    };
    for (const char *detail : { "full", "graphlet" }) {
        // support kmer, and annotate mode: whatever the seed
        Json::Value r = with_coordinates(request_of({ seed_of(ix.S), seed_of(ix.orphan()) },
                                                    R"({"bounds": {"max_extension_bp": 30}})", detail));
        Json::Value out = process_traverse_request(r, *ix.anno, "");
        EXPECT_EQ("support kmer", reason(out["results"][0])) << detail;
        EXPECT_EQ("support kmer", reason(out["results"][1])) << detail;   // failed too
        r["strategy"]["labels"]["mode"] = "annotate";
        r["strategy"]["branching"]["on_reconverge"] = "keep";
        out = process_traverse_request(r, *ix.anno, "");
        EXPECT_EQ("support kmer", reason(out["results"][0])) << detail;
        // trace: a walked seed has the block, a failed derivation none
        r = with_coordinates(request_of({ seed_of(ix.S), seed_of(ix.orphan()) }, kTrace, detail));
        out = process_traverse_request(r, *ix.anno, "");
        EXPECT_TRUE(out["results"][0]["coordinates"].isObject()) << detail;
        EXPECT_FALSE(out["results"][0].isMember("coordinates_reason")) << detail;
        EXPECT_EQ("failed", out["results"][1]["outcome"]["walks"].asString());
        EXPECT_EQ("no traversal", reason(out["results"][1])) << detail;
        // a budget that does not hold the seed (its depth-0 state at 1 MiB in detail full,
        // with 64 labels)
        // ... and the time budget spent after the derivation's first k-mer (D3): under trace
        // its set is not followed over the seed, no occurrences are known
        TraverseLimits late;
        late.max_time_ms = 1e-9;
        r = with_coordinates(request_of({ seed_of(ix.S) }, kTrace, detail));
        out = process_traverse_request(r, *ix.anno, "", late);
        EXPECT_EQ("partial", out["results"][0]["outcome"]["walks"].asString()) << detail;
        EXPECT_EQ("partial derivation", reason(out["results"][0])) << detail;
    }
    // the index's reason comes first, under every support and mode
    auto plain = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, { ix.S + ix.P + ix.S + ix.Q, ix.S + ix.P, ix.M + ix.N, ix.N + ix.O },
            { "C", "D", "E", "G" }, DeBruijnGraph::BASIC);
    for (const char *strategy : { R"({"bounds": {"max_extension_bp": 30}})",
                                  R"({"labels": {"mode": "annotate"}, "bounds": {"max_extension_bp": 30}})" }) {
        const Json::Value out = process_traverse_request(
                with_coordinates(request_of({ seed_of(ix.S), seed_of(ix.orphan()) }, strategy)),
                *plain, "");
        EXPECT_EQ("index has no coordinates", reason(out["results"][0])) << strategy;
        EXPECT_EQ("index has no coordinates", reason(out["results"][1])) << strategy;
    }
    // a seed the memory budget does not hold: 200 labels at the seed's depth-0 state
    {
        const std::string X = random_seq(30, 2501), Y = random_seq(30, 2502);
        std::vector<std::string> seqs, labels;
        for (size_t i = 0; i < 200; ++i) {
            seqs.push_back(X + Y);
            labels.push_back("L" + std::to_string(i));
        }
        auto wide = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                11, seqs, labels, DeBruijnGraph::BASIC, true, std::vector<uint64_t>(200, 0));
        Json::Value r = with_coordinates(request_of({ seed_of(X) }, kTrace, "full"));
        r["strategy"]["bounds"]["max_memory_mb"] = 1;
        const Json::Value out = process_traverse_request(r, *wide, "");
        ASSERT_TRUE(out["results"][0].isMember("error")) << compact_json(out["results"][0]);
        EXPECT_EQ("memory", out["results"][0]["resource_stop"]["resource"].asString());
        EXPECT_EQ("no traversal", reason(out["results"][0]));
        // the depth-0 state holds the seed's and the roots' coordinates: dropping them helps
        std::set<std::string> actions;
        for (const Json::Value &a : out["results"][0]["resource_stop"]["actions"]) actions.insert(a.asString());
        EXPECT_TRUE(actions.count("drop_coordinates")) << compact_json(out["results"][0]["resource_stop"]);
    }
    // a refused dictionary (a name that is not UTF-8)
    {
        const std::string R = random_seq(40, 2601), T = random_seq(30, 2602);
        auto bad = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
                11, { R + T }, { "bad\xFF" }, DeBruijnGraph::BASIC, true, { 0 });
        const Json::Value out = process_traverse_request(
                with_coordinates(request_of({ seed_of(R) }, kTrace, "full")), *bad, "");
        EXPECT_EQ("unrepresentable_label_name", out["results"][0]["limitations"][0]["cause"].asString());
        EXPECT_EQ("no traversal", reason(out["results"][0]));
    }
    // seeds never started (the attempt cancelled between seeds)
    {
        const Json::Value one = with_coordinates(request_of({ seed_of(ix.S, { "C" }) }, kTrace));
        CancelAtPoll single(std::numeric_limits<uint64_t>::max());
        process_traverse_request(one, *ix.anno, "", {}, nullptr, single.attempt.get());
        Json::Value two = one;
        two["seeds"].append(one["seeds"][0]);
        CancelAtPoll cancel(single.calls + 1);
        const Json::Value out = process_traverse_request(two, *ix.anno, "", {}, nullptr,
                                                         cancel.attempt.get());
        ASSERT_EQ(2u, out["results"].size());
        EXPECT_TRUE(out["results"][0]["coordinates"].isObject());
        EXPECT_EQ("not_started", out["results"][1]["resource_stop"]["phase"].asString());
        EXPECT_EQ("no traversal", reason(out["results"][1]));
    }
}

// One block in every detail (the envelope carries it; the graphlet's summary too)
TEST(GraphletCoordinates, SameBlockInEveryDetail) {
    const CoordIndex ix;
    for (const Json::Value &cap : { Json::Value(1), Json::Value(16), Json::Value("unlimited") }) {
        std::string block;
        for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
            const Json::Value out = process_traverse_request(
                    with_coordinates(request_of({ seed_of(ix.S, { "C", "D" }) }, kTrace, detail), cap),
                    *ix.anno, "");
            const Json::Value &c = out["results"][0]["coordinates"];
            ASSERT_TRUE(c.isObject()) << detail;
            if (block.empty()) {
                block = compact_json(c);
            } else {
                EXPECT_EQ(block, compact_json(c)) << detail << " cap " << compact_json(cap);
            }
        }
        EXPECT_NE(std::string::npos, block.find("\"kind\":\"column\""));
    }
}

// The MGT body changes only where a list was cut, by one K record of kind coordinates (C12's free
// token, extra field lists_cut), which the reader keeps; otherwise it is byte-identical to the
// opt-out body. Stripped of the block, its reason, the limitation and the echo, an opt-in response
// is the opt-out response, in every detail
TEST(GraphletCoordinates, KRecordOnlyWhenCutAndStrippedEqualsOptOut) {
    const CoordIndex ix;
    for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
        for (const std::string &strategy : { std::string(kTrace),
                                             std::string(R"({"bounds": {"max_extension_bp": 30}})") }) {
            const Json::Value plain = request_of({ seed_of(ix.S, { "C", "D" }),
                                                   seed_of(ix.orphan()) }, strategy, detail);
            const std::string off = compact_json(process_traverse_request(plain, *ix.anno, ""));
            for (const Json::Value &cap : { Json::Value(1), Json::Value(16), Json::Value("unlimited") }) {
                const Json::Value on = process_traverse_request(with_coordinates(plain, cap), *ix.anno, "");
                size_t k_records = 0;
                EXPECT_EQ(off, compact_json(strip_coordinates(on, &k_records)))
                    << detail << " " << strategy << " cap " << compact_json(cap);
                const Json::Value &first = on["results"][0];
                const bool cut = first.isMember("coordinates") && first["coordinates"].isObject()
                              && !first["coordinates"]["complete"].asBool();
                const Json::Value *lim = limitation_of(first["limitations"], "coordinates");
                // S occurs twice in C: only a cap of 1 cuts (the seed's list and C's root runs)
                const bool expect_cut = strategy == kTrace && cap == Json::Value(1);
                EXPECT_EQ(expect_cut, cut) << detail << " cap " << compact_json(cap);
                EXPECT_EQ(expect_cut, lim != nullptr) << detail;
                EXPECT_EQ(expect_cut && std::string(detail) == "graphlet" ? 1u : 0u, k_records);
                if (lim) {
                    EXPECT_EQ("output.max_coordinate_occurrences", (*lim)["knob"].asString());
                    EXPECT_EQ(1u, (*lim)["limit"].asUInt64());
                    EXPECT_EQ(2u, (*lim)["observed"].asUInt64());       // the largest true count
                    EXPECT_GE((*lim)["lists_cut"].asUInt64(), 1u);
                    // in no outcome class (decision C2): the outcome is the opt-out's (the
                    // stripped response above equals it, outcome included)
                    if (std::string(detail) == "graphlet") {
                        const std::string text = first["graphlet"].asString();
                        const size_t at = text.find("\nK * coordinates output.max_coordinate_occurrences i:1 i:2 * lists_cut=i:");
                        EXPECT_NE(std::string::npos, at) << text;
                    }
                }
            }
        }
    }
}

// The delivery model bounds the block too (DeliveryCostsBoundTheOutput's check): many occurrences
// with "unlimited", positions of 20 digits, every detail; the walker charges at least the model
namespace {

std::shared_ptr<graph::AnnotatedDBG> repeat_coordinate_index(const std::string &S, uint64_t start) {
    std::string rec;
    for (size_t i = 0; i < 40; ++i) {
        rec += S + random_seq(12, 2700 + i);
    }
    return test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, { rec, S + random_seq(30, 2799) }, { "F", "H" }, DeBruijnGraph::BASIC, true,
            { start, start });
}

} // namespace

TEST(GraphletCoordinates, DeliveryCostsBoundTheBlock) {
    const std::string S = random_seq(30, 2700);
    size_t checked = 0;
    for (uint64_t start : { uint64_t(0), uint64_t(10'000'000'000'000'000'000ull) }) {
        auto anno = repeat_coordinate_index(S, start);
        LabelOracle oracle(*anno);
        for (size_t cap : { size_t(16), Strategy::kUnlimited }) {
            for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
                const std::string what = std::string(detail) + " cap " + std::to_string(cap)
                                       + " start " + std::to_string(start);
                Strategy st = strategy_of(kTrace);
                st.coordinates = true;
                st.max_coordinate_occurrences = cap;
                st.delivery = delivery_costs(detail, st.sequences, kMgtFloatWidth,
                                             CoordinatesOutput::BLOCK);
                const SeedResult r = traverse_seed(oracle, seed_of(S, { "F", "H" }), st,
                                                   LabelChangeCost::forbid());
                ASSERT_TRUE(r.coordinates_recorded) << what;
                const DeliveryCosts &d = st.delivery;
                const uint64_t model = modelled_delivery(r, d);
                Json::Value j = seed_result_to_json(r, st, detail, false);
                ASSERT_TRUE(j["coordinates"].isObject()) << what;
                uint64_t before = 0;
                if (std::string(detail) == "graphlet") {
                    const std::string text = graphlet_text(r, seed_of(S, { "F", "H" }), st,
                                                           context_of(oracle), j);
                    uint64_t index = 0;
                    for (const ArmResult &arm : r.arms) {
                        index += 40 * arm.segments.size() + 8 * arm.runs.size();
                        for (const Segment &s : arm.segments) {
                            for (const Event &ev : s.events) index += 96 * (ev.type == EventType::LABEL_END);
                        }
                    }
                    before = 3 * text.size() + index + json_tree_bytes(j);
                    j["graphlet"] = text;
                }
                check_widths(j, 2, "", what);
                EXPECT_LE(j["limitations"].size(), 6u) << what;
                const uint64_t actual = delivery_peak(j, before);
                EXPECT_LE(actual, model) << what;
                EXPECT_LE(model, r.account.memory_final) << what;
                EXPECT_EQ(40u, r.seed_coordinates[0].total) << what;
                checked++;
            }
        }
    }
    EXPECT_EQ(16u, checked);
}

// A memory stop of a request asking for coordinates offers drop_coordinates (C12): the recorded
// ones are part of the account (4-20 times a run's own cost), and so is the null form with its
// reason where none can be recorded (support kmer), which moved such stops shallower with no
// lever naming the cause (review of W1, finding 2); dropping them is then a stop no shallower.
// Without coordinates it is not offered
TEST(GraphletCoordinates, MemoryStopOffersDropCoordinates) {
    const std::string S = random_seq(30, 2800);
    std::vector<std::string> seqs, labels;
    std::vector<uint64_t> starts;
    uint64_t at = 0;
    for (size_t i = 0; i < 128; ++i) {
        seqs.push_back(S + random_seq(200, 2900 + i));
        labels.push_back("F");
        starts.push_back(at);
        at += seqs.back().size() - 11 + 1;
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, seqs, labels, DeBruijnGraph::BASIC, true, starts);
    size_t offered[2] = { 0, 0 };
    for (uint64_t mb : { 1, 2, 3, 4 }) {
        for (bool trace : { true, false }) {
            Json::Value r = request_of({ seed_of(S, { "F" }) },
                                       trace ? R"({"support": "trace", "direction": "right",
                                                   "branching": {"on_reconverge": "keep",
                                                                 "max_label_branches": "unlimited",
                                                                 "max_splits_per_path": "unlimited"},
                                                   "bounds": {"max_extension_bp": 300}})"
                                             : R"({"direction": "right",
                                                   "branching": {"on_reconverge": "keep",
                                                                 "max_label_branches": "unlimited",
                                                                 "max_splits_per_path": "unlimited"},
                                                   "bounds": {"max_extension_bp": 300}})", "full");
            r["strategy"]["bounds"]["max_memory_mb"] = Json::UInt64(mb);
            uint64_t depth[2] = { 0, 0 };
            for (bool coordinates : { true, false }) {
                const Json::Value res = process_traverse_request(
                        coordinates ? with_coordinates(r) : r, *anno, "")["results"][0];
                depth[coordinates] = res["arms"]["right"]["complete_to_bp"].asUInt64();
                const Json::Value &q = res["resource_stop"];
                if (!q.isObject() || q["resource"].asString() != "memory")
                    continue;
                std::set<std::string> actions;
                for (const Json::Value &a : q["actions"]) actions.insert(a.asString());
                EXPECT_EQ(coordinates, actions.count("drop_coordinates") == 1)
                    << mb << " MiB trace " << trace << " coordinates " << coordinates;
                offered[trace] += actions.count("drop_coordinates");
                if (coordinates && !trace) {
                    // inert: the null form and its reason are what the opt-in adds
                    EXPECT_TRUE(res["coordinates"].isNull()) << mb;
                    EXPECT_EQ("support kmer", res["coordinates_reason"].asString()) << mb;
                }
            }
            // the lever holds: without coordinates the walk stops no shallower
            EXPECT_GE(depth[false], depth[true]) << mb << " MiB trace " << trace;
        }
    }
    EXPECT_GT(offered[true], 0u);
    EXPECT_GT(offered[false], 0u);
}

namespace {

// One result as the server builds it from a walk (seed_result_to_json, and in a graphlet the
// body with its byte and line counts): |rj|'s graphlet members from |r| and |rj| itself
void attach_graphlet(Json::Value *rj, const SeedResult &r, const Seed &seed, const Strategy &st,
                     const LabelOracle &oracle) {
    size_t lines = 0;
    const std::string text = graphlet_text(r, seed, st, context_of(oracle), *rj, &lines);
    (*rj)["graphlet_bytes"] = Json::UInt64(text.size());
    (*rj)["graphlet_lines"] = Json::UInt64(lines);
    (*rj)["graphlet"] = text;
}

// |rj| without what record coordinates add to it: the block (or null and reason), the cut
// list's limitation and the drop_coordinates action; a graphlet's body written again from the
// stripped JSON (its K and Q records follow it) — the text the same walk gives without them
Json::Value strip_result(Json::Value rj, const SeedResult &r, const Seed &seed,
                         const Strategy &st, const LabelOracle &oracle) {
    rj.removeMember("coordinates");
    rj.removeMember("coordinates_reason");
    Json::Value kept(Json::arrayValue);
    for (const Json::Value &l : rj["limitations"]) {
        if (l["kind"].asString() != "coordinates")
            kept.append(l);
    }
    rj["limitations"] = kept;
    if (rj.isMember("resource_stop")) {
        Json::Value actions(Json::arrayValue);
        for (const Json::Value &a : rj["resource_stop"]["actions"]) {
            if (a.asString() != "drop_coordinates")
                actions.append(a);
        }
        rj["resource_stop"]["actions"] = actions;
    }
    if (rj.isMember("graphlet"))
        attach_graphlet(&rj, r, seed, st, oracle);
    return rj;
}

} // namespace

// Plan revision 3: the text record coordinates add to a seed's result is counted exactly from
// digits and punctuation (coordinates_text_bytes), whatever the detail, the cap, the digits of
// the positions (20 here), a cut list, a memory stop offering drop_coordinates or the null
// form: the result's text less it is the text of the same walk without coordinates. And the
// account's coordinate share (ResourceAccount::coordinates: DeliveryCosts::coordinate_fixed
// and the entries) prices that text at kCoordinateAccountPerTextByte at least, which the
// delivery reserve relies on; without coordinates both are 0. compact_json_size is the
// writer's length of every output here
TEST(GraphletCoordinates, CoordinateTextIsExactAndBoundedByItsAccount) {
    const std::string S = random_seq(30, 2700);
    size_t checked = 0, cut = 0, dropped = 0, nulls = 0;
    for (uint64_t start : { uint64_t(0), uint64_t(10'000'000'000'000'000'000ull) }) {
        auto anno = repeat_coordinate_index(S, start);
        for (size_t cap : { size_t(1), size_t(16), Strategy::kUnlimited }) {
            for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
                for (uint64_t mb : { uint64_t(0), uint64_t(1), uint64_t(2) }) {
                    for (bool trace : { true, false }) {
                        // fresh per walk: the oracle's counters are the result's
                        LabelOracle oracle(*anno);
                        const std::string what = std::string(detail) + " cap " + std::to_string(cap)
                                               + " start " + std::to_string(start) + " "
                                               + std::to_string(mb) + " MiB"
                                               + (trace ? " trace" : " kmer");
                        Strategy st = strategy_of(trace ? kTrace : R"({"branching": {"on_reconverge": "keep", "max_label_branches": 2}, "bounds": {"max_extension_bp": 60}})");
                        st.coordinates = true;
                        st.max_coordinate_occurrences = cap;
                        st.max_memory_bytes = mb << 20;
                        st.delivery = delivery_costs(detail, st.sequences, kMgtFloatWidth,
                                                     trace ? CoordinatesOutput::BLOCK
                                                           : CoordinatesOutput::REASON);
                        const Seed seed = seed_of(S, { "F", "H" });
                        SeedResult r;
                        try {
                            r = traverse_seed(oracle, seed, st, LabelChangeCost::forbid());
                        } catch (const SeedBudgetError &) {
                            continue;       // its failed result has no walk: priced apart
                        }
                        Json::Value rj = seed_result_to_json(r, st, detail, false);
                        if (std::string(detail) == "graphlet")
                            attach_graphlet(&rj, r, seed, st, oracle);
                        const std::string text = json_text(rj, true);
                        EXPECT_EQ(text.size(), compact_json_size(rj)) << what;
                        const uint64_t coordinate_text = coordinates_text_bytes(rj);
                        const Json::Value stripped = strip_result(rj, r, seed, st, oracle);
                        EXPECT_EQ(text.size() - coordinate_text, json_text(stripped, true).size())
                            << what;
                        EXPECT_GT(r.account.coordinates, 0u) << what;
                        EXPECT_GE(r.account.coordinates, st.delivery.coordinate_fixed) << what;
                        EXPECT_LE(coordinate_text * kCoordinateAccountPerTextByte,
                                  r.account.coordinates) << what;
                        // and the opt-out walk's account is the rest, unbudgeted (the same walk)
                        if (!mb) {
                            Strategy off = st;
                            off.coordinates = false;
                            off.max_coordinate_occurrences = 16;
                            off.delivery = delivery_costs(detail, off.sequences);
                            LabelOracle fresh(*anno);
                            const SeedResult o = traverse_seed(fresh, seed, off,
                                                               LabelChangeCost::forbid());
                            EXPECT_EQ(0u, o.account.coordinates) << what;
                            EXPECT_EQ(o.account.memory_final,
                                      r.account.memory_final - r.account.coordinates) << what;
                            Json::Value oj = seed_result_to_json(o, off, detail, false);
                            if (std::string(detail) == "graphlet")
                                attach_graphlet(&oj, o, seed, off, fresh);
                            EXPECT_EQ(json_text(oj, true), json_text(stripped, true)) << what;
                            EXPECT_EQ(0u, coordinates_text_bytes(oj)) << what;
                        }
                        cut += limitation_of(rj["limitations"], "coordinates") != nullptr;
                        nulls += rj["coordinates"].isNull();
                        if (rj.isMember("resource_stop")) {
                            for (const Json::Value &a : rj["resource_stop"]["actions"]) {
                                dropped += a.asString() == "drop_coordinates";
                            }
                        }
                        checked++;
                    }
                }
            }
        }
    }
    std::cerr << checked << " results, " << cut << " cut, " << dropped << " offering "
              << "drop_coordinates, " << nulls << " null" << std::endl;
    EXPECT_GT(checked, 100u);
    EXPECT_GT(cut, 0u);
    EXPECT_GT(dropped, 0u);
    EXPECT_GT(nulls, 0u);
}

// kCoordinateAccountPerTextByte is a bound, not an estimate: each part of the coordinate share
// is priced at least 12 times the most text it can write — an occurrence of two 20-digit
// numbers, a run's entry and a seed label's with every optional member at its widest, the
// block's skeleton with a cut list's limitation (a 20-digit count, "unlimited"), the K record
// the writer makes of it, the drop_coordinates tokens and the digits they add to the counts,
// and the null form with its longest reason — in every detail
TEST(GraphletCoordinates, CoordinateAccountBoundsItsText) {
    const uint64_t kMax = std::numeric_limits<uint64_t>::max();
    // what a member or an element of a list adds: itself and its comma
    auto element = [](const Json::Value &v) { return compact_json_size(v) + 1; };
    auto member = [](const std::string &key, const Json::Value &v) {
        return key.size() + 2 + 1 + compact_json_size(v) + 1;
    };
    Json::Value occurrence(Json::arrayValue);
    occurrence.append(Json::UInt64(kMax));
    occurrence.append(Json::UInt64(kMax));
    Json::Value run;
    run["run"] = Json::UInt64(kMax);
    run["label"] = Json::UInt(std::numeric_limits<uint32_t>::max());
    run["from_bp"] = Json::UInt64(kMax);
    run["to_bp"] = Json::UInt64(kMax);
    run["occurrences"] = Json::Value(Json::arrayValue);
    run["occurrences_total"] = Json::UInt64(kMax);
    run["chains_ended"] = Json::UInt64(kMax);
    run["lower_bound"] = true;
    Json::Value seed;
    seed["label"] = Json::UInt(std::numeric_limits<uint32_t>::max());
    seed["occurrences"] = Json::Value(Json::arrayValue);
    seed["occurrences_total"] = Json::UInt64(kMax);
    Json::Value block;
    block["kind"] = "column";
    block["k"] = Json::UInt64(kMax);
    block["max_occurrences"] = "unlimited";
    block["complete"] = false;
    block["runs_lower_bound"] = Json::UInt64(kMax);
    block["seed"] = Json::Value(Json::arrayValue);
    block["arms"]["left"] = Json::Value(Json::arrayValue);
    block["arms"]["right"] = Json::Value(Json::arrayValue);
    // the limitation of a real cut list, then at its widest
    const std::string S = random_seq(30, 2700);
    auto anno = repeat_coordinate_index(S, 0);
    LabelOracle oracle(*anno);
    Strategy st = strategy_of(kTrace);
    st.coordinates = true;
    st.max_coordinate_occurrences = 1;
    st.delivery = delivery_costs("graphlet", st.sequences, kMgtFloatWidth, CoordinatesOutput::BLOCK);
    const Seed s = seed_of(S, { "F", "H" });
    const SeedResult r = traverse_seed(oracle, s, st, LabelChangeCost::forbid());
    Json::Value rj = seed_result_to_json(r, st, "graphlet", false);
    Json::Value *lim = nullptr;
    for (Json::Value &l : rj["limitations"]) {
        if (l["kind"].asString() == "coordinates")
            lim = &l;
    }
    ASSERT_TRUE(lim);
    const std::string effect = (*lim)["effect"].asString();
    (*lim)["effect"] = std::to_string(kMax) + effect.substr(effect.find(' '));
    (*lim)["limit"] = "unlimited";
    (*lim)["observed"] = Json::UInt64(kMax);
    (*lim)["lists_cut"] = Json::UInt64(kMax);
    const Json::Value widest = *lim;
    const std::string body = graphlet_text(r, s, st, context_of(oracle), rj);
    const size_t at = body.find("\nK * coordinates ");
    ASSERT_NE(std::string::npos, at);
    const std::string k_line = body.substr(at + 1, body.find('\n', at + 1) - at);
    EXPECT_NE(std::string::npos, k_line.find(std::to_string(kMax)));
    const uint64_t drop_json = 2 + std::strlen("drop_coordinates") + 1;
    const uint64_t drop_mgt = 1 + std::strlen("drop_coordinates");
    // the counts the K record and the token can widen by a digit: Z, graphlet_lines,
    // graphlet_bytes
    const uint64_t digits = 3;
    const uint64_t block_json = member("coordinates", block) + element(widest) + drop_json;
    const uint64_t block_graphlet = block_json + json_escaped_size(k_line) + drop_mgt + digits;
    const uint64_t null_form = member("coordinates", Json::Value())
                             + member("coordinates_reason", "index has no coordinates")
                             + drop_json;
    for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
        const bool graphlet = std::string(detail) == "graphlet";
        const DeliveryCosts d = delivery_costs(detail, true, kMgtFloatWidth, CoordinatesOutput::BLOCK);
        EXPECT_GE(2 * sizeof(Coord) + d.occurrence, kCoordinateAccountPerTextByte * element(occurrence))
            << detail;
        EXPECT_GE(2 * sizeof(RunCoordinates) + d.coordinate_run, kCoordinateAccountPerTextByte * element(run))
            << detail;
        EXPECT_GE(2 * sizeof(SeedCoordinates) + d.coordinate_seed, kCoordinateAccountPerTextByte * element(seed))
            << detail;
        EXPECT_GE(d.coordinate_fixed,
                  kCoordinateAccountPerTextByte * (graphlet ? block_graphlet : block_json)) << detail;
        const DeliveryCosts n = delivery_costs(detail, true, kMgtFloatWidth, CoordinatesOutput::REASON);
        EXPECT_GE(n.coordinate_fixed,
                  kCoordinateAccountPerTextByte * (null_form + (graphlet ? drop_mgt + digits : 0)))
            << detail;
        // without coordinates: no share, and every other price as before
        const DeliveryCosts none = delivery_costs(detail, true);
        EXPECT_EQ(0u, none.coordinate_fixed + none.coordinate_run + none.coordinate_seed + none.occurrence);
        EXPECT_EQ(none.fixed + d.coordinate_fixed, d.fixed) << detail;
        EXPECT_EQ(none.fixed + n.coordinate_fixed, n.fixed) << detail;
        if (std::string(detail) == "full") {
            std::cerr << "full: occurrence " << 2 * sizeof(Coord) + d.occurrence << " / "
                      << element(occurrence) << ", run " << 2 * sizeof(RunCoordinates) + d.coordinate_run
                      << " / " << element(run) << ", seed " << 2 * sizeof(SeedCoordinates) + d.coordinate_seed
                      << " / " << element(seed) << ", block " << d.coordinate_fixed << " / " << block_json
                      << ", null " << n.coordinate_fixed << " / " << null_form << std::endl;
        } else if (graphlet) {
            std::cerr << "graphlet: block " << d.coordinate_fixed << " / " << block_graphlet
                      << ", null " << n.coordinate_fixed << " / " << null_form + drop_mgt + digits
                      << std::endl;
        }
    }
}

// Plan revision 3, end to end (the server's path: each seed's text written once built, its
// account and its coordinate share passed to the attempt): an attempt asking for record
// coordinates measures the very account per text byte the same attempt without them measures
// — their share of the account and the exact text they wrote left out — on seeds of more than
// measured_text_bytes, at every cap, so the server's measured_account_per_text_byte, and the
// walk-until of any later attempt, does not depend on whether one with coordinates came first
TEST(GraphletCoordinates, AttemptsMeasureTheSameRatioWithCoordinates) {
    // one seed S, then 800 distinct tails of 1,300 bp under one column: a fan of 800 paths,
    // about 1 MB of bases
    const std::string S = random_seq(40, 3100);
    std::vector<std::string> seqs, labels;
    std::vector<uint64_t> starts;
    uint64_t at = 0;
    for (size_t i = 0; i < 800; ++i) {
        seqs.push_back(S + random_seq(1300, 3200 + i));
        labels.push_back("F");
        starts.push_back(at);
        at += seqs.back().size() - 31 + 1;
    }
    auto anno = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            31, seqs, labels, DeBruijnGraph::BASIC, true, starts);
    size_t measured = 0;
    // a cut list (S's 800 occurrences at cap 1: the limitation, and its K record in a graphlet)
    // and none
    for (const auto &[detail, cap] : { std::make_pair("full", Json::Value(1)),
                                       std::make_pair("graphlet", Json::Value(1)),
                                       std::make_pair("full", Json::Value("unlimited")) }) {
        {
            const Json::Value plain = request_of({ seed_of(S, { "F" }) },
                    R"({"support": "trace", "direction": "right",
                        "branching": {"on_reconverge": "keep", "max_label_branches": "unlimited",
                                      "max_splits_per_path": "unlimited"},
                        "bounds": {"max_extension_bp": 1400, "max_live_paths": 100000,
                                   "max_paths": 100000, "max_steps": 100000000,
                                   "max_output_bp": 100000000}})", detail);
            double ratio[2] = { 0, 0 };
            uint64_t text[2] = { 0, 0 };
            for (bool coordinates : { false, true }) {
                AttemptSettings settings;
                AttemptRegistry registry(settings);
                auto attempt = std::make_shared<Attempt>(0, std::chrono::system_clock::time_point(),
                                                         registry.settings(),
                                                         registry.server_instance(), true);
                ResultTexts texts;
                process_traverse_request(coordinates ? with_coordinates(plain, cap) : plain,
                                         *anno, "", {}, nullptr, attempt.get(), &texts);
                ASSERT_EQ(1u, texts.texts.size());
                text[coordinates] = texts.texts[0].size();
                ratio[coordinates] = attempt->own_account_per_text_byte();
                // as the server does once the response is written
                registry.note_account_per_text_byte(detail, ratio[coordinates]);
                if (coordinates) {
                    EXPECT_NE(std::string::npos, texts.texts[0].find("\"coordinates\":{"));
                } else {
                    EXPECT_EQ(std::string::npos, texts.texts[0].find("\"coordinates\""));
                }
            }
            const std::string what = std::string(detail) + " cap " + compact_json(cap);
            EXPECT_GT(text[false], kMeasuredTextBytes) << what;
            EXPECT_GT(text[true], text[false]) << what;
            EXPECT_GT(ratio[false], 0.0) << what;
            EXPECT_EQ(ratio[false], ratio[true]) << what;
            measured += ratio[true] > 0;
        }
    }
    EXPECT_EQ(3u, measured);
}

// The probe's coordinates block (feature level 6, plan revisions 7 and 8) follows the index it
// describes: supported exactly where trace is (supports_trace: coordinates on a basic graph), the
// kinds that index can report (record and mixed only with a CoordToHeader), the cap's default
// the request's; and it is in GET /traverse/capabilities only, never in the per-request
// capabilities, whose level-6 difference is the digit alone
TEST(GraphletCoordinates, CapabilitiesBlockFollowsTheIndex) {
    const std::string S = random_seq(40, 3300);
    const std::vector<std::string> seqs { S + random_seq(20, 3301), S + random_seq(20, 3302) };
    auto plain = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, seqs, { "F", "H" });
    auto coords = test::build_anno_graph<DBGSuccinct, annot::ColumnCompressed<>>(
            11, seqs, { "F", "H" }, DeBruijnGraph::BASIC, true, { 0, 0 });
    const HeaderIndexCase headers = header_index(11, { { "r1", seqs[0] }, { "r2", seqs[1] } });
    const LabelOracle without(*plain), with(*coords), with_headers(*headers.anno, headers.cth.get());
    struct Case {
        const char *name;
        const LabelOracle *oracle;
        bool supported;
        std::vector<std::string> kinds;
    };
    const Case cases[] = {
        { "no coordinates", &without, false, {} },
        { "coordinates", &with, true, { "column" } },
        { "coordinates and headers", &with_headers, true, { "record", "column", "mixed" } },
    };
    for (const Case &c : cases) {
        const Json::Value block = coordinates_capabilities_json(*c.oracle);
        const Json::Value per_request = capabilities_to_json(*c.oracle, "", nullptr);
        EXPECT_EQ(c.supported, block["supported"].asBool()) << c.name;
        EXPECT_EQ(per_request["supports_trace"].asBool(), block["supported"].asBool()) << c.name;
        std::vector<std::string> kinds;
        for (const Json::Value &k : block["kinds"]) {
            kinds.push_back(k.asString());
        }
        EXPECT_EQ(c.kinds, kinds) << c.name;
        EXPECT_FALSE(per_request.isMember("coordinates")) << c.name;
        EXPECT_EQ(6, per_request["feature_level"].asInt()) << c.name;
        // the same knobs, default and texts whatever the index: what a client reads before
        // asking
        EXPECT_EQ("output.coordinates", block["knob"].asString());
        EXPECT_EQ("output.max_coordinate_occurrences", block["cap_knob"].asString());
        ASSERT_TRUE(block["max_occurrences_default"].isUInt64());
        EXPECT_EQ(Strategy().max_coordinate_occurrences, block["max_occurrences_default"].asUInt64());
        EXPECT_EQ(16u, block["max_occurrences_default"].asUInt64());
        EXPECT_EQ("coordinates", block["limitation"].asString());
        EXPECT_EQ("drop_coordinates", block["action"].asString());
        // the column record-end numbering (review of W1, finding 6) and the true bound, which
        // names the probe's maxima rather than a value it could contradict
        EXPECT_NE(std::string::npos, block["rule"].asString().find(
                "a column label's interval in a record's last k - 1 bases shares its numbers with "
                "the next record's first k - 1 positions")) << c.name;
        EXPECT_NE(std::string::npos, block["output_bound"].asString().find(
                "this server's max_memory_mb")) << c.name;
    }
}

// D3: a seed delivered with the set derived from part of it states a `derivation` limitation the
// fixed part does not hold; the walker charges it (DeliveryCosts::extra_limitation) and the model
// still bounds what the serialisers hold, in every detail, with and without a memory budget
TEST(Graphlet, PartialDerivationIsPricedAndBound) {
    size_t checked = 0;
    for (BudgetCase c : budget_cases()) {
        // the seed's labels derived (forbid cost, no extra labels: a derived set is then valid)
        if (c.st.label_mode == LabelMode::ANNOTATE || c.cost.finite() || !c.st.extra.empty())
            continue;
        c.seed.labels.clear();
        LabelOracle oracle(*c.anno);
        // at least two k-mers: after the last one the derivation is complete, not partial
        if (c.seed.sequence.size() < oracle.get_k() + 1)
            continue;
        for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
            for (uint64_t mb : { uint64_t(0), uint64_t(64) }) {
                const std::string what = c.name + " " + detail + " " + std::to_string(mb) + " MiB";
                Strategy st = c.st;
                st.time_budget_ms = 1e-9;            // spent after the first k-mer (read whole)
                st.max_memory_bytes = mb << 20;
                st.delivery = delivery_costs(detail, st.sequences);
                const SeedResult r = run_case(c, st);
                ASSERT_TRUE(r.derivation_partial) << what;
                const uint64_t model = modelled_delivery(r, st.delivery);
                Json::Value j = seed_result_to_json(r, st, detail, false);
                ASSERT_EQ("derivation", j["limitations"][0]["kind"].asString()) << what;
                EXPECT_EQ("qualified", j["outcome"]["label_evidence"].asString()) << what;
                uint64_t before = 0;
                if (std::string(detail) == "graphlet") {
                    const std::string text = graphlet_text(r, c.seed, st, context_of(oracle), j);
                    EXPECT_NE(std::string::npos, text.find("\nK * derivation bounds.time_budget_ms "))
                        << what;
                    // the record's extra fields are a failed derivation's (cause, server_limit):
                    // n is the S record's num_kmers (review of W1, finding 5)
                    EXPECT_NE(std::string::npos, text.find(" * cause=s:time_budget ")) << what;
                    EXPECT_EQ(std::string::npos, text.find("num_kmers=")) << what;
                    EXPECT_NE(std::string::npos, text.find("\nO p c q i\n")) << what;
                    before = 3 * text.size() + json_tree_bytes(j);
                    j["graphlet"] = text;
                }
                check_widths(j, 2, "", what);
                // what the serialisers hold is within the model and within what the walker
                // charged. (This test's model counts a continuation label as a whole end label
                // more, which delivery_costs prices within the end label's own (its member for
                // end_reasons): at a walk censored at depth 0 no reservation's slack covers that
                // overcount, so the model is compared with the output, not with the charge)
                const uint64_t actual = delivery_peak(j, before);
                EXPECT_LE(actual, model) << what;
                EXPECT_LE(actual, r.account.memory_final) << what;
                checked++;
            }
        }
    }
    EXPECT_GT(checked, 0u);
}
