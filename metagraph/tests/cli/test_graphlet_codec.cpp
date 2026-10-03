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
#include "cli/traverse.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
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
            if (s.partition != "*")
                EXPECT_FALSE(annotate || s.parents.size() < 2) << "G.partition explicit where * holds";
            s.total = s.entry_total == "*" ? s.entry_set.size() : u(s.entry_total);
            if (s.entry_total != "*") EXPECT_NE(s.total, s.entry_set.size());
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
            if (!sequences_absent) EXPECT_EQ(s.length, bases_of(s).size());
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
    EXPECT_NE(fp, index_manifest_fingerprint(manifest("m2.json", "annotation B"), {}));
    // a stated index_fp must be the computed one
    EXPECT_EQ(fp, index_manifest_fingerprint(manifest("m3.json", "annotation A", fp), { graph, anno }));
    EXPECT_THROW(index_manifest_fingerprint(manifest("m4.json", "annotation A", std::string(64, '0')), {}),
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

void count_json(const Json::Value &v, uint64_t *members, uint64_t *elements, uint64_t *strings) {
    if (v.isObject()) {
        for (const std::string &name : v.getMemberNames()) {
            (*members)++;
            *strings += name.size();
            count_json(v[name], members, elements, strings);
        }
    } else if (v.isArray()) {
        for (const Json::Value &x : v) {
            (*elements)++;
            count_json(x, members, elements, strings);
        }
    } else if (v.isString()) {
        *strings += v.asString().size();
    }
}

// what the delivery model charges for |r| in |detail| (DeliveryCosts per object), from
// the result's own objects: a lower bound of the walker's charge, which prices lists at
// their bounds
uint64_t modelled_delivery(const SeedResult &r, const DeliveryCosts &d) {
    uint64_t m = d.fixed;
    for (const LabelRef &l : r.label_dict) m += d.label + l.name.size() * d.label_name;
    for (const ArmResult &arm : r.arms) {
        for (const Segment &s : arm.segments) {
            m += d.segment + s.length_bp * d.base;
            size_t labels = s.labels_start.size() + s.labels_end.size() + s.parents.size();
            for (const auto &via : s.labels_via_parent) labels += via.size() + 1;
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
            for (const auto &rf : be.refused) m += (1 + rf.labels.size()) * d.branch_event_entry;
        }
        m += arm.growth.size() * d.bin;
    }
    return m;
}

} // namespace

// The delivery model bounds what the serialisers hold (§14: every committed prefix stays
// deliverable within what was reserved for it): per detail, the MGT body three times
// (text, its JSON copy, the response), or the JSON tree (jsoncpp's map nodes, keys and
// strings) with its indented text — and stays within a small factor of it, or the
// budget would refuse far too early.
TEST(Graphlet, DeliveryCostsBoundTheOutput) {
    size_t checked = 0;
    for (const BudgetCase &c : budget_cases()) {
        LabelOracle oracle(*c.anno);
        const SeedResult r = run_case(c, c.st);
        for (const char *detail : { "summary", "tree", "full", "graphlet" }) {
            const DeliveryCosts d = delivery_costs(detail, c.st.sequences);
            const uint64_t model = modelled_delivery(r, d);
            uint64_t actual = 0;
            const Json::Value j = seed_result_to_json(r, c.st, detail, false);
            if (std::string(detail) == "graphlet") {
                actual = 3 * graphlet_text(r, c.seed, c.st, context_of(oracle), j).size();
            } else {
                uint64_t members = 0, elements = 0, strings = 0;
                count_json(j, &members, &elements, &strings);
                actual = 128 * members + 96 * elements + 2 * strings + j.toStyledString().size();
            }
            const std::string what = c.name + " " + detail;
            EXPECT_LE(actual, model) << what;
            // An MGT body writes derivable fields as '*' and sets as deltas, which a bound
            // fixed per object cannot know before the walk, so its model is looser (and
            // matters less: the walker's own state is several times the body)
            const uint64_t factor = std::string(detail) == "graphlet" ? 10 : 6;
            EXPECT_LE(model - d.fixed, factor * actual) << what << ": the model is not useful";
            checked++;
        }
    }
    EXPECT_GT(checked, 40u);
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
    // the lever works: the default cap records 64 labels per node, which fit
    const Json::Value ok = request(64);
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

    // a long seed, explicit and derived: cut within about one interval
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
            // checked once the interval passed, within one more chunk or row of it
            EXPECT_LE(e.stop().used, 2.0 * kWorkCheckInterval) << what;
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
