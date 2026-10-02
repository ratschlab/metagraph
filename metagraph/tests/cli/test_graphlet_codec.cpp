#include "gtest/gtest.h"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <unistd.h>

#include <json/json.h>

#include "tests/annotation/test_annotated_dbg_helpers.hpp"
#include "cli/traverse.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"


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


// Label names come from FASTA headers and file names, which need not be UTF-8. They are
// made valid once (every maximal ill-formed subsequence -> one U+FFFD), so that detail:
// full, the summary and the MGT body carry the same bytes. The expectations are Python's
// bytes.decode('utf-8', 'replace'), the reader's side of the same rule.
TEST(GraphletCodec, NamesAreMadeValidUtf8) {
    const std::string R = "\xEF\xBF\xBD";   // U+FFFD
    const std::vector<std::pair<std::string, std::string>> cases = {
        { "abc", "abc" },
        { "caf\xC3\xA9", "caf\xC3\xA9" },
        { "\xF0\x9F\x98\x80" "z", "\xF0\x9F\x98\x80" "z" },
        { "\xC0" "abc", R + "abc" },                 // C0: never a lead byte
        { "\xE9" "ab", R + "ab" },                   // Latin-1: a lead without its tail
        { "\xED\xA0\x80", R + R + R },              // a surrogate
        { "\xE0\x80\x80", R + R + R },              // overlong
        { "\xF4\x90\x80\x80", R + R + R + R },      // beyond U+10FFFF
        { "\xF0\x9F\x98", R },                      // truncated: one maximal subpart
        { "\xE2\x82" "x", R + "x" },
        { "\xFF", R },
        { "\x80\x80", R + R },
    };
    for (const auto &[raw, want] : cases) {
        EXPECT_EQ(want, to_valid_utf8(raw)) << raw;
        // front coding then works on the bytes that are transported
        EXPECT_EQ(want, front_decode("", front_code("", want)));
    }
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

const char kCodes[] = "DLBRUVJTXSPNOMW";

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
            auto f = fields(line, line[0] == 'L' || line[0] == 'X' ? 4 : line[0] == 'K' ? 9 : SIZE_MAX);
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
    for (const char *kind : { "H", "S", "L", "O", "K", "A", "B", "V", "G", "P", "E", "T", "C", "R", "Z",
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
