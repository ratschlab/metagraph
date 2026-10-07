#include <algorithm>
#include <chrono>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <json/json.h>
#include "gtest/gtest.h"

#include "../annotation/test_annotated_dbg_helpers.hpp"

#include "annotation/coord_to_header.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "cli/pattern.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


// output.labels "all" of POST /pattern (src/cli/pattern_retrieval.cpp, increment 3 of
// docs/DESIGN-pattern-search.md): the labels of every released context (LabelRecorder
// discovery), their placement (LabelQuery with coordinates, the record mapping first and the
// offset after), the deduplication of placed occurrences, the modes, and every refusal,
// truncation, cut and stop stated. Tiny graphs built from explicit records in labelled columns
// with their coordinates and a record mapping; the expectations come from a RECORD-SCAN
// ORACLE over those records (every occurrence of the pattern and of its reverse complement,
// 1-based, per (column, seq_id, strand)), never from the engine.

namespace {

using namespace mtg;
using namespace mtg::graph;
using namespace mtg::cli;
using Clock = pattern::Deadline::Clock;

struct Record {
    std::string column;
    std::string header;
    std::string seq;
};

std::string rc(std::string s) {
    reverse_complement(s.begin(), s.end());
    return s;
}

/**
 * An index of |records| (each a column label and a record of it, in seq_id order per column):
 * a BASIC (or |mode|) succinct graph at |k| with the dummy-edge mask, an annotation of type
 * |Annotation| (with coordinates when |coordinates|: each record's k-mers numbered from its
 * column's running count, as `annotate --coordinates` numbers them), and the record mapping
 * (the .seqs of `annotate --index-header-coords`) built from the same records.
 */
struct Index {
    size_t k = 0;
    std::vector<Record> records;
    std::unique_ptr<AnnotatedDBG> anno;
    std::unique_ptr<annot::CoordToHeader> cth;
    // the first column coordinate of each record
    std::vector<uint64_t> starts;
};

template <class Annotation>
Index build(size_t k, const std::vector<Record> &records, bool coordinates,
            DeBruijnGraph::Mode mode = DeBruijnGraph::BASIC) {
    Index idx;
    idx.k = k;
    idx.records = records;
    std::vector<std::string> seqs, labels;
    std::map<std::string, uint64_t> next;
    for (const Record &r : records) {
        EXPECT_GE(r.seq.size(), k);
        seqs.push_back(r.seq);
        labels.push_back(r.column);
        idx.starts.push_back(next[r.column]);
        next[r.column] += r.seq.size() - k + 1;
    }
    idx.anno = test::build_anno_graph<DBGSuccinct, Annotation>(
            k, seqs, labels, mode, coordinates,
            coordinates ? idx.starts : std::vector<uint64_t>{});
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

// k = 7. The pattern AC (and GT, its reverse complement) in:
//  - c1/r0: AC at 0-based 2 and 5, both inside the k-mer GACGACT (offsets 1 and 4): two
//    offsets in one k-mer, two occurrences;
//  - c1/r1: GT at 4 only (the - strand);
//  - c2/r0 and c2/r1: equal lengths, AC at local 4 in both: told apart by seq_id;
//  - c3/r0: the k-mer ACGTACG three times in one record (locals 0, 4, 8): one context, three
//    occurrences;
//  - c4/r0: c1/r0 again, so that its rows carry two labels (the per-anchor cap).
const size_t kK = 7;
const std::vector<Record> kRecords = {
    { "c1", "c1_r0", "GGACGACTTTG" },
    { "c1", "c1_r1", "CCCTGTCCCCA" },
    { "c2", "c2_r0", "TTTTACTTTTG" },
    { "c2", "c2_r1", "GGGGACGGGGC" },
    { "c3", "c3_r0", "ACGTACGTACGTACG" },
    { "c4", "c4_r0", "GGACGACTTTG" },
};

PatternLimits limits() {
    PatternLimits l;
    // the 2-mers of the tiny graphs are below any useful floor
    l.min_information_bits = 0;
    return l;
}

Json::Value run(const Index &idx, const std::string &body, RetrievalHooks hooks = {},
                bool records = true,
                const std::function<Clock::time_point()> &clock = nullptr) {
    if (records)
        hooks.coord_to_header = idx.cth.get();
    return process_pattern_request(parse_pattern_body(body), *idx.anno, limits(), "rel",
                                   nullptr, nullptr, clock, &hooks);
}

std::pair<int, std::string> refusal(const Index &idx, const std::string &body) {
    try {
        run(idx, body);
    } catch (const PatternRefusal &e) {
        return { e.status(), e.code() };
    }
    return { 200, "" };
}

std::vector<std::pair<std::string, std::string>> oriented(const std::string &p) {
    if (rc(p) == p)
        return { { "=", p } };
    return { { "+", p }, { "-", rc(p) } };
}

using Context = std::tuple<std::string, std::string, uint64_t>;        // strand, k-mer, offset
using Placed = std::tuple<std::string, uint64_t, uint64_t, std::string>; // column, seq_id, start, strand

// the record-scan oracle: every occurrence of the oriented pattern in every record,
// (column, seq_id, 1-based start, strand)
std::set<Placed> placed_oracle(const Index &idx, const std::string &p) {
    std::set<Placed> out;
    std::map<std::string, uint64_t> seq_id;
    for (const Record &r : idx.records) {
        const uint64_t id = seq_id[r.column]++;
        for (const auto &[strand, q] : oriented(p)) {
            for (size_t i = 0; i + q.size() <= r.seq.size(); ++i) {
                if (r.seq.compare(i, q.size(), q) == 0)
                    out.emplace(r.column, id, i + 1, strand);
            }
        }
    }
    return out;
}

// the graph contexts of |p| carried by each column (BASIC: the records' k-mers as deposited)
std::map<std::string, std::set<Context>> contexts_by_column(const Index &idx,
                                                            const std::string &p) {
    std::map<std::string, std::set<Context>> out;
    for (const Record &r : idx.records) {
        for (size_t s = 0; s + idx.k <= r.seq.size(); ++s) {
            const std::string kmer = r.seq.substr(s, idx.k);
            for (const auto &[strand, q] : oriented(p)) {
                for (size_t o = 0; o + q.size() <= idx.k; ++o) {
                    if (kmer.compare(o, q.size(), q) == 0)
                        out[r.column].emplace(strand, kmer, o);
                }
            }
        }
    }
    return out;
}

const Record& record_of(const Index &idx, const std::string &column, uint64_t seq_id) {
    uint64_t i = 0;
    for (const Record &r : idx.records) {
        if (r.column == column && i++ == seq_id)
            return r;
    }
    throw std::out_of_range("no record " + column + "/" + std::to_string(seq_id));
}

std::pair<uint64_t, uint64_t> nt_coords(const std::string &s) {
    const size_t dash = s.find('-');
    return { std::stoull(s.substr(0, dash)), std::stoull(s.substr(dash + 1)) };
}

/**
 * A complete answer with record placement against both oracles: every context with exactly
 * the columns carrying its k-mer, every placed occurrence the record-scan oracle's, inside
 * the context's own k-mer of its record (the mapping first, the offset after), the per-label
 * summary, and every count exact.
 */
void check_complete_record_answer(const Index &idx, const Json::Value &e, const std::string &p) {
    const size_t L = p.size();
    const auto by_column = contexts_by_column(idx, p);
    std::set<Context> all;
    for (const auto &[column, set] : by_column) {
        all.insert(set.begin(), set.end());
    }
    EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
    EXPECT_TRUE(e["withheld"].isNull());
    EXPECT_TRUE(e["stop"].isNull());
    EXPECT_EQ("record", e["placement"].asString());
    EXPECT_EQ(0u, e["rows_refused"].size());
    EXPECT_EQ(0u, e["anchors_truncated"].size());
    ASSERT_EQ(all.size(), e["results"].size());

    std::map<std::string, std::set<std::tuple<uint64_t, uint64_t, std::string>>> unions;
    std::set<Context> got;
    for (const Json::Value &r : e["results"]) {
        const Context c { r["strand"].asString(), r["kmer"].asString(), r["offset"].asUInt64() };
        got.insert(c);
        std::set<std::string> expected_columns;
        for (const auto &[column, set] : by_column) {
            if (set.count(c))
                expected_columns.insert(column);
        }
        EXPECT_EQ("kmer", r["support"].asString());
        EXPECT_EQ("complete", r["labels_status"].asString());
        EXPECT_EQ(expected_columns.size(), r["labels_total"].asUInt64());
        std::set<std::string> columns;
        for (const Json::Value &l : r["labels"]) {
            const std::string column = l["column"].asString();
            columns.insert(column);
            EXPECT_EQ("kmer", l["support"].asString());
            const Json::Value &list = l["occurrence_list"];
            EXPECT_EQ(list.size(), l["occurrences"]["value"].asUInt64());
            EXPECT_EQ("exact", l["occurrences"]["relation"].asString());
            EXPECT_EQ("placed_occurrences", l["occurrences"]["unit"].asString());
            EXPECT_GT(list.size(), 0u) << r;
            for (const Json::Value &o : list) {
                const uint64_t seq_id = o["seq_id"].asUInt64();
                const Record &rec = record_of(idx, column, seq_id);
                EXPECT_EQ(rec.header, o["record"].asString());
                EXPECT_EQ(rec.seq.size(), o["nt_length"].asUInt64());
                EXPECT_EQ(std::get<0>(c), o["strand"].asString());
                const auto [start, end] = nt_coords(o["nt_coords"].asString());
                EXPECT_EQ(start + L - 1, end);
                // the context's own k-mer sits where the occurrence lies, at its offset
                ASSERT_GE(start - 1, std::get<2>(c));
                EXPECT_EQ(std::get<1>(c), rec.seq.substr(start - 1 - std::get<2>(c), idx.k))
                        << column << "/" << seq_id << " " << o;
                unions[column].emplace(seq_id, start, o["strand"].asString());
            }
        }
        EXPECT_EQ(expected_columns, columns) << r;
    }
    EXPECT_EQ(all, got);

    // the occurrences: the record scan's, each once (the contexts of one occurrence collapse)
    std::set<Placed> placed;
    uint64_t total = 0;
    for (const auto &[column, set] : unions) {
        for (const auto &[seq_id, start, strand] : set) {
            placed.emplace(column, seq_id, start, strand);
        }
        total += set.size();
    }
    EXPECT_EQ(placed_oracle(idx, p), placed);

    // by_label: (contexts desc, column asc), counts exact
    const Json::Value &by_label = e["by_label"];
    ASSERT_EQ(by_column.size(), by_label.size());
    std::vector<std::pair<int64_t, std::string>> order;
    for (const Json::Value &b : by_label) {
        const std::string column = b["column"].asString();
        order.emplace_back(-b["contexts"]["value"].asInt64(), column);
        EXPECT_EQ(by_column.at(column).size(), b["contexts"]["value"].asUInt64());
        EXPECT_EQ("exact", b["contexts"]["relation"].asString());
        uint64_t suffix = 0;
        for (const Context &c : by_column.at(column)) {
            suffix += std::get<2>(c) == idx.k - L;
        }
        EXPECT_EQ(suffix, b["contexts_suffix"]["value"].asUInt64());
        EXPECT_EQ(unions[column].size(), b["occurrences"]["value"].asUInt64()) << column;
        EXPECT_EQ("exact", b["occurrences"]["relation"].asString());
        EXPECT_TRUE(b["graph"].isNull());
    }
    EXPECT_TRUE(std::is_sorted(order.begin(), order.end()));
    EXPECT_EQ(by_column.size(), e["counts"]["labels"]["value"].asUInt64());
    EXPECT_EQ("exact", e["counts"]["labels"]["relation"].asString());
    EXPECT_EQ("labels", e["counts"]["labels"]["unit"].asString());
    EXPECT_EQ(total, e["counts"]["occurrences"]["value"].asUInt64());
    EXPECT_EQ("exact", e["counts"]["occurrences"]["relation"].asString());
}

std::string body(const std::string &patterns, const std::string &rest = "") {
    return "{\"patterns\": [" + patterns + "], \"output\": {\"labels\": \"all\"}"
           + (rest.empty() ? "" : ", " + rest) + "}";
}


TEST(PatternRetrieval, RecordPlacementMatchesTheRecordScan) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}, {\"dna\": \"ACGTA\"}"));
    EXPECT_EQ("all", out["output"]["labels"].asString());
    EXPECT_TRUE(out["output"]["occurrences"].asBool());
    const Json::Value &l = out["limits"];
    EXPECT_EQ(64u, l["max_labels_per_anchor"].asUInt64());
    EXPECT_EQ(256u, l["max_memory_mb"].asUInt64());
    EXPECT_EQ(1000u, l["max_labels"].asUInt64());
    EXPECT_EQ(16u, l["max_occurrences_per_label"].asUInt64());
    EXPECT_EQ(100000000u, l["max_annotation_work"].asUInt64());
    EXPECT_FALSE(l["allow_unbudgeted_annotation"].asBool());
    ASSERT_EQ(2u, out["patterns"].size());
    for (size_t i = 0; i < 2; ++i) {
        const Json::Value &e = out["patterns"][static_cast<int>(i)];
        EXPECT_EQ("budgeted", e["annotation"].asString());
        EXPECT_EQ(0u, e["notes"].size()) << e["notes"];
        EXPECT_GT(e["work"]["annotation_rows"].asUInt64(), 0u);
        EXPECT_GT(e["work"]["annotation_units"].asUInt64(), 0u);
        EXPECT_GT(e["work"]["memory_bytes"].asUInt64(), 0u);
        EXPECT_TRUE(e["labels_cut"].isNull());
        EXPECT_TRUE(e["occurrences_cut"].isNull());
        check_complete_record_answer(idx, e, i ? "ACGTA" : "AC");
    }
}

TEST(PatternRetrieval, TwoOffsetsInOneKmerAreTwoOccurrences) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"strands\": \"forward\""));
    const Json::Value &e = out["patterns"][0];
    // GACGACT: c1/r0 local 1 (and c4/r0), AC at offsets 1 and 4 -> starts 3 and 6
    std::map<uint64_t, std::set<std::pair<uint64_t, std::string>>> by_offset;
    for (const Json::Value &r : e["results"]) {
        if (r["kmer"].asString() != "GACGACT")
            continue;
        for (const Json::Value &label : r["labels"]) {
            if (label["column"].asString() != "c1")
                continue;
            for (const Json::Value &o : label["occurrence_list"]) {
                by_offset[r["offset"].asUInt64()].emplace(o["seq_id"].asUInt64(),
                                                          o["nt_coords"].asString());
            }
        }
    }
    const std::map<uint64_t, std::set<std::pair<uint64_t, std::string>>> expected {
        { 1, { { 0, "3-4" } } }, { 4, { { 0, "6-7" } } },
    };
    EXPECT_EQ(expected, by_offset);
}

TEST(PatternRetrieval, RepeatedKmerInOneRecord) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"strands\": \"forward\""));
    // ACGTACG at locals 0, 4, 8 of c3/r0: one context at offset 0, three occurrences
    bool seen = false;
    for (const Json::Value &r : out["patterns"][0]["results"]) {
        if (r["kmer"].asString() != "ACGTACG" || r["offset"].asUInt64() != 0)
            continue;
        seen = true;
        ASSERT_EQ(1u, r["labels"].size());
        const Json::Value &label = r["labels"][0];
        EXPECT_EQ("c3", label["column"].asString());
        std::vector<std::string> coords;
        for (const Json::Value &o : label["occurrence_list"]) {
            EXPECT_EQ(0u, o["seq_id"].asUInt64());
            EXPECT_EQ("c3_r0", o["record"].asString());
            coords.push_back(o["nt_coords"].asString());
        }
        EXPECT_EQ(std::vector<std::string>({ "1-2", "5-6", "9-10" }), coords);
        EXPECT_EQ(3u, label["occurrences"]["value"].asUInt64());
    }
    EXPECT_TRUE(seen);
}

TEST(PatternRetrieval, EqualLengthRecordsToldApartBySeqId) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"strands\": \"forward\""));
    std::set<std::tuple<uint64_t, std::string, std::string>> c2;
    for (const Json::Value &r : out["patterns"][0]["results"]) {
        for (const Json::Value &label : r["labels"]) {
            if (label["column"].asString() != "c2")
                continue;
            for (const Json::Value &o : label["occurrence_list"]) {
                c2.emplace(o["seq_id"].asUInt64(), o["record"].asString(),
                           o["nt_coords"].asString());
            }
        }
    }
    // the same local coordinates in two records of one column: two occurrences
    const std::set<std::tuple<uint64_t, std::string, std::string>> expected {
        { 0, "c2_r0", "5-6" }, { 1, "c2_r1", "5-6" },
    };
    EXPECT_EQ(expected, c2);
    for (const Json::Value &b : out["patterns"][0]["by_label"]) {
        if (b["column"].asString() == "c2") {
            EXPECT_EQ(2u, b["occurrences"]["value"].asUInt64());
        }
    }
}

TEST(PatternRetrieval, OffsetAddedAfterTheRecordMapping) {
    // k = 5: the column's records ACGTA, CCCCC, CCCCC, CCCCC hold one k-mer each, at column
    // coordinates 0, 1, 2, 3. GTA sits at offset 2 of ACGTA: record 0, local 0, start 3 — not
    // the record of column coordinate 0 + 2
    const std::vector<Record> records = {
        { "m", "m0", "ACGTA" }, { "m", "m1", "CCCCC" }, { "m", "m2", "CCCCC" },
        { "m", "m3", "CCCCC" },
    };
    Index idx = build<annot::RowDiffColumnAnnotator>(5, records, true);
    Json::Value out = run(idx, body("{\"dna\": \"GTA\"}", "\"strands\": \"forward\""));
    const Json::Value &e = out["patterns"][0];
    ASSERT_EQ(1u, e["results"].size()) << e;
    const Json::Value &r = e["results"][0];
    EXPECT_EQ("ACGTA", r["kmer"].asString());
    EXPECT_EQ(2u, r["offset"].asUInt64());
    ASSERT_EQ(1u, r["labels"].size());
    ASSERT_EQ(1u, r["labels"][0]["occurrence_list"].size());
    const Json::Value &o = r["labels"][0]["occurrence_list"][0];
    EXPECT_EQ(0u, o["seq_id"].asUInt64());
    EXPECT_EQ("m0", o["record"].asString());
    EXPECT_EQ("3-5", o["nt_coords"].asString());
    EXPECT_EQ(5u, o["nt_length"].asUInt64());
    EXPECT_EQ("+", o["strand"].asString());
    check_complete_record_answer(idx, out["patterns"][0], "GTA");
}

TEST(PatternRetrieval, RefusedRow) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    RetrievalHooks hooks;
    // every charge of every read refused: no row fits
    hooks.deny_decode = [](uint64_t) { return true; };
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"mode\": \"partial\""), hooks);
    const Json::Value &e = out["patterns"][0];
    EXPECT_FALSE(e["retrieval_complete"].asBool());
    EXPECT_TRUE(e["withheld"].isNull());
    EXPECT_TRUE(e["stop"].isNull());
    // every distinct row stated once, with its k-mer and the phase
    std::set<uint64_t> rows;
    for (const Json::Value &r : e["results"]) {
        EXPECT_EQ("refused", r["labels_status"].asString());
        EXPECT_TRUE(r["labels"].isNull());
        EXPECT_TRUE(r["labels_total"].isNull());
        rows.insert(r["row"].asUInt64());
    }
    ASSERT_EQ(rows.size(), e["rows_refused"].size());
    for (const Json::Value &x : e["rows_refused"]) {
        EXPECT_TRUE(rows.count(x["row"].asUInt64()));
        EXPECT_EQ("label_discovery", x["phase"].asString());
        EXPECT_EQ("max_memory", x["reason"].asString());
        EXPECT_EQ(7u, x["kmer"].asString().size());
    }
    EXPECT_EQ("at_least", e["counts"]["labels"]["relation"].asString());
    EXPECT_EQ(0u, e["by_label"].size());

    // all_or_count: withheld at the first refused row, the row stated
    out = run(idx, body("{\"dna\": \"AC\"}"), hooks);
    const Json::Value &w = out["patterns"][0];
    EXPECT_EQ("annotation_budget", w["withheld"]["reason"].asString());
    EXPECT_EQ(1u, w["rows_refused"].size());
    EXPECT_EQ(0u, w["results"].size());
    EXPECT_EQ(0u, w["returned"].asUInt64());
    EXPECT_FALSE(w["retrieval_complete"].asBool());
    EXPECT_TRUE(w["by_label"].isNull());
    EXPECT_EQ("unknown", w["counts"]["labels"]["relation"].asString());
    EXPECT_EQ("unknown", w["counts"]["occurrences"]["relation"].asString());
    // the graph count is unchanged: exact
    EXPECT_EQ("exact", w["counts"]["contexts"]["relation"].asString());
}

TEST(PatternRetrieval, TruncatedAnchor) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    // the k-mers of c1/r0 are c4/r0's too: two labels each, one kept
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"max_labels_per_anchor\": 1"));
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("anchor_labels_truncated", e["withheld"]["reason"].asString());
    EXPECT_EQ(0u, e["results"].size());
    std::set<std::string> truncated;
    for (const Json::Value &t : e["anchors_truncated"]) {
        EXPECT_EQ(1u, t["cap"].asUInt64());
        EXPECT_EQ(2u, t["total"].asUInt64());
        truncated.insert(t["kmer"].asString());
    }
    // every k-mer of GGACGACTTTG contains AC
    const std::set<std::string> expected { "GGACGAC", "GACGACT", "ACGACTT", "CGACTTT",
                                           "GACTTTG" };
    EXPECT_EQ(expected, truncated);

    // partial: the cut stated per context, the first label (ascending column) kept
    out = run(idx, body("{\"dna\": \"AC\"}", "\"max_labels_per_anchor\": 1, \"mode\": \"partial\""));
    const Json::Value &p = out["patterns"][0];
    EXPECT_TRUE(p["withheld"].isNull());
    EXPECT_FALSE(p["retrieval_complete"].asBool());
    EXPECT_EQ(5u, p["anchors_truncated"].size());
    size_t cut = 0;
    for (const Json::Value &r : p["results"]) {
        if (!expected.count(r["kmer"].asString())) {
            EXPECT_EQ("complete", r["labels_status"].asString());
            continue;
        }
        cut++;
        EXPECT_EQ("truncated", r["labels_status"].asString());
        EXPECT_EQ(2u, r["labels_total"].asUInt64());
        ASSERT_EQ(1u, r["labels"].size());
        EXPECT_EQ("c1", r["labels"][0]["column"].asString());
    }
    EXPECT_GT(cut, 0u);
    EXPECT_EQ("at_least", p["counts"]["labels"]["relation"].asString());
    EXPECT_EQ("at_least", p["counts"]["occurrences"]["relation"].asString());
    for (const Json::Value &b : p["by_label"]) {
        EXPECT_EQ("at_least", b["contexts"]["relation"].asString());
    }
}

TEST(PatternRetrieval, MemoryAccountStopsTheRun) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    // the descriptors of six contexts (the model: 512 + 2k bytes each) and nothing more:
    // partial's descriptors take at most half of the account, three of them
    RetrievalHooks hooks;
    hooks.max_memory_bytes = 6 * (512 + 2 * kK);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"mode\": \"partial\""), hooks);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ(3u, e["returned"].asUInt64());
    EXPECT_EQ(3u, e["results"].size());
    EXPECT_EQ("max_memory", e["cut"]["reason"].asString());
    EXPECT_EQ("output", e["stop"]["phase"].asString());
    EXPECT_EQ("max_memory", e["stop"]["reason"].asString());
    EXPECT_FALSE(e["retrieval_complete"].asBool());
    // the graph count stays what discovery found
    EXPECT_EQ("exact", e["counts"]["contexts"]["relation"].asString());
    EXPECT_GT(e["counts"]["contexts"]["value"].asUInt64(), 3u);

    out = run(idx, body("{\"dna\": \"AC\"}"), hooks);
    EXPECT_EQ("output_budget", out["patterns"][0]["withheld"]["reason"].asString());
    EXPECT_EQ(0u, out["patterns"][0]["results"].size());

    // one byte under what a complete run held at its peak: something is refused or stopped,
    // and said
    Json::Value full = run(idx, body("{\"dna\": \"AC\"}"));
    const uint64_t peak = full["patterns"][0]["work"]["memory_bytes"].asUInt64();
    ASSERT_TRUE(full["patterns"][0]["retrieval_complete"].asBool());
    hooks.max_memory_bytes = peak - 1;
    out = run(idx, body("{\"dna\": \"AC\"}"), hooks);
    const Json::Value &w = out["patterns"][0];
    EXPECT_FALSE(w["retrieval_complete"].asBool()) << w;
    EXPECT_TRUE(w["withheld"]["reason"].asString() == "annotation_budget"
                || w["withheld"]["reason"].asString() == "output_budget") << w["withheld"];
    out = run(idx, body("{\"dna\": \"AC\"}", "\"mode\": \"partial\""), hooks);
    const Json::Value &q = out["patterns"][0];
    EXPECT_FALSE(q["retrieval_complete"].asBool());
    EXPECT_TRUE(q["rows_refused"].size() || !q["stop"].isNull() || !q["cut"].isNull()) << q;
    EXPECT_LE(q["work"]["memory_bytes"].asUInt64(), peak - 1);
}

TEST(PatternRetrieval, WorkBudgetStopsTheReads) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}, {\"dna\": \"GACG\"}",
                                    "\"mode\": \"partial\", \"max_annotation_work\": 1"));
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("label_discovery", e["stop"]["phase"].asString());
    EXPECT_EQ("max_annotation_work", e["stop"]["reason"].asString());
    EXPECT_FALSE(e["retrieval_complete"].asBool());
    size_t read = 0, not_read = 0;
    for (const Json::Value &r : e["results"]) {
        read += r["labels_status"].asString() == "complete";
        not_read += r["labels_status"].asString() == "not_read";
    }
    // one row read (its read overshot the budget), every other stated as not read
    EXPECT_GE(read, 1u);
    EXPECT_GT(not_read, 0u);
    EXPECT_EQ(e["results"].size(), read + not_read);
    // the request's budget: the next pattern reads nothing
    const Json::Value &next = out["patterns"][1];
    EXPECT_EQ("max_annotation_work", next["stop"]["reason"].asString());
    EXPECT_EQ(0u, next["work"]["annotation_rows"].asUInt64());

    out = run(idx, body("{\"dna\": \"AC\"}", "\"max_annotation_work\": 1"));
    EXPECT_EQ("annotation_budget", out["patterns"][0]["withheld"]["reason"].asString());
}

TEST(PatternRetrieval, DeadlineInTheReads) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    // a virtual clock every read advances by far more than the budget
    const Clock::time_point start = Clock::now();
    auto virtual_ms = std::make_shared<double>(0);
    auto clock = [start, virtual_ms]() {
        return start + std::chrono::duration_cast<Clock::duration>(
                std::chrono::duration<double, std::milli>(*virtual_ms));
    };
    RetrievalHooks hooks;
    hooks.read_hook = [virtual_ms](size_t) { *virtual_ms += 1e6; };
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}, {\"dna\": \"GACG\"}",
                                    "\"mode\": \"partial\""), hooks, true, clock);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("label_discovery", e["stop"]["phase"].asString());
    EXPECT_EQ("time", e["stop"]["reason"].asString());
    EXPECT_EQ("time_limited", e["determinism"].asString());
    EXPECT_FALSE(e["retrieval_complete"].asBool());
    // the graph discovery completed before the clock moved
    EXPECT_EQ("exact", e["counts"]["contexts"]["relation"].asString());
    // the next pattern: after a time stop, as after any
    const Json::Value &next = out["patterns"][1];
    EXPECT_EQ("discovery", next["stop"]["phase"].asString());
    EXPECT_EQ("time", next["stop"]["reason"].asString());
    EXPECT_EQ("unknown", next["counts"]["contexts"]["relation"].asString());

    *virtual_ms = 0;
    out = run(idx, body("{\"dna\": \"AC\"}"), hooks, true, clock);
    EXPECT_EQ("deadline", out["patterns"][0]["withheld"]["reason"].asString());
    EXPECT_EQ("time_limited", out["patterns"][0]["determinism"].asString());
}

TEST(PatternRetrieval, UnbudgetedBackend) {
    Index budgeted = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Index column = build<annot::ColumnCompressed<>>(kK, kRecords, true);
    EXPECT_EQ(std::make_pair(400, std::string("annotation_unbudgeted")),
              refusal(column, body("{\"dna\": \"AC\"}")));
    // mode count reads nothing: no refusal
    EXPECT_EQ(std::make_pair(200, std::string()),
              refusal(column, body("{\"dna\": \"AC\"}", "\"mode\": \"count\"")));
    Json::Value a = run(budgeted, body("{\"dna\": \"AC\"}"));
    Json::Value b = run(column, body("{\"dna\": \"AC\"}", "\"allow_unbudgeted_annotation\": true"));
    const Json::Value &e = b["patterns"][0];
    EXPECT_EQ("unbudgeted", e["annotation"].asString());
    ASSERT_EQ(1u, e["notes"].size());
    EXPECT_EQ("annotation_unbudgeted", e["notes"][0].asString());
    EXPECT_TRUE(b["limits"]["allow_unbudgeted_annotation"].asBool());
    check_complete_record_answer(column, e, "AC");
    // the same labels and occurrences as on the budgeted backend (node ids may differ)
    auto labels_of = [](const Json::Value &entry) {
        std::map<std::tuple<std::string, std::string, uint64_t>, Json::Value> m;
        for (const Json::Value &r : entry["results"]) {
            m[{ r["strand"].asString(), r["kmer"].asString(), r["offset"].asUInt64() }]
                    = r["labels"];
        }
        return m;
    };
    EXPECT_EQ(labels_of(a["patterns"][0]), labels_of(e));
    EXPECT_EQ(a["patterns"][0]["by_label"], e["by_label"]);
    EXPECT_EQ(a["patterns"][0]["counts"], e["counts"]);
}

TEST(PatternRetrieval, GlobalPlacementWithoutRecordMapping) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}"), {}, false);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("global", e["placement"].asString());
    ASSERT_EQ(1u, e["notes"].size());
    EXPECT_EQ("record_bounds_unknown", e["notes"][0].asString());
    // nothing is placed in a record: no occurrence count is claimed
    EXPECT_EQ("unknown", e["counts"]["occurrences"]["relation"].asString());
    EXPECT_EQ("exact", e["counts"]["labels"]["relation"].asString());
    EXPECT_TRUE(e["retrieval_complete"].asBool());
    size_t entries = 0;
    for (const Json::Value &r : e["results"]) {
        for (const Json::Value &l : r["labels"]) {
            EXPECT_FALSE(l.isMember("occurrences"));
            const std::string column = l["column"].asString();
            for (const Json::Value &o : l["occurrence_list"]) {
                entries++;
                EXPECT_EQ(r["offset"].asUInt64(), o["offset"].asUInt64());
                EXPECT_EQ(r["strand"].asString(), o["strand"].asString());
                // kmer_coord is the column coordinate of the context's own k-mer
                const uint64_t c = o["kmer_coord"].asUInt64();
                bool found = false;
                for (size_t i = 0; i < idx.records.size(); ++i) {
                    const Record &rec = idx.records[i];
                    const uint64_t n = rec.seq.size() - kK + 1;
                    if (rec.column == column && idx.starts[i] <= c && c < idx.starts[i] + n) {
                        EXPECT_EQ(r["kmer"].asString(), rec.seq.substr(c - idx.starts[i], kK));
                        found = true;
                    }
                }
                EXPECT_TRUE(found) << column << " " << c;
            }
        }
    }
    EXPECT_GT(entries, 0u);
    for (const Json::Value &b : e["by_label"]) {
        EXPECT_EQ("unknown", b["occurrences"]["relation"].asString());
    }
}

TEST(PatternRetrieval, NoCoordinatesNoPlacement) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, false);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}"), {}, false);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("none", e["placement"].asString());
    EXPECT_TRUE(e["retrieval_complete"].asBool());
    EXPECT_EQ("unknown", e["counts"]["occurrences"]["relation"].asString());
    const auto by_column = contexts_by_column(idx, "AC");
    EXPECT_EQ(by_column.size(), e["counts"]["labels"]["value"].asUInt64());
    for (const Json::Value &r : e["results"]) {
        const Context c { r["strand"].asString(), r["kmer"].asString(), r["offset"].asUInt64() };
        std::set<std::string> expected, got;
        for (const auto &[column, set] : by_column) {
            if (set.count(c))
                expected.insert(column);
        }
        for (const Json::Value &l : r["labels"]) {
            got.insert(l["column"].asString());
            EXPECT_FALSE(l.isMember("occurrence_list"));
            EXPECT_FALSE(l.isMember("occurrences"));
        }
        EXPECT_EQ(expected, got);
    }
}

TEST(PatternRetrieval, OccurrencesNotRequested) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, "{\"patterns\": [{\"dna\": \"AC\"}], "
                               "\"output\": {\"labels\": \"all\", \"occurrences\": false}}");
    EXPECT_FALSE(out["output"]["occurrences"].asBool());
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("not_requested", e["placement"].asString());
    EXPECT_TRUE(e["retrieval_complete"].asBool());
    EXPECT_EQ("unknown", e["counts"]["occurrences"]["relation"].asString());
    for (const Json::Value &r : e["results"]) {
        for (const Json::Value &l : r["labels"]) {
            EXPECT_FALSE(l.isMember("occurrence_list"));
        }
    }
}

TEST(PatternRetrieval, PartialCutsLabelsAndOccurrences) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value full = run(idx, body("{\"dna\": \"AC\"}"));
    const Json::Value &f = full["patterns"][0];
    ASSERT_GT(f["by_label"].size(), 1u);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"mode\": \"partial\", \"max_labels\": 1, "
                                                     "\"max_occurrences_per_label\": 1"));
    const Json::Value &e = out["patterns"][0];
    EXPECT_FALSE(e["retrieval_complete"].asBool());
    EXPECT_EQ("max_labels", e["labels_cut"]["reason"].asString());
    EXPECT_EQ(1u, e["labels_cut"]["returned"].asUInt64());
    ASSERT_EQ(1u, e["by_label"].size());
    // the first label of the complete answer's order, with its counts: the cut lists, the
    // counts stay whole
    EXPECT_EQ(f["by_label"][0], e["by_label"][0]);
    EXPECT_EQ(f["counts"], e["counts"]);
    EXPECT_EQ("max_occurrences_per_label", e["occurrences_cut"]["reason"].asString());
    const std::string kept = e["by_label"][0]["column"].asString();
    std::set<std::string> listed;
    for (const Json::Value &r : e["results"]) {
        for (const Json::Value &l : r["labels"]) {
            EXPECT_EQ(kept, l["column"].asString());
            for (const Json::Value &o : l["occurrence_list"]) {
                listed.insert(o["seq_id"].asString() + ":" + o["nt_coords"].asString()
                              + o["strand"].asString());
            }
        }
    }
    // the label's union lists one occurrence (each context that holds it lists it)
    EXPECT_EQ(1u, listed.size());
}

TEST(PatternRetrieval, CountModeAndProjectionNoneReadNothing) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    size_t reads = 0;
    RetrievalHooks hooks;
    hooks.read_hook = [&reads](size_t) { reads++; };
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"mode\": \"count\""), hooks);
    EXPECT_TRUE(out["output"].isNull());
    const Json::Value &e = out["patterns"][0];
    ASSERT_EQ(1u, e["notes"].size());
    EXPECT_EQ("annotation_not_read", e["notes"][0].asString());
    EXPECT_FALSE(e.isMember("placement"));
    EXPECT_FALSE(e.isMember("by_label"));
    EXPECT_FALSE(out["limits"].isMember("max_labels"));
    EXPECT_EQ("unknown", e["counts"]["labels"]["relation"].asString());

    out = run(idx, "{\"patterns\": [{\"dna\": \"AC\"}], \"max_labels\": 5}", hooks);
    EXPECT_EQ("none", out["output"]["labels"].asString());
    EXPECT_EQ(1u, out["output"].size());
    ASSERT_EQ(1u, out["patterns"][0]["notes"].size());
    EXPECT_EQ("annotation_not_read", out["patterns"][0]["notes"][0].asString());
    EXPECT_FALSE(out["patterns"][0]["results"][0].isMember("labels"));

    // a request without any annotation field: no note
    out = run(idx, "{\"patterns\": [{\"dna\": \"AC\"}]}", hooks);
    EXPECT_EQ(0u, out["patterns"][0]["notes"].size());
    EXPECT_EQ(0u, reads);
}

TEST(PatternRetrieval, ClampsAndRefusals) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"max_memory_mb\": 100000, "
                                    "\"max_labels_per_anchor\": 1000"));
    const Json::Value &clamped = out["limits"]["clamped"];
    ASSERT_EQ(2u, clamped.size());
    EXPECT_EQ("max_labels_per_anchor", clamped[0]["field"].asString());
    EXPECT_EQ(64u, clamped[0]["effective"].asUInt64());
    EXPECT_EQ("max_memory_mb", clamped[1]["field"].asString());
    EXPECT_EQ(256u, clamped[1]["effective"].asUInt64());
    EXPECT_EQ(256u, out["limits"]["max_memory_mb"].asUInt64());

    EXPECT_EQ(std::make_pair(400, std::string("later_increment")),
              refusal(idx, "{\"patterns\": [{\"dna\": \"AC\"}], "
                           "\"output\": {\"labels\": \"predicate_only\"}}"));
    EXPECT_EQ(std::make_pair(400, std::string("invalid_request")),
              refusal(idx, "{\"patterns\": [{\"dna\": \"AC\"}], "
                           "\"output\": {\"labels\": \"none\", \"occurrences\": true}}"));
    EXPECT_EQ(std::make_pair(400, std::string("invalid_request")),
              refusal(idx, body("{\"dna\": \"AC\"}", "\"max_labels_per_anchor\": 0")));
}

TEST(PatternRetrieval, WithheldByTheCountReadsNothing) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    size_t reads = 0;
    RetrievalHooks hooks;
    hooks.read_hook = [&reads](size_t) { reads++; };
    // more contexts than max_contexts: the count is not admitted, no annotation is read
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"max_contexts\": 2"), hooks);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("count_above_threshold", e["withheld"]["reason"].asString());
    EXPECT_EQ(0u, reads);
    EXPECT_TRUE(e["by_label"].isNull());
    EXPECT_EQ("unknown", e["counts"]["labels"]["relation"].asString());
    EXPECT_EQ(0u, e["work"]["annotation_rows"].asUInt64());
    // an absent pattern: complete, exact zeros
    out = run(idx, body("{\"dna\": \"TTTTTTT\"}"), hooks);
    const Json::Value &z = out["patterns"][0];
    EXPECT_TRUE(z["retrieval_complete"].asBool());
    EXPECT_EQ(0u, z["counts"]["labels"]["value"].asUInt64());
    EXPECT_EQ("exact", z["counts"]["labels"]["relation"].asString());
    EXPECT_EQ("exact", z["counts"]["occurrences"]["relation"].asString());
    EXPECT_EQ(0u, z["by_label"].size());
    EXPECT_TRUE(z["by_label"].isArray());
}

TEST(PatternRetrieval, CanonicalGraphLabelsWithoutPlacement) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, false,
                                                     DeBruijnGraph::CANONICAL);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}"), {}, false);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("none_canonical", e["placement"].asString());
    EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
    // a stored k-mer carries the columns of both its orientations
    for (const Json::Value &r : e["results"]) {
        const std::string kmer = r["kmer"].asString();
        std::set<std::string> expected, got;
        for (const Record &rec : idx.records) {
            if (rec.seq.find(kmer) != std::string::npos
                    || rec.seq.find(rc(kmer)) != std::string::npos)
                expected.insert(rec.column);
        }
        for (const Json::Value &l : r["labels"]) {
            got.insert(l["column"].asString());
        }
        EXPECT_EQ(expected, got) << kmer;
    }
}

TEST(PatternRetrieval, PrimaryGraphLabels) {
    Index idx = build<annot::ColumnCompressed<>>(kK, kRecords, false, DeBruijnGraph::PRIMARY);
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"allow_unbudgeted_annotation\": true"),
                          {}, false);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("primary", out["index"]["graph_mode"].asString());
    EXPECT_EQ("none_canonical", e["placement"].asString());
    EXPECT_TRUE(e["retrieval_complete"].asBool()) << e;
    ASSERT_GT(e["results"].size(), 0u);
    for (const Json::Value &r : e["results"]) {
        const std::string kmer = r["kmer"].asString();
        std::set<std::string> expected, got;
        for (const Record &rec : idx.records) {
            if (rec.seq.find(kmer) != std::string::npos
                    || rec.seq.find(rc(kmer)) != std::string::npos)
                expected.insert(rec.column);
        }
        for (const Json::Value &l : r["labels"]) {
            got.insert(l["column"].asString());
        }
        EXPECT_EQ(expected, got) << kmer;
    }
}

// The answer states why an entry is incomplete (the owner's guarantee rule)
void expect_incomplete_stated(const Json::Value &e) {
    if (e.isMember("error") || e["retrieval_complete"].asBool())
        return;
    EXPECT_TRUE(!e["withheld"].isNull() || !e["cut"].isNull() || !e["stop"].isNull()
                || e["rows_refused"].size() || e["anchors_truncated"].size()
                || !e["labels_cut"].isNull() || !e["occurrences_cut"].isNull()) << e;
}

const uint64_t kDescriptor = 512 + 2 * kK;
const uint64_t kStatement = 384 + kK;
// a label name of the dictionary: 192 + 3 x its length ("c1")
const uint64_t kName = 192 + 3 * 2;

// Review finding (unsigned wrap): an unbudgeted read's label names are held whether or not
// they fit; once they took the account past its maximum, `max - held` wrapped and every later
// charge succeeded, so the account bounded nothing (in the review's probe the next pattern
// listed 35 contexts with labels at 48,647 bytes against a maximum of 626). Now nothing fits
// while it is past the maximum, and the reads stop, stated.
TEST(PatternRetrieval, UnbudgetedNamesPastTheMaximumStopTheReads) {
    Index idx = build<annot::ColumnCompressed<>>(kK, kRecords, true);
    // GACGACT: one context, its row carries c1 and c4 (two names, 2 x 198 bytes). The account
    // holds its descriptor and leaves 393 bytes: enough to reserve one statement (391) for the
    // read, not for the two names the read returns, which take it 3 bytes past its maximum.
    RetrievalHooks hooks;
    hooks.max_memory_bytes = kDescriptor + kStatement + 2;
    ASSERT_LT(hooks.max_memory_bytes - kDescriptor, 2 * kName);
    Json::Value out = run(idx, body("{\"dna\": \"GACGACT\"}, {\"dna\": \"AC\"}",
                                    "\"allow_unbudgeted_annotation\": true"), hooks);
    const Json::Value &first = out["patterns"][0];
    EXPECT_EQ("label_discovery", first["stop"]["phase"].asString()) << first;
    EXPECT_EQ("max_memory", first["stop"]["reason"].asString());
    EXPECT_EQ("annotation_budget", first["withheld"]["reason"].asString());
    EXPECT_EQ(0u, first["results"].size());
    // past the maximum by the forced names, no more
    EXPECT_EQ(kDescriptor + 2 * kName, first["work"]["memory_bytes"].asUInt64());
    // the next pattern: the names stay with the dictionary, its first descriptor does not fit
    // what is left, and it is said
    const Json::Value &next = out["patterns"][1];
    EXPECT_FALSE(next["retrieval_complete"].asBool());
    EXPECT_EQ("output_budget", next["withheld"]["reason"].asString()) << next;
    EXPECT_EQ("output", next["stop"]["phase"].asString());
    EXPECT_EQ("max_memory", next["stop"]["reason"].asString());
    EXPECT_EQ(0u, next["work"]["annotation_rows"].asUInt64());
    EXPECT_EQ(first["work"]["memory_bytes"], next["work"]["memory_bytes"]);
}

// The memory account bounds the request (§5.3): over a sweep of maxima, both modes, a
// budgeted backend (with every read admitted, or every read refused) and an unbudgeted one, no
// entry's peak passes the maximum (the unbudgeted reads: by the names of one read at most,
// bounded here by the whole dictionary's), and every incomplete entry says why.
TEST(PatternRetrieval, MemoryAccountBoundsEveryMaximum) {
    Index budgeted = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Index column = build<annot::ColumnCompressed<>>(kK, kRecords, true);
    const uint64_t names = 4 * kName;
    enum Backend { BUDGETED, REFUSING, UNBUDGETED };
    for (Backend backend : { BUDGETED, REFUSING, UNBUDGETED }) {
        const Index &idx = backend == UNBUDGETED ? column : budgeted;
        for (const char *mode : { "all_or_count", "partial" }) {
            SCOPED_TRACE(std::string(mode) + " backend " + std::to_string(backend));
            const std::string b = body("{\"dna\": \"GACGACT\"}, {\"dna\": \"AC\"}, "
                                       "{\"dna\": \"GACG\"}",
                                       "\"allow_unbudgeted_annotation\": true, \"mode\": \""
                                       + std::string(mode) + "\"");
            uint64_t complete = 0, stated = 0;
            for (uint64_t max = 500; max < 40000; max += max < 2500 ? 1 : 97) {
                RetrievalHooks hooks;
                hooks.max_memory_bytes = max;
                if (backend == REFUSING)
                    hooks.deny_decode = [](uint64_t) { return true; };
                Json::Value out = run(idx, b, hooks);
                uint64_t peak = 0;
                for (const Json::Value &e : out["patterns"]) {
                    // the account's peak is the request's: non-decreasing over its patterns
                    EXPECT_GE(e["work"]["memory_bytes"].asUInt64(), peak);
                    peak = e["work"]["memory_bytes"].asUInt64();
                    expect_incomplete_stated(e);
                    complete += e["retrieval_complete"].asBool();
                    stated += !e["retrieval_complete"].asBool();
                }
                ASSERT_LE(peak, max + (backend == UNBUDGETED ? names : 0))
                        << "max " << max << "\n" << out;
            }
            EXPECT_GT(stated, 0u);
            if (backend != REFUSING) {
                EXPECT_GT(complete, 0u);
            }
        }
    }
}

// Review finding (descriptors): partial filled the whole account with the contexts'
// descriptors before reading any label, so every row was refused and a memory cut returned
// contexts without a single label. Now the descriptors take at most half of what the account
// has left, a context's descriptor is charged before its result object is built (none built
// past the cut), and the rows have the other half.
TEST(PatternRetrieval, PartialMemoryCutKeepsLabelledContexts) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    Json::Value full = run(idx, body("{\"dna\": \"AC\"}", "\"mode\": \"partial\""));
    const uint64_t all = full["patterns"][0]["returned"].asUInt64();
    ASSERT_GT(all, 4u);
    // the account holds the descriptors of every context, no more
    RetrievalHooks hooks;
    hooks.max_memory_bytes = all * kDescriptor;
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"mode\": \"partial\""), hooks);
    const Json::Value &e = out["patterns"][0];
    EXPECT_EQ("max_memory", e["cut"]["reason"].asString()) << e;
    EXPECT_EQ(all / 2, e["returned"].asUInt64());
    EXPECT_EQ(e["returned"].asUInt64(), e["results"].size());
    EXPECT_LE(e["work"]["memory_bytes"].asUInt64(), hooks.max_memory_bytes);
    uint64_t labelled = 0;
    for (const Json::Value &r : e["results"]) {
        labelled += r["labels_status"].asString() == "complete" && r["labels"].size() > 0;
    }
    EXPECT_GT(labelled, 0u) << e;
    // the counts over the graph stay the engine's
    EXPECT_EQ(full["patterns"][0]["counts"]["contexts"], e["counts"]["contexts"]);
}

// The statements of refused rows are in the account (priced, reserved before the read that
// can refuse): with every read refused it holds the descriptors and one statement per row.
TEST(PatternRetrieval, RefusedRowStatementsAreCharged) {
    Index idx = build<annot::RowDiffColumnAnnotator>(kK, kRecords, true);
    RetrievalHooks hooks;
    hooks.deny_decode = [](uint64_t) { return true; };
    Json::Value out = run(idx, body("{\"dna\": \"AC\"}", "\"mode\": \"partial\""), hooks);
    const Json::Value &e = out["patterns"][0];
    const uint64_t contexts = e["returned"].asUInt64();
    const uint64_t rows = e["rows_refused"].size();
    ASSERT_GT(rows, 1u);
    EXPECT_EQ(contexts * kDescriptor + rows * kStatement, e["work"]["memory_bytes"].asUInt64());
    // all_or_count: the first refusal ends the reads; its statement is held with the
    // descriptors until the pattern is withheld (the statement stays with the answer)
    out = run(idx, body("{\"dna\": \"AC\"}"), hooks);
    const Json::Value &w = out["patterns"][0];
    ASSERT_EQ(1u, w["rows_refused"].size());
    EXPECT_EQ(contexts * kDescriptor + kStatement, w["work"]["memory_bytes"].asUInt64());
}

} // namespace
