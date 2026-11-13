#ifndef __SEARCH_ANSWER_HPP__
#define __SEARCH_ANSWER_HPP__

#include <cstdint>
#include <string>

#include <json/json.h>

#include "graph/annotated_dbg.hpp"
#include "graph/search_budget.hpp"


namespace mtg {
namespace cli {

class SeqSearchResult;

/**
 * The /search answer budget (SPEC-multigraph-envelope.md section 3; the server's
 * --search-max-answer-memory-mb): the prices of what a /search answer holds, charged to the
 * request's accounts (graph::SearchBudget) before it is built. A model of the deployed
 * build's heap (libstdc++ and jemalloc, alloc_bytes; jsoncpp 1.9.3's nodes and strings as
 * /resolve's answer account prices them, cli/traverse.cpp), never a measurement.
 */

// the rule a client reads the budget by (GET /capabilities' search.rule)
constexpr char kSearchBudgetRule[]
        = "SPEC-multigraph-envelope.md section 3, the /search answer budget";

// What one FASTA record's result holds beyond its values' buffers (which the values account
// was charged with as they grew): its place in the pair's results, its name, sequence and
// alignment, its rows' entries and labels (the pair's values account)
uint64_t search_record_rows_bytes(const SeqSearchResult &result);

// What one FASTA record's object of the answer costs (the answer account): |tree|, its JSON
// tree as SeqSearchResult::to_json and the server build it (its rows, the pair tags on each
// row and on its cut, the cut, its node in the merged answer on a multi-graph server);
// |transient|, the most that building one of its rows holds beside the tree (the run strings
// before they are copied into it)
struct SearchRecordPrice {
    uint64_t tree = 0;
    uint64_t transient = 0;
};
// |tagged|: the rows and the cut carry the pair's tags (a multi-graph server); |index_fp|
// empty: none (stated null)
SearchRecordPrice search_record_price(const SeqSearchResult &result, bool verbose_output,
                                      bool tagged, const std::string &graph,
                                      const std::string &graph_path,
                                      const std::string &annotation_path,
                                      const std::string &index_fp);

// The bytes of the per-position strings SeqSearchResult::to_json writes for a row (with
// |verbose| as --verbose-output), computed without writing them: its signature's presence_mask,
// its kmer_abundances, its kmer_coords
uint64_t presence_mask_text_length(const graph::PresenceRuns &mask, bool verbose);
uint64_t count_runs_text_length(const graph::CountRuns &counts, bool verbose);
uint64_t coords_text_length(const graph::PositionCoords &coords, bool verbose);

// The text a writer of |builder| writes for |root| (as Json::writeString(builder, root)):
// its length, counted without storing it; and the text, written into a string reserved to
// |length| bytes (the bytes are Json::writeString's)
uint64_t json_text_length(const Json::Value &root, const Json::StreamWriterBuilder &builder);
std::string json_text_of_length(const Json::Value &root,
                                const Json::StreamWriterBuilder &builder, uint64_t length);

// The /search answer's text (Json::writeString(Json::StreamWriterBuilder(), |answer|),
// byte for byte), charged to |account| (the request's answer account) before it is written:
// its buffer beside the tree, then, where that is more than the tree held, what is held once
// the tree is freed: the text beside the transport's copy (or gzip's compressed copy and its
// state) and, with the server's cache on (|cached|), the cache's copy
std::string search_answer_text(const Json::Value &answer, graph::SearchAccount &account,
                               bool cached);

// The refusal of a /search priced past the bound: {"error": the message, "code":
// "search_answer_too_large", "max_memory_mb": the bound}, which the server writes error first
// (error_first_json_text) with status 400
Json::Value search_answer_too_large_json(uint64_t max_mb);

// GET /capabilities' `search` block: {"max_memory_mb": the bound (null when 0, none),
// "rule": kSearchBudgetRule, "too_large": {"status": 400, "code": "search_answer_too_large",
// "members": ["code", "error", "max_memory_mb"]}}
Json::Value search_capabilities_json(uint64_t max_mb);

} // namespace cli
} // namespace mtg

#endif // __SEARCH_ANSWER_HPP__
