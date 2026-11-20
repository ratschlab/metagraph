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
// the code of a /search refused whole (400) or of a unit refused alone (its cut entry's
// `refused`)
constexpr char kSearchAnswerTooLarge[] = "search_answer_too_large";

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
// |tagged|: the rows and the cut carry the pair's tags, |graph| and |pair_id| (a multi-graph
// server); |index_fp| empty: none (stated null)
SearchRecordPrice search_record_price(const SeqSearchResult &result, bool verbose_output,
                                      bool tagged, const std::string &graph,
                                      const std::string &pair_id,
                                      const std::string &index_fp);

// What measuring one record's price (search_record_price) holds while it runs: the coordinate
// ranges collapsed_coords_text keeps open for a query_coords row (not |verbose_output|), at
// most coord_ranges_bytes of the record's longest coordinate list; the other measures hold
// nothing (REVIEW 9726cb420 [P2], lane T6-W2)
uint64_t search_record_scratch(const SeqSearchResult &result, bool verbose_output);

// Charges |answer| (the request's answer account) for one record's object before the server
// builds it, in two steps: measure_search_record charges the scratch measuring its price holds
// (search_record_scratch), before it is measured, and returns the price; charge_search_price
// charges its tree and its transient, the transient (and the scratch, which it covers) as the
// most any record of the pair took so far, |*transient_charged|, which both raise (the
// objects of a pair are built one at a time). charge_search_record makes both. Past the bound
// they throw graph::SearchAnswerTooLarge
SearchRecordPrice measure_search_record(graph::SearchAccount &answer, uint64_t *transient_charged,
                                        const SeqSearchResult &result, bool verbose_output,
                                        bool tagged, const std::string &graph,
                                        const std::string &pair_id,
                                        const std::string &index_fp);
void charge_search_price(graph::SearchAccount &answer, uint64_t *transient_charged,
                         const SearchRecordPrice &price);
SearchRecordPrice charge_search_record(graph::SearchAccount &answer, uint64_t *transient_charged,
                                       const SeqSearchResult &result, bool verbose_output,
                                       bool tagged, const std::string &graph,
                                       const std::string &pair_id,
                                       const std::string &index_fp);

// One unit's object (a FASTA record on a pair; lane T7-P1), as process_search_request charges
// each in turn before it builds it: the unit's object check, its own price against the bound
// (first the scratch of measuring it alone, search_record_scratch, before it is measured or
// charged; then measure_search_record and its tree and transient, saturating): past it, with
// the request's `partial_refusal` (|budget|) the unit is refused alone (|*result| becomes its
// stub, SeqSearchResult::refuse with the account "answer" and the price in MiB rounded up, the
// scratch's when that alone passed, and the refusal is noted for the log) and the stub's price
// is charged; without it the request is refused on the answer side (SearchBudget::
// note_object_passed) and false is returned, no object to build. Otherwise (and for a stub
// refused on its values) the price is charged (charge_search_price) and the unit counted.
// |pair|: the pair as the log names it; |pair_index|: its place in the answer (the log's
// order). Past the request's bound throws graph::SearchAnswerTooLarge
bool charge_search_unit(graph::SearchBudget &budget, uint64_t *transient_charged,
                        SeqSearchResult *result, bool verbose_output, bool tagged,
                        const std::string &graph, const std::string &pair_id,
                        const std::string &index_fp, const std::string &pair,
                        size_t pair_index = 0);

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

// The refusal of a /search priced past the bound: {"error": the message, "account": |account|
// ("values" or "answer", SearchBudget::account()), "code": "search_answer_too_large",
// "max_memory_mb": the bound}, which the server writes error first (error_first_json_text)
// with status 400
Json::Value search_answer_too_large_json(uint64_t max_mb, const char *account);

// GET /capabilities' `search` block (lane T7-P1): {"handler_threads": the HTTP handler threads,
// "max_memory_mb": the bound (null when 0, none), "pair_searches_at_once": the pair searches
// the server runs at once over every request, "records_at_once": the records one pair search
// runs at once, "rule": kSearchBudgetRule, "too_large": {"status": 400, "code":
// "search_answer_too_large", "accounts": ["answer", "values"], "members": ["account", "code",
// "error", "max_memory_mb"]}, "unit_refused": {"field": "partial_refusal", "code":
// "search_answer_too_large", "accounts": ["answer", "values"], "members": ["account", "code",
// "max_memory_mb", "priced_mb"]}}
Json::Value search_capabilities_json(uint64_t max_mb, uint64_t pair_searches_at_once,
                                     uint64_t records_at_once, uint64_t handler_threads);

} // namespace cli
} // namespace mtg

#endif // __SEARCH_ANSWER_HPP__
