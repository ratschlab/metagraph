#ifndef __STATS_HPP__
#define __STATS_HPP__

#include <cstdlib>
#include <string>
#include <vector>

#include <json/json.h>

namespace mtg {

namespace graph {
    class DeBruijnGraph;
    namespace boss {
        class BOSS;
    } // namespace boss
} // namespace graph

namespace annot {
    class CoordToHeader;
} // namespace annot

namespace cli {

class Config;


/**
 * The names of an index as the serving routes need them (`stats --check-names`): every label
 * name valid UTF-8 — a route refuses to list one that is not (unrepresentable_label_name),
 * since jsoncpp would replace its bytes by U+FFFD and a valid name holding U+FFFD could be
 * another label's — and, with the record mapping, every record header valid UTF-8 and distinct
 * within its column (CoordToHeader::find_header answers a header with its first record only).
 * One JSON object: labels, labels_invalid_utf8, headers_checked (a mapping was given), headers,
 * headers_invalid_utf8, duplicate_record_headers, first_invalid (at most 10 of {kind: label |
 * header, problem: invalid_utf8 | duplicate, column, record (0-based in its column; null for a
 * label), hex[, duplicate_of]}: the labels, then per column the headers that are not UTF-8,
 * then its repeats) and names_ok. Offending bytes are spelled in hex, never raw. One pass
 * over the headers with a per-column hash set of record indices.
 */
Json::Value index_names_report(const std::vector<std::string> &labels,
                               const annot::CoordToHeader *mapping);

void print_stats(const graph::DeBruijnGraph &graph, bool print_counts_hist = false);

void print_boss_stats(const graph::boss::BOSS &boss_graph,
                      bool count_dummy = false,
                      size_t num_threads = 0);

int print_stats(Config *config);

int compare(Config *config);

} // namespace cli
} // namespace mtg

#endif // __STATS_HPP__
