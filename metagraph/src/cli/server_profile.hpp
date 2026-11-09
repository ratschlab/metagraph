#ifndef __METAGRAPH_SERVER_PROFILE_HPP__
#define __METAGRAPH_SERVER_PROFILE_HPP__

#include <string>
#include <vector>

#include <json/json.h>

#include "graph/annotated_dbg.hpp"


namespace mtg {
namespace cli {

struct IndexIdentity;

// The capability profile of a (graph, annotation) pair (E17, SPEC-multigraph-envelope.md
// section 8): GET /capabilities' graph_summary.pairs[] entry of a multi-graph server, and the
// top-level pair_summary of a single-graph server, with the fingerprint of the pair's
// capability fields (support_fp)

// The annotation's type as `--anno-type` names it (Config::annotype_to_string), read from the
// extension of the annotation class the loader constructed (its file's extension)
std::string annotation_form(const graph::AnnotatedDBG::Annotator &annotation);

// The graph's representation as `--graph` names it, read from the extension of the graph class
// the loader constructed (under a CanonicalDBG, the graph it wraps): succinct, hash (an ordered
// hash graph, packed or not: one extension), hashstr, hashfast, sshash or bitmap
std::string graph_representation(const graph::DeBruijnGraph &graph);

// The first line of the text support_fp hashes, the version of its rule: a field added to the
// profile later changes it, and with it every support_fp (2: seqs_order, H4)
constexpr const char kSupportVersion[] = "mgt-support 2";

// The capability fields support_fp covers, by their names in a pair entry (a member of
// `traversal` as "traversal.<name>"), in byte order of the names: the order of the text
const std::vector<std::string>& support_fields();

// The text support_fp hashes: kSupportVersion, then per field of support_fields() the line
// "<name>=<value>", the value the entry's member as compact JSON (null when absent), each line
// ending with "\n"
std::string support_text(const Json::Value &pair);

// FNV-1a-64 of support_text(|pair|), 16 lowercase hex digits (as index_meta_fp)
std::string support_fingerprint(const Json::Value &pair);

/**
 * The profile of the pair |index| (SPEC-pattern-search.md §24.4, graph_summary_pair): its tags
 * (|graph|, |graph_path|, |annotation_path|: null on a single-graph server), its identity
 * (index_ns, index_fp, index_meta_fp from |identity|, whose meta fingerprint is set), k,
 * alphabet, graph_representation, annotation_form and decode_charged; seqs_order, what its load
 * verified of its record mapping's column order ("verified", "count_only", null without a
 * .seqs; SPEC-multigraph-envelope.md section 9.1); what the pattern search
 * makes of it (graph_mode, available, unavailable_reason, mask, counting, as its
 * GET /pattern/capabilities block states them); what a traversal reads (`traversal`: regime,
 * num_labels, has_coordinates, has_coord_to_header, supports_trace, counts, as every
 * /traverse response states them); and support_fp over the capability fields
 */
Json::Value pair_summary_json(const graph::AnnotatedDBG &index, const IndexIdentity &identity,
                              const Json::Value &graph, const Json::Value &graph_path,
                              const Json::Value &annotation_path);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_PROFILE_HPP__
