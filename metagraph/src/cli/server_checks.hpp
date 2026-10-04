#ifndef __METAGRAPH_SERVER_CHECKS_HPP__
#define __METAGRAPH_SERVER_CHECKS_HPP__

#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <json/json.h>


namespace mtg {
namespace cli {

// The parts of the server's request handling that need no HTTP library (server_utils.cpp),
// declared apart so that they are tested on their own

// Whether the peer of the connected TCP socket |fd| is gone: a non-blocking peek (consuming
// nothing) finds the connection closed — an orderly close or a half-close (the peer will send
// nothing more; nginx's default treats it the same) — or reset, or |fd| is no socket. Data
// waiting (a pipelined request, a trailing CRLF) hides a close behind it from the peek, so then
// the kernel's TCP state is read (Linux TCP_INFO, macOS TCP_CONNECTION_INFO): CLOSE_WAIT (the
// FIN arrived) or CLOSED (reset) is gone, otherwise connected. On another platform, or for a
// socket that is not TCP, data waiting reads "connected" (a close behind it is then not seen
// until the data is read). Nothing to read yet: connected. An error that says nothing about
// the peer reads "connected".
bool peer_closed(int fd);

// |value| as the server writes it: compact (no indentation) or the default writer's
// indentation, byte for byte what Json::writeString writes; with |check|, called every 64 KiB
// of text, whose exception stops the writing and reaches the caller
std::string json_text(const Json::Value &value, bool compact,
                      const std::function<void()> &check = nullptr);

// The compact text of a /traverse response whose results were written as text one by one
// (ResultTexts): |envelope|'s members before "results" and after it (jsoncpp writes an
// object's members in byte order of their names), with "results":[t0,t1,...] between them —
// byte for byte json_text(envelope with results, true); |check| as json_text's, also called
// every 64 KiB of the results' texts
std::string assemble_traverse_response(const Json::Value &envelope,
                                       const std::vector<std::string> &results,
                                       const std::function<void()> &check = nullptr);

/**
 * One line of a multi-graph list (`server_query GRAPHS.csv`):
 *     name,graph_path,annotation_path[,manifest_path[,index_ns]]
 * split on every comma, the three first columns as before (a line of three columns is read as
 * it always was), the two last optional — the per-graph identity of
 * DESIGN-traverse-graphlet.md §21: the manifest of the pair's bundle, checked at start-up as
 * --index-manifest is, and the name the pair's responses state. An empty optional column means
 * none (index_fp, index_ns null). Paths are as the server opens them (relative to its working
 * directory).
 */
struct GraphListEntry {
    size_t line = 0;                // 1-based, for messages
    std::string name;
    std::string graph_path;
    std::string annotation_path;
    std::string manifest_path;      // "" none
    std::string index_ns;           // "" none
};
// Throws std::invalid_argument naming the line's problem: fewer than three columns, more than
// five (a fourth column used to be dropped silently; it now names a manifest, so a column
// nothing reads is refused rather than ignored), or an index_ns that is not [A-Za-z0-9._-]+
GraphListEntry parse_graph_list_line(const std::string &text, size_t line);

/**
 * The identity each (graph, annotation) pair of a graph list states: (index_ns, index_fp),
 * "" for none. |fingerprint| returns the digest of an entry's manifest checked against the
 * entry's files (index_manifest_fingerprint), called once per (pair, manifest). One index,
 * one identity: lines naming the same pair must agree on what they state — an empty column
 * states nothing and takes what another line states — and a conflict throws
 * std::invalid_argument naming both lines.
 */
std::map<std::pair<std::string, std::string>, std::pair<std::string, std::string>>
graph_list_identities(const std::vector<GraphListEntry> &entries,
                      const std::function<std::string(const GraphListEntry &)> &fingerprint);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_CHECKS_HPP__
