#ifndef __METAGRAPH_SERVER_CHECKS_HPP__
#define __METAGRAPH_SERVER_CHECKS_HPP__

#include <functional>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <json/json.h>


namespace mtg {
namespace cli {

// The parts of the server's request handling that need no HTTP library (server_utils.cpp),
// declared apart so that they are tested on their own

// Answers a request with |status| and the JSON |body| (an error that states more than its
// message: the usage of a ledger-managed /traverse, the state of a conflicting attempt)
class HttpError : public std::runtime_error {
  public:
    HttpError(int status, Json::Value body)
          : std::runtime_error(body.get("error", "").asString()), status_(status),
            body_(std::move(body)) {}
    int status() const { return status_; }
    const Json::Value& body() const { return body_; }

  private:
    int status_;
    Json::Value body_;
};

// The client of the request is gone: nothing is written, and the connection is closed
class ClientGone : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

class CurrentlyInitializingError : public std::runtime_error {
  public:
    CurrentlyInitializingError()
        : std::runtime_error("Server is currently initializing") {}
};

// What process_request does besides writing the result (null: nothing)
struct ResponseControl {
    // called while the body is written (every 64 KiB of JSON text) and compressed (every
    // block): throws to stop — ClientGone (nothing is written) or HttpError (that answer)
    std::function<void()> check;
    // after the handler is done: the status written and the body's size in bytes, or 0 and
    // nullopt when nothing was written (the client is gone)
    std::function<void(int status, std::optional<size_t> bytes)> on_written;
    // the zlib level of a compressed body (1-9; 9, the best compression, for every route
    // that does not choose: the traversal routes choose a faster one, see
    // Config::traverse_compression_level)
    int compression_level = 9;
    // after a body of |text_bytes| was compressed in |seconds| (a delivery rate measured)
    std::function<void(size_t text_bytes, double seconds)> on_compressed;
    // writes the result as text instead of json_text(result, compact, check) — a route that
    // wrote parts of it already (the /traverse results, each written once built) assembles
    // them here; must give the same bytes json_text would
    std::function<std::string(const Json::Value &result, const std::function<void()> &check)> write;
    // A route that answers nobody who left (/pattern, SPEC-pattern-search.md §3: a client
    // that is gone — a half-close counts — or a server that stops gets nothing written):
    // asked once every answer is built, an error's included (a refusal, a 400 of a malformed
    // body, a 503 while the index loads or at the deadline, an unexpected failure), and
    // before a byte of it is written; true: nothing is written and the connection is closed,
    // as for ClientGone. Not the deadline, which |check| reads: a 503 at the deadline is
    // still written to a client that is there. Unset (every other route): every error is
    // written; a check reaches only the success path
    std::function<bool()> gone;
};

// The response of one request as process_request writes it: |status| (0: nothing is
// written, the connection is closed — a ClientGone, or ResponseControl::gone), the header
// fields after the Content-Type ("application/json", every response's first) in the order
// they are added, and the body
struct RequestAnswer {
    int status = 0;
    std::vector<std::pair<std::string, std::string>> header;
    std::string body;
};

/**
 * What process_request answers, without the HTTP library (server_utils.hpp): runs |process|
 * on |content| and writes its JSON result (|control|->write, else json_text(result, |compact|)
 * under |control|->check), compressed when |encoding| is "gzip" or "deflate" (the request's
 * requested_encoding; "" none), then asks the check once more. Its failures are answered:
 * ClientGone with nothing (status 0); HttpError with its status and body (uncompressed, as
 * every error); CurrentlyInitializingError with 503 and Retry-After; any other exception with
 * 400 {"error": its message}, anything else with 500. Last, |control|->gone, when set: true
 * answers nothing (status 0) whatever was built. |request_id| names the request in the log.
 */
RequestAnswer answer_request(const std::string &content, const std::string &encoding,
                             size_t request_id,
                             const std::function<Json::Value(const std::string &)> &process,
                             bool compact = false,
                             const ResponseControl *control = nullptr);

// Whether the peer of the connected TCP socket |fd| is gone: a non-blocking peek (consuming
// nothing) finds the connection closed — an orderly close or a half-close (the peer will send
// nothing more; nginx's default treats it the same) — or reset, or |fd| holds no connection (a
// negative or closed descriptor, or one that is no socket). Data waiting (a pipelined request,
// a trailing CRLF) hides a close behind it from the peek, so then the kernel's TCP state is
// read (Linux TCP_INFO, macOS TCP_CONNECTION_INFO): any state past ESTABLISHED is gone — the
// peer's FIN (CLOSE_WAIT) or reset (CLOSED), and the server's own shutdown at its content
// timeout (FIN_WAIT1/2 and the states after it), in which nothing can be delivered either —,
// ESTABLISHED (and SYN_RECV, a TCP Fast Open connection's, which this server does not enable)
// connected. On another platform, or for a socket that is not TCP, data waiting reads
// "connected" (a close behind it is then not seen until the data is read). Nothing to read
// yet: connected. An error that says nothing about the peer reads "connected".
bool peer_closed(int fd);

// |value| as the server writes it: compact (no indentation) or the default writer's
// indentation, byte for byte what Json::writeString writes; with |check|, called every 64 KiB
// of text (however large the pieces the writer hands over: they are copied in pieces up to the
// next check), whose exception stops the writing and reaches the
// caller. What stays uninterruptible is what the writer does between two pieces — preparing
// one token, e.g. escaping one string value of 16 MiB, before it is copied — and |max_gap_ms|,
// when given, receives the longest time between two checks (and from the start to the first,
// and from the last to the end) for observed_max_uninterruptible_ms
std::string json_text(const Json::Value &value, bool compact,
                      const std::function<void()> &check = nullptr,
                      double *max_gap_ms = nullptr);

// The compact text of a /traverse response whose results were written as text one by one
// (ResultTexts): |envelope|'s members before "results" and after it (jsoncpp writes an
// object's members in byte order of their names), with "results":[t0,t1,...] between them —
// byte for byte json_text(envelope with results, true); |check| and |max_gap_ms| as
// json_text's, the results' texts copied in pieces of at most 64 KiB with a check between
// them. |results| is taken by value (the server moves its texts in): the response's exact
// size is reserved before anything is copied, and each text is freed once copied, so that
// assembling holds about the response once, not the texts beside a copy that grows by
// doubling (up to 3x the text before compression or the transport's copy begins)
std::string assemble_traverse_response(const Json::Value &envelope,
                                       std::vector<std::string> results,
                                       const std::function<void()> &check = nullptr,
                                       double *max_gap_ms = nullptr);

// |text| compressed by zlib at |level| (1-9), in a gzip container when |gzip|, else a zlib
// stream; |check| called before each 32 KiB block of output, its exception (after the stream
// is released) reaching the caller. zlib counts its input in 32 bits: the text is handed over
// in pieces of at most |max_piece| bytes (0: the most zlib takes, 2^32 - 1), so a text of any
// size is compressed whole; a text of one piece is compressed in one call (|max_piece|
// below that is for tests)
std::string compress_string(const std::string &text, int level, bool gzip,
                            const std::function<void()> &check = nullptr,
                            size_t max_piece = 0);

// The interval of the checks of json_text and assemble_traverse_response (bytes of text)
constexpr size_t kDeliveryCheckBytes = size_t(1) << 16;

/**
 * One line of a multi-graph list (`server_query GRAPHS.csv`):
 *     name,graph_path,annotation_path[,manifest_path[,index_ns]]
 * split on every comma; the first three columns are required, the last two optional — the
 * per-graph identity of
 * DESIGN-traverse-graphlet.md §16.1: the manifest of the pair's bundle, checked at start-up as
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
// five (a column nothing reads is refused rather than ignored), or an index_ns that is not
// [A-Za-z0-9._-]+
GraphListEntry parse_graph_list_line(const std::string &text, size_t line);

/**
 * The identity each (graph, annotation) pair of a graph list states: (index_ns, index_fp),
 * "" for none. |fingerprint| returns the digest of an entry's manifest checked against the
 * entry's files (index_manifest_fingerprint), called once per (bundle, manifest). One index,
 * one identity: lines naming the same index must agree on what they state — an empty column
 * states nothing and takes what another line states — and a conflict throws
 * std::invalid_argument naming both lines. Two lines name the same index when the complete
 * loader inventories of their pairs (|inventory|, by default index_bundle_files of the listed
 * spellings: the identity files, never the graph's derived mask or Bloom filter) resolve to the
 * same real paths in the same roles; two different indexes stating
 * one index_fp are refused too.
 */
std::map<std::pair<std::string, std::string>, std::pair<std::string, std::string>>
graph_list_identities(const std::vector<GraphListEntry> &entries,
                      const std::function<std::string(const GraphListEntry &)> &fingerprint,
                      const std::function<std::vector<std::string>(const GraphListEntry &)>
                              &inventory = nullptr);

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_SERVER_CHECKS_HPP__
