#include <atomic>
#include <chrono>
#include <condition_variable>
#include <csignal>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <mutex>
#include <set>
#include <thread>

#include <unistd.h>

#include <json/json.h>
#include <tsl/hopscotch_map.h>
#include <server_http.hpp>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "query.hpp"
#include "align.hpp"
#include "traverse.hpp"
#include "traverse_attempts.hpp"
#include "pattern.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "server_utils.hpp"
#include "cli/json_helpers.hpp"
#include "cli/load/load_annotation.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using namespace mtg::graph;

using HttpServer = SimpleWeb::Server<SimpleWeb::HTTP>;

// How long the HTTP server gives a request after its header was read (the body, the handler
// and sending the response) before it shuts the connection: also the cap on the duration
// bound of a /traverse attempt (traverse_attempts.hpp), less a second for the response
constexpr long kContentTimeoutS = static_cast<long>(kServerContentTimeoutS);

// The rules of GET /traverse/capabilities (and GET /capabilities' deadline_check), stated by
// reference to the SPEC section that holds each: the documents a service returns in one piece
// have a ceiling of 32 KiB, and the numbers a client computes with are fields beside them
// (chunk_target_ms, poll_stride, work_check_interval, path_cache_mb). ASCII only: the writers
// escape any other byte as \uXXXX
constexpr char kDeadlineCheckRule[]
        = "SPEC-labeled-traversal-core.md section 6.8, chunked deadlines";
constexpr char kWorkBoundRule[] = "SPEC-labeled-traversal-core.md section 6.8, work units";
constexpr char kDecodeCacheRule[]
        = "SPEC-labeled-traversal-core.md section 8.4, the row-diff path cache";


Json::Value process_search_request(const Json::Value &json,
                                   const graph::AnnotatedDBG &anno_graph,
                                   const Config &config_orig) {
    const auto &fasta = json["FASTA"];
    if (fasta.isNull())
        throw std::invalid_argument("No input sequences received from client");

    Config config(config_orig);
    // discovery_fraction a proxy of 1 - %similarity
    config.discovery_fraction
            = json.get("discovery_fraction", config.discovery_fraction).asDouble();

    config.alignment_min_exact_match
            = json.get("min_exact_match",
                       config.alignment_min_exact_match).asDouble();

    config.alignment_max_nodes_per_seq_char = json.get(
        "max_num_nodes_per_seq_char",
        config.alignment_max_nodes_per_seq_char).asDouble();

    if (config.discovery_fraction < 0.0 || config.discovery_fraction > 1.0) {
        throw std::invalid_argument(
                "Discovery fraction should be within [0, 1.0]. Instead got "
                + std::to_string(config.discovery_fraction));
    }

    if (config.alignment_min_exact_match < 0.0
            || config.alignment_min_exact_match > 1.0) {
        throw std::invalid_argument(
                "Minimum exact match should be within [0, 1.0]. Instead got "
                + std::to_string(config.alignment_min_exact_match));
    }

    config.num_top_labels = json.get("top_labels", config.num_top_labels).asInt();

    if (json.get("query_coords", false).asBool()) {
        config.query_mode = COORDS;
    } else if (json.get("query_counts", false).asBool()) {
        config.query_mode = COUNTS;
    } else if (json.get("with_signature", false).asBool()) {
        config.query_mode = SIGNATURE;
    } else if (json.get("abundance_sum", false).asBool()) {
        config.query_mode = COUNTS_SUM;
    } else {
        config.query_mode = MATCHES;
    }

    // Throw client an error if they try to query coordinates/kmer-counts on unsupported indexes
    if ((config.query_mode == COUNTS || config.query_mode == COUNTS_SUM)
            && !dynamic_cast<const annot::matrix::IntMatrix *>(
                            &anno_graph.get_annotator().get_matrix())) {
        throw std::invalid_argument("Annotation does not support k-mer count queries");
    }

    if (config.query_mode == COORDS
            && !dynamic_cast<const annot::matrix::MultiIntMatrix *>(
                            &anno_graph.get_annotator().get_matrix())) {
        throw std::invalid_argument("Annotation does not support k-mer coordinate queries");
    }

    std::unique_ptr<align::DBGAlignerConfig> aligner_config;
    if (json.get("align", false).asBool()) {
        aligner_config.reset(new align::DBGAlignerConfig(
            initialize_aligner_config(config, anno_graph.get_graph())
        ));
    }

    // Need mutex while appending to vector
    std::vector<SeqSearchResult> search_results;
    std::mutex result_mutex;

    // writing to temporary file in order to reuse query code. This is not optimal and
    // may turn out to be an issue in production. However, adapting FastaParser to
    // work on strings seems non-trivial. An alternative would be to use
    // read_fasta_from_string.
    utils::TempFile tf(config.tmp_dir);
    tf.ofstream() << fasta.asString();
    tf.ofstream().close();

    // Query sequences and callback by appending result to vector with mutex for thread safety
    query_fasta(tf.name(),
        [&](const SeqSearchResult &result) {
            std::lock_guard<std::mutex> lock(result_mutex);
            search_results.emplace_back(std::move(result));
        },
        config, anno_graph, aligner_config.get()
    );

    // Ensure JSON results are sorted by their ID
    std::sort(search_results.begin(), search_results.end(),
              [](const SeqSearchResult &lhs, const SeqSearchResult &rhs) {
                  return lhs.get_sequence().id < rhs.get_sequence().id;
              });

    // Create full JSON object
    Json::Value search_response(Json::arrayValue);
    for (const auto &seq_result : search_results) {
        search_response.append(seq_result.to_json(config.verbose_output, anno_graph.get_graph().get_k()));
    }

    return search_response;
}

// TODO: implement alignment_result.to_json as in process_search_request
Json::Value process_align_request(const std::string &received_message,
                                  const graph::DeBruijnGraph &graph,
                                  const Config &config_orig) {
    Json::Value json = parse_json_string(received_message);

    const auto &fasta = json["FASTA"];

    Json::Value root = Json::Value(Json::arrayValue);

    Config config(config_orig);

    config.alignment_num_alternative_paths = json.get(
        "max_alternative_alignments",
        (uint64_t)config.alignment_num_alternative_paths).asInt();

    if (!config.alignment_num_alternative_paths) {
        // TODO: better throw an exception and send an error response to the client
        logger->warn("[Server] Got invalid value of alignment_num_alternative_paths = {}."
                     " The default value of 1 will be used instead...", config.alignment_num_alternative_paths);
        config.alignment_num_alternative_paths = 1;
    }

    config.alignment_min_exact_match
            = json.get("min_exact_match",
                       config.alignment_min_exact_match).asDouble();

    if (config.alignment_min_exact_match < 0.0
            || config.alignment_min_exact_match > 1.0) {
        throw std::invalid_argument(
                "Minimum exact match should be within [0, 1.0]. Instead got "
                + std::to_string(config.alignment_min_exact_match));
    }

    config.alignment_max_nodes_per_seq_char = json.get(
        "max_num_nodes_per_seq_char",
        config.alignment_max_nodes_per_seq_char).asDouble();

    align::DBGAligner aligner(graph, initialize_aligner_config(config, graph));
    const align::DBGAlignerConfig &aligner_config = aligner.get_config();

    // TODO: make parallel?
    seq_io::read_fasta_from_string(fasta.asString(),
                                   [&](seq_io::kseq_t *read_stream) {
        Json::Value align_entry;
        align_entry[SeqSearchResult::SEQ_DESCRIPTION_JSON_FIELD] = read_stream->name.s;

        // not supporting reverse complement yet
        Json::Value alignments = Json::Value(Json::arrayValue);

        for (const auto &path : aligner.align(read_stream->seq.s)) {
            Json::Value a;
            a[SeqSearchResult::SCORE_JSON_FIELD] = path.get_score();
            a[SeqSearchResult::MAX_SCORE_JSON_FIELD] = aligner_config.match_score(read_stream->seq.s)
                + aligner_config.left_end_bonus + aligner_config.right_end_bonus;
            aligner.get_config().match_score(read_stream->seq.s);
            a[SeqSearchResult::SEQUENCE_JSON_FIELD] = std::string(path.get_sequence());
            a[SeqSearchResult::CIGAR_JSON_FIELD] = path.get_cigar().to_string();
            a[SeqSearchResult::ORIENTATION_JSON_FIELD] = path.get_orientation();

            alignments.append(a);
        }

        align_entry[SeqSearchResult::ALIGNMENT_JSON_FIELD] = alignments;

        root.append(align_entry);
    });

    return root;
}

// ---------------------------------------------------------------- shutdown (SIGTERM, SIGINT)

namespace {

/**
 * The server's shutdown on SIGTERM or SIGINT. Without a handler the process kept running: in a
 * container it is PID 1, for which the kernel drops a signal with the default action, so every
 * `docker stop` of a deploy waited its whole timeout (120 s) and then killed it. Now:
 *  - the signal handler only writes a byte to a pipe (async-signal-safe); a thread reading it
 *    does the rest. A second signal ends the process at once (_exit);
 *  - nothing new is served: a /traverse or /resolve request already reading its request or
 *    arriving now finds its client "gone" at its first check;
 *  - every traversal in flight (each /traverse request has an Attempt, with attempt_id or not)
 *    is stopped at its next poll as if its client had gone away — nothing is written for it
 *    and its connection is closed. A partial result written now would state a cause that is
 *    not true ("cancelled": POST /traverse/cancel), and a ledger treats an attempt without an
 *    answer as unanswered (not_after_ms, bound_ms), which is the truth;
 *  - once they returned, or after kDrainMs, the HTTP server stops (its acceptor and every
 *    connection closed, a response still being sent cut: incomplete, never a shorter complete
 *    one, its Content-Length unmet) and the process exits from run_server (_Exit: returning
 *    would join an index load still in progress in the loader's destructor);
 *  - a handler that does not return within kExitMs more (a long /search, /align, or one
 *    uninterruptible annotation read) does not hold the process: the log is flushed and the
 *    process exits (_Exit, the state is in memory only).
 * Before the server accepts connections the process exits at once. In single-index mode the
 * server accepts them at once and answers 503 while its index loads: a signal then stops the
 * server like any other, and the exit does not wait for the load.
 */
class Shutdown {
  public:
    static constexpr uint64_t kDrainMs = 3000;
    static constexpr uint64_t kExitMs = 3000;

    bool stopping() const { return stopping_.load(std::memory_order_acquire); }

    // a /traverse request's attempt, for the duration of its handler
    void enter(const std::shared_ptr<Attempt> &attempt) {
        std::lock_guard<std::mutex> lock(mutex_);
        in_flight_.insert(attempt);
        if (stopping())
            attempt->request_stop(graph::traversal::ExternalStop::CLIENT_GONE);
    }
    void leave(const std::shared_ptr<Attempt> &attempt) {
        std::lock_guard<std::mutex> lock(mutex_);
        in_flight_.erase(attempt);
        if (in_flight_.empty())
            changed_.notify_all();
    }
    // the server accepts connections: from now on a shutdown stops it (stopping it before
    // start() would be lost: start() opens the acceptor anew)
    void started(HttpServer *server) {
        std::lock_guard<std::mutex> lock(mutex_);
        server_ = server;
    }
    // run_server is exiting (the server stopped): the shutdown thread need not force the exit
    void finished() {
        std::lock_guard<std::mutex> lock(mutex_);
        finished_ = true;
        changed_.notify_all();
    }

    // installs the handlers and starts the thread that waits for them (once, at the start of
    // run_server, before any request can arrive)
    void install() {
        if (::pipe(pipe_) != 0) {
            logger->warn("[Server] No shutdown pipe ({}): SIGTERM and SIGINT keep their "
                         "default action", std::strerror(errno));
            return;
        }
        write_end_.store(pipe_[1]);
        struct sigaction action;
        std::memset(&action, 0, sizeof(action));
        action.sa_handler = &Shutdown::on_signal;
        sigemptyset(&action.sa_mask);
        // the threads a signal interrupts carry on (reads, waits) as if it had not come
        action.sa_flags = SA_RESTART;
        sigaction(SIGTERM, &action, nullptr);
        sigaction(SIGINT, &action, nullptr);
        std::thread([this]() { wait_and_stop(); }).detach();
    }

  private:
    static void on_signal(int sig) {
        if (signalled_.exchange(true)) {
            // a second signal: the first one's shutdown is not awaited
            _exit(128 + sig);
        }
        const unsigned char byte = static_cast<unsigned char>(sig);
        const ssize_t written = ::write(write_end_.load(), &byte, 1);
        (void)written;
    }

    void wait_and_stop() {
        unsigned char sig = 0;
        while (::read(pipe_[0], &sig, 1) < 0 && errno == EINTR) {}
        std::unique_lock<std::mutex> lock(mutex_);
        stopping_.store(true, std::memory_order_release);
        if (!server_) {
            logger->info("[Server] Signal {} before the server accepted connections: exiting", sig);
            logger->flush();
            std::_Exit(0);
        }
        logger->info("[Server] Signal {}: shutting down; {} traversal(s) in flight stopped as "
                     "if their clients had gone (nothing written for them), the server stops "
                     "within {} ms", sig, in_flight_.size(), kDrainMs);
        for (const auto &attempt : in_flight_) {
            attempt->request_stop(graph::traversal::ExternalStop::CLIENT_GONE);
        }
        changed_.wait_for(lock, std::chrono::milliseconds(kDrainMs),
                          [&]() { return in_flight_.empty(); });
        HttpServer *server = server_;
        const size_t left = in_flight_.size();
        lock.unlock();
        if (left) {
            logger->warn("[Server] {} traversal(s) still running: their connections are closed",
                         left);
        }
        server->stop();
        lock.lock();
        if (!changed_.wait_for(lock, std::chrono::milliseconds(kExitMs),
                               [&]() { return finished_; })) {
            logger->warn("[Server] A request handler did not return within {} ms of the stop: "
                         "exiting without it", kExitMs);
            logger->flush();
            std::_Exit(0);
        }
    }

    static inline std::atomic<bool> signalled_ { false };
    static inline std::atomic<int> write_end_ { -1 };
    int pipe_[2] = { -1, -1 };
    std::atomic<bool> stopping_ { false };
    std::mutex mutex_;
    std::condition_variable changed_;
    std::set<std::shared_ptr<Attempt>> in_flight_;
    HttpServer *server_ = nullptr;
    bool finished_ = false;
};

// a /traverse handler's attempt is in flight while this lives
class InFlight {
  public:
    InFlight(Shutdown &shutdown, std::shared_ptr<Attempt> attempt)
          : shutdown_(shutdown), attempt_(std::move(attempt)) { shutdown_.enter(attempt_); }
    ~InFlight() { shutdown_.leave(attempt_); }
    InFlight(const InFlight&) = delete;
    InFlight& operator=(const InFlight&) = delete;

  private:
    Shutdown &shutdown_;
    std::shared_ptr<Attempt> attempt_;
};

} // namespace

std::thread start_server(HttpServer &server_startup, Config &config, size_t num_threads,
                         std::function<void()> on_accepting = nullptr) {
    server_startup.config.thread_pool_size = num_threads;

    if (config.host_address != "") {
        server_startup.config.address = config.host_address;
    }
    server_startup.config.port = config.port;
    server_startup.config.timeout_request = 30;    // 30 sec to finish headers
    server_startup.config.timeout_content = kContentTimeoutS;   // 15 minutes for body/compute (per request) max
    if (config.max_request_body_mb) {
        // the library's bound on a request's headers and body (size_t max by default): a
        // Content-Length above it, or a chunked body growing past it, is reported to on_error as
        // message_size and the connection dropped — the 413 the library builds there is never
        // sent (server_http.hpp, read_request_and_content)
        server_startup.config.max_request_streambuf_size
                = static_cast<size_t>(config.max_request_body_mb) << 20;
    }

    logger->info("[Server] Will listen on {} port {}",
                 server_startup.config.address, server_startup.config.port);
    logger->info("[Server] Maximum connections: {}", num_threads);
    if (config.max_request_body_mb)
        logger->info("[Server] Maximum request body: {} MiB", config.max_request_body_mb);
    return std::thread([&server_startup, on_accepting]() {
        server_startup.start([on_accepting](unsigned short) {
            if (on_accepting)
                on_accepting();
        });
    });
}

using GraphPair = std::pair<std::string, std::string>;
using GraphIndexes = tsl::hopscotch_map<std::string, std::vector<GraphPair>>;

/**
 * Address exactly one physical (graph, annotation) pair for a traversal request or the probe.
 * A traversal must stay inside one graph, so a name covering several shards is rejected
 * unless the request disambiguates it with |graph_path|; a graph_path naming one graph with
 * several annotations under the name is rejected too — a traversal reads ONE annotation, and
 * answering from the first of them would leave the others unread without saying so.
 */
const GraphPair& select_traverse_pair(const std::string &name, const std::string *graph_path,
                                      const GraphIndexes &indexes) {
    auto it = indexes.find(name);
    if (it == indexes.end())
        throw InvalidRequest("Bad request: unknown graph '" + name + "'");

    // the same pair listed twice under the name is one index (/search keeps its list as is)
    std::vector<const GraphPair*> pairs;
    for (const auto &pair : it->second) {
        if (std::none_of(pairs.begin(), pairs.end(),
                         [&](const GraphPair *p) { return *p == pair; })) {
            pairs.push_back(&pair);
        }
    }
    if (pairs.size() == 1)
        return *pairs[0];

    // a name whose pairs share one graph has several annotations, which graph_path cannot
    // choose between
    auto several_annotations = [&](const std::string &graph) {
        std::string annotations;
        for (const GraphPair *pair : pairs) {
            if (pair->first == graph)
                annotations += (annotations.empty() ? "" : ", ") + pair->second;
        }
        return InvalidRequest("Bad request: " + (graph_path
                                  ? "graph_path '" + graph + "' names a graph"
                                  : "index '" + name + "' lists one graph (" + graph + ")")
                              + " with several annotations" + (graph_path ? " under index '"
                                  + name + "'" : std::string()) + " (" + annotations
                              + "): a traversal reads one annotation, which graph_path cannot "
                              "choose; list each annotation under a name of its own");
    };
    if (graph_path) {
        const GraphPair *found = nullptr;
        size_t matching = 0;
        for (const GraphPair *pair : pairs) {
            if (pair->first != *graph_path)
                continue;
            matching++;
            if (!found)
                found = pair;
        }
        if (!found) {
            throw InvalidRequest("Bad request: graph_path '" + *graph_path + "' is not part of "
                                 "index '" + name + "'");
        }
        if (matching > 1)
            throw several_annotations(*graph_path);
        return *found;
    }
    std::vector<std::string> graphs;
    for (const GraphPair *pair : pairs) {
        if (std::find(graphs.begin(), graphs.end(), pair->first) == graphs.end())
            graphs.push_back(pair->first);
    }
    if (graphs.size() == 1)
        throw several_annotations(graphs[0]);
    std::string candidates;
    for (const std::string &graph : graphs) {
        candidates += (candidates.empty() ? "" : ", ") + graph;
    }
    throw InvalidRequest("Bad request: index '" + name + "' spans several graphs (" + candidates
                         + "); pass 'graph_path' to pick one");
}

std::vector<std::string> filter_graphs_from_list(
        const GraphIndexes &indexes,
        const Json::Value &content_json,
        size_t request_id,
        size_t max_names_without_filtering = 10) {
    std::vector<std::string> graphs_to_query;

    if (content_json.isMember("graphs") && content_json["graphs"].isArray()) {
        for (const auto &item : content_json["graphs"]) {
            graphs_to_query.push_back(item.asString());
            if (!indexes.count(graphs_to_query.back()))
                throw std::invalid_argument("Request with an uninitialized graph " + graphs_to_query.back());
        }
        // deduplicate
        std::sort(graphs_to_query.begin(), graphs_to_query.end());
        graphs_to_query.erase(std::unique(graphs_to_query.begin(), graphs_to_query.end()), graphs_to_query.end());
    } else {
        if (indexes.size() > max_names_without_filtering) {
            throw std::invalid_argument(
                fmt::format("Bad request: requests without names (no \"graphs\" field) are "
                            "only supported for small indexes (<={} names)",
                            max_names_without_filtering));
        }
        // query all graphs from list `config->fnames`
        for (const auto &[name, _] : indexes) {
            graphs_to_query.push_back(name);
        }
    }
    return graphs_to_query;
}


namespace {

// A request abandoned while it waited for the memory to load its index (`in_ram`): its client
// is gone, or the server stops; nothing was loaded, nothing is written
class LoadAbandoned : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

/**
 * The index one request of a traversal or pattern route reads on a multi-graph server: the
 * pair's resident index, or (`in_ram` on a server on mmap, the pair within --mem-cap-gb) its
 * own copy loaded into RAM for it under a reservation of the server's LoadReservations (the
 * rule and the wait of /search), freed and released when the lease ends. |load_ms|: what the
 * request spent before its work could begin, waiting for the memory and loading (0 without a
 * load); the request's budgets start after it.
 */
class IndexLease {
  public:
    explicit IndexLease(const AnnotatedDBG &resident, InRamPlan plan)
          : index_(&resident), plan_(plan) {}
    IndexLease(std::unique_ptr<AnnotatedDBG> loaded, LoadReservations *reservations,
               size_t bytes, double load_ms)
          : loaded_(std::move(loaded)), index_(loaded_.get()), plan_(InRamPlan::LOAD),
            reservations_(reservations), bytes_(bytes), load_ms_(load_ms) {}
    ~IndexLease() {
        loaded_.reset();
        if (reservations_)
            reservations_->release(bytes_);
    }
    IndexLease(const IndexLease&) = delete;
    IndexLease& operator=(const IndexLease&) = delete;

    const AnnotatedDBG& index() const { return *index_; }
    InRamPlan plan() const { return plan_; }
    double load_ms() const { return load_ms_; }

  private:
    std::unique_ptr<AnnotatedDBG> loaded_;
    const AnnotatedDBG *index_;
    InRamPlan plan_;
    LoadReservations *reservations_ = nullptr;
    size_t bytes_ = 0;
    double load_ms_ = 0;
};

// The parameters of a GET probe of one pair on a multi-graph server (?graph=<name>
// [&graph_path=<path>]), refused as GET /traverse/capabilities refuses them (|route| names the
// route in the messages)
std::pair<std::string, std::optional<std::string>>
probe_parameters(const SimpleWeb::CaseInsensitiveMultimap &query, const std::string &route) {
    std::optional<std::string> name, graph_path;
    for (const auto &[key, value] : query) {
        std::optional<std::string> *slot = key == "graph" ? &name
                                         : key == "graph_path" ? &graph_path : nullptr;
        if (!slot) {
            throw InvalidRequest("Bad request: unknown parameter '" + key + "' of GET "
                                 + route + " (graph, graph_path)");
        }
        if (*slot)
            throw InvalidRequest("Bad request: the parameter '" + key + "' is given twice");
        *slot = value;
    }
    if (!name) {
        throw InvalidRequest("Bad request: in multi-graph mode GET " + route + " "
                             "needs ?graph=<name> (and graph_path=<path> when the name "
                             "spans several graphs); GET /capabilities lists the graphs");
    }
    return { *name, graph_path };
}

} // namespace

// The reverse index of the sequence headers (CoordToHeader::find_header), built while the
// index loads rather than by the first request naming a header: on refseq33m (33M headers)
// that request would spend seconds in it (9.3 s), inside its time budget and outside every
// read timer
static void build_header_index(const graph::AnnotatedDBG &index) {
    const auto *coord_to_header = index.get_coord_to_header();
    if (!coord_to_header)
        return;
    Timer timer;
    const size_t headers = coord_to_header->build_header_index();
    logger->info("[Server] Sequence header index built: {} headers in {:.1f} s", headers,
                 timer.elapsed());
}

// The pairs of a multi-graph server in the order its answers list them: the names in byte
// order, each name's pairs in the graph list's order, a pair listed twice under a name once
static std::vector<std::pair<std::string, GraphPair>>
ordered_pairs(const GraphIndexes &indexes, const std::vector<std::string> &names) {
    std::vector<std::pair<std::string, GraphPair>> pairs;
    for (const std::string &name : names) {
        const size_t first = pairs.size();
        for (const GraphPair &pair : indexes.at(name)) {
            if (std::none_of(pairs.begin() + first, pairs.end(),
                             [&](const auto &p) { return p.second == pair; })) {
                pairs.emplace_back(name, pair);
            }
        }
    }
    return pairs;
}

static std::vector<std::string> sorted_names(const GraphIndexes &indexes) {
    std::vector<std::string> names;
    for (const auto &[name, _] : indexes) {
        names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    return names;
}

/**
 * GET /capabilities' graph_summary of a multi-graph server, computed once at start-up, so that
 * a service learns every pair with one probe: per (name, pair), in ordered_pairs' order, the
 * pair (graph, graph_path, annotation_path), its identity (index_ns, index_fp), k, and what the
 * pattern search makes of it — graph_mode, available and unavailable_reason, mask and counting,
 * as its GET /pattern/capabilities?graph= block states them (the dummy fraction of a pair
 * without its mask, sampled at load, is its block's only) — and what a traversal reads
 * (regime, num_labels, has_coordinates, has_coord_to_header, supports_trace, as every
 * /traverse response states them); and whether the pairs' columns are disjoint
 * (column_overlap), which licenses summing a label's counts and occurrences over the pairs
 */
static Json::Value multi_graph_summary(
        const GraphIndexes &indexes,
        const VectorMap<GraphPair, std::unique_ptr<AnnotatedDBG>> &graphs_cache,
        const std::map<GraphPair, IndexIdentity> &identities) {
    Timer timer;
    Json::Value summary;
    Json::Value pairs(Json::arrayValue);
    uint64_t masked = 0;
    for (const auto &[name, pair] : ordered_pairs(indexes, sorted_names(indexes))) {
        const AnnotatedDBG &index = *graphs_cache.at(pair);
        const DeBruijnGraph &graph = index.get_graph();
        const IndexIdentity &id = identities.at(pair);
        Json::Value p;
        p["graph"] = name;
        p["graph_path"] = pair.first;
        p["annotation_path"] = pair.second;
        p["index_ns"] = string_or_null(id.name);
        p["index_fp"] = string_or_null(id.fp);
        p["k"] = uint_json(graph.get_k());
        // the pattern search, as the pair's block states it (pattern_capabilities_json)
        const graph::pattern::GraphSupport support = route_support(graph);
        const bool recognised = support.supported || support.reason == "alphabet_unsupported"
                                    || support.reason == "alphabet_untested"
                                    || support.reason == "mask_invalid";
        p["available"] = support.supported;
        p["unavailable_reason"] = support.supported ? Json::Value()
                                                    : Json::Value(support.reason);
        p["graph_mode"] = recognised ? Json::Value(to_string(support.mode)) : Json::Value();
        p["mask"] = !recognised ? Json::Value()
                  : !support.mask_present ? Json::Value("absent")
                  : mask_built_at_load(graph) ? Json::Value("built_at_load")
                                              : Json::Value("file");
        p["counting"] = !support.supported ? Json::Value()
                      : support.mask_present ? Json::Value("exact")
                                             : Json::Value("upper_bound");
        masked += recognised && support.mask_present;
        // what a traversal reads, as the capabilities of every /traverse response state it
        const graph::traversal::LabelOracle oracle(index);
        Json::Value t;
        t["regime"] = to_string(oracle.regime());
        t["num_labels"] = uint_json(oracle.num_columns());
        t["has_coordinates"] = oracle.has_coordinates();
        t["has_coord_to_header"] = oracle.coord_to_header() != nullptr;
        t["supports_trace"] = oracle.has_coordinates()
                && oracle.regime() == graph::traversal::Regime::BASIC;
        p["traversal"] = std::move(t);
        pairs.append(std::move(p));
    }
    // the columns of every pair (each pair once, whatever names list it)
    std::vector<std::vector<std::string>> columns;
    for (const auto &[pair, index] : graphs_cache) {
        columns.push_back(index->get_annotator().get_label_encoder().get_labels());
    }
    const ColumnOverlap overlap = column_overlap(columns);
    summary["columns_disjoint"] = overlap.disjoint;
    summary["shared_columns"] = uint_json(overlap.shared);
    summary["pairs"] = std::move(pairs);
    if (overlap.disjoint) {
        logger->info("[Server] The {} columns of the {} (graph, annotation) pairs are disjoint: "
                     "a label's counts can be summed over them", overlap.columns,
                     graphs_cache.size());
    } else {
        logger->info("[Server] {} of the {} column names are columns of more than one (graph, "
                     "annotation) pair (one of them: '{}'): columns_disjoint false, a label's "
                     "counts must not be summed over the pairs", overlap.shared,
                     overlap.columns, overlap.example);
    }
    logger->info("[Server] Pattern search: {} of the {} pairs with a dummy-edge mask (exact "
                 "counts); the others count upper bounds with estimates, their dummy fraction "
                 "sampled at load. Summary computed in {:.3f} s", masked,
                 summary["pairs"].size(), timer.elapsed());
    return summary;
}

int run_server(Config *config) {
    assert(config);
    // SIGTERM and SIGINT shut the server down promptly (class Shutdown); allocated for the
    // process's lifetime, which its detached thread shares
    Shutdown &shutdown = *new Shutdown();
    shutdown.install();
    std::atomic<size_t> num_requests = 0;

    ThreadPool graph_loader(1, 1);
    std::shared_future<std::unique_ptr<AnnotatedDBG>> anno_graph;

    GraphIndexes indexes;

    ThreadPool graphs_pool(get_num_threads(), 1000 /* max_num_tasks */);
    size_t num_server_threads = std::max(1u, get_num_threads());
    set_num_threads(std::max(1u, config->parallel_each));
    config->parallel_each = 1;  // query one batch at a time
    logger->info("[Server] Threads per graph: {}", get_num_threads());

    VectorMap<std::pair<std::string, std::string>, std::unique_ptr<AnnotatedDBG>> graphs_cache;
    bool loaded_with_mmap = utils::with_mmap();

    // The index identity stated by /traverse and /resolve (DESIGN-traverse-graphlet.md
    // §3.1). Single-index mode: --index-name and --index-manifest, checked against the
    // loaded files before the index is served; the meta fingerprint once the index is
    // loaded. Multi-index mode: per (graph, annotation) pair, the name and the manifest of the
    // graph list's optional columns, checked before loading, and the meta fingerprint computed
    // on first use (a pair without them states null).
    IndexIdentity single_identity;
    // the graph's derived data the loader reads (index_derived_files: its mask and Bloom
    // filter), named beside a checked manifest so that an operator sees what is loaded
    // although index_fp does not cover it
    auto log_derived_data = [&](const std::string &graph) {
        std::vector<std::string> loaded;
        for (const IndexDerivedFile &f : index_derived_files(graph)) {
            if (f.loaded)
                loaded.push_back(f.path);
        }
        if (!loaded.empty()) {
            logger->info("[Server] Derived data of {} loaded beside it, not part of index_fp: {}",
                         graph, fmt::join(loaded, ", "));
        }
    };
    // per (graph, annotation) pair: its name and manifest digest from the graph list, its meta
    // fingerprint computed on first use, from the resident index or a copy of it loaded for a
    // request (`in_ram`: the same files, the same identity)
    std::mutex identities_mutex;
    std::map<GraphPair, IndexIdentity> pair_identities;
    std::map<const AnnotatedDBG*, GraphPair> resident_pairs;
    auto identity_of_pair = [&](const GraphPair &pair, const AnnotatedDBG &index) {
        std::lock_guard<std::mutex> lock(identities_mutex);
        IndexIdentity &id = pair_identities[pair];
        if (id.meta_fp.empty())
            id.meta_fp = index_meta_fingerprint(graph::traversal::LabelOracle(index));
        return id;
    };
    auto identity_of = [&](const AnnotatedDBG &index) -> IndexIdentity {
        if (config->fnames.empty())
            return single_identity;
        return identity_of_pair(resident_pairs.at(&index), index);
    };
    // the pair's index_fp alone (the graph list's manifest digest, known at start-up; "" without
    // a manifest): what the envelope entry of a /pattern pair states when its index could not
    // be leased for the request (the meta fingerprint needs the index, index_fp does not)
    auto fp_of_pair = [&](const GraphPair &pair) {
        std::lock_guard<std::mutex> lock(identities_mutex);
        return pair_identities[pair].fp;
    };
    // the multi-graph server's per-pair summary of GET /capabilities and whether its pairs'
    // columns are disjoint, computed once its indexes are loaded
    Json::Value graph_summary;

    if (config->infbase_annotators.size() == 1) {
        assert(config->fnames.empty());
        // a manifest of another bundle is refused before hours of loading, not after
        try {
            single_identity.name = config->index_name;
            if (!config->index_manifest.empty()) {
                single_identity.fp = index_manifest_fingerprint(
                        config->index_manifest,
                        index_bundle_files(config->infbase, config->infbase_annotators[0],
                                           !config->no_coord_mapping),
                        nullptr,
                        index_unloaded_optional_files(config->infbase,
                                                      config->infbase_annotators[0],
                                                      !config->no_coord_mapping));
                logger->info("[Server] Index manifest {}: index_fp {}", config->index_manifest,
                             single_identity.fp);
                log_derived_data(config->infbase);
            }
        } catch (const std::exception &e) {
            logger->error("[Server] {}", e.what());
            std::exit(1);
        }
        anno_graph = graph_loader.enqueue([&]() {
            logger->info("[Server] Loading graph and annotation in parallel...");
            auto anno_graph = initialize_annotated_dbg(*config);
            logger->info("[Server] Annotated graph loaded. Current mem usage: {} MiB", get_curr_RSS() >> 20);
            // set before the future is ready, so every request sees it complete
            single_identity.meta_fp = index_meta_fingerprint(graph::traversal::LabelOracle(*anno_graph));
            build_header_index(*anno_graph);
            return anno_graph;
        });
    } else {
        assert(config->fnames.size() == 1);

        std::ifstream file(config->fnames[0]);

        if (!file.is_open()) {
            logger->error("[Server] Could not open file {} for reading", config->fnames[0]);
            std::exit(1);
        }

        size_t num_indexes = 0;
        std::string line;
        std::vector<GraphListEntry> entries;
        for (size_t line_no = 1; std::getline(file, line); ++line_no) {
            if (line.empty())
                continue; // skip empty lines
            try {
                entries.push_back(parse_graph_list_line(line, line_no));
            } catch (const std::exception &e) {
                logger->error("[Server] Invalid line in the csv file: {}", e.what());
                std::exit(1);
            }
            const GraphListEntry &e = entries.back();
            indexes[e.name].emplace_back(e.graph_path, e.annotation_path);
            num_indexes++;
        }
        // The per-graph identity (DESIGN-traverse-graphlet.md §16.1): every listed manifest is
        // checked against the files its pair loads — sizes and the digest of its list, no
        // re-hashing, as --index-manifest — before hours of loading, and a mismatch, or two
        // lines stating different identities for one pair, refuses to start
        std::map<GraphPair, std::pair<std::string, std::string>> pair_identity;
        try {
            // what each pair loads (the loader dependency inventory of its listed spellings)
            auto inventory = [&](const GraphListEntry &e) {
                return index_bundle_files(e.graph_path, e.annotation_path,
                                          !config->no_coord_mapping);
            };
            pair_identity = graph_list_identities(entries, [&](const GraphListEntry &e) {
                std::string stated;
                std::string fp;
                try {
                    fp = index_manifest_fingerprint(
                            e.manifest_path, inventory(e), &stated,
                            index_unloaded_optional_files(e.graph_path, e.annotation_path,
                                                          !config->no_coord_mapping));
                } catch (const std::exception &ex) {
                    throw std::invalid_argument("line " + std::to_string(e.line) + " of the "
                                                "graph list: " + ex.what());
                }
                if (!stated.empty() && stated != e.index_ns) {
                    logger->warn("[Server] Index manifest {} names the index '{}'; line {} of "
                                 "the graph list names it '{}', which its responses state",
                                 e.manifest_path, stated, e.line, e.index_ns);
                }
                return fp;
            }, inventory);
        } catch (const std::exception &e) {
            logger->error("[Server] {}", e.what());
            std::exit(1);
        }
        std::set<std::string> derived_logged;
        for (const auto &[pair, id] : pair_identity) {
            if (!id.first.empty() || !id.second.empty()) {
                logger->info("[Server] Index ({}, {}): index_ns {}, index_fp {}", pair.first,
                             pair.second, id.first.empty() ? "null" : id.first,
                             id.second.empty() ? "null (no manifest)" : id.second);
            }
            if (!id.second.empty() && derived_logged.insert(pair.first).second)
                log_derived_data(pair.first);
        }
        std::vector<std::string> names;
        for (const auto &[name, _] : indexes) {
            names.push_back(name);
        }
        logger->info("[Server] Loaded a list of {} graphs for {} names: {}",
                     num_indexes, indexes.size(), fmt::join(names, ", "));
        if (loaded_with_mmap) {
            logger->info("[Server] Graphs will be loaded with mmap (--mmap set)."
                         " Make sure they're on a fast drive.");
        } else {
            logger->info("[Server] Graphs will be loaded into RAM (--mmap not set)."
                         " Total memory ≈ sum of all graph+annotation file sizes.");
        }

        // Deduplicate graph paths so the same underlying graph file is loaded only
        // once, even if it appears in multiple (graph, annotation) entries.
        std::vector<std::string> unique_graph_paths;
        tsl::hopscotch_map<std::string, size_t> graph_path_to_idx;
        for (const auto &[name, graphs] : indexes) {
            for (const auto &[graph_fname, anno_fname] : graphs) {
                if (graph_path_to_idx.emplace(graph_fname, unique_graph_paths.size()).second)
                    unique_graph_paths.push_back(graph_fname);
                graphs_cache[{ graph_fname, anno_fname }] = nullptr;
            }
        }
        if (graphs_cache.empty()) {
            logger->error("[Server] No graphs to serve. Exiting.");
            exit(1);
        }
        if (unique_graph_paths.size() < graphs_cache.size()) {
            logger->info("[Server] Deduplicated {} graph references down to {} unique graph files.",
                         graphs_cache.size(), unique_graph_paths.size());
        }

        logger->info("[Server] Loading {} unique graph(s)...", unique_graph_paths.size());
        std::vector<std::shared_ptr<DeBruijnGraph>> loaded_graphs(unique_graph_paths.size());
        // each graph prepared for the pattern search in its loading thread, as a single-graph
        // server prepares its graph: its mask checked on a sample of its W = $ edges
        // (sub-second), and the dummy fraction of a graph without a mask sampled (10,000
        // entries), so that no request pays for it inside its deadline
        PatternPreparation list_preparation;
        list_preparation.check_mask = true;
        list_preparation.sample_fraction = true;
        list_preparation.progress = false;
        #pragma omp parallel for num_threads(get_num_threads() * num_server_threads) schedule(dynamic)
        for (size_t i = 0; i < unique_graph_paths.size(); ++i) {
            loaded_graphs[i] = load_critical_dbg(unique_graph_paths[i]);
            prepare_pattern_graph(loaded_graphs[i], unique_graph_paths[i], list_preparation);
        }

        logger->info("[Server] Loading {} annotation(s)...", graphs_cache.size());
        #pragma omp parallel for num_threads(get_num_threads() * num_server_threads) schedule(dynamic)
        for (size_t i = 0; i < graphs_cache.size(); ++i) {
            Config config_copy = *config;
            auto it = graphs_cache.nth(i);
            const auto &[graph_fname, anno_fname] = it.key();
            config_copy.infbase = graph_fname;
            config_copy.infbase_annotators = { anno_fname };
            it.value() = initialize_annotated_dbg(loaded_graphs[graph_path_to_idx.at(graph_fname)],
                                                  config_copy);
        }
        for (const auto &[pair, index] : graphs_cache) {
            build_header_index(*index);
        }
        // every response computed on a pair states its identity
        for (const auto &[pair, index] : graphs_cache) {
            const auto &[name, fp] = pair_identity.at(pair);
            IndexIdentity &id = pair_identities[pair];
            id.name = name;
            id.fp = fp;
            resident_pairs[index.get()] = pair;
        }
        graph_summary = multi_graph_summary(indexes, graphs_cache, pair_identities);
        logger->info("[Server] All graphs were loaded ({}). Ready to serve queries.",
                     loaded_with_mmap ? "with mmap" : "into RAM");
        // Dynamic per-request loads (in_ram path) should always be in RAM.
        utils::set_mmap(false);
    }

    size_t memory_all = config->memory_available * 1e9;
    // the memory lent to per-request loads (`in_ram`, every route that reads an index): a load
    // waits until its files' size is free (--mem-cap-gb)
    LoadReservations reservations(memory_all);
    std::atomic<size_t> graphs_being_queried = 0;

    // the actual server
    HttpServer server;
    server.resource["^/search"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                              shared_ptr<HttpServer::Request> request) {
        size_t request_id = num_requests++;
        process_request(response, request, request_id, [&](const std::string& content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError();  // the index is not loaded yet, so we can't process the request

            Json::Value content_json = parse_json_string(content);
            logger->info("[Server] Request {}: {}", request_id, content_json.toStyledString());
            Json::Value result;

            // simple case with a single graph pair
            if (!config->fnames.size()) {
                if (content_json.isMember("graphs"))
                    throw std::invalid_argument("Bad request: no support for filtering graphs on this server");
                logger->trace("Request {}: Started querying graph {}, in total graphs being queried at the moment: {}",
                              request_id, config->infbase, graphs_being_queried.fetch_add(1) + 1);
                try {
                    result = process_search_request(content_json, *anno_graph.get(), *config);
                } catch (...) {
                    graphs_being_queried--;
                    throw;
                }
                graphs_being_queried--;
            } else {
                std::vector<std::string> graphs_to_query
                        = filter_graphs_from_list(indexes, content_json, request_id,
                                                  config->max_graphs_without_selection);
                std::mutex mu;
                std::vector<std::shared_future<std::exception_ptr>> futures;
                for (const auto &name : graphs_to_query) {
                    for (const auto &[graph_fname, anno_fname] : indexes[name]) {
                        futures.push_back(graphs_pool.enqueue([&,config,graph_fname=graph_fname,anno_fname=anno_fname]() {
                            logger->trace("Request {}: Started querying graph {}. In total graphs being queried at the moment: {}",
                                          request_id, graph_fname, graphs_being_queried.fetch_add(1) + 1);
                            size_t index_size_reserved = 0;
                            auto release_memory = [&]() {
                                if (!index_size_reserved)
                                    return;
                                reservations.release(index_size_reserved);
                                index_size_reserved = 0;
                            };
                            try {
                                std::unique_ptr<AnnotatedDBG> index_loaded;
                                const AnnotatedDBG *index;
                                bool in_ram = content_json.isMember("in_ram") && content_json["in_ram"].asBool();
                                if (in_ram && !loaded_with_mmap)
                                    in_ram = false;  // already in RAM, no need to re-load
                                size_t index_size = in_ram ? std::filesystem::file_size(graph_fname)
                                                                + std::filesystem::file_size(anno_fname)
                                                           : -1;
                                if (in_ram && index_size > memory_all) {
                                    logger->warn("Request {}: Graph of size {} GB is too large to fit into "
                                                 "RAM (reserved memory: {} GB). It will be queried with mmap",
                                                 request_id, index_size / 1e9, memory_all / 1e9);
                                    in_ram = false;
                                }
                                Timer timer;
                                if (in_ram) {
                                    reservations.reserve(index_size);
                                    index_size_reserved = index_size;
                                    logger->trace("Request {}: Loading graph {} of size {} GB to RAM...",
                                                  request_id, graph_fname, index_size / 1e9);
                                    timer.reset();
                                    Config config_copy = *config;
                                    config_copy.infbase = graph_fname;
                                    config_copy.infbase_annotators = { anno_fname };
                                    index_loaded = initialize_annotated_dbg(config_copy);
                                    index = index_loaded.get();
                                } else {
                                    index = graphs_cache.at({ graph_fname, anno_fname }).get();
                                }

                                auto json = process_search_request(content_json, *index, *config);
                                logger->trace("Request {}: {} graph {} {} in {} sec",
                                              request_id, in_ram ? "Loaded and searched" : "Searched",
                                              graph_fname, anno_fname, timer.elapsed());

                                index_loaded.reset();
                                release_memory();

                                std::lock_guard<std::mutex> lock(mu);
                                if (result.empty()) {
                                    result = std::move(json);
                                } else {
                                    assert(json.size() == result.size());
                                    for (Json::ArrayIndex i = 0; i < result.size(); ++i) {
                                        if (result[i][SeqSearchResult::SEQ_DESCRIPTION_JSON_FIELD]
                                                != json[i][SeqSearchResult::SEQ_DESCRIPTION_JSON_FIELD]) {
                                            throw std::logic_error("ERROR: Results for different sequences can't be merged");
                                        }
                                        for (auto&& value : json[i]["results"]) {
                                            result[i]["results"].append(std::move(value));
                                        }
                                    }
                                }
                            } catch (...) {
                                release_memory();
                                graphs_being_queried--;
                                return std::current_exception();
                            }
                            graphs_being_queried--;
                            return std::exception_ptr();
                        }));
                    }
                }
                for (auto &future : futures) {
                    if (auto ex = future.get())
                        std::rethrow_exception(ex);
                }
            }
            return result;
        });
    };

    server.resource["^/align"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                             shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string &content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError(); // the index is not loaded yet, so we can't process the request

            if (!config->fnames.size())
                return process_align_request(content, anno_graph.get()->get_graph(), *config);

            throw std::invalid_argument("Bad request: alignment requests are not yet supported for "
                                        "servers with multiple graphs");
        });
    };

    // Whether the client of |request| is gone or the server stops (a shutdown is treated as a
    // client gone, class Shutdown)
    auto gone_of = [&shutdown](const shared_ptr<HttpServer::Request> &request) {
        return [r = request.get(), &shutdown]() { return shutdown.stopping() || client_gone(*r); };
    };

    // The index of a /resolve, /traverse or /pattern request on |pair| of a multi-graph server:
    // the resident one, or with `in_ram` (in_ram_plan, /search's rule: the server on mmap, the
    // pair's files within --mem-cap-gb) a copy loaded into RAM for the request once its memory
    // is reserved (/search's wait, abandoned when |gone| answers true: LoadAbandoned), prepared
    // as the resident one was — for /pattern its mask checked and its dummy fraction sampled —
    // so that the request's work, and its budgets, start on a ready index after the load. The
    // copy's record mapping (the sequence headers and their reverse index, built at start-up)
    // is the resident pair's own object: the copy is of the same files, the same .seqs, so a
    // mapping loaded and indexed again would be identical, for seconds and gigabytes per
    // request on a chunk with millions of records
    auto lease_index = [&](const GraphPair &pair, bool in_ram, size_t request_id,
                           const std::function<bool()> &gone, bool for_pattern) {
        const AnnotatedDBG &resident = *graphs_cache.at(pair);
        const size_t bytes = in_ram && loaded_with_mmap
                ? std::filesystem::file_size(pair.first) + std::filesystem::file_size(pair.second)
                : 0;
        const InRamPlan plan = in_ram_plan(in_ram, loaded_with_mmap, bytes, memory_all);
        if (plan == InRamPlan::RESIDENT_TOO_LARGE) {
            logger->warn("Request {}: Graph of size {} GB is too large to fit into RAM (reserved "
                         "memory: {} GB). It will be queried with mmap", request_id,
                         bytes / 1e9, memory_all / 1e9);
        }
        if (plan != InRamPlan::LOAD)
            return std::make_unique<IndexLease>(resident, plan);
        Timer timer;
        if (!reservations.reserve(bytes, gone)) {
            throw LoadAbandoned("the client closed its connection (or the server stops) while "
                                "the request waited for the memory to load its index");
        }
        std::unique_ptr<AnnotatedDBG> loaded;
        try {
            logger->trace("Request {}: Loading graph {} of size {} GB to RAM...", request_id,
                          pair.first, bytes / 1e9);
            Config config_copy = *config;
            config_copy.infbase = pair.first;
            config_copy.infbase_annotators = { pair.second };
            PatternPreparation preparation;
            preparation.check_mask = for_pattern;
            preparation.sample_fraction = for_pattern;
            preparation.progress = false;
            // the resident index of this very pair (graphs_cache's key is the pair, its
            // annotation path the one loaded here), so its mapping is this annotation's .seqs
            loaded = initialize_annotated_dbg(config_copy, preparation,
                                              resident.share_coord_to_header());
        } catch (...) {
            reservations.release(bytes);
            throw;
        }
        const double load_ms = timer.elapsed() * 1000;
        logger->info("[Server] Request {}: ({}, {}) loaded into RAM for it in {:.0f} ms (the "
                     "wait for its {:.3f} GB included)", request_id, pair.first, pair.second,
                     load_ms, bytes / 1e9);
        return std::make_unique<IndexLease>(std::move(loaded), &reservations, bytes, load_ms);
    };

    /**
     * POST /pattern on a multi-graph server, as /search: the request's `graphs` (every name
     * when there are at most 10 and it names none), each pair of each name answered as one
     * single-graph request would be (its own deadline, caps and memory account, starting after
     * its index is loaded for it with `in_ram`), in parallel on the graphs' pool; the envelope
     * (pattern_envelope) carries one entry per pair in ordered_pairs' order, tagged with its
     * pair (graph, graph_path, annotation_path) and its index_fp: the pair's answer, or the
     * refusal the pair alone would have been answered (SPEC §24.1) -- its graph's support,
     * its annotation, its deadline, a failure while it was processed. What no pair decides
     * (the body, `graphs`, `in_ram`, the request's fields) is checked before any pair and
     * refuses the whole request as before. The envelope is written by the latest deadline of
     * the pairs that answered, past it 503 deadline (|delivery|); a refused pair binds
     * nothing (one refused for its deadline has passed it).
     */
    auto multi_graph_pattern = [&](const std::string &content, size_t request_id,
                                   const std::function<bool()> &gone,
                                   PatternDelivery *delivery) {
        Timer timer;
        Json::Value json = parse_pattern_body(content);
        if (!json.isObject())
            throw PatternRefusal(400, "invalid_request", "request: expected an object");
        std::vector<std::string> names;
        std::optional<bool> in_ram;
        try {
            names = pattern_graph_names(json, sorted_names(indexes),
                                        config->max_graphs_without_selection);
            in_ram = in_ram_field(json);
        } catch (const std::invalid_argument &e) {
            throw PatternRefusal(400, "invalid_request", e.what());
        }
        // each pair's request: the body without the selection, its fields checked once here
        // (the same for every pair), so that their refusal is the request's and not a pair's
        json.removeMember("graphs");
        validate_pattern_request(json, pattern_limits(*config));
        const auto targets = ordered_pairs(indexes, names);
        logger->info("[Server] Request {}: /pattern on {} pair(s) of {} graph name(s){}",
                     request_id, targets.size(), names.size(),
                     in_ram.value_or(false) ? ", in_ram" : "");

        struct Slot {
            Json::Value entry;
            bool answered = false;
            bool aborted = false;
            // when answered: the pair's deadline, from which the envelope's is taken
            graph::pattern::Deadline::Clock::time_point start;
            double budget_ms = 0;
        };
        std::vector<Slot> slots(targets.size());
        std::vector<std::shared_future<void>> futures;
        for (size_t i = 0; i < targets.size(); ++i) {
            futures.push_back(graphs_pool.enqueue([&, i]() {
                Slot &slot = slots[i];
                const std::string &name = targets[i].first;
                const GraphPair &pair = targets[i].second;
                // the pair's entry when it is refused: what it alone would have been answered
                auto refused = [&](int http_status, const Json::Value &body) {
                    slot.entry = pattern_pair_refused(name, pair.first, pair.second,
                                                      fp_of_pair(pair), http_status, body);
                };
                try {
                    if (gone()) {
                        slot.aborted = true;
                        return;
                    }
                    const auto lease = lease_index(pair, in_ram.value_or(false), request_id,
                                                   gone, /* for_pattern */ true);
                    const IndexIdentity identity = identity_of_pair(pair, lease->index());
                    PatternDelivery own;
                    own.set_abort(gone);
                    slot.start = graph::pattern::Deadline::Clock::now();
                    Json::Value answer = process_pattern_request(json, lease->index(), *config,
                                                                 &identity, &own);
                    slot.budget_ms = answer["limits"]["time_budget_ms"].asDouble();
                    if (in_ram)
                        answer["timing"]["load_ms"] = lease->load_ms();
                    slot.entry = pattern_pair_answered(name, pair.first, pair.second,
                                                       identity.fp, std::move(answer));
                    slot.answered = true;
                } catch (const PatternRefusal &e) {
                    logger->warn("[Server] Request {}: pair ({}, {}) of '{}' refused ({} {}): {}",
                                 request_id, pair.first, pair.second, name, e.status(),
                                 e.code(), e.what());
                    refused(e.status(), e.body());
                } catch (const graph::pattern::Aborted &) {
                    slot.aborted = true;
                } catch (const LoadAbandoned &) {
                    slot.aborted = true;
                } catch (const std::exception &e) {
                    // a failure of this pair, as the route answers one for a single graph
                    // (answer_request): 400 with its text, no code
                    logger->warn("[Server] Request {}: pair ({}, {}) of '{}' failed: {}",
                                 request_id, pair.first, pair.second, name, e.what());
                    Json::Value body;
                    body["error"] = e.what();
                    refused(400, body);
                } catch (...) {
                    logger->warn("[Server] Request {}: pair ({}, {}) of '{}' failed",
                                 request_id, pair.first, pair.second, name);
                    Json::Value body;
                    body["error"] = "Internal server error";
                    refused(500, body);
                }
            }));
        }
        for (auto &future : futures) {
            future.wait();
        }
        for (const Slot &slot : slots) {
            if (slot.aborted)
                throw graph::pattern::Aborted();
        }
        // the deadline of the writing: the latest of the answered pairs'
        const Slot *latest = nullptr;
        for (const Slot &slot : slots) {
            if (!slot.answered)
                continue;
            if (!latest || slot.start + std::chrono::duration<double, std::milli>(slot.budget_ms)
                    > latest->start + std::chrono::duration<double, std::milli>(latest->budget_ms)) {
                latest = &slot;
            }
        }
        if (latest) {
            delivery->set_deadline(graph::pattern::Deadline(
                    latest->start, latest->budget_ms,
                    static_cast<double>(config->pattern_finalize_ms)));
        }
        std::vector<Json::Value> entries;
        entries.reserve(slots.size());
        for (Slot &slot : slots) {
            entries.push_back(std::move(slot.entry));
        }
        Json::Value out = pattern_envelope(names, std::move(entries), timer.elapsed() * 1000);
        logger->info("[Server] Request {}: /pattern envelope: {} pair(s) answered, {} refused",
                     request_id, out["answered"].asUInt(), out["refused"].asUInt());
        delivery->check();
        return out;
    };

    // Count, and extract without reading annotation, the graph contexts of short motifs and
    // IUPAC patterns (DESIGN-pattern-search.md; pattern.hpp): on the single graph, or on a
    // multi-graph server on the graphs a request selects (multi_graph_pattern). Every refusal
    // is {"error", "code"}; an answer is written under the request's own deadline (its
    // finalisation reserve), past which it is 503 "deadline", never a partial answer.
    server.resource["^/pattern$"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                                shared_ptr<HttpServer::Request> request) {
        // what the route throws, in the server's terms: a refusal as its HTTP answer, an
        // abandoned request as a client gone
        auto translated = [](const auto &f) {
            try {
                return f();
            } catch (const PatternRefusal &e) {
                throw HttpError(e.status(), e.body());
            } catch (const graph::pattern::Aborted &e) {
                throw ClientGone(e.what());
            }
        };
        // the request's deadline, set once its body is parsed; the writing and the
        // compression of the answer are checked against it
        PatternDelivery delivery;
        // a client that is gone, or a shutdown, is not answered: the work ends at its next
        // clock reading, the writing at its next check, and nothing is written (as /resolve's and
        // /traverse's)
        auto gone = gone_of(request);
        delivery.set_abort(gone);
        ResponseControl control;
        // nor with an error: a refusal, a 400 of a malformed body, a 503 (the index loading,
        // the deadline) or an unexpected failure is not written to a client that left — a
        // half-close counts — or during a shutdown (SPEC §3).
        // Asked apart from the deadline: a 503 at the deadline reaches a client still there
        control.gone = gone;
        control.check = [&]() { translated([&]() { delivery.check(); }); };
        // the time to compress is inside the reserve: the traversal routes' faster level
        control.compression_level = config->traverse_compression_level;
        const size_t request_id = num_requests++;
        process_request(response, request, request_id, [&](const std::string &content) {
            return translated([&]() {
                if (config->fnames.size())
                    return multi_graph_pattern(content, request_id, gone, &delivery);
                if (anno_graph.wait_for(0s) != std::future_status::ready)
                    throw CurrentlyInitializingError();
                const AnnotatedDBG &index = *anno_graph.get();
                const IndexIdentity identity = identity_of(index);
                const Json::Value json = parse_pattern_body(content);
                Json::Value answer = process_pattern_request(json, index, *config, &identity,
                                                             &delivery);
                // in_ram: a single-graph server's index is the one it holds, nothing is loaded
                if (json.isObject() && json.isMember("in_ram"))
                    answer["timing"]["load_ms"] = 0.0;
                return answer;
            });
        }, /* compact */ true, &control);
    };

    // The ledger-managed /traverse attempts of this process (requests with attempt_id): the
    // backend half of the attempt ledger of DESIGN-traverse-graphlet.md §14
    // (traverse_attempts.hpp)
    AttemptSettings attempt_settings;
    attempt_settings.allowance_ms = config->traverse_attempt_allowance_ms;
    attempt_settings.hard_cap_ms = kContentTimeoutS * 1000.0 - 1000;
    attempt_settings.retention_s = config->traverse_attempt_retention_s;
    attempt_settings.retention_count = config->traverse_attempt_retention;
    attempt_settings.tombstone_max_s = config->traverse_attempt_tombstone_max_s;
    attempt_settings.clock_skew_ms = config->traverse_clock_skew_ms;
    attempt_settings.content_timeout_s = kContentTimeoutS;
    attempt_settings.delivery_compress_mbps = config->traverse_delivery_compress_mbps;
    attempt_settings.delivery_build_mbps = config->traverse_delivery_build_mbps;
    // until measured longer, the walk is taken to end at most a chunk of a read and 950 ms (the
    // heads between two readings of the clock, the stopped seed's finalisation) after its
    // walk-until (calibrated: 352-1,001 ms measured on SRA, 1,699 once under load)
    attempt_settings.delivery_stop_ms = static_cast<double>(config->traverse_chunk_target_ms) + 950;
    AttemptRegistry attempts(attempt_settings);
    logger->info("[Server] Traverse attempts: server_instance {}, allowance {} ms, {}",
                 attempts.server_instance(), attempt_settings.allowance_ms,
                 attempts.retention_text());

    // The traversal routes' transport: compact JSON, compressed (when the client accepts it)
    // at a faster zlib level than the other routes' 9 — the time to build a response is what
    // an attempt's delivery window bounds, and level 1 writes 3-4 times faster than 9 at
    // about 1.8 times the bytes (measured on real responses); the decompressed bytes are the
    // same at every level
    ResponseControl traversal_io;
    traversal_io.compression_level = config->traverse_compression_level;

    // Report where a query is supported and which labels carry which blocks, and
    // optionally freeze seeds for /traverse. No graph traversal. A request with
    // bounds.time_budget_ms runs under that deadline (capped by --traverse-max-time-ms, as
    // /traverse's), its answer written under it: past it, 503 "deadline", never a partial
    // answer; without the field no deadline is set or read
    // the cap of the opt-in deadline: --traverse-max-time-ms, and never above what the
    // transport can honour (with the flag 0, uncapped, or above it, a budget past the content
    // timeout would be accepted and the connection closed before any answer or 503); lowered
    // requests are stated in limits.clamped
    const ResolveTimeLimits resolve_time {
        config->traverse_max_time_ms > 0
                && config->traverse_max_time_ms < static_cast<double>(kServerMaxDeadlineMs)
            ? config->traverse_max_time_ms : static_cast<double>(kServerMaxDeadlineMs),
        kResolveFinalizeMs };
    server.resource["^/resolve$"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                                shared_ptr<HttpServer::Request> request) {
        auto as_http = [](const ResolveDeadline &e) {
            return HttpError(503, resolve_deadline_body(e));
        };
        // the request's deadline, set once its body is parsed (none without the field: the
        // check then does nothing, and the text and its compression are the same bytes)
        ResolveDelivery delivery;
        ResponseControl control;
        control.compression_level = traversal_io.compression_level;
        control.check = [&]() {
            try {
                delivery.check();
            } catch (const ResolveDeadline &e) {
                throw as_http(e);
            }
        };
        const size_t request_id = num_requests++;
        process_request(response, request, request_id, [&](const std::string &content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError();

            Json::Value json = parse_json_string(content);
            const bool multi = !config->fnames.empty();
            const GraphSelection selection = traverse_graph_selection(json, multi);
            const std::optional<bool> in_ram = in_ram_field(json);
            std::unique_ptr<IndexLease> lease;
            IndexIdentity identity = single_identity;
            if (multi) {
                const GraphPair &pair = select_traverse_pair(
                        selection.name, selection.graph_path ? &*selection.graph_path : nullptr,
                        indexes);
                // a malformed request is refused before its index is loaded for it
                if (in_ram.value_or(false) && loaded_with_mmap)
                    parse_resolve_request(json);
                try {
                    lease = lease_index(pair, in_ram.value_or(false), request_id,
                                        gone_of(request), /* for_pattern */ false);
                } catch (const LoadAbandoned &e) {
                    throw ClientGone(e.what());
                }
                identity = identity_of_pair(pair, lease->index());
            }
            const AnnotatedDBG &index = lease ? lease->index() : *anno_graph.get();
            // a client that is gone is not answered: abandoned between the request's phases
            try {
                Json::Value out = process_resolve_request(json, index, config->index_release,
                                                          config->resolve_max_query_bp,
                                                          &identity, gone_of(request),
                                                          resolve_time, &delivery);
                // in_ram: the time before the work began, waiting for and loading its index
                if (in_ram)
                    out["timing"]["load_ms"] = lease ? lease->load_ms() : 0.0;
                return out;
            } catch (const graph::traversal::AttemptAborted &e) {
                throw ClientGone(e.what());
            } catch (const ResolveDeadline &e) {
                throw as_http(e);
            }
        }, /* compact */ true, &control);
    };

    // Extend frozen seeds along consistent annotation labels.
    server.resource["^/traverse$"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                                 shared_ptr<HttpServer::Request> request) {
        const size_t request_id = num_requests++;
        // Every traversal is stopped when its client is gone: its walk polls the client's
        // connection, and nothing is written. A request with attempt_id is also an attempt of
        // the service's ledger: registered, cancellable by id, bounded in duration by the
        // server itself, and every response to it states its usage (traverse_attempts.hpp)
        // (the client check holds the request weakly: a finished attempt is retained for its
        // state, not with the request's body)
        auto attempt = std::make_shared<Attempt>(
                request_id, request->header_read_time, attempts.settings(),
                attempts.server_instance(), /* enforced */ true,
                [weak = std::weak_ptr<HttpServer::Request>(request), &shutdown]() {
                    // a shutdown is treated as a client gone (class Shutdown)
                    if (shutdown.stopping())
                        return true;
                    auto r = weak.lock();
                    return !r || client_gone(*r);
                });
        // stopped with the others if the server shuts down while it runs
        const InFlight in_flight(shutdown, attempt);
        // what this server measured of its deliveries replaces the configured rates and ratios
        // in the attempt's reserve (the slowest, the smallest, of its recent responses)
        attempt->set_measured(attempts.measured());
        bool registered = false;
        // what an error states besides its message once the attempt is registered: its usage
        auto with_usage = [&](int status, const std::string &what, const std::string &reason) {
            Json::Value body;
            body["error"] = what;
            body["usage"] = attempt->usage_json(reason);
            return HttpError(status, std::move(body));
        };
        ResponseControl control;
        control.compression_level = config->traverse_compression_level;
        control.on_compressed = [&](size_t text_bytes, double seconds) {
            if (text_bytes >= kMeasuredTextBytes && seconds > 0)
                attempts.note_compress_rate(static_cast<double>(text_bytes) / seconds / 1e6);
        };
        // each seed's result is written as text once built (its tree freed at once, its bytes
        // known to the attempt's delivery reserve); the response is assembled from them, byte
        // for byte the text of the whole tree. The texts are moved into the assembly, which
        // frees each once copied: nothing reads them after it, and kept here they would live
        // until the handler returned, through the compression and the transport's copy
        ResultTexts texts;
        control.write = [&texts, attempt](const Json::Value &envelope,
                                          const std::function<void()> &check) {
            // the longest stretch between two checks (deadline_check)
            double gap = 0;
            std::string text = texts.active
                ? assemble_traverse_response(envelope, std::move(texts.texts), check, &gap)
                : json_text(envelope, true, check, &gap);
            attempt->note_delivery_gap_ms(gap);
            return text;
        };
        // the attempt's client and bound while the response is written and compressed
        control.check = [&]() {
            try {
                attempt->check_delivery();
            } catch (const graph::traversal::AttemptAborted &e) {
                throw ClientGone(e.what());
            } catch (const AttemptAtBound &e) {
                throw with_usage(503, e.what(), "deadline");
            }
        };
        control.on_written = [&](int status, std::optional<size_t> bytes) {
            // the longest read or head piece of any /traverse walk (deadline_check), and the
            // slowest build rate and smallest account per text byte measured on its large seeds
            // (the delivery reserve; the latter without the record coordinates' share, so
            // that an output with them measures what the same output without them would)
            attempts.note_uninterruptible(attempt->max_uninterruptible_ms());
            attempts.note_uninterruptible(attempt->max_delivery_gap_ms());
            attempts.note_build_rate(attempt->own_build_mbps());
            attempts.note_account_per_text_byte(attempt->delivery_detail(),
                                                attempt->own_account_per_text_byte());
            attempts.note_stop_latency(attempt->own_stop_ms());
            if (!registered) {
                if (!status) {
                    logger->info("[Server] Request {}: client gone, {}; no response written",
                                 request_id, attempt->stop_summary());
                }
                return;
            }
            const std::string reason = !status ? "client_gone"
                                     : status == 503 ? "deadline"
                                     : status >= 400 ? "error"
                                     : attempt->walk_reason();
            attempts.finish(attempt, reason, status, bytes);
            logger->info("[Server] Attempt {} (request {}) finished ({}): {}; {}",
                         attempt->ids().attempt_id, request_id, reason, attempt->stop_summary(),
                         status ? fmt::format("response {}, {} bytes, {:.0f} ms", status,
                                              bytes.value_or(0), attempt->elapsed_ms())
                                : std::string("no response written (the client is gone)"));
        };
        process_request(response, request, request_id, [&](const std::string &content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError();

            Json::Value json = parse_json_string(content);
            // a malformed id or not_after_ms is refused before anything is registered (400, no
            // usage)
            attempt->set_ids(attempt_ids(json));
            if (attempt->managed()) {
                // an id does not run twice at once, nor again while it is retained or held, on
                // this server process: a second request with a running, retained, held or
                // tombstoned id is refused, without usage (it would be reconciled against the
                // other attempt); and one whose not_after_ms has passed is not started at all
                // (its ledger may already have released it), also without usage. Once its
                // retention and hold are over the id runs again
                if (auto refused = attempts.start(attempt)) {
                    if (refused->instance_mismatch || refused->expired) {
                        logger->info("[Server] Attempt {} (request {}): not started, {}",
                                     attempt->ids().attempt_id, request_id,
                                     refused->body["error"].asString());
                        throw HttpError(409, std::move(refused->body));
                    }
                    Json::Value body;
                    if (refused->body["tombstone"].asBool()) {
                        // its attempt object states the suppression, judged against this
                        // request's not_after_ms (the tombstone now covers it when it can)
                        body["error"] = "attempt_id '" + attempt->ids().attempt_id + "' was "
                                        "cancelled before this request arrived (tombstoned: "
                                      + attempts.retention_text() + "): it is not run";
                    } else {
                        // a finished id is also refused while held past its retention, and runs
                        // again once neither keeps it
                        body["error"] = "attempt_id '" + attempt->ids().attempt_id + "' is "
                                        "running, or retained or held after it finished, on "
                                        "this server (" + attempts.retention_text()
                                      + "): it is not run while so; its state is in attempt";
                    }
                    body["attempt"] = std::move(refused->body);
                    throw HttpError(409, std::move(body));
                }
                registered = true;
                logger->info("[Server] Attempt {} (request {}): registered",
                             attempt->ids().attempt_id, request_id);
            } else if (const uint64_t now = attempts.now_ms();
                       not_after_passed(attempt->ids(), now)) {
                // not_after_ms without attempt_id is honoured the same way: refused at
                // handler start once passed, nothing run
                Json::Value body = expired_json(attempt->ids(), now, attempts.server_instance());
                logger->info("[Server] Request {}: not started, {}", request_id,
                             body["error"].asString());
                throw HttpError(409, std::move(body));
            }
            try {
                const bool multi = !config->fnames.empty();
                const GraphSelection selection = traverse_graph_selection(json, multi);
                const std::optional<bool> in_ram = in_ram_field(json);
                std::unique_ptr<IndexLease> lease;
                IndexIdentity identity = single_identity;
                if (multi) {
                    const GraphPair &pair = select_traverse_pair(
                            selection.name,
                            selection.graph_path ? &*selection.graph_path : nullptr, indexes);
                    // a malformed request is refused before its index is loaded for it
                    if (in_ram.value_or(false) && loaded_with_mmap)
                        parse_traverse_request(json);
                    // the wait for the memory is abandoned when the client is gone or the server
                    // stops; a cancel takes effect when the work begins (its seeds not started)
                    lease = lease_index(pair, in_ram.value_or(false), request_id,
                                        gone_of(request), /* for_pattern */ false);
                    identity = identity_of_pair(pair, lease->index());
                }
                const AnnotatedDBG &index = lease ? lease->index() : *anno_graph.get();
                TraverseLimits limits;
                limits.max_time_ms = config->traverse_max_time_ms;
                limits.max_seeds = config->traverse_max_seeds;
                limits.max_seed_bp = config->traverse_max_seed_bp;
                limits.max_seed_labels = config->traverse_max_seed_labels;
                limits.max_memory_mb = config->traverse_max_memory_mb;
                limits.max_work_units = config->traverse_max_work_units;
                limits.chunk_target_ms = static_cast<double>(config->traverse_chunk_target_ms);
                limits.path_cache_bytes = config->traverse_path_cache_mb << 20;
                // a multi-graph server's graphs are the chunks of an index, which the fan-out
                // of /search sends each seed to whether or not they hold it: a seed a chunk does
                // not hold is its result, whichever spelling (`graph`, `graphs`) named the chunk
                limits.not_in_graph_per_seed = multi;
                if (in_ram)
                    limits.load_ms = lease ? lease->load_ms() : 0.0;
                return process_traverse_request(json, index, config->index_release, limits,
                                                &identity, attempt.get(), &texts);
            } catch (const LoadAbandoned &e) {
                throw ClientGone(e.what());
            } catch (const graph::traversal::AttemptAborted &e) {
                throw ClientGone(e.what());
            } catch (const AttemptAtBound &e) {
                throw with_usage(503, e.what(), "deadline");
            } catch (const std::exception &e) {
                if (!attempt->managed())
                    throw;
                throw with_usage(400, e.what(), "error");
            } catch (...) {
                if (!attempt->managed())
                    throw;
                throw with_usage(500, "Internal server error", "error");
            }
        }, /* compact */ true, &control);
    };

    // Cancel a running attempt by its id (POST {"attempt_id", "wait_ms"?, "not_after_ms"?}):
    // 200 when it was asked to stop (state stopping, or finished within wait_ms), 404 when it
    // has finished or is unknown (an unknown id tombstoned, held to the not_after_ms given +
    // the clock skew allowance), 429 when an unknown id was not tombstoned. Not refused while
    // the index loads: an attempt may be cancelled at any time
    server.resource["^/traverse/cancel$"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                                        shared_ptr<HttpServer::Request> request) {
        const size_t request_id = num_requests++;
        process_request(response, request, request_id, [&](const std::string &content) {
            Json::Value json = parse_json_string(content);
            if (!json.isObject())
                throw InvalidRequest("request: expected an object");
            for (const std::string &name : json.getMemberNames()) {
                if (name != "attempt_id" && name != "wait_ms" && name != "not_after_ms")
                    throw InvalidRequest("request: unknown field '" + name + "'");
            }
            if (!json["attempt_id"].isString() || !valid_attempt_id(json["attempt_id"].asString())) {
                throw InvalidRequest("request.attempt_id: expected a string of 1 to 128 "
                                     "characters from [A-Za-z0-9._:-]");
            }
            uint64_t wait_ms = 0;
            if (json.isMember("wait_ms")) {
                const Json::Value &w = json["wait_ms"];
                if (!w.isIntegral() || (w.isInt64() && w.asInt64() < 0) || w.asUInt64() > 10'000)
                    throw InvalidRequest("request.wait_ms: expected an integer in [0, 10000]");
                wait_ms = w.asUInt64();
            }
            // the not_after_ms of the request being cancelled (the same rule as /traverse's): a
            // tombstone of an unknown id is held until it has passed + the skew allowance
            const std::optional<uint64_t> not_after_ms = read_not_after_ms(json, "request");
            const std::string id = json["attempt_id"].asString();
            auto [status, body] = attempts.cancel(id, wait_ms, not_after_ms);
            logger->info("[Server] Attempt {}: cancel requested (request {}): {} {}", id,
                         request_id, status, body.get("state", "").asString());
            if (status != 200)
                throw HttpError(status, std::move(body));
            return body;
        }, /* compact */ true, &traversal_io);
    };

    // The state of an attempt (running | stopping | finished, with the reason and when it
    // stopped), kept for the retention period after it finished; 404 after that
    server.resource["^/traverse/attempt/([^/]+)$"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                               shared_ptr<HttpServer::Request> request) {
        const std::string id = request->path_match[1].str();
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!valid_attempt_id(id)) {
                throw InvalidRequest("attempt_id: expected 1 to 128 characters from "
                                     "[A-Za-z0-9._:-]");
            }
            auto [status, body] = attempts.state(id);
            if (status != 200)
                throw HttpError(status, std::move(body));
            return body;
        }, /* compact */ true, &traversal_io);
    };

    // The content encodings of the traversal routes (compact JSON, Accept-Encoding honoured:
    // gzip preferred, deflate accepted)
    auto encodings_json = []() { return strings_json({ "gzip", "deflate" }); };

    // How a deadline reaches the walk (both capabilities routes): the time-sized chunks of
    // the annotation reads it may fall into, what stays uninterruptible, the longest single
    // piece seen, and how often an attempt's poll reads the clock; the rule is the SPEC's
    // (kDeadlineCheckRule)
    auto deadline_check_json = [&]() {
        Json::Value d;
        d["chunk_target_ms"] = static_cast<Json::UInt64>(config->traverse_chunk_target_ms);
        // No time bound on one piece is stated, and none will be: checkpoints bound the index
        // operations of a piece, not its wall time, which page faults and scheduling leave
        // open
        d["max_uninterruptible_ms"] = Json::Value();
        d["observed_max_uninterruptible_ms"]
            = static_cast<Json::UInt64>(attempts.observed_max_uninterruptible_ms());
        // one in poll_stride of a walk's polls (the one before a head) reads the clock for an
        // attempt's walk-until and bound, so a walk-until is seen up to poll_stride - 1 heads
        // late (the rule's number, stated as a number)
        d["poll_stride"] = static_cast<Json::UInt64>(attempts.settings().poll_stride);
        d["rule"] = kDeadlineCheckRule;
        return d;
    };

    // `in_ram` (both capabilities routes): the routes that accept it (/search's field: the
    // index loaded into RAM for the request), whether this server loads one for such a request
    // (a multi-graph server on mmap with --mem-cap-gb above 0; the others serve the index they
    // hold), the memory it lends to such loads, and that a request's budgets start after its
    // load (timing.load_ms states the load)
    auto in_ram_json = [&]() {
        Json::Value r;
        r["routes"] = strings_json({ "pattern", "resolve", "search", "traverse" });
        r["loads"] = !config->fnames.empty() && loaded_with_mmap && memory_all > 0;
        r["mem_cap_gb"] = number_json(config->memory_available);
        r["budgets_start"] = "after_load";
        return r;
    };

    // The keys both capabilities routes state alike: which walk a response gives (every
    // /traverse response names it too); the ledger-managed attempts (requests with
    // attempt_id): how they are named, cancelled, queried, kept and bounded
    // (traverse_attempts.hpp); the transport: the traversal routes write compact JSON and
    // honour Accept-Encoding (gzip preferred, deflate accepted), at this zlib level; how a
    // deadline reaches the walk; and /resolve's deadline (bounds.time_budget_ms; index-free, so
    // stated while the single index loads too). The deadline is not a feature_level bump, which
    // every /resolve and /traverse response states: a client gates on the block's presence
    auto put_contract = [&](Json::Value *c) {
        (*c)["algorithm_version"] = kTraverseAlgorithmVersion;
        (*c)["attempts"] = attempts.capabilities_json();
        (*c)["content_encodings"] = encodings_json();
        (*c)["compression_level"] = config->traverse_compression_level;
        (*c)["deadline_check"] = deadline_check_json();
        (*c)["in_ram"] = in_ram_json();
        (*c)["resolve"] = resolve_capabilities_json(resolve_time);
    };

    // What one index of this deployment supports, so a client can pick a strategy before
    // asking: GET /traverse/capabilities (per graph in multi-graph mode)
    auto probe_json = [&](const AnnotatedDBG &index, const IndexIdentity &identity) {
        graph::traversal::LabelOracle oracle(index);
        Json::Value caps = capabilities_to_json(oracle, config->index_release, &identity);
        caps["max_time_ms"] = config->traverse_max_time_ms;
        caps["max_seeds"] = static_cast<Json::UInt64>(config->traverse_max_seeds);
        caps["max_seed_bp"] = static_cast<Json::UInt64>(config->traverse_max_seed_bp);
        caps["max_seed_labels"] = static_cast<Json::UInt64>(config->traverse_max_seed_labels);
        caps["max_query_bp"] = static_cast<Json::UInt64>(config->resolve_max_query_bp);
        // the request budgets of DESIGN-traverse-graphlet.md §14 (bounds.max_memory_mb,
        // bounds.max_work_units) and W, the interval in charged work units at which the
        // walker reads the clock at the latest; stated here, not in the per-request
        // capabilities of every response. How far a work stop can exceed its budget is
        // stated with what bounds it, not as a fixed maximum: a fetch call's rows are
        // decoded and charged whole, and each stop states the most its seed charged
        // between two comparisons (no fixed kind of charge bounds them all)
        caps["budgets"] = strings_json({ "max_memory_mb", "max_work_units" });
        // the server's maxima of those budgets (feature level 4; 0: off): a larger
        // budget is lowered to it and an omitted one set to it, echoed in strategy.clamped
        caps["max_memory_mb"] = static_cast<Json::UInt64>(config->traverse_max_memory_mb);
        caps["max_work_units"] = static_cast<Json::UInt64>(config->traverse_max_work_units);
        caps["work_check_interval"]
            = static_cast<Json::UInt64>(graph::traversal::kWorkCheckInterval);
        // what bounds.max_work_units counts and how far a work stop can exceed it: the rule is
        // the SPEC's (work is deterministic logical work, not measured decode effort; the
        // physical decode counters are in each response's timing)
        caps["work_bound"] = kWorkBoundRule;
        caps["memory_bound"] = "soft";
        // the row-diff path cache of the reads (feature level 4)
        Json::Value decode_cache;
        decode_cache["path_cache_mb"] = static_cast<Json::UInt64>(config->traverse_path_cache_mb);
        decode_cache["rule"] = kDecodeCacheRule;
        caps["decode_cache"] = std::move(decode_cache);
        put_contract(&caps);
        // record coordinates (feature level 6): only here, the per-request capabilities
        // change only in their feature_level
        caps["coordinates"] = coordinates_capabilities_json(oracle);
        // the pattern search (SPEC-pattern-search.md §23): here because this is the document a
        // service's probe reads; the full block with `details`, the route of the full block
        caps["pattern"] = pattern_traverse_block(pattern_capabilities_json(
                &index, pattern_limits(*config), /* multi_graph */ false));
        return caps;
    };

    // What this deployment supports, so a client can pick a strategy before asking.
    server.resource["^/traverse/capabilities$"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                            shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError();
            const auto query = request->parse_query_string();
            if (config->fnames.empty()) {
                // a single index: the parameters that would select one are refused (they name
                // what this server does not have), any other is ignored, as it always was
                for (const auto &[key, value] : query) {
                    if (key == "graph" || key == "graph_path") {
                        throw InvalidRequest("Bad request: this server hosts a single graph; "
                                             "remove the 'graph' / 'graph_path' parameter");
                    }
                }
                return probe_json(*anno_graph.get(), single_identity);
            }
            // Multi-graph mode: the probe of one (graph, annotation) pair, selected by the
            // request fields' rules (?graph=<name>, and graph_path=<path> when the name spans
            // several graphs); it states which pair it describes
            const auto [name, graph_path] = probe_parameters(query, "/traverse/capabilities");
            const GraphPair &pair = select_traverse_pair(name, graph_path ? &*graph_path : nullptr,
                                                         indexes);
            const AnnotatedDBG &index = *graphs_cache.at(pair);
            Json::Value caps = probe_json(index, identity_of(index));
            caps["graph"] = name;
            caps["graph_path"] = pair.first;
            return caps;
        }, /* compact */ true, &traversal_io);
    };

    // The full pattern block (SPEC-pattern-search.md §23): answered while the single index
    // loads (available null, the graph fields null), so that a client learns the contract and
    // the caps before the index is ready. The parameters that would select a graph are refused
    // on a single-graph server, as on /traverse/capabilities; any other is ignored. A
    // multi-graph server answers the block of one pair, selected as /traverse/capabilities
    // selects it (?graph=<name>[&graph_path=<path>]), naming it (graph, graph_path)
    server.resource["^/pattern/capabilities$"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                           shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            const bool multi = !config->fnames.empty();
            const auto query = request->parse_query_string();
            if (multi) {
                const auto [name, graph_path] = probe_parameters(query, "/pattern/capabilities");
                const GraphPair &pair = select_traverse_pair(
                        name, graph_path ? &*graph_path : nullptr, indexes);
                Json::Value block = pattern_capabilities_json(
                        graphs_cache.at(pair).get(), pattern_limits(*config), false);
                block["graph"] = name;
                block["graph_path"] = pair.first;
                return block;
            }
            for (const auto &[key, value] : query) {
                if (key == "graph" || key == "graph_path") {
                    throw InvalidRequest("Bad request: this server hosts a single graph; "
                                         "remove the 'graph' / 'graph_path' parameter");
                }
            }
            const bool ready = anno_graph.wait_for(0s) == std::future_status::ready;
            return pattern_capabilities_json(ready ? anno_graph.get().get() : nullptr,
                                             pattern_limits(*config), false);
        }, /* compact */ true, &traversal_io);
    };

    // The server-wide capabilities (DESIGN-traverse-graphlet.md §17.1): the
    // routes and features this server offers, its mode and graphs, the attempts and how
    // deadlines are checked — answered while the single index loads (ready: false), so that a
    // service can learn the server's instance and contract before it routes anything to it
    server.resource["^/capabilities$"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                   shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            const bool multi = !config->fnames.empty();
            Json::Value c;
            put_contract(&c);
            c["feature_level"] = kTraverseFeatureLevel;
            Json::Value features(Json::arrayValue);
            Json::Value routes;
            routes["capabilities"] = "GET /capabilities";
            routes["search"] = "POST /search";
            features.append("search");
            if (!multi) {
                // a multi-graph server answers /align with 400
                routes["align"] = "POST /align";
                features.append("align");
            }
            routes["resolve"] = "POST /resolve";
            routes["traverse"] = "POST /traverse";
            routes["traverse_capabilities"] = multi
                ? "GET /traverse/capabilities?graph={name}[&graph_path={path}]"
                : "GET /traverse/capabilities";
            routes["cancel"] = "POST /traverse/cancel";
            routes["attempt"] = "GET /traverse/attempt/{attempt_id}";
            routes["column_labels"] = "GET /column_labels";
            routes["stats"] = "GET /stats";
            for (const char *f : { "resolve", "traverse", "attempts" }) {
                features.append(f);
            }
            // whether a graph can be searched is pattern.available (single-graph) or its
            // graph_summary entry's available (multi-graph); a client gates on both. The full
            // pattern block's own route beside it (SPEC-pattern-search.md §23), per pair on a
            // multi-graph server
            routes["pattern"] = "POST /pattern";
            routes["pattern_capabilities"] = multi
                ? "GET /pattern/capabilities?graph={name}[&graph_path={path}]"
                : "GET /pattern/capabilities";
            features.append("pattern");
            c["features"] = std::move(features);
            if (multi) {
                std::vector<std::string> names;
                for (const auto &[name, _] : indexes) {
                    names.push_back(name);
                }
                std::sort(names.begin(), names.end());
                Json::Value graphs(Json::arrayValue);
                for (const std::string &name : names) {
                    graphs.append(name);
                }
                c["graphs"] = std::move(graphs);
            } else {
                c["graphs"] = Json::Value();
            }
            c["mode"] = multi ? "multi" : "single";
            // a multi-graph server loads every index before it listens
            c["ready"] = multi || anno_graph.wait_for(0s) == std::future_status::ready;
            // the pattern search (DESIGN-pattern-search.md §7.3): its graph fields are null
            // while the single index loads. On a multi-graph server they are per pair
            // (graph_summary; GET /pattern/capabilities?graph=), and the block states the
            // contract and the caps with `available` true when /pattern serves one pair at least
            c["pattern"] = pattern_capabilities_json(
                    !multi && c["ready"].asBool() ? anno_graph.get().get() : nullptr,
                    pattern_limits(*config), false);
            if (multi) {
                // available when a pair is served; else the first pair's reason (each pair's is
                // in graph_summary)
                Json::Value reason;
                bool any = false;
                for (const Json::Value &pair : graph_summary["pairs"]) {
                    any |= pair["available"].asBool();
                    if (reason.isNull())
                        reason = pair["unavailable_reason"];
                }
                c["pattern"]["available"] = any;
                c["pattern"]["unavailable_reason"] = any ? Json::Value() : reason;
            }
            // per pair, so that a service learns a multi-graph server with one probe; null on a
            // single-graph server (pattern and GET /traverse/capabilities describe its graph)
            c["graph_summary"] = multi ? graph_summary : Json::Value();
            // the threshold of /search's and /pattern's selection rule (a request without
            // `graphs`), so that a client sees it before a longer list refuses its requests
            c["max_graphs_without_selection"]
                    = max_graphs_without_selection_json(multi, config->max_graphs_without_selection);
            // the largest request body the server reads (a longer one is dropped), so that a
            // client can check it before a request is cut off without a response
            c["max_request_body_mb"] = max_request_body_mb_json(config->max_request_body_mb);
            c["release"] = config->index_release;
            c["routes"] = std::move(routes);
            c["schema_version"] = 1;
            c["server_instance"] = attempts.server_instance();
            return c;
        }, /* compact */ true, &traversal_io);
    };

    server.resource["^/column_labels"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                    shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError(); // the index is not loaded yet, so we can't process the request

            Json::Value root(Json::arrayValue);
            if (!config->fnames.size()) {
                auto labels = anno_graph.get()->get_annotator().get_label_encoder().get_labels();
                for (const std::string &label : labels) {
                    root.append(label);
                }
            } else {
                for (const auto &[name, graphs] : indexes) {
                    for (const auto &[graph_fname, anno_fname] : graphs) {
                        const auto &labels = graphs_cache.at({ graph_fname, anno_fname })->get_annotator().get_label_encoder().get_labels();
                        for (const std::string &label : labels) {
                            root.append(label);
                        }
                    }
                }
            }
            return root;
        });
    };

    server.resource["^/stats"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                            shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError(); // the index is not loaded yet, so we can't process the request

            auto get_num_labels = [](const AnnotatedDBG &anno_dbg) {
                uint64_t num_labels = 0;
                if (const auto *coord_to_header = anno_dbg.get_coord_to_header()) {
                    for (uint64_t col = 0; col < coord_to_header->num_columns(); ++col) {
                        num_labels += coord_to_header->num_sequences(col);
                    }
                } else {
                    num_labels = anno_dbg.get_annotator().num_labels();
                }
                return num_labels;
            };

            Json::Value root;
            if (config->fnames.size()) {
                // for scenarios with multiple graphs
                const auto &graph = graphs_cache.begin()->second->get_graph();
                size_t k = graph.get_k();
                bool is_consistent_k = true;
                bool is_canonical = (graph.get_mode() == graph::DeBruijnGraph::CANONICAL);
                bool is_consistent_canonical = true;
                for (const auto &[graph_anno, anno_dbg] : graphs_cache) {
                    const auto &graph = anno_dbg->get_graph();
                    if (k != graph.get_k())
                        is_consistent_k = false;
                    if (is_canonical != (graph.get_mode() == graph::DeBruijnGraph::CANONICAL))
                        is_consistent_canonical = false;
                }
                if (is_consistent_k)
                    root["graph"]["k"] = static_cast<uint64_t>(k);
                if (is_consistent_canonical)
                    root["graph"]["is_canonical_mode"] = is_canonical;
                uint64_t num_labels = 0;
                for (const auto &[name, graphs] : indexes) {
                    for (const auto &[graph_fname, anno_fname] : graphs) {
                        num_labels += get_num_labels(*graphs_cache.at({ graph_fname, anno_fname }));
                    }
                }
                root["annotation"]["labels"] = num_labels;
            } else {
                root["graph"]["filename"] = std::filesystem::path(config->infbase).filename().string();
                root["graph"]["k"] = static_cast<uint64_t>(anno_graph.get()->get_graph().get_k());
                root["graph"]["nodes"] = anno_graph.get()->get_graph().num_nodes();
                root["graph"]["is_canonical_mode"] = (anno_graph.get()->get_graph().get_mode()
                                                        == graph::DeBruijnGraph::CANONICAL);
                const auto &annotation = anno_graph.get()->get_annotator();
                root["annotation"]["filename"] = std::filesystem::path(config->infbase_annotators.front()).filename().string();
                root["annotation"]["labels"] = get_num_labels(*anno_graph.get());
                root["annotation"]["objects"] = static_cast<uint64_t>(annotation.num_objects());
            }
            return root;
        });
    };

    server.default_resource["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                        shared_ptr<HttpServer::Request> request) {
        size_t request_id = num_requests++;
        logger->warn("[Server] Not found {} for {} request {} from {}",
                     request->path, request->method, request_id,
                     request->remote_endpoint().address().to_string());
        response->write(SimpleWeb::StatusCode::client_error_not_found,
                        "Could not find path " + request->path);
    };
    server.default_resource["POST"] = server.default_resource["GET"];

    server.on_error = [config](shared_ptr<HttpServer::Request> request,
                               const SimpleWeb::error_code &ec) {
        if (ec == SimpleWeb::make_error_code::make_error_code(SimpleWeb::errc::message_size)) {
            // a body over --max-request-body-mb: the library closes the connection without a
            // response (it builds a 413 it never sends), so this log line is the only trace
            const auto it = request->header.find("Content-Length");
            logger->warn("[Server] {} {} from {} dropped: its body ({}) exceeds "
                         "--max-request-body-mb {} MiB; connection closed without a response",
                         request->method, request->path,
                         request->remote_endpoint().address().to_string(),
                         it != request->header.end() ? it->second + " bytes" : "chunked",
                         config->max_request_body_mb);
            return;
        }
        // Handle errors here, ignoring a few trivial ones.
        if (ec.value() != asio::stream_errc::eof
                && ec.value() != asio::error::operation_aborted) {
            logger->warn("[Server] Got error {} {} {}",
                         ec.message(), ec.category().name(), ec.value());
        }
    };

    std::thread server_thread = start_server(server, *config, num_server_threads,
                                             [&]() { shutdown.started(&server); });
    server_thread.join();
    // The server stopped (only a shutdown stops it): the process exits here, without returning
    // through this function's destructors. They join what may still run: graph_loader an index
    // load in progress (single-index mode serves 503 while it loads, so a stop can come before
    // it ends — on a cold disk minutes later, the whole docker-stop timeout this shutdown is
    // there to avoid), graphs_pool its tasks. Nothing of the process's state
    // outlives it (an index is only read), so nothing is lost by not unwinding
    logger->info("[Server] Stopped");
    shutdown.finished();
    logger->flush();
    std::fflush(nullptr);
    std::_Exit(0);
}

} // namespace cli
} // namespace mtg
