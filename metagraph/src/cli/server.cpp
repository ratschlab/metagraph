#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <csignal>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <thread>

#include <unistd.h>

#include <json/json.h>
#include <tsl/hopscotch_map.h>
#include <server_http.hpp>
#include <csignal>

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
#include "server_cache.hpp"
#include "server_request.hpp"
#include "server_selection.hpp"
#include "server_envelope.hpp"
#include "server_fanout.hpp"
#include "server_profile.hpp"
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


/**
 * The pair a /search answer is computed on, as the answer states it
 * (SPEC-multigraph-envelope.md section 3): on a multi-graph server (|tagged|) its name and
 * paths, written on every row and on each query sequence's cut; on a single-graph server none,
 * the rows untagged and the cut's tags null. |index_fp|: the pair's identity (the server's on a
 * single-graph server), "" for none, stated null.
 */
struct SearchSource {
    bool tagged = false;
    std::string graph;
    std::string graph_path;
    std::string annotation_path;
    std::string index_fp;
};

Json::Value process_search_request(const SearchRequest &request,
                                   const graph::AnnotatedDBG &anno_graph,
                                   const SearchSource &source) {
    const auto &config = request.config;

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
    if (request.align) {
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
    tf.ofstream() << request.fasta;
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

    // Create full JSON object: per query sequence its rows, tagged with their pair on a
    // multi-graph server, and the cut its top N was taken with (one entry, for this pair)
    auto string_or_null = [](bool stated, const std::string &value) {
        return stated ? Json::Value(value) : Json::Value();
    };
    Json::Value search_response(Json::arrayValue);
    for (const auto &seq_result : search_results) {
        Json::Value json = seq_result.to_json(config.verbose_output, anno_graph.get_graph().get_k());
        if (source.tagged) {
            for (Json::Value &row : json["results"]) {
                row["graph"] = source.graph;
                row["graph_path"] = source.graph_path;
                row["annotation_path"] = source.annotation_path;
            }
        }
        Json::Value cut = seq_result.cut_to_json();
        cut["graph"] = string_or_null(source.tagged, source.graph);
        cut["graph_path"] = string_or_null(source.tagged, source.graph_path);
        cut["annotation_path"] = string_or_null(source.tagged, source.annotation_path);
        cut["index_fp"] = string_or_null(!source.index_fp.empty(), source.index_fp);
        json["cuts"] = Json::Value(Json::arrayValue);
        json["cuts"].append(std::move(cut));
        search_response.append(std::move(json));
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

    const auto &alignment_count = json.get(
        "max_alternative_alignments",
        (uint64_t)config.alignment_num_alternative_paths);
    if (!alignment_count.isInt() || alignment_count.asInt() <= 0) {
        throw std::invalid_argument("'max_alternative_alignments' must be a positive integer.");
    }
    config.alignment_num_alternative_paths = alignment_count.asInt();

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
 *  - once they returned, or after kDrainMs, the HTTP server's workers stop, and once they have
 *    returned its acceptor and every connection are closed (a response still being sent cut:
 *    incomplete, never a shorter complete one, its Content-Length unmet) and the process exits
 *    from run_server (_Exit: returning would join an index load still in progress in the
 *    loader's destructor);
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
        logger->info("[Server] Shutting down...");
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
        // Only the HTTP workers are stopped here (their io_service): closing the acceptor and
        // the connections from this thread would race with the accept and read handlers the
        // workers run, which crashed in asio (master's #681). run_server closes them once the
        // workers have returned (HttpServer::stop after the server thread's join)
        server->io_service->stop();
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
        // message_size, which answers 413 body_too_large (refuse_body_too_large; the library
        // prepares a 413 of its own there but does not send it: server_http.hpp, read() and
        // read_chunked_transfer_encoded()). A header block above it is closed without an answer
        server_startup.config.max_request_streambuf_size
                = static_cast<size_t>(config.max_request_body_mb) << 20;
    }

    logger->info("[Server] Will listen on {} port {}",
                 server_startup.config.address, server_startup.config.port);
    // the handler threads (this line said "Maximum connections" before H3), and the true cap
    // on the connections (--max-connections, CappedHttpServer)
    logger->info("[Server] HTTP handler threads: {}", num_threads);
    if (config.max_connections) {
        logger->info("[Server] Maximum connections: {} (--max-connections; a further connection "
                     "waits in the listen queue until one closes)", config.max_connections);
    } else {
        logger->info("[Server] Maximum connections: unlimited (--max-connections 0)");
    }
    if (config.max_request_body_mb)
        logger->info("[Server] Maximum request body: {} MiB", config.max_request_body_mb);
    return std::thread([&server_startup, on_accepting]() {
        server_startup.start([on_accepting](unsigned short port) {
            logger->info("[Server] Listening on port {}", port);
            if (on_accepting)
                on_accepting();
        });
    });
}

namespace {

// A request abandoned while it waited for the memory to load its index (`in_ram`): its client
// is gone, or the server stops; nothing was loaded, nothing is written
class LoadAbandoned : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

// A request whose completion deadline (finish_by_ms less the clock skew allowance) passed
// before, or while, it waited for the memory to load its index (SPEC-multigraph-envelope.md
// section 4: the wait is not entered then, and given up when it passes); nothing was loaded
class LoadExpired : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

// A pair of a multi-graph fan-out that was not started (PairControl refused its start: the
// request's work cutoff or its not_after_ms passed while it waited); nothing was loaded
class PairNotStarted : public std::runtime_error {
  public:
    PairNotStarted() : std::runtime_error("the pair was not started") {}
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
    IndexLease(std::unique_ptr<AnnotatedDBG> loaded, LoadReservation &&reservation, double load_ms)
          : reservation_(std::move(reservation)), loaded_(std::move(loaded)),
            index_(loaded_.get()), plan_(InRamPlan::LOAD), load_ms_(load_ms) {}
    IndexLease(const IndexLease&) = delete;
    IndexLease& operator=(const IndexLease&) = delete;

    const AnnotatedDBG& index() const { return *index_; }
    InRamPlan plan() const { return plan_; }
    double load_ms() const { return load_ms_; }

  private:
    // (in this order: the copy is freed before its bytes are given back)
    std::optional<LoadReservation> reservation_;
    std::unique_ptr<AnnotatedDBG> loaded_;
    const AnnotatedDBG *index_;
    InRamPlan plan_;
    double load_ms_ = 0;
};

// A body with more JSON values than --max-request-json-values, as every strict-JSON route
// refuses it (H1): 400 too_many_json_values, its members written in the order error, code, the
// cap (error_first_json_text). Nothing of the request ran
HttpError json_values_refusal(const TooManyJsonValues &e) {
    Json::Value body = too_many_json_values_json(e.max_values());
    std::string text = error_first_json_text(body);
    return HttpError(400, std::move(body), std::move(text));
}

// A /resolve listing more explicit seeds than --resolve-max-explicit-seeds (GPT review 7,
// finding 1): 400 too_many_explicit_seeds, its members written in the order error, code, the cap
// (error_first_json_text), as the intake refusals. Nothing of the request ran or was registered
HttpError explicit_seeds_refusal(const TooManyExplicitSeeds &e) {
    Json::Value body = too_many_explicit_seeds_json(e);
    std::string text = error_first_json_text(body);
    return HttpError(400, std::move(body), std::move(text));
}

// A /resolve answer priced past --resolve-max-discover-memory-mb (SPEC-labeled-traversal-core.md
// §4.6, GPT review 7, finding 2): 400 resolve_answer_too_large, its members written in the order
// error, code, the bound (error_first_json_text). The work ran; nothing partial is sent
HttpError answer_refusal(const graph::traversal::AnswerTooLarge &e) {
    Json::Value body = resolve_answer_too_large_json(e);
    std::string text = error_first_json_text(body);
    return HttpError(400, std::move(body), std::move(text));
}

// The body of /resolve, /traverse or /traverse/cancel, read as /pattern reads its body
// (parse_strict_json: one RFC 8259 JSON text, member names unique); a fault is 400
// {"error": "Bad json received: ..."} without a code, as /search's reader answers a body that
// is no JSON. Its values are counted against |max_values| (--max-request-json-values; 0:
// none) before it is parsed: one more is json_values_refusal
Json::Value parse_traversal_body(const std::string &content, uint64_t max_values) {
    try {
        return parse_strict_json<RefuseRequest>(content, "Bad json received: ", max_values);
    } catch (const TooManyJsonValues &e) {
        throw json_values_refusal(e);
    }
}

// The body of /pattern (parse_pattern_body: its faults are PatternRefusal 400
// invalid_request), its values counted as parse_traversal_body counts them
Json::Value parse_pattern_route_body(const std::string &content, uint64_t max_values) {
    try {
        return parse_pattern_body(content, max_values);
    } catch (const TooManyJsonValues &e) {
        throw json_values_refusal(e);
    }
}

// The caps of POST /traverse this server applies (Config::traverse_max_*, the H1 counts
// included), without the fields of one request
TraverseLimits traverse_limits(const Config &config) {
    TraverseLimits limits;
    limits.max_time_ms = config.traverse_max_time_ms;
    limits.max_seeds = config.traverse_max_seeds;
    limits.max_seed_bp = config.traverse_max_seed_bp;
    limits.max_seed_labels = config.traverse_max_seed_labels;
    limits.max_explicit_labels = config.traverse_max_explicit_labels;
    limits.max_extra_labels = config.traverse_max_extra_labels;
    limits.max_cost_entries = config.traverse_max_cost_entries;
    limits.max_memory_mb = config.traverse_max_memory_mb;
    limits.max_work_units = config.traverse_max_work_units;
    limits.max_request_memory_mb = config.traverse_max_request_memory_mb;
    limits.chunk_target_ms = static_cast<double>(config.traverse_chunk_target_ms);
    limits.path_cache_bytes = config.traverse_path_cache_mb << 20;
    return limits;
}

// The deadline block of a multi-graph server's /traverse envelope (SPEC-multigraph-envelope.md
// section 4): the request's not_after_ms and finish_by_ms, as its attempt holds them, and when
// the request gave either of them (null otherwise, as on /resolve and /pattern: RequestDeadlines
// ::json) the effective deadline and the work cutoff. /traverse's request budget is its
// attempt's bound, which finish_by_ms less the clock skew allowance and the transport's cap
// already bound, and its finalisation reserve the delivery reserve: the effective deadline is
// the receipt plus the bound, the cutoff the walk-until (Attempt::envelope_deadlines)
Json::Value traverse_envelope_deadline(const Attempt &attempt) {
    std::optional<std::chrono::system_clock::time_point> effective, cutoff;
    if (attempt.ids().not_after_ms || attempt.finish_by_ms()) {
        if (const auto deadlines = attempt.envelope_deadlines()) {
            effective = deadlines->effective;
            cutoff = deadlines->cutoff;
        }
    }
    return deadline_json(attempt.ids().not_after_ms, attempt.finish_by_ms(), effective, cutoff);
}

} // namespace

// The reverse index of the sequence headers (CoordToHeader::find_header), built while the
// index loads rather than by the first request naming a header: on refseq33m (33M headers)
// that request would spend seconds in it (9.3 s), inside its time budget and outside every
// read timer. Its sort runs on |num_threads| threads. An index whose columns or records exceed
// the header index's 32-bit fields is not served: the server exits naming the limit.
static void build_header_index(const graph::AnnotatedDBG &index, size_t num_threads) {
    const auto *coord_to_header = index.get_coord_to_header();
    if (!coord_to_header)
        return;
    Timer timer;
    try {
        const size_t distinct = coord_to_header->build_header_index(num_threads);
        // its entries (one per record) and the distinct headers among them: fewer when a
        // header names several records (headers_unique, E22)
        uint64_t entries = 0;
        for (uint64_t col = 0; col < coord_to_header->num_columns(); ++col) {
            entries += coord_to_header->num_sequences(col);
        }
        logger->info("[Server] Sequence header index built: {} entries, {} distinct headers, "
                     "in {:.1f} s ({} MiB)", entries, distinct, timer.elapsed(),
                     coord_to_header->header_index_bytes() >> 20);
    } catch (const std::length_error &e) {
        logger->error("[Server] {}", e.what());
        std::exit(1);
    }
}

// GET /capabilities' headers_unique (headers_unique_json) over |pairs|, logged with the time it
// took: a merge of the pairs' header indexes, which the load has built
static Json::Value headers_unique(const std::vector<HeaderPair> &pairs, size_t num_threads) {
    Timer timer;
    Json::Value certificate = headers_unique_json(pairs, num_threads);
    const uint64_t collisions = certificate["collisions"].asUInt64();
    if (certificate["value"].isNull()) {
        logger->info("[Server] headers_unique: null, a pair has no .seqs ({} repeated headers "
                     "among the pairs that have one), computed in {:.1f} s", collisions,
                     timer.elapsed());
    } else if (collisions) {
        logger->warn("[Server] headers_unique: false, {} headers occur more than once (one of "
                     "them: '{}'); a header label answers its first occurrence only. Computed "
                     "in {:.1f} s", collisions, certificate["examples"][0]["header"].asString(),
                     timer.elapsed());
    } else {
        logger->info("[Server] headers_unique: true, computed in {:.1f} s", timer.elapsed());
    }
    return certificate;
}

// GET /capabilities' labels_disjoint (labels_disjoint_json) over |pairs|, logged with the time it
// took (beside headers_unique's): the pairs' column names looked up in their header indexes
static Json::Value labels_disjoint(const std::vector<HeaderPair> &pairs, size_t num_threads) {
    Timer timer;
    Json::Value certificate = labels_disjoint_json(pairs, num_threads);
    const uint64_t collisions = certificate["collisions"].asUInt64();
    if (certificate["value"].isNull()) {
        logger->info("[Server] labels_disjoint: null, a pair has no .seqs ({} headers that are "
                     "also column names among the pairs that have one), computed in {:.1f} s",
                     collisions, timer.elapsed());
    } else if (collisions) {
        logger->warn("[Server] labels_disjoint: false, {} headers are also column names (one of "
                     "them: '{}'); a label with such a name selects the column where one has it. "
                     "Computed in {:.1f} s", collisions,
                     certificate["examples"][0]["header"].asString(), timer.elapsed());
    } else {
        logger->info("[Server] labels_disjoint: true, computed in {:.1f} s", timer.elapsed());
    }
    return certificate;
}

/**
 * GET /capabilities' graph_summary of a multi-graph server, computed once at start-up, so that
 * a service learns every pair with one probe: per (name, pair), in ordered_pairs' order, the
 * pair's profile (pair_summary_json: its tags and identity, k, alphabet, representation,
 * annotation form, decode_charged, what the pattern search makes of it — graph_mode, available
 * and unavailable_reason, mask and counting, as its GET /pattern/capabilities?graph= block
 * states them (the dummy fraction of a pair without its mask, sampled at load, is its block's
 * only) — what a traversal reads, and support_fp); and whether the pairs' columns are disjoint
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
    std::set<std::string> profiles;
    for (const auto &[name, pair] : ordered_pairs(indexes, sorted_names(indexes))) {
        Json::Value p = pair_summary_json(*graphs_cache.at(pair), identities.at(pair),
                                          name, pair.first, pair.second);
        masked += p["mask"].isString() && p["mask"].asString() != "absent";
        profiles.insert(p["support_fp"].asString());
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
                 "sampled at load. {} capability profile(s) (support_fp). Summary computed in "
                 "{:.3f} s", masked, summary["pairs"].size(), profiles.size(), timer.elapsed());
    return summary;
}

// The pairs of a multi-graph server's headers_unique and labels_disjoint: each once, tagged
// with the first name that lists it in ordered_pairs' order
static std::vector<HeaderPair> multi_graph_header_pairs(
        const GraphIndexes &indexes,
        const VectorMap<GraphPair, std::unique_ptr<AnnotatedDBG>> &graphs_cache) {
    std::vector<HeaderPair> pairs;
    std::set<GraphPair> seen;
    for (const auto &[name, pair] : ordered_pairs(indexes, sorted_names(indexes))) {
        if (!seen.insert(pair).second)
            continue;
        const AnnotatedDBG &index = *graphs_cache.at(pair);
        pairs.push_back({ name, pair.first, pair.second, index.get_coord_to_header(),
                          &index.get_annotator().get_label_encoder().get_labels() });
    }
    return pairs;
}

int run_server(Config *config) {
    assert(config);
    // SIGTERM and SIGINT shut the server down promptly (class Shutdown); allocated for the
    // process's lifetime, which its detached thread shares
    Shutdown &shutdown = *new Shutdown();
    shutdown.install();
    std::atomic<size_t> num_requests = 0;

    std::shared_future<std::unique_ptr<AnnotatedDBG>> anno_graph;
    // Destroy the loader before the future holding its result.
    ThreadPool graph_loader(1, 1);

    GraphIndexes indexes;

    ThreadPool graphs_pool(get_num_threads(), 1000 /* max_num_tasks */);
    size_t num_server_threads = std::max(1u, get_num_threads());
    // the HTTP server's handler threads, on which every route's handler runs (E19):
    // --http-threads, else -p's value, as before E19; -p keeps sizing the graphs' pool (above)
    // and the threads that load a graph list (num_server_threads)
    const size_t http_threads = config->http_threads ? config->http_threads : num_server_threads;
    set_num_threads(std::max(1u, config->parallel_each));
    config->parallel_each = 1;  // query one batch at a time
    logger->info("[Server] Threads per graph: {}", get_num_threads());

    VectorMap<std::pair<std::string, std::string>, std::unique_ptr<AnnotatedDBG>> graphs_cache;
    bool loaded_with_mmap = utils::with_mmap();

    // The index identity stated by /traverse and /resolve (DESIGN-traverse-graphlet.md
    // §3.1). Single-index mode: --index-name and --index-manifest, checked against the
    // loaded files before the index is served; the meta fingerprint once the index is
    // loaded. Multi-index mode: per (graph, annotation) pair, the name and the manifest of the
    // graph list's optional columns, checked before loading (a pair without them states null),
    // and the meta fingerprint computed once the pairs are loaded, before the server listens
    // (E17: GET /capabilities states it per pair).
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
    // fingerprint computed at start-up from the resident index (a copy of it loaded for a
    // request, `in_ram`, has the same files, the same identity; computed here only if missing)
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
    // and the identity of the indexes it serves together (GET /capabilities' composite_fp)
    Json::Value composite_fp;
    // The stamps of every file a copy of a pair loaded into RAM for one request reads
    // (pair_file_stamps: the loader inventory and the graph's derived data), taken at start-up
    // before any file of any pair was read (a multi-graph server)
    std::map<GraphPair, std::vector<PairFileStamp>> pair_stamps;
    // A request's load of |pair| into RAM is refused once a file of the pair changed on disk
    // since start-up (GPT review 7 finding 4): its answer would carry the start-up identity
    // (index_fp, index_meta_fp), the start-up certificates (headers_unique, labels_disjoint,
    // seqs_order, composite_fp) and the resident record mapping, and a replaced graph answered
    // silently wrong. Thrown as IndexChanged: the route's 400 with the refusal's text and the
    // code index_changed (the owner's decision of 2026-11-11), logged for the operator
    auto refuse_changed_pair = [&](const GraphPair &pair) {
        if (std::string why = pair_files_changed(pair_stamps.at(pair), pair.first, pair.second);
                !why.empty()) {
            logger->warn("[Server] in_ram refused: {}", why);
            throw IndexChanged(why);
        }
    };
    // the refusal of a copy whose annotation has other column labels than the resident's,
    // logged for the operator as refuse_changed_pair's (T4 review R3, N4)
    auto other_labels = [&](const GraphPair &pair) {
        const std::string why = fmt::format(
                "The annotation '{}' loaded into RAM for the request has other column labels "
                "than the one the server started with: restart onto the new files, or send the "
                "request without in_ram", pair.second);
        logger->warn("[Server] in_ram refused: {}", why);
        return IndexChanged(why);
    };
    // A copy of |pair| loaded into RAM for one request: with the resident pair's record mapping
    // (the .seqs the stamps vouch for, checked at start-up), so that no request reads a .seqs
    // or builds its header index, while every file of the pair is the one stamped at start-up,
    // before the load and after it (a file replaced during the load), and the copy's labels
    // the resident's (an annotation of another column count is refused by its load, without
    // reading its .seqs). Else refused (refuse_changed_pair, other_labels): no copy is served
    // from other files.
    // The status change time in the stamps makes a file swapped in and back within one load
    // differ too (it takes a rename or a link, which sets it); a directory or a symlink on the
    // path swapped and back does not. |load|(mapping) loads it
    auto load_copy = [&](const GraphPair &pair, const AnnotatedDBG &resident,
                         const std::function<std::unique_ptr<AnnotatedDBG>(
                                 std::shared_ptr<const annot::CoordToHeader>)> &load) {
        refuse_changed_pair(pair);
        std::unique_ptr<AnnotatedDBG> copy;
        try {
            copy = load(resident.share_coord_to_header());
        } catch (const RecordMappingMismatch &) {
            // the annotation loaded has another column count than the resident's mapping (it
            // was replaced during the load): refused, its .seqs never read (T4 review R3, S1)
            refuse_changed_pair(pair);
            throw other_labels(pair);
        }
        refuse_changed_pair(pair);
        if (copy->get_annotator().get_label_encoder().get_labels()
                != resident.get_annotator().get_label_encoder().get_labels())
            throw other_labels(pair);
        return copy;
    };
    // the bytes of |pair|'s graph and annotation, which an in_ram load reserves; a file that
    // cannot be stat'ed is refused as changed when it was there at start-up (it is gone)
    auto pair_bytes = [&](const GraphPair &pair) -> size_t {
        std::error_code graph_error;
        std::error_code annotation_error;
        const uintmax_t graph = std::filesystem::file_size(pair.first, graph_error);
        const uintmax_t annotation = std::filesystem::file_size(pair.second, annotation_error);
        if (graph_error || annotation_error) {
            refuse_changed_pair(pair);
            throw std::filesystem::filesystem_error(
                    "Cannot read the size of the index file",
                    graph_error ? pair.first : pair.second,
                    graph_error ? graph_error : annotation_error);
        }
        return graph + annotation;
    };
    // a single-graph server's headers_unique and labels_disjoint of GET /capabilities, set before
    // its index is ready
    Json::Value single_headers_unique;
    Json::Value single_labels_disjoint;
    // and its pair's profile (pair_summary of GET /capabilities, E17), set then too
    Json::Value single_pair_summary;
    // the list GET /column_labels answers on a multi-graph server (H3): each distinct pair
    // once, at its first place in the answers' order, its labels in column order; set once its
    // indexes are loaded (a single-graph server's list is its one annotation's)
    ColumnLabelList multi_column_labels;

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
            build_header_index(*anno_graph, get_num_threads() * num_server_threads);
            // the one pair, with null tags as a single-graph server's /search cuts
            const std::vector<HeaderPair> header_pair { { Json::Value(), Json::Value(),
                    Json::Value(), anno_graph->get_coord_to_header(),
                    &anno_graph->get_annotator().get_label_encoder().get_labels() } };
            single_headers_unique = headers_unique(header_pair,
                                                   get_num_threads() * num_server_threads);
            single_labels_disjoint = labels_disjoint(header_pair,
                                                     get_num_threads() * num_server_threads);
            // its profile, with null tags as its headers_unique
            single_pair_summary = pair_summary_json(*anno_graph, single_identity, Json::Value(),
                                                    Json::Value(), Json::Value());
            logger->info("[Server] Capability profile: {} {} annotation, support_fp {}, "
                         "index_meta_fp {}", single_pair_summary["graph_representation"].asString(),
                         single_pair_summary["annotation_form"].asString(),
                         single_pair_summary["support_fp"].asString(), single_identity.meta_fp);
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

        // every file a request's copy of a pair reads, stamped before any is read (load_copy)
        for (const auto &[pair, index] : graphs_cache) {
            pair_stamps[pair] = pair_file_stamps(pair.first, pair.second,
                                                 !config->no_coord_mapping);
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
            build_header_index(*index, get_num_threads() * num_server_threads);
        }
        // every response computed on a pair states its identity
        for (const auto &[pair, index] : graphs_cache) {
            const auto &[name, fp] = pair_identity.at(pair);
            IndexIdentity &id = pair_identities[pair];
            id.name = name;
            id.fp = fp;
            resident_pairs[index.get()] = pair;
        }
        // and its meta fingerprint, computed here for every pair (one pass over its column
        // names, the pairs in parallel), so that GET /capabilities states it per pair before
        // any request has read the pair (E17)
        {
            Timer timer;
            std::vector<std::string> meta_fps(graphs_cache.size());
            std::vector<uint64_t> num_columns(graphs_cache.size());
            #pragma omp parallel for num_threads(get_num_threads() * num_server_threads) schedule(dynamic)
            for (size_t i = 0; i < graphs_cache.size(); ++i) {
                const graph::traversal::LabelOracle oracle(*graphs_cache.nth(i).value());
                meta_fps[i] = index_meta_fingerprint(oracle);
                num_columns[i] = oracle.num_columns();
            }
            uint64_t columns = 0;
            for (size_t i = 0; i < graphs_cache.size(); ++i) {
                pair_identities[graphs_cache.nth(i).key()].meta_fp = meta_fps[i];
                columns += num_columns[i];
            }
            logger->info("[Server] index_meta_fp of the {} pairs ({} columns) computed in "
                         "{:.3f} s", graphs_cache.size(), columns, timer.elapsed());
        }
        graph_summary = multi_graph_summary(indexes, graphs_cache, pair_identities);
        // over each (graph, annotation) pair once, as loaded: a pair listed under several
        // names is one index
        std::vector<std::string> index_fps;
        for (const auto &[pair, index] : graphs_cache) {
            index_fps.push_back(pair_identities.at(pair).fp);
        }
        composite_fp = composite_fp_json(std::move(index_fps));
        logger->info("[Server] Composite fingerprint of the {} pairs: {}", graphs_cache.size(),
                     composite_fp.isNull() ? std::string("null (a pair has no manifest)")
                                           : composite_fp.asString());
        const std::vector<HeaderPair> header_pairs
                = multi_graph_header_pairs(indexes, graphs_cache);
        graph_summary["headers_unique"]
                = headers_unique(header_pairs, get_num_threads() * num_server_threads);
        graph_summary["labels_disjoint"]
                = labels_disjoint(header_pairs, get_num_threads() * num_server_threads);
        {
            // GET /column_labels' list (H3): the rule of headers_unique, each pair once
            std::vector<const std::vector<std::string> *> parts;
            std::set<GraphPair> seen;
            for (const auto &[name, pair] : ordered_pairs(indexes, sorted_names(indexes))) {
                if (seen.insert(pair).second) {
                    parts.push_back(&graphs_cache.at(pair)->get_annotator().get_label_encoder()
                                                          .get_labels());
                }
            }
            multi_column_labels = ColumnLabelList(std::move(parts));
        }
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

    size_t cache_size_bytes = static_cast<size_t>(config->server_cache_size * 1e9);
    ServerQueryCache search_cache(cache_size_bytes);
    if (cache_size_bytes) {
        logger->info("[Server] /search result cache enabled, max size: {} GB",
                     config->server_cache_size);
    } else {
        logger->info("[Server] /search result cache disabled (--cache-size 0)");
    }

    // The caps on the requests of /traverse, /resolve and /pattern whose handlers run at once
    // (E19, SPEC-labeled-traversal-core.md §10.3, "Handler capacity"; 0: none). A request of a
    // route at its cap is refused busy when its handler starts (refuse_busy), so that the other
    // routes keep a handler thread while the long calls run: the caps set sum to less than the
    // handler threads (Config refuses them otherwise)
    const RouteCaps route_caps { config->traverse_max_concurrent, config->resolve_max_concurrent,
                                 config->pattern_max_concurrent };
    RouteSlots traverse_slots("/traverse", "--traverse-max-concurrent", route_caps.traverse);
    RouteSlots resolve_slots("/resolve", "--resolve-max-concurrent", route_caps.resolve);
    RouteSlots pattern_slots("/pattern", "--pattern-max-concurrent", route_caps.pattern);
    if (config->http_threads || route_caps.traverse || route_caps.resolve || route_caps.pattern) {
        logger->info("[Server] Handler threads: {} ({}); requests run at once at most: /traverse "
                     "{}, /resolve {}, /pattern {} (0: no cap; above a cap a request is refused "
                     "503 busy, Retry-After {} s); the graphs' pool: {} thread(s) (-p)",
                     http_threads, config->http_threads ? "--http-threads" : "-p's value",
                     route_caps.traverse, route_caps.resolve, route_caps.pattern,
                     kBusyRetryAfterS, num_server_threads);
    }
    // the HTTP server, with the cap on the connections it holds open (H3: --max-connections,
    // 0 = none; a further connection waits in the listen queue until one closes)
    CappedHttpServer server;
    server.set_max_connections(config->max_connections);

    // Keep the cache handle alive until the asynchronous response write completes.
    server.resource["^/search"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                              shared_ptr<HttpServer::Request> request) {
        size_t request_id = num_requests++;
        auto handle_holder = cache_size_bytes
                ? std::make_shared<ServerQueryCache::Handle>() : nullptr;
        process_serialized_request(response, request, request_id,
            [&, handle_holder](const std::string &content) -> std::string {
                if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                    throw CurrentlyInitializingError();  // graph still loading

                Json::Value content_json = parse_json_string(content);
                if (logger->should_log(spdlog::level::trace))
                    logger->trace("[Server] Request {}: {}", request_id, content_json.toStyledString());

                const auto search_request = parse_search_request(content_json, *config);
                // Validate fields omitted from the key before returning a cached result.
                // the pairs searched, in the order the answer lists their rows and cuts
                // (ordered_pairs: the names in byte order, each name's pairs in the list's
                // order, a pair listed twice under a name once), each with its name
                std::vector<std::pair<std::string, GraphPair>> targets;
                bool load_in_ram = false;
                if (config->fnames.empty()) {
                    if (content_json.isMember("graphs"))
                        throw std::invalid_argument("Bad request: no support for filtering graphs on this server");
                    // nor for `pairs` (SPEC-multigraph-envelope.md section 6), as on a
                    // multi-graph server; before the cache, so that no cached answer is given
                    // for it (before release 10 the field was ignored here)
                    refuse_search_pairs(content_json);
                } else {
                    // the pairs of the selection (SPEC-multigraph-envelope.md section 6:
                    // `graphs` only, `pairs` refused with pairs_unsupported, an unknown name
                    // with graph_unknown) in the index's pair order (section 3: the names in byte
                    // order, each name's pairs in the list's order, a pair listed twice under a
                    // name once), refused before any pair is queried
                    for (SelectedPair &selected : select_pairs(
                                SelectionRoute::SEARCH, content_json, indexes,
                                config->max_graphs_without_selection)) {
                        targets.emplace_back(std::move(selected.name), std::move(selected.pair));
                    }
                    load_in_ram = content_json.get("in_ram", false).asBool() && loaded_with_mmap;
                }

                if (handle_holder) {
                    // the rows are tagged with the name a pair was selected under, so the key
                    // holds the names beside their pairs: a pair listed under two names is
                    // never answered with the other name's rows
                    const std::string graph_identity = config->fnames.empty()
                            ? config->infbase : search_pairs_identity(targets);
                    *handle_holder = search_cache.acquire(
                        make_search_cache_key(search_request, graph_identity));
                    if (!handle_holder->is_miss()) {
                        auto cached = handle_holder->get();
                        logger->info("[Server] Request {}: cache HIT  {:.1f} KB  "
                                     "(occupancy {:.1f}/{:.0f} MB, pinned {:.1f} MB, {} "
                                     "entries)",
                                     request_id, cached->size() / 1e3,
                                     search_cache.size_bytes() / 1e6,
                                     cache_size_bytes / 1e6,
                                     search_cache.pinned_bytes() / 1e6,
                                     search_cache.entry_count());
                        return *cached;
                    }
                }

                // Cache miss: run the search, publish the result. If the
                // compute throws, propagate the exception to any duplicate
                // waiters via set_exception() before re-throwing.
                Json::Value result;
                try {
                    if (!config->fnames.size()) {
                        logger->trace("Request {}: Started querying graph {}, in total graphs being queried at the moment: {}",
                                      request_id, config->infbase, graphs_being_queried.fetch_add(1) + 1);
                        try {
                            SearchSource source;
                            source.index_fp = single_identity.fp;
                            result = process_search_request(search_request, *anno_graph.get(),
                                                            source);
                        } catch (...) {
                            graphs_being_queried--;
                            throw;
                        }
                        graphs_being_queried--;
                    } else {
                        // each pair's answer in its place of |targets|, merged in that order
                        // once every pair has answered (not in the order the pairs finish)
                        std::vector<Json::Value> answers(targets.size());
                        std::vector<SearchSource> sources(targets.size());
                        for (size_t t = 0; t < targets.size(); ++t) {
                            const auto &[name, pair] = targets[t];
                            sources[t].tagged = true;
                            sources[t].graph = name;
                            sources[t].graph_path = pair.first;
                            sources[t].annotation_path = pair.second;
                            sources[t].index_fp = fp_of_pair(pair);
                        }
                        struct PendingQueries {
                            std::vector<std::shared_future<std::exception_ptr>> futures;
                            ~PendingQueries() {
                                // Tasks capture request-local state, including answers.
                                for (auto &future : futures)
                                    future.wait();
                            }
                        } pending;
                        auto &futures = pending.futures;
                        // Once a task is launched, recording its future must not allocate.
                        futures.reserve(targets.size());
                        for (size_t t = 0; t < targets.size(); ++t) {
                            {
                                const std::string &graph_fname = targets[t].second.first;
                                const std::string &anno_fname = targets[t].second.second;
                                futures.push_back(graphs_pool.enqueue([&,config,t,graph_fname=graph_fname,anno_fname=anno_fname]() {
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
                                        bool in_ram = load_in_ram;
                                        size_t index_size = in_ram ? pair_bytes({ graph_fname, anno_fname })
                                                                   : -1;
                                        if (in_ram && index_size > memory_all) {
                                            logger->warn("Request {}: Graph of size {} GB is too large to fit into "
                                                         "RAM (reserved memory: {} GB). It will be queried with mmap",
                                                         request_id, index_size / 1e9, memory_all / 1e9);
                                            in_ram = false;
                                        }
                                        Timer timer;
                                        if (in_ram) {
                                            // a pair whose files changed since start-up is
                                            // refused before its memory is waited for
                                            refuse_changed_pair({ graph_fname, anno_fname });
                                            reservations.reserve(index_size);
                                            index_size_reserved = index_size;
                                            logger->trace("Request {}: Loading graph {} of size {} GB to RAM...",
                                                          request_id, graph_fname, index_size / 1e9);
                                            timer.reset();
                                            Config config_copy = *config;
                                            config_copy.infbase = graph_fname;
                                            config_copy.infbase_annotators = { anno_fname };
                                            // as /traverse, /resolve and /pattern load a pair
                                            // into RAM: with the resident pair's mapping while
                                            // the pair's files are unchanged (load_copy), and a
                                            // load that fails or is refused throws, so that it
                                            // fails this request alone (the route's 400 with
                                            // its message), never the server (H4 review S2)
                                            const AnnotatedDBG &resident
                                                    = *graphs_cache.at({ graph_fname, anno_fname });
                                            config_copy.no_coord_mapping
                                                    = !resident.get_coord_to_header();
                                            PatternPreparation preparation;
                                            preparation.progress = false;
                                            index_loaded = load_copy(
                                                    { graph_fname, anno_fname }, resident,
                                                    [&](auto mapping) {
                                                return initialize_annotated_dbg(
                                                        config_copy, preparation, std::move(mapping));
                                            });
                                            index = index_loaded.get();
                                        } else {
                                            index = graphs_cache.at({ graph_fname, anno_fname }).get();
                                        }

                                        answers[t] = process_search_request(search_request, *index,
                                                                            sources[t]);
                                        logger->trace("Request {}: {} graph {} {} in {} sec",
                                                      request_id, in_ram ? "Loaded and searched" : "Searched",
                                                      graph_fname, anno_fname, timer.elapsed());

                                        index_loaded.reset();
                                        release_memory();
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
                        // the rows and the cuts of every query sequence in the pairs' order
                        for (Json::Value &json : answers) {
                            if (result.empty()) {
                                result = std::move(json);
                                continue;
                            }
                            if (json.size() != result.size())
                                throw std::logic_error("ERROR: Results for different sequences can't be merged");
                            for (Json::ArrayIndex i = 0; i < result.size(); ++i) {
                                if (result[i][SeqSearchResult::SEQ_DESCRIPTION_JSON_FIELD]
                                        != json[i][SeqSearchResult::SEQ_DESCRIPTION_JSON_FIELD]) {
                                    throw std::logic_error("ERROR: Results for different sequences can't be merged");
                                }
                                for (auto&& value : json[i]["results"]) {
                                    result[i]["results"].append(std::move(value));
                                }
                                for (auto&& value : json[i]["cuts"]) {
                                    result[i]["cuts"].append(std::move(value));
                                }
                            }
                        }
                    }
                } catch (...) {
                    if (handle_holder)
                        handle_holder->set_exception(std::current_exception());
                    throw;
                }

                // Serialize once: the cached form *is* the response
                // body. Cache hits then skip serialization entirely.
                std::string serialized = Json::writeString(Json::StreamWriterBuilder(), result);
                if (handle_holder) {
                    bool retained = handle_holder->set_result(serialized);
                    // BYPASS: larger than the budget, or than what the entries still being
                    // delivered (pinned) leave of it: delivered, not retained (H3)
                    logger->info("[Server] Request {}: cache {} {:.1f} KB  "
                                 "(occupancy {:.1f}/{:.0f} MB, pinned {:.1f} MB, {} entries)",
                                 request_id, retained ? "STORE" : "BYPASS", serialized.size() / 1e3,
                                 search_cache.size_bytes() / 1e6,
                                 cache_size_bytes / 1e6,
                                 search_cache.pinned_bytes() / 1e6,
                                 search_cache.entry_count());
                }
                return serialized;
            },
            [handle_holder, request_id](const SimpleWeb::error_code &ec) {
                if (!handle_holder || !*handle_holder)
                    return;  // never acquired (e.g. early throw before acquire)
                if (ec) {
                    logger->info("[Server] Request {}: delivery FAILED ({})", request_id, ec.message());
                    handle_holder->mark_protected();
                } else {
                    handle_holder->mark_delivered();
                }
            });
    };

    server.resource["^/align"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                             shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string &content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError(); // the index is not loaded yet, so we can't process the request

            if (config->fnames.size()) {
                throw std::invalid_argument("Bad request: alignment requests are not yet supported for "
                                            "servers with multiple graphs");
            }
            return process_align_request(content, anno_graph.get()->get_graph(), *config);
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
    //
    // A pair of a fan-out (|control|) states where it is: waiting for the memory (not started),
    // loading, running; it starts only when its control lets it (else PairNotStarted, nothing
    // held), and its wait polls the control's predicate (the request's stop, its not_after_ms,
    // its wait limit). The lease, and the reservation in it, belong to the caller: a fan-out's
    // worker holds them until it ends, whenever its request was answered. |expired| (/traverse
    // with finish_by_ms): the wait is not entered once it answers true and given up when it
    // does (LoadExpired)
    auto lease_index = [&](const GraphPair &pair, bool in_ram, size_t request_id,
                           const std::function<bool()> &gone, bool for_pattern,
                           PairControl *control, const std::function<bool()> &expired) {
        const AnnotatedDBG &resident = *graphs_cache.at(pair);
        const size_t bytes = in_ram && loaded_with_mmap ? pair_bytes(pair) : 0;
        const InRamPlan plan = in_ram_plan(in_ram, loaded_with_mmap, bytes, memory_all);
        if (plan == InRamPlan::RESIDENT_TOO_LARGE) {
            logger->warn("Request {}: Graph of size {} GB is too large to fit into RAM (reserved "
                         "memory: {} GB). It will be queried with mmap", request_id,
                         bytes / 1e9, memory_all / 1e9);
        }
        if (plan != InRamPlan::LOAD) {
            if (control && !control->start_running())
                throw PairNotStarted();
            return std::make_unique<IndexLease>(resident, plan);
        }
        // a pair whose files changed since start-up is refused before its memory is waited
        // for (load_copy refuses it again after the wait, before and after the load)
        refuse_changed_pair(pair);
        const char *expired_text = "the request's finish_by_ms (less the clock skew allowance) "
                                   "passed while it waited for the memory to load its index";
        if (expired && expired())
            throw LoadExpired(expired_text);
        // a pair of a fan-out does not enter the wait once its control's predicate answers
        // (the request stopped, its client gone, not_after_ms, its wait limit: finish_by_ms
        // inline, the cutoff pooled), as /traverse does not with LoadExpired: the memory may be
        // free, and a load then would run for an answer already due
        if (control && control->wait_predicate()())
            throw PairNotStarted();
        Timer timer;
        if (control)
            control->waiting_memory();
        const std::function<bool()> wait
                = control ? control->wait_predicate()
                : expired ? std::function<bool()>([&]() { return gone() || expired(); })
                : gone;
        if (!reservations.reserve(bytes, wait)) {
            if (control)
                throw PairNotStarted();
            if (expired && expired())
                throw LoadExpired(expired_text);
            throw LoadAbandoned("the client closed its connection (or the server stops) while "
                                "the request waited for the memory to load its index");
        }
        auto lease = load_under_reservation(LoadReservation(&reservations, bytes), [&]() {
            return !control || control->start_loading(bytes);
        }, [&]() {
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
            // annotation path the one loaded here), so its mapping is this annotation's
            // .seqs while the pair's files are the ones stamped at start-up (load_copy; else
            // the load is refused), and a resident without one (no .seqs beside the annotation,
            // or --no-coord-mapping) has none for its copy either: the .seqs is not looked for
            // again, so the operator's notice of a missing one is logged at start-up, not by
            // every request
            config_copy.no_coord_mapping = !resident.get_coord_to_header();
            return load_copy(pair, resident, [&](auto mapping) {
                return initialize_annotated_dbg(config_copy, preparation, std::move(mapping));
            });
        }, [&](std::unique_ptr<AnnotatedDBG> &&loaded, LoadReservation &&reservation) {
            const double load_ms = timer.elapsed() * 1000;
            logger->info("[Server] Request {}: ({}, {}) loaded into RAM for it in {:.0f} ms (the "
                         "wait for its {:.3f} GB included)", request_id, pair.first, pair.second,
                         load_ms, bytes / 1e9);
            return std::make_unique<IndexLease>(std::move(loaded), std::move(reservation),
                                                load_ms);
        });
        if (!lease)
            throw PairNotStarted();
        // the load ended: the pair's work begins
        if (control)
            control->start_running();
        return lease;
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
    attempt_settings.tombstone_hold_s = config->traverse_attempt_tombstone_hold_s;
    attempt_settings.clock_skew_ms = config->traverse_clock_skew_ms;
    attempt_settings.content_timeout_s = kContentTimeoutS;
    attempt_settings.observed_window_s = config->traverse_observed_window_s;
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

    // The requests of /resolve and /pattern by request_id (SPEC-multigraph-envelope.md section
    // 5): running, and ended ones for the retention period, answered by GET /requests/{id}
    RequestSettings request_settings;
    request_settings.retention_s = config->request_retention_s;
    request_settings.retention_count = config->request_retention;
    request_settings.server_instance = attempts.server_instance();
    // the clock of every admission's not_after_ms check (both routes' registries), which the
    // registry's 404 states as server_time_ms (section 5)
    request_settings.now_ms = [&attempts]() { return attempts.now_ms(); };
    RequestRegistry requests(request_settings);
    logger->info("[Server] Requests by request_id: {}", requests.retention_text());
    // a request_id and an attempt_id are one namespace (GET /requests/{id} reads both)
    IdAdmission id_admission(requests, attempts);

    // Whether a single index is still loading: every route that reads it answers 503 then
    auto index_loading = [&]() {
        return config->fnames.empty() && anno_graph.wait_for(0s) != std::future_status::ready;
    };
    // The busy refusal (E19) of a request whose route is at its cap, written by its handler when it
    // starts, before anything of the request is parsed or registered: 503 with Retry-After and
    // {code: "busy", error, server_instance} (busy_answer), one log line. While the server shuts
    // down, or to a client that is gone, nothing is written and the connection is closed (nothing
    // new is served, class Shutdown)
    auto refuse_busy = [&](const RouteSlots &slots, shared_ptr<HttpServer::Response> &response,
                           const shared_ptr<HttpServer::Request> &request, size_t request_id) {
        const bool silent = shutdown.stopping() || client_gone(*request);
        logger->info("[Server] {} request {} from {}: busy, {} of {} running ({}); {}",
                     slots.route(), request_id, request->remote_endpoint().address().to_string(),
                     slots.running(), slots.max(), slots.flag(),
                     silent ? "nothing written (the server stops or the client is gone)" : "503");
        if (silent) {
            response->close_connection_after_response = true;
            return;
        }
        write_answer(response, busy_answer(ServerBusy(slots.busy_message(),
                                                      attempts.server_instance())));
    };
    // The slot of a request whose handler refused nothing at its start because the single index
    // was loading, taken once the index is ready (inside process_request, after the loading 503):
    // at the route's cap the request is refused busy there (ServerBusy), or, while the server
    // stops, not answered (ClientGone)
    auto take_slot = [&](RouteSlots &slots) {
        RouteSlots::Slot slot = slots.try_acquire();
        if (!slot) {
            if (shutdown.stopping()) {
                throw ClientGone("the server stops, and the route runs its most requests at "
                                 "once: not run");
            }
            throw ServerBusy(slots.busy_message(), attempts.server_instance());
        }
        return slot;
    };

    // The deadlines of a /resolve or /pattern request whose handler starts now (section 4):
    // |budget_ms| the route's request budget, |reserve_ms| its finalisation reserve
    auto deadlines_of = [&](const shared_ptr<HttpServer::Request> &request, double budget_ms,
                            const RequestHandles &handles, double reserve_ms) {
        return request_deadlines(request->header_read_time, std::chrono::system_clock::now(),
                                 std::chrono::steady_clock::now(), budget_ms, handles,
                                 attempt_settings.clock_skew_ms, reserve_ms,
                                 static_cast<double>(kServerMaxDeadlineMs));
    };

    /**
     * The admission of a /resolve or /pattern request (sections 4 and 5), before anything of it
     * runs (IdAdmission::admit): one meant for another server_instance is refused (409
     * instance_mismatch, F3b), then a request_id that is running or retained here, on any
     * route (409, its entry and an error: no second worker), then one whose not_after_ms has
     * passed (409 expired, as /traverse refuses it); otherwise |entry| is registered under its
     * request_id (a request without one is not registered). A malformed field was refused
     * before (400), so that nothing a client must correct holds its id. On a multi-graph
     * server it is the request's: no pair starts before it
     */
    auto admit = [&](const std::shared_ptr<RequestEntry> &entry, const RequestHandles &handles,
                     size_t request_id) {
        std::optional<Json::Value> refused = id_admission.admit(entry, handles);
        if (refused) {
            logger->info("[Server] Request {}: not started, {}", request_id,
                         refused->get("error", "").asString());
            throw HttpError(409, std::move(*refused));
        }
    };

    // What a /resolve or /pattern request's response records in its registry entry (section 5,
    // none before it was admitted): a 400 found after the admission (only the work could find
    // it) frees the request_id as if it had been refused before it, any other status ends the
    // entry and retains it
    auto record_response = [&requests](const std::shared_ptr<RequestEntry> &entry, int status,
                                       std::optional<size_t> bytes) {
        if (!entry)
            return;
        if (status == 400)
            requests.forget(entry);
        entry->response(status, bytes);
    };

    // The pairs' refs of a request's registry entry: the selected pairs, or on a single-graph
    // server one pair with null tags
    auto pair_refs = [](const std::vector<SelectedPair> &pairs) {
        std::vector<Json::Value> refs;
        for (const SelectedPair &pair : pairs) {
            refs.push_back(envelope_pair_ref(pair));
        }
        if (refs.empty()) {
            Json::Value ref(Json::objectValue);
            for (const char *key : { "graph", "graph_path", "annotation_path" }) {
                ref[key] = Json::Value();
            }
            refs.push_back(std::move(ref));
        }
        return refs;
    };

    // The selected names of |pairs| (in byte order, each once): the envelope's `graphs`
    auto names_of = [](const std::vector<SelectedPair> &pairs) {
        std::vector<std::string> names;
        for (const SelectedPair &pair : pairs) {
            if (names.empty() || names.back() != pair.name)
                names.push_back(pair.name);
        }
        std::sort(names.begin(), names.end());
        names.erase(std::unique(names.begin(), names.end()), names.end());
        return names;
    };

    // two or more pairs run on the graphs' pool (one runs on its handler's thread); its queue is
    // bounded, and a handler that waits for room there must still see its cutoff and its client
    const Launch on_pool = [&graphs_pool](std::function<void()> task,
                                          SteadyClock::time_point until) {
        return graphs_pool.try_enqueue_until(until, std::move(task));
    };

    /**
     * POST /pattern on a multi-graph server: the pairs the request selects (select_pairs:
     * `graphs` or `pairs` with the coded refusals of SPEC-multigraph-envelope.md section 6,
     * every name when there are at most --max-graphs-without-selection and it names none),
     * each answered as one single-graph request would be (its own deadline, caps and
     * memory account, starting after its index is loaded for it with `in_ram`) by the fan-out
     * (run_fanout: one pair on this thread, two or more on the graphs' pool); the envelope
     * (post_envelope, SPEC-multigraph-envelope.md section 1) carries one entry per pair in the
     * index's pair order: the pair's answer, or the refusal the pair alone would have been
     * answered (SPEC §24.1), or its not_reached / unfinished outcome at the work cutoff
     * (section 4). What no pair decides (the body, the selection, `in_ram`, the request's
     * fields) is checked before any pair and refuses the whole request as before. Each pair
     * keeps its own deadline (from its start, capped by finish_by_ms and, pooled, by the
     * work cutoff, that cap stated only when it binds); two or more pairs
     * are written by the request's effective deadline, one pair by its own deadline, as a
     * single-graph request
     */
    auto multi_graph_pattern = [&](const std::string &content, size_t request_id,
                                   const shared_ptr<HttpServer::Request> &request,
                                   const std::function<bool()> &gone,
                                   PatternDelivery *delivery,
                                   std::shared_ptr<RequestEntry> *tracked) {
        Timer timer;
        // the request's CPU time and page faults: this thread's own part (the parse, the
        // checks, the envelope, an inline pair) and each pooled pair's work on its thread
        const common::ThreadUsageMeter own_usage;
        Json::Value json = parse_pattern_route_body(content, config->max_request_json_values);
        if (!json.isObject())
            throw PatternRefusal(400, "invalid_request", "request: expected an object");
        // the checks in the release-9 order (the selection, in_ram, then the fields), so that
        // a request it refused keeps its words; the request's handles (new in release 10, an
        // unknown field before) after its fields (validate_pattern_request reads them last)
        std::vector<SelectedPair> pairs;
        std::optional<bool> in_ram;
        try {
            pairs = select_pairs(SelectionRoute::PATTERN, json, indexes,
                                 config->max_graphs_without_selection);
            in_ram = in_ram_field(json);
        } catch (const std::invalid_argument &e) {
            throw PatternRefusal(400, "invalid_request", e.what());
        }
        const std::vector<std::string> names = names_of(pairs);
        // each pair's request: the body without the selection (`graphs` or `pairs`, SPEC-
        // multigraph-envelope.md section 6), its fields checked once here (the same for every
        // pair), so that their refusal is the request's and not a pair's
        json.removeMember("graphs");
        json.removeMember("pairs");
        const PatternLimits limits = pattern_limits(*config);
        const double budget_ms = validate_pattern_request(json, limits);
        RequestHandles handles;
        try {
            handles = read_request_handles(json);
        } catch (const std::invalid_argument &e) {
            throw PatternRefusal(400, "invalid_request", e.what());
        }
        const RequestDeadlines deadlines = deadlines_of(
                request, budget_ms, handles, static_cast<double>(config->pattern_finalize_ms));
        auto entry = std::make_shared<RequestEntry>(handles.request_id, "/pattern",
                                                    attempts.server_instance(),
                                                    deadlines.arrival_wall, pair_refs(pairs));
        admit(entry, handles, request_id);
        *tracked = entry;
        logger->info("[Server] Request {}: /pattern on {} pair(s) of {} graph name(s){}{}",
                     request_id, pairs.size(), names.size(),
                     in_ram.value_or(false) ? ", in_ram" : "",
                     handles.request_id ? ", request_id " + *handles.request_id : std::string());

        // A pair's work, on the thread that took it. It holds what it reads (the request's
        // body, in_ram, the request's number) and run_server's state, which lives as long as
        // the process; |delivery| and |gone| are this handler's and read by an inline pair only
        PairJob job = [&, json, in_ram, request_id, delivery](PairControl &control) {
            PairResult result;
            const SelectedPair &selected = control.pair();
            const GraphPair &pair = selected.pair;
            // the pair's entry when it is refused: what it alone would have been answered
            auto refused = [&](int http_status, const Json::Value &body) {
                result.entry = envelope_entry_refused(selected, fp_of_pair(pair), http_status,
                                                      body);
            };
            try {
                const auto lease = lease_index(pair, in_ram.value_or(false), request_id,
                                               control.stop_predicate(), /* for_pattern */ true,
                                               &control, nullptr);
                const IndexIdentity identity = identity_of_pair(pair, lease->index());
                // an inline pair is written under this handler's delivery (its own deadline,
                // as a single-graph request); a pooled one under its own, the envelope's
                // writing being the request's
                PatternDelivery own;
                PatternDelivery *pair_delivery = control.inline_run() ? delivery : &own;
                if (!control.inline_run())
                    own.set_abort(control.stop_predicate());
                if (const auto cap = control.deadline_cap())
                    pair_delivery->set_deadline_cap(*cap);
                // a pooled pair's budget is capped at the work cutoff, stated when it binds
                if (const auto cut = control.cutoff_cap())
                    pair_delivery->set_cutoff_cap(*cut);
                Json::Value answer = process_pattern_request(json, lease->index(), *config,
                                                             &identity, pair_delivery);
                if (in_ram)
                    answer["timing"]["load_ms"] = lease->load_ms();
                result.memory_bytes = answer["work"]["memory_bytes"].asUInt64();
                result.entry = envelope_entry_answered(selected, identity.fp, std::move(answer));
                result.answered = true;
            } catch (const PairNotStarted &) {
                result.not_started = true;
            } catch (const IndexChanged &e) {
                // a pair whose files changed since start-up (in_ram): its refusal block is the
                // refusal's body, code index_changed beside the text
                refused(e.status(), e.body());
            } catch (const PatternRefusal &e) {
                logger->warn("[Server] Request {}: pair ({}, {}) of '{}' refused ({} {}): {}",
                             request_id, pair.first, pair.second, selected.name, e.status(),
                             e.code(), e.what());
                refused(e.status(), e.body());
            } catch (const graph::pattern::Aborted &) {
                result.aborted = true;
            } catch (const LoadAbandoned &) {
                result.aborted = true;
            } catch (const std::exception &e) {
                // a failure of this pair, as the route answers one for a single graph
                // (answer_request): 400 with its text, no code
                logger->warn("[Server] Request {}: pair ({}, {}) of '{}' failed: {}",
                             request_id, pair.first, pair.second, selected.name, e.what());
                Json::Value body;
                body["error"] = e.what();
                refused(400, body);
            } catch (...) {
                logger->warn("[Server] Request {}: pair ({}, {}) of '{}' failed",
                             request_id, pair.first, pair.second, selected.name);
                Json::Value body;
                body["error"] = "Internal server error";
                refused(500, body);
            }
            return result;
        };
        FanoutSettings settings;
        settings.route = "/pattern";
        settings.pairs = pairs;
        for (const SelectedPair &pair : pairs) {
            settings.index_fps.push_back(fp_of_pair(pair.pair));
        }
        settings.deadlines = deadlines;
        settings.gone = gone;
        FanoutOutcome outcome = run_fanout(settings, entry, std::move(job), on_pool);
        if (outcome.client_gone)
            throw graph::pattern::Aborted();
        const bool answered_inline = pairs.size() == 1
                && outcome.entries[0]["outcome"].asString() == "answered";
        if (!answered_inline) {
            // the writing's deadline: the request's effective one (an inline pair that answered
            // is written by its own, as a single-graph request); a refused inline pair binds
            // nothing beyond the transport
            delivery->set_deadline(pairs.size() > 1
                ? graph::pattern::Deadline(deadlines.arrival,
                                           std::chrono::duration<double, std::milli>(
                                                   deadlines.effective - deadlines.arrival).count(),
                                           static_cast<double>(config->pattern_finalize_ms))
                : graph::pattern::Deadline(deadlines.arrival,
                                           static_cast<double>(kServerMaxDeadlineMs),
                                           static_cast<double>(config->pattern_finalize_ms)));
        }
        EnvelopeContext context;
        context.request_id = handles.request_id;
        context.server_instance = attempts.server_instance();
        context.deadline = deadlines.json();
        Json::Value out = post_envelope("/pattern", names, std::move(outcome.entries), context,
                                        timer.elapsed() * 1000);
        put_thread_usage(&out["timing"], pairs.size() > 1
                                             ? own_usage.elapsed().add(outcome.pairs_usage)
                                             : own_usage.elapsed());
        logger->info("[Server] Request {}: /pattern envelope: {} pair(s) answered, {} refused "
                     "({} not reached, {} unfinished)", request_id, out["answered"].asUInt(),
                     out["refused"].asUInt(), out["not_reached"].size(),
                     out["unfinished"].size());
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
        const size_t request_id = num_requests++;
        // the route's slot (E19), from the handler's start to its return; at the route's cap the
        // request is refused here, before anything of it is read (after the loading 503: a
        // request arriving while the single index loads takes its slot once it is ready)
        RouteSlots::Slot slot = pattern_slots.try_acquire();
        if (!slot && !index_loading()) {
            refuse_busy(pattern_slots, response, request, request_id);
            return;
        }
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
        // clock reading, the writing at its next check, and nothing is written (as on /resolve
        // and /traverse, which write their errors to it, /resolve until it is tracked)
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
        // the request's registry entry once it was admitted (every request on a multi-graph
        // server, one with request_id on a single-graph one): its response is recorded when
        // written, or that nothing was (SPEC-multigraph-envelope.md section 5)
        std::shared_ptr<RequestEntry> tracked;
        control.on_written = [&](int status, std::optional<size_t> bytes) {
            record_response(tracked, status, bytes);
        };
        process_request(response, request, request_id, [&](const std::string &content) {
            return translated([&]() {
                if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                    throw CurrentlyInitializingError();
                if (!slot)
                    slot = take_slot(pattern_slots);
                if (config->fnames.size()) {
                    return multi_graph_pattern(content, request_id, request, gone, &delivery,
                                               &tracked);
                }
                const AnnotatedDBG &index = *anno_graph.get();
                const IndexIdentity identity = identity_of(index);
                const Json::Value json = parse_pattern_route_body(
                        content, config->max_request_json_values);
                // the request's handles (sections 4 and 5): its not_after_ms admits it, its
                // finish_by_ms caps its deadline, its request_id registers it — once its
                // fields are valid, with or without a request_id: a request refused 400 holds
                // no id, and a malformed one is refused 400 also when its not_after_ms has
                // passed (the request's own 400s come before the expired 409)
                // (the fields first, the handles last: a request the release-9 server refused
                // keeps its words, as their fields were unknown to it)
                validate_pattern_request(json, pattern_limits(*config));
                RequestHandles handles;
                try {
                    handles = read_request_handles(json);
                } catch (const std::invalid_argument &e) {
                    throw PatternRefusal(400, "invalid_request", e.what());
                }
                std::shared_ptr<RequestEntry> entry;
                if (handles.request_id) {
                    entry = std::make_shared<RequestEntry>(
                            handles.request_id, "/pattern", attempts.server_instance(),
                            request->header_read_time, pair_refs({}));
                }
                admit(entry, handles, request_id);
                tracked = entry;
                if (handles.finish_by_ms) {
                    const RequestDeadlines d = deadlines_of(
                            request, 0, handles, static_cast<double>(config->pattern_finalize_ms));
                    delivery.set_deadline_cap(*d.finish);
                }
                if (entry)
                    entry->pair_running(0);
                Json::Value answer;
                try {
                    answer = process_pattern_request(json, index, *config, &identity, &delivery);
                } catch (...) {
                    if (entry) {
                        entry->pair_outcome(0, "refused");
                        entry->pair_ended(0);
                    }
                    throw;
                }
                if (entry) {
                    entry->pair_outcome(0, "answered");
                    entry->pair_ended(0, answer["work"]["memory_bytes"].asUInt64(),
                                      answer["timing"].get("cpu_ms", Json::Value()).isNumeric()
                                          ? std::optional<double>(
                                                  answer["timing"]["cpu_ms"].asDouble())
                                          : std::nullopt);
                }
                // in_ram: a single-graph server's index is the one it holds, nothing is loaded
                if (json.isObject() && json.isMember("in_ram"))
                    answer["timing"]["load_ms"] = 0.0;
                return answer;
            });
        }, /* compact */ true, &control);
    };

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
    // (and a discovery's bounds, SPEC-labeled-traversal-core.md §4.6)
    const ResolveTimeLimits resolve_time {
        config->traverse_max_time_ms > 0
                && config->traverse_max_time_ms < static_cast<double>(kServerMaxDeadlineMs)
            ? config->traverse_max_time_ms : static_cast<double>(kServerMaxDeadlineMs),
        kResolveFinalizeMs, config->resolve_max_discover_labels,
        config->resolve_max_discover_memory_mb, config->resolve_max_explicit_labels,
        config->resolve_max_discover_answer_labels, config->resolve_max_explicit_seeds };

    /**
     * POST /resolve on a multi-graph server (SPEC-multigraph-envelope.md section 1, E2): the
     * pairs the request selects (select_pairs), each answered as one single-graph /resolve
     * would be (its own caps and deadline, its identity, its in_ram load) by the fan-out, in
     * the envelope; a named graph answers one entry, on this thread. The request's own fields
     * (the body, the selection, in_ram, the query's length and budget) are checked before any
     * pair and refuse the whole request, as on a single-graph server
     */
    auto multi_graph_resolve = [&](Json::Value json, size_t request_id,
                                   const shared_ptr<HttpServer::Request> &request,
                                   ResolveWriting *writing,
                                   std::shared_ptr<RequestEntry> *tracked) {
        Timer timer;
        const common::ThreadUsageMeter own_usage;
        // the checks in the release-9 order, so that a request it refused keeps its words:
        // the selection's shape, in_ram, the pairs the names select, then the fields
        check_selection_shape(SelectionRoute::RESOLVE, json, indexes,
                              config->max_graphs_without_selection);
        const std::optional<bool> in_ram = in_ram_field(json);
        const std::vector<SelectedPair> pairs = select_pairs(
                SelectionRoute::RESOLVE, json, indexes, config->max_graphs_without_selection);
        // `pairs` is read: the pairs' requests are parsed without it (no field of the parser)
        json.removeMember("pairs");
        const std::vector<std::string> names = names_of(pairs);
        // what no pair decides refuses the whole request: its fields (the request's handles,
        // sections 4 and 5, read last, where they were unknown fields before release 10), the
        // query's length (the server's limit and bounds.max_query_bp), a budget not above the
        // reserve
        const ResolveRequest req = [&]() {
            try {
                return parse_resolve_request(json, resolve_time.max_explicit_labels,
                                             resolve_time.max_explicit_seeds);
            } catch (const TooManyExplicitSeeds &e) {
                throw explicit_seeds_refusal(e);
            }
        }();
        check_resolve_fields(req, config->resolve_max_query_bp, resolve_time);
        const RequestHandles handles = read_request_handles(json);
        // the request budget: bounds.time_budget_ms (lowered to the server's cap), else the
        // server's maximum
        const double budget_ms = req.time_budget_ms
                ? std::min(*req.time_budget_ms, resolve_time.max_time_ms)
                : resolve_time.max_time_ms;
        const RequestDeadlines deadlines = deadlines_of(request, budget_ms, handles,
                                                        resolve_time.finalize_ms);
        auto entry = std::make_shared<RequestEntry>(handles.request_id, "/resolve",
                                                    attempts.server_instance(),
                                                    deadlines.arrival_wall, pair_refs(pairs));
        admit(entry, handles, request_id);
        *tracked = entry;
        logger->info("[Server] Request {}: /resolve on {} pair(s) of {} graph name(s){}{}",
                     request_id, pairs.size(), names.size(),
                     in_ram.value_or(false) ? ", in_ram" : "",
                     handles.request_id ? ", request_id " + *handles.request_id : std::string());
        // A pair's own deadline is capped by finish_by_ms (a budget-less pair then gets the
        // server's maximum, process_resolve_request) and, pooled, by the work cutoff, stated
        // only when that cap stopped its work (the owner's decision of 2026-11-01); a pair
        // without a deadline (no budget, no finish_by_ms) has no budget to cap and is told to
        // stop at the cutoff (at the server's maximum)
        PairJob job = [&, json, in_ram, request_id, writing](PairControl &control) {
            PairResult result;
            const SelectedPair &selected = control.pair();
            const GraphPair &pair = selected.pair;
            auto refused = [&](int http_status, const Json::Value &body) {
                result.entry = envelope_entry_refused(selected, fp_of_pair(pair), http_status,
                                                      body);
            };
            try {
                const auto lease = lease_index(pair, in_ram.value_or(false), request_id,
                                               control.stop_predicate(), /* for_pattern */ false,
                                               &control, nullptr);
                const IndexIdentity identity = identity_of_pair(pair, lease->index());
                // an inline pair is written under this handler's delivery (its own deadline);
                // a pooled one's check is its own, the envelope's writing the request's
                ResolveDelivery own;
                ResolveDelivery *pair_delivery = control.inline_run() ? &writing->delivery : &own;
                if (const auto cap = control.deadline_cap())
                    pair_delivery->set_deadline_cap(*cap);
                // a pooled pair's budget is capped at the work cutoff, stated when it binds
                if (const auto cut = control.cutoff_cap())
                    pair_delivery->set_cutoff_cap(*cut);
                // the pair answers for the labels it holds and lists the others
                // (unknown_labels, SPEC-multigraph-envelope.md section 2)
                Json::Value out = process_resolve_request(
                        json, lease->index(), config->index_release, config->resolve_max_query_bp,
                        &identity, control.stop_predicate(), resolve_time, pair_delivery,
                        nullptr, /* labels_not_in_graph_dropped */ true);
                // in_ram: the time before the work began, waiting for and loading its index
                if (in_ram)
                    out["timing"]["load_ms"] = lease->load_ms();
                result.memory_bytes = out["work"]["memory_bytes"].asUInt64();
                result.entry = envelope_entry_answered(selected, identity.fp, std::move(out));
                result.answered = true;
            } catch (const PairNotStarted &) {
                result.not_started = true;
            } catch (const LoadAbandoned &) {
                result.aborted = true;
            } catch (const graph::traversal::AttemptAborted &) {
                result.aborted = true;
            } catch (const ResolveDeadline &e) {
                refused(503, resolve_deadline_body(e));
            } catch (const graph::traversal::AnswerTooLarge &e) {
                // the pair's answer priced past the bound: that pair refused, the others answered
                refused(400, resolve_answer_too_large_json(e));
            } catch (const HttpError &e) {
                refused(e.status(), e.body());
            } catch (const std::exception &e) {
                // what the pair alone would have been answered (answer_request): 400 with its
                // text, no code
                logger->warn("[Server] Request {}: pair ({}, {}) of '{}' failed: {}",
                             request_id, pair.first, pair.second, selected.name, e.what());
                Json::Value body;
                body["error"] = e.what();
                refused(400, body);
            } catch (...) {
                Json::Value body;
                body["error"] = "Internal server error";
                refused(500, body);
            }
            return result;
        };
        FanoutSettings settings;
        settings.route = "/resolve";
        settings.pairs = pairs;
        for (const SelectedPair &pair : pairs) {
            settings.index_fps.push_back(fp_of_pair(pair.pair));
        }
        settings.deadlines = deadlines;
        settings.gone = gone_of(request);
        FanoutOutcome outcome = run_fanout(settings, entry, std::move(job), on_pool);
        if (outcome.client_gone) {
            throw ClientGone("the client closed its connection (or the server stops) before "
                             "the answer was built");
        }
        const bool answered_inline = pairs.size() == 1
                && outcome.entries[0]["outcome"].asString() == "answered";
        if (pairs.size() > 1) {
            // the writing's deadline: the request's effective one
            const auto effective = deadlines.effective;
            const std::string late = "resolve: the envelope could not be written by the "
                    "request's effective deadline (" + iso_utc(deadlines.effective_wall)
                    + "): nothing partial is sent";
            writing->delivery.set_check([effective, late]() {
                if (std::chrono::steady_clock::now() > effective)
                    throw ResolveDeadline(late);
            });
        } else if (!answered_inline) {
            // a refused inline pair binds nothing beyond the transport
            writing->delivery.set_check(nullptr);
        }
        EnvelopeContext context;
        context.request_id = handles.request_id;
        context.server_instance = attempts.server_instance();
        context.deadline = deadlines.json();
        Json::Value out = post_envelope("/resolve", names, std::move(outcome.entries), context,
                                        timer.elapsed() * 1000);
        put_thread_usage(&out["timing"], pairs.size() > 1
                                             ? own_usage.elapsed().add(outcome.pairs_usage)
                                             : own_usage.elapsed());
        logger->info("[Server] Request {}: /resolve envelope: {} pair(s) answered, {} refused "
                     "({} not reached, {} unfinished)", request_id, out["answered"].asUInt(),
                     out["refused"].asUInt(), out["not_reached"].size(),
                     out["unfinished"].size());
        return out;
    };

    server.resource["^/resolve$"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                                shared_ptr<HttpServer::Request> request) {
        const size_t request_id = num_requests++;
        // the route's slot (E19), from the handler's start to its return; at the route's cap the
        // request is refused here, before anything of it is read (after the loading 503)
        RouteSlots::Slot slot = resolve_slots.try_acquire();
        if (!slot && !index_loading()) {
            refuse_busy(resolve_slots, response, request, request_id);
            return;
        }
        auto as_http = [](const ResolveDeadline &e) {
            return HttpError(503, resolve_deadline_body(e));
        };
        // a client that is gone, or a shutdown, is not answered: the work is abandoned at its
        // next phase (process_resolve_request), the writing and the compression of the answer
        // at the next check, the connection asked at most every client_check_ms, and nothing
        // is written (as on /pattern and /traverse); then the request's deadline, set once its
        // body is parsed (none without the field: the check then does nothing, and the text
        // and its compression are the same bytes). Its errors are written to it until the
        // request is tracked, its registry entry then recording client_gone (ResolveWriting)
        ResolveWriting writing(gone_of(request),
                               std::chrono::milliseconds(attempt_settings.client_check_ms),
                               traversal_io.compression_level);
        // the request's registry entry once it was admitted (every request on a multi-graph
        // server, one with request_id on a single-graph one): its response is recorded when
        // written, or that nothing was (SPEC-multigraph-envelope.md section 5)
        std::shared_ptr<RequestEntry> tracked;
        writing.tracked = [&tracked]() { return tracked != nullptr; };
        writing.control.on_written = [&](int status, std::optional<size_t> bytes) {
            record_response(tracked, status, bytes);
        };
        process_request(response, request, request_id, [&](const std::string &content) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError();
            if (!slot)
                slot = take_slot(resolve_slots);

            Json::Value json = parse_traversal_body(content, config->max_request_json_values);
            const bool multi = !config->fnames.empty();
            if (multi) {
                try {
                    return multi_graph_resolve(std::move(json), request_id, request, &writing,
                                               &tracked);
                } catch (const ResolveDeadline &e) {
                    throw as_http(e);
                }
            }
            traverse_graph_selection(json, multi);
            const std::optional<bool> in_ram = in_ram_field(json);
            // the request's handles (sections 4 and 5), checked before anything is registered:
            // by the parser of process_resolve_request, after every other field (they were
            // unknown fields before release 10, refused last), so that a request the
            // release-9 server refused keeps its words. Malformed, they are none here: the
            // parser refuses the request before it is admitted
            RequestHandles handles;
            try {
                handles = read_request_handles(json);
            } catch (const std::invalid_argument &) {}
            // its not_after_ms admits it, its finish_by_ms caps its deadline (a request without
            // bounds.time_budget_ms then gets the server's maximum), its request_id registers
            // it — once every refusal of its own is past (its fields, the query's length and
            // budget, its labels: process_resolve_request admits it then), so that a request
            // refused 400 holds no id and a malformed one is refused 400 also when expired
            std::shared_ptr<RequestEntry> entry;
            if (handles.request_id) {
                entry = std::make_shared<RequestEntry>(handles.request_id, "/resolve",
                                                       attempts.server_instance(),
                                                       request->header_read_time,
                                                       pair_refs({}));
            }
            const auto admitted = [&]() {
                admit(entry, handles, request_id);
                tracked = entry;
                if (entry)
                    entry->pair_running(0);
            };
            if (handles.finish_by_ms) {
                writing.delivery.set_deadline_cap(
                        *deadlines_of(request, 0, handles, resolve_time.finalize_ms).finish);
            }
            const IndexIdentity identity = single_identity;
            try {
                // a single-graph server keeps the 400 for a label it does not hold
                // (SPEC-multigraph-envelope.md section 2): multi is false here. A client that
                // is gone is not answered: abandoned between the request's phases
                Json::Value out = process_resolve_request(json, *anno_graph.get(),
                                                          config->index_release,
                                                          config->resolve_max_query_bp,
                                                          &identity, gone_of(request),
                                                          resolve_time, &writing.delivery,
                                                          nullptr, multi, admitted);
                // in_ram: a single-graph server's index is the one it holds, nothing is loaded
                if (in_ram)
                    out["timing"]["load_ms"] = 0.0;
                if (entry) {
                    entry->pair_outcome(0, "answered");
                    entry->pair_ended(0, out["work"]["memory_bytes"].asUInt64(),
                                      out["timing"]["cpu_ms"].isNumeric()
                                          ? std::optional<double>(out["timing"]["cpu_ms"].asDouble())
                                          : std::nullopt);
                }
                return out;
            } catch (...) {
                // (an entry refused before its admission was never registered)
                if (tracked) {
                    tracked->pair_outcome(0, "refused");
                    tracked->pair_ended(0);
                }
                try {
                    throw;
                } catch (const graph::traversal::AttemptAborted &e) {
                    throw ClientGone(e.what());
                } catch (const ResolveDeadline &e) {
                    throw as_http(e);
                } catch (const TooManyExplicitSeeds &e) {
                    throw explicit_seeds_refusal(e);
                } catch (const graph::traversal::AnswerTooLarge &e) {
                    throw answer_refusal(e);
                }
            }
        }, /* compact */ true, &writing.control);
    };

    // A multi-graph server's /traverse answer (SPEC-multigraph-envelope.md section 1.1): the
    // envelope with one entry, the |pair|'s single-graph answer |out| tagged with the pair and
    // |identity|'s index_fp (the identity its capabilities state), without its usage: the usage
    // is the attempt's, and the attempt is the request's, so the envelope states it at its top
    // level (for a ledger-managed request only, as a single-graph answer states it), the
    // attempt_id as its request_id (null without one) and the attempt registry's
    // server_instance. Only an answer of HTTP 200 is wrapped: every other status of the route
    // (409, 503 with usage at the bound, 400) stays HTTP-level and unwrapped, as the service's
    // attempt ledger reads it
    auto traverse_envelope = [&](const SelectedPair &pair, const IndexIdentity &identity,
                                 const Attempt &attempt, Json::Value out) {
        EnvelopeContext context;
        out.removeMember("usage", &context.usage);
        if (attempt.managed())
            context.request_id = attempt.ids().attempt_id;
        context.server_instance = attempts.server_instance();
        context.deadline = traverse_envelope_deadline(attempt);
        std::vector<Json::Value> entries;
        entries.push_back(envelope_entry_answered(pair, identity.fp, std::move(out)));
        return post_envelope("/traverse", { pair.name }, std::move(entries), context,
                             attempt.elapsed_ms());
    };

    // Extend frozen seeds along consistent annotation labels.
    server.resource["^/traverse$"]["POST"] = [&](shared_ptr<HttpServer::Response> response,
                                                 shared_ptr<HttpServer::Request> request) {
        const size_t request_id = num_requests++;
        // the route's slot (E19), from the handler's start to its return: released after
        // everything below, the attempt's response recorded. At the route's cap the request is
        // refused here, before anything of it is read and before its attempt exists (after the
        // loading 503)
        RouteSlots::Slot slot = traverse_slots.try_acquire();
        if (!slot && !index_loading()) {
            refuse_busy(traverse_slots, response, request, request_id);
            return;
        }
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
        // (a request without attempt_id states none: a 503 of its finish_by_ms is its message
        // and code alone)
        auto with_usage = [&](int status, const std::string &what, const std::string &reason) {
            Json::Value body;
            body["error"] = what;
            if (attempt->managed()) {
                body["usage"] = attempt->usage_json(reason);
            } else if (status == 503) {
                body["code"] = "deadline";
            }
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
        // until the handler returned, through the compression and the transport's copy. On a
        // multi-graph server they are the results of the envelope's one entry
        ResultTexts texts;
        const bool enveloped = !config->fnames.empty();
        control.write = [&texts, attempt, enveloped](const Json::Value &envelope,
                                                     const std::function<void()> &check) {
            // the longest stretch between two checks (deadline_check)
            double gap = 0;
            std::string text = !texts.active ? json_text(envelope, true, check, &gap)
                : enveloped
                    ? assemble_traverse_envelope(envelope, std::move(texts.texts), check, &gap)
                    : assemble_traverse_response(envelope, std::move(texts.texts), check, &gap);
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
            // (the delivery reserve; the latter without the record coordinates' share and the
            // names share, each with its exact text, so that an output with them measures what
            // the same output without them would: Attempt::note_delivered)
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
            if (!slot)
                slot = take_slot(traverse_slots);

            Json::Value json = parse_traversal_body(content, config->max_request_json_values);
            // a malformed id or not_after_ms is refused before anything is registered (400, no
            // usage)
            attempt->set_ids(attempt_ids(json));
            // a not_after_ms beyond the horizon (hard_cap_ms + clock skew ahead: H3, GPT review
            // 6 #6a) is refused here too, before anything is registered or held (400, no
            // usage): once run, its id would be held until it. A past one is never beyond it
            // (the expired 409 below stays)
            if (auto beyond = attempts.beyond_horizon_json(attempt->ids().not_after_ms,
                                                           "request")) {
                std::string text = error_first_json_text(*beyond);
                throw HttpError(400, std::move(*beyond), std::move(text));
            }
            // finish_by_ms, the completion deadline (SPEC-multigraph-envelope.md section 4):
            // the attempt's bound and walk-until become the minimum of their own and it less
            // the clock skew allowance. Malformed, it is refused before anything is registered
            // (400, no usage) — after the checks the release-9 server made before its parser
            // refused it as an unknown field (the selection, in_ram, the other fields), so that
            // a request it refused for another fault keeps its words
            try {
                attempt->set_finish_by(read_finish_by_ms(json, "request"));
            } catch (const std::invalid_argument &e) {
                if (!config->fnames.empty()) {
                    check_selection_shape(SelectionRoute::TRAVERSE, json, indexes,
                                          config->max_graphs_without_selection);
                    in_ram_field(json);
                    select_pairs(SelectionRoute::TRAVERSE, json, indexes,
                                 config->max_graphs_without_selection);
                    json.removeMember("pairs");
                } else {
                    traverse_graph_selection(json, false);
                    in_ram_field(json);
                }
                // its fields, finish_by_ms last (the counts first, before anything is copied)
                parse_traverse_request(json, traverse_limits(*config));
                throw;
            }
            if (attempt->managed()) {
                // an id does not run twice at once, nor again while it is retained or held, on
                // this server process: a second request with a running, retained, held or
                // tombstoned id is refused, without usage (it would be reconciled against the
                // other attempt); and one whose not_after_ms has passed is not started at all
                // (its ledger may already have released it), also without usage. Once its
                // retention and hold are over the id runs again
                if (auto refused = id_admission.start_attempt(attempt)) {
                    // the refusal's body is read with get(): operator[] of a non-const value
                    // adds a missing member (as null) to the body served
                    if (refused->instance_mismatch || refused->expired) {
                        logger->info("[Server] Attempt {} (request {}): not started, {}",
                                     attempt->ids().attempt_id, request_id,
                                     refused->body.get("error", "").asString());
                        throw HttpError(409, std::move(refused->body));
                    }
                    Json::Value body;
                    if (refused->body.get("tombstone", false).asBool()) {
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
                // a single-graph server refuses the selectors; a multi-graph one checks their
                // shape before in_ram and resolves the names after it, as the release-9 server
                // did (a request it refused keeps its words)
                if (!multi) {
                    traverse_graph_selection(json, false);
                } else {
                    check_selection_shape(SelectionRoute::TRAVERSE, json, indexes,
                                          config->max_graphs_without_selection);
                }
                const std::optional<bool> in_ram = in_ram_field(json);
                std::unique_ptr<IndexLease> lease;
                IndexIdentity identity = single_identity;
                // the pair a multi-graph server's envelope is answered for
                std::optional<SelectedPair> selected;
                if (multi) {
                    // exactly one pair (SPEC-multigraph-envelope.md section 6): `graph` or
                    // `graphs` resolving to one, or one object of `pairs`; several are refused
                    // with the candidates (traverse_single_pair)
                    selected = select_pairs(SelectionRoute::TRAVERSE, json, indexes,
                                            config->max_graphs_without_selection).front();
                    // `pairs` is read: the request is parsed without it (no field of the
                    // parser, which refuses it on a single-graph server as before release 10)
                    json.removeMember("pairs");
                    const GraphPair &pair = selected->pair;
                    // a malformed request is refused before its index is loaded for it
                    if (in_ram.value_or(false) && loaded_with_mmap)
                        parse_traverse_request(json, traverse_limits(*config));
                    // GET /requests/{attempt_id} names the pair the attempt ran on, for as
                    // long as the attempt is answered
                    attempt->set_pair_ref(pair_ref_json(*selected));
                    // the wait for the memory is abandoned when the client is gone or the server
                    // stops, and not entered or given up once finish_by_ms (less the clock skew
                    // allowance) has passed; a cancel takes effect when the work begins (its
                    // seeds not started)
                    std::function<bool()> expired;
                    if (const auto finish_by = attempt->finish_by_ms()) {
                        const uint64_t limit = *finish_by > attempt_settings.clock_skew_ms
                                ? *finish_by - attempt_settings.clock_skew_ms : 0;
                        expired = [&attempts, limit]() { return attempts.now_ms() > limit; };
                    }
                    lease = lease_index(pair, in_ram.value_or(false), request_id,
                                        gone_of(request), /* for_pattern */ false, nullptr,
                                        expired);
                    identity = identity_of_pair(pair, lease->index());
                }
                const AnnotatedDBG &index = lease ? lease->index() : *anno_graph.get();
                TraverseLimits limits = traverse_limits(*config);
                // a multi-graph server's graphs are the chunks of an index, which the fan-out
                // of /search sends each seed to whether or not they hold it: a seed a chunk does
                // not hold is its result, whichever spelling (`graph`, `graphs`) named the chunk
                limits.not_in_graph_per_seed = multi;
                // and its pairs need not hold every label a request names: the absent ones are
                // dropped and stated (SPEC-multigraph-envelope.md section 2)
                limits.labels_not_in_graph_dropped = multi;
                if (in_ram)
                    limits.load_ms = lease ? lease->load_ms() : 0.0;
                Json::Value out = process_traverse_request(json, index, config->index_release,
                                                           limits, &identity, attempt.get(),
                                                           &texts);
                if (!selected)
                    return out;
                return traverse_envelope(*selected, identity, *attempt, std::move(out));
            } catch (const LoadAbandoned &e) {
                throw ClientGone(e.what());
            } catch (const LoadExpired &e) {
                throw with_usage(503, e.what(), "deadline");
            } catch (const graph::traversal::AttemptAborted &e) {
                throw ClientGone(e.what());
            } catch (const AttemptAtBound &e) {
                throw with_usage(503, e.what(), "deadline");
            } catch (const SelectionRefusal &e) {
                // its code and candidates kept, with the usage of a managed request
                if (!attempt->managed())
                    throw;
                Json::Value body = e.body();
                body["usage"] = attempt->usage_json("error");
                throw HttpError(400, std::move(body));
            } catch (const IndexChanged &e) {
                // a pair whose files changed since start-up (in_ram): its code kept, written
                // error first, with the usage of a managed request
                if (!attempt->managed())
                    throw;
                Json::Value body = e.body();
                body["usage"] = attempt->usage_json("error");
                std::string text = error_first_json_text(body);
                throw HttpError(400, std::move(body), std::move(text));
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
            const CancelRequest cancel = parse_cancel_request(
                    parse_traversal_body(content, config->max_request_json_values));
            // the horizon of /traverse (H3, GPT review 6 #6a), before anything is looked up or
            // tombstoned: a cancel names the not_after_ms of a request /traverse accepts
            if (auto beyond = attempts.beyond_horizon_json(cancel.not_after_ms, "request")) {
                std::string text = error_first_json_text(*beyond);
                throw HttpError(400, std::move(*beyond), std::move(text));
            }
            auto [status, body] = attempts.cancel(cancel.attempt_id, cancel.wait_ms,
                                                  cancel.not_after_ms);
            logger->info("[Server] Attempt {}: cancel requested (request {}): {} {}",
                         cancel.attempt_id, request_id, status,
                         body.get("state", "").asString());
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

    // The state of a request by its handle (SPEC-multigraph-envelope.md section 5): a /resolve
    // or /pattern request's request_id (its pairs' states and usage, its response), kept for
    // the retention period after it ended, then 404; for /traverse the attempt_id is the handle
    // and the attempt registry answers (attempt: the object of GET /traverse/attempt/{id}), one
    // polling path for every route. The 404 states the server's clock read in the lookup
    // (server_time_ms; IdAdmission::state)
    server.resource["^/requests/([^/]+)$"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                       shared_ptr<HttpServer::Request> request) {
        const std::string id = request->path_match[1].str();
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!valid_attempt_id(id)) {
                throw InvalidRequest("request_id: expected 1 to 128 characters from "
                                     "[A-Za-z0-9._:-]");
            }
            auto [status, body] = id_admission.state(id);
            if (status == 200)
                return body;
            throw HttpError(404, std::move(body));
        }, /* compact */ true, &traversal_io);
    };

    // The content encodings of the traversal routes (compact JSON, Accept-Encoding honoured:
    // gzip preferred, deflate accepted)
    auto encodings_json = []() { return strings_json({ "gzip", "deflate" }); };

    // How a deadline reaches the walk (both capabilities routes): the time-sized chunks of
    // the annotation reads it may fall into, what stays uninterruptible, the longest single
    // piece seen within the last observed_window_s seconds, and how often an attempt's poll
    // reads the clock; the rule is the SPEC's (kDeadlineCheckRule)
    auto deadline_check_json = [&]() {
        Json::Value d;
        d["chunk_target_ms"] = static_cast<Json::UInt64>(config->traverse_chunk_target_ms);
        // No time bound on one piece is stated, and none will be: checkpoints bound the index
        // operations of a piece, not its wall time, which page faults and scheduling leave
        // open
        d["max_uninterruptible_ms"] = Json::Value();
        d["observed_max_uninterruptible_ms"]
            = static_cast<Json::UInt64>(attempts.observed_max_uninterruptible_ms());
        d["observed_window_s"] = static_cast<Json::UInt64>(attempts.settings().observed_window_s);
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
        // the counts refused before anything of a request is copied (H1; 0: none): one seed's
        // explicit labels, strategy.labels.extra, a change_cost table's entries
        caps["max_explicit_labels"]
            = static_cast<Json::UInt64>(config->traverse_max_explicit_labels);
        caps["max_extra_labels"] = static_cast<Json::UInt64>(config->traverse_max_extra_labels);
        caps["max_cost_entries"] = static_cast<Json::UInt64>(config->traverse_max_cost_entries);
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
        // the request-wide memory bound (E5b; 0: off): a seed after the first is not started
        // when the results before it and its budget would pass it (resource_stop scope
        // "request"); its presence says that this server can leave seeds unstarted for it
        caps["max_request_memory_mb"]
            = static_cast<Json::UInt64>(config->traverse_max_request_memory_mb);
        // the caps on the requests each long route runs at once (E19; 0: none), as GET
        // /capabilities' multigraph block states them: above one a request is refused 503 busy
        caps["max_concurrent"] = max_concurrent_json(route_caps);
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
        // service's probe reads; the gate block (the full block's gate fields and caps) with
        // `details`, the route of the full block
        caps["pattern"] = pattern_traverse_block(pattern_capabilities_json(
                &index, pattern_limits(*config), /* multi_graph */ false));
        return caps;
    };

    // a single-graph server refuses the parameters that would select a pair (they name what it
    // does not have); any other is ignored, as it always was (annotation_path too: a single-graph
    // server answers as before SPEC-multigraph-envelope.md)
    auto refuse_probe_selection = [](const SimpleWeb::CaseInsensitiveMultimap &query) {
        for (const auto &[key, value] : query) {
            if (key == "graph" || key == "graph_path") {
                throw InvalidRequest("Bad request: this server hosts a single graph; "
                                     "remove the 'graph' / 'graph_path' parameter");
            }
        }
    };

    // What this deployment supports, so a client can pick a strategy before asking.
    server.resource["^/traverse/capabilities$"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                            shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError();
            const auto query = request->parse_query_string();
            if (config->fnames.empty()) {
                refuse_probe_selection(query);
                return probe_json(*anno_graph.get(), single_identity);
            }
            // Multi-graph mode: the envelope (SPEC-multigraph-envelope.md section 1) of the
            // pairs the parameters select, each entry the probe of its (graph, annotation) pair
            // tagged with the pair and its index_fp: one pair, selected by the request fields'
            // rules (?graph=<name>, and graph_path=<path> when the name spans several graphs),
            // or without parameters every pair (up to the names cap)
            std::vector<Json::Value> entries;
            for (const SelectedPair &selected : select_pairs_from_query(
                     SelectionRoute::TRAVERSE, query, indexes,
                     config->max_graphs_without_selection)) {
                const AnnotatedDBG &index = *graphs_cache.at(selected.pair);
                const IndexIdentity identity = identity_of_pair(selected.pair, index);
                entries.push_back(envelope_entry_answered(selected, identity.fp,
                                                          probe_json(index, identity)));
            }
            return get_envelope("/traverse/capabilities", std::move(entries),
                                attempts.server_instance());
        }, /* compact */ true, &traversal_io);
    };

    // The full pattern block (SPEC-pattern-search.md §23): answered while the single index
    // loads (available null, the graph fields null), so that a client learns the contract and
    // the caps before the index is ready. The parameters that would select a graph are refused
    // on a single-graph server, as on /traverse/capabilities; any other is ignored. A
    // multi-graph server answers the envelope of the pairs selected as /traverse/capabilities
    // selects them (?graph=<name>[&graph_path=<path>], or every pair without parameters), each
    // entry the pair's block tagged with the pair and its index_fp
    server.resource["^/pattern/capabilities$"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                           shared_ptr<HttpServer::Request> request) {
        process_request(response, request, num_requests++, [&](const std::string&) {
            const bool multi = !config->fnames.empty();
            const auto query = request->parse_query_string();
            if (multi) {
                std::vector<Json::Value> entries;
                for (const SelectedPair &selected : select_pairs_from_query(
                         SelectionRoute::PATTERN, query, indexes,
                         config->max_graphs_without_selection)) {
                    entries.push_back(envelope_entry_answered(
                            selected, fp_of_pair(selected.pair),
                            pattern_capabilities_json(graphs_cache.at(selected.pair).get(),
                                                      pattern_limits(*config), false)));
                }
                return get_envelope("/pattern/capabilities", std::move(entries),
                                    attempts.server_instance());
            }
            refuse_probe_selection(query);
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
            // whether a header label names one record (E22): in graph_summary on a multi-graph
            // server, here on a single-graph one (null while its index loads)
            if (!multi)
                c["headers_unique"] = c["ready"].asBool() ? single_headers_unique : Json::Value();
            // and whether no header is also a column name (H4: a header label selects one
            // record only when both hold), placed as headers_unique
            if (!multi)
                c["labels_disjoint"] = c["ready"].asBool() ? single_labels_disjoint : Json::Value();
            // the capability profile of a single-graph server's pair (E17): a graph_summary
            // entry with null tags, null while its index loads (each pair's is in graph_summary
            // on a multi-graph server)
            if (!multi)
                c["pair_summary"] = c["ready"].asBool() ? single_pair_summary : Json::Value();
            // the threshold of /search's and /pattern's selection rule (a request without
            // `graphs`), so that a client sees it before a longer list refuses its requests
            c["max_graphs_without_selection"]
                    = max_graphs_without_selection_json(multi, config->max_graphs_without_selection);
            // the largest request body the server reads (a longer one is refused 413
            // body_too_large before it is read), so that a client can check it before it sends
            c["max_request_body_mb"] = max_request_body_mb_json(config->max_request_body_mb);
            // what the server reads of a request (H1): that cap and its 413, the value count of
            // the strict-JSON routes and its 400; in the multigraph block on a multi-graph
            // server, here on a single-graph one (one reader for both, as
            // expect_server_instance)
            const Json::Value intake = request_intake_json(config->max_request_body_mb,
                                                           config->max_request_json_values,
                                                           config->max_connections);
            if (!multi)
                c["request_intake"] = intake;
            // GET /column_labels (H3): its byte cap, its paging, its 400 and the list's length
            // (null while a single index loads); placed as request_intake
            std::optional<uint64_t> labels_total;
            if (multi) {
                labels_total = multi_column_labels.total();
            } else if (c["ready"].asBool()) {
                labels_total = anno_graph.get()->get_annotator().get_label_encoder().size();
            }
            const Json::Value column_labels = column_labels_capabilities_json(
                    config->max_column_labels_bytes, labels_total);
            if (!multi)
                c["column_labels"] = column_labels;
            c["release"] = config->index_release;
            if (multi) {
                // the identity of the indexes served together (SPEC-multigraph-envelope.md
                // section 8), the rule it is computed by, and what the multi-graph routes do
                // in this release: the routes answering in the envelope, the selectors, the
                // deadlines and the request registry
                c["composite_fp"] = composite_fp;
                c["composite_rule"] = kCompositeRule;
                c["multigraph"] = multigraph_capabilities_json(
                        config->traverse_max_request_memory_mb, route_caps, intake,
                        column_labels,
                        { resolve_time.max_discover_labels,
                          resolve_time.max_discover_memory_mb,
                          resolve_time.max_explicit_labels,
                          resolve_time.max_discover_answer_labels,
                          resolve_time.max_explicit_seeds });
            }
            // the request registry (section 5) on every server: a single-graph /resolve or
            // /pattern with request_id registers too, and GET /requests/{attempt_id} reads a
            // /traverse attempt
            routes["requests"] = kRequestsRoute;
            // the routes that take expect_server_instance (section 5, F3b): in the multigraph
            // block on a multi-graph server, here on a single-graph one (one reader for both)
            if (!multi)
                c["expect_server_instance"] = expect_server_instance_routes_json();
            c["routes"] = std::move(routes);
            c["schema_version"] = 1;
            c["server_instance"] = attempts.server_instance();
            return c;
        }, /* compact */ true, &traversal_io);
    };

    // The column labels (H3, review 5 #4; SPEC-labeled-traversal-core.md §10.3): a single-graph
    // server's annotation in column order; a multi-graph server's pairs each once, at their
    // first place in the answers' order (multi_column_labels). Without parameters the array as
    // before, refused 400 column_labels_too_large when its labels' JSON strings pass
    // --max-column-labels-bytes; with ?offset=&limit= a page. Written as text directly: no copy
    // of a label vector, no JSON tree of the whole list
    server.resource["^/column_labels"]["GET"] = [&](shared_ptr<HttpServer::Response> response,
                                                    shared_ptr<HttpServer::Request> request) {
        std::string text;
        ResponseControl control;
        control.write = [&](const Json::Value &, const std::function<void()> &) {
            return std::move(text);
        };
        process_request(response, request, num_requests++, [&](const std::string&) {
            if (!config->fnames.size() && anno_graph.wait_for(0s) != std::future_status::ready)
                throw CurrentlyInitializingError(); // the index is not loaded yet, so we can't process the request

            const auto query = request->parse_query_string();
            const ColumnLabelsPage page = column_labels_parameters(
                    std::vector<std::pair<std::string, std::string>>(query.begin(), query.end()));
            const ColumnLabelList single_list = config->fnames.empty()
                ? ColumnLabelList({ &anno_graph.get()->get_annotator().get_label_encoder()
                                                       .get_labels() })
                : ColumnLabelList();
            const ColumnLabelList &list = config->fnames.empty() ? single_list
                                                                 : multi_column_labels;
            const uint64_t max_bytes = config->max_column_labels_bytes;
            if (page.paged) {
                text = column_labels_page_text(list, page, max_bytes);
                return Json::Value();
            }
            std::optional<std::string> all = column_labels_text(list, max_bytes);
            if (!all) {
                Json::Value body = column_labels_too_large_json(list.total(), max_bytes);
                std::string refusal = error_first_json_text(body);
                throw HttpError(400, std::move(body), std::move(refusal));
            }
            text = std::move(*all);
            return Json::Value();
        }, /* compact */ false, &control);
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
            // a body over --max-request-body-mb, before it is read (a Content-Length above it)
            // or at the chunk that passes it: the library prepares a 413 but does not send it
            // and drops the request, so the 413 body_too_large is written here, on its
            // connection (refuse_body_too_large: the connection closed after it, what the
            // client still sends discarded meanwhile), or the connection is closed when the
            // Content-Length is no byte count. Nothing of the request ran
            const auto it = request->header.find("Content-Length");
            const std::string address = request->remote_endpoint().address().to_string();
            const bool answered = refuse_body_too_large(request, config->max_request_body_mb);
            const std::string length = it != request->header.end()
                    ? ::utils::clip_for_message(it->second, 64) : std::string();
            if (answered) {
                // the count as the 413 states it (content_length_value: " 3145728 " is 3145728)
                const auto count = it != request->header.end()
                        ? content_length_value(it->second) : std::nullopt;
                logger->warn("[Server] {} {} from {} refused 413: its body ({}) exceeds "
                             "--max-request-body-mb {} MiB",
                             request->method, ::utils::clip_for_message(request->path, 256),
                             address, count ? std::to_string(*count) + " bytes" : "chunked",
                             config->max_request_body_mb);
            } else {
                logger->warn("[Server] {} {} from {} dropped without an answer: its "
                             "Content-Length '{}' is no byte count",
                             request->method, ::utils::clip_for_message(request->path, 256),
                             address, length);
            }
            return;
        }
        // Handle errors here, ignoring a few trivial ones.
        if (ec.value() != asio::stream_errc::eof
                && ec.value() != asio::error::operation_aborted) {
            logger->warn("[Server] Got error {} {} {}",
                         ec.message(), ec.category().name(), ec.value());
        }
    };

    std::thread server_thread = start_server(server, *config, http_threads,
                                             [&]() { shutdown.started(&server); });
    server_thread.join();
    // No HTTP worker is left: close the acceptor and the open connections (a response still
    // being sent is cut)
    server.stop();
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
