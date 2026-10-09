#ifndef __LOAD_ANNOTATED_GRAPH_HPP__
#define __LOAD_ANNOTATED_GRAPH_HPP__

#include <cstdint>
#include <future>
#include <memory>
#include <string>

namespace mtg {

namespace graph {
class DeBruijnGraph;
class DBGSuccinct;
class AnnotatedDBG;
} // namespace graph

namespace cli {

class Config;

inline constexpr size_t kDefaultMaxChunksOpen = 2000;

std::unique_ptr<graph::AnnotatedDBG>
initialize_annotated_dbg(std::shared_ptr<graph::DeBruijnGraph> graph,
                         const Config &config,
                         size_t max_chunks_open = kDefaultMaxChunksOpen);

std::unique_ptr<graph::AnnotatedDBG> initialize_annotated_dbg(const Config &config);

/**
 * What the loading thread does for the pattern search to a graph it loads, before the graph is
 * shared with anyone (prepare_pattern_graph):
 *   check_mask       the sampled check of the graph's mask (check_mask_at_load)
 *   build_mask       --pattern-build-mask: the mask built in memory when there is none
 *   notices          the start-up log of a graph without its mask (or with one it could not
 *                    open), for its operator; a graph list logs one line for all its graphs
 *   sample_fraction  the dummy fraction of a graph without its mask sampled now
 *                    (sample_dummy_fraction_at_load); a graph prepared without it is sampled
 *                    at its first use (dummy_fraction)
 *   stdout_reserved  the CLI's stdout carries its answers: the progress is logged at trace
 *   progress         the check's and the sample's time logged at info level (else trace)
 */
struct PatternPreparation {
    bool check_mask = false;
    bool build_mask = false;
    bool notices = false;
    bool sample_fraction = false;
    bool stdout_reserved = false;
    bool progress = true;
};

// The preparation of a graph |config| loads: the CLI's and a single-graph server's graph
// (-i / -a) are served to the pattern search, with every step; any other none
PatternPreparation pattern_preparation(const Config &config);

// Does |prep| to |graph| (loaded from |path|), in this order: the mask check, the mask built
// (build_mask), the notices, the dummy fraction sampled
void prepare_pattern_graph(const std::shared_ptr<graph::DeBruijnGraph> &graph,
                           const std::string &path, const PatternPreparation &prep);

// Kick off graph loading on a worker thread, the graph prepared there for the pattern search
// (pattern_preparation(config), or |prep|) before the future is ready.
std::shared_future<std::shared_ptr<graph::DeBruijnGraph>>
async_load_critical_dbg(const Config &config);
std::shared_future<std::shared_ptr<graph::DeBruijnGraph>>
async_load_critical_dbg(const std::string &path, const PatternPreparation &prep);

// The graph of |config| (-i) loaded and prepared as |prep| says, its annotation (-a) loaded
// beside it in parallel (a server's per-request load of one pair, `in_ram`)
std::unique_ptr<graph::AnnotatedDBG>
initialize_annotated_dbg(const Config &config, const PatternPreparation &prep);

/**
 * --pattern-build-mask (DESIGN-pattern-search.md §4, mask: built_at_load): when |graph| is a
 * succinct graph loaded without its dummy-edge mask, build the mask in memory as
 * `transform --mask-dummy` and `build --mask-dummy` do (mask_dummy_edges), log the counts and
 * the time, and record |graph| as one whose mask was built at load. A graph that has its mask
 * (read from its .edgemask) keeps it, and one that is not succinct has none to build (the
 * pattern search refuses it for that). Exits the process when the mask cannot be built: the
 * operator asked for it (exact counts), and a server without it would answer upper bounds
 * throughout.
 * |stdout_reserved|: the caller's stdout carries its answers (`metagraph pattern`), where the
 * logger writes its info lines, so the progress is logged at trace level (stderr, with -v).
 */
void build_mask_at_load(const std::shared_ptr<graph::DeBruijnGraph> &graph,
                        bool stdout_reserved = false);

// Whether this process built |graph|'s mask at load (build_mask_at_load) rather than reading
// it from its file. |graph| is the served graph: a DBGSuccinct, or the CanonicalDBG that wraps
// a PRIMARY one.
bool mask_built_at_load(const graph::DeBruijnGraph &graph);

// The edges with W = $ the check of a mask looks up at load (sample_mask_sentinels)
inline constexpr uint64_t kMaskCheckSamples = 128;

// What sample_mask_sentinels looked up: the graph's edges with W = $ (plain or marked), how
// many of them it looked up (the main dummy edge 1 and the drawn ones: all of them when there
// are at most kMaskCheckSamples) and how many of those the mask marks valid
struct MaskSentinelSample {
    uint64_t sentinel_edges = 0;
    uint64_t looked_up = 0;
    uint64_t valid = 0;
};

/**
 * The sampled check of a dummy-edge mask: the main dummy edge 1 and |samples| edges with W = $
 * (plain or marked) drawn uniformly with replacement — by std::mt19937_64 seeded with the
 * graph's number of edges, as the dummy fraction is drawn, so the same edges in every process
 * of every platform — each looked up in the mask; every edge with W = $ when there are at most
 * |samples| of them (the full check then). A mask_dummy_kmers mask marks none of them valid. A
 * few dozen random reads of W and the mask: sub-second on a cold, busy host, where the full
 * count (DBGSuccinct::count_valid_sentinel_edges, a select on W per such edge) took 25 minutes
 * once. A mask that marks only a small fraction of those edges valid can pass it: the full
 * check is made where a mask is written (`transform --mask-dummy` refuses to write one that
 * fails it). Nothing is looked up without a mask (all zero).
 */
MaskSentinelSample sample_mask_sentinels(const graph::DBGSuccinct &graph,
                                         uint64_t samples = kMaskCheckSamples);

/**
 * The check of a mask the pattern search is served on (mask_invalid, SPEC-pattern-search.md
 * §6), at every load: sample_mask_sentinels, and when an edge it looked up is marked valid,
 * logs the remedy and records |graph| as one whose mask is invalid (mask_invalid_at_load): the
 * pattern search then refuses it (mask_invalid), since its counts would take those dummies for
 * k-mers and could be claimed exact while too large. Done once per load, in the loading thread,
 * for a graph the pattern search serves (prepare_pattern_graph); a mask built at load is
 * correct by construction and not checked. Returns the edges found marked valid (0: none, or no
 * mask, or not succinct). |stdout_reserved| as for build_mask_at_load; |progress|: whether the
 * check's time is logged at info level (else trace: a graph list's graphs, one line each)
 */
uint64_t check_mask_at_load(const std::shared_ptr<graph::DeBruijnGraph> &graph,
                            bool stdout_reserved = false, bool progress = true);

// Whether check_mask_at_load found |graph|'s mask invalid. |graph| is the served graph: a
// DBGSuccinct, or the CanonicalDBG that wraps a PRIMARY one.
bool mask_invalid_at_load(const graph::DeBruijnGraph &graph);

/**
 * Start loading the graph and the AnnotatedDBG in parallel. Returns futures
 * for both. The annotated DBG future resolves to nullptr when no annotation
 * is configured.
 */
std::pair<std::shared_future<std::shared_ptr<graph::DeBruijnGraph>>,
          std::future<std::unique_ptr<graph::AnnotatedDBG>>>
load_graph_with_async_annotation(const Config &config);

} // namespace cli
} // namespace mtg

#endif // __LOAD_ANNOTATED_GRAPH_HPP__
