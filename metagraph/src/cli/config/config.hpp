#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <filesystem>
#include <string>
#include <vector>

#include "kmer/kmer_collector_config.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/representation/base/sequence_graph.hpp"
#include "cli/query.hpp"


namespace mtg {
namespace cli {

// The HTTP server's content timeout (s): the body, the handler and the response of every
// request must fit in it (Simple-Web-Server's timeout_content, server.cpp), and the second a
// route's own deadline leaves under it for the transport (/traverse's attempts' hard cap; the
// cap of /pattern's time_budget_ms and of /resolve's opt-in deadline, review of 2026-10-07,
// R1-05)
constexpr uint64_t kServerContentTimeoutS = 900;
constexpr uint64_t kServerTransportMarginMs = 1000;
// the longest route deadline the transport can honour (ms)
constexpr uint64_t kServerMaxDeadlineMs = kServerContentTimeoutS * 1000 - kServerTransportMarginMs;

class Config {
  public:
    Config(int argc, char *argv[]);

    static constexpr auto UNINITIALIZED_STR = "\0";

    bool print_graph = false;
    bool print_graph_internal_repr = false;
    bool print_column_names = false;
    bool print_counts_hist = false;
    bool forward_and_reverse = false;
    bool complete = false;
    bool dynamic = false;
    bool mark_dummy_kmers = false;
    bool filename_anno = false;
    bool annotate_sequence_headers = false;
    bool to_adj_list = false;
    bool to_fasta = false;
    bool enumerate_out_sequences = false;
    bool to_gfa = false;
    bool output_compacted = false;
    bool unitigs = false;
    bool kmers_in_single_form = false;
    bool initialize_bloom = false;
    bool count_kmers = false;
    bool query_presence = false;
    bool verbose_output = false;
    bool filter_present = false;
    bool dump_text_anno = false;
    bool sparse = false;
    bool subsample_rows = false;
    bool batch_align = false;
    bool suppress_unlabeled = false;
    bool in_ram = false;
    bool clear_dummy = false;
    bool count_dummy = false;
    bool greedy_brwt = false;
    bool cluster_linkage = false;
    bool separately = false;
    bool map_sequences = false;
    bool align_sequences = false;
    bool align_only_forwards = false;
    bool filter_by_kmer = false;
    bool output_json = false;
    bool aggregate_columns = false;
    bool coordinates = false;
    bool index_header_coords = false;
    bool no_coord_mapping = false;
    // Opt-in: include the bulky VG-style `path.mapping[]` object in
    // `metagraph align --json` output. Off by default since nothing in-tree
    // consumes it and it dominates the JSON size for long alignments.
    bool align_output_path = false;
    bool advanced = false;

    unsigned int k = 3;

    // Cache ranges of nodes in succinct graphs to search faster.
    // For DNA4, index nodes for all possible suffixes of length 12.
    // In general, the default value is: log_{|Sigma|}(2^24)
    static const size_t kDefaultIndexSuffixLen;
    unsigned int node_suffix_length = kDefaultIndexSuffixLen;
    unsigned int distance = 0;
    unsigned int parallel_each = 1;
    unsigned int parallel_nodes = -1;  // if not set, redefined by |parallel|
    unsigned int num_bins_per_thread = 1;
    unsigned int parts_total = 1;
    unsigned int part_idx = 0;
    unsigned int suffix_len = 0;
    unsigned int frequency = 1;
    unsigned int alignment_length = 0;
    double memory_available = 1;
    unsigned int min_count = 1;
    unsigned int max_count = std::numeric_limits<unsigned int>::max();
    unsigned int min_value = 1;
    unsigned int max_value = std::numeric_limits<unsigned int>::max();
    unsigned int num_top_labels = -1;
    unsigned int genome_binsize_anno = 1000;
    unsigned int arity_brwt = 2;
    unsigned int relax_arity_brwt = 10;
    unsigned long long RA_ivbuffer_size = 16'384; // in B
    unsigned int min_tip_size = 1;
    unsigned int min_unitig_median_kmer_abundance = 1;
    int fallback_abundance_cutoff = 1;
    unsigned int port = 5555;
    unsigned int bloom_max_num_hash_functions = 10;
    unsigned int num_columns_cached = 10;
    unsigned int max_hull_forks = 4;
    unsigned int row_diff_stage = 0;

    // traversal (see docs/SPEC-labeled-traversal-core.md)
    bool traverse_resolve = false;          // run resolve/select instead of traversal
    // print the loader dependency inventory of -i / -a (index_load_inventory) and exit
    bool traverse_index_inventory = false;
    std::string index_release;              // echoed in responses; requests may pin it
    // the index identity of DESIGN-traverse-graphlet.md §3.1 (single-index mode): a
    // name for humans and routing, and the bundle manifest whose digest is the identity
    std::string index_name;
    std::string index_manifest;
    // Server-side caps for POST /traverse (0 = unlimited). They are deliberately
    // non-zero by default: a request that names NO labels makes the server derive the
    // permitted set, whose cost is set by the seed and the index rather than by
    // anything the caller declared, so an unconfigured deployment must not be the
    // unlimited one. The CLI, which is operator-run, applies no caps.
    double traverse_max_time_ms = 30'000;
    size_t traverse_max_seeds = 64;
    uint64_t traverse_max_seed_bp = 100'000;
    size_t traverse_max_seed_labels = 10'000;
    uint64_t resolve_max_query_bp = 0;
    // The server's maxima of a request's budgets (bounds.max_memory_mb, bounds.max_work_units;
    // owner decision R16): 0 = off (a request's budgets as given, an omitted one none). Set,
    // a larger budget is lowered to it and an omitted one set to it, echoed as clamped like
    // the time cap. Off by default: a budget changes how a walk stops, so a deployment that
    // did not choose one keeps the results it always gave
    uint64_t traverse_max_memory_mb = 0;
    uint64_t traverse_max_work_units = 0;
    // Ledger-managed /traverse attempts (requests with attempt_id, traverse_attempts.hpp):
    // the allowance added to n_seeds x the per-seed time budget in the duration bound the
    // server enforces, and how long (and how many) finished attempts stay queryable
    double traverse_attempt_allowance_ms = 10'000;
    uint64_t traverse_attempt_retention_s = 3600;
    size_t traverse_attempt_retention = 10'000;
    // The longest a tombstone is held to cover the not_after_ms a cancel (or a refused copy of
    // the request) names, plus the clock skew allowance (s; never less than retention_s, which
    // every tombstone is held at least): what bounds how long one cancel can keep an id
    // refused (review of pass 5, finding 1)
    uint64_t traverse_attempt_tombstone_max_s = 86'400;
    // the accepted ranges of the retention settings (refused at start-up beyond them): a year,
    // ten million attempts
    static constexpr uint64_t kMaxAttemptRetentionS = 31'536'000;
    static constexpr uint64_t kMaxAttemptRetentionCount = 10'000'000;
    // what a ledger adds to a request's not_after_ms before it treats an unanswered attempt as
    // one that cannot start subsequently (stated in the capabilities; the server's own check is
    // strict), ms
    uint64_t traverse_clock_skew_ms = 2000;
    // the expected duration of one uninterruptible piece of a /traverse annotation read under
    // a deadline (the reads are decoded in chunks sized from the observed per-row time, the
    // deadline checked between them); 0: one piece per read, as before (ms)
    uint64_t traverse_chunk_target_ms = 50;
    // the bound of a /traverse request's row-diff path cache (MiB; the rows its reads
    // reconstruct, kept so that a later read's row-diff path stops at a cached row): within
    // it without a memory budget, within the label cache's allotment under one; 0 = off
    uint64_t traverse_path_cache_mb = 128;
    // the zlib level of the traversal routes' compressed bodies (the other routes keep 9):
    // measured on real responses, level 1 writes 3-4 times faster than 9 at 1.8 times the
    // bytes, and the time to build a response is what a delivery window bounds
    int traverse_compression_level = 1;
    // the throughputs the delivery reserve of an attempt assumes (MB/s): compressing at
    // compression_level, and building the JSON text of a seed's result (replaced by the
    // slowest rate measured on the attempt's own seeds of 1 MB or more)
    double traverse_delivery_compress_mbps = 50;
    double traverse_delivery_build_mbps = 10;

    // POST /pattern and `metagraph pattern` (docs/DESIGN-pattern-search.md §5.3): each cap is a
    // request field's maximum (a larger request value is lowered to it and stated as clamped)
    // and, but for the time budget, its default; the CLI applies the same ones, so that both
    // answer alike
    double pattern_min_information_bits = 24;
    // owner decision #24 of 2026-10-08: on a graph without its dummy-edge mask, a pattern whose
    // unchecked candidates number at most this has each of them tested at query time (k - 1
    // steps each, charged to max_steps), its counts then exact; 0 tests none. Capped low (the
    // owner: per-query checking of large blocks is too expensive): at most
    // kMaxPatternCheckedEntries. Not a request field (stated in the capabilities' caps)
    uint64_t pattern_max_checked_entries = 50;
    static constexpr uint64_t kMaxPatternCheckedEntries = 1'000;
    uint64_t pattern_max_contexts = 10'000;
    uint64_t pattern_max_anchors = 1'000;
    // increment 4 (long_search "paths", §4.2): the default and maximum of max_paths, the
    // retrieval threshold on the completed paths of a pattern longer than k (all_or_count
    // releases them only when their exact count is at most this; partial's cap)
    uint64_t pattern_max_paths = 1'000;
    // the step cap of a request (§5.3). Not calibrated against the time budget on a deployed
    // index (review of 2026-10-07, X-EFFICIENCY-05): measured in RAM at 0.2-0.7 us a step
    // (the mini index, random graphs of 0.5 and 2 billion edges; M5 Max, shared), so 1e8 steps
    // take 20-75 s, about the default 60 s budget: which of the two stops a heavy request
    // first depends on the machine and its load (a time stop is time_limited). The rate on
    // a large mmapped index is unmeasured (the milestone-6 benchmark, DESIGN §13); a request
    // cannot raise this cap, only the operator can
    uint64_t pattern_max_steps = 100'000'000;
    // the time budget of a request that names none, and the most one may name (the owner,
    // 2026-10-07: 60 s by default, capped under the 900 s content timeout with room for the
    // answer's serialisation and compression)
    uint64_t pattern_default_time_ms = 60'000;
    uint64_t pattern_max_time_ms = 600'000;
    // the finalisation reserve inside the time budget: work stops at least this long before
    // the deadline so that the answer can still be written by it (ms); longer by the
    // estimated time to write what the answer buffers (review of 2026-10-07,
    // X-EFFICIENCY-04): at the rates below (MB/s), building and writing its JSON text and
    // compressing it, with a margin of 1.25 (pattern_retrieval.hpp, AnswerVolume). Starting
    // estimates as /traverse's delivery rates: conservative (16 x 10,000 results, 18.4 MB of
    // text, were written and gzipped in 0.3-0.45 s on an M-series Mac; the model gives 2.8 s)
    uint64_t pattern_finalize_ms = 250;
    double pattern_delivery_build_mbps = 10;
    double pattern_delivery_compress_mbps = 50;
    // patterns per request (a longer list is refused, not cut)
    uint64_t pattern_max_patterns = 16;
    // output.labels "all" (increment 3, §5.3): the labels kept per row (more are stated as a
    // truncated anchor), the annotation work per request (the oracle's units: 8 per row, 1
    // per entry and coordinate, and the rows' row-diff dependencies), the request's memory
    // account (MiB), and partial's lists: labels per pattern, occurrences per label
    uint64_t pattern_max_labels_per_anchor = 64;
    uint64_t pattern_max_annotation_work = 100'000'000;
    uint64_t pattern_max_memory_mb = 256;
    uint64_t pattern_max_labels = 1'000;
    uint64_t pattern_max_occurrences = 16;
    // increment 5b, a predicate's selection (SPEC-pattern-search.md §19.2, §19.3; owner
    // decision P3 of 2026-10-08): the default and maximum of max_predicate_contexts (the raw
    // contexts a pattern's selection may test) and of max_predicate_work (the selection's work
    // per request, the oracle's units, a budget of its own beside max_annotation_work), and the
    // names a predicate may list (not a request field: a larger predicate is refused,
    // predicate_too_large). The last is at most kMaxPatternPredicateLabels (refused at
    // start-up above): the bound predicate and its echo are linear in it
    uint64_t pattern_max_predicate_contexts = 100'000;
    uint64_t pattern_max_predicate_work = 100'000'000;
    uint64_t pattern_max_predicate_labels = 10'000;
    static constexpr uint64_t kMaxPatternPredicateLabels = 1'000'000;
    // server_query (-i / -a) and pattern: a succinct graph loaded without its .edgemask gets
    // the same dummy-edge mask built in memory before it is served (§4, mask: built_at_load),
    // for small indexes and tests; a large one is given the file once by transform --mask-dummy
    bool pattern_build_mask = false;
    // transform --mask-dummy: replace an existing .edgemask (refused without, since a mask
    // written by build or by another run would be overwritten unseen)
    bool force = false;

    unsigned int max_path_length = 100;
    unsigned int smoothing_window = 1;  // no smoothing by default
    unsigned int num_kmers_in_seq = 0;  // assume all input reads have this length

    unsigned long long int query_batch_size = 100'000'000;
    unsigned long long int num_rows_subsampled = 1'000'000;
    unsigned long long int num_singleton_kmers = 0;
    unsigned long long int max_hull_depth = -1;  // the default is a function of input
    unsigned long long int num_chars = 0;

    uint8_t count_width = 8;

    // Alignment options
    bool alignment_edit_distance = false;
    bool alignment_chain = false;
    bool alignment_post_chain = false;
    bool alignment_seed_complexity_filter = true;

    int8_t alignment_match_score = 2;
    int8_t alignment_mm_transition_score = 3;
    int8_t alignment_mm_transversion_score = 3;
    int8_t alignment_gap_opening_penalty = 6;
    int8_t alignment_gap_extension_penalty = 2;
    int8_t alignment_end_bonus = 5;

    int32_t alignment_min_path_score = 0;
    int32_t alignment_xdrop = 27;

    size_t alignment_num_alternative_paths = 1;
    size_t alignment_min_seed_length = 19;
    size_t alignment_max_seed_length = std::numeric_limits<size_t>::max();
    size_t alignment_max_num_seeds_per_locus = 1000;

    double alignment_rel_score_cutoff = 0.95;

    double discovery_fraction = 0.7;
    double presence_fraction = 0.0;
    double min_count_quantile = 0.0;
    double max_count_quantile = 1.0;
    double bloom_fpp = 1.0;
    double bloom_bpk = 4.0;
    double alignment_max_nodes_per_seq_char = 5.0;
    double alignment_max_ram = 200;
    // TODO: rename to min_covered_by_seeds
    double alignment_min_exact_match = 0.7;
    double min_fraction = 0.0;
    double max_fraction = 1.0;
    double cleaning_threshold_percentile = 0.001;
    std::vector<double> count_slice_quantiles;
    std::vector<double> count_quantiles;

    std::vector<std::string> fnames;
    std::vector<std::string> anno_labels;
    std::vector<std::string> infbase_annotators;
    std::string outfbase;
    std::string infbase;
    std::string rename_instructions_file;
    std::string refpath;
    std::string suffix;
    std::string fasta_header_delimiter;
    std::string anno_labels_delimiter = ":";
    std::string fasta_anno_comment_delim = UNINITIALIZED_STR;
    std::string header = "";
    std::string host_address;
    std::string assembly_config_file;
    std::string linkage_file;
    std::string intersected_columns;

    std::filesystem::path tmp_dir;

    size_t disk_cap_bytes = -1;

    enum IdentityType {
        NO_IDENTITY = -1,
        BUILD = 1,
        CLEAN,
        EXTEND,
        MERGE,
        CONCATENATE,
        COMPARE,
        ALIGN,
        STATS,
        ANNOTATE,
        MERGE_ANNOTATIONS,
        TRANSFORM,
        TRANSFORM_ANNOTATION,
        ASSEMBLE,
        RELAX_BRWT,
        QUERY,
        SERVER_QUERY,
        TRAVERSE,
        PATTERN,
    };
    IdentityType identity = NO_IDENTITY;

    graph::boss::BOSS::State state = graph::boss::BOSS::State::STAT;

    static std::string state_to_string(graph::boss::BOSS::State state);
    static graph::boss::BOSS::State string_to_state(const std::string &string);

    enum AnnotationType {
        ColumnCompressed = 1,
        RowCompressed,
        BRWT,
        BinRelWT,
        RowDiff,
        RowDiffBRWT,
        RowDiffRowFlat,
        RowDiffRowSparse,
        RowDiffDisk,
        RowFlat,
        RowSparse,
        RBFish,
        RbBRWT,
        IntBRWT,
        IntRowDiffBRWT,
        IntRowDiffDisk,
        ColumnCoord,
        BRWTCoord,
        RowDiffCoord,
        RowDiffBRWTCoord,
        RowDiffDiskCoord,
    };

    enum GraphType {
        INVALID = -1,
        SUCCINCT = 1,
        HASH,
        HASH_PACKED,
        HASH_STR,
        HASH_FAST,
        SSHASH,
        BITMAP,
    };

    AnnotationType anno_type = ColumnCompressed;
    static std::string annotype_to_string(AnnotationType state);
    static AnnotationType string_to_annotype(const std::string &string);

    GraphType graph_type = SUCCINCT;
    static GraphType string_to_graphtype(const std::string &string);

    graph::DeBruijnGraph::Mode graph_mode = graph::DeBruijnGraph::BASIC;
    static std::string graphmode_to_string(graph::DeBruijnGraph::Mode mode);
    static graph::DeBruijnGraph::Mode string_to_graphmode(const std::string &string);

    QueryMode query_mode = MATCHES;
    static std::string querymode_to_string(QueryMode mode);
    static QueryMode string_to_querymode(const std::string &string);

    void print_usage(const std::string &prog_name,
                     IdentityType identity = NO_IDENTITY);
};

} // namespace cli
} // namespace mtg

#endif // __CONFIG_HPP__
