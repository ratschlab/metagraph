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
