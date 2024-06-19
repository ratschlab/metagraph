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
    bool inplace = false;
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

    QueryMode query_mode = LABELS;
    static std::string querymode_to_string(QueryMode mode);
    static QueryMode string_to_querymode(const std::string &string);

    void print_usage(const std::string &prog_name,
                     IdentityType identity = NO_IDENTITY);
};

} // namespace cli
} // namespace mtg

#endif // __CONFIG_HPP__
