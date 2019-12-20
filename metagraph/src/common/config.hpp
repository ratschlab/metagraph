#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <string>
#include <vector>

#include "kmer/kmer_collector_config.hpp"


class Config {
  public:
    Config(int argc, char *argv[]);

    static constexpr auto UNINITIALIZED_STR = "\0";

    bool print_graph = false;
    bool print_graph_internal_repr = false;
    bool print_column_names = false;
    bool forward_and_reverse = false;
    bool canonical = false;
    bool complete = false;
    bool dynamic = false;
    bool mark_dummy_kmers = true;
    bool filename_anno = false;
    bool annotate_sequence_headers = false;
    bool to_adj_list = false;
    bool to_fasta = false;
    bool to_gfa = false;
    bool unitigs = false;
    bool kmers_in_single_form = false;
    bool initialize_bloom = false;
    bool count_kmers = false;
    bool query_presence = false;
    bool filter_present = false;
    bool dump_text_anno = false;
    bool sparse = false;
    bool fast = false;
    bool count_labels = false;
    bool suppress_unlabeled = false;
    bool clear_dummy = false;
    bool count_dummy = false;
    bool canonical_mode = false;
    bool greedy_brwt = false;
    bool separately = false;
    bool files_sequentially = false;
    bool call_bubbles = false;
    bool call_breakpoints = false;
    bool map_sequences = false;
    bool align_sequences = false;
    bool filter_by_kmer = false;
    bool output_json = false;

    unsigned int k = 3;
    unsigned int distance = 0;
    unsigned int parallel_nodes = -1;  // if not set, redefined by |parallel|
    unsigned int num_bins_per_thread = 1;
    unsigned int parts_total = 1;
    unsigned int part_idx = 0;
    unsigned int suffix_len = 0;
    unsigned int frequency = 1;
    unsigned int alignment_length = 0;
    unsigned int memory_available = 1;
    unsigned int min_count = 1;
    unsigned int max_count = std::numeric_limits<unsigned int>::max();
    unsigned int num_top_labels = -1;
    unsigned int genome_binsize_anno = 1000;
    unsigned int arity_brwt = 2;
    unsigned int relax_arity_brwt = 10;
    unsigned int row_cache_size = 0;
    unsigned int min_tip_size = 1;
    unsigned int min_unitig_median_kmer_abundance = 1;
    unsigned int fallback_abundance_cutoff = 1;
    unsigned int port = 5555;
    unsigned int bloom_max_num_hash_functions = 10;

    // Alignment options
    bool alignment_seed_unimems = false;
    bool alignment_edit_distance = false;

    int8_t alignment_match_score = 2;
    int8_t alignment_mm_transition = 1;
    int8_t alignment_mm_transversion = 2;
    int8_t alignment_gap_opening_penalty = 3;
    int8_t alignment_gap_extension_penalty = 1;

    int32_t alignment_min_cell_score = 0;
    int32_t alignment_min_path_score = 0;

    size_t alignment_queue_size = 20;
    size_t alignment_vertical_bandwidth = 17;
    size_t alignment_num_alternative_paths = 1;
    size_t alignment_min_seed_length = 0;
    size_t alignment_max_seed_length = std::numeric_limits<size_t>::max();
    size_t alignment_max_num_seeds_per_locus = std::numeric_limits<size_t>::max();

    double discovery_fraction = 1.0;
    double label_mask_in_fraction = 1.0;
    double label_mask_out_fraction = 0.0;
    double label_other_fraction = 1.0;
    double min_count_quantile = 0.;
    double max_count_quantile = 1.;
    double bloom_fpp = 1.0;
    double bloom_bpk = 4.0;
    std::vector<double> count_slice_quantiles;

    std::vector<std::string> fname;
    std::vector<std::string> anno_labels;
    std::vector<std::string> infbase_annotators;
    std::vector<std::string> label_mask_in;
    std::vector<std::string> label_mask_out;
    std::vector<std::string> label_filter;
    std::string outfbase;
    std::string infbase;
    std::string rename_instructions_file;
    std::string dbpath;
    std::string refpath;
    std::string suffix;
    std::string fasta_header_delimiter;
    std::string anno_labels_delimiter = ":";
    std::string fasta_anno_comment_delim = UNINITIALIZED_STR;
    std::string annotation_label = "";
    std::string header = "";
    std::string accession2taxid;
    std::string taxonomy_nodes;
    std::string taxonomy_map;

    enum IdentityType {
        NO_IDENTITY = -1,
        BUILD = 1,
        CLEAN,
        EXTEND,
        EXPERIMENT,
        MERGE,
        CONCATENATE,
        COMPARE,
        ALIGN,
        STATS,
        ANNOTATE,
        ANNOTATE_COORDINATES,
        MERGE_ANNOTATIONS,
        TRANSFORM,
        TRANSFORM_ANNOTATION,
        ASSEMBLE,
        RELAX_BRWT,
        QUERY,
        SERVER_QUERY,
        CALL_VARIANTS,
        PARSE_TAXONOMY
    };
    IdentityType identity = NO_IDENTITY;

    enum StateType { STAT = 1, DYN, SMALL, FAST };
    StateType state = STAT;

    static std::string state_to_string(StateType state);
    static StateType string_to_state(const std::string &string);

    enum AnnotationType {
        ColumnCompressed = 1,
        RowCompressed,
        BRWT,
        BinRelWT_sdsl,
        BinRelWT,
        RowFlat,
        RBFish,
    };

    enum GraphType {
        INVALID = -1,
        SUCCINCT = 1,
        HASH,
        HASH_PACKED,
        HASH_STR,
        HASH_FAST,
        BITMAP,
    };

    AnnotationType anno_type = ColumnCompressed;
    GraphType graph_type = SUCCINCT;

    mg::kmer::ContainerType container;
    static mg::kmer::ContainerType string_to_container(const std::string &string);

    static std::string annotype_to_string(AnnotationType state);
    static AnnotationType string_to_annotype(const std::string &string);
    static GraphType string_to_graphtype(const std::string &string);

    void print_usage(const std::string &prog_name,
                     IdentityType identity = NO_IDENTITY);
};

#endif // __CONFIG_HPP__
