#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <string>
#include <vector>


class Config {
  public:
    Config(int argc, const char *argv[]);

    bool verbose = false;
    bool print_graph = false;
    bool print_graph_internal_repr = false;
    bool print_column_names = false;
    bool reverse = false;
    bool canonical = false;
    bool complete = false;
    bool dynamic = false;
    bool mark_dummy_kmers = true;
    bool filename_anno = false;
    bool fasta_anno = false;
    bool to_adj_list = false;
    bool unitigs = false;
    bool count_kmers = false;
    bool query_presence = false;
    bool filter_present = false;
    bool dump_raw_anno = false;
    bool sparse = false;
    bool fast = false;
    bool count_labels = false;
    bool suppress_unlabeled = false;
    bool generate_filtered_fasta = false;
    bool generate_filtered_fastq = false;
    bool internal = false;
    bool clear_dummy = false;
    bool count_dummy = false;
    bool canonical_mode = false;
    bool greedy_brwt = false;
    bool separately = false;
    bool call_bubbles = false;

    unsigned int k = 3;
    unsigned int filter_k = 3;
    unsigned int distance = 0;
    unsigned int parallel = 1;
    unsigned int num_bins_per_thread = 1;
    unsigned int parts_total = 1;
    unsigned int part_idx = 0;
    unsigned int suffix_len = 0;
    unsigned int frequency = 1;
    unsigned int alignment_length = 0;
    unsigned int alignment_num_top_paths = 10;
    unsigned int memory_available = 0;
    unsigned int min_count = 1;
    unsigned int max_count = std::numeric_limits<unsigned int>::max();
    unsigned int unreliable_kmers_threshold = 0;
    unsigned int num_top_labels = -1;
    unsigned int genome_binsize_anno = 1000;
    unsigned int kmer_mapping_mode = 0;
    unsigned int arity_brwt = 2;
    unsigned int relax_arity_brwt = 10;
    unsigned int pruned_dead_end_size = 0;
    unsigned int port = 5555;

    double discovery_fraction = 1.0;
    double label_mask_out_fraction = 0.0;
    double alignment_sw_threshold = 0.8;

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
    std::string annotation_label = "";
    std::string header = "";
    std::string accession2taxid;
    std::string taxonomy_nodes;
    std::string taxonomy_map;

    enum IdentityType {
        NO_IDENTITY = -1,
        BUILD = 1,
        EXTEND,
        FILTER,
        FILTER_STATS,
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
        HASH_STR,
        BITMAP,
    };

    AnnotationType anno_type = ColumnCompressed;
    GraphType graph_type = SUCCINCT;

    static std::string annotype_to_string(AnnotationType state);
    static AnnotationType string_to_annotype(const std::string &string);
    static GraphType string_to_graphtype(const std::string &string);

    void print_usage(const std::string &prog_name,
                     IdentityType identity = NO_IDENTITY);
};

#endif // __CONFIG_HPP__
