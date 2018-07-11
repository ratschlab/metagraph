#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <string>
#include <vector>


class Config {
  public:
    Config(int argc, const char *argv[]);

    bool verbose = false;
    bool quiet = false;
    bool print_graph_succ = false;
    bool reverse = false;
    bool dynamic = false;
    bool fast = false;
    bool filename_anno = false;
    bool fasta_anno = false;
    bool to_adj_list = false;
    bool to_sequences = false;
    bool to_row_annotator = false;
    bool count_kmers_query = false;
    bool query_presence = false;
    bool filter_present = false;
    bool dump_raw_anno = false;
    bool use_row_annotator = false;
    bool sparse = false;
    bool count_labels = false;
    bool suppress_unlabeled = false;
    bool generate_filtered_fasta = false;
    bool generate_filtered_fastq = false;
    bool use_kmc = false;

    unsigned int k = 3;
    unsigned int filter_k = 3;
    unsigned int distance = 0;
    unsigned int parallel = 1;
    unsigned int num_bins_per_thread = 1;
    unsigned int parts_total = 1;
    unsigned int part_idx = 0;
    unsigned int suffix_len = 0;
    unsigned int frequency = 1;
    unsigned int nsplits = 1;
    unsigned int alignment_length = 0;
    unsigned int memory_available = 0;
    unsigned int bloom_num_hash_functions = 0;
    unsigned int bloom_test_num_kmers = 0;
    unsigned int max_unreliable_abundance = 0;
    unsigned int unreliable_kmers_threshold = 0;
    unsigned int num_top_labels = -1;
    unsigned int genome_binsize_anno = 1000;

    double bloom_fpp = -1;
    double bloom_bits_per_edge = -1;
    double discovery_fraction = 1.0;

    std::vector<std::string> fname;
    std::vector<std::string> anno_labels;
    std::vector<std::string> infbase_annotators;
    std::string outfbase;
    std::string infbase;
    std::string sqlfbase;
    std::string dbpath;
    std::string refpath;
    std::string suffix;
    std::string fasta_header_delimiter;
    std::string anno_labels_delimiter = ":";
    std::string annotation_label = "";

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
        ANNOTATE_BLOOM,
        ANNOTATE_COORDINATES,
        MERGE_ANNOTATIONS,
        TRANSFORM,
        TRANSFORM_ANNOTATION,
        CLASSIFY
    };
    IdentityType identity = NO_IDENTITY;

    enum StateType { STAT = 1, DYN, SMALL };
    StateType state = STAT;

    static std::string state_to_string(StateType state);
    static StateType string_to_state(const std::string &string);

    void print_usage(const std::string &prog_name,
                     IdentityType identity = NO_IDENTITY);
};

#endif // __CONFIG_HPP__
