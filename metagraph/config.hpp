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
    bool fast = false;
    bool fasta_anno = false;
    bool to_adj_list = false;
    bool count_kmers_query = false;
    bool query_presence = false;
    bool traversal_merge = false;
    bool dump_raw_anno = false;
    bool use_row_annotator = false;
    bool count_labels = false;

    unsigned int k = 3;
    unsigned int distance = 0;
    unsigned int parallel = 1;
    unsigned int num_bins_per_thread = 1;
    unsigned int parts_total = 1;
    unsigned int part_idx = 0;
    unsigned int collect = 1;
    unsigned int frequency = 1;
    unsigned int nsplits = 1;
    unsigned int alignment_length = 0;
    unsigned int memory_available = 0;
    unsigned int bloom_num_hash_functions = 0;
    unsigned int bloom_test_num_kmers = 0;
    unsigned int noise_kmer_frequency = 0;
    unsigned int num_top_labels = -1;

    double bloom_fpp = -1;
    double bloom_bits_per_edge = -1;
    double discovery_fraction = 1.0;

    std::vector<std::string> fname;
    std::string outfbase;
    std::string infbase;
    std::string sqlfbase;
    std::string dbpath;
    std::string refpath;
    std::string suffix;
    std::string fasta_header_delimiter;
    std::string anno_labels_delimiter = ":";

    enum IdentityType {
        NO_IDENTITY = -1,
        BUILD = 1,
        FILTER,
        EXPERIMENT,
        MERGE,
        COMPARE,
        ALIGN,
        STATS,
        ANNOTATE,
        ANNOTATE_BLOOM,
        TRANSFORM,
        CLASSIFY
    };
    IdentityType identity = NO_IDENTITY;

    enum StateType { STAT = 1, DYN };
    StateType state = STAT;

    void print_usage(const std::string &prog_name,
                     IdentityType identity = NO_IDENTITY);
};

#endif // __CONFIG_HPP__
