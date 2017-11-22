#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <string>
#include <vector>


class Config {
  public:
    Config(int argc, const char *argv[]);

    bool verbose = false;
    bool quiet = false;
    bool print_graph = false;
    bool print_graph_succ = false;
    bool reverse = false;
    bool fast = false;
    bool add_anno = false;

    unsigned int k = 3;
    unsigned int distance = 0;
    unsigned int parallel = 1;
    unsigned int bins_per_thread = 1;
    unsigned int parts_total = 1;
    unsigned int part_idx = 0;
    unsigned int collect = 1;
    unsigned int frequency = 1;
    unsigned int nsplits = 1;

    std::vector<std::string> fname;
    std::string outfbase;
    std::string infbase;
    std::string sqlfbase;
    std::string dbpath;
    std::string refpath;
    std::string suffix;

    enum IdentityType {
        NO_IDENTITY = -1,
        BUILD = 1,
        MERGE,
        COMPARE,
        ALIGN,
        STATS,
        ANNOTATE,
        DUMP,
        CLASSIFY
    };
    IdentityType identity = NO_IDENTITY;

    enum StateType { STAT = 1, DYN, CSTR };
    StateType state = STAT;

    void print_usage(const std::string &prog_name, IdentityType identity = NO_IDENTITY);
};

#endif // __CONFIG_HPP__
