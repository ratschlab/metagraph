//#include <cstdlib>
//#include <cstdio>
//#include <cstring>
#include "config.hpp"

// default constructor
Config::Config() {
    init();
}

Config::Config(int argc, const char *argv[]) {

    init();
    // provide help overview if no identity was given
    if (argc == 1) {
        print_usage(std::string(argv[0]));
        exit(-1);
    }

    // parse identity from first command line argument
    if (!strcmp(argv[1], "merge")) {
        identity = merge;
    } else if (!strcmp(argv[1], "compare")) {
        identity = compare;
    } else if (!strcmp(argv[1], "align")) {
        identity = align;
    } else if (!strcmp(argv[1], "build")) {
        identity = build;
    } else if (!strcmp(argv[1], "stats")) {
        identity = stats;
    }
    // provide help screen for chosen identity
    if (argc == 2) {
        print_usage(std::string(argv[0]), identity);
        exit(-1);
    }
    // parse remaining command line items
    int i = 2;
    while (i < argc) {
        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
            verbose = true;
        } else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--print-graph")) {
            print_graph = true;
        } else if (!strcmp(argv[i], "-k") || !strcmp(argv[i], "--kmer-length")) {
            k = atoi(argv[i++]);
        } else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--distance")) {
            distance = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-O") || !strcmp(argv[i], "--outfile-base")) {
            outfbase = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-S") || !strcmp(argv[i], "--sql-base")) {
            sqlfbase = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-I") || !strcmp(argv[i], "--infile-base")) {
            infbase = std::string(argv[++i]);
        //} else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--threads")) {
        //    num_threads = atoi(argv[++i]);
        //} else if (!strcmp(argv[i], "--debug")) {
        //    debug = true;
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(std::string(argv[0]), identity);
            exit(0);
        } else {
            if (argv[i][0] == '-') {
                fprintf(stderr, "\nERROR: Unknown option %s\n", argv[i]);
                print_usage(std::string(argv[0], identity));
                exit(-1);
            }
            else {
                fname.push_back(std::string(argv[i]));
            }
        }
        i++;
    }
}

// nothing to destro
Config::~Config() {}

void Config::print_usage(std::string prog_name, int identity) {
    fprintf(stderr, "Comprehensive metagenome graph representation -- Version 0.1\n\n");
    //fprintf(stderr, "This program is the first implementation of a\n");
    //fprintf(stderr, "meta-metagenome graph for identification and annotation\n");
    //fprintf(stderr, "purposes.\n\n");

    switch (identity) {
        case noidentity: {
            fprintf(stderr, "Usage: %s <command> [command specific options]\n\n", prog_name.c_str());
            fprintf(stderr, "Available commands:\n");
            fprintf(stderr, "\tbuild\t\tconstruct a graph object from input sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats or integrate sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats into a given graph\n\n");

            fprintf(stderr, "\tmerge\t\tintegrate a given set of graph structures\n");
            fprintf(stderr, "\t\t\tand output a new graph structure\n\n");

            fprintf(stderr, "\tcompare\t\tcheck whether two given graphs are identical\n\n");

            fprintf(stderr, "\talign\t\talign the reads provided in files in fast[a|q]\n");
            fprintf(stderr, "\t\t\tformats to the graph\n\n");

            fprintf(stderr, "\tstats\t\tprint graph statistics for given graph(s)\n\n");
        } break;
        case build: {
            fprintf(stderr, "Usage: %s build [options] FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());
            fprintf(stderr, "Available options for build:\n");
            fprintf(stderr, "\t-O --outfile-base [STR] \tbasename of output file []\n");
            fprintf(stderr, "\t-S --sql-base [STR] \tbasename for SQL output file\n");
            fprintf(stderr, "\t-I --infile-base [STR] \tbasename for loading graph input file\n");
            fprintf(stderr, "\t-k --kmer-length [INT] \tlength of the k-mer to use [3]\n");
            fprintf(stderr, "\t-p --print-graph \tprint graph table to the screen [off]\n");
        } break;
        case align: {
            fprintf(stderr, "Usage: %s align [options] FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());
            fprintf(stderr, "Available options for align:\n");
            fprintf(stderr, "\t-d --distance [INT] \tMax allowed alignment distance [0]\n");
        } break;
        case compare: {
            fprintf(stderr, "Usage: %s compare [options] GRAPH1 [[GRAPH2] ...]\n\n", prog_name.c_str());
            fprintf(stderr, "Available options for compare:\n");
            fprintf(stderr, "\t-I --infile-base [STR] \tbasename for loading graph input file\n");
        } break;
        case merge: {
            fprintf(stderr, "Usage: %s merge [options] GRAPH1 [[GRAPH2] ...]\n\n", prog_name.c_str());
            fprintf(stderr, "Available options for merge:\n");
            fprintf(stderr, "\t-O --outfile-base [STR] \tbasename of output file []\n");
            fprintf(stderr, "\t-p --print-graph \tprint graph table to the screen [off]\n");
        } break;
        case stats: {
            fprintf(stderr, "Usage: %s stats [options] GRAPH1 [[GRAPH2] ...]\n\n", prog_name.c_str());
            fprintf(stderr, "Available options for stats:\n");
            fprintf(stderr, "\t-I --infile-base [STR] \tbasename for loading graph input file\n");
            fprintf(stderr, "\t-O --outfile-base [STR] \tbasename of output file []\n");
            fprintf(stderr, "\t-p --print-graph \tprint graph table to the screen [off]\n");
        } break;

    }
    if (identity != noidentity) {
        fprintf(stderr, "\n\tGeneral options:\n");
        fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
        fprintf(stderr, "\t-h --help \t\tprint usage info\n");
        fprintf(stderr, "\n");
    }
}

/*void Config::print_call(std::string prog_name) {
    fprintf(stdout, "%s has been started with the following parameters:\n\n", prog_name.c_str());
    fprintf(stdout, "\t input file:           %s\n", infile.c_str()); 
    fprintf(stdout, "\t output file:          %s\n", outfile.c_str()); 
    fprintf(stdout, "\t strand specific:      %s\n", strand_specific?"yes":"no");
    fprintf(stdout, "\t pre filter:           %s\n", pre_filter?"on":"off");
    fprintf(stdout, "\t init on secondary:    %s\n", take_non_secondary_only?"off":"on");
    if (pre_filter) {
        fprintf(stdout, "\t filter dist:          %i\n", filter_distance);
        fprintf(stdout, "\t use variants:         %s\n", use_variants?"on":"off");
    }
    fprintf(stdout, "\t max list length:      %i\n", max_list_length);
    fprintf(stdout, "\t pair usage:           %s\n", use_pair_info?"on":"off");
    if (use_pair_info) {
        fprintf(stdout, "\t max frag size size:   %i\n", max_gen_frag_size);
        fprintf(stdout, "\t max pair list length: %i\n", max_pair_list_length);
    }
    if (trim_id > 0)
        fprintf(stdout, "\t trim read id by:      %i\n", trim_id);
    fprintf(stdout, "\t print best only:      %s\n", print_best_only?"on":"off");
    fprintf(stdout, "\t print unmapped:       %s\n", print_unmapped?"on":"off");
    fprintf(stdout, "\t iterations:           %i\n", iterations);
    fprintf(stdout, "\t 1 iteration burn in:  %s\n", burn_in?"on":"off");
    if (use_brkpts) 
        fprintf(stdout, "\t annotation file:      %s\n", annotation.c_str()); 
    fprintf(stdout, "\t threads:              %i\n", num_threads);
    if (use_mip_variance || ! use_mip_objective) {
        fprintf(stdout, "\t window size:          %i\n", window_size);
    } else {
        fprintf(stdout, "\t use MiTie objective:    %s\n", use_mip_objective?"on":"off");
        fprintf(stdout, "\t segment file:         %s\n", segmentfile.c_str()); 
        fprintf(stdout, "\t loss file:            %s\n", lossfile.c_str()); 
        fprintf(stdout, "\t read length:          %i\n", read_len); 
        fprintf(stdout, "\t zero for unpred seg:  %s\n", zero_unpred?"yes":"no");
        fprintf(stdout, "\t use variance if no MiTie-segment overlaps: %s\n", use_mip_variance?"on":"off");
    }
}*/

// PRIVATE
void Config::init() {
    verbose = false;
    print_graph = false;
    distance = 0;
    k = 3;
    identity = noidentity;
}

