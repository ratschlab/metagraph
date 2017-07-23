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
    } else if (!strcmp(argv[1], "annotate")) {
        identity = annotate;
    } else if (!strcmp(argv[1], "classify")) {
        identity = classify;
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
        } else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--quiet")) {
            quiet = true;
        } else if (!strcmp(argv[i], "-P") || !strcmp(argv[i], "--print")) {
            print_graph = true;
        } else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--reverse")) {
            reverse = true;
        } else if (!strcmp(argv[i], "--fast")) {
            fast = true;
        } else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--parallel")) {
            parallel = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--parts-total")) {
            parts_total = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--part-idx")) {
            part_idx = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bins-per-thread")) {
            bins_per_thread = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-k") || !strcmp(argv[i], "--kmer-length")) {
            k = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "--frequency")) {
            frequency = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--distance")) {
            distance = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-O") || !strcmp(argv[i], "--outfile-base")) {
            outfbase = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-R") || !strcmp(argv[i], "--reference")) {
            refpath = std::string(argv[++i]);
        //} else if (!strcmp(argv[i], "-D") || !strcmp(argv[i], "--db-path")) {
        //    dbpath = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-S") || !strcmp(argv[i], "--sql-base")) {
            sqlfbase = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-I") || !strcmp(argv[i], "--infile-base")) {
            infbase = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-C") || !strcmp(argv[i], "--collect")) {
            collect = atoi(argv[++i]);
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

            fprintf(stderr, "\tannotate\tgiven a graph and a fast[a|q] file, annotate\n");
            fprintf(stderr, "\t\t\tthe respective kmers\n\n");
        } break;
        case build: {
            fprintf(stderr, "Usage: %s build [options] FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());
            fprintf(stderr, "Available options for build:\n");
            fprintf(stderr, "\t-O --outfile-base [STR] \tbasename of output file []\n");
            fprintf(stderr, "\t-S --sql-base [STR] \tbasename for SQL output file\n");
            fprintf(stderr, "\t-I --infile-base [STR] \tbasename for loading graph input file\n");
            fprintf(stderr, "\t-k --kmer-length [INT] \tlength of the k-mer to use [3]\n");
            fprintf(stderr, "\t-r --reverse \tbuild graph from reverse complement of input [off]\n");
            fprintf(stderr, "\t-P --print \tprint graph table to the screen [off]\n");
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
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
            fprintf(stderr, "\t-b --bins-per-thread [INT] \tnumber of bins each thread computes on average [1]\n");
            fprintf(stderr, "\t-P --print \t\tprint graph table to the screen [off]\n");
            fprintf(stderr, "\t   --part-idx [INT] \t\tidx to use when doing external merge []\n");
            fprintf(stderr, "\t   --parts-total [INT] \t\ttotal number of parts in external merge[]\n");
            fprintf(stderr, "\t-C --collect [INT] \t\tinitiate collection of external merge, provide total number of splits [1]\n");
        } break;
        case stats: {
            fprintf(stderr, "Usage: %s stats [options] GRAPH1 [[GRAPH2] ...]\n\n", prog_name.c_str());
            fprintf(stderr, "Available options for stats:\n");
            fprintf(stderr, "\t-O --outfile-base [STR] \tbasename of output file []\n");
            fprintf(stderr, "\t-P --print \tprint graph table to the screen [off]\n");
        } break;
        case annotate: {
            fprintf(stderr, "Usage: %s annotate [options] PATH1 [[PATH2] ...]\n\tEach path is given as file in fasta or fastq format.\n\n", prog_name.c_str());
            fprintf(stderr, "Available options for annotate:\n");
            //fprintf(stderr, "\t-D --db-path \tpath that is used to store the annotations database []\n");
            fprintf(stderr, "\t-I --infile-base [STR] \tbasename for loading graph to be annotated\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
            fprintf(stderr, "\t-b --bins-per-thread [INT] \tnumber of bins each thread computes on average [1]\n");
            fprintf(stderr, "\t-f --frequency [INT] \t\twhen a, annotate only every a-th kmer [1]\n");
        } break;
        case classify: {
            fprintf(stderr, "Usage: %s classify [options] FILE1 [[FILE2] ...]\n\tEach read file is given in fasta or fastq format.\n\n", prog_name.c_str());
            fprintf(stderr, "\t-I --infile-base [STR] \tbasename for graph with annotation used for classifying\n");
            fprintf(stderr, "\t-d --distance [INT] \tMax allowed alignment distance [0]\n");
        } break;

    }
    if (identity != noidentity) {
        fprintf(stderr, "\n\tGeneral options:\n");
        fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
        fprintf(stderr, "\t-q --quiet \t\tproduce as little log output as posible [off]\n");
        fprintf(stderr, "\t-h --help \t\tprint usage info\n");
        fprintf(stderr, "\n");
    }
}

// PRIVATE
void Config::init() {
    verbose = false;
    quiet = false;
    print_graph = false;
    reverse = false;
    fast = false;
    distance = 0;
    parallel = 1;
    bins_per_thread = 1;
    parts_total = 1;
    part_idx = 0;
    collect = 1;
    frequency = 1;
    k = 3;
    identity = noidentity;
}

