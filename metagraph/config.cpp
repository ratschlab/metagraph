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
    int i = 1;
    if (argc == 1) {
        print_usage(std::string(argv[0]));
        exit(-1);
    }
    while (i < argc) {
        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
            verbose = true;
        } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--integrate")) {
            integrate = true;
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
        } else if (!strcmp(argv[i], "-m") || !strcmp(argv[i], "--mergelist")) {
            merge = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--compare")) {
            compare = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-a") || !strcmp(argv[i], "--align")) {
            align = std::string(argv[++i]);
        //} else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--threads")) {
        //    num_threads = atoi(argv[++i]);
        //} else if (!strcmp(argv[i], "--debug")) {
        //    debug = true;
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(std::string(argv[0]));
            exit(0);
        } else {
            if (argv[i][0] == '-') {
                fprintf(stderr, "\nERROR: Unknown option %s\n", argv[i]);
                print_usage(std::string(argv[0]));
                exit(-1);
            }
            else {
                fname.push_back(std::string(argv[i]));
            }
        }
        i++;
    }
}

// nothing to destroy
Config::~Config() {}

void Config::print_usage(std::string prog_name) {
    fprintf(stderr, "A comprehensive graph represenatation of metagenome information\nVersion 0.1\n\n");
    fprintf(stderr, "This program is the first implementation of a\n");
    fprintf(stderr, "meta-metagenome graph for identifciation and annotation\n");
    fprintf(stderr, "purposes.\n\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "\t%s -O OUTBASE [options] FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());
    fprintf(stderr, "Available Options:\n\n");
    // Input and threads
    //fprintf(stderr, "\n\tInput handling and paralellization:\n");
    fprintf(stderr, "\t-S --sql-base [STR] \tbasename for SQL output file\n");
    fprintf(stderr, "\t-I --infile-base [STR] \tbasename for loading graph input files\n");
    fprintf(stderr, "\t-i --integrate \t\tintegrates fasta into given (-I) graph [off]\n");
    fprintf(stderr, "\t-k --kmer-length [INT] \tLength of the k-mer to use [3]\n");
    // Filter options
    //fprintf(stderr, "\n\tInput file filtering:\n");
    fprintf(stderr, "\t-d --distance [INT] \tMax allowed alignment distance [0]\n");
    fprintf(stderr, "\t-m --merge \t\tlist of graph file basenames used as input for merging [off]]\n");
    fprintf(stderr, "\t-c --compare \t\tlist of graph file basenames used as input for merging [off]\n");
    fprintf(stderr, "\t-a --align\t\trun in alignment mode [off]\n");
    // Output Options
    fprintf(stderr, "\n\tGeneral:\n");
    fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
    fprintf(stderr, "\t-h --help \t\tprint usage info\n");
    fprintf(stderr, "\n");
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
    integrate = false;
    distance = 0;
    k = 3;
}

