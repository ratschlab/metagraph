#include "config.hpp"

#include <cstring>
#include <iostream>

#include "utils.hpp"


Config::Config(int argc, const char *argv[]) {
    // provide help overview if no identity was given
    if (argc == 1) {
        print_usage(argv[0]);
        exit(-1);
    }

    // parse identity from first command line argument
    if (!strcmp(argv[1], "merge")) {
        identity = MERGE;
    } else if (!strcmp(argv[1], "extend")) {
        identity = EXTEND;
    } else if (!strcmp(argv[1], "concatenate")) {
        identity = CONCATENATE;
    } else if (!strcmp(argv[1], "compare")) {
        identity = COMPARE;
    } else if (!strcmp(argv[1], "align")) {
        identity = ALIGN;
    } else if (!strcmp(argv[1], "build")) {
        identity = BUILD;
    } else if (!strcmp(argv[1], "filter")) {
        identity = FILTER;
    } else if (!strcmp(argv[1], "filter_stats")) {
        identity = FILTER_STATS;
    } else if (!strcmp(argv[1], "experiment")) {
        identity = EXPERIMENT;
    } else if (!strcmp(argv[1], "stats")) {
        identity = STATS;
    } else if (!strcmp(argv[1], "annotate")) {
        identity = ANNOTATE;
    } else if (!strcmp(argv[1], "bloom")) {
        identity = ANNOTATE_BLOOM;
    } else if (!strcmp(argv[1], "coordinate")) {
        identity = ANNOTATE_COORDINATES;
    } else if (!strcmp(argv[1], "merge_anno")) {
        identity = MERGE_ANNOTATIONS;
    } else if (!strcmp(argv[1], "classify")) {
        identity = CLASSIFY;
    } else if (!strcmp(argv[1], "transform")) {
        identity = TRANSFORM;
    } else if (!strcmp(argv[1], "transform_anno")) {
        identity = TRANSFORM_ANNOTATION;
    } else {
        print_usage(argv[0]);
        exit(-1);
    }

    // provide help screen for chosen identity
    if (argc == 2) {
        print_usage(argv[0], identity);
        exit(-1);
    }

    // parse remaining command line items
    for (int i = 2; i < argc; ++i) {
        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
            verbose = true;
        } else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--quiet")) {
            quiet = true;
        } else if (!strcmp(argv[i], "--print")) {
            print_graph_succ = true;
        } else if (!strcmp(argv[i], "--print-internal")) {
            print_graph_internal_repr = true;
        } else if (!strcmp(argv[i], "--count-kmers")) {
            count_kmers_query = true;
        } else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--reverse")) {
            reverse = true;
        } else if (!strcmp(argv[i], "--fast")) {
            fast = true;
        } else if (!strcmp(argv[i], "--dynamic")) {
            dynamic = true;
        } else if (!strcmp(argv[i], "--anno-filename")) {
            filename_anno = true;
        } else if (!strcmp(argv[i], "--anno-header")) {
            fasta_anno = true;
        } else if (!strcmp(argv[i], "--anno-label")) {
            anno_labels.emplace_back(argv[++i]);
        } else if (!strcmp(argv[i], "--coord-binsize")) {
            genome_binsize_anno = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--suppress-unlabeled")) {
            suppress_unlabeled = true;
        } else if (!strcmp(argv[i], "--row-annotator")) {
            use_row_annotator = true;
        } else if (!strcmp(argv[i], "--sparse")) {
            sparse = true;
        } else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--parallel")) {
            parallel = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--parts-total")) {
            parts_total = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--part-idx")) {
            part_idx = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bins-per-thread")) {
            num_bins_per_thread = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-k") || !strcmp(argv[i], "--kmer-length")) {
            k = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--filter-abund")) {
            max_unreliable_abundance = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--filter-thres")) {
            unreliable_kmers_threshold = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--filter-k")) {
            filter_k = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--mem-cap-gb")) {
            memory_available = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--dump-raw-anno")) {
            dump_raw_anno = true;
        } else if (!strcmp(argv[i], "--generate-fasta")) {
            generate_filtered_fasta = true;
        } else if (!strcmp(argv[i], "--generate-fastq")) {
            generate_filtered_fastq = true;
        } else if (!strcmp(argv[i], "--kmc")) {
            use_kmc = true;
            //TODO: add into some USAGE description
        } else if (!strcmp(argv[i], "--bloom-false-pos-prob")) {
            bloom_fpp = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--bloom-bits-per-edge")) {
            bloom_bits_per_edge = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--discovery-fraction")) {
            discovery_fraction = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--query-presence")) {
            query_presence = true;
        } else if (!strcmp(argv[i], "--filter-present")) {
            filter_present = true;
        } else if (!strcmp(argv[i], "--count-labels")) {
            count_labels = true;
        } else if (!strcmp(argv[i], "--bloom-hash-functions")) {
            bloom_num_hash_functions = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--bloom-test-num-kmers")) {
            bloom_test_num_kmers = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--align-length")) {
            alignment_length = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "--frequency")) {
            frequency = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--distance")) {
            distance = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--outfile-base")) {
            outfbase = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--reference")) {
            refpath = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--header-delimiter")) {
            fasta_header_delimiter = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--labels-delimiter")) {
            anno_labels_delimiter = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "--num-splits")) {
            nsplits = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--kmer-mapping-mode")) {
            kmer_mapping_mode = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--num-top-labels")) {
            num_top_labels = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--suffix")) {
            suffix = argv[++i];
        } else if (!strcmp(argv[i], "--state")) {
            state = string_to_state(argv[++i]);
        } else if (!strcmp(argv[i], "--rename-cols")) {
            rename_instructions_file = std::string(argv[++i]);
        //} else if (!strcmp(argv[i], "--db-path")) {
        //    dbpath = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-a") || !strcmp(argv[i], "--annotator")) {
            infbase_annotators.emplace_back(argv[++i]);
        } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--infile-base")) {
            infbase = std::string(argv[++i]);
            infbase = utils::remove_suffix(infbase, ".dbg");
        } else if (!strcmp(argv[i], "--to-adj-list")) {
            to_adj_list = true;
        } else if (!strcmp(argv[i], "--to-fasta")) {
            to_fasta = true;
        } else if (!strcmp(argv[i], "--contigs")) {
            to_fasta = true;
            contigs = true;
        } else if (!strcmp(argv[i], "--to-row-format")) {
            to_row_annotator = true;
        } else if (!strcmp(argv[i], "--clear-dummy")) {
            clear_dummy = true;
        } else if (!strcmp(argv[i], "--internal")) {
            internal = true;
        } else if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--len-suffix")) {
            suffix_len = atoi(argv[++i]);
        //} else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--threads")) {
        //    num_threads = atoi(argv[++i]);
        //} else if (!strcmp(argv[i], "--debug")) {
        //    debug = true;
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(argv[0], identity);
            exit(0);
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "\nERROR: Unknown option %s\n\n", argv[i]);
            print_usage(argv[0], identity);
            exit(-1);
        } else {
            fname.push_back(argv[i]);
        }
    }

    if (!fname.size() && identity != STATS
                      && (identity != CONCATENATE || infbase.empty())) {
        std::string line;
        while (std::getline(std::cin, line)) {
            if (line.size())
                fname.push_back(line);
        }
    }

    bool print_usage_and_exit = false;

    if (outfbase.size() && !utils::check_if_writable(outfbase)) {
        std::cerr << "Error: Can't write to " << outfbase << std::endl
                  << "Check if the path is correct" << std::endl;
        exit(1);
    }

    if (nsplits == 0) {
        std::cerr << "Error: Invalid number of splits" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity != CONCATENATE && identity != STATS && !fname.size())
        print_usage_and_exit = true;

    if (identity == CONCATENATE && !(fname.empty() ^ infbase.empty())) {
        std::cerr << "Error: Either set all chunk filenames"
                  << " or use the -i and -l options" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == FILTER)
        filter_k = k;

    if (identity == ALIGN && infbase.empty())
        print_usage_and_exit = true;

    if (identity == CLASSIFY && infbase.empty())
        print_usage_and_exit = true;

    if (identity == ANNOTATE && infbase.empty())
        print_usage_and_exit = true;

    if (identity == ANNOTATE_COORDINATES && infbase.empty())
        print_usage_and_exit = true;

    if ((identity == ANNOTATE || identity == ANNOTATE_BLOOM)
            && !filename_anno && !fasta_anno && !anno_labels.size()) {
        std::cerr << "Error: No annotation to add" << std::endl;
        print_usage_and_exit = true;
    }

    if ((identity == ANNOTATE || identity == ANNOTATE_BLOOM
                              || identity == ANNOTATE_COORDINATES)
            && outfbase.empty())
        outfbase = infbase;

    if (identity == EXTEND && (outfbase.empty() || infbase.empty()))
        print_usage_and_exit = true;

    if (identity == MERGE_ANNOTATIONS && outfbase.empty())
        print_usage_and_exit = true;

    if (identity == CLASSIFY && infbase_annotators.empty())
        infbase_annotators.push_back(infbase);

    if (identity == ANNOTATE_BLOOM && infbase.empty())
        print_usage_and_exit = true;

    if (identity == ANNOTATE_BLOOM && bloom_fpp < 0
                                   && bloom_bits_per_edge < 0) {
        std::cerr << "ERROR: please specify either a false positive probability"
                  << " or the number of bits per edge." << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == TRANSFORM && fname.size() != 1)
        print_usage_and_exit = true;

    if (identity == TRANSFORM_ANNOTATION && outfbase.empty())
        print_usage_and_exit = true;

    if (identity == MERGE && fname.size() < 2)
        print_usage_and_exit = true;

    if (identity == COMPARE && fname.size() != 2)
        print_usage_and_exit = true;

    if (discovery_fraction < 0 || discovery_fraction > 1)
        print_usage_and_exit = true;

    // if misused, provide help screen for chosen identity and exit
    if (print_usage_and_exit) {
        print_usage(argv[0], identity);
        exit(-1);
    }
}

std::string Config::state_to_string(StateType state) {
    switch (state) {
        case STAT:
            return "fast";
        case DYN:
            return "dynamic";
        case SMALL:
            return "small";
        default:
            assert(false);
            return "Never happens";
    }
}

Config::StateType Config::string_to_state(const std::string &string) {
    if (string == "fast") {
        return StateType::STAT;
    } else if (string == "dynamic") {
        return StateType::DYN;
    } else if (string == "small") {
        return StateType::SMALL;
    } else {
        throw std::runtime_error("Error: unknown graph state");
    }
}

void Config::print_usage(const std::string &prog_name, IdentityType identity) {
    fprintf(stderr, "Comprehensive metagenome graph representation -- Version 0.1\n\n");
    //fprintf(stderr, "This program is the first implementation of a\n");
    //fprintf(stderr, "meta-metagenome graph for identification and annotation\n");
    //fprintf(stderr, "purposes.\n\n");

    switch (identity) {
        case NO_IDENTITY: {
            fprintf(stderr, "Usage: %s <command> [command specific options]\n\n", prog_name.c_str());

            fprintf(stderr, "Available commands:\n");

            fprintf(stderr, "\tfilter\t\tfilter out reads with rare k-mers and dump\n");
            fprintf(stderr, "\t\t\tfilters to disk\n\n");

            fprintf(stderr, "\tbuild\t\tconstruct a graph object from input sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats or integrate sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats into a given graph\n\n");

            fprintf(stderr, "\textend\t\textend an existing graph with new sequences from\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats or integrate sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats\n\n");

            fprintf(stderr, "\tmerge\t\tintegrate a given set of graph structures\n");
            fprintf(stderr, "\t\t\tand output a new graph structure\n\n");

            fprintf(stderr, "\tconcatenate\tcombine the results of the external merge or\n");
            fprintf(stderr, "\t\t\tconstruction and output the resulting graph structure\n\n");

            fprintf(stderr, "\tcompare\t\tcheck whether two given graphs are identical\n\n");

            fprintf(stderr, "\talign\t\talign sequences provided in fast[a|q] files\n");
            fprintf(stderr, "\t\t\tto graph\n\n");

            fprintf(stderr, "\tstats\t\tprint graph statistics for given graph(s)\n\n");

            fprintf(stderr, "\tfilter_stats\tget statistics for filters\n\n");

            fprintf(stderr, "\tannotate\tgiven a graph and a fast[a|q] file, annotate\n");
            fprintf(stderr, "\t\t\tthe respective kmers\n\n");

            fprintf(stderr, "\tcoordinate\tgiven a graph and a fast[a|q] file, annotate\n");
            fprintf(stderr, "\t\t\tkmers with their respective coordinates in genomes\n\n");

            fprintf(stderr, "\tbloom\t\tgiven a graph and a fast[a|q] file, annotate\n");
            fprintf(stderr, "\t\t\tthe respective kmers using Bloom filters\n\n");

            fprintf(stderr, "\tmerge_anno\tmerge annotation columns\n\n");

            fprintf(stderr, "\tclassify\tannotate sequences from fast[a|q] files\n\n");

            fprintf(stderr, "\ttransform\tgiven a graph, transform it to other formats\n\n");

            fprintf(stderr, "\ttransform_anno\tchange representation of the graph annotation\n\n");

            return;
        }
        case EXPERIMENT: {
            fprintf(stderr, "Usage: %s experiment ???\n\n", prog_name.c_str());
        } break;
        case BUILD: {
            fprintf(stderr, "Usage: %s build [options] FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for build:\n");
            fprintf(stderr, "\t   --reference [STR] \tbasename of reference sequence []\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\tbasename of output file []\n");
            fprintf(stderr, "\t   --mem-cap-gb [INT] \tmaximum memory available, in Gb [inf]\n");
            fprintf(stderr, "\t-k --kmer-length [INT] \tlength of the k-mer to use [3]\n");
            fprintf(stderr, "\t-r --reverse \t\tadd reverse complement reads [off]\n");
            fprintf(stderr, "\t   --filter-abund [INT] threshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] max allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --dynamic \t\tuse dynamic build method [off]\n");
            fprintf(stderr, "\t   --print \t\tprint graph table to the screen [off]\n");
            fprintf(stderr, "\t   --suffix \t\tbuild graph chunk only for k-mers with the suffix given [off]\n");
            fprintf(stderr, "\t-s --num-splits \tdefine the minimum number of bins to split kmers into [1]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case EXTEND: {
            fprintf(stderr, "Usage: %s extend -i <graph_basename> -o <extended_graph_basename> [options] <FASTQ1> [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for extend:\n");
            fprintf(stderr, "\t   --reference [STR] \tbasename of reference sequence []\n");
            fprintf(stderr, "\t-a --annotator [STR] \tbasename of annotator to extend []\n");
            fprintf(stderr, "\t   --row-annotator \tuse row based annotator instead of column compressor [off]\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\tbasename of output file []\n");
            fprintf(stderr, "\t   --state [STR] \tstate of the extended graph (either 'fast' or 'dynamic') [fast]\n");
            fprintf(stderr, "\t-r --reverse \t\tadd reverse complement reads [off]\n");
            fprintf(stderr, "\t   --filter-abund [INT] threshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] max allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --print \t\tprint graph table to the screen [off]\n");
            // fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case FILTER: {
            fprintf(stderr, "Usage: %s filter [options] --filter-abund <cutoff> FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for filter:\n");
            fprintf(stderr, "\t-k --kmer-length [INT] \tlength of the k-mer to use [3]\n");
            fprintf(stderr, "\t   --kmc \t\tuse precomputed KMC database with k-mer counts\n");
            fprintf(stderr, "\t-r --reverse \t\tadd reverse complement reads [off]\n");
            fprintf(stderr, "\t   --filter-abund [INT] threshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] max allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --generate-fasta \twrite filtered reads to disk in fasta format [off]\n");
            fprintf(stderr, "\t   --generate-fastq \twrite filtered reads to disk in fastq format [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case FILTER_STATS: {
            fprintf(stderr, "Usage: %s filter_stats [options] --filter-abund <cutoff> FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for filter:\n");
            fprintf(stderr, "\t   --filter-abund [INT] \tthreshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] \tmax allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \t\tlength of k-mers used for counting and filtering [3]\n");
        } break;
        case ALIGN: {
            fprintf(stderr, "Usage: %s align -i <graph_basename> [options] <FASTQ1> [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for align:\n");
            fprintf(stderr, "\t-r --reverse \t\t\talso annotate reverse complement reads [off]\n");
            fprintf(stderr, "\t   --filter-abund [INT] \tthreshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] \tmax allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \t\tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --query-presence\t\ttest sequences for presence [off]\n");
            fprintf(stderr, "\t   --kmer-mapping-mode\t\tlevel of heuristics to use for unmapped k-mers (0, 1, or 2) [0]\n");
            fprintf(stderr, "\t   --discovery-fraction [FLOAT]\tfraction of k-mers required to count sequence [1.0]\n");
            fprintf(stderr, "\t   --filter-present\t\treport only present input sequences [off]\n");
            fprintf(stderr, "\t   --count-kmers \t\tquery the number of k-mers discovered [off]\n");
            fprintf(stderr, "\t   --align-length [INT]\t\tlength of subsequences to align [k]\n");
            fprintf(stderr, "\t-d --distance [INT] \t\tmax allowed alignment distance [0]\n");
        } break;
        case COMPARE: {
            fprintf(stderr, "Usage: %s compare [options] GRAPH1 GRAPH2\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for compare:\n");
            fprintf(stderr, "\t   --internal \t\tcompare internal graph representations\n");
        } break;
        case MERGE: {
            fprintf(stderr, "Usage: %s merge [options] GRAPH1 GRAPH2 [[GRAPH3] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for merge:\n");
            fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file []\n");
            fprintf(stderr, "\t   --print \t\t\tprint graph table to the screen [off]\n");
            fprintf(stderr, "\t-b --bins-per-thread [INT] \tnumber of bins each thread computes on average [1]\n");
            fprintf(stderr, "\t   --dynamic \t\t\tdynamic merge by adding traversed paths [off]\n");
            fprintf(stderr, "\t   --part-idx [INT] \t\tidx to use when doing external merge []\n");
            fprintf(stderr, "\t   --parts-total [INT] \t\ttotal number of parts in external merge[]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case CONCATENATE: {
            fprintf(stderr, "Usage: %s concatenate [options] [[CHUNK] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for merge:\n");
            fprintf(stderr, "\t-i --infile-base [STR] \t\tload graph chunks from files '<infile-base>.<suffix>.dbgchunk' []\n");
            fprintf(stderr, "\t-l --len-suffix [INT] \t\titerate all possible suffices of the length given [0]\n");
            fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file []\n");
            fprintf(stderr, "\t   --print \t\t\tprint graph table to the screen [off]\n");
            // fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case STATS: {
            fprintf(stderr, "Usage: %s stats [options] GRAPH1 [[GRAPH2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for stats:\n");
            // fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t-a --annotator [STR] \tbasename of annotator to update []\n");
            fprintf(stderr, "\t   --row-annotator \tuse row based annotator instead of column compressor [off]\n");
            fprintf(stderr, "\t   --print \t\tprint graph table to the screen [off]\n");
            fprintf(stderr, "\t   --print-internal \tprint internal graph representation to the screen [off]\n");
        } break;
        case ANNOTATE: {
            fprintf(stderr, "Usage: %s annotate -i <graph_basename> [options] <PATH1> [[PATH2] ...]\n"
                            "\tEach path is given as file in fasta or fastq format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t   --reference [STR] \t\tbasename of reference sequence []\n");
            fprintf(stderr, "\t-a --annotator [STR] \t\tbasename of annotator to update []\n");
            fprintf(stderr, "\t   --row-annotator \t\tuse row based annotator instead of column compressor [off]\n");
            fprintf(stderr, "\t   --sparse \t\t\tuse the row-major sparse matrix to annotate graph [off]\n");
            fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file [<graph_basename>]\n");
            fprintf(stderr, "\t-r --reverse \t\t\talso annotate reverse complement reads [off]\n");
            fprintf(stderr, "\t   --filter-abund [INT] \tthreshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] \tmax allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \t\tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --anno-filename \t\tinclude filenames as annotation labels [off]\n");
            fprintf(stderr, "\t   --anno-header \t\textract annotation labels from headers of sequences in files [off]\n");
            fprintf(stderr, "\t   --header-delimiter [STR]\tdelimiter for splitting annotation header into multiple labels [off]\n");
            fprintf(stderr, "\t   --anno-label [STR]\t\tadd label to annotation for all sequences from the files passed []\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case ANNOTATE_COORDINATES: {
            fprintf(stderr, "Usage: %s coordinate -i <graph_basename> [options] <PATH1> [[PATH2] ...]\n"
                            "\tEach path is given as file in fasta or fastq format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t-a --annotator [STR] \t\tbasename of annotator to update []\n");
            fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file [<graph_basename>]\n");
            fprintf(stderr, "\t-r --reverse \t\t\talso annotate reverse complement reads [off]\n");
            fprintf(stderr, "\t   --filter-abund [INT] \tthreshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] \tmax allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \t\tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --coord-binsize [INT]\tstepsize for k-mer coordinates in input sequences from the fasta files [1000]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case MERGE_ANNOTATIONS: {
            fprintf(stderr, "Usage: %s merge_anno [options] -o <annotator_basename> <ANNOT1> [[ANNOT2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t   --row-annotator \tuse row based annotator instead of column compressor [off]\n");
            fprintf(stderr, "\t   --sparse \t\tuse the row-major sparse matrix to annotate graph [off]\n");
            // fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case ANNOTATE_BLOOM: {
            fprintf(stderr, "Usage: %s bloom -i <graph_basename> [options] <PATH1> [[PATH2] ...]\n"
                            "\tEach path is given as file in fasta or fastq format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for bloom:\n");
            fprintf(stderr, "\t   --reference [STR] \t\t\tbasename of reference sequence []\n");
            fprintf(stderr, "\t-r --reverse \t\t\t\talso annotate reverse complement reads [off]\n");
            fprintf(stderr, "\t   --filter-abund [INT] \t\tthreshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] \t\tmax allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \t\t\tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --anno-filename \t\t\tinclude filenames as annotation labels [off]\n");
            fprintf(stderr, "\t   --anno-header \t\t\textract annotation labels from headers of sequences in files [off]\n");
            fprintf(stderr, "\t   --header-delimiter [STR]\t\tdelimiter for splitting annotation header into multiple labels [off]\n");
            fprintf(stderr, "\t   --anno-label [STR]\t\t\tadd label to annotation for all sequences from the files passed []\n");
            // fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
            // fprintf(stderr, "\t-b --bins-per-thread [INT] \tnumber of bins each thread computes on average [1]\n");
            // fprintf(stderr, "\t-f --frequency [INT] \t\twhen a, annotate only every a-th kmer [1]\n");
            // fprintf(stderr, "\t   --db-path \tpath that is used to store the annotations database []\n");
            fprintf(stderr, "\t   --bloom-false-pos-prob [FLOAT]\tFalse positive probability in bloom filter [-1]\n");
            fprintf(stderr, "\t   --bloom-bits-per-edge [FLOAT] \tBits per edge used in bloom filter annotator [0.4]\n");
            fprintf(stderr, "\t   --bloom-hash-functions [INT] \tNumber of hash functions used in bloom filter [off]\n");
            fprintf(stderr, "\t   --bloom-test-num-kmers \t\tEstimate false positive rate for every n k-mers [0]\n");
        } break;
        case CLASSIFY: {
            fprintf(stderr, "Usage: %s classify -i <graph_basename> [options] <FILE1> [[FILE2] ...]\n"
                            "\tEach file is given in fasta or fastq format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for classify:\n");
            fprintf(stderr, "\t-r --reverse \t\t\tclassify reverse complement sequences [off]\n");
            fprintf(stderr, "\t   --filter-abund [INT] \tthreshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] \tmax allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \t\tlength of k-mers used for counting and filtering [3]\n");
            // fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file []\n");
            fprintf(stderr, "\t-a --annotator [STR] \t\tbasename of annotator [<graph_basename>]\n");
            fprintf(stderr, "\t   --fast \t\t\tuse fast column based annotator (with auxiliary index) [off]\n");
            fprintf(stderr, "\t   --row-annotator \t\tuse row based annotator instead of column compressor [off]\n");
            fprintf(stderr, "\t   --sparse \t\t\tuse the row-major sparse matrix to annotate graph [off]\n");
            fprintf(stderr, "\t   --suppress-unlabeled \tdo not show results for sequences missing in graph [off]\n");
            fprintf(stderr, "\t   --count-labels \t\tcount labels for k-mers from querying sequences [off]\n");
            fprintf(stderr, "\t   --num-top-labels \t\tmaximum number of frequent labels to print [off]\n");
            fprintf(stderr, "\t   --discovery-fraction \tfraction of labeled k-mers required for annotation [1.0]\n");
            fprintf(stderr, "\t   --labels-delimiter [STR]\tdelimiter for annotation labels [\":\"]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
            // fprintf(stderr, "\t-d --distance [INT] \tMax allowed alignment distance [0]\n");
        } break;
        case TRANSFORM: {
            fprintf(stderr, "Usage: %s transform [options] GRAPH\n\n", prog_name.c_str());

            fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file []\n");
            fprintf(stderr, "\t   --clear-dummy \t\terase all redundant dummy edges [off]\n");
            fprintf(stderr, "\t   --state [STR] \t\tchange the graph state (either 'fast' or 'dynamic') [fast]\n");
            fprintf(stderr, "\t   --to-adj-list \t\twrite the adjacency list to file [off]\n");
            // fprintf(stderr, "\t   --to-sequences \t\twrite contig sequences to file [off]\n");
            fprintf(stderr, "\t   --to-fasta \t\t\textract sequences from graph and write to compressed fasta file [off]\n");
            fprintf(stderr, "\t   --contigs \t\t\textract all simple paths from graph and write to compressed fasta file [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case TRANSFORM_ANNOTATION: {
            fprintf(stderr, "Usage: %s transform_anno [options] -o <annotator_basename> ANNOTATOR\n\n", prog_name.c_str());

            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --rename-cols [STR]\tfile with rules for renaming annotation labels []\n");
            fprintf(stderr, "\t                      \texample: 'L_1 L_1_renamed\n");
            fprintf(stderr, "\t                      \t          L_2 L_2_renamed\n");
            fprintf(stderr, "\t                      \t          L_2 L_2_renamed\n");
            fprintf(stderr, "\t                      \t          ... ...........'\n");
            fprintf(stderr, "\t   --to-row-format \ttransform annotations to the row-wise format [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
    }

    fprintf(stderr, "\n\tGeneral options:\n");
    fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
    fprintf(stderr, "\t-q --quiet \t\tproduce as little log output as posible [off]\n");
    fprintf(stderr, "\t-h --help \t\tprint usage info\n");
    fprintf(stderr, "\n");
}
