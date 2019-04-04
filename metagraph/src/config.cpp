#include "config.hpp"

#include <cstring>
#include <iostream>
#include <unordered_set>

#include "utils.hpp"
#include "threading.hpp"


Config::Config(int argc, const char *argv[]) {
    // provide help overview if no identity was given
    if (argc == 1) {
        print_usage(argv[0]);
        exit(-1);
    }

    // parse identity from first command line argument
    if (!strcmp(argv[1], "build")) {
        identity = BUILD;
    } else if (!strcmp(argv[1], "merge")) {
        identity = MERGE;
    } else if (!strcmp(argv[1], "extend")) {
        identity = EXTEND;
    } else if (!strcmp(argv[1], "concatenate")) {
        identity = CONCATENATE;
    } else if (!strcmp(argv[1], "compare")) {
        identity = COMPARE;
    } else if (!strcmp(argv[1], "align")) {
        identity = ALIGN;
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
    } else if (!strcmp(argv[1], "coordinate")) {
        identity = ANNOTATE_COORDINATES;
    } else if (!strcmp(argv[1], "merge_anno")) {
        identity = MERGE_ANNOTATIONS;
    } else if (!strcmp(argv[1], "query")) {
        identity = QUERY;
    } else if (!strcmp(argv[1], "server_query")) {
        identity = SERVER_QUERY;
    } else if (!strcmp(argv[1], "transform")) {
        identity = TRANSFORM;
    } else if (!strcmp(argv[1], "transform_anno")) {
        identity = TRANSFORM_ANNOTATION;
    } else if (!strcmp(argv[1], "assemble")) {
        identity = ASSEMBLE;
    } else if (!strcmp(argv[1], "relax_brwt")) {
        identity = RELAX_BRWT;
    } else if (!strcmp(argv[1], "call_variants")) {
        identity = CALL_VARIANTS;
    } else if (!strcmp(argv[1], "parse_taxonomy")) {
        identity = PARSE_TAXONOMY;
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
        } else if (!strcmp(argv[i], "--print")) {
            print_graph = true;
        } else if (!strcmp(argv[i], "--print-col-names")) {
            print_column_names = true;
        } else if (!strcmp(argv[i], "--print-internal")) {
            print_graph_internal_repr = true;
        } else if (!strcmp(argv[i], "--count-kmers")) {
            count_kmers = true;
        } else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--reverse")) {
            reverse = true;
        } else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--canonical")) {
            canonical = true;
        } else if (!strcmp(argv[i], "--complete")) {
            complete = true;
        } else if (!strcmp(argv[i], "--dynamic")) {
            dynamic = true;
        } else if (!strcmp(argv[i], "--no-shrink")) {
            mark_dummy_kmers = false;
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
        } else if (!strcmp(argv[i], "--sparse")) {
            sparse = true;
        } else if (!strcmp(argv[i], "--fast")) {
            fast = true;
        } else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--parallel")) {
            parallel = atoi(argv[++i]);
            set_num_threads(parallel);
        } else if (!strcmp(argv[i], "--parts-total")) {
            parts_total = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--part-idx")) {
            part_idx = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bins-per-thread")) {
            num_bins_per_thread = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-k") || !strcmp(argv[i], "--kmer-length")) {
            k = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--min-count")) {
            min_count = std::max(atoi(argv[++i]), 1);
        } else if (!strcmp(argv[i], "--max-count")) {
            max_count = atoi(argv[++i]);
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
        } else if (!strcmp(argv[i], "--discovery-fraction")) {
            discovery_fraction = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--query-presence")) {
            query_presence = true;
        } else if (!strcmp(argv[i], "--filter-present")) {
            filter_present = true;
        } else if (!strcmp(argv[i], "--count-labels")) {
            count_labels = true;
        } else if (!strcmp(argv[i], "--align-length")) {
            alignment_length = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--align-num-paths")) {
            alignment_num_top_paths = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--align-sw-threshold")) {
            alignment_sw_threshold = std::stof(argv[++i]);
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
        } else if (!strcmp(argv[i], "--separately")) {
            separately = true;
        } else if (!strcmp(argv[i], "--kmer-mapping-mode")) {
            kmer_mapping_mode = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--num-top-labels")) {
            num_top_labels = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--port")) {
            port = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--suffix")) {
            suffix = argv[++i];
        } else if (!strcmp(argv[i], "--state")) {
            state = string_to_state(argv[++i]);

        } else if (!strcmp(argv[i], "--anno-type")) {
            anno_type = string_to_annotype(argv[++i]);
        } else if (!strcmp(argv[i], "--graph")) {
            graph_type = string_to_graphtype(argv[++i]);
        } else if (!strcmp(argv[i], "--rename-cols")) {
            rename_instructions_file = std::string(argv[++i]);
        //} else if (!strcmp(argv[i], "--db-path")) {
        //    dbpath = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-a") || !strcmp(argv[i], "--annotator")) {
            infbase_annotators.emplace_back(argv[++i]);
        } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--infile-base")) {
            infbase = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--to-adj-list")) {
            to_adj_list = true;
        } else if (!strcmp(argv[i], "--unitigs")) {
            unitigs = true;
        } else if (!strcmp(argv[i], "--header")) {
            header = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--prune-end")) {
            pruned_dead_end_size = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--count-dummy")) {
            count_dummy = true;
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
        } else if (!strcmp(argv[i], "--greedy")) {
            greedy_brwt = true;
        } else if (!strcmp(argv[i], "--arity")) {
            arity_brwt = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--relax-arity")) {
            relax_arity_brwt = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(argv[0], identity);
            exit(0);
        } else if (!strcmp(argv[i], "--label-mask-in")) {
            label_mask_in.emplace_back(argv[++i]);
        } else if (!strcmp(argv[i], "--label-mask-out")) {
            label_mask_out.emplace_back(argv[++i]);
        } else if (!strcmp(argv[i], "--label-mask-out-fraction")) {
            label_mask_out_fraction = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--label-filter")) {
            label_filter.emplace_back(argv[++i]);
        } else if (!strcmp(argv[i], "--call-bubbles")) {
            call_bubbles = true;
        } else if (!strcmp(argv[i], "--accession")) {
            accession2taxid = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--taxonomy")) {
            taxonomy_nodes = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--taxonomy-map")) {
            taxonomy_map = std::string(argv[++i]);
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "\nERROR: Unknown option %s\n\n", argv[i]);
            print_usage(argv[0], identity);
            exit(-1);
        } else {
            fname.push_back(argv[i]);
        }
    }

    // given kmc_pre and kmc_suf pair, only include one
    // this still allows for the same file to be included multiple times
    std::unordered_set<std::string> kmc_file_set;

    for (auto it = fname.begin(); it != fname.end(); ++it) {
        if (utils::get_filetype(*it) == "KMC"
                && !kmc_file_set.insert(utils::remove_suffix(*it, ".kmc_pre", ".kmc_suf")).second)
            fname.erase(it--);
    }

    utils::set_verbose(verbose);

    if (!fname.size() && identity != STATS
                      && identity != SERVER_QUERY
                      && !(identity == BUILD && complete)
                      && !(identity == CALL_VARIANTS)
                      && !(identity == PARSE_TAXONOMY)
                      && !(identity == CONCATENATE && !infbase.empty())) {
        std::string line;
        while (std::getline(std::cin, line)) {
            if (line.size())
                fname.push_back(line);
        }
    }

    bool print_usage_and_exit = false;

    if (identity != CONCATENATE
            && identity != STATS
            && identity != SERVER_QUERY
            && !(identity == BUILD && complete)
            && !(identity == CALL_VARIANTS)
            && !(identity == PARSE_TAXONOMY)
            && !fname.size())
        print_usage_and_exit = true;

    if (identity == CONCATENATE && !(fname.empty() ^ infbase.empty())) {
        std::cerr << "Error: Either set all chunk filenames"
                  << " or use the -i and -l options" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == CONCATENATE && outfbase.empty())
        print_usage_and_exit = true;

    if (identity == FILTER)
        filter_k = k;

    if (identity == ALIGN && infbase.empty())
        print_usage_and_exit = true;

    if ((identity == QUERY || identity == SERVER_QUERY) && infbase.empty())
        print_usage_and_exit = true;

    if (identity == ANNOTATE && infbase.empty())
        print_usage_and_exit = true;

    if (identity == ANNOTATE_COORDINATES && infbase.empty())
        print_usage_and_exit = true;

    if ((identity == ANNOTATE_COORDINATES
            || identity == ANNOTATE
            || identity == EXTEND) && infbase_annotators.size() > 1) {
        std::cerr << "Error: one annotator at most is allowed for extension." << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == ANNOTATE
            && !filename_anno && !fasta_anno && !anno_labels.size()) {
        std::cerr << "Error: No annotation to add" << std::endl;
        print_usage_and_exit = true;
    }

    if ((identity == ANNOTATE || identity == ANNOTATE_COORDINATES)
            && outfbase.empty())
        outfbase = utils::remove_suffix(infbase, ".orhashdbg", ".bitmapdbg", ".dbg");

    if (identity == EXTEND && (outfbase.empty() || infbase.empty()))
        print_usage_and_exit = true;

    if (identity == MERGE_ANNOTATIONS && outfbase.empty())
        print_usage_and_exit = true;

    if ((identity == QUERY || identity == SERVER_QUERY) && infbase_annotators.size() != 1)
        print_usage_and_exit = true;

    if ((identity == TRANSFORM
            || identity == TRANSFORM_ANNOTATION
            || identity == ASSEMBLE
            || identity == RELAX_BRWT)
                    && fname.size() != 1)
        print_usage_and_exit = true;

    if (identity == PARSE_TAXONOMY &&
            ((accession2taxid == "" && taxonomy_nodes == "") || outfbase == ""))
        print_usage_and_exit = true;

    if (identity == TRANSFORM_ANNOTATION && outfbase.empty())
        print_usage_and_exit = true;

    if (identity == MERGE && fname.size() < 2)
        print_usage_and_exit = true;

    if (identity == MERGE && outfbase.empty())
        print_usage_and_exit = true;

    if (identity == COMPARE && fname.size() != 2)
        print_usage_and_exit = true;

    if (discovery_fraction < 0 || discovery_fraction > 1)
        print_usage_and_exit = true;

    if (min_count >= max_count) {
        std::cerr << "Error: max-count must be greater than min-count" << std::endl;
        print_usage(argv[0], identity);
    }
    if (alignment_sw_threshold < 0 || alignment_sw_threshold > 1)
        print_usage_and_exit = true;

    if (outfbase.size() && !utils::check_if_writable(outfbase)) {
        std::cerr << "Error: Can't write to " << outfbase << std::endl
                  << "Check if the path is correct" << std::endl;
        exit(1);
    }

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
        case FAST:
            return "faster";
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
    } else if (string == "faster") {
        return StateType::FAST;
    } else {
        throw std::runtime_error("Error: unknown graph state");
    }
}

std::string Config::annotype_to_string(AnnotationType state) {
    switch (state) {
        case ColumnCompressed:
            return "column";
        case RowCompressed:
            return "row";
        case BRWT:
            return "brwt";
        case BinRelWT_sdsl:
            return "bin_rel_wt_sdsl";
        case BinRelWT:
            return "bin_rel_wt";
        case RowFlat:
            return "flat";
        case RBFish:
            return "rbfish";
        default:
            assert(false);
            return "Never happens";
    }
}

Config::AnnotationType Config::string_to_annotype(const std::string &string) {
    if (string == "column") {
        return AnnotationType::ColumnCompressed;
    } else if (string == "row") {
        return AnnotationType::RowCompressed;
    } else if (string == "brwt") {
        return AnnotationType::BRWT;
    } else if (string == "bin_rel_wt_sdsl") {
        return AnnotationType::BinRelWT_sdsl;
    } else if (string == "bin_rel_wt") {
        return AnnotationType::BinRelWT;
    } else if (string == "flat") {
        return AnnotationType::RowFlat;
    } else if (string == "rbfish") {
        return AnnotationType::RBFish;
    } else {
        std::cerr << "Error: unknown annotation representation" << std::endl;
        exit(1);
    }
}

Config::GraphType Config::string_to_graphtype(const std::string &string) {
    if (string == "succinct") {
        return GraphType::SUCCINCT;
    } else if (string == "hash") {
        return GraphType::HASH;
    } else if (string == "hashstr") {
        return GraphType::HASH_STR;
    } else if (string == "bitmap") {
        return GraphType::BITMAP;
    } else {
        std::cerr << "Error: unknown graph representation" << std::endl;
        exit(1);
    }
}

void Config::print_usage(const std::string &prog_name, IdentityType identity) {
    fprintf(stderr, "Metagraph: comprehensive metagenome graph representation -- Version 0.1\n\n");

    const char annotation_list[] = "('column', 'row', 'bin_rel_wt_sdsl', 'bin_rel_wt', 'flat', 'rbfish', 'brwt')";

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

            fprintf(stderr, "\tmerge_anno\tmerge annotation columns\n\n");

            fprintf(stderr, "\ttransform\tgiven a graph, transform it to other formats\n\n");

            fprintf(stderr, "\ttransform_anno\tchange representation of the graph annotation\n\n");

            fprintf(stderr, "\tassemble\tgiven a graph, extract sequences from it\n\n");

            fprintf(stderr, "\trelax_brwt\toptimize the tree structure in brwt annotator\n\n");

            fprintf(stderr, "\tquery\t\tannotate sequences from fast[a|q] files\n\n");
            fprintf(stderr, "\tserver_query\tannotate received sequences and send annotations back\n\n");

            fprintf(stderr, "\tcall_variants\tgenerate a masked annotated graph and call variants\n");
            fprintf(stderr, "\t\t\trelative to unmasked graph\n\n");

            fprintf(stderr, "\tparse_taxonomy\tgenerate NCBI Accession ID to Taxonomy ID mapper\n\n");

            return;
        }
        case EXPERIMENT: {
            fprintf(stderr, "Usage: %s experiment ???\n\n", prog_name.c_str());
        } break;
        case FILTER: {
            fprintf(stderr, "Usage: %s filter [options] --min-count <cutoff> FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for filter:\n");
            fprintf(stderr, "\t   --max-count [INT] \tmax k-mer abundance, excluding [inf]\n");
            fprintf(stderr, "\t   --filter-thres [INT] max allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t-r --reverse \t\tprocess reverse complement sequences as well [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-k --kmer-length [INT] \tlength of the k-mer to use [3]\n");
            fprintf(stderr, "\t   --generate-fasta \twrite filtered reads to disk in FASTA format [off]\n");
            fprintf(stderr, "\t   --generate-fastq \twrite filtered reads to disk in FASTQ format [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case FILTER_STATS: {
            fprintf(stderr, "Usage: %s filter_stats [options] --min-count <cutoff> FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for filter:\n");
            fprintf(stderr, "\t   --filter-k [INT] \tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --filter-thres [INT] max allowed number of unreliable kmers in reliable reads [0]\n");
        } break;
        case BUILD: {
            fprintf(stderr, "Usage: %s build [options] FILE1 [[FILE2] ...]\n"
                            "\tEach input file is given in FASTA, FASTQ, or VCF format.\n"
                            "\tNote that VCF files must be in plain text or bgzip format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for build:\n");
            fprintf(stderr, "\t   --min-count [INT] \tmin k-mer abundance, including [1]\n");
            fprintf(stderr, "\t   --max-count [INT] \tmax k-mer abundance, excluding [inf]\n");
            fprintf(stderr, "\t   --filter-k [INT] \tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --filter-thres [INT] max allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --reference [STR] \tbasename of reference sequence (for parsing VCF files) []\n");
            fprintf(stderr, "\t-r --reverse \t\tprocess reverse complement sequences as well [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --graph [STR] \tgraph representation: succinct / bitmap / hash / hashstr [succinct]\n");
            fprintf(stderr, "\t   --count-kmers \tcount k-mers and build weighted graph [off]\n");
            fprintf(stderr, "\t-k --kmer-length [INT] \tlength of the k-mer to use [3]\n");
            fprintf(stderr, "\t-c --canonical \t\tindex only canonical k-mers (e.g. for read sets) [off]\n");
            fprintf(stderr, "\t   --complete \t\tconstruct a complete graph (only for Bitmap graph) [off]\n");
            fprintf(stderr, "\t   --mem-cap-gb [INT] \tpreallocated buffer size in Gb [0]\n");
            fprintf(stderr, "\t   --dynamic \t\tuse dynamic build method [off]\n");
            fprintf(stderr, "\t-l --len-suffix [INT] \tk-mer suffix length for building graph from chunks [0]\n");
            fprintf(stderr, "\t   --suffix \t\tbuild graph chunk only for k-mers with the suffix given [off]\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\tbasename of output file []\n");
            fprintf(stderr, "\t   --no-shrink \t\tdo not build mask for dummy k-mers (only for Succinct graph) [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case EXTEND: {
            fprintf(stderr, "Usage: %s extend -i <GRAPH> -o <extended_graph_basename> [options] FILE1 [[FILE2] ...]\n"
                            "\tEach input file is given in FASTA, FASTQ, or VCF format.\n"
                            "\tNote that VCF files must be in plain text or bgzip format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for extend:\n");
            fprintf(stderr, "\t   --min-count [INT] \tmin k-mer abundance, including [1]\n");
            fprintf(stderr, "\t   --max-count [INT] \tmax k-mer abundance, excluding [inf]\n");
            fprintf(stderr, "\t   --filter-k [INT] \tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --filter-thres [INT] max allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --reference [STR] \tbasename of reference sequence (for parsing VCF files) []\n");
            fprintf(stderr, "\t-r --reverse \t\tprocess reverse complement sequences as well [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-a --annotator [STR] \tannotator to extend []\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\tbasename of output file []\n");
            // fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case ALIGN: {
            fprintf(stderr, "Usage: %s align -i <GRAPH> [options] FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for align:\n");
            fprintf(stderr, "\t-r --reverse \t\t\talign reverse complement sequences as well [off]\n");
            fprintf(stderr, "\t   --query-presence \t\ttest sequences for presence [off]\n");
            fprintf(stderr, "\t   --kmer-mapping-mode \t\tlevel of heuristics to use for unmapped k-mers (0, 1, or 2) [0]\n");
            fprintf(stderr, "\t   --discovery-fraction [FLOAT] fraction of k-mers required to count sequence [1.0]\n");
            fprintf(stderr, "\t   --filter-present \t\treport only present input sequences [off]\n");
            fprintf(stderr, "\t   --count-kmers \t\tquery the number of k-mers discovered [off]\n");
            fprintf(stderr, "\t   --align-length [INT]\t\tlength of subsequences to align [k]\n");
            fprintf(stderr, "\t   --align-num-paths [INT]\t\tnumber of parallel paths to explore at any point [10]\n");
            fprintf(stderr, "\t   --align-sw-threshold [FLOAT]\t\tthreshold proportion to the path length to determine if Smith Waterman should be computed for a path. [0.1]\n");
            fprintf(stderr, "\t-d --distance [INT] \t\tmax allowed alignment distance [0]\n");
        } break;
        case COMPARE: {
            fprintf(stderr, "Usage: %s compare [options] GRAPH1 GRAPH2\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for compare:\n");
            fprintf(stderr, "\t   --internal \t\tcompare internal graph representations\n");
        } break;
        case MERGE: {
            fprintf(stderr, "Usage: %s merge -o <graph_basename> [options] GRAPH1 GRAPH2 [[GRAPH3] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for merge:\n");
            fprintf(stderr, "\t-b --bins-per-thread [INT] \tnumber of bins each thread computes on average [1]\n");
            fprintf(stderr, "\t   --dynamic \t\t\tdynamic merge by adding traversed paths [off]\n");
            fprintf(stderr, "\t   --part-idx [INT] \t\tidx to use when doing external merge []\n");
            fprintf(stderr, "\t   --parts-total [INT] \t\ttotal number of parts in external merge[]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case CONCATENATE: {
            fprintf(stderr, "Usage: %s concatenate -o <graph_basename> [options] [[CHUNK] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for merge:\n");
            fprintf(stderr, "\t   --graph [STR] \tgraph representation: succinct / bitmap [succinct]\n");
            fprintf(stderr, "\t-i --infile-base [STR] \tload graph chunks from files '<infile-base>.<suffix>.<type>.chunk' []\n");
            fprintf(stderr, "\t-l --len-suffix [INT] \titerate all possible suffices of the length given [0]\n");
            fprintf(stderr, "\t-c --canonical \t\tcanonical graph mode (e.g. for read sets) [off]\n");
            // fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case TRANSFORM: {
            fprintf(stderr, "Usage: %s transform [options] GRAPH\n\n", prog_name.c_str());

            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --clear-dummy \terase all redundant dummy edges [off]\n");
            fprintf(stderr, "\t   --prune-end [INT] \tprune all dead ends of this length and shorter [0]\n");
            fprintf(stderr, "\t   --state [STR] \tchange state of succinct graph: fast / faster / dynamic / small [fast]\n");
            fprintf(stderr, "\t   --to-adj-list \twrite adjacency list to file [off]\n");
            fprintf(stderr, "\t   --header [STR] \theader for sequences in FASTA output []\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case ASSEMBLE: {
            fprintf(stderr, "Usage: %s assemble [options] GRAPH\n"
                            "\tAssemble contigs from de Bruijn graph and dump to compressed FASTA file.\n\n", prog_name.c_str());

            fprintf(stderr, "\t-o --outfile-base [STR] \t\tbasename of output file []\n");
            fprintf(stderr, "\t   --prune-end [INT] \t\t\tprune all dead ends of this length and shorter [0]\n");
            fprintf(stderr, "\t   --unitigs \t\t\t\textract unitigs [off]\n");
            fprintf(stderr, "\t   --header [STR] \t\t\theader for sequences in FASTA output []\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\t\tuse multiple threads for computation [1]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-a --annotator [STR] \t\t\tannotator to load []\n");
            fprintf(stderr, "\t   --label-mask-in [STR] \t\tlabel to include in masked graph\n");
            fprintf(stderr, "\t   --label-mask-out [STR] \t\tlabel to exclude from masked graph\n");
            fprintf(stderr, "\t   --label-mask-out-fraction [FLOAT] \tmaximum fraction of mask-out labels among the set of\n");
            fprintf(stderr, "\t                                     \tall matching mask-in and mask-out labels [0.0]\n");
        } break;
        case STATS: {
            fprintf(stderr, "Usage: %s stats [options] GRAPH1 [[GRAPH2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for stats:\n");
            fprintf(stderr, "\t   --print \t\tprint graph table to the screen [off]\n");
            fprintf(stderr, "\t   --print-internal \tprint internal graph representation to screen [off]\n");
            fprintf(stderr, "\t   --count-dummy \tshow number of dummy source and sink edges [off]\n");
            fprintf(stderr, "\t-a --annotator [STR] \tannotation []\n");
            fprintf(stderr, "\t   --print-col-names \tprint names of the columns in annotation to screen [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case ANNOTATE: {
            fprintf(stderr, "Usage: %s annotate -i <GRAPH> [options] FILE1 [[FILE2] ...]\n"
                            "\tEach file is given in FASTA, FASTQ, or VCF format.\n"
                            "\tNote that VCF files must be in plain text or bgzip format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t   --min-count [INT] \tmin k-mer abundance, including [1]\n");
            fprintf(stderr, "\t   --max-count [INT] \tmax k-mer abundance, excluding [inf]\n");
            fprintf(stderr, "\t   --filter-k [INT] \tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --filter-thres [INT] max allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --reference [STR] \tbasename of reference sequence (for parsing VCF files) []\n");
            fprintf(stderr, "\t-r --reverse \t\tprocess reverse complement sequences as well [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --anno-type [STR] \ttarget annotation representation: column / row [column]\n");
            fprintf(stderr, "\t-a --annotator [STR] \tannotator to update []\n");
            fprintf(stderr, "\t   --sparse \t\tuse the row-major sparse matrix to annotate graph [off]\n");
            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file [<GRAPH>]\n");
            fprintf(stderr, "\t   --separately \tannotate each file independently and dump to the same directory [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --anno-filename \t\tinclude filenames as annotation labels [off]\n");
            fprintf(stderr, "\t   --anno-header \t\textract annotation labels from headers of sequences in files [off]\n");
            fprintf(stderr, "\t   --header-delimiter [STR]\tdelimiter for splitting annotation header into multiple labels [off]\n");
            fprintf(stderr, "\t   --anno-label [STR]\t\tadd label to annotation for all sequences from the files passed []\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
            // fprintf(stderr, "\t   --fast \t\t\tannotate in fast regime (leads to repeated labels and bigger annotation) [off]\n");
        } break;
        case ANNOTATE_COORDINATES: {
            fprintf(stderr, "Usage: %s coordinate -i <GRAPH> [options] FASTA1 [[FASTA2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t-r --reverse \t\t\tprocess reverse complement sequences as well [off]\n");
            fprintf(stderr, "\t-a --annotator [STR] \t\tannotator to update []\n");
            fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file [<GRAPH>]\n");
            fprintf(stderr, "\t   --coord-binsize [INT]\tstepsize for k-mer coordinates in input sequences from the fasta files [1000]\n");
            fprintf(stderr, "\t   --fast \t\t\tannotate in fast regime [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case MERGE_ANNOTATIONS: {
            fprintf(stderr, "Usage: %s merge_anno [options] -o <annotator_basename> ANNOT1 [[ANNOT2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t   --anno-type [STR] \ttarget annotation representation [column]\n");
            fprintf(stderr, "\t\t"); fprintf(stderr, annotation_list); fprintf(stderr, "\n");
            // fprintf(stderr, "\t   --sparse \t\tuse the row-major sparse matrix to annotate graph [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case TRANSFORM_ANNOTATION: {
            fprintf(stderr, "Usage: %s transform_anno [options] -o <annotator_basename> ANNOTATOR\n\n", prog_name.c_str());

            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --rename-cols [STR]\tfile with rules for renaming annotation labels []\n");
            fprintf(stderr, "\t                      \texample: 'L_1 L_1_renamed\n");
            fprintf(stderr, "\t                      \t          L_2 L_2_renamed\n");
            fprintf(stderr, "\t                      \t          L_2 L_2_renamed\n");
            fprintf(stderr, "\t                      \t          ... ...........'\n");
            fprintf(stderr, "\t   --anno-type [STR] \ttarget annotation format [column]\n");
            fprintf(stderr, "\t\t"); fprintf(stderr, annotation_list); fprintf(stderr, "\n");
            fprintf(stderr, "\t   --arity  \t\tarity in the brwt tree [2]\n");
            fprintf(stderr, "\t   --greedy  \t\tuse greedy column partitioning in brwt construction [off]\n");
            fprintf(stderr, "\t   --fast  \t\ttransform annotation in memory without streaming [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case RELAX_BRWT: {
            fprintf(stderr, "Usage: %s relax_brwt [options] -o <annotator_basename> ANNOTATOR\n\n", prog_name.c_str());

            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --relax-arity [INT] \trelax brwt tree to optimize arity limited to this number [10]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case QUERY: {
            fprintf(stderr, "Usage: %s query -i <GRAPH> -a <ANNOTATION> [options] FILE1 [[FILE2] ...]\n"
                            "\tEach input file is given in FASTA or FASTQ format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for query:\n");
            fprintf(stderr, "\t-r --reverse \t\tquery reverse complement sequences as well [off]\n");
            fprintf(stderr, "\t   --sparse \t\tuse row-major sparse matrix for row annotation [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --count-labels \t\tcount labels for k-mers from querying sequences [off]\n");
            fprintf(stderr, "\t   --num-top-labels \t\tmaximum number of frequent labels to print [off]\n");
            fprintf(stderr, "\t   --discovery-fraction [FLOAT] fraction of labeled k-mers required for annotation [1.0]\n");
            fprintf(stderr, "\t   --labels-delimiter [STR]\tdelimiter for annotation labels [\":\"]\n");
            fprintf(stderr, "\t   --suppress-unlabeled \tdo not show results for sequences missing in graph [off]\n");
            // fprintf(stderr, "\t-d --distance [INT] \tmax allowed alignment distance [0]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case SERVER_QUERY: {
            fprintf(stderr, "Usage: %s server_query -i <GRAPH> -a <ANNOTATION> [options]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for server_query:\n");
            fprintf(stderr, "\t   --port [INT] \tTCP port for incoming connections [5555]\n");
            fprintf(stderr, "\t   --sparse \t\tuse the row-major sparse matrix to annotate graph [off]\n");
            // fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file []\n");
            // fprintf(stderr, "\t-d --distance [INT] \tmax allowed alignment distance [0]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tmaximum number of parallel connections [1]\n");
        } break;
        case CALL_VARIANTS: {
            fprintf(stderr, "Usage: %s call_variants -a <annotation> [options]\n", prog_name.c_str());

            fprintf(stderr, "Available options for call_variants:\n");
            fprintf(stderr, "\t-a --annotator [STR] \t\t\tannotator to load []\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --label-mask-in [STR] \t\tlabel to include in masked graph []\n");
            fprintf(stderr, "\t   --label-mask-out [STR] \t\tlabel to exclude from masked graph []\n");
            fprintf(stderr, "\t   --label-mask-out-fraction [FLOAT] \tmaximum fraction of mask-out labels among the set of\n");
            fprintf(stderr, "\t                                     \tall matching mask-in and mask-out labels [0.0]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --call-bubbles \t\t\tcall labels from bubbles\n");
            fprintf(stderr, "\t   --label-filter [STR] \t\tdiscard variants with this label []\n");
            fprintf(stderr, "\t   --taxonomy-map [STR] \t\tfilename of taxonomy map file []\n");
        } break;
        case PARSE_TAXONOMY: {
            fprintf(stderr, "Usage: %s parse_taxonomy -o <OUTBASE> [options]\n", prog_name.c_str());

            fprintf(stderr, "Available options for parse_taxonomy:\n");
            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --accession [STR] \tfilename of the accession2taxid.gz file []\n");
            fprintf(stderr, "\t   --taxonomy [STR] \tfilename of the nodes.dmp file []\n");
        } break;
    }

    fprintf(stderr, "\n\tGeneral options:\n");
    fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
    fprintf(stderr, "\t-h --help \t\tprint usage info\n");
    fprintf(stderr, "\n");
}
