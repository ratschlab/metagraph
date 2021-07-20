#include "config.hpp"

#include <cstring>
#include <iostream>
#include <unordered_set>
#include <filesystem>

#include "common/threads/threading.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "seq_io/formats.hpp"
#include "kmer/kmer_extractor.hpp"


namespace mtg {
namespace cli {

using mtg::graph::boss::BOSS;
using mtg::graph::DeBruijnGraph;


const size_t Config::kDefaultIndexSuffixLen
    = 20 / std::log2(kmer::KmerExtractor2Bit().alphabet.size());

void print_welcome_message() {
    fprintf(stderr, "#############################\n");
    fprintf(stderr, "### Welcome to MetaGraph! ###\n");
    fprintf(stderr, "#############################\n\n");
    fprintf(stderr, "Metagraph: comprehensive metagenome graph representation -- Version 0.1\n\n");
}

Config::Config(int argc, char *argv[]) {
    // provide help overview if no identity was given
    if (argc == 1) {
        print_usage(argv[0]);
        exit(-1);
    }

    // parse identity from first command line argument
    if (!strcmp(argv[1], "build")) {
        identity = BUILD;
    } else if (!strcmp(argv[1], "clean")) {
        identity = CLEAN;
    } else if (!strcmp(argv[1], "merge")) {
        identity = MERGE;
    } else if (!strcmp(argv[1], "extend")) {
        identity = EXTEND;
    } else if (!strcmp(argv[1], "concatenate")) {
        identity = CONCATENATE;
        clear_dummy = true;
    } else if (!strcmp(argv[1], "compare")) {
        identity = COMPARE;
    } else if (!strcmp(argv[1], "align")) {
        identity = ALIGN;
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
        tmp_dir = "OUTFBASE_TEMP_DIR";
    } else if (!strcmp(argv[1], "assemble")) {
        identity = ASSEMBLE;
    } else if (!strcmp(argv[1], "relax_brwt")) {
        identity = RELAX_BRWT;
    } else if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        print_welcome_message();
        print_usage(argv[0]);
        exit(0);
    } else {
        print_usage(argv[0]);
        exit(-1);
    }

    // provide help screen for chosen identity
    if (argc == 2) {
        print_usage(argv[0], identity);
        exit(-1);
    }

    const auto get_value = [&](int i) {
        assert(i > 0);
        assert(i < argc);

        if (i + 1 == argc) {
            std::cerr << "Error: no value provided for option "
                      << argv[i] << std::endl;
            print_usage(argv[0], identity);
            exit(-1);
        }
        return argv[i + 1];
    };

    bool print_usage_and_exit = false;

    // parse remaining command line items
    for (int i = 2; i < argc; ++i) {
        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
            common::set_verbose(true);
        } else if (!strcmp(argv[i], "--print")) {
            print_graph = true;
        } else if (!strcmp(argv[i], "--print-col-names")) {
            print_column_names = true;
        } else if (!strcmp(argv[i], "--print-internal")) {
            print_graph_internal_repr = true;
        } else if (!strcmp(argv[i], "--coordinates")) {
            coordinates = true;
        } else if (!strcmp(argv[i], "--num-kmers-in-seq")) {
            // FYI: experimental
            std::cerr << "WARNING: Flag --num-kmers-in-seq is experimental and"
                         " should only be used for experimental purposes" << std::endl;
            num_kmers_in_seq = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--count-kmers")) {
            count_kmers = true;
        } else if (!strcmp(argv[i], "--count-width")) {
            count_width = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--fwd-and-reverse")) {
            forward_and_reverse = true;
        } else if (!strcmp(argv[i], "--mode")) {
            graph_mode = string_to_graphmode(get_value(i++));
        } else if (!strcmp(argv[i], "--complete")) {
            complete = true;
        } else if (!strcmp(argv[i], "--dynamic")) {
            dynamic = true;
        } else if (!strcmp(argv[i], "--mask-dummy")) {
            mark_dummy_kmers = true;
        } else if (!strcmp(argv[i], "--anno-filename")) {
            filename_anno = true;
        } else if (!strcmp(argv[i], "--anno-header")) {
            annotate_sequence_headers = true;
        } else if (!strcmp(argv[i], "--header-comment-delim")) {
            fasta_anno_comment_delim = std::string(get_value(i++));
        } else if (!strcmp(argv[i], "--anno-label")) {
            anno_labels.emplace_back(get_value(i++));
        } else if (!strcmp(argv[i], "--coord-binsize")) {
            genome_binsize_anno = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--suppress-unlabeled")) {
            suppress_unlabeled = true;
        } else if (!strcmp(argv[i], "--sparse")) {
            sparse = true;
        } else if (!strcmp(argv[i], "--cache")) {
            num_columns_cached = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--fast")) {
            fast = true;
        } else if (!strcmp(argv[i], "--batch-size")) {
            query_batch_size_in_bytes = atoll(get_value(i++));
        } else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--parallel")) {
            set_num_threads(atoi(get_value(i++)));
        } else if (!strcmp(argv[i], "--parallel-nodes")) {
            parallel_nodes = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--max-path-length")) {
            max_path_length = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--parts-total")) {
            parts_total = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--part-idx")) {
            part_idx = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "--bins-per-thread")) {
            num_bins_per_thread = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "-k") || !strcmp(argv[i], "--kmer-length")) {
            k = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--min-count")) {
            min_count = std::max(atoi(get_value(i++)), 1);
        } else if (!strcmp(argv[i], "--max-count")) {
            max_count = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--min-count-q")) {
            min_count_quantile = std::max(std::stod(get_value(i++)), 0.);
        } else if (!strcmp(argv[i], "--max-count-q")) {
            max_count_quantile = std::min(std::stod(get_value(i++)), 1.);
        } else if (!strcmp(argv[i], "--count-bins-q")) {
            for (const auto &border : utils::split_string(get_value(i++), " ")) {
                count_slice_quantiles.push_back(std::stod(border));
            }
        } else if (!strcmp(argv[i], "--count-quantiles")) {
            for (const auto &p : utils::split_string(get_value(i++), " ")) {
                count_quantiles.push_back(std::stod(p));
            }
        } else if (!strcmp(argv[i], "--aggregate-columns")) {
            aggregate_columns = true;
        } else if (!strcmp(argv[i], "--intersected-anno")) {
            intersected_columns = get_value(i++);
        } else if (!strcmp(argv[i], "--min-fraction")) {
            min_fraction = std::stod(get_value(i++));
        } else if (!strcmp(argv[i], "--max-fraction")) {
            max_fraction = std::stod(get_value(i++));
        } else if (!strcmp(argv[i], "--mem-cap-gb")) {
            memory_available = atof(get_value(i++));
        } else if (!strcmp(argv[i], "--dump-text-anno")) {
            dump_text_anno = true;
        } else if (!strcmp(argv[i], "--discovery-fraction")) {
            discovery_fraction = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--align-rel-score-cutoff")) {
            alignment_rel_score_cutoff = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--query-presence")) {
            query_presence = true;
        } else if (!strcmp(argv[i], "--query-coords")) {
            query_coords = true;
        } else if (!strcmp(argv[i], "--filter-present")) {
            filter_present = true;
        } else if (!strcmp(argv[i], "--count-labels")) {
            count_labels = true;
        } else if (!strcmp(argv[i], "--print-signature")) {
            print_signature = true;
        } else if (!strcmp(argv[i], "--map")) {
            map_sequences = true;
        } else if (!strcmp(argv[i], "--align")) {
            align_sequences = true;
        } else if (!strcmp(argv[i], "--align-one-strand")) {
            align_one_strand = true;
        } else if (!strcmp(argv[i], "--align-edit-distance")) {
            alignment_edit_distance = true;
        } else if (!strcmp(argv[i], "--max-hull-depth")) {
            max_hull_depth = atoll(get_value(i++));
        } else if (!strcmp(argv[i], "--batch-align")) {
            batch_align = true;
        } else if (!strcmp(argv[i], "--align-length")) {
            alignment_length = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-match-score")) {
            alignment_match_score = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-mm-transition-penalty")) {
            alignment_mm_transition_score = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-mm-transversion-penalty")) {
            alignment_mm_transversion_score = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-gap-open-penalty")) {
            alignment_gap_opening_penalty = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-gap-extension-penalty")) {
            alignment_gap_extension_penalty = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-alternative-alignments")) {
            alignment_num_alternative_paths = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-min-path-score")) {
            alignment_min_path_score = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-xdrop")) {
            alignment_xdrop = atol(get_value(i++));
        } else if (!strcmp(argv[i], "--align-min-seed-length")) {
            alignment_min_seed_length = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-max-seed-length")) {
            alignment_max_seed_length = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-max-num-seeds-per-locus")) {
            alignment_max_num_seeds_per_locus = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-max-nodes-per-seq-char")) {
            alignment_max_nodes_per_seq_char = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--align-min-exact-match")) {
            alignment_min_exact_match = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--max-hull-forks")) {
            max_hull_forks = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-max-ram")) {
            alignment_max_ram = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "--frequency")) {
            frequency = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--distance")) {
            distance = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--outfile-base")) {
            outfbase = std::string(get_value(i++));
        } else if (!strcmp(argv[i], "--reference")) {
            refpath = std::string(get_value(i++));
        } else if (!strcmp(argv[i], "--header-delimiter")) {
            fasta_header_delimiter = std::string(get_value(i++));
        } else if (!strcmp(argv[i], "--labels-delimiter")) {
            anno_labels_delimiter = std::string(get_value(i++));
        } else if (!strcmp(argv[i], "--separately")) {
            separately = true;
        } else if (!strcmp(argv[i], "--sequentially")) {
            files_sequentially = true;
        } else if (!strcmp(argv[i], "--num-top-labels")) {
            num_top_labels = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--port")) {
            port = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--address")) {
            host_address = get_value(i++);
        }else if (!strcmp(argv[i], "--suffix")) {
            suffix = get_value(i++);
        } else if (!strcmp(argv[i], "--initialize-bloom")) {
            initialize_bloom = true;
        } else if (!strcmp(argv[i], "--bloom-fpp")) {
            bloom_fpp = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--bloom-bpk")) {
            bloom_bpk = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--bloom-max-num-hash-functions")) {
            bloom_max_num_hash_functions = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--state")) {
            state = string_to_state(get_value(i++));

        } else if (!strcmp(argv[i], "--anno-type")) {
            anno_type = string_to_annotype(get_value(i++));
        } else if (!strcmp(argv[i], "--graph")) {
            graph_type = string_to_graphtype(get_value(i++));
        } else if (!strcmp(argv[i], "--rename-cols")) {
            rename_instructions_file = std::string(get_value(i++));
        } else if (!strcmp(argv[i], "-a") || !strcmp(argv[i], "--annotator")) {
            infbase_annotators.emplace_back(get_value(i++));
        } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--infile-base")) {
            infbase = std::string(get_value(i++));
        } else if (!strcmp(argv[i], "--to-adj-list")) {
            to_adj_list = true;
        } else if (!strcmp(argv[i], "--to-fasta")) {
            to_fasta = true;
        } else if (!strcmp(argv[i], "--enumerate")) {
            enumerate_out_sequences = true;
        } else if (!strcmp(argv[i], "--to-gfa")) {
            to_gfa = true;
        } else if (!strcmp(argv[i], "--compacted")) {
            output_compacted = true;
        } else if (!strcmp(argv[i], "--json")) {
            output_json = true;
        } else if (!strcmp(argv[i], "--unitigs")) {
            to_fasta = true;
            unitigs = true;
        } else if (!strcmp(argv[i], "--primary-kmers")) {
            kmers_in_single_form = true;
        } else if (!strcmp(argv[i], "--header")) {
            header = std::string(get_value(i++));
        } else if (!strcmp(argv[i], "--prune-tips")) {
            min_tip_size = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--prune-unitigs")) {
            min_unitig_median_kmer_abundance = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--fallback")) {
            fallback_abundance_cutoff = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--smoothing-window")) {
            smoothing_window = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--num-singletons")) {
            num_singleton_kmers = atoll(get_value(i++));
        } else if (!strcmp(argv[i], "--count-dummy")) {
            count_dummy = true;
        } else if (!strcmp(argv[i], "--clear-dummy")) {
            clear_dummy = true;
        } else if (!strcmp(argv[i], "--index-ranges")) {
            node_suffix_length = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--no-postprocessing")) {
            clear_dummy = false;
        } else if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--len-suffix")) {
            suffix_len = atoi(get_value(i++));
        //} else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--threads")) {
        //    num_threads = atoi(get_value(i++));
        //} else if (!strcmp(argv[i], "--debug")) {
        //    debug = true;
        } else if (!strcmp(argv[i], "--greedy")) {
            greedy_brwt = true;
        } else if (!strcmp(argv[i], "--row-diff-stage")) {
            row_diff_stage = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--linkage")) {
            cluster_linkage = true;
        } else if (!strcmp(argv[i], "--subsample")) {
            num_rows_subsampled = atoll(get_value(i++));
        } else if (!strcmp(argv[i], "--linkage-file")) {
            linkage_file = get_value(i++);
        } else if (!strcmp(argv[i], "--arity")) {
            arity_brwt = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--relax-arity")) {
            relax_arity_brwt = atoi(get_value(i++));
        // } else if (!strcmp(argv[i], "--cache-size")) {
        //     row_cache_size = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_welcome_message();
            print_usage(argv[0], identity);
            exit(0);
        } else if (!strcmp(argv[i], "--label-mask-in")) {
            label_mask_in.emplace_back(get_value(i++));
        } else if (!strcmp(argv[i], "--label-mask-out")) {
            label_mask_out.emplace_back(get_value(i++));
        } else if (!strcmp(argv[i], "--label-mask-in-fraction")) {
            label_mask_in_fraction = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--label-mask-out-fraction")) {
            label_mask_out_fraction = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--label-other-fraction")) {
            label_other_fraction = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--filter-by-kmer")) {
            filter_by_kmer = true;
        } else if (!strcmp(argv[i], "--disk-swap")) {
            tmp_dir = get_value(i++);
        } else if (!strcmp(argv[i], "--disk-cap-gb")) {
            disk_cap_bytes = atoi(get_value(i++)) * 1e9;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "\nERROR: Unknown option %s\n\n", argv[i]);
            print_usage(argv[0], identity);
            exit(-1);
        } else {
            fnames.push_back(argv[i]);
        }
    }

    if (parallel_nodes == static_cast<unsigned int>(-1))
        parallel_nodes = get_num_threads();

    if (identity == TRANSFORM && to_fasta)
        identity = ASSEMBLE;

    // given kmc_pre and kmc_suf pair, only include one
    // this still allows for the same file to be included multiple times
    std::unordered_set<std::string> kmc_file_set;

    for (auto it = fnames.begin(); it != fnames.end(); ++it) {
        if (seq_io::file_format(*it) == "KMC"
                && !kmc_file_set.insert(utils::remove_suffix(*it, ".kmc_pre", ".kmc_suf")).second)
            fnames.erase(it--);
    }

    if (!fnames.size() && identity != STATS
                      && identity != SERVER_QUERY
                      && !(identity == BUILD && complete)
                      && !(identity == CONCATENATE && !infbase.empty())) {
        std::string line;
        while (std::getline(std::cin, line)) {
            if (line.size())
                fnames.push_back(line);
        }
    }

    if (!count_slice_quantiles.size()) {
        count_slice_quantiles.push_back(0);
        count_slice_quantiles.push_back(1);
    }

    if (count_width <= 1) {
        std::cerr << "Error: bad value for count-width, need at least 2 bits"
                     " to represent k-mer abundance" << std::endl;
        print_usage_and_exit = true;
    }
    if (!count_kmers)
        count_width = 0;

    if (count_width > 32) {
        std::cerr << "Error: bad value for count-width, can use maximum 32 bits"
                     " to represent k-mer abundance" << std::endl;
        print_usage_and_exit = true;
    }

    for (size_t i = 1; i < count_slice_quantiles.size(); ++i) {
        if (count_slice_quantiles[i - 1] >= count_slice_quantiles[i]) {
            std::cerr << "Error: bin count quantiles must be provided in strictly increasing order"
                      << std::endl;
            print_usage_and_exit = true;
        }
    }
    if (count_slice_quantiles.front() < 0 || count_slice_quantiles.back() > 1) {
        std::cerr << "Error: bin count quantiles must be in range [0, 1]"
                  << std::endl;
        print_usage_and_exit = true;
    }
    if (count_slice_quantiles.size() == 1) {
        std::cerr << "Error: provide at least two bin count borders"
                  << std::endl;
        print_usage_and_exit = true;
    }

    if (min_fraction < 0 || min_fraction > 1 || max_fraction < 0 || max_fraction > 1) {
        std::cerr << "Error: min_fraction and max_fraction must be in range [0, 1]"
                  << std::endl;
        print_usage_and_exit = true;
    }

#if _PROTEIN_GRAPH
    if (graph_mode != DeBruijnGraph::BASIC || forward_and_reverse) {
        std::cerr << "Error: reverse complement not defined for protein alphabets"
                  << std::endl;
        print_usage_and_exit = true;
    }
#endif

    if (tmp_dir == "OUTFBASE_TEMP_DIR") {
        tmp_dir = std::filesystem::path(outfbase).remove_filename();
    }
    utils::set_swap_path(tmp_dir);

    if (identity != CONCATENATE
            && identity != STATS
            && identity != SERVER_QUERY
            && !(identity == TRANSFORM_ANNOTATION && anno_type == Config::RowDiff)
            && !(identity == BUILD && complete)
            && !fnames.size()) {
        std::cerr << "Error: No input file(s) passed" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == CONCATENATE && !(fnames.empty() ^ infbase.empty())) {
        std::cerr << "Error: Either set all chunk filenames"
                  << " or use the -i and -l options" << std::endl;
        print_usage_and_exit = true;
    }

    if (alignment_min_seed_length > alignment_max_seed_length) {
        std::cerr << "Error: min_seed_length must be <= max_seed_length" << std::endl;
        print_usage_and_exit = true;
    }

    // only the best alignment is used in query
    // |alignment_num_alternative_paths| must be set to 1
    if (identity == QUERY && align_sequences
                          && alignment_num_alternative_paths != 1)
        print_usage_and_exit = true;

    if (identity == ALIGN && infbase.empty())
        print_usage_and_exit = true;

    if (identity == ALIGN &&
            (alignment_mm_transition_score < 0
            || alignment_mm_transversion_score < 0
            || alignment_gap_opening_penalty < 0
            || alignment_gap_extension_penalty < 0)) {
        std::cerr << "Error: alignment penalties should be given as positive integers"
                  << std::endl;
        print_usage_and_exit = true;
    }

    if (count_kmers || query_presence)
        map_sequences = true;

    if ((identity == QUERY || identity == SERVER_QUERY) && infbase.empty())
        print_usage_and_exit = true;

    if ((identity == QUERY || identity == SERVER_QUERY || identity == ALIGN)
            && alignment_num_alternative_paths == 0) {
        std::cerr << "Error: align-alternative-alignments must be > 0" << std::endl;
        print_usage_and_exit = true;
    }

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
            && !filename_anno && !annotate_sequence_headers && !anno_labels.size()) {
        std::cerr << "Error: No annotation to add" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == ANNOTATE_COORDINATES && outfbase.empty())
        outfbase = utils::remove_suffix(infbase, ".dbg",
                                                 ".orhashdbg",
                                                 ".hashstrdbg",
                                                 ".hashfastdbg",
                                                 ".bitmapdbg");
    if (identity == EXTEND && infbase.empty())
        print_usage_and_exit = true;

    if ((identity == QUERY || identity == SERVER_QUERY) && infbase_annotators.size() != 1)
        print_usage_and_exit = true;

    if ((identity == TRANSFORM
            || identity == CLEAN
            || identity == ASSEMBLE
            || identity == RELAX_BRWT)
                    && fnames.size() != 1)
        print_usage_and_exit = true;

    if ((identity == TRANSFORM
            || identity == BUILD
            || identity == ANNOTATE
            || identity == CONCATENATE
            || identity == EXTEND
            || identity == MERGE
            || identity == CLEAN
            || identity == TRANSFORM_ANNOTATION
            || identity == MERGE_ANNOTATIONS
            || identity == ASSEMBLE
            || identity == RELAX_BRWT)
                    && outfbase.empty())
        print_usage_and_exit = true;

    if (identity == TRANSFORM_ANNOTATION) {
        const bool to_row_diff = anno_type == RowDiff
                                    || anno_type == RowDiffBRWT
                                    || anno_type == IntRowDiffBRWT
                                    || anno_type == RowDiffRowSparse
                                    || anno_type == RowDiffCoord;
        if (to_row_diff && !infbase.size()) {
            std::cerr << "Path to graph must be passed with '-i <GRAPH>'" << std::endl;
            print_usage_and_exit = true;
        } else if (!to_row_diff && infbase.size()) {
            std::cerr << "Graph is only required for transform to row_diff types" << std::endl;
            print_usage_and_exit = true;
        }
    }

    if (identity == MERGE && fnames.size() < 2)
        print_usage_and_exit = true;

    if (identity == COMPARE && fnames.size() != 2)
        print_usage_and_exit = true;

    if (discovery_fraction < 0 || discovery_fraction > 1)
        print_usage_and_exit = true;

    if (min_count >= max_count) {
        std::cerr << "Error: max-count must be greater than min-count" << std::endl;
        print_usage(argv[0], identity);
    }

    if (alignment_max_seed_length < alignment_min_seed_length) {
        std::cerr << "Error: align-max-seed-length has to be at least align-min-seed-length" << std::endl;
        print_usage_and_exit = true;
    }

    if (bloom_fpp <= 0.0 || bloom_fpp > 1.0) {
        std::cerr << "Error: bloom-fpp must be > 0.0 and < 1.0" << std::endl;
        print_usage_and_exit = true;
    }

    if (bloom_bpk <= 0.0) {
        std::cerr << "Error: bloom-bpk must > 0.0" << std::endl;
        print_usage_and_exit = true;
    }

    if (initialize_bloom && bloom_bpk == 0.0 && bloom_fpp == 1.0) {
        std::cerr << "Error: at least one of 0.0 < bloom_fpp < 1.0 or 0.0 < bloom_bpk must be true" << std::endl;
        print_usage_and_exit = true;
    }

    if (outfbase.size()
            && !(utils::check_if_writable(outfbase)
                    || (separately
                        && std::filesystem::is_directory(std::filesystem::status(outfbase))))) {
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


std::string Config::state_to_string(BOSS::State state) {
    switch (state) {
        case BOSS::State::STAT:
            return "stat";
        case BOSS::State::DYN:
            return "dynamic";
        case BOSS::State::SMALL:
            return "small";
        case BOSS::State::FAST:
            return "fast";
    }
    throw std::runtime_error("Never happens");
}

BOSS::State Config::string_to_state(const std::string &string) {
    if (string == "stat") {
        return BOSS::State::STAT;
    } else if (string == "dynamic") {
        return BOSS::State::DYN;
    } else if (string == "small") {
        return BOSS::State::SMALL;
    } else if (string == "fast") {
        return BOSS::State::FAST;
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
        case RbBRWT:
            return "rb_brwt";
        case RowDiff:
            return "row_diff";
        case RowDiffBRWT:
            return "row_diff_brwt";
        case RowDiffRowSparse:
            return "row_diff_sparse";
        case RowSparse:
            return "row_sparse";
        case IntBRWT:
            return "int_brwt";
        case IntRowDiffBRWT:
            return "row_diff_int_brwt";
        case ColumnCoord:
            return "column_coord";
        case RowDiffCoord:
            return "row_diff_coord";
    }
    throw std::runtime_error("Never happens");
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
    } else if (string == "rb_brwt") {
        return AnnotationType::RbBRWT;
    } else if (string == "row_diff") {
        return AnnotationType::RowDiff;
    } else if (string == "row_diff_brwt") {
        return AnnotationType::RowDiffBRWT;
    } else if (string == "row_diff_sparse") {
        return AnnotationType::RowDiffRowSparse;
    } else if (string == "row_sparse") {
        return AnnotationType::RowSparse;
    } else if (string == "int_brwt") {
        return AnnotationType::IntBRWT;
    } else if (string == "row_diff_int_brwt") {
        return AnnotationType::IntRowDiffBRWT;
    } else if (string == "column_coord") {
        return AnnotationType::ColumnCoord;
    } else if (string == "row_diff_coord") {
        return AnnotationType::RowDiffCoord;
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

    } else if (string == "hashpacked") {
        return GraphType::HASH_PACKED;

    } else if (string == "hashstr") {
        return GraphType::HASH_STR;

    } else if (string == "hashfast") {
        return GraphType::HASH_FAST;

    } else if (string == "bitmap") {
        return GraphType::BITMAP;

    } else {
        std::cerr << "Error: unknown graph representation" << std::endl;
        exit(1);
    }
}

std::string Config::graphmode_to_string(DeBruijnGraph::Mode mode) {
    switch (mode) {
        case DeBruijnGraph::BASIC:
            return "basic";
        case DeBruijnGraph::CANONICAL:
            return "canonical";
        case DeBruijnGraph::PRIMARY:
            return "primary";
    }
    throw std::runtime_error("Never happens");
}

DeBruijnGraph::Mode Config::string_to_graphmode(const std::string &string) {
    if (string == "basic") {
        return DeBruijnGraph::BASIC;

    } else if (string == "canonical") {
        return DeBruijnGraph::CANONICAL;

    } else if (string == "primary") {
        return DeBruijnGraph::PRIMARY;

    } else {
        std::cerr << "Error: unknown graph mode" << std::endl;
        exit(1);
    }
}


void Config::print_usage(const std::string &prog_name, IdentityType identity) {
    const char annotation_list[] = "\t\t( column, brwt, rb_brwt, int_brwt,\n"
                                   "\t\t  column_coord,\n"
                                   "\t\t  row_diff, row_diff_brwt, row_diff_sparse, row_diff_int_brwt,\n"
                                   "\t\t  row, flat, row_sparse, rbfish, bin_rel_wt, bin_rel_wt_sdsl )";

    switch (identity) {
        case NO_IDENTITY: {
            fprintf(stderr, "Usage: %s <command> [command specific options]\n\n", prog_name.c_str());

            fprintf(stderr, "Available commands:\n");

            fprintf(stderr, "\tbuild\t\tconstruct a graph object from input sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats into a given graph\n\n");

            fprintf(stderr, "\tclean\t\tclean an existing graph and extract sequences from it\n");
            fprintf(stderr, "\t\t\tin fast[a|q] formats\n\n");

            fprintf(stderr, "\textend\t\textend an existing graph with new sequences from\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats\n\n");

            fprintf(stderr, "\tmerge\t\tintegrate a given set of graph structures\n");
            fprintf(stderr, "\t\t\tand output a new graph structure\n\n");

            fprintf(stderr, "\tconcatenate\tcombine the results of the external merge or\n");
            fprintf(stderr, "\t\t\tconstruction and output the resulting graph structure\n\n");

            fprintf(stderr, "\tcompare\t\tcheck whether two given graphs are identical\n\n");

            fprintf(stderr, "\talign\t\talign sequences provided in fast[a|q] files to graph\n\n");

            fprintf(stderr, "\tstats\t\tprint graph statistics for given graph(s)\n\n");

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

            return;
        }
        case BUILD: {
            fprintf(stderr, "Usage: %s build [options] -o <outfile-base> FILE1 [[FILE2] ...]\n"
                            "\tEach input file is given in FASTA, FASTQ, VCF, or KMC format.\n"
                            "\tNote that VCF files must be in plain text or bgzip format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for build:\n");
            fprintf(stderr, "\t   --min-count [INT] \tmin k-mer abundance, including [1]\n");
            fprintf(stderr, "\t   --max-count [INT] \tmax k-mer abundance, excluding [inf]\n");
            fprintf(stderr, "\t   --min-count-q [INT] \tmin k-mer abundance quantile (min-count is used by default) [0.0]\n");
            fprintf(stderr, "\t   --max-count-q [INT] \tmax k-mer abundance quantile (max-count is used by default) [1.0]\n");
            fprintf(stderr, "\t   --reference [STR] \tbasename of reference sequence (for parsing VCF files) []\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --graph [STR] \tgraph representation: succinct / bitmap / hash / hashstr / hashfast [succinct]\n");
            fprintf(stderr, "\t   --count-kmers \tcount k-mers and build weighted graph [off]\n");
            fprintf(stderr, "\t   --count-width \tnumber of bits used to represent k-mer abundance [8]\n");
            fprintf(stderr, "\t   --index-ranges [INT]\tindex all node ranges in BOSS for suffixes of given length [%zu]\n", kDefaultIndexSuffixLen);
            fprintf(stderr, "\t-k --kmer-length [INT] \tlength of the k-mer to use [3]\n");
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --mode \t\tk-mer indexing mode: basic / canonical / primary [basic]\n");
#endif
            fprintf(stderr, "\t   --complete \t\tconstruct a complete graph (only for Bitmap graph) [off]\n");
            fprintf(stderr, "\t   --mem-cap-gb [INT] \tpreallocated buffer size in GB [1]\n");
            fprintf(stderr, "\t   --dynamic \t\tuse dynamic build method [off]\n");
            fprintf(stderr, "\t-l --len-suffix [INT] \tk-mer suffix length for building graph from chunks [0]\n");
            fprintf(stderr, "\t   --suffix \t\tbuild graph chunk only for k-mers with the suffix given [off]\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\tbasename of output file []\n");
            fprintf(stderr, "\t   --mask-dummy \tbuild mask for dummy k-mers (only for Succinct graph) [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
            fprintf(stderr, "\t   --disk-swap [STR] \tdirectory to use for temporary files [off]\n");
            fprintf(stderr, "\t   --disk-cap-gb [INT] \tmax temp disk space to use before forcing a merge, in GB [inf]\n");
        } break;
        case CLEAN: {
            fprintf(stderr, "Usage: %s clean -o <outfile-base> [options] GRAPH\n\n", prog_name.c_str());
            fprintf(stderr, "Available options for clean:\n");
            fprintf(stderr, "\t   --min-count [INT] \t\tmin k-mer abundance, including [1]\n");
            fprintf(stderr, "\t   --max-count [INT] \t\tmax k-mer abundance, excluding [inf]\n");
            fprintf(stderr, "\t   --num-singletons [INT] \treset the number of count 1 k-mers in histogram (0: off) [0]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --prune-tips [INT] \t\tprune all dead ends shorter than this value [1]\n");
            fprintf(stderr, "\t   --prune-unitigs [INT] \tprune all unitigs with median k-mer counts smaller\n"
                            "\t                         \t\tthan this value (0: auto) [1]\n");
            fprintf(stderr, "\t   --fallback [INT] \t\tfallback threshold if the automatic one cannot be\n"
                            "\t                         \t\tdetermined (-1: disables fallback) [1]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --smoothing-window [INT] \twindow size for smoothing k-mer counts in unitigs [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --count-bins-q [FLOAT ...] \tbinning quantiles for partitioning k-mers with\n"
                            "\t                              \t\tdifferent abundance levels ['0 1']\n"
                            "\t                              \t\tExample: --count-bins-q '0 0.33 0.66 1'\n");
            // fprintf(stderr, "\n");
            // fprintf(stderr, "\t-o --outfile-base [STR]\tbasename of output file []\n");
            fprintf(stderr, "\t   --unitigs \t\t\textract unitigs instead of contigs [off]\n");
            fprintf(stderr, "\t   --to-fasta \t\t\tdump clean sequences to compressed FASTA file [off]\n");
            fprintf(stderr, "\t   --enumerate \t\t\tenumerate sequences in FASTA [off]\n");
            // fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case EXTEND: {
            fprintf(stderr, "Usage: %s extend -i <GRAPH> -o <extended_graph_basename> [options] FILE1 [[FILE2] ...]\n"
                            "\tEach input file is given in FASTA, FASTQ, VCF, or KMC format.\n"
                            "\tNote that VCF files must be in plain text or bgzip format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for extend:\n");
            fprintf(stderr, "\t   --min-count [INT] \tmin k-mer abundance, including [1]\n");
            fprintf(stderr, "\t   --max-count [INT] \tmax k-mer abundance, excluding [inf]\n");
            fprintf(stderr, "\t   --reference [STR] \tbasename of reference sequence (for parsing VCF files) []\n");
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --fwd-and-reverse \tadd both forward and reverse complement sequences [off]\n");
#endif
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-a --annotator [STR] \tannotator to extend []\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\tbasename of output file []\n");
            // fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case ALIGN: {
            fprintf(stderr, "Usage: %s align -i <GRAPH> [options] FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());

#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --fwd-and-reverse \t\tfor each input sequence, align its reverse complement as well [off]\n");
#endif
            fprintf(stderr, "\t   --header-comment-delim [STR]\tdelimiter for joining fasta header with comment [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --map \t\t\tmap k-mers to graph exactly instead of aligning.\n");
            fprintf(stderr, "\t         \t\t\t\tTurned on if --count-kmers or --query-presence are set [off]\n");
            fprintf(stderr, "\t   --compacted\t\t\tdump the GFA's 'P' lines in a compacted mode [off]\n");
            fprintf(stderr, "\t-k --kmer-length [INT]\t\tlength of mapped k-mers (at most graph's k) [k]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --count-kmers \t\tfor each sequence, report the number of k-mers discovered in graph [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --query-presence \t\ttest sequences for presence, report as 0 or 1 [off]\n");
            fprintf(stderr, "\t   --filter-present \t\treport only present input sequences as FASTA [off]\n");
            fprintf(stderr, "\t   --batch-size \t\tquery batch size (number of base pairs) [100000000]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Available options for alignment:\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\t\t\t\tbasename of output file []\n");
            fprintf(stderr, "\t   --json \t\t\t\t\toutput alignment in JSON format [off]\n");
            fprintf(stderr, "\t   --align-one-strand \t\t\t\tdo not align backwards from a seed on basic-mode graphs [off]\n");
            fprintf(stderr, "\t   --align-alternative-alignments \t\tthe number of alternative paths to report per seed [1]\n");
            fprintf(stderr, "\t   --align-min-path-score [INT]\t\t\tthe minimum score that a reported path can have [0]\n");
            fprintf(stderr, "\t   --align-max-nodes-per-seq-char [FLOAT]\tmaximum number of nodes to consider per sequence character [12.0]\n");
            fprintf(stderr, "\t   --align-max-ram [FLOAT]\t\t\tmaximum amount of RAM used per alignment in MB [200.0]\n");
            fprintf(stderr, "\t   --align-xdrop [INT]\t\t\t\tthe maximum difference between the current and the best alignment [27]\n");
            fprintf(stderr, "\t   --align-rel-score-cutoff [FLOAT]\t\tmin score relative to the current best alignment to use as a lower bound for subsequent extensions [0.8]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Advanced options for scoring:\n");
            fprintf(stderr, "\t   --align-match-score [INT]\t\t\tpositive match score [2]\n");
            fprintf(stderr, "\t   --align-mm-transition-penalty [INT]\t\tpositive transition penalty (DNA only) [3]\n");
            fprintf(stderr, "\t   --align-mm-transversion-penalty [INT]\tpositive transversion penalty (DNA only) [3]\n");
            fprintf(stderr, "\t   --align-gap-open-penalty [INT]\t\tpositive gap opening penalty [5]\n");
            fprintf(stderr, "\t   --align-gap-extension-penalty [INT]\t\tpositive gap extension penalty [2]\n");
            fprintf(stderr, "\t   --align-edit-distance \t\t\tuse unit costs for scoring matrix [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Advanced options for seeding:\n");
            fprintf(stderr, "\t   --align-min-seed-length [INT]\t\tthe minimum length of a seed [graph k]\n");
            fprintf(stderr, "\t   --align-max-seed-length [INT]\t\tthe maximum length of a seed [graph k]\n");
            fprintf(stderr, "\t   --align-min-exact-match [FLOAT] \t\tfraction of matching nucleotides required to align sequence [0.0]\n");
            fprintf(stderr, "\t   --align-max-num-seeds-per-locus [INT]\tthe maximum number of allowed inexact seeds per locus [inf]\n");
        } break;
        case COMPARE: {
            fprintf(stderr, "Usage: %s compare [options] GRAPH1 GRAPH2\n\n", prog_name.c_str());

            // fprintf(stderr, "Available options for compare:\n");
            // fprintf(stderr, "\t   --internal \t\tcompare internal graph representations\n");
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
            fprintf(stderr, "\t-l --len-suffix [INT] \titerate all possible suffixes of the length given [0]\n");
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --mode \t\tk-mer indexing mode: basic / canonical / primary [basic]\n");
#endif
            fprintf(stderr, "\t   --no-postprocessing \tdo not erase redundant dummy edges after concatenation [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case TRANSFORM: {
            fprintf(stderr, "Usage: %s transform -o <outfile-base> [options] GRAPH\n\n", prog_name.c_str());

            // fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --index-ranges [INT]\tindex all node ranges in BOSS for suffixes of given length [%zu]\n", kDefaultIndexSuffixLen);
            fprintf(stderr, "\t   --clear-dummy \terase all redundant dummy edges and build an edgemask for non-redundant [off]\n");
            fprintf(stderr, "\t   --prune-tips [INT] \tprune all dead ends of this length and shorter [0]\n");
            fprintf(stderr, "\t   --state [STR] \tchange state of succinct graph: small / dynamic / fast [stat]\n");
            fprintf(stderr, "\t   --to-adj-list \twrite adjacency list to file [off]\n");
            fprintf(stderr, "\t   --to-fasta \t\textract sequences from graph and dump to compressed FASTA file [off]\n");
            fprintf(stderr, "\t   --enumerate \t\tenumerate sequences in FASTA [off]\n");
            fprintf(stderr, "\t   --initialize-bloom \tconstruct a Bloom filter for faster detection of non-existing k-mers [off]\n");
            fprintf(stderr, "\t   --unitigs \t\textract all unitigs from graph and dump to compressed FASTA file [off]\n");
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --primary-kmers \toutput each k-mer only in one if its forms (canonical/non-canonical) [off]\n");
#endif
            fprintf(stderr, "\t   --to-gfa \t\tdump graph layout to GFA [off]\n");
            fprintf(stderr, "\t   --compacted \t\tdump compacted de Bruijn graph to GFA [off]\n");
            fprintf(stderr, "\t   --header [STR] \theader for sequences in FASTA output []\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Advanced options for --initialize-bloom. bloom-fpp, when < 1, overrides bloom-bpk.\n");
            fprintf(stderr, "\t   --bloom-fpp [FLOAT] \t\t\t\texpected false positive rate [1.0]\n");
            fprintf(stderr, "\t   --bloom-bpk [FLOAT] \t\t\t\tnumber of bits per kmer [4.0]\n");
            fprintf(stderr, "\t   --bloom-max-num-hash-functions [INT] \tmaximum number of hash functions [10]\n");
        } break;
        case ASSEMBLE: {
            fprintf(stderr, "Usage: %s assemble -o <outfile-base> [options] GRAPH\n"
                            "\tAssemble contigs from de Bruijn graph and dump to compressed FASTA file.\n\n", prog_name.c_str());

            // fprintf(stderr, "\t-o --outfile-base [STR] \t\tbasename of output file []\n");
            fprintf(stderr, "\t   --prune-tips [INT] \tprune all dead ends of this length and shorter [0]\n");
            fprintf(stderr, "\t   --unitigs \t\textract unitigs [off]\n");
            fprintf(stderr, "\t   --enumerate \t\tenumerate sequences assembled and dumped to FASTA [off]\n");
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --primary-kmers \toutput each k-mer only in one if its forms (canonical/non-canonical) [off]\n");
#endif
            fprintf(stderr, "\t   --to-gfa \t\tdump graph layout to GFA [off]\n");
            fprintf(stderr, "\t   --compacted \t\tdump compacted de Bruijn graph to GFA [off]\n");
            fprintf(stderr, "\t   --header [STR] \theader for sequences in FASTA output []\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-a --annotator [STR] \t\t\tannotator to load []\n");
            fprintf(stderr, "\t   --label-mask-in [STR] \t\tlabel to include in masked graph\n");
            fprintf(stderr, "\t   --label-mask-out [STR] \t\tlabel to exclude from masked graph\n");
            fprintf(stderr, "\t   --label-mask-in-fraction [FLOAT] \tminimum fraction of mask-in labels among the set of masked labels [1.0]\n");
            fprintf(stderr, "\t   --label-mask-out-fraction [FLOAT] \tmaximum fraction of mask-out labels among the set of masked labels [0.0]\n");
            fprintf(stderr, "\t   --label-other-fraction [FLOAT] \tmaximum fraction of other labels allowed [1.0]\n");
            fprintf(stderr, "\t   --filter-by-kmer \t\t\tmask out graph k-mers individually [off]\n");
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
            fprintf(stderr, "Usage: %s annotate -i <GRAPH> -o <annotation-basename> [options] FILE1 [[FILE2] ...]\n"
                            "\tEach file is given in FASTA, FASTQ, VCF, or KMC format.\n"
                            "\tNote that VCF files must be in plain text or bgzip format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t   --min-count [INT] \tmin k-mer abundance, including [1]\n");
            fprintf(stderr, "\t   --max-count [INT] \tmax k-mer abundance, excluding [inf]\n");
            fprintf(stderr, "\t   --reference [STR] \tbasename of reference sequence (for parsing VCF files) []\n");
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --fwd-and-reverse \tprocess both forward and reverse complement sequences [off]\n");
#endif
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --anno-type [STR] \ttarget annotation representation: column / row [column]\n");
            fprintf(stderr, "\t-a --annotator [STR] \tannotator to update []\n");
            fprintf(stderr, "\t   --sparse \t\tuse the row-major sparse matrix to annotate graph [off]\n");
            fprintf(stderr, "\t   --cache \t\tnumber of columns in cache (for column representation only) [10]\n");
            fprintf(stderr, "\t   --disk-swap [STR] \tdirectory to use for temporary files [off]\n");
            fprintf(stderr, "\t   --mem-cap-gb [FLOAT]\tbuffer size in GB (per column in construction) [1]\n");
            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file (or directory, for --separately) []\n");
            fprintf(stderr, "\t   --separately \tannotate each file independently and dump to the same directory [off]\n");
            fprintf(stderr, "\t   --sequentially \tannotate files sequentially (each may use multiple threads) [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --anno-filename \t\tinclude filenames as annotation labels [off]\n");
            fprintf(stderr, "\t   --anno-header \t\textract annotation labels from headers of sequences in files [off]\n");
            fprintf(stderr, "\t   --header-comment-delim [STR]\tdelimiter for joining fasta header with comment [off]\n");
            fprintf(stderr, "\t   --header-delimiter [STR]\tdelimiter for splitting annotation header into multiple labels [off]\n");
            fprintf(stderr, "\t   --anno-label [STR]\t\tadd label to annotation for all sequences from the files passed []\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --count-kmers \tadd k-mer counts to the annotation [off]\n");
            fprintf(stderr, "\t   --count-width \tnumber of bits used to represent k-mer abundance [8]\n");
            fprintf(stderr, "\t   --coordinates \tannotate coordinates as multi-integer attributes [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
            // fprintf(stderr, "\t   --fast \t\t\tannotate in fast regime (leads to repeated labels and bigger annotation) [off]\n");
        } break;
        case ANNOTATE_COORDINATES: {
            fprintf(stderr, "Usage: %s coordinate -i <GRAPH> [options] FASTA1 [[FASTA2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --fwd-and-reverse \t\tprocess both forward and reverse complement sequences [off]\n");
#endif
            fprintf(stderr, "\t-a --annotator [STR] \t\tannotator to update []\n");
            fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file [<GRAPH>]\n");
            fprintf(stderr, "\t   --coord-binsize [INT]\tstepsize for k-mer coordinates in input sequences from the fasta files [1000]\n");
            fprintf(stderr, "\t   --fast \t\t\tannotate in fast regime [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case MERGE_ANNOTATIONS: {
            fprintf(stderr, "Usage: %s merge_anno -o <annotation-basename> [options] ANNOT1 [[ANNOT2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t   --anno-type [STR] \ttarget annotation representation [column]\n");
            fprintf(stderr, "%s\n", annotation_list);
            // fprintf(stderr, "\t   --sparse \t\tuse the row-major sparse matrix to annotate graph [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case TRANSFORM_ANNOTATION: {
            fprintf(stderr, "Usage: %s transform_anno -o <annotation-basename> [options] ANNOTATOR\n\n", prog_name.c_str());

            // fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --aggregate-columns \t\taggregate annotation columns into a bitmask (new column) [off]\n");
            fprintf(stderr, "\t   --anno-label [STR]\t\tname of the aggregated output column [mask]\n");
            fprintf(stderr, "\t   --min-count [INT] \t\texclude k-mers appearing in fewer than this number of columns [1]\n");
            fprintf(stderr, "\t   --min-fraction [FLOAT] \texclude k-mers appearing in fewer than this fraction of columns [0.0]\n");
            fprintf(stderr, "\t   --max-count [INT] \t\texclude k-mers appearing in more than this number of columns [inf]\n");
            fprintf(stderr, "\t   --max-fraction [FLOAT] \texclude k-mers appearing in more than this fraction of columns [1.0]\n");
            fprintf(stderr, "\t   --intersected-anno [STR] \tannotation with columns to intersect with ANNOTATOR (count shared bits) [off]\n");
            fprintf(stderr, "\t   --rename-cols [STR] \tfile with rules for renaming annotation labels []\n");
            fprintf(stderr, "\t                       \texample: 'L_1 L_1_renamed\n");
            fprintf(stderr, "\t                       \t          L_2 L_2_renamed\n");
            fprintf(stderr, "\t                       \t          L_2 L_2_renamed\n");
            fprintf(stderr, "\t                       \t          ... ...........'\n");
            fprintf(stderr, "\t   --anno-type [STR] \ttarget annotation format [column]\n");
            fprintf(stderr, "%s\n", annotation_list);
            fprintf(stderr, "\t   --arity \t\tarity in the brwt tree [2]\n");
            fprintf(stderr, "\t   --greedy \t\tuse greedy column partitioning in brwt construction [off]\n");
            fprintf(stderr, "\t   --linkage \t\tcluster columns and construct linkage matrix [off]\n");
            fprintf(stderr, "\t   --linkage-file [STR]\tlinkage matrix specifying brwt tree structure []\n");
            fprintf(stderr, "\t                       \texample: '0 1 <dist> 4\n");
            fprintf(stderr, "\t                       \t          2 3 <dist> 5\n");
            fprintf(stderr, "\t                       \t          4 5 <dist> 6'\n");
            fprintf(stderr, "\t   --fast \t\ttransform annotation in memory without streaming / sparse subsampling [off]\n");
            fprintf(stderr, "\t   --subsample [INT] \tnumber of rows subsampled for distance estimation in column clustering [1000000]\n");
            fprintf(stderr, "\t   --dump-text-anno \tdump the columns of the annotator as separate text files [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --row-diff-stage [0|1|2] \tstage of the row_diff construction [0]\n");
            fprintf(stderr, "\t   --max-path-length [INT] \tmaximum path length in row_diff annotation [100]\n");
            fprintf(stderr, "\t-i --infile-base [STR] \t\tgraph for generating succ/pred/anchors (for row_diff types) []\n");
            fprintf(stderr, "\t   --count-kmers \t\tadd k-mer counts to the row_diff annotation [off]\n");
            fprintf(stderr, "\t   --coordinates \t\tadd k-mer coordinates to the row_diff annotation [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --parallel-nodes [INT] \tnumber of nodes processed in parallel in brwt tree [n_threads]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --disk-swap [STR] \tdirectory for temporary files [OUT_BASEDIR]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case RELAX_BRWT: {
            fprintf(stderr, "Usage: %s relax_brwt -o <annotation-basename> [options] ANNOTATOR\n\n", prog_name.c_str());

            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --relax-arity [INT] \trelax brwt tree to optimize arity limited to this number [10]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case QUERY: {
            fprintf(stderr, "Usage: %s query -i <GRAPH> -a <ANNOTATION> [options] FILE1 [[FILE2] ...]\n"
                            "\tEach input file is given in FASTA or FASTQ format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for query:\n");
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --fwd-and-reverse \tfor each input sequence, query its reverse complement as well [off]\n");
#endif
            fprintf(stderr, "\t   --align \t\talign sequences instead of mapping k-mers [off]\n");
            fprintf(stderr, "\t   --sparse \t\tuse row-major sparse matrix for row annotation [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --count-labels \t\tcount labels for k-mers from querying sequences [off]\n");
            fprintf(stderr, "\t   --count-kmers \t\tweight k-mers with their annotated counts (requires count annotation) [off]\n");
            fprintf(stderr, "\t   --count-quantiles [FLOAT ...] \tk-mer count quantiles to compute for each label [off]\n"
                            "\t                                 \t\tExample: --count-quantiles '0 0.33 0.5 0.66 1'\n"
                            "\t                                 \t\t(0 corresponds to MIN, 1 corresponds to MAX)\n");
            fprintf(stderr, "\t   --query-coords \t\tquery k-mer coordinates (requires coord annotation) [off]\n");
            fprintf(stderr, "\t   --print-signature \t\tprint vectors indicating present/absent k-mers [off]\n");
            fprintf(stderr, "\t   --num-top-labels \t\tmaximum number of frequent labels to print [off]\n");
            fprintf(stderr, "\t   --discovery-fraction [FLOAT] fraction of labeled k-mers required for annotation [0.7]\n");
            fprintf(stderr, "\t   --labels-delimiter [STR]\tdelimiter for annotation labels [\":\"]\n");
            fprintf(stderr, "\t   --suppress-unlabeled \tdo not show results for sequences missing in graph [off]\n");
            // fprintf(stderr, "\t-d --distance [INT] \tmax allowed alignment distance [0]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
            // fprintf(stderr, "\t   --cache-size [INT] \tnumber of uncompressed rows to store in the cache [0]\n");
            fprintf(stderr, "\t   --fast \t\tquery in batches [off]\n");
            fprintf(stderr, "\t   --batch-size \tquery batch size (number of base pairs) [100000000]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Available options for --align:\n");
            fprintf(stderr, "\t   --align-one-strand \t\t\tdo not align backwards from a seed on basic-mode graphs [off]\n");
            // fprintf(stderr, "\t   --align-alternative-alignments \tthe number of alternative paths to report per seed [1]\n");
            fprintf(stderr, "\t   --align-min-path-score [INT]\t\t\tthe minimum score that a reported path can have [0]\n");
            fprintf(stderr, "\t   --align-max-nodes-per-seq-char [FLOAT]\tmaximum number of nodes to consider per sequence character [12.0]\n");
            fprintf(stderr, "\t   --align-max-ram [FLOAT]\t\t\tmaximum amount of RAM used per alignment in MB [200.0]\n");
            fprintf(stderr, "\t   --align-xdrop [INT]\t\t\t\tthe maximum difference between the current and the best alignment [27]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --batch-align \t\talign against query graph [off]\n");
            fprintf(stderr, "\t   --max-hull-forks [INT]\tmaximum number of forks to take when expanding query graph [4]\n");
            fprintf(stderr, "\t   --max-hull-depth [INT]\tmaximum number of steps to traverse when expanding query graph [max_nodes_per_seq_char * max_seq_len]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Advanced options for scoring:\n");
            fprintf(stderr, "\t   --align-match-score [INT]\t\t\tpositive match score [2]\n");
            fprintf(stderr, "\t   --align-mm-transition-penalty [INT]\t\tpositive transition penalty (DNA only) [3]\n");
            fprintf(stderr, "\t   --align-mm-transversion-penalty [INT]\tpositive transversion penalty (DNA only) [3]\n");
            fprintf(stderr, "\t   --align-gap-open-penalty [INT]\t\tpositive gap opening penalty [5]\n");
            fprintf(stderr, "\t   --align-gap-extension-penalty [INT]\t\tpositive gap extension penalty [2]\n");
            fprintf(stderr, "\t   --align-edit-distance \t\t\tuse unit costs for scoring matrix [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Advanced options for seeding:\n");
            fprintf(stderr, "\t   --align-min-seed-length [INT]\t\tthe minimum length of a seed [graph k]\n");
            fprintf(stderr, "\t   --align-max-seed-length [INT]\t\tthe maximum length of a seed [graph k]\n");
            fprintf(stderr, "\t   --align-min-exact-match [FLOAT] fraction of matching nucleotides required to align sequence [0.0]\n");
            fprintf(stderr, "\t   --align-max-num-seeds-per-locus [INT]\tthe maximum number of allowed inexact seeds per locus [inf]\n");
        } break;
        case SERVER_QUERY: {
            fprintf(stderr, "Usage: %s server_query -i <GRAPH> -a <ANNOTATION> [options]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for server_query:\n");
            fprintf(stderr, "\t   --port [INT] \tTCP port for incoming connections [5555]\n");
            fprintf(stderr, "\t   --address \t\tinterface for incoming connections (default: all)\n");
            fprintf(stderr, "\t   --sparse \t\tuse the row-major sparse matrix to annotate graph [off]\n");
            // fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file []\n");
            // fprintf(stderr, "\t-d --distance [INT] \tmax allowed alignment distance [0]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tmaximum number of parallel connections [1]\n");
            // fprintf(stderr, "\t   --cache-size [INT] \tnumber of uncompressed rows to store in the cache [0]\n");
        } break;
    }

    fprintf(stderr, "\n\tGeneral options:\n");
    fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
    fprintf(stderr, "\t-h --help \t\tprint usage info\n");
    fprintf(stderr, "\n");
}

} // namespace cli
} // namespace mtg
