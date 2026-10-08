#include "config.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <optional>
#include <unordered_set>
#include <filesystem>

#include "common/threads/threading.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/logger.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "seq_io/formats.hpp"
#include "kmer/kmer_extractor.hpp"


namespace mtg {
namespace cli {

using mtg::graph::boss::BOSS;
using mtg::graph::DeBruijnGraph;


const size_t Config::kDefaultIndexSuffixLen
    = 24 / std::log2(kmer::KmerExtractor2Bit().alphabet.size());

void print_welcome_message() {
    fprintf(stderr, "#############################\n");
    fprintf(stderr, "### Welcome to MetaGraph! ###\n");
    fprintf(stderr, "#############################\n\n");
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
    } else if (!strcmp(argv[1], "merge_anno")) {
        identity = MERGE_ANNOTATIONS;
    } else if (!strcmp(argv[1], "query")) {
        identity = QUERY;
    } else if (!strcmp(argv[1], "traverse")) {
        identity = TRAVERSE;
    } else if (!strcmp(argv[1], "pattern")) {
        identity = PATTERN;
    } else if (!strcmp(argv[1], "server_query")) {
        identity = SERVER_QUERY;
        num_top_labels = 10'000;
        memory_available = 0;
    } else if (!strcmp(argv[1], "transform")) {
        identity = TRANSFORM;
    } else if (!strcmp(argv[1], "transform_anno")) {
        identity = TRANSFORM_ANNOTATION;
        tmp_dir = "OUTFBASE_TEMP_DIR";
        memory_available = 1000; // 1 TB
    } else if (!strcmp(argv[1], "assemble")) {
        identity = ASSEMBLE;
    } else if (!strcmp(argv[1], "relax_brwt")) {
        identity = RELAX_BRWT;
    } else if (!strcmp(argv[1], "--version")) {
        std::cout << "Version: " VERSION << std::endl;
        exit(0);
    } else if (!strcmp(argv[1], "--advanced")) {
        advanced = true;
        print_welcome_message();
        print_usage(argv[0]);
        exit(0);
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
    bool xdrop_override = false;
    // An integer (ms) the traversal capabilities state: non-negative and exactly representable
    // as a double, so that a JSON client reads it as written and a ledger can add to it (atoll
    // wrapped -1 to 2^64 - 1, stated as such; review of pass 5). The attempts' allowance too: a
    // fraction would be stated rounded and enforced unrounded
    // the value of |text| when it is an integer in [0, 2^53 - 1] (none otherwise)
    const auto exact_integer = [](const char *text) -> std::optional<uint64_t> {
        char *end = nullptr;
        const double v = std::strtod(text, &end);
        if (end == text || *end != '\0' || !(v >= 0) || v != std::floor(v)
                || v > 9007199254740991.0) {
            return std::nullopt;
        }
        return static_cast<uint64_t>(v);
    };
    const auto exact_ms = [&](const char *option, const char *text, uint64_t *out) {
        const std::optional<uint64_t> v = exact_integer(text);
        if (!v) {
            std::cerr << "Error: " << option << " must be an integer in [0, 2^53 - 1], got '"
                      << text << "'" << std::endl;
            print_usage_and_exit = true;
            return;
        }
        *out = *v;
    };
    // An integer in [0, |max|]: the attempts' retention settings (review of pass 5, finding 4:
    // atoll read -1 as 2^64 - 1 seconds, which the capabilities stated while the conversion to
    // a signed std::chrono::seconds expired every tombstone at once). |max| keeps every
    // conversion the registry makes exact (seconds to steady-clock nanoseconds, to Unix-epoch
    // milliseconds) far from overflowing, so a value is refused at start-up rather than wrapped.
    // Every refusal names the option's own range (a negative or non-numeric value too)
    const auto bounded = [&](const char *option, const char *text, uint64_t max, uint64_t *out) {
        const std::optional<uint64_t> v = exact_integer(text);
        if (!v || *v > max) {
            std::cerr << "Error: " << option << " must be an integer in [0, " << max
                      << "], got '" << text << "'" << std::endl;
            print_usage_and_exit = true;
            return;
        }
        *out = *v;
    };

    // parse remaining command line items
    for (int i = 2; i < argc; ++i) {
        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
            common::set_verbose(true);
        } else if (!strcmp(argv[i], "--mmap")) {
            utils::set_mmap(true);
        } else if (!strcmp(argv[i], "--madv-random")) {
            utils::set_madvise(true);
        } else if (!strcmp(argv[i], "--one-pass-brwt")) {
            annot::matrix::set_one_pass_brwt(true);
        } else if (!strcmp(argv[i], "--print")) {
            print_graph = true;
        } else if (!strcmp(argv[i], "--advanced")) {
            advanced = true;
            if (argc == 3)
                print_usage_and_exit = true;
        } else if (!strcmp(argv[i], "--print-col-names")) {
            print_column_names = true;
        } else if (!strcmp(argv[i], "--print-internal")) {
            print_graph_internal_repr = true;
        } else if (!strcmp(argv[i], "--print-counts-hist")) {
            print_counts_hist = true;
        } else if (!strcmp(argv[i], "--coordinates")) {
            coordinates = true;
        } else if (!strcmp(argv[i], "--index-header-coords")) {
            index_header_coords = true;
        } else if (!strcmp(argv[i], "--no-coord-mapping")) {
            no_coord_mapping = true;
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
        } else if (!strcmp(argv[i], "--query-mode")) {
            query_mode = string_to_querymode(get_value(i++));
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
        } else if (!strcmp(argv[i], "--batch-size")) {
            query_batch_size = atoll(get_value(i++));
        } else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--parallel")) {
            set_num_threads(atoi(get_value(i++)));
        } else if (!strcmp(argv[i], "--parallel-nodes")) {
            parallel_nodes = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--threads-each")) {
            parallel_each = atoi(get_value(i++));
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
        } else if (!strcmp(argv[i], "--compute-overlap")) {
            intersected_columns = get_value(i++);
        } else if (!strcmp(argv[i], "--min-fraction")) {
            min_fraction = std::stod(get_value(i++));
        } else if (!strcmp(argv[i], "--max-fraction")) {
            max_fraction = std::stod(get_value(i++));
        } else if (!strcmp(argv[i], "--min-value")) {
            min_value = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--max-value")) {
            max_value = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--mem-cap-gb")) {
            memory_available = atof(get_value(i++));
        } else if (!strcmp(argv[i], "--dump-text-anno")) {
            dump_text_anno = true;
        } else if (!strcmp(argv[i], "--min-kmers-fraction-label")) {
            discovery_fraction = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--align-rel-score-cutoff")) {
            alignment_rel_score_cutoff = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--min-kmers-fraction-graph")) {
            presence_fraction = std::stof(get_value(i++));
        } else if (!strcmp(argv[i], "--query-presence")) {
            query_presence = true;
        } else if (!strcmp(argv[i], "--verbose-output")) {
            verbose_output = true;
        } else if (!strcmp(argv[i], "--filter-present")) {
            filter_present = true;
        } else if (!strcmp(argv[i], "--map")) {
            map_sequences = true;
        } else if (!strcmp(argv[i], "--align")) {
            align_sequences = true;
        } else if (!strcmp(argv[i], "--align-only-forwards")) {
            align_only_forwards = true;
        } else if (!strcmp(argv[i], "--align-output-path")) {
            align_output_path = true;
        } else if (!strcmp(argv[i], "--align-edit-distance")) {
            alignment_edit_distance = true;
        } else if (!strcmp(argv[i], "--align-chain")) {
            alignment_chain = true;
        } else if (!strcmp(argv[i], "--align-post-chain")) {
            alignment_post_chain = true;
        } else if (!strcmp(argv[i], "--align-no-seed-complexity-filter")) {
            alignment_seed_complexity_filter = false;
        } else if (!strcmp(argv[i], "--num-chars")) {
            num_chars = atoll(get_value(i++));
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
        } else if (!strcmp(argv[i], "--align-end-bonus")) {
            alignment_end_bonus = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-alternative-alignments")) {
            alignment_num_alternative_paths = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-min-path-score")) {
            alignment_min_path_score = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--align-xdrop")) {
            alignment_xdrop = atol(get_value(i++));
            xdrop_override = true;
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
        } else if (!strcmp(argv[i], "--num-top-labels")) {
            num_top_labels = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--port")) {
            port = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--address")) {
            host_address = get_value(i++);
        }else if (!strcmp(argv[i], "--suffix")) {
            suffix = get_value(i++);
        } else if (!strcmp(argv[i], "--diff-assembly-rules")) {
            assembly_config_file = get_value(i++);
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
        } else if (!strcmp(argv[i], "--resolve")) {
            traverse_resolve = true;
        } else if (!strcmp(argv[i], "--index-inventory")) {
            traverse_index_inventory = true;
        } else if (!strcmp(argv[i], "--index-release")) {
            index_release = get_value(i++);
        } else if (!strcmp(argv[i], "--index-name")) {
            index_name = get_value(i++);
        } else if (!strcmp(argv[i], "--index-manifest")) {
            index_manifest = get_value(i++);
        } else if (!strcmp(argv[i], "--traverse-max-time-ms")) {
            traverse_max_time_ms = atof(get_value(i++));
        } else if (!strcmp(argv[i], "--traverse-max-seeds")) {
            traverse_max_seeds = atoll(get_value(i++));
        } else if (!strcmp(argv[i], "--traverse-max-seed-bp")) {
            traverse_max_seed_bp = atoll(get_value(i++));
        } else if (!strcmp(argv[i], "--traverse-max-seed-labels")) {
            traverse_max_seed_labels = atoll(get_value(i++));
        } else if (!strcmp(argv[i], "--resolve-max-query-bp")) {
            resolve_max_query_bp = atoll(get_value(i++));
        } else if (!strcmp(argv[i], "--traverse-max-memory-mb")) {
            exact_ms(argv[i], get_value(i), &traverse_max_memory_mb);
            i++;
        } else if (!strcmp(argv[i], "--traverse-max-work-units")) {
            exact_ms(argv[i], get_value(i), &traverse_max_work_units);
            i++;
        } else if (!strcmp(argv[i], "--traverse-attempt-allowance-ms")) {
            uint64_t ms = 0;
            exact_ms(argv[i], get_value(i), &ms);
            traverse_attempt_allowance_ms = static_cast<double>(ms);
            i++;
        } else if (!strcmp(argv[i], "--traverse-clock-skew-ms")) {
            exact_ms(argv[i], get_value(i), &traverse_clock_skew_ms);
            i++;
        } else if (!strcmp(argv[i], "--traverse-chunk-target-ms")) {
            exact_ms(argv[i], get_value(i), &traverse_chunk_target_ms);
            i++;
        } else if (!strcmp(argv[i], "--traverse-path-cache-mb")) {
            exact_ms(argv[i], get_value(i), &traverse_path_cache_mb);
            i++;
        } else if (!strcmp(argv[i], "--traverse-compression-level")) {
            traverse_compression_level = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--traverse-delivery-compress-mbps")) {
            traverse_delivery_compress_mbps = atof(get_value(i++));
        } else if (!strcmp(argv[i], "--traverse-delivery-build-mbps")) {
            traverse_delivery_build_mbps = atof(get_value(i++));
        } else if (!strcmp(argv[i], "--traverse-attempt-retention-s")) {
            bounded(argv[i], get_value(i), kMaxAttemptRetentionS, &traverse_attempt_retention_s);
            i++;
        } else if (!strcmp(argv[i], "--traverse-attempt-retention")) {
            uint64_t count = 0;
            bounded(argv[i], get_value(i), kMaxAttemptRetentionCount, &count);
            traverse_attempt_retention = count;
            i++;
        } else if (!strcmp(argv[i], "--traverse-attempt-tombstone-max-s")) {
            bounded(argv[i], get_value(i), kMaxAttemptRetentionS,
                    &traverse_attempt_tombstone_max_s);
            i++;
        } else if (!strcmp(argv[i], "--pattern-min-information-bits")) {
            // a number of bits: stated in the capabilities as given
            char *end = nullptr;
            const char *text = get_value(i++);
            pattern_min_information_bits = std::strtod(text, &end);
            if (end == text || *end != '\0' || !(pattern_min_information_bits >= 0)
                    || !std::isfinite(pattern_min_information_bits)) {
                std::cerr << "Error: --pattern-min-information-bits must be a number >= 0, got '"
                          << text << "'" << std::endl;
                print_usage_and_exit = true;
            }
        } else if (!strcmp(argv[i], "--pattern-max-contexts")) {
            exact_ms(argv[i], get_value(i), &pattern_max_contexts);
            i++;
        } else if (!strcmp(argv[i], "--pattern-max-anchors")) {
            exact_ms(argv[i], get_value(i), &pattern_max_anchors);
            i++;
        } else if (!strcmp(argv[i], "--pattern-max-paths")) {
            exact_ms(argv[i], get_value(i), &pattern_max_paths);
            i++;
        } else if (!strcmp(argv[i], "--pattern-max-steps")) {
            exact_ms(argv[i], get_value(i), &pattern_max_steps);
            i++;
        } else if (!strcmp(argv[i], "--pattern-max-time-ms")) {
            exact_ms(argv[i], get_value(i), &pattern_max_time_ms);
            i++;
        } else if (!strcmp(argv[i], "--pattern-default-time-ms")) {
            exact_ms(argv[i], get_value(i), &pattern_default_time_ms);
            i++;
        } else if (!strcmp(argv[i], "--pattern-finalize-ms")) {
            exact_ms(argv[i], get_value(i), &pattern_finalize_ms);
            i++;
        } else if (!strcmp(argv[i], "--pattern-delivery-build-mbps")) {
            pattern_delivery_build_mbps = atof(get_value(i++));
        } else if (!strcmp(argv[i], "--pattern-delivery-compress-mbps")) {
            pattern_delivery_compress_mbps = atof(get_value(i++));
        } else if (!strcmp(argv[i], "--pattern-max-patterns")) {
            exact_ms(argv[i], get_value(i), &pattern_max_patterns);
            i++;
        } else if (!strcmp(argv[i], "--pattern-max-labels-per-anchor")) {
            exact_ms(argv[i], get_value(i), &pattern_max_labels_per_anchor);
            i++;
        } else if (!strcmp(argv[i], "--pattern-max-annotation-work")) {
            exact_ms(argv[i], get_value(i), &pattern_max_annotation_work);
            i++;
        } else if (!strcmp(argv[i], "--pattern-max-memory-mb")) {
            exact_ms(argv[i], get_value(i), &pattern_max_memory_mb);
            i++;
        } else if (!strcmp(argv[i], "--pattern-max-labels")) {
            exact_ms(argv[i], get_value(i), &pattern_max_labels);
            i++;
        } else if (!strcmp(argv[i], "--pattern-max-occurrences")) {
            exact_ms(argv[i], get_value(i), &pattern_max_occurrences);
            i++;
        } else if (!strcmp(argv[i], "--pattern-build-mask")) {
            pattern_build_mask = true;
        } else if (!strcmp(argv[i], "--force")) {
            force = true;
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
        } else if (!strcmp(argv[i], "--cleaning-threshold-percentile")) {
            cleaning_threshold_percentile = std::stod(get_value(i++));
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
        } else if (!strcmp(argv[i], "--in-ram")) {
            in_ram = true;
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
        } else if (!strcmp(argv[i], "--subsample-rows")) {
            subsample_rows = true;
        } else if (!strcmp(argv[i], "--linkage-file")) {
            linkage_file = get_value(i++);
        } else if (!strcmp(argv[i], "--arity")) {
            arity_brwt = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--relax-arity")) {
            relax_arity_brwt = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "--RA-ivbuff-size")) {
            RA_ivbuffer_size = atoll(get_value(i++));
        // } else if (!strcmp(argv[i], "--cache-size")) {
        //     row_cache_size = atoi(get_value(i++));
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_welcome_message();
            print_usage(argv[0], identity);
            exit(0);
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

    // transform --mask-dummy (DESIGN-pattern-search.md §4) writes the dummy-edge mask beside
    // the graph and nothing else, so another transformation asked with it would be dropped
    // unseen: refused instead (checked before --to-fasta turns the identity into clean)
    if (identity == TRANSFORM && mark_dummy_kmers) {
        if (clear_dummy || initialize_bloom || to_adj_list || to_fasta
                || graph_mode != DeBruijnGraph::BASIC || state != BOSS::State::STAT
                || node_suffix_length != kDefaultIndexSuffixLen) {
            std::cerr << "Error: --mask-dummy writes only the dummy-edge mask "
                         "(<graph>.edgemask) and leaves the graph as it is; it cannot be "
                         "combined with --clear-dummy, --initialize-bloom, --to-adj-list, "
                         "--to-fasta, --unitigs, --mode, --state or --index-ranges" << std::endl;
            print_usage_and_exit = true;
        }
        // the loader reads the mask beside the graph it loads, under the graph's own name:
        // a mask written under another name would never be read
        if (outfbase.size() && fnames.size() == 1
                && utils::remove_suffix(outfbase, ".dbg")
                        != utils::remove_suffix(fnames[0], ".dbg")) {
            std::cerr << "Error: --mask-dummy writes the mask beside the graph, where the "
                         "loader reads it (" << utils::remove_suffix(fnames[0], ".dbg")
                      << ".edgemask): omit -o, or name the graph itself" << std::endl;
            print_usage_and_exit = true;
        }
    }
    if (force && !(identity == TRANSFORM && mark_dummy_kmers)) {
        std::cerr << "Error: --force applies only to transform --mask-dummy" << std::endl;
        print_usage_and_exit = true;
    }
    // the mask is built where the graph is loaded at start-up: one graph (-i / -a) of
    // server_query or pattern. A graph list's graphs are loaded by requests, and /pattern is
    // not served on them yet (a later increment), so the flag would do nothing there
    if (pattern_build_mask
            && !(identity == PATTERN || (identity == SERVER_QUERY && fnames.empty()))) {
        std::cerr << "Error: --pattern-build-mask applies to server_query with one graph "
                     "(-i / -a) and to pattern" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == TRANSFORM && to_fasta)
        identity = CLEAN;

    if (!xdrop_override && alignment_chain)
        alignment_xdrop = 100;

    // given kmc_pre and kmc_suf pair, only include one
    // this still allows for the same file to be included multiple times
    std::unordered_set<std::string> kmc_file_set;

    for (auto it = fnames.begin(); it != fnames.end(); ++it) {
        if (seq_io::file_format(*it) == "KMC"
                && !kmc_file_set.insert(utils::remove_suffix(*it, ".kmc_pre", ".kmc_suf")).second)
            fnames.erase(it--);
    }

    if (!print_usage_and_exit && !fnames.size()
                      && identity != STATS
                      && identity != SERVER_QUERY
                      && !(identity == BUILD && complete)
                      && !(identity == CONCATENATE && !infbase.empty())
                      // the inventory names its files with -i/-a; reading stdin would wait for EOF
                      && !(identity == TRAVERSE && traverse_index_inventory)) {
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
            && !(identity == TRAVERSE && traverse_index_inventory)
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

    if (identity == QUERY && infbase.empty())
        print_usage_and_exit = true;

    if ((identity == QUERY || identity == SERVER_QUERY || identity == ALIGN)
            && alignment_num_alternative_paths == 0) {
        std::cerr << "Error: align-alternative-alignments must be > 0" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == ANNOTATE && infbase.empty())
        print_usage_and_exit = true;

    if ((identity == ANNOTATE || identity == EXTEND) && infbase_annotators.size() > 1) {
        std::cerr << "Error: one annotator at most is allowed for extension." << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == ANNOTATE
            && !filename_anno && !annotate_sequence_headers && !anno_labels.size()) {
        std::cerr << "Error: no annotation labels passed (see flags --anno-filename --anno-header --anno-label)" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == ASSEMBLE
            && (infbase_annotators.size() && assembly_config_file.empty())) {
        std::cerr << "Error: annotator passed, but no differential assembly rule config file provided" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == EXTEND && infbase.empty())
        print_usage_and_exit = true;

    if (identity == QUERY && infbase_annotators.size() != 1)
        print_usage_and_exit = true;

    if (identity == TRAVERSE && (infbase.empty() || infbase_annotators.size() != 1)) {
        std::cerr << "Error: traverse requires a graph (-i) and exactly one annotation (-a)" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == PATTERN && (infbase.empty() || infbase_annotators.size() != 1)) {
        std::cerr << "Error: pattern requires a graph (-i) and exactly one annotation (-a)" << std::endl;
        print_usage_and_exit = true;
    }

    // the /pattern caps (DESIGN-pattern-search.md §5.3): a request needs at least one step, a
    // pattern, and time beyond the finalisation reserve, which lies inside the time budget;
    // the default budget is one a request may name
    if ((identity == PATTERN || identity == SERVER_QUERY)
            && (pattern_max_steps < 1 || pattern_max_patterns < 1
                || pattern_default_time_ms <= pattern_finalize_ms
                || pattern_default_time_ms > pattern_max_time_ms)) {
        std::cerr << "Error: --pattern-max-steps and --pattern-max-patterns must be at least 1, "
                     "and --pattern-default-time-ms above --pattern-finalize-ms and at most "
                     "--pattern-max-time-ms" << std::endl;
        print_usage_and_exit = true;
    }
    // a /pattern deadline the transport cannot honour would be accepted and echoed, and the
    // connection closed at the content timeout before any answer or 503 (review of
    // 2026-10-07, R1-05): the cap stays under it (the CLI has no transport)
    if (identity == SERVER_QUERY && pattern_max_time_ms > kServerMaxDeadlineMs) {
        std::cerr << "Error: --pattern-max-time-ms must be at most " << kServerMaxDeadlineMs
                  << " on server_query (the " << kServerContentTimeoutS << " s content timeout "
                  "less " << kServerTransportMarginMs << " ms for the transport)" << std::endl;
        print_usage_and_exit = true;
    }
    // the rates the finalisation estimate of an answer assumes (review of 2026-10-07,
    // X-EFFICIENCY-04): a rate that is not a positive number would make the estimate meaningless
    if ((identity == PATTERN || identity == SERVER_QUERY)
            && (!(pattern_delivery_build_mbps > 0) || !std::isfinite(pattern_delivery_build_mbps)
                || !(pattern_delivery_compress_mbps > 0)
                || !std::isfinite(pattern_delivery_compress_mbps))) {
        std::cerr << "Error: --pattern-delivery-build-mbps and --pattern-delivery-compress-mbps "
                     "must be positive numbers" << std::endl;
        print_usage_and_exit = true;
    }
    // the labelled retrieval's caps (increment 3): a row keeps at least one label, a read
    // needs work and memory; the memory account in bytes must fit 64 bits (2^40 MiB)
    if ((identity == PATTERN || identity == SERVER_QUERY)
            && (pattern_max_labels_per_anchor < 1 || pattern_max_annotation_work < 1
                || pattern_max_memory_mb < 1 || pattern_max_memory_mb > (uint64_t(1) << 40))) {
        std::cerr << "Error: --pattern-max-labels-per-anchor, --pattern-max-annotation-work and "
                     "--pattern-max-memory-mb must be at least 1, and --pattern-max-memory-mb "
                     "at most 2^40" << std::endl;
        print_usage_and_exit = true;
    }

    if ((identity == TRAVERSE || identity == SERVER_QUERY || identity == PATTERN)
            && (!index_name.empty() || !index_manifest.empty())) {
        if (fnames.size() && identity == SERVER_QUERY) {
            // one name and one manifest describe one index, not a list of them
            std::cerr << "Error: --index-name and --index-manifest describe a single index "
                         "(-i / -a), not a graph list" << std::endl;
            print_usage_and_exit = true;
        }
        for (char c : index_name) {
            if (!std::isalnum(static_cast<unsigned char>(c)) && c != '.' && c != '_' && c != '-') {
                std::cerr << "Error: --index-name must match [A-Za-z0-9._-]+" << std::endl;
                print_usage_and_exit = true;
                break;
            }
        }
    }

    // a request's bounds.max_memory_mb is at most 1,048,576 (1 TiB): a larger maximum could
    // not be echoed as a budget a request may give
    if (traverse_path_cache_mb > 1'048'576) {
        std::cerr << "Error: --traverse-path-cache-mb must be in [0, 1048576]" << std::endl;
        print_usage_and_exit = true;
    }
    if (traverse_max_memory_mb > 1'048'576) {
        std::cerr << "Error: --traverse-max-memory-mb must be in [0, 1048576]" << std::endl;
        print_usage_and_exit = true;
    }
    if (traverse_compression_level < 1 || traverse_compression_level > 9) {
        std::cerr << "Error: --traverse-compression-level must be in [1, 9]" << std::endl;
        print_usage_and_exit = true;
    }
    if (!(traverse_delivery_compress_mbps > 0) || !(traverse_delivery_build_mbps > 0)) {
        std::cerr << "Error: the delivery rates (--traverse-delivery-compress-mbps, "
                     "--traverse-delivery-build-mbps) must be positive" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == SERVER_QUERY
            && (fnames.size() > 1
                || (fnames.size() && (infbase.size() || infbase_annotators.size()))
                || (fnames.empty() && (infbase.empty() || infbase_annotators.size() != 1))))
        print_usage_and_exit = true;  // only one of fnames or (infbase & annotator) must be used

    if ((identity == TRANSFORM
            || identity == CLEAN
            || identity == ASSEMBLE
            || identity == RELAX_BRWT)
                    && fnames.size() != 1) {
        std::cerr << "Error: exactly one graph must be provided for this mode" << std::endl;
        print_usage_and_exit = true;
    }

    // (transform --mask-dummy writes beside its input graph and needs no -o)
    if (((identity == TRANSFORM && !mark_dummy_kmers)
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
                                    || anno_type == RowDiffDisk
                                    || anno_type == IntRowDiffBRWT
                                    || anno_type == IntRowDiffDisk
                                    || anno_type == RowDiffRowFlat
                                    || anno_type == RowDiffRowSparse
                                    || anno_type == RowDiffDiskCoord
                                    || anno_type == RowDiffBRWTCoord
                                    || anno_type == RowDiffCoord;
        // For row_diff_brwt and row_diff_flat inputs, the anchor and fork-succ
        // bitmaps are embedded in the .annodbg, so no graph is needed when
        // re-formatting to another row_diff representation.
        const bool input_has_embedded_anchor = fnames.size()
                && (utils::ends_with(fnames.front(), annot::RowDiffBRWTAnnotator::kExtension)
                        || utils::ends_with(fnames.front(), annot::RowDiffRowFlatAnnotator::kExtension));
        if (input_has_embedded_anchor && fnames.size() != 1) {
            std::cerr << "Can only convert " << fnames.front()
                      << " annotations one at a time" << std::endl;
            print_usage_and_exit = true;
        }
        if (to_row_diff && !infbase.size() && !input_has_embedded_anchor) {
            std::cerr << "Path to graph must be passed with '-i <GRAPH>'" << std::endl;
            print_usage_and_exit = true;
        } else if (!to_row_diff && infbase.size()) {
            std::cerr << "Graph is only required for transform to row_diff types"
                      << std::endl;
            print_usage_and_exit = true;
        } else if (input_has_embedded_anchor && infbase.size()) {
            common::logger->warn(
                    "Graph is not needed when the input annotation has embedded "
                    "anchors (row_diff_brwt, row_diff_flat); ignoring '-i {}'",
                    infbase);
            infbase.clear();
        }
    }

    if (identity == MERGE && fnames.size() < 2)
        print_usage_and_exit = true;

    if (identity == COMPARE && fnames.size() != 2)
        print_usage_and_exit = true;

    if (discovery_fraction < 0 || discovery_fraction > 1)
        print_usage_and_exit = true;

    if (presence_fraction < 0 || presence_fraction > 1)
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
        case RowDiffRowFlat:
            return "row_diff_flat";
        case RowDiffRowSparse:
            return "row_diff_sparse";
        case RowSparse:
            return "row_sparse";
        case RowDiffDisk:
            return "row_diff_disk";
        case IntRowDiffDisk:
            return "row_diff_int_disk";
        case RowDiffDiskCoord:
            return "row_diff_disk_coord";
        case IntBRWT:
            return "int_brwt";
        case IntRowDiffBRWT:
            return "row_diff_int_brwt";
        case ColumnCoord:
            return "column_coord";
        case BRWTCoord:
            return "brwt_coord";
        case RowDiffCoord:
            return "row_diff_coord";
        case RowDiffBRWTCoord:
            return "row_diff_brwt_coord";
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
    } else if (string == "row_diff_flat") {
        return AnnotationType::RowDiffRowFlat;
    } else if (string == "row_diff_sparse") {
        return AnnotationType::RowDiffRowSparse;
    } else if (string == "row_sparse") {
        return AnnotationType::RowSparse;
    } else if (string == "row_diff_disk") {
        return AnnotationType::RowDiffDisk;
    } else if (string == "row_diff_int_disk") {
        return AnnotationType::IntRowDiffDisk;
    } else if (string == "row_diff_disk_coord") {
        return AnnotationType::RowDiffDiskCoord;
    } else if (string == "int_brwt") {
        return AnnotationType::IntBRWT;
    } else if (string == "row_diff_int_brwt") {
        return AnnotationType::IntRowDiffBRWT;
    } else if (string == "column_coord") {
        return AnnotationType::ColumnCoord;
    } else if (string == "brwt_coord") {
        return AnnotationType::BRWTCoord;
    } else if (string == "row_diff_coord") {
        return AnnotationType::RowDiffCoord;
    } else if (string == "row_diff_brwt_coord") {
        return AnnotationType::RowDiffBRWTCoord;
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

    } else if (string == "sshash") {
        return GraphType::SSHASH;
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

std::string Config::querymode_to_string(QueryMode mode) {
    switch (mode) {
        case QueryMode::LABELS:
                return "labels";
        case QueryMode::MATCHES:
                return "matches";
        case QueryMode::COUNTS_SUM:
                return "counts-sum";
        case QueryMode::COUNTS:
                return "counts";
        case QueryMode::COORDS:
                return "coords";
        case QueryMode::SIGNATURE:
                return "signature";
    }
    throw std::runtime_error("Never happens");
}

QueryMode Config::string_to_querymode(const std::string &string) {
    if (string == "labels") {
        return QueryMode::LABELS;
    } else if (string == "matches") {
        return QueryMode::MATCHES;
    } else if (string == "counts-sum") {
        return QueryMode::COUNTS_SUM;
    } else if (string == "counts") {
        return QueryMode::COUNTS;
    } else if (string == "coords") {
        return QueryMode::COORDS;
    } else if (string == "signature") {
        return QueryMode::SIGNATURE;
    } else {
        std::cerr << "Error: unknown query mode. Check value passed with flag '--query-mode'." << std::endl;
        exit(1);
    }
}


void Config::print_usage(const std::string &prog_name, IdentityType identity) {
    const char annotation_list[] = "\t\t( column, brwt, rb_brwt, int_brwt,\n"
                                   "\t\t  column_coord, brwt_coord, row_diff_coord, row_diff_brwt_coord,\n"
                                   "\t\t  row_diff, row_diff_brwt, row_diff_flat, row_diff_sparse, row_diff_int_brwt,\n"
                                   "\t\t  row_diff_disk, row_diff_int_disk, row_diff_disk_coord,\n"
                                   "\t\t  row, flat, row_sparse, rbfish, bin_rel_wt )";

    switch (identity) {
        case NO_IDENTITY: {
            fprintf(stderr, "Usage: %s <command> [command specific options]\n\n", prog_name.c_str());

            fprintf(stderr, "Available commands:\n");

            fprintf(stderr, "\tbuild\t\tconstruct a graph object from input sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats into a given graph\n\n");

            fprintf(stderr, "\tclean\t\tclean an existing graph and extract sequences from it\n");
            fprintf(stderr, "\t\t\tin fast[a|q] formats\n\n");

            fprintf(stderr, "\ttransform\tgiven a graph, transform it to other formats\n\n");

if (advanced) {
            fprintf(stderr, "\textend\t\textend an existing graph with new sequences from\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats\n\n");

            fprintf(stderr, "\tmerge\t\tintegrate a given set of graph structures\n");
            fprintf(stderr, "\t\t\tand output a new graph structure\n\n");

            fprintf(stderr, "\tconcatenate\tcombine the results of the external merge or\n");
            fprintf(stderr, "\t\t\tconstruction and output the resulting graph structure\n\n");

            fprintf(stderr, "\tcompare\t\tcheck whether two given graphs are identical\n\n");
}
            fprintf(stderr, "\talign\t\talign sequences provided in fast[a|q] files to graph\n\n");

            fprintf(stderr, "\tannotate\tgiven a graph and a fast[a|q] file, annotate\n");
            fprintf(stderr, "\t\t\tthe respective kmers\n\n");
if (advanced) {
            fprintf(stderr, "\tmerge_anno\tmerge annotations\n\n");
}
            fprintf(stderr, "\trelax_brwt\toptimize the tree structure in brwt annotator\n\n");

            fprintf(stderr, "\ttransform_anno\tchange representation of the graph annotation\n\n");

            fprintf(stderr, "\tassemble\tgiven a graph, extract sequences from it\n\n");

            fprintf(stderr, "\tquery\t\tannotate sequences from fast[a|q] files\n\n");
            fprintf(stderr, "\tserver_query\tannotate received sequences and send annotations back\n\n");

            fprintf(stderr, "\ttraverse\textend sequence seeds through the graph along consistent\n");
            fprintf(stderr, "\t\t\tannotation labels (JSON request files)\n\n");

            fprintf(stderr, "\tpattern\t\tcount and extract the graph contexts of short DNA or IUPAC\n");
            fprintf(stderr, "\t\t\tpatterns, and read their labels (JSON request files)\n\n");

            fprintf(stderr, "\tstats\t\tprint graph statistics for given graph(s) or annotation\n\n");

            fprintf(stderr, "General options:\n");
            fprintf(stderr, "\t--advanced \tshow other advanced and legacy options [off]\n");
            fprintf(stderr, "\t--version \tprint version\n");
            fprintf(stderr, "\n");
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
            fprintf(stderr, "\t   --graph [STR] \tgraph representation: succinct / bitmap / hash / hashstr / hashfast [succinct] / sshash\n");
            fprintf(stderr, "\t   --state [STR] \tstate of succinct graph: small / dynamic / stat / fast [stat]\n");
            fprintf(stderr, "\t   --in-ram \t\tconstruct succinct graph in RAM instead of inplace [off]\n");
            fprintf(stderr, "\t   --count-kmers \tcount k-mers and build weighted graph [off]\n");
            fprintf(stderr, "\t   --count-width \tnumber of bits used to represent k-mer abundance [8]\n");
            fprintf(stderr, "\t   --index-ranges [INT]\tindex all node ranges in BOSS for suffixes of given length [%zu]\n", kDefaultIndexSuffixLen);
            fprintf(stderr, "\t   --num-chars [INT]\tif the number of characters is known beforehand, enter it here [0]\n");
            fprintf(stderr, "\t-k --kmer-length [INT] \tlength of the k-mer to use [3]\n");
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --mode \t\tk-mer indexing mode: basic / canonical / primary [basic]\n");
#endif
            fprintf(stderr, "\t   --complete \t\tconstruct a complete graph (only for Bitmap graph) [off]\n");
            fprintf(stderr, "\t   --mem-cap-gb [INT] \tpreallocated buffer size in GB [1]\n");
if (advanced) {
            fprintf(stderr, "\t   --dynamic \t\tuse dynamic build method [off]\n");
            fprintf(stderr, "\t-l --len-suffix [INT] \tk-mer suffix length for building graph from chunks [0]\n");
            fprintf(stderr, "\t   --suffix \t\tbuild graph chunk only for k-mers with the suffix given [off]\n");
}
            fprintf(stderr, "\t-o --outfile-base [STR]\tbasename of output file []\n");
if (advanced) {
            fprintf(stderr, "\t   --mask-dummy \tbuild mask for dummy k-mers (only for Succinct graph; requires --in-ram) [off]\n");
}
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
            fprintf(stderr, "\t   --disk-swap [STR] \tdirectory to use for temporary files [off]\n");
if (advanced) {
            fprintf(stderr, "\t   --disk-cap-gb [INT] \tmax temp disk space to use before forcing a merge, in GB [inf]\n");
}
        } break;
        case CLEAN: {
            fprintf(stderr, "Usage: %s clean -o <outfile-base> [options] GRAPH\n\n", prog_name.c_str());
            fprintf(stderr, "Available options for clean:\n");
            fprintf(stderr, "\t   --min-count [INT] \t\tmin k-mer abundance, including [1]\n");
            fprintf(stderr, "\t   --max-count [INT] \t\tmax k-mer abundance, excluding [inf]\n");
if (advanced) {
            fprintf(stderr, "\t   --num-singletons [INT] \treset the number of count 1 k-mers in histogram (0: off) [0]\n");
}
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --prune-tips [INT] \t\tprune all dead ends shorter than this value [1]\n");
            fprintf(stderr, "\t   --prune-unitigs [INT] \tprune all unitigs with median k-mer counts smaller\n"
                            "\t                         \t\tthan this value (0: auto) [1]\n");
            fprintf(stderr, "\t   --cleaning-threshold-percentile [FLOAT] the percentile of the k-mer count distribution to set as the cleaning threshold [0.001]\n");
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
if (advanced) {
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --fwd-and-reverse \t\tfor each input sequence, report a separate alignment for its reverse complement as well [off]\n");
#endif
            fprintf(stderr, "\t   --header-comment-delim [STR]\tdelimiter for joining fasta header with comment [off]\n");
}
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
            fprintf(stderr, "\t   --batch-size [INT] \t\tquery batch size (number of base pairs) [100'000'000]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Available options for alignment:\n");
            fprintf(stderr, "\t-a --annotator [STR] \t\t\t\tannotator to load for label/trace-consistent alignment []\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\t\t\t\tbasename of output file []\n");
            fprintf(stderr, "\t   --json \t\t\t\t\toutput alignment in JSON format [off]\n");
if (advanced) {
            fprintf(stderr, "\t   --align-only-forwards \t\t\tdo not align backwards from a seed on basic-mode graphs [off]\n");
            fprintf(stderr, "\t   --align-no-seed-complexity-filter \t\t\t\tdisable the filter for low-complexity seeds. [off]\n");
            fprintf(stderr, "\t   --align-output-path \t\t\t\twith --json, also emit the VG-style path.mapping object [off]\n");
}
            fprintf(stderr, "\t   --align-alternative-alignments \t\tthe number of alternative paths to report per seed [1]\n");
            fprintf(stderr, "\t   --align-chain \t\t\t\tconstruct seed chains before alignment. Useful for long error-prone reads. [off]\n");
            fprintf(stderr, "\t   --align-post-chain \t\t\tperform multiple local alignments and chain them together into a single alignment. Useful for long error-prone reads. [off]\n");
            fprintf(stderr, "\t         \t\t\t\t\t\tA '$' inserted into the reference sequence indicates a jump in the graph.\n");
            fprintf(stderr, "\t         \t\t\t\t\t\tA 'G' in the reported CIGAR string indicates inserted graph nodes.\n");
if (advanced) {
            fprintf(stderr, "\t   --align-min-path-score [INT]\t\t\tmin score that a reported path can have [0]\n");
            fprintf(stderr, "\t   --align-max-nodes-per-seq-char [FLOAT]\tmaximum number of nodes to consider per sequence character [5.0]\n");
            fprintf(stderr, "\t   --align-max-ram [FLOAT]\t\t\tmaximum amount of RAM used per alignment in MB [200.0]\n");
}
            fprintf(stderr, "\t   --align-xdrop [INT]\t\t\t\tmaximum difference between the current score and the best alignment score [27, 100 if chaining is enabled]\n");
            fprintf(stderr, "\t   \t\t\t\t\t\t\tNote that this parameter should be scaled accordingly when changing the default scoring parameters.\n");
            fprintf(stderr, "\t   --align-rel-score-cutoff [FLOAT]\t\tmin score relative to the current best alignment to use as a lower bound for subsequent extensions [0.95]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Advanced options for scoring:\n");
            fprintf(stderr, "\t   --align-match-score [INT]\t\t\tpositive match score [2]\n");
            fprintf(stderr, "\t   --align-mm-transition-penalty [INT]\t\tpositive transition penalty (DNA only) [3]\n");
            fprintf(stderr, "\t   --align-mm-transversion-penalty [INT]\tpositive transversion penalty (DNA only) [3]\n");
            fprintf(stderr, "\t   --align-gap-open-penalty [INT]\t\tpositive gap opening penalty [6]\n");
            fprintf(stderr, "\t   --align-gap-extension-penalty [INT]\t\tpositive gap extension penalty [2]\n");
            fprintf(stderr, "\t   --align-end-bonus [INT]\t\tscore bonus for each endpoint of the query covered by an alignment [5]\n");
            fprintf(stderr, "\t   --align-edit-distance \t\t\tuse unit costs for scoring matrix [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Advanced options for seeding:\n");
            fprintf(stderr, "\t   --align-min-seed-length [INT]\t\tmin length of a seed [19]\n");
            fprintf(stderr, "\t   --align-max-seed-length [INT]\t\tmax length of a seed [inf]\n");
if (advanced) {
            fprintf(stderr, "\t   --align-min-exact-match [FLOAT] \t\tfraction of matching nucleotides required to align sequence [0.7]\n");
            fprintf(stderr, "\t   --align-max-num-seeds-per-locus [INT]\tmaximum number of allowed inexact seeds per locus [1000]\n");
}
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
            fprintf(stderr, "Usage: %s transform -o <outfile-base> [options] GRAPH\n"
                            "       %s transform --mask-dummy [--force] [-p INT] GRAPH.dbg\n\n", prog_name.c_str(), prog_name.c_str());

            // fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --index-ranges [INT]\tindex all node ranges in BOSS for suffixes of given length [%zu]\n", kDefaultIndexSuffixLen);
            fprintf(stderr, "\t   --clear-dummy \terase all redundant dummy edges and build an edgemask for non-redundant [off]\n");
            fprintf(stderr, "\t   --mask-dummy \twrite only the mask of dummy edges, <GRAPH without .dbg>.edgemask, as build --mask-dummy\n"
                            "\t                \twould: nothing pruned, the graph file untouched, node ids and annotation valid [off]\n");
            fprintf(stderr, "\t   --force \t\twith --mask-dummy: replace an existing .edgemask [off]\n");
            fprintf(stderr, "\t   --prune-tips [INT] \tprune all dead ends of this length and shorter [0]\n");
            fprintf(stderr, "\t   --state [STR] \tchange state of succinct graph: small / dynamic / stat / fast [stat]\n");
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
            fprintf(stderr, "\t-a --annotator [STR] \t\tannotator to load []\n");
            fprintf(stderr, "\t   --diff-assembly-rules [STR] \tJSON file describing labels to mask in and out and their relative fractions []\n");
            fprintf(stderr, "\t                       \t\tSee the manual for the specification.\n");
        } break;
        case STATS: {
            fprintf(stderr, "Usage: %s stats [options] GRAPH1 [[GRAPH2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for stats:\n");
            fprintf(stderr, "\t   --print \t\tprint graph table to the screen [off]\n");
            fprintf(stderr, "\t   --print-internal \tprint internal graph representation to screen [off]\n");
            fprintf(stderr, "\t   --count-quantiles [FLOAT ...] \tk-mer count quantiles to compute for each label [off]\n"
                            "\t                                 \t\tExample: --count-quantiles '0 0.33 0.5 0.66 1'\n"
                            "\t                                 \t\t(0 corresponds to MIN, 1 corresponds to MAX)\n");
            fprintf(stderr, "\t   --print-counts-hist \tprint histogram of k-mer weights as pairs (weight: num_kmers) [off]\n");
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
if (advanced) {
            fprintf(stderr, "\t   --sparse \t\tuse the row-major sparse matrix to annotate graph [off]\n");
}
            fprintf(stderr, "\t   --cache \t\tnumber of columns in cache (for column representation only) [10]\n");
            fprintf(stderr, "\t   --disk-swap [STR] \tdirectory to use for temporary files [off]\n");
            fprintf(stderr, "\t   --mem-cap-gb [FLOAT]\tbuffer size in GB (per column in construction) [1]\n");
            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file (or directory, for --separately) []\n");
            fprintf(stderr, "\t   --separately \tannotate each file independently and dump to the same directory [off]\n");
            fprintf(stderr, "\t   --threads-each [INT]\tnumber of threads to use when annotating each file with --separately [1]\n");
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
            fprintf(stderr, "\t   --index-header-coords \tgenerate a CoordToHeader mapping (.seqs file) from input FASTA files [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case MERGE_ANNOTATIONS: {
            fprintf(stderr, "Usage: %s merge_anno -o <annotation-basename> [options] ANNOT1 [[ANNOT2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case TRANSFORM_ANNOTATION: {
            fprintf(stderr, "Usage: %s transform_anno -o <annotation-basename> [options] ANNOTATOR\n\n", prog_name.c_str());

            // fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --aggregate-columns \t\taggregate annotation columns into a bitmask (new column) [off]\n");
            fprintf(stderr, "\t                       \t\t\tFormula: min-count <= \\sum_i 1{min-value <= c_i <= max-value} <= max-count\n");
            fprintf(stderr, "\t                       \t\t\tWith --count-kmers: min-count <= \\sum_i c_i 1{min-value <= c_i <= max-value} <= max-count\n");
            fprintf(stderr, "\t   --anno-label [STR]\t\tname of the aggregated output column [mask]\n");
            fprintf(stderr, "\t   --count-kmers \t\tsum up k-mer counts across columns [off]\n");
            fprintf(stderr, "\t   --count-width [INT] \t\tnumber of bits for aggregated k-mer counts (--count-kmers only) [8]\n");
            fprintf(stderr, "\t                       \t\tvalid range: [2, 32]. Values saturate at 2^W - 1.\n");
            fprintf(stderr, "\t   --min-value [INT] \t\tignore pre-aggregation counts smaller than this [1]\n");
            fprintf(stderr, "\t   --min-count [INT] \t\texclude k-mers with aggregated count smaller than this [1]\n");
            fprintf(stderr, "\t   --min-fraction [FLOAT] \texclude k-mers appearing in fewer than this fraction of columns [0.0]\n");
            fprintf(stderr, "\t                          \t\tignored in --count-kmers mode, use --min-count instead\n");
            fprintf(stderr, "\t   --max-value [INT] \t\tignore pre-aggregation counts larger than this [inf]\n");
            fprintf(stderr, "\t   --max-count [INT] \t\texclude k-mers with aggregated count larger than this [inf]\n");
            fprintf(stderr, "\t                       \t\tin --count-kmers mode: max_count+1 must fit in count_width bits (unless inf)\n");
            fprintf(stderr, "\t   --max-fraction [FLOAT] \texclude k-mers appearing in more than this fraction of columns [1.0]\n");
            fprintf(stderr, "\t                          \t\tignored in --count-kmers mode, use --max-count instead\n");
            fprintf(stderr, "\t   --compute-overlap [STR] \tcompute the number of shared bits in columns of this annotation and ANNOTATOR [off]\n");
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
            fprintf(stderr, "\t   --subsample [INT] \tnumber of bits subsampled for distance estimation in column clustering [1'000'000]\n");
            fprintf(stderr, "\t   --subsample-rows \tsubsample rows (the same positions in all columns) instead of only set bits [off]\n");
            fprintf(stderr, "\t   --dump-text-anno \tdump the columns of the annotator as separate text files [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --row-diff-stage [0|1|2] \tstage of the row_diff construction [0]\n");
            fprintf(stderr, "\t   --max-path-length [INT] \tmaximum path length in row_diff annotation [100]\n");
            fprintf(stderr, "\t   --mem-cap-gb [FLOAT]\tmemory in GB available for the transform [1000]\n");
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
                            "\tEach input file is given in FASTA or FASTQ format.\n"
                            "\tOutput format: tsv with rows '<query id>\t<query name>\t<results ...>'.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for query:\n");
#if ! _PROTEIN_GRAPH
            fprintf(stderr, "\t   --fwd-and-reverse \tfor each input sequence, query its reverse complement as well [off]\n");
#endif
            fprintf(stderr, "\t   --align \t\talign sequences instead of mapping k-mers [off]\n");
if (advanced) {
            fprintf(stderr, "\t   --sparse \t\tuse row-major sparse matrix for row annotation [off]\n");
}
            fprintf(stderr, "\t   --json \t\toutput query results in JSON format [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t   --query-mode \tquery mode (only labels with enough k-mer matches are reported) [%s]\n", querymode_to_string(MATCHES).c_str());
            fprintf(stderr, "\t       Available modes:\n");
            fprintf(stderr, "\t                %s \t\tprint labels (with enough k-mer matches)\n", querymode_to_string(LABELS).c_str());
            fprintf(stderr, "\t                %s \tprint labels and the number of k-mer matches (for every label with enough k-mer matches)\n", querymode_to_string(MATCHES).c_str());
            fprintf(stderr, "\t                %s \tprint masks indicating present/absent k-mers\n", querymode_to_string(SIGNATURE).c_str());
            fprintf(stderr, "\t                \t\t\t\tOutput format: run-length encoding 'x<N>o<N>...' where x=present, o=absent\n"
                            "\t                \t\t\t\t    e.g. 'x3o2x2o1' for bitmask '11100110'\n"
                            "\t                \t\t\t\tWith --verbose-output: full binary string, e.g. '11100110'\n");
if (advanced) {
            fprintf(stderr, "\t                %s \tprint sum of counts for the matched k-mers, requires count or coord annotation (...)\n", querymode_to_string(COUNTS_SUM).c_str());
}
            fprintf(stderr, "\t                %s \t\tprint k-mer counts, requires count or coord annotation (...)\n", querymode_to_string(COUNTS).c_str());
            fprintf(stderr, "\t                \t\t\t\tOutput format: '<pos in query>=<abundance>' (single k-mer match)\n"
                            "\t                \t\t\t\t    or '<first pos>-<last pos>=<abundance>' (segment match)\n"
                            "\t                \t\t\t\tAll positions start with 0\n");
            fprintf(stderr, "\t                %s \t\tprint k-mer coordinates, requires coord annotation (...)\n", querymode_to_string(COORDS).c_str());
            fprintf(stderr, "\t                \t\t\t\tOutput format: '<pos in query>-<pos in sample>' (single k-mer match)\n"
                            "\t                \t\t\t\t    or '<start pos in query>-<first pos in sample>-<last pos in sample>' (segment match)\n"
                            "\t                \t\t\t\tAll positions start with 0\n");
if (advanced) {
            fprintf(stderr, "\t   --verbose-output \t\tfor coords/counts: do not collapse ranges; for signature: print full mask instead of RLE [off]\n");
}
            fprintf(stderr, "\t   --num-top-labels [INT] \t\tmaximum number of top labels to output [inf]\n");
            fprintf(stderr, "\t   --min-kmers-fraction-label [FLOAT] \tmin fraction of k-mers from the query required to be present in a label [0.7]\n");
            fprintf(stderr, "\t   --min-kmers-fraction-graph [FLOAT] \tmin fraction of k-mers from the query required to be present in the graph [0.0]\n");
            fprintf(stderr, "\t   --no-coord-mapping \t\t\tquery without mapping coords to sequence headers even if the .seqs index exists [off]\n");
if (advanced) {
            fprintf(stderr, "\t   --labels-delimiter [STR]\tdelimiter for annotation labels [\":\"]\n");
            fprintf(stderr, "\t   --suppress-unlabeled \tdo not show results for sequences missing in graph [off]\n");
}
            // fprintf(stderr, "\t-d --distance [INT] \tmax allowed alignment distance [0]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
            // fprintf(stderr, "\t   --cache-size [INT] \tnumber of uncompressed rows to store in the cache [0]\n");
            fprintf(stderr, "\t   --batch-size [INT] \tquery batch size in bp (0 to disable batch query) [100'000'000]\n");
if (advanced) {
            fprintf(stderr, "\t   --threads-each [INT]\tnumber of parallel batches [1]\n");
            fprintf(stderr, "\t   --RA-ivbuff-size [INT] \tsize (in bytes) of int_vector_buffer used in random access mode (e.g. by row disk annotator) [16384]\n");
}
            fprintf(stderr, "\t   --one-pass-brwt \tuse one-pass parallel BRWT traversal for queries [off]\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "Available options for --align:\n");
if (advanced) {
            fprintf(stderr, "\t   --align-only-forwards \t\t\tdo not align backwards from a seed on basic-mode graphs [off]\n");
}
            // fprintf(stderr, "\t   --align-alternative-alignments \tthe number of alternative paths to report per seed [1]\n");
            fprintf(stderr, "\t   --align-min-path-score [INT]\t\t\tmin score that a reported path can have [0]\n");
if (advanced) {
            fprintf(stderr, "\t   --align-max-nodes-per-seq-char [FLOAT]\tmaximum number of nodes to consider per sequence character [5.0]\n");
            fprintf(stderr, "\t   --align-max-ram [FLOAT]\t\t\tmaximum amount of RAM used per alignment in MB [200.0]\n");
}
            fprintf(stderr, "\t   --align-xdrop [INT]\t\t\t\tmaximum difference between the current score and the best alignment score [27, 100 if chaining is enabled]\n");
            fprintf(stderr, "\t   \t\t\t\t\t\t\tNote that this parameter should be scaled accordingly when changing the default scoring parameters.\n");
            fprintf(stderr, "\n");
if (advanced) {
            fprintf(stderr, "\t   --batch-align \t\talign against query graph [off]\n");
            fprintf(stderr, "\t   --max-hull-forks [INT]\tmaximum number of forks to take when expanding query graph [4]\n");
            fprintf(stderr, "\t   --max-hull-depth [INT]\tmaximum number of steps to traverse when expanding query graph [max_nodes_per_seq_char * max_seq_len]\n");
            fprintf(stderr, "\n");
}
            fprintf(stderr, "\tAdvanced options for scoring:\n");
            fprintf(stderr, "\t   --align-match-score [INT]\t\t\tpositive match score [2]\n");
            fprintf(stderr, "\t   --align-mm-transition-penalty [INT]\t\tpositive transition penalty (DNA only) [3]\n");
            fprintf(stderr, "\t   --align-mm-transversion-penalty [INT]\tpositive transversion penalty (DNA only) [3]\n");
            fprintf(stderr, "\t   --align-gap-open-penalty [INT]\t\tpositive gap opening penalty [6]\n");
            fprintf(stderr, "\t   --align-gap-extension-penalty [INT]\t\tpositive gap extension penalty [2]\n");
if (advanced) {
            fprintf(stderr, "\t   --align-end-bonus [INT]\t\tscore bonus for each endpoint of the query covered by an alignment [5]\n");
            fprintf(stderr, "\t   --align-edit-distance \t\t\tuse unit costs for scoring matrix [off]\n");
}
            fprintf(stderr, "\n");
            fprintf(stderr, "\tAdvanced options for seeding:\n");
            fprintf(stderr, "\t   --align-min-seed-length [INT]\t\tmin length of a seed [19]\n");
            fprintf(stderr, "\t   --align-max-seed-length [INT]\t\tmax length of a seed [inf]\n");
            fprintf(stderr, "\t   --align-min-exact-match [FLOAT]\t\tfraction of matching nucleotides required to align sequence [0.7]\n");
if (advanced) {
            fprintf(stderr, "\t   --align-max-num-seeds-per-locus [INT]\tmaximum number of allowed inexact seeds per locus [1000]\n");
}
        } break;
        case TRAVERSE: {
            fprintf(stderr, "Usage: %s traverse [options] -i <GRAPH> -a <ANNOTATION> REQUEST.json [[REQUEST2.json] ...]\n"
                            "\tEach request is a JSON object as documented in\n"
                            "\tdocs/SPEC-labeled-traversal-core.md (§5 for traverse, §4.1 for --resolve).\n"
                            "\tOne JSON result is written to stdout per request.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for traverse:\n");
            fprintf(stderr, "\t   --resolve \t\t\treport label support and seed candidates instead of traversing [off]\n");
            fprintf(stderr, "\t   --index-release [STR]\trelease id echoed in results; requests may pin it []\n");
            fprintf(stderr, "\t   --index-name [STR]\t\tname of the index in capabilities and graphlets, [A-Za-z0-9._-]+ []\n");
            fprintf(stderr, "\t   --index-manifest [FILE]\tmanifest of the index bundle (files with size and sha256); its digest is the index identity []\n");
            fprintf(stderr, "\t   --index-inventory \t\tprint the files the loaders open for -i / -a (JSON: the identity files a manifest lists, and apart the graph's derived mask and Bloom filter, which it never lists) and exit, loading nothing [off]\n");
            fprintf(stderr, "\t   --traverse-chunk-target-ms [INT]\tdecode an annotation read a deadline may fall into in chunks of about this duration, the deadline checked between them; 0 = one piece per read [50]\n");
            fprintf(stderr, "\t   --traverse-path-cache-mb [INT]\tbound of a request's row-diff path cache (decoded rows kept so that later reads stop their row-diff paths at them), within the label cache's allotment under a memory budget; 0 = off [128]\n");
            fprintf(stderr, "\t   --json \t\t\tprint compact JSON (one line per request) [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for loading [1]\n");
            fprintf(stderr, "\n");
            return;
        }
        case PATTERN: {
            fprintf(stderr, "Usage: %s pattern [options] -i <GRAPH> -a <ANNOTATION> REQUEST.json [[REQUEST2.json] ...]\n"
                            "\tEach request is the JSON body of POST /pattern (docs/DESIGN-pattern-search.md §7);\n"
                            "\tone JSON answer is written to stdout per request, the server's, under the same caps.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for pattern:\n");
            fprintf(stderr, "\t   --index-release [STR]\trelease id echoed in answers []\n");
            fprintf(stderr, "\t   --index-name [STR]\t\tname of the index in answers, [A-Za-z0-9._-]+ []\n");
            fprintf(stderr, "\t   --index-manifest [FILE]\tmanifest of the index bundle; its digest is the index identity []\n");
            fprintf(stderr, "\t   --pattern-min-information-bits [FLOAT] \tinformation floor of a pattern (bits; an exact pattern in suffix scope is exempt) [24]\n");
            fprintf(stderr, "\t   --pattern-max-contexts [INT] \tdefault and maximum of max_contexts per pattern (retrieval threshold, partial's cap) [10000]\n");
            fprintf(stderr, "\t   --pattern-max-anchors [INT] \tdefault and maximum of max_anchors per pattern longer than k [1000]\n");
            fprintf(stderr, "\t   --pattern-max-paths [INT] \tdefault and maximum of max_paths per pattern longer than k with long_search paths (retrieval threshold, partial's cap) [1000]\n");
            fprintf(stderr, "\t   --pattern-max-steps [INT] \tdefault and maximum of max_steps per request (range and mask-scan steps) [100000000]\n");
            fprintf(stderr, "\t   --pattern-default-time-ms [INT] \ttime_budget_ms of a request that names none [60000]\n");
            fprintf(stderr, "\t   --pattern-max-time-ms [INT] \tmaximum of time_budget_ms per request [600000]\n");
            fprintf(stderr, "\t   --pattern-finalize-ms [INT] \tfinalisation reserve inside the time budget: work stops at least this long before the deadline, longer by the estimated time to write what the answer holds (raise it with --pattern-max-contexts or --pattern-max-patterns, or on a slow host) [250]\n");
            fprintf(stderr, "\t   --pattern-delivery-build-mbps [FLOAT] \trate (MB/s) at which an answer's JSON text is assumed to be built and written, for that estimate [10]\n");
            fprintf(stderr, "\t   --pattern-delivery-compress-mbps [FLOAT] \trate (MB/s) at which it is assumed to be compressed [50]\n");
            fprintf(stderr, "\t   --pattern-max-patterns [INT] \tpatterns per request (a longer list is refused) [16]\n");
            fprintf(stderr, "\t   --pattern-max-labels-per-anchor [INT] \tdefault and maximum of max_labels_per_anchor: labels kept per row with output.labels all (more: a truncated anchor) [64]\n");
            fprintf(stderr, "\t   --pattern-max-annotation-work [INT] \tdefault and maximum of max_annotation_work per request (annotation work units) [100000000]\n");
            fprintf(stderr, "\t   --pattern-max-memory-mb [INT] \tdefault and maximum of max_memory_mb: the memory account of a request reading labels [256]\n");
            fprintf(stderr, "\t   --pattern-max-labels [INT] \tdefault and maximum of max_labels: labels listed per pattern in mode partial [1000]\n");
            fprintf(stderr, "\t   --pattern-max-occurrences [INT] \tdefault and maximum of max_occurrences_per_label in mode partial [16]\n");
            fprintf(stderr, "\t   --pattern-build-mask \tbuild the dummy-edge mask in memory at load when the graph has no .edgemask, with -p threads, for exact counts (without a mask: upper bounds with estimates; small graphs; else transform --mask-dummy once) [off]\n");
            fprintf(stderr, "\t   --json \t\t\tprint compact JSON (one line per request) [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for loading [1]\n");
            fprintf(stderr, "\n");
            return;
        }
        case SERVER_QUERY: {
            fprintf(stderr, "Usage: %s server_query (-i <GRAPH> -a <ANNOTATION> | <GRAPHS.csv>) [options]\n\n"
                            "\tThe index must be passed with flags -i -a or with a file GRAPHS.csv listing one\n"
                            "\tor more indexes, a file with rows: '<name>,<graph_path>,<annotation_path>\\n'.\n"
                            "\t(If multiple rows have the same name, all those graphs will be queried for that name.)\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for server_query:\n");
            fprintf(stderr, "\t   --port [INT] \tTCP port for incoming connections [5555]\n");
            fprintf(stderr, "\t   --address \t\tinterface for incoming connections (default: all)\n");
if (advanced) {
            fprintf(stderr, "\t   --sparse \t\tuse the row-major sparse matrix to annotate graph [off]\n");
}
            // fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file []\n");
            // fprintf(stderr, "\t-d --distance [INT] \tmax allowed alignment distance [0]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tmaximum number of parallel connections [1]\n");
            fprintf(stderr, "\t   --threads-each [INT] \tnumber of threads per graph [1]\n");
            fprintf(stderr, "\t   --one-pass-brwt \tuse one-pass parallel BRWT traversal for queries [off]\n");
            // fprintf(stderr, "\t   --cache-size [INT] \tnumber of uncompressed rows to store in the cache [0]\n");
            fprintf(stderr, "\n\t   --num-top-labels [INT] \tmaximum number of top labels per query by default [10'000]\n");
            fprintf(stderr, "\t   --no-coord-mapping \t\tquery without mapping coords to sequence headers even if the .seqs index exists [off]\n");
            fprintf(stderr, "\t   --mem-cap-gb [FLOAT] \tmemory in GB available for the server to load graphs for queries into RAM [0]\n");
            fprintf(stderr, "\n\t   --index-release [STR] \trelease id echoed by /traverse and /resolve; requests may pin it []\n");
            fprintf(stderr, "\t   --index-name [STR] \t\tname of the index (-i / -a only) in capabilities and graphlets, [A-Za-z0-9._-]+ []\n");
            fprintf(stderr, "\t   --index-manifest [FILE] \tmanifest of the index bundle (-i / -a only); its digest is the index identity []\n");
            fprintf(stderr, "\t   --traverse-max-time-ms [FLOAT] \tcap on bounds.time_budget_ms per /traverse seed, 0 = unlimited [30000]\n");
            fprintf(stderr, "\t   --traverse-max-seeds [INT] \t\tcap on seeds per /traverse request, 0 = unlimited [64]\n");
            fprintf(stderr, "\t   --traverse-max-seed-bp [INT] \t\tcap on the length of a /traverse seed, 0 = unlimited [100000]\n");
            fprintf(stderr, "\t   --traverse-max-seed-labels [INT] \tcap on labels DERIVED from a seed, 0 = unlimited [10000]\n");
            fprintf(stderr, "\t   --resolve-max-query-bp [INT] \t\tcap on the /resolve query length, 0 = unlimited [0]\n");
            fprintf(stderr, "\t   --traverse-max-memory-mb [INT] \tmaximum of bounds.max_memory_mb per /traverse seed: a larger one is lowered to it, an omitted one set to it; 0 = off [0]\n");
            fprintf(stderr, "\t   --traverse-max-work-units [INT] \tmaximum of bounds.max_work_units per /traverse seed, likewise; 0 = off [0]\n");
            fprintf(stderr, "\t   --traverse-attempt-allowance-ms [INT] \tadded to seeds x time budget in the bound enforced on a /traverse with attempt_id [10000]\n");
            fprintf(stderr, "\t   --traverse-attempt-retention-s [INT] \tfinished attempts stay queryable this long (and their ids refused; one sent with not_after_ms is held until not_after_ms + the clock skew allowance, within the tombstone cap), and a cancel of an unknown id tombstones it at least this long, in [0, 31536000]; 0 = nothing kept, no tombstones (such a cancel is refused, 429) [3600]\n");
            fprintf(stderr, "\t   --traverse-attempt-retention [INT] \tat most this many finished attempts are kept (oldest dropped first; one sent with not_after_ms stays held, its id refused, as --traverse-attempt-retention-s says), and a cancel of an unknown id is tombstoned only while fewer than this many tombstones and held attempts are held (else it is refused, 429), in [0, 10000000] [10000]\n");
            fprintf(stderr, "\t   --traverse-attempt-tombstone-max-s [INT] \tthe longest a tombstone, or a finished attempt's id, is held to cover a not_after_ms (+ the clock skew allowance), from the latest cancel, finish or refused copy of the request, in [0, 31536000]; never less than --traverse-attempt-retention-s [86400]\n");
            fprintf(stderr, "\t   --traverse-clock-skew-ms [INT] \tclock skew a ledger adds to not_after_ms, stated in the capabilities [2000]\n");
            fprintf(stderr, "\t   --traverse-chunk-target-ms [INT] \tdecode an annotation read a deadline may fall into in chunks of about this duration, the deadline checked between them; 0 = one piece per read [50]\n");
            fprintf(stderr, "\t   --traverse-path-cache-mb [INT] \tbound of a /traverse request's row-diff path cache (decoded rows kept so that later reads stop their row-diff paths at them), within the label cache's allotment under a memory budget; 0 = off [128]\n");
            fprintf(stderr, "\t   --traverse-compression-level [INT] \tzlib level (1-9) of the traversal routes' compressed bodies; the other routes use 9 [1]\n");
            fprintf(stderr, "\t   --traverse-delivery-compress-mbps [FLOAT] \tcompression rate the delivery reserve of an attempt assumes [50]\n");
            fprintf(stderr, "\t   --traverse-delivery-build-mbps [FLOAT] \tresponse-building rate it assumes until the attempt measures its own [10]\n");
            fprintf(stderr, "\t   --pattern-min-information-bits [FLOAT] \tinformation floor of a pattern (bits; an exact pattern in suffix scope is exempt) [24]\n");
            fprintf(stderr, "\t   --pattern-max-contexts [INT] \tdefault and maximum of max_contexts per pattern (retrieval threshold, partial's cap) [10000]\n");
            fprintf(stderr, "\t   --pattern-max-anchors [INT] \tdefault and maximum of max_anchors per pattern longer than k [1000]\n");
            fprintf(stderr, "\t   --pattern-max-paths [INT] \tdefault and maximum of max_paths per pattern longer than k with long_search paths (retrieval threshold, partial's cap) [1000]\n");
            fprintf(stderr, "\t   --pattern-max-steps [INT] \tdefault and maximum of max_steps per request (range and mask-scan steps) [100000000]\n");
            fprintf(stderr, "\t   --pattern-default-time-ms [INT] \ttime_budget_ms of a request that names none [60000]\n");
            fprintf(stderr, "\t   --pattern-max-time-ms [INT] \tmaximum of time_budget_ms per request, at most 899000 (the 900 s content timeout less 1 s for the transport) [600000]\n");
            fprintf(stderr, "\t   --pattern-finalize-ms [INT] \tfinalisation reserve inside the time budget: work stops at least this long before the deadline, longer by the estimated time to write what the answer holds (raise it with --pattern-max-contexts or --pattern-max-patterns, or on a slow host) [250]\n");
            fprintf(stderr, "\t   --pattern-delivery-build-mbps [FLOAT] \trate (MB/s) at which an answer's JSON text is assumed to be built and written, for that estimate [10]\n");
            fprintf(stderr, "\t   --pattern-delivery-compress-mbps [FLOAT] \trate (MB/s) at which it is assumed to be compressed [50]\n");
            fprintf(stderr, "\t   --pattern-max-patterns [INT] \tpatterns per request (a longer list is refused) [16]\n");
            fprintf(stderr, "\t   --pattern-max-labels-per-anchor [INT] \tdefault and maximum of max_labels_per_anchor: labels kept per row with output.labels all (more: a truncated anchor) [64]\n");
            fprintf(stderr, "\t   --pattern-max-annotation-work [INT] \tdefault and maximum of max_annotation_work per request (annotation work units) [100000000]\n");
            fprintf(stderr, "\t   --pattern-max-memory-mb [INT] \tdefault and maximum of max_memory_mb: the memory account of a request reading labels [256]\n");
            fprintf(stderr, "\t   --pattern-max-labels [INT] \tdefault and maximum of max_labels: labels listed per pattern in mode partial [1000]\n");
            fprintf(stderr, "\t   --pattern-max-occurrences [INT] \tdefault and maximum of max_occurrences_per_label in mode partial [16]\n");
            fprintf(stderr, "\t   --pattern-build-mask \tbuild the dummy-edge mask in memory at load when the graph has no .edgemask, with --threads-each threads, for exact counts (without a mask: upper bounds with estimates; small graphs; else transform --mask-dummy once; a mask written later is read only after a restart) [off]\n");
        } break;
    }

    fprintf(stderr, "\nGeneral options:\n");
    fprintf(stderr, "\t   --mmap \t\tuse memory mapping when loading to reduce RAM [off]\n");
    if (identity == SERVER_QUERY)
        fprintf(stderr, "\t   --madv-random \tenable MADV_RANDOM hints for graphs loaded with mmap (speeds up tiny queries, may slow down large ones) [off]\n");
    fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
    fprintf(stderr, "\t   --advanced \t\tshow other advanced and legacy options [off]\n");
    fprintf(stderr, "\t-h --help \t\tprint usage info\n");
    fprintf(stderr, "\n");
}

} // namespace cli
} // namespace mtg
