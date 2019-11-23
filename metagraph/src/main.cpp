#include <filesystem>
#include <json/json.h>
#include <ips4o.hpp>
#include <fmt/format.h>

#include "unix_tools.hpp"
#include "config.hpp"
#include "sequence_io.hpp"
#include "boss_construct.hpp"
#include "boss_chunk.hpp"
#include "boss_merge.hpp"
#include "annotated_dbg.hpp"
#include "annotate_row_compressed.hpp"
#include "annotate_column_compressed.hpp"
#include "serialization.hpp"
#include "algorithms.hpp"
#include "string_utils.hpp"
#include "file_utils.hpp"
#include "threading.hpp"
#include "reverse_complement.hpp"
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"
#include "kmc_parser.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_bitmap.hpp"
#include "dbg_bitmap_construct.hpp"
#include "dbg_succinct.hpp"
#include "graph_cleaning.hpp"
#include "dbg_aligner.hpp"
#include "aligner_methods.hpp"
#include "server.hpp"
#include "node_weights.hpp"
#include "masked_graph.hpp"
#include "annotated_graph_algorithm.hpp"
#include "taxid_mapper.hpp"
#include <typeinfo>

typedef annotate::MultiLabelEncoded<uint64_t, std::string> Annotator;

const size_t kNumCachedColumns = 10;
const size_t kBitsPerCount = 8;
static const size_t kRowBatchSize = 100'000;


Config::GraphType parse_graph_extension(const std::string &filename) {
    if (utils::ends_with(filename, ".dbg")) {
        return Config::GraphType::SUCCINCT;

    } else if (utils::ends_with(filename, ".orhashdbg")) {
        return Config::GraphType::HASH;

    } else if (utils::ends_with(filename, ".hashstrdbg")) {
        return Config::GraphType::HASH_STR;

    } else if (utils::ends_with(filename, ".bitmapdbg")) {
        return Config::GraphType::BITMAP;

    } else {
        return Config::GraphType::INVALID;
    }
}

Config::AnnotationType parse_annotation_type(const std::string &filename) {
    if (utils::ends_with(filename, annotate::kColumnAnnotatorExtension)) {
        return Config::AnnotationType::ColumnCompressed;

    } else if (utils::ends_with(filename, annotate::kRowAnnotatorExtension)) {
        return Config::AnnotationType::RowCompressed;

    } else if (utils::ends_with(filename, annotate::kBRWTExtension)) {
        return Config::AnnotationType::BRWT;

    } else if (utils::ends_with(filename, annotate::kBinRelWT_sdslExtension)) {
        return Config::AnnotationType::BinRelWT_sdsl;

    } else if (utils::ends_with(filename, annotate::kBinRelWTExtension)) {
        return Config::AnnotationType::BinRelWT;

    } else if (utils::ends_with(filename, annotate::kRowPackedExtension)) {
        return Config::AnnotationType::RowFlat;

    } else if (utils::ends_with(filename, annotate::kRainbowfishExtension)) {
        return Config::AnnotationType::RBFish;

    } else {
        std::cerr << "Error: unknown annotation format in "
                  << filename << std::endl;
        exit(1);
    }
}

std::string remove_graph_extension(const std::string &filename) {
    return utils::remove_suffix(filename, ".dbg",
                                          ".orhashdbg",
                                          ".hashstrdbg",
                                          ".bitmapdbg");
}

template <class Graph = BOSS>
std::shared_ptr<Graph> load_critical_graph_from_file(const std::string &filename) {
    auto graph = std::make_shared<Graph>(2);
    if (!graph->load(filename)) {
        std::cerr << "ERROR: can't load graph from file " << filename << std::endl;
        exit(1);
    }
    return graph;
}

std::shared_ptr<DeBruijnGraph> load_critical_dbg(const std::string &filename) {
    auto graph_type = parse_graph_extension(filename);
    switch (graph_type) {
        case Config::GraphType::SUCCINCT:
            return load_critical_graph_from_file<DBGSuccinct>(filename);

        case Config::GraphType::HASH:
            return load_critical_graph_from_file<DBGHashOrdered>(filename);

        case Config::GraphType::HASH_PACKED:
            return load_critical_graph_from_file<DBGHashOrdered>(filename);

        case Config::GraphType::HASH_STR:
            return load_critical_graph_from_file<DBGHashString>(filename);

        case Config::GraphType::BITMAP:
            return load_critical_graph_from_file<DBGBitmap>(filename);

        case Config::GraphType::INVALID:
            std::cerr << "ERROR: can't load graph from file '"
                      << filename
                      << "', needs valid file extension" << std::endl;
            exit(1);
    }
    assert(false);
    exit(1);
}

void annotate_data(const std::vector<std::string> &files,
                   const std::string &ref_sequence_path,
                   AnnotatedDBG *anno_graph,
                   bool forward_and_reverse,
                   size_t min_count,
                   size_t max_count,
                   bool filename_anno,
                   bool fasta_anno,
                   const std::string &fasta_anno_comment_delim,
                   const std::string &fasta_header_delimiter,
                   const std::vector<std::string> &anno_labels,
                   bool verbose) {
    size_t total_seqs = 0;

    Timer timer;

    // iterate over input files
    for (const auto &file : files) {
        Timer data_reading_timer;

        if (verbose) {
            std::cout << std::endl << "Parsing " << file << std::endl;
        }
        // read files
        if (utils::get_filetype(file) == "VCF") {
            read_vcf_file_with_annotations_critical(
                file,
                ref_sequence_path,
                dynamic_cast<const DeBruijnGraph &>(anno_graph->get_graph()).get_k(),
                [&](auto&& seq, const auto &variant_labels) {
                    std::vector<std::string> labels(variant_labels.begin(),
                                                    variant_labels.end());

                    if (filename_anno)
                        labels.push_back(file);

                    for (const auto &label : anno_labels) {
                        labels.push_back(label);
                    }

                    anno_graph->annotate_sequence(seq, labels);
                },
                forward_and_reverse
            );
        } else if (utils::get_filetype(file) == "KMC") {
            std::vector<std::string> labels;

            if (filename_anno) {
                labels.push_back(file);
            }

            for (const auto &label : anno_labels) {
                labels.push_back(label);
            }

            kmc::read_kmers(
                file,
                [&](std::string&& sequence) {
                    anno_graph->annotate_sequence(std::move(sequence), labels);

                    total_seqs += 1;
                    if (verbose && total_seqs % 10000 == 0) {
                        std::cout << "processed " << total_seqs << " sequences"
                                  << ", trying to annotate as ";
                        for (const auto &label : labels) {
                            std::cout << "<" << label << ">";
                        }
                        std::cout << ", " << timer.elapsed() << "sec" << std::endl;
                    }
                },
                !dynamic_cast<const DeBruijnGraph&>(anno_graph->get_graph()).is_canonical_mode(),
                min_count,
                max_count
            );
        } else if (utils::get_filetype(file) == "FASTA"
                    || utils::get_filetype(file) == "FASTQ") {
            read_fasta_file_critical(file,
                [&](kseq_t *read_stream) {
                    std::vector<std::string> labels;

                    if (fasta_anno) {
                        labels = utils::split_string(
                            fasta_anno_comment_delim != Config::UNINITIALIZED_STR
                                && read_stream->comment.l
                                    ? utils::join_strings(
                                        { read_stream->name.s, read_stream->comment.s },
                                        fasta_anno_comment_delim,
                                        true)
                                    : read_stream->name.s,
                            fasta_header_delimiter
                        );
                    }
                    if (filename_anno) {
                        labels.push_back(file);
                    }

                    for (const auto &label : anno_labels) {
                        labels.push_back(label);
                    }

                    anno_graph->annotate_sequence(read_stream->seq.s, labels);

                    total_seqs += 1;
                    if (verbose && total_seqs % 10000 == 0) {
                        std::cout << "processed " << total_seqs << " sequences"
                                  << ", last was " << read_stream->name.s
                                  << ", trying to annotate as ";
                        for (const auto &label : labels) {
                            std::cout << "<" << label << ">";
                        }
                        std::cout << ", " << timer.elapsed() << "sec" << std::endl;
                    }
                },
                forward_and_reverse
            );
        } else {
            std::cerr << "ERROR: Filetype unknown for file "
                      << file << std::endl;
            exit(1);
        }

        if (verbose) {
            std::cout << "File processed in "
                      << data_reading_timer.elapsed()
                      << "sec, current mem usage: "
                      << (get_curr_RSS() >> 20) << " MiB"
                      << ", total time: " << timer.elapsed()
                      << "sec" << std::endl;
        }
    }

    // join threads if any were initialized
    anno_graph->join();
}


void annotate_coordinates(const std::vector<std::string> &files,
                          AnnotatedDBG *anno_graph,
                          bool forward_and_reverse,
                          size_t genome_bin_size,
                          bool verbose) {
    size_t total_seqs = 0;

    Timer timer;

    const size_t k = dynamic_cast<const DeBruijnGraph &>(anno_graph->get_graph()).get_k();

    // iterate over input files
    for (const auto &file : files) {
        Timer data_reading_timer;

        if (verbose)
            std::cout << std::endl << "Parsing " << file << std::endl;

        // open stream
        if (utils::get_filetype(file) == "FASTA"
                    || utils::get_filetype(file) == "FASTQ") {

            bool forward_strand = true;

            read_fasta_file_critical(file,
                [&](kseq_t *read_stream) {
                    std::vector<std::string> labels {
                        file,
                        read_stream->name.s,
                        std::to_string(forward_and_reverse && (total_seqs % 2)), // whether the read is reverse
                        "",
                    };

                    const std::string sequence(read_stream->seq.s);
                    for (size_t i = 0; i < sequence.size(); i += genome_bin_size) {
                        labels.back() = fmt::format_int(i).c_str();

                        // forward: |0 =>  |6 =>  |12=>  |18=>  |24=>  |30=>|
                        // reverse: |<=30|  <=24|  <=18|  <=12|  <= 6|  <= 0|
                        const size_t bin_size = std::min(
                            static_cast<size_t>(sequence.size() - i),
                            static_cast<size_t>(genome_bin_size + k - 1)
                        );
                        anno_graph->annotate_sequence(
                            forward_strand
                                ? sequence.substr(i, bin_size)
                                : sequence.substr(sequence.size() - i - bin_size, bin_size),
                            { utils::join_strings(labels, "\1"), }
                        );
                    }

                    total_seqs += 1;
                    if (verbose && total_seqs % 10000 == 0) {
                        std::cout << "processed " << total_seqs << " sequences"
                                  << ", last was " << read_stream->name.s
                                  << ", trying to annotate as ";
                        for (const auto &label : labels) {
                            std::cout << "<" << label << ">";
                        }
                        std::cout << ", " << timer.elapsed() << "sec" << std::endl;
                    }

                    // If we read both strands, the next sequence is
                    // either reverse (if the current one is forward)
                    // or new (if the current one is reverse), and therefore forward
                    if (forward_and_reverse)
                        forward_strand = !forward_strand;
                },
                forward_and_reverse
            );
        } else {
            std::cerr << "ERROR: the type of file "
                      << file << " is not supported" << std::endl;
            exit(1);
        }

        if (verbose) {
            std::cout << "File processed in "
                      << data_reading_timer.elapsed()
                      << "sec, current mem usage: "
                      << (get_curr_RSS() >> 20) << " MiB"
                      << ", total time: " << timer.elapsed()
                      << "sec" << std::endl;
        }
    }

    // join threads if any were initialized
    anno_graph->join();
}

void execute_query(const std::string &seq_name,
                   const std::string &sequence,
                   bool count_labels,
                   bool suppress_unlabeled,
                   size_t num_top_labels,
                   double discovery_fraction,
                   std::string anno_labels_delimiter,
                   const AnnotatedDBG &anno_graph,
                   std::ostream &output_stream,
                   IDBGAligner *aligner = nullptr) {
    std::vector<std::string> sequences;

    std::vector<double> weights;

    if (aligner) {
        auto alignments = aligner->align(sequence);
        sequences.reserve(alignments.size());
        weights.reserve(alignments.size());

        for (const auto &alignment : alignments) {
            sequences.emplace_back(alignment.get_sequence());
            weights.emplace_back(std::exp(alignment.get_score()
                - aligner->get_config().match_score(sequences.back().begin(),
                                                    sequences.back().end())));
        }
    }

    assert(sequences.size() == weights.size());
    assert(!aligner || sequences.size());

    std::string output;
    output.reserve(1'000);

    if (count_labels) {
        auto top_labels = aligner
            ? anno_graph.get_top_labels(sequences,
                                        weights,
                                        num_top_labels,
                                        discovery_fraction)
            : anno_graph.get_top_labels(sequence, num_top_labels, discovery_fraction);

        if (!top_labels.size() && suppress_unlabeled)
            return;

        output += seq_name;

        for (const auto &[label, count] : top_labels) {
            output += "\t<";
            output += label;
            output += ">:";
            output += fmt::format_int(count).c_str();
        }

        output += '\n';

    } else {
        auto labels_discovered = aligner
            ? anno_graph.get_labels(sequences, weights, discovery_fraction)
            : anno_graph.get_labels(sequence, discovery_fraction);

        if (!labels_discovered.size() && suppress_unlabeled)
            return;

        output += seq_name;
        output += '\t';
        output += utils::join_strings(labels_discovered,
                                      anno_labels_delimiter);
        output += '\n';
    }

    output_stream << output;
}

std::unique_ptr<Annotator> initialize_annotation(Config::AnnotationType anno_type,
                                                 const Config &config,
                                                 uint64_t num_rows) {
    std::unique_ptr<Annotator> annotation;

    switch (anno_type) {
        case Config::ColumnCompressed: {
            annotation.reset(
                new annotate::ColumnCompressed<>(
                    num_rows, kNumCachedColumns, config.verbose
                )
            );
            break;
        }
        case Config::RowCompressed: {
            annotation.reset(new annotate::RowCompressed<>(num_rows, config.sparse));
            break;
        }
        case Config::BRWT: {
            annotation.reset(new annotate::BRWTCompressed<>(config.row_cache_size));
            break;
        }
        case Config::BinRelWT_sdsl: {
            annotation.reset(new annotate::BinRelWT_sdslAnnotator(config.row_cache_size));
            break;
        }
        case Config::BinRelWT: {
            annotation.reset(new annotate::BinRelWTAnnotator(config.row_cache_size));
            break;
        }
        case Config::RowFlat: {
            annotation.reset(new annotate::RowFlatAnnotator(config.row_cache_size));
            break;
        }
        case Config::RBFish: {
            annotation.reset(new annotate::RainbowfishAnnotator(config.row_cache_size));
            break;
        }
    }

    return annotation;
}

std::unique_ptr<Annotator> initialize_annotation(const std::string &filename,
                                                 const Config &config) {
    return initialize_annotation(parse_annotation_type(filename), config, 0);
}

std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(std::shared_ptr<DeBruijnGraph> graph,
                                                       const Config &config) {
    // TODO: introduce something like graph->max_node_index() to replace num_nodes() here
    auto annotation_temp = config.infbase_annotators.size()
            ? initialize_annotation(parse_annotation_type(config.infbase_annotators.at(0)), config, 0)
            : initialize_annotation(config.anno_type, config, graph->num_nodes());

    if (config.infbase_annotators.size()
            && !annotation_temp->load(config.infbase_annotators.at(0))) {
        std::cerr << "ERROR: can't load annotations for graph "
                  << config.infbase
                  << ", file corrupted" << std::endl;
        exit(1);
    }

    // load graph
    auto anno_graph = std::make_unique<AnnotatedDBG>(std::move(graph),
                                                     std::move(annotation_temp),
                                                     config.parallel);

    if (!anno_graph->check_compatibility()) {
        std::cerr << "Error: graph and annotation are not compatible."
                  << std::endl;
        exit(1);
    }

    return anno_graph;
}

std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(const Config &config) {
    return initialize_annotated_dbg(load_critical_dbg(config.infbase), config);
}

std::unique_ptr<MaskedDeBruijnGraph>
mask_graph(const AnnotatedDBG &anno_graph, Config *config) {
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());

    if (!graph.get())
        throw std::runtime_error("Masking only supported for DeBruijnGraph");

    // Remove non-present labels
    config->label_mask_in.erase(
        std::remove_if(config->label_mask_in.begin(),
                       config->label_mask_in.end(),
                       [&](const auto &label) {
                           bool exists = anno_graph.label_exists(label);
                           if (!exists && config->verbose)
                               std::cout << "Removing mask-in label " << label << std::endl;

                           return !exists;
                       }),
        config->label_mask_in.end()
    );

    config->label_mask_out.erase(
        std::remove_if(config->label_mask_out.begin(),
                       config->label_mask_out.end(),
                       [&](const auto &label) {
                           bool exists = anno_graph.label_exists(label);
                           if (!exists && config->verbose)
                               std::cout << "Removing mask-out label " << label << std::endl;

                           return !exists;
                       }),
        config->label_mask_out.end()
    );

    if (config->verbose) {
        std::cout << "Masked in:";
        for (const auto &in : config->label_mask_in) {
            std::cout << " " << in;
        }
        std::cout << std::endl;

        std::cout << "Masked out:";
        for (const auto &out : config->label_mask_out) {
            std::cout << " " << out;
        }
        std::cout << std::endl;
    }

    if (!config->filter_by_kmer) {
        return std::make_unique<MaskedDeBruijnGraph>(
            graph,
            annotated_graph_algorithm::mask_nodes_by_unitig_labels(
                anno_graph,
                config->label_mask_in,
                config->label_mask_out,
                config->label_mask_in_fraction,
                config->label_mask_out_fraction,
                config->label_other_fraction
            )
        );
    }

    return std::make_unique<MaskedDeBruijnGraph>(
        graph,
        annotated_graph_algorithm::mask_nodes_by_node_label(
            anno_graph,
            config->label_mask_in,
            config->label_mask_out,
            [config,&anno_graph](auto index,
                                 auto get_num_in_labels,
                                 auto get_num_out_labels) {
                assert(index != DeBruijnGraph::npos);

                size_t num_in_labels = get_num_in_labels();

                if (num_in_labels < config->label_mask_in_fraction
                                        * config->label_mask_in.size())
                    return false;

                size_t num_out_labels = get_num_out_labels();

                if (num_out_labels < config->label_mask_out_fraction
                                        * config->label_mask_out.size())
                    return false;

                size_t num_total_labels = anno_graph.get_labels(index).size();

                return num_total_labels - num_in_labels - num_out_labels
                            <= config->label_other_fraction * num_total_labels;
            }
        )
    );
}


template <class AnnotatorTo, class AnnotatorFrom>
void convert(std::unique_ptr<AnnotatorFrom> annotator,
             const Config &config,
             const Timer &timer) {
    if (config.verbose)
        std::cout << "Converting to " << Config::annotype_to_string(config.anno_type)
                  << " annotator...\t" << std::flush;

    auto target_annotator = annotate::convert<AnnotatorTo>(std::move(*annotator));
    annotator.reset();
    if (config.verbose)
        std::cout << timer.elapsed() << "sec" << std::endl;

    if (config.verbose)
        std::cout << "Serializing to " << config.outfbase
                  << "...\t" << std::flush;
    target_annotator->serialize(config.outfbase);
    if (config.verbose)
        std::cout << timer.elapsed() << "sec" << std::endl;
}


void set_aligner_parameters(const DeBruijnGraph &graph, Config &config) {
    // fix seed length bounds
    if (!config.alignment_min_seed_length || config.alignment_seed_unimems)
        config.alignment_min_seed_length = graph.get_k();

    if (config.alignment_max_seed_length == std::numeric_limits<size_t>::max()
            && !config.alignment_seed_unimems)
        config.alignment_max_seed_length = graph.get_k();

    if (config.verbose) {
        std::cout << "Alignment settings:" << "\n"
                  << "\t Seeding: " << (config.alignment_seed_unimems ? "unimems" : "nodes") << "\n"
                  << "\t Alignments to report: " << config.alignment_num_alternative_paths << "\n"
                  << "\t Priority queue size: " << config.alignment_queue_size << "\n"
                  << "\t Min seed length: " << config.alignment_min_seed_length << "\n"
                  << "\t Max seed length: " << config.alignment_max_seed_length << "\n"
                  << "\t Max num seeds per locus: " << config.alignment_max_num_seeds_per_locus << "\n"
                  << "\t Scoring matrix: " << (config.alignment_edit_distance ? "unit costs" : "matrix") << "\n"
                  << "\t Gap opening penalty: " << int64_t(config.alignment_gap_opening_penalty) << "\n"
                  << "\t Gap extension penalty: " << int64_t(config.alignment_gap_extension_penalty) << "\n"
                  << "\t Min DP table cell score: " << int64_t(config.alignment_min_cell_score) << "\n"
                  << "\t Min alignment score: " << config.alignment_min_path_score << std::endl;

        if (!config.alignment_edit_distance)
            std::cout << "\t Match score: " << int64_t(config.alignment_match_score) << "\n"
                      << "\t (DNA) Transition score: " << int64_t(config.alignment_mm_transition) << "\n"
                      << "\t (DNA) Transversion score: " << int64_t(config.alignment_mm_transversion) << "\n";

        std::cout << std::endl;
    }
}

std::unique_ptr<IDBGAligner> build_aligner(const DeBruijnGraph &graph, Config &config) {
    set_aligner_parameters(graph, config);

    // TODO: fix this when alphabets are no longer set at compile time
    #if _PROTEIN_GRAPH
        const auto *alphabet = alphabets::kAlphabetProtein;
        const auto *alphabet_encoding = alphabets::kCharToProtein;
    #elif _DNA_CASE_SENSITIVE_GRAPH
        const auto *alphabet = alphabets::kAlphabetDNA;
        const auto *alphabet_encoding = alphabets::kCharToDNA;
    #elif _DNA5_GRAPH
        const auto *alphabet = alphabets::kAlphabetDNA;
        const auto *alphabet_encoding = alphabets::kCharToDNA;
    #elif _DNA_GRAPH
        const auto *alphabet = alphabets::kAlphabetDNA;
        const auto *alphabet_encoding = alphabets::kCharToDNA;
    #else
        static_assert(false,
            "Define an alphabet: either "
            "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
        );
    #endif

    Cigar::initialize_opt_table(alphabet, alphabet_encoding);

    if (config.alignment_seed_unimems) {
        return std::make_unique<DBGAligner<UniMEMSeeder<>>>(graph, DBGAlignerConfig(config));

    } else if (config.alignment_min_seed_length < graph.get_k()) {
        if (!dynamic_cast<const DBGSuccinct*>(&graph)) {
            std::cerr << "ERROR: SuffixSeeder can be used only with succinct graph representation"
                      << std::endl;
            exit(1);
        }

        // Use the seeder that seeds to node suffixes
        return std::make_unique<DBGAligner<SuffixSeeder<>>>(graph, DBGAlignerConfig(config));

    } else {
        return std::make_unique<DBGAligner<>>(graph, DBGAlignerConfig(config));
    }
}

void map_sequences_in_file(const std::string &file,
                           const DeBruijnGraph &graph,
                           std::shared_ptr<DBGSuccinct> dbg,
                           const Config &config,
                           const Timer &timer,
                           ThreadPool *thread_pool = nullptr,
                           std::mutex *print_mutex = nullptr) {
    // TODO: multithreaded
    std::ignore = std::tie(thread_pool, print_mutex);

    Timer data_reading_timer;

    std::function<bool(std::string)> map_kmers;
    if (dbg) {
        map_kmers = [&](std::string sequence) {
            return dbg->get_boss().find(sequence,
                                        config.discovery_fraction,
                                        config.kmer_mapping_mode);
        };
    } else {
        map_kmers = [&](std::string sequence) {
            return graph.find(sequence, config.discovery_fraction);
        };
    }

    read_fasta_file_critical(file, [&](kseq_t *read_stream) {
        if (config.verbose)
            std::cout << "Sequence: " << read_stream->seq.s << "\n";

        if (config.query_presence
                && config.alignment_length == graph.get_k()) {
            if (config.filter_present) {
                if (map_kmers(std::string(read_stream->seq.s)))
                    std::cout << ">" << read_stream->name.s << "\n"
                                     << read_stream->seq.s << "\n";
            } else {
                std::cout << map_kmers(std::string(read_stream->seq.s)) << "\n";
            }
            return;
        }

        assert(config.alignment_length <= graph.get_k());

        std::vector<DeBruijnGraph::node_index> graphindices;
        if (config.alignment_length == graph.get_k()) {
            graph.map_to_nodes(read_stream->seq.s,
                               [&](const auto &node) {
                                   graphindices.emplace_back(node);
                               });
        } else if (config.query_presence || config.count_kmers) {
            // TODO: make more efficient
            for (size_t i = 0; i + graph.get_k() <= read_stream->seq.l; ++i) {
                dbg->call_nodes_with_suffix(
                    read_stream->seq.s + i,
                    read_stream->seq.s + i + config.alignment_length,
                    [&](auto node, auto) {
                        if (graphindices.empty())
                            graphindices.emplace_back(node);
                    },
                    config.alignment_length
                );
            }
        }

        size_t num_discovered = std::count_if(graphindices.begin(),
                                              graphindices.end(),
                                              [](const auto &x) { return x > 0; });

        const size_t num_kmers = graphindices.size();

        if (config.query_presence) {
            const size_t min_kmers_discovered =
                num_kmers - num_kmers * (1 - config.discovery_fraction);
            if (config.filter_present) {
                if (num_discovered >= min_kmers_discovered)
                    std::cout << ">" << read_stream->name.s << "\n"
                                     << read_stream->seq.s << "\n";
            } else {
                std::cout << (num_discovered >= min_kmers_discovered) << "\n";
            }
            return;
        }

        if (config.count_kmers) {
            std::cout << "Kmers matched (discovered/total): "
                      << num_discovered << "/"
                      << num_kmers << "\n";
            return;
        }

        if (config.alignment_length == graph.get_k()) {
            for (size_t i = 0; i < graphindices.size(); ++i) {
                assert(i + config.alignment_length <= read_stream->seq.l);
                std::cout << std::string(read_stream->seq.s + i, config.alignment_length)
                          << ": " << graphindices[i] << "\n";
            }
        } else {
            // map input subsequences to multiple nodes
            for (size_t i = 0; i + graph.get_k() <= read_stream->seq.l; ++i) {
                // TODO: make more efficient
                std::string subseq(read_stream->seq.s + i,
                                   read_stream->seq.s + i + config.alignment_length);

                dbg->call_nodes_with_suffix(subseq.begin(),
                                            subseq.end(),
                                            [&](auto node, auto) {
                                                std::cout << subseq << ": "
                                                          << node
                                                          << "\n";
                                            },
                                            config.alignment_length);
            }
        }

    }, config.forward_and_reverse);

    if (config.verbose) {
        std::cout << "File processed in "
                  << data_reading_timer.elapsed()
                  << "sec, current mem usage: "
                  << (get_curr_RSS() >> 20) << " MiB"
                  << ", total time: " << timer.elapsed()
                  << "sec" << std::endl;
    }
}

typedef std::function<void(const std::string&)> SequenceCallback;

std::unique_ptr<AnnotatedDBG>
construct_query_graph(const AnnotatedDBG &anno_graph,
                      std::function<void(SequenceCallback)> call_sequences,
                      double discovery_fraction,
                      size_t num_threads) {
    const auto *full_dbg = dynamic_cast<const DeBruijnGraph*>(&anno_graph.get_graph());
    if (!full_dbg)
        throw std::runtime_error("Error: batch queries are supported only for de Bruijn graphs");

    const auto &full_annotation = anno_graph.get_annotation();

    Timer timer;

    // construct graph storing all k-mers in query
    std::shared_ptr<DeBruijnGraph> graph
        = std::make_shared<DBGHashOrdered>(full_dbg->get_k(), false);

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(full_dbg);
    const auto *bloom_filter = dbg_succ ? dbg_succ->get_bloom_filter() : nullptr;

    call_sequences([&](const std::string &sequence) {
        if (!bloom_filter || sdsl::util::cnt_one_bits(bloom_filter->check_kmer_presence(sequence)))
            graph->add_sequence(sequence);
    });

    if (utils::get_verbose()) {
        std::cout << "Query graph --- k-mers indexed: "
                  << timer.elapsed() << " sec" << std::endl;
        timer.reset();
    }

    exit(1);

    // pull contigs from query graph
    std::vector<std::pair<std::string, std::vector<DeBruijnGraph::node_index>>> contigs;
    graph->call_sequences(
        [&](const std::string &contig, const auto &path) { contigs.emplace_back(contig, path); },
        full_dbg->is_canonical_mode()
    );

    if (utils::get_verbose()) {
        std::cout << "Query graph --- contigs extracted: "
                  << timer.elapsed() << " sec" << std::endl;
        timer.reset();
    }

    if (full_dbg->is_canonical_mode()) {
        // construct graph storing all distinct k-mers in query
        graph = std::make_shared<DBGHashOrdered>(full_dbg->get_k(), true);

        for (const auto &pair : contigs) {
            graph->add_sequence(pair.first);
        }

        if (utils::get_verbose()) {
            std::cout << "Query graph --- reindexed k-mers in canonical mode: "
                      << timer.elapsed() << " sec" << std::endl;
            timer.reset();
        }
    }

    // map contigs onto the full graph
    auto index_in_full_graph
        = std::make_shared<std::vector<uint64_t>>(graph->num_nodes() + 1, 0);

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 10)
    for (size_t i = 0; i < contigs.size(); ++i) {

        auto contig = std::move(contigs[i].first);
        auto path = std::move(contigs[i].second);

        if (graph->is_canonical_mode()) {
            size_t j = 0;
            graph->map_to_nodes(contig, [&](auto node) { path[j++] = node; });
            assert(j == path.size());
        }

        size_t j = 0;

        full_dbg->map_to_nodes(contig,
            [&](auto node_in_full) { (*index_in_full_graph)[path[j++]] = node_in_full; }
        );

        assert(j == path.size());
    }

    if (utils::get_verbose()) {
        std::cout << "Query graph --- contigs mapped to graph: "
                  << timer.elapsed() << " sec" << std::endl;
        timer.reset();
    }

    contigs.clear();

    assert(!(*index_in_full_graph)[0]);

    if (discovery_fraction > 0) {
        sdsl::bit_vector mask(graph->num_nodes() + 1, false);

        call_sequences([&](const std::string &sequence) {
            const size_t num_kmers = sequence.length() - graph->get_k() + 1;
            const size_t max_kmers_missing = num_kmers * (1 - discovery_fraction);
            const size_t min_kmers_discovered = num_kmers - max_kmers_missing;
            size_t num_kmers_discovered = 0;
            size_t num_kmers_missing = 0;

            std::vector<DeBruijnGraph::node_index> nodes;
            nodes.reserve(num_kmers);

            graph->map_to_nodes(sequence,
                [&](auto node) {
                    if ((*index_in_full_graph)[node]) {
                        num_kmers_discovered++;
                        nodes.push_back(node);
                    } else {
                        num_kmers_missing++;
                    }
                },
                [&]() { return num_kmers_missing > max_kmers_missing
                                || num_kmers_discovered >= min_kmers_discovered; }
            );

            if (num_kmers_missing <= max_kmers_missing) {
                for (auto node : nodes) { mask[node] = true; }
            }
        });

        // correcting the mask
        call_zeros(mask, [&](auto i) { (*index_in_full_graph)[i] = 0; });

        if (utils::get_verbose()) {
            std::cout << "Query graph --- reduced k-mer dictionary: "
                      << timer.elapsed() << " sec" << std::endl;
            timer.reset();
        }
    }

    assert(index_in_full_graph.get());

    std::vector<std::pair<uint64_t, uint64_t>> from_full_to_query;
    from_full_to_query.reserve(index_in_full_graph->size());

    for (uint64_t node = 0; node < index_in_full_graph->size(); ++node) {
        if ((*index_in_full_graph)[node]) {
            from_full_to_query.emplace_back(
                AnnotatedDBG::graph_to_anno_index((*index_in_full_graph)[node]),
                AnnotatedDBG::graph_to_anno_index(node)
            );
        }
    }

    ips4o::parallel::sort(from_full_to_query.begin(), from_full_to_query.end(),
        [](const auto &first, const auto &second) { return first.first < second.first; },
        num_threads
    );

    // initialize fast query annotation
    // copy annotations from the full graph to the query graph
    auto annotation = std::make_unique<annotate::RowCompressed<>>(
        graph->num_nodes(),
        full_annotation.get_label_encoder().get_labels(),
        [&](annotate::RowCompressed<>::CallRow call_row) {

            #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
            for (uint64_t batch_begin = 0;
                                batch_begin < from_full_to_query.size();
                                                batch_begin += kRowBatchSize) {

                const uint64_t batch_end
                    = std::min(batch_begin + kRowBatchSize,
                               static_cast<uint64_t>(from_full_to_query.size()));

                std::vector<uint64_t> row_indexes;
                row_indexes.reserve(batch_end - batch_begin);

                for (uint64_t i = batch_begin; i < batch_end; ++i) {
                    assert(from_full_to_query[i].first < full_annotation.num_objects());

                    row_indexes.push_back(from_full_to_query[i].first);
                }

                auto rows = full_annotation.get_label_codes(row_indexes);

                assert(rows.size() == batch_end - batch_begin);

                for (uint64_t i = batch_begin; i < batch_end; ++i) {
                    call_row(from_full_to_query[i].second,
                             std::move(rows[i - batch_begin]));
                }
            }
        }
    );

    if (utils::get_verbose()) {
        std::cout << "Query graph --- constructed query annotation: "
                  << timer.elapsed() << " sec" << std::endl;
        timer.reset();
    }

    graph = std::make_shared<MaskedDeBruijnGraph>(graph,
        [=](auto i) -> bool { return (*index_in_full_graph)[i]; }
    );

    // build annotated graph from the query graph and copied annotations
    return std::make_unique<AnnotatedDBG>(graph, std::move(annotation));
}


void print_boss_stats(const BOSS &boss_graph,
                      bool count_dummy = false,
                      size_t num_threads = 0,
                      bool verbose = false) {
    std::cout << "====================== BOSS STATS ======================" << std::endl;
    std::cout << "k: " << boss_graph.get_k() + 1 << std::endl;
    std::cout << "nodes (k-1): " << boss_graph.num_nodes() << std::endl;
    std::cout << "edges ( k ): " << boss_graph.num_edges() << std::endl;
    std::cout << "state: " << Config::state_to_string(boss_graph.get_state()) << std::endl;

    assert(boss_graph.rank_W(boss_graph.num_edges(), boss_graph.alph_size) == 0);
    std::cout << "W stats: {'" << boss_graph.decode(0) << "': "
              << boss_graph.rank_W(boss_graph.num_edges(), 0);
    for (int i = 1; i < boss_graph.alph_size; ++i) {
        std::cout << ", '" << boss_graph.decode(i) << "': "
                  << boss_graph.rank_W(boss_graph.num_edges(), i)
                        + boss_graph.rank_W(boss_graph.num_edges(), i + boss_graph.alph_size);
    }
    std::cout << "}" << std::endl;

    assert(boss_graph.get_F(0) == 0);
    std::cout << "F stats: {'";
    for (int i = 1; i < boss_graph.alph_size; ++i) {
        std::cout << boss_graph.decode(i - 1) << "': "
                  << boss_graph.get_F(i) - boss_graph.get_F(i - 1)
                  << ", '";
    }
    std::cout << boss_graph.decode(boss_graph.alph_size - 1) << "': "
              << boss_graph.num_edges() - boss_graph.get_F(boss_graph.alph_size - 1)
              << "}" << std::endl;

    if (count_dummy) {
        std::cout << "dummy source edges: "
                  << boss_graph.mark_source_dummy_edges(NULL, num_threads, verbose)
                  << std::endl;
        std::cout << "dummy sink edges: "
                  << boss_graph.mark_sink_dummy_edges()
                  << std::endl;
    }
    std::cout << "========================================================" << std::endl;
}

void print_stats(const DeBruijnGraph &graph) {
    std::cout << "====================== GRAPH STATS =====================" << std::endl;
    std::cout << "k: " << graph.get_k() << std::endl;
    std::cout << "nodes (k): " << graph.num_nodes() << std::endl;
    std::cout << "canonical mode: " << (graph.is_canonical_mode() ? "yes" : "no") << std::endl;

    if (auto weights = graph.get_extension<NodeWeights>()) {
        double sum_weights = 0;
        uint64_t num_non_zero_weights = 0;
        if (const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph)) {
            // In DBGSuccinct some of the nodes may be masked out
            // TODO: Fix this by using non-contiguous indexing in graph
            //       so that mask of dummy edges does not change indexes.
            for (uint64_t i = 1; i <= dbg_succ->get_boss().num_edges(); ++i) {
                if (uint64_t weight = (*weights)[i]) {
                    sum_weights += weight;
                    num_non_zero_weights++;
                }
            }
        } else {
            if (!weights->is_compatible(graph)) {
                std::cerr << "ERROR: node weights are not compatible with graph" << std::endl;
                exit(1);
            }
            graph.call_nodes([&](auto i) {
                if (uint64_t weight = (*weights)[i]) {
                    sum_weights += weight;
                    num_non_zero_weights++;
                }
            });
        }
        std::cout << "nnz weights: " << num_non_zero_weights << std::endl;
        std::cout << "avg weight: " << static_cast<double>(sum_weights) / num_non_zero_weights << std::endl;

        if (utils::get_verbose()) {
            if (const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(&graph)) {
                // In DBGSuccinct some of the nodes may be masked out
                // TODO: Fix this by using non-contiguous indexing in graph
                //       so that mask of dummy edges does not change indexes.
                for (uint64_t i = 1; i <= dbg_succ->get_boss().num_edges(); ++i) {
                    if (uint64_t weight = (*weights)[i])
                        std::cout << weight << " ";
                }
            } else {
                graph.call_nodes([&](auto i) { std::cout << (*weights)[i] << " "; });
            }
            std::cout << std::endl;
        }
    }

    std::cout << "========================================================" << std::endl;
}

void print_stats(const Annotator &annotation) {
    std::cout << "=================== ANNOTATION STATS ===================" << std::endl;
    std::cout << "labels:  " << annotation.num_labels() << std::endl;
    std::cout << "objects: " << annotation.num_objects() << std::endl;
    std::cout << "density: " << static_cast<double>(annotation.num_relations())
                                    / annotation.num_objects()
                                    / annotation.num_labels() << std::endl;
    std::cout << "representation: ";

    if (dynamic_cast<const annotate::ColumnCompressed<std::string> *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::ColumnCompressed) << std::endl;

    } else if (dynamic_cast<const annotate::RowCompressed<std::string> *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::RowCompressed) << std::endl;

    } else if (dynamic_cast<const annotate::BRWTCompressed<std::string> *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::BRWT) << std::endl;
        const auto &brwt = dynamic_cast<const annotate::BRWTCompressed<std::string> &>(annotation).data();
        std::cout << "=================== Multi-BRWT STATS ===================" << std::endl;
        std::cout << "num nodes: " << brwt.num_nodes() << std::endl;
        std::cout << "avg arity: " << brwt.avg_arity() << std::endl;
        std::cout << "shrinkage: " << brwt.shrinking_rate() << std::endl;
        if (utils::get_verbose()) {
            std::cout << "==================== Multi-BRWT TREE ===================" << std::endl;
            brwt.print_tree_structure(std::cout);
        }

    } else if (dynamic_cast<const annotate::BinRelWT_sdslAnnotator *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::BinRelWT_sdsl) << std::endl;

    } else if (dynamic_cast<const annotate::BinRelWTAnnotator *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::BinRelWT) << std::endl;

    } else if (dynamic_cast<const annotate::RowFlatAnnotator *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::RowFlat) << std::endl;

    } else if (dynamic_cast<const annotate::RainbowfishAnnotator *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::RBFish) << std::endl;

    } else {
        assert(false);
        throw std::runtime_error("Unknown annotator");
    }
    std::cout << "========================================================" << std::endl;
}

template <class Callback, class CountedKmer, class Loop>
void parse_sequences(const std::vector<std::string> &files,
                     const Config &config,
                     const Timer &timer,
                     Callback call_sequence,
                     CountedKmer call_kmer,
                     Loop call_sequences) {
    // iterate over input files
    for (const auto &file : files) {
        if (config.verbose) {
            std::cout << std::endl << "Parsing " << file << std::endl;
        }

        Timer data_reading_timer;

        if (utils::get_filetype(file) == "VCF") {
            read_vcf_file_critical(file,
                                   config.refpath,
                                   config.k,
                                   [&](std::string&& sequence) {
                                       call_sequence(std::move(sequence));
                                   },
                                   config.forward_and_reverse);

        } else if (utils::get_filetype(file) == "KMC") {
            bool warning_different_k = false;

            auto min_count = config.min_count;
            auto max_count = config.max_count;

            if (config.min_count_quantile > 0 || config.max_count_quantile < 1) {
                std::unordered_map<uint64_t, uint64_t> count_hist;
                kmc::read_kmers(
                    file,
                    [&](std::string&&, uint32_t count) { count_hist[count]++; },
                    !config.canonical
                );

                if (count_hist.size()) {
                    std::vector<std::pair<uint64_t, uint64_t>> count_hist_v(count_hist.begin(),
                                                                            count_hist.end());

                    ips4o::parallel::sort(count_hist_v.begin(), count_hist_v.end(),
                        [](const auto &first, const auto &second) {
                            return first.first < second.first;
                        },
                        config.parallel
                    );

                    if (config.min_count_quantile > 0)
                        min_count = utils::get_quantile(count_hist_v, config.min_count_quantile);
                    if (config.max_count_quantile < 1)
                        max_count = utils::get_quantile(count_hist_v, config.max_count_quantile);

                    std::cout << "Used k-mer count thresholds:\n"
                              << "min (including): " << min_count << "\n"
                              << "max (excluding): " << max_count << std::endl;
                }
            }

            kmc::read_kmers(
                file,
                [&](std::string&& sequence, uint32_t count) {
                    if (!warning_different_k && sequence.size() != config.k) {
                        std::cerr << "Warning: k-mers parsed from KMC database "
                                  << file << " have length " << sequence.size()
                                  << " but graph is constructed for k=" << config.k
                                  << std::endl;
                        warning_different_k = true;
                    }
                    call_kmer(std::move(sequence), count);
                },
                !config.canonical,
                min_count,
                max_count
            );

        } else if (utils::get_filetype(file) == "FASTA"
                    || utils::get_filetype(file) == "FASTQ") {
            if (files.size() >= config.parallel) {
                auto forward_and_reverse = config.forward_and_reverse;

                // capture all required values by copying to be able
                // to run task from other threads
                call_sequences([=](auto callback) {
                    read_fasta_file_critical(file, [=](kseq_t *read_stream) {
                        // add read to the graph constructor as a callback
                        callback(read_stream->seq.s);
                    }, forward_and_reverse);
                });
            } else {
                read_fasta_file_critical(file, [&](kseq_t *read_stream) {
                    // add read to the graph constructor as a callback
                    call_sequence(read_stream->seq.s);
                }, config.forward_and_reverse);
            }
        } else {
            std::cerr << "ERROR: Filetype unknown for file "
                      << file << std::endl;
            exit(1);
        }

        if (config.verbose) {
            std::cout << "Finished extracting sequences from file " << file
                      << " in " << timer.elapsed() << "sec" << std::endl;
        }
        if (config.verbose) {
            std::cout << "File processed in "
                      << data_reading_timer.elapsed()
                      << "sec, current mem usage: "
                      << (get_curr_RSS() >> 20) << " MiB"
                      << ", total time: " << timer.elapsed()
                      << "sec" << std::endl;
        }
    }

    if (config.verbose) {
        std::cout << std::endl;
    }
}

std::string form_client_reply(const std::string &received_message,
                              const AnnotatedDBG &anno_graph,
                              const Config &config,
                              IDBGAligner *aligner = nullptr) {
    try {
        Json::Value json;

        {
            Json::CharReaderBuilder rbuilder;
            std::unique_ptr<Json::CharReader> reader { rbuilder.newCharReader() };
            std::string errors;

            if (!reader->parse(received_message.data(),
                               received_message.data() + received_message.size(),
                               &json,
                               &errors)) {
                std::cerr << "Error: bad json file:\n" << errors << std::endl;
                //TODO: send error message back in a json file
                throw std::domain_error("bad json received");
            }
        }

        const auto &fasta = json["FASTA"];
        const auto &seq = json["SEQ"];

        // discovery_fraction a proxy of 1 - %similarity
        auto discovery_fraction = json.get("discovery_fraction",
                                           config.discovery_fraction).asDouble();
        auto count_labels = json.get("count_labels", config.count_labels).asBool();
        auto num_top_labels = json.get("num_labels", config.num_top_labels).asInt();

        std::ostringstream oss;

        // query callback shared by FASTA and sequence modes
        auto execute_server_query = [&](const std::string &name,
                                        const std::string &sequence) {
            execute_query(name,
                          sequence,
                          count_labels,
                          config.suppress_unlabeled,
                          num_top_labels,
                          discovery_fraction,
                          config.anno_labels_delimiter,
                          anno_graph,
                          oss,
                          aligner);
        };

        if (!seq.isNull()) {
            // input is plain sequence
            execute_server_query(seq.asString(), seq.asString());
        } else if (!fasta.isNull()) {
            // input is a FASTA sequence
            read_fasta_from_string(
                fasta.asString(),
                [&](kseq_t *read_stream) {
                    execute_server_query(read_stream->name.s,
                                         read_stream->seq.s);
                }
            );
        } else {
            std::cerr << "Error: no input sequences received from client" << std::endl;
            // TODO: no input sequences -> form an error message for the client
            throw std::domain_error("No input sequences");
        }

        return oss.str();

    } catch (const Json::LogicError &e) {
        std::cerr << "Error: bad json file: " << e.what() << std::endl;
        //TODO: send errors in a json file
        throw;
    } catch (const std::exception &e) {
        std::cerr << "Error: processing request error: " << e.what() << std::endl;
        //TODO: send errors in a json file
        throw;
    } catch (...) {
        std::cerr << "Error: processing request error" << std::endl;
        //TODO: send errors in a json file
        throw;
    }
}


int main(int argc, const char *argv[]) {
    auto config = std::make_unique<Config>(argc, argv);

    if (config->verbose) {
        std::cout << "#############################\n"
                  << "### Welcome to MetaGraph! ###\n"
                  << "#############################\n" << std::endl;
    }

    const auto &files = config->fname;

    switch (config->identity) {
        case Config::EXPERIMENT: {
            break;
        }
        case Config::BUILD: {
            std::unique_ptr<DeBruijnGraph> graph;

            if (config->verbose)
                std::cout << "Build De Bruijn Graph with k-mer size k="
                          << config->k << std::endl;

            Timer timer;

            if (config->canonical)
                config->forward_and_reverse = false;

            if (config->complete) {
                if (config->graph_type != Config::GraphType::BITMAP) {
                    std::cerr << "Error: Only bitmap-graph can be built"
                              << " in complete mode" << std::endl;
                    exit(1);
                }

                graph.reset(new DBGBitmap(config->k, config->canonical));

            } else if (config->graph_type == Config::GraphType::SUCCINCT && !config->dynamic) {
                auto boss_graph = std::make_unique<BOSS>(config->k - 1);

                if (config->verbose) {
                    std::cout << "Start reading data and extracting k-mers" << std::endl;
                }
                //enumerate all suffixes
                assert(boss_graph->alph_size > 1);
                std::vector<std::string> suffixes;
                if (config->suffix.size()) {
                    suffixes = { config->suffix };
                } else {
                    suffixes = KmerExtractorBOSS::generate_suffixes(config->suffix_len);
                }

                BOSS::Chunk graph_data(KmerExtractorBOSS::alphabet.size(),
                                       boss_graph->get_k(),
                                       config->canonical);

                //one pass per suffix
                for (const std::string &suffix : suffixes) {
                    timer.reset();

                    if (suffix.size() > 0 || suffixes.size() > 1) {
                        std::cout << "\nSuffix: " << suffix << std::endl;
                    }

                    auto constructor = IBOSSChunkConstructor::initialize(
                        boss_graph->get_k(),
                        config->canonical,
                        config->count_kmers,
                        suffix,
                        config->parallel,
                        static_cast<uint64_t>(config->memory_available) << 30,
                        config->verbose
                    );

                    parse_sequences(files, *config, timer,
                        [&](std::string&& read) { constructor->add_sequence(std::move(read)); },
                        [&](std::string&& kmer, uint32_t count) { constructor->add_sequence(std::move(kmer), count); },
                        [&](const auto &loop) { constructor->add_sequences(loop); }
                    );

                    auto next_block = constructor->build_chunk();
                    if (config->verbose) {
                        std::cout << "Graph chunk with " << next_block->size()
                                  << " k-mers was built in "
                                  << timer.elapsed() << "sec" << std::endl;
                    }

                    if (config->outfbase.size() && config->suffix.size()) {
                        std::cout << "Serialize the graph chunk for suffix '"
                                  << suffix << "'...\t" << std::flush;
                        timer.reset();
                        next_block->serialize(config->outfbase + "." + suffix);
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    }

                    if (config->suffix.size())
                        return 0;

                    graph_data.extend(*next_block);
                    delete next_block;
                }

                if (config->count_kmers) {
                    sdsl::int_vector<> kmer_counts(0, 0, kBitsPerCount);
                    graph_data.initialize_boss(boss_graph.get(), &kmer_counts);
                    graph.reset(new DBGSuccinct(boss_graph.release(), config->canonical));
                    graph->add_extension(std::make_shared<NodeWeights>(std::move(kmer_counts)));
                    assert(graph->get_extension<NodeWeights>()->is_compatible(*graph));
                } else {
                    graph_data.initialize_boss(boss_graph.get());
                    graph.reset(new DBGSuccinct(boss_graph.release(), config->canonical));
                }

            } else if (config->graph_type == Config::GraphType::BITMAP && !config->dynamic) {

                if (!config->outfbase.size()) {
                    std::cerr << "Error: No output file provided" << std::endl;
                    exit(1);
                }

                if (config->verbose) {
                    std::cout << "Start reading data and extracting k-mers" << std::endl;
                }
                // enumerate all suffixes
                std::vector<std::string> suffixes;
                if (config->suffix.size()) {
                    suffixes = { config->suffix };
                } else {
                    suffixes = KmerExtractor2Bit().generate_suffixes(config->suffix_len);
                }

                std::unique_ptr<DBGBitmapConstructor> constructor;
                std::vector<std::string> chunk_filenames;

                //one pass per suffix
                for (const std::string &suffix : suffixes) {
                    timer.reset();

                    if (config->verbose && (suffix.size() > 0 || suffixes.size() > 1)) {
                        std::cout << "\nSuffix: " << suffix << std::endl;
                    }

                    constructor.reset(
                        new DBGBitmapConstructor(
                            config->k,
                            config->canonical,
                            config->count_kmers ? kBitsPerCount : 0,
                            suffix,
                            config->parallel,
                            static_cast<uint64_t>(config->memory_available) << 30,
                            config->verbose
                        )
                    );

                    parse_sequences(files, *config, timer,
                        [&](std::string&& read) { constructor->add_sequence(std::move(read)); },
                        [&](std::string&& kmer, uint32_t count) { constructor->add_sequence(std::move(kmer), count); },
                        [&](const auto &loop) { constructor->add_sequences(loop); }
                    );

                    if (!suffix.size()) {
                        assert(suffixes.size() == 1);

                        auto *bitmap_graph = new DBGBitmap(config->k);
                        constructor->build_graph(bitmap_graph);
                        graph.reset(bitmap_graph);

                    } else {
                        std::unique_ptr<DBGBitmap::Chunk> chunk { constructor->build_chunk() };
                        if (config->verbose) {
                            std::cout << "Graph chunk with " << chunk->num_set_bits()
                                      << " k-mers was built in "
                                      << timer.elapsed() << "sec" << std::endl;

                            std::cout << "Serialize the graph chunk for suffix '"
                                      << suffix << "'...\t" << std::flush;
                        }

                        chunk_filenames.push_back(
                            utils::join_strings({ config->outfbase, suffix }, ".")
                                + DBGBitmap::kChunkFileExtension
                        );
                        std::ofstream out(chunk_filenames.back(), std::ios::binary);
                        chunk->serialize(out);
                        if (config->verbose)
                            std::cout << timer.elapsed() << "sec" << std::endl;
                    }

                    // only one chunk had to be constructed
                    if (config->suffix.size())
                        return 0;
                }

                if (suffixes.size() > 1) {
                    assert(chunk_filenames.size());
                    timer.reset();
                    graph.reset(constructor->build_graph_from_chunks(chunk_filenames,
                                                                     config->canonical,
                                                                     config->verbose));
                }

            } else {
                //slower method
                switch (config->graph_type) {

                    case Config::GraphType::SUCCINCT:
                        graph.reset(new DBGSuccinct(config->k, config->canonical));
                        break;

                    case Config::GraphType::HASH:
                        graph.reset(new DBGHashOrdered(config->k, config->canonical));
                        break;

                    case Config::GraphType::HASH_PACKED:
                        graph.reset(new DBGHashOrdered(config->k, config->canonical, true));
                        break;

                    case Config::GraphType::HASH_STR:
                        if (config->canonical) {
                            std::cerr << "Warning: string hash-based de Bruijn graph"
                                      << " does not support canonical mode."
                                      << " Normal mode will be used instead." << std::endl;
                        }
                        // TODO: implement canonical mode
                        graph.reset(new DBGHashString(config->k/*, config->canonical*/));
                        break;

                    case Config::GraphType::BITMAP:
                        std::cerr << "Error: Bitmap-graph construction"
                                  << " in dynamic regime is not supported" << std::endl;
                        exit(1);

                    case Config::GraphType::INVALID:
                        assert(false);
                }
                assert(graph);

                parse_sequences(files, *config, timer,
                    [&graph](std::string&& seq) {
                        graph->add_sequence(std::move(seq));
                    },
                    [&graph](std::string&& kmer, uint32_t /*count*/) {
                        graph->add_sequence(std::move(kmer));
                    },
                    [&graph](const auto &loop) {
                        loop([&graph](const char *seq) { graph->add_sequence(seq); });
                    }
                );

                if (config->count_kmers) {
                    graph->add_extension(std::make_shared<NodeWeights>(graph->num_nodes() + 1, kBitsPerCount));
                    auto node_weights = graph->get_extension<NodeWeights>();
                    assert(node_weights->is_compatible(*graph));

                    if (graph->is_canonical_mode())
                        config->forward_and_reverse = true;

                    parse_sequences(files, *config, timer,
                        [&graph,&node_weights](std::string&& seq) {
                            graph->map_to_nodes_sequentially(seq.begin(), seq.end(),
                                [&](auto node) { node_weights->add_weight(node, 1); }
                            );
                        },
                        [&graph,&node_weights](std::string&& kmer, uint32_t count) {
                            node_weights->add_weight(graph->kmer_to_node(kmer), count);
                        },
                        [&graph,&node_weights](const auto &loop) {
                            loop([&graph,&node_weights](const char *seq) {
                                std::string seq_str(seq);
                                graph->map_to_nodes_sequentially(seq_str.begin(), seq_str.end(),
                                    [&](auto node) { node_weights->add_weight(node, 1); }
                                );
                            });
                        }
                    );
                }
            }

            if (config->verbose)
                std::cout << "Graph construction finished in "
                          << timer.elapsed() << "sec" << std::endl;

            if (!config->outfbase.empty()) {
                if (dynamic_cast<DBGSuccinct*>(graph.get()) && config->mark_dummy_kmers) {
                    if (config->verbose)
                        std::cout << "Detecting all dummy k-mers..." << std::flush;

                    timer.reset();
                    dynamic_cast<DBGSuccinct&>(*graph).mask_dummy_kmers(config->parallel, false);

                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                }

                graph->serialize(config->outfbase);
                graph->serialize_extensions(config->outfbase);
            }

            return 0;
        }
        case Config::EXTEND: {
            assert(config->infbase_annotators.size() <= 1);

            Timer timer;

            // load graph
            auto graph = load_critical_dbg(config->infbase);

            auto node_weights = graph->load_extension<NodeWeights>(config->infbase);
            // TODO: fix extension of DBGSuccinct with k-mer counts
            //       DBGSuccinct with mask of dummy edges initialized uses
            //       contiguous indexes that are not compatible with node weights,
            //       which are indexed by the rows of the BOSS table.
            //       This can be fixed by using the same indexes in all cases
            //       (non-contiguous indexing)
            if (!node_weights->is_compatible(*graph)) {
                std::cerr << "Error: node weights are not compatible with graph "
                          << config->infbase
                          << " and will not be updated." << std::endl;
                node_weights.reset();
            }

            if (config->verbose) {
                std::cout << "De Bruijn graph with k-mer size k="
                          << graph->get_k() << " has been loaded in "
                          << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();

            if (dynamic_cast<DBGSuccinct*>(graph.get())) {
                auto &succinct_graph = dynamic_cast<DBGSuccinct&>(*graph);

                if (succinct_graph.get_state() != Config::DYN) {
                    if (config->verbose)
                        std::cout << "Switching state of succinct graph to dynamic..." << std::flush;

                    succinct_graph.switch_state(Config::DYN);

                    if (config->verbose)
                        std::cout << "\tdone in " << timer.elapsed() << "sec" << std::endl;
                }
            }

            std::unique_ptr<bit_vector_dyn> inserted_edges;
            if (config->infbase_annotators.size() || node_weights)
                inserted_edges.reset(new bit_vector_dyn(graph->num_nodes() + 1, 0));

            timer.reset();

            if (config->verbose)
                std::cout << "Start graph extension" << std::endl;

            if (graph->is_canonical_mode())
                config->forward_and_reverse = false;

            config->canonical = graph->is_canonical_mode();

            parse_sequences(files, *config, timer,
                [&graph,&inserted_edges](std::string&& seq) {
                    graph->add_sequence(seq, inserted_edges.get());
                },
                [&graph,&inserted_edges](std::string&& kmer, uint32_t /*count*/) {
                    graph->add_sequence(kmer, inserted_edges.get());
                },
                [&graph,&inserted_edges](const auto &loop) {
                    loop([&graph,&inserted_edges](const char *seq) {
                        graph->add_sequence(seq, inserted_edges.get());
                    });
                }
            );

            if (config->verbose)
                std::cout << "Graph extension finished in "
                          << timer.elapsed() << "sec" << std::endl;
            timer.reset();

            if (node_weights) {
                node_weights->insert_nodes(*inserted_edges);

                assert(node_weights->is_compatible(*graph));

                if (graph->is_canonical_mode())
                    config->forward_and_reverse = true;

                parse_sequences(files, *config, timer,
                    [&graph,&node_weights](std::string&& seq) {
                        graph->map_to_nodes_sequentially(seq.begin(), seq.end(),
                            [&](auto node) { node_weights->add_weight(node, 1); }
                        );
                    },
                    [&graph,&node_weights](std::string&& kmer, uint32_t count) {
                        node_weights->add_weight(graph->kmer_to_node(kmer), count);
                    },
                    [&graph,&node_weights](const auto &loop) {
                        loop([&graph,&node_weights](const char *seq) {
                            std::string seq_str(seq);
                            graph->map_to_nodes_sequentially(seq_str.begin(), seq_str.end(),
                                [&](auto node) { node_weights->add_weight(node, 1); }
                            );
                        });
                    }
                );
            }

            if (config->verbose)
                std::cout << "Node weights updated in "
                          << timer.elapsed() << "sec" << std::endl;

            assert(config->outfbase.size());

            // serialize graph
            timer.reset();

            graph->serialize(config->outfbase);
            graph->serialize_extensions(config->outfbase);
            graph.reset();

            if (config->verbose)
                std::cout << "Serialized in " << timer.elapsed() << "sec" << std::endl;

            timer.reset();

            if (!config->infbase_annotators.size())
                return 0;

            auto annotation = initialize_annotation(config->infbase_annotators.at(0), *config);

            if (!annotation->load(config->infbase_annotators.at(0))) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            } else if (config->verbose) {
                std::cout << "Annotation was loaded in "
                          << timer.elapsed() << "sec" << std::endl;
            }

            timer.reset();

            assert(inserted_edges);

            if (annotation->num_objects() + 1 != inserted_edges->size()
                                                - inserted_edges->num_set_bits()) {
                std::cerr << "ERROR: incompatible graph and annotation." << std::endl;
                exit(1);
            }

            if (config->verbose)
                std::cout << "Insert empty rows to the annotation matrix..." << std::flush;

            AnnotatedDBG::insert_zero_rows(annotation.get(), *inserted_edges);

            if (config->verbose)
                std::cout << "\tdone in " << timer.elapsed() << "sec" << std::endl;

            annotation->serialize(config->outfbase);

            return 0;
        }
        case Config::ANNOTATE: {
            assert(config->infbase_annotators.size() <= 1);

            const auto graph = load_critical_dbg(config->infbase);

            if (graph->is_canonical_mode())
                config->forward_and_reverse = false;

            if (!config->separately) {
                auto anno_graph = initialize_annotated_dbg(graph, *config);

                annotate_data(files,
                              config->refpath,
                              anno_graph.get(),
                              config->forward_and_reverse,
                              config->min_count,
                              config->max_count,
                              config->filename_anno,
                              config->fasta_anno,
                              config->fasta_anno_comment_delim,
                              config->fasta_header_delimiter,
                              config->anno_labels,
                              config->verbose);

                anno_graph->get_annotation().serialize(config->outfbase);

            } else {
                size_t num_threads = config->parallel;
                // annotate multiple columns in parallel, each in a single thread
                config->parallel = 1;

                #pragma omp parallel for num_threads(num_threads) default(shared) schedule(dynamic, 1)
                for (size_t i = 0; i < files.size(); ++i) {
                    auto anno_graph = initialize_annotated_dbg(graph, *config);

                    annotate_data({ files[i] },
                                  config->refpath,
                                  anno_graph.get(),
                                  config->forward_and_reverse,
                                  config->min_count,
                                  config->max_count,
                                  config->filename_anno,
                                  config->fasta_anno,
                                  config->fasta_anno_comment_delim,
                                  config->fasta_header_delimiter,
                                  config->anno_labels,
                                  config->verbose);

                    anno_graph->get_annotation().serialize(
                        config->outfbase.size()
                            ? config->outfbase + "/" + utils::split_string(files[i], "/").back()
                            : files[i]
                    );
                }
            }

            return 0;
        }
        case Config::ANNOTATE_COORDINATES: {
            assert(config->infbase_annotators.size() <= 1);

            auto graph_temp = load_critical_dbg(config->infbase);

            auto annotation_temp
                = std::make_unique<annotate::RowCompressed<>>(graph_temp->num_nodes());

            if (config->infbase_annotators.size()
                    && !annotation_temp->load(config->infbase_annotators.at(0))) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }

            // load graph
            AnnotatedDBG anno_graph(graph_temp,
                                    std::move(annotation_temp),
                                    config->parallel,
                                    config->fast);

            if (!anno_graph.check_compatibility()) {
                std::cerr << "Error: graph and annotation are not compatible."
                          << std::endl;
                exit(1);
            }

            annotate_coordinates(files,
                                 &anno_graph,
                                 config->forward_and_reverse,
                                 config->genome_binsize_anno,
                                 config->verbose);

            anno_graph.get_annotation().serialize(config->outfbase);

            return 0;
        }
        case Config::MERGE_ANNOTATIONS: {
            if (config->anno_type == Config::ColumnCompressed) {
                annotate::ColumnCompressed<> annotation(0, kNumCachedColumns, config->verbose);
                if (!annotation.merge_load(files)) {
                    std::cerr << "ERROR: can't load annotations" << std::endl;
                    exit(1);
                }
                annotation.serialize(config->outfbase);
                return 0;
            }

            std::vector<std::unique_ptr<Annotator>> annotators;
            std::vector<std::string> stream_files;

            for (const auto &filename : files) {
                auto anno_file_type = parse_annotation_type(filename);
                if (anno_file_type == Config::AnnotationType::RowCompressed) {
                    stream_files.push_back(filename);
                } else {
                    auto annotator = initialize_annotation(filename, *config);
                    if (!annotator->load(filename)) {
                        std::cerr << "ERROR: can't load annotation from file "
                                  << filename << std::endl;
                        exit(1);
                    }
                    annotators.push_back(std::move(annotator));
                }
            }

            if (config->anno_type == Config::RowCompressed) {
                annotate::merge<annotate::RowCompressed<>>(std::move(annotators), stream_files, config->outfbase);
            } else if (config->anno_type == Config::RowFlat) {
                annotate::merge<annotate::RowFlatAnnotator>(std::move(annotators), stream_files, config->outfbase);
            } else if (config->anno_type == Config::RBFish) {
                annotate::merge<annotate::RainbowfishAnnotator>(std::move(annotators), stream_files, config->outfbase);
            } else if (config->anno_type == Config::BinRelWT_sdsl) {
                annotate::merge<annotate::BinRelWT_sdslAnnotator>(std::move(annotators), stream_files, config->outfbase);
            } else if (config->anno_type == Config::BinRelWT) {
                annotate::merge<annotate::BinRelWTAnnotator>(std::move(annotators), stream_files, config->outfbase);
            } else if (config->anno_type == Config::BRWT) {
                annotate::merge<annotate::BRWTCompressed<>>(std::move(annotators), stream_files, config->outfbase);
            } else {
                std::cerr << "ERROR: Merging of annotations to '"
                          << config->annotype_to_string(config->anno_type)
                          << "' is not implemented." << std::endl;
                exit(1);
            }

            return 0;
        }
        case Config::QUERY: {
            assert(config->infbase_annotators.size() == 1);

            auto graph = load_critical_dbg(config->infbase);
            auto anno_graph = initialize_annotated_dbg(graph, *config);

            ThreadPool thread_pool(std::max(1u, config->parallel) - 1);

            Timer timer;

            std::unique_ptr<IDBGAligner> aligner;
            // TODO: make aligner work with batch querying
            if (config->align_sequences && !config->fast)
                aligner.reset(build_aligner(*graph, *config).release());

            // iterate over input files
            for (const auto &file : files) {
                if (config->verbose) {
                    std::cout << "\nParsing sequences from " + file + '\n' << std::flush;
                }

                Timer curr_timer;

                size_t seq_count = 0;

                const auto *graph_to_query = anno_graph.get();

                // Graph constructed from a batch of queried sequences
                // Used only in fast mode
                std::unique_ptr<AnnotatedDBG> query_graph;
                if (config->fast) {
                    query_graph = construct_query_graph(*anno_graph,
                        [&](auto call_sequence) {
                            read_fasta_file_critical(file,
                                [&](kseq_t *seq) { call_sequence(seq->seq.s); },
                                config->forward_and_reverse
                            );
                        },
                        config->count_labels ? 0 : config->discovery_fraction,
                        config->parallel
                    );

                    graph_to_query = query_graph.get();

                    if (config->verbose) {
                        std::cout << "Query graph constructed for "
                                        + file + " in "
                                        + std::to_string(curr_timer.elapsed())
                                        + " sec\n" << std::flush;
                    }
                }

                read_fasta_file_critical(file,
                    [&](kseq_t *read_stream) {
                        thread_pool.enqueue(execute_query,
                            fmt::format_int(seq_count++).str() + "\t"
                                + read_stream->name.s,
                            std::string(read_stream->seq.s),
                            config->count_labels,
                            config->suppress_unlabeled,
                            config->num_top_labels,
                            config->discovery_fraction,
                            config->anno_labels_delimiter,
                            std::ref(*graph_to_query),
                            std::ref(std::cout),
                            aligner.get()
                        );
                    },
                    config->forward_and_reverse
                );

                // wait while all threads finish processing the current file
                thread_pool.join();

                if (config->verbose) {
                    std::cout << "File " + file + " was processed in "
                                    + std::to_string(curr_timer.elapsed())
                                    + " sec, total time: "
                                    + std::to_string(timer.elapsed())
                                    + " sec\n" << std::flush;
                }
            }

            return 0;
        }
        case Config::SERVER_QUERY: {
            assert(config->infbase_annotators.size() == 1);

            Timer timer;

            std::cout << "Loading graph..." << std::endl;

            auto graph = load_critical_dbg(config->infbase);
            auto anno_graph = initialize_annotated_dbg(graph, *config);

            std::cout << "Graph loaded in "
                      << timer.elapsed() << "sec, current mem usage: "
                      << (get_curr_RSS() >> 20) << " MiB" << std::endl;

            std::unique_ptr<IDBGAligner> aligner;
            // TODO: make aligner work with batch querying
            if (config->align_sequences && !config->fast)
                aligner.reset(build_aligner(*graph, *config).release());

            const size_t num_threads = std::max(1u, config->parallel);

            std::cout << "Initializing tcp service with "
                      << num_threads << " threads, listening port "
                      << config->port << std::endl;

            try {
                asio::io_context io_context;

                asio::signal_set signals(io_context, SIGINT, SIGTERM);
                signals.async_wait([&](auto, auto) { io_context.stop(); });

                Server server(io_context, config->port,
                    [&](const std::string &received_message) {
                        return form_client_reply(
                            received_message,
                            *anno_graph,
                            *config,
                            aligner.get()
                        );
                    }
                );

                std::vector<std::thread> workers;
                for (size_t i = 0; i < std::max(1u, config->parallel); ++i) {
                    workers.emplace_back([&io_context]() { io_context.run(); });
                }
                for (auto &thread : workers) {
                    thread.join();
                }
            } catch (const std::exception &e) {
                std::cerr << "Exception: " << e.what() << std::endl;
            } catch (...) {
                std::cerr << "Error: Unknown exception" << std::endl;
            }

            return 0;
        }
        case Config::COMPARE: {
            assert(files.size());

            std::cout << "Loading graph                " << files.at(0) << std::endl;
            auto graph = load_critical_dbg(files.at(0));

            for (size_t f = 1; f < files.size(); ++f) {
                std::cout << "Loading graph for comparison " << files[f] << std::endl;
                auto second = load_critical_dbg(files[f]);
                if (*graph == *second) {
                    std::cout << "Graphs are identical" << std::endl;
                } else {
                    std::cout << "Graphs are not identical" << std::endl;
                }
            }

            return 0;
        }
        case Config::CONCATENATE: {
            assert(config->outfbase.size());

            auto chunk_files = files;

            Timer timer;

            if (!files.size()) {
                assert(config->infbase.size());

                const auto sorted_suffixes = config->graph_type == Config::GraphType::SUCCINCT
                        ? KmerExtractorBOSS().generate_suffixes(config->suffix_len)
                        : KmerExtractor2Bit().generate_suffixes(config->suffix_len);

                for (const std::string &suffix : sorted_suffixes) {
                    assert(suffix.size() == config->suffix_len);
                    chunk_files.push_back(config->infbase + "." + suffix);
                }
            }

            if (!chunk_files.size()) {
                std::cerr << "Error: no input files provided, nothing to concatenate" << std::endl;
                exit(1);
            }

            for (auto &filename : chunk_files) {
                filename = utils::remove_suffix(filename,
                                                BOSS::Chunk::kFileExtension,
                                                DBGBitmap::kChunkFileExtension);
            }

            // collect results on an external merge or construction
            std::unique_ptr<DeBruijnGraph> graph;
            switch (config->graph_type) {
                case Config::GraphType::SUCCINCT: {
                    auto p = BOSS::Chunk::build_boss_from_chunks(chunk_files, config->verbose);
                    auto dbg_succ = std::make_unique<DBGSuccinct>(p.first, p.second);

                    if (config->verbose) {
                        std::cout << "Chunks concatenated in "
                                  << timer.elapsed() << "sec" << std::endl;
                    }

                    if (config->clear_dummy) {
                        if (config->verbose) {
                            std::cout << "Traverse source dummy edges,"
                                      << " remove redundant ones, and mark"
                                      << " those that cannot be removed."
                                      << std::endl;
                        }
                        dbg_succ->mask_dummy_kmers(config->parallel, true);
                    }
                    graph = std::move(dbg_succ);
                    break;
                }
                case Config::GraphType::BITMAP: {
                    graph.reset(DBGBitmapConstructor::build_graph_from_chunks(
                        chunk_files, config->canonical, config->verbose
                    ));
                    break;
                }
                default:
                    std::cout << "ERROR: Cannot concatenate chunks for "
                              << "this graph representation" << std::endl;
                    exit(1);
            }
            assert(graph);

            if (config->verbose) {
                std::cout << "Graph was assembled in "
                          << timer.elapsed() << "sec" << std::endl;
                print_stats(*graph);
                if (config->graph_type == Config::GraphType::SUCCINCT) {
                    print_boss_stats(
                        dynamic_cast<DBGSuccinct*>(graph.get())->get_boss()
                    );
                }
            }

            // graph output
            graph->serialize(config->outfbase);

            return 0;
        }
        case Config::MERGE: {
            BOSS *graph = NULL;

            Timer timer;

            std::vector<std::shared_ptr<DBGSuccinct>> dbg_graphs;
            std::vector<const BOSS*> graphs;

            config->canonical = true;

            for (const auto &file : files) {
                std::cout << "Opening file " << file << std::endl;

                dbg_graphs.emplace_back(load_critical_graph_from_file<DBGSuccinct>(file));

                graphs.push_back(&dbg_graphs.back()->get_boss());

                if (config->verbose)
                    print_boss_stats(*graphs.back());

                config->canonical &= dbg_graphs.back()->is_canonical_mode();
            }

            std::cout << "Graphs are loaded in " << timer.elapsed()
                                                 << "sec" << std::endl;

            if (config->dynamic) {
                std::cout << "Start merging traversal" << std::endl;
                timer.reset();

                graph = dbg_graphs.at(0)->release_boss();

                if (graph->get_state() != Config::DYN) {
                    if (config->verbose)
                        std::cout << "Switching state of succinct graph to dynamic..." << std::flush;

                    graph->switch_state(Config::DYN);

                    if (config->verbose)
                        std::cout << "\tdone in " << timer.elapsed() << "sec" << std::endl;
                }

                for (size_t i = 1; i < graphs.size(); ++i) {
                    graph->merge(dbg_graphs.at(i)->get_boss());

                    std::cout << "traversal " << files[i] << " done\t"
                              << timer.elapsed() << "sec" << std::endl;

                    dbg_graphs.at(i).reset();
                }
            } else if (config->parallel > 1 || config->parts_total > 1) {
                std::cout << "Start merging blocks" << std::endl;
                timer.reset();

                auto *chunk = merge::merge_blocks_to_chunk(
                    graphs,
                    config->part_idx,
                    config->parts_total,
                    config->parallel,
                    config->num_bins_per_thread,
                    config->verbose
                );
                if (!chunk) {
                    std::cerr << "ERROR when building chunk "
                              << config->part_idx << std::endl;
                    exit(1);
                }
                std::cout << "Blocks merged\t" << timer.elapsed()
                          << "sec" << std::endl;

                if (config->parts_total > 1) {
                    chunk->serialize(config->outfbase
                                      + "." + std::to_string(config->part_idx)
                                      + "_" + std::to_string(config->parts_total));
                } else {
                    graph = new BOSS(graphs[0]->get_k());
                    chunk->initialize_boss(graph);
                }
                delete chunk;
            } else {
                std::cout << "Start merging graphs" << std::endl;
                timer.reset();

                graph = merge::merge(graphs, config->verbose);
            }
            dbg_graphs.clear();

            assert(graph);

            std::cout << "Graphs merged in " << timer.elapsed() << "sec" << std::endl;

            // graph output
            DBGSuccinct(graph, config->canonical).serialize(config->outfbase);

            return 0;
        }
        case Config::CLEAN: {
            assert(files.size() == 1);
            assert(config->outfbase.size());

            config->min_count = std::max(1u, config->min_count);

            if (!config->to_fasta) {
                std::cerr << "Error: Clean graph can be serialized only in"
                          << " form of contigs/unitigs, add flag --to-fasta" << std::endl;
                exit(1);
            }

            Timer timer;
            if (config->verbose)
                std::cout << "Graph loading...\t" << std::flush;

            auto graph = load_critical_dbg(files.at(0));

            if (config->min_count > 1
                    || config->max_count < std::numeric_limits<unsigned int>::max()
                    || config->min_unitig_median_kmer_abundance != 1
                    || config->count_slice_quantiles[0] != 0
                    || config->count_slice_quantiles[1] != 1) {
                // load k-mer counts
                auto node_weights = graph->load_extension<NodeWeights>(files.at(0));

                if (!(node_weights)) {
                    std::cerr << "ERROR: Cannot load k-mer counts from file "
                              << files.at(0) << std::endl;
                    exit(1);
                }

                if (auto *dbg_succ = dynamic_cast<DBGSuccinct*>(graph.get()))
                    dbg_succ->reset_mask();

                if (!node_weights->is_compatible(*graph)) {
                    std::cerr << "Error: k-mer counts are not compatible with graph "
                              << files.at(0) << std::endl;
                    exit(1);
                }

                if (config->min_count > 1
                        || config->max_count < std::numeric_limits<unsigned int>::max()) {
                    const auto &weights = *graph->get_extension<NodeWeights>();

                    graph = std::make_shared<MaskedDeBruijnGraph>(graph,
                        [&](auto i) { return weights[i] >= config->min_count
                                            && weights[i] <= config->max_count; });
                    graph->add_extension(node_weights);

                    assert(node_weights->is_compatible(*graph));
                }

                if (config->min_unitig_median_kmer_abundance == 0) {
                    // skip zero k-mer counts for dummy k-mers in DBGSuccinct
                    const auto _graph = dynamic_cast<DBGSuccinct*>(graph.get())
                            ? std::make_shared<MaskedDeBruijnGraph>(graph, [&](auto i) { return (*node_weights)[i] > 0; })
                            : graph;

                    config->min_unitig_median_kmer_abundance
                        = estimate_min_kmer_abundance(*_graph, *node_weights,
                                                      config->fallback_abundance_cutoff);
                }
            }

            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            if (config->verbose) {
                if (dynamic_cast<const MaskedDeBruijnGraph*>(graph.get())) {
                    std::cout << "Extracting sequences from subgraph..." << std::endl;
                } else {
                    std::cout << "Extracting sequences from graph..." << std::endl;
                }
            }

            timer.reset();

            auto call_clean_contigs = [&](auto callback) {
                if (config->min_unitig_median_kmer_abundance != 1) {
                    auto node_weights = graph->get_extension<NodeWeights>();
                    assert(node_weights);
                    if (!node_weights->is_compatible(*graph)) {
                        std::cerr << "Error: k-mer counts are not compatible with subgraph" << std::endl;
                        exit(1);
                    }

                    std::cout << "Threshold for median k-mer abundance in unitigs: "
                              << config->min_unitig_median_kmer_abundance << std::endl;

                    graph->call_unitigs([&](const std::string &unitig, const auto &path) {
                        if (!is_unreliable_unitig(path,
                                                  *node_weights,
                                                  config->min_unitig_median_kmer_abundance))
                            callback(unitig, path);
                    }, config->min_tip_size);

                } else if (config->unitigs || config->min_tip_size > 1) {
                    graph->call_unitigs(callback, config->min_tip_size);

                } else {
                    graph->call_sequences(callback);
                }
            };

            assert(config->count_slice_quantiles.size() >= 2);

            if (config->count_slice_quantiles[0] == 0
                    && config->count_slice_quantiles[1] == 1) {
                FastaWriter writer(utils::remove_suffix(config->outfbase, ".gz", ".fasta") + ".fasta.gz",
                                   config->header, true);

                call_clean_contigs([&](const std::string &contig, const auto &) {
                    writer.write(contig);
                });

            } else {
                auto node_weights = graph->get_extension<NodeWeights>();
                if (!node_weights) {
                    std::cerr << "Error: need k-mer counts for binning k-mers by abundance"
                              << std::endl;
                    exit(1);
                }
                assert(node_weights->is_compatible(*graph));

                auto &weights = node_weights->get_data();

                assert(graph->num_nodes() + 1 == weights.size());

                // compute clean count histogram
                std::unordered_map<uint64_t, uint64_t> count_hist;

                if (config->min_unitig_median_kmer_abundance != 1 || config->min_tip_size > 1) {
                    // cleaning required
                    sdsl::bit_vector removed_nodes(weights.size(), 1);

                    call_clean_contigs([&](const std::string&, const auto &path) {
                        for (auto i : path) {
                            assert(weights[i]);
                            count_hist[weights[i]]++;
                            removed_nodes[i] = 0;
                        }
                    });

                    call_ones(removed_nodes, [&weights](auto i) { weights[i] = 0; });

                } else if (auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(graph)) {
                    // use entire graph without dummy BOSS edges
                    graph->call_nodes([&](auto i) {
                        if (uint64_t count = weights[i])
                            count_hist[count]++;
                    });
                } else {
                    // use entire graph
                    graph->call_nodes([&](auto i) {
                        assert(weights[i]);
                        count_hist[weights[i]]++;
                    });
                }
                // must not have any zero weights
                assert(!count_hist.count(0));

                std::vector<std::pair<uint64_t, uint64_t>> count_hist_v(count_hist.begin(),
                                                                        count_hist.end());
                count_hist.clear();

                ips4o::parallel::sort(count_hist_v.begin(), count_hist_v.end(),
                    [](const auto &first, const auto &second) {
                        return first.first < second.first;
                    },
                    config->parallel
                );

                #pragma omp parallel for num_threads(config->parallel)
                for (size_t i = 1; i < config->count_slice_quantiles.size(); ++i) {
                    // extract sequences for k-mer counts bin |i|
                    assert(config->count_slice_quantiles[i - 1] < config->count_slice_quantiles[i]);

                    FastaWriter writer(utils::remove_suffix(config->outfbase, ".gz", ".fasta")
                                        + "." + std::to_string(config->count_slice_quantiles[i - 1])
                                        + "." + std::to_string(config->count_slice_quantiles[i]) + ".fasta.gz",
                                       config->header, true);

                    if (!count_hist_v.size())
                        continue;

                    uint64_t min_count = config->count_slice_quantiles[i - 1] > 0
                        ? utils::get_quantile(count_hist_v, config->count_slice_quantiles[i - 1])
                        : 1;
                    uint64_t max_count = config->count_slice_quantiles[i] < 1
                        ? utils::get_quantile(count_hist_v, config->count_slice_quantiles[i])
                        : std::numeric_limits<uint64_t>::max();

                    std::cout << "Used k-mer count thresholds:\n"
                              << "min (including): " << min_count << "\n"
                              << "max (excluding): " << max_count << std::endl;

                    assert(node_weights->is_compatible(*graph));

                    MaskedDeBruijnGraph graph_slice(graph,
                        [&](auto i) { return weights[i] >= min_count && weights[i] < max_count; });

                    graph_slice.call_unitigs([&](const auto &contig, auto&&) { writer.write(contig); });
                }
            }

            if (config->verbose)
                std::cout << "Graph cleaning finished in "
                          << timer.elapsed() << "sec" << std::endl;

            return 0;
        }
        case Config::STATS: {
            for (const auto &file : files) {
                std::shared_ptr<DeBruijnGraph> graph;

                graph = load_critical_dbg(file);
                graph->load_extension<NodeWeights>(file);

                std::cout << "Statistics for graph " << file << std::endl;

                print_stats(*graph);

                if (dynamic_cast<DBGSuccinct*>(graph.get())) {
                    const auto &boss_graph = dynamic_cast<DBGSuccinct&>(*graph).get_boss();

                    print_boss_stats(boss_graph,
                                     config->count_dummy,
                                     config->parallel,
                                     config->verbose);

                    if (config->print_graph_internal_repr)
                        boss_graph.print_internal_representation(std::cout);
                }

                if (config->print_graph)
                    std::cout << *graph;
            }

            for (const auto &file : config->infbase_annotators) {
                auto annotation = initialize_annotation(file, *config);

                if (config->print_column_names) {
                    annotate::LabelEncoder<std::string> label_encoder;

                    std::cout << "INFO: Scanning annotation " << file << std::endl;

                    try {
                        std::ifstream instream(file, std::ios::binary);

                        // TODO: make this more reliable
                        if (dynamic_cast<const annotate::ColumnCompressed<> *>(annotation.get())) {
                            // Column compressed dumps the number of rows first
                            // skipping it...
                            load_number(instream);
                        }

                        if (!label_encoder.load(instream))
                            throw std::ios_base::failure("");

                    } catch (...) {
                        std::cerr << "Error: Can't read label encoder from file "
                                  << file << std::endl;
                        exit(1);
                    }

                    std::cout << "INFO: Number of columns: " << label_encoder.size() << std::endl;
                    for (size_t c = 0; c < label_encoder.size(); ++c) {
                        std::cout << label_encoder.decode(c) << std::endl;
                    }

                    continue;
                }

                if (!annotation->load(file)) {
                    std::cerr << "ERROR: can't load annotation from file "
                              << file << std::endl;
                    exit(1);
                }

                std::cout << "Statistics for annotation " << file << std::endl;
                print_stats(*annotation);
            }

            return 0;
        }
        case Config::TRANSFORM_ANNOTATION: {
            assert(files.size() == 1);

            Timer timer;

            /********************************************************/
            /***************** dump labels to text ******************/
            /********************************************************/

            if (config->dump_raw_anno || config->dump_text_anno) {
                const Config::AnnotationType input_anno_type
                    = parse_annotation_type(files.at(0));

                auto annotation = initialize_annotation(files.at(0), *config);

                if (config->verbose)
                    std::cout << "Loading annotation..." << std::endl;

                if (config->anno_type == Config::ColumnCompressed) {
                    if (!annotation->merge_load(files)) {
                        std::cerr << "ERROR: can't load annotations" << std::endl;
                        exit(1);
                    }
                } else {
                    // Load annotation from disk
                    if (!annotation->load(files.at(0))) {
                        std::cerr << "ERROR: can't load annotation from file "
                                  << files.at(0) << std::endl;
                        exit(1);
                    }
                }

                if (config->verbose) {
                    std::cout << "Annotation loaded in "
                              << timer.elapsed() << "sec" << std::endl;

                    std::cout << "Dumping annotators...\t" << std::flush;
                }

                if (input_anno_type == Config::ColumnCompressed) {
                    assert(dynamic_cast<annotate::ColumnCompressed<>*>(annotation.get()));
                    if (config->dump_raw_anno) {
                        dynamic_cast<annotate::ColumnCompressed<>*>(
                            annotation.get()
                        )->dump_columns(config->outfbase, true, get_num_threads());
                    }

                    if (config->dump_text_anno) {
                        dynamic_cast<annotate::ColumnCompressed<>*>(
                            annotation.get()
                        )->dump_columns(config->outfbase, false, get_num_threads());
                    }
                } else if (input_anno_type == Config::BRWT) {
                    assert(dynamic_cast<annotate::BRWTCompressed<>*>(annotation.get()));
                    if (config->dump_raw_anno) {
                        dynamic_cast<annotate::BRWTCompressed<>*>(
                            annotation.get()
                        )->dump_columns(config->outfbase, true, get_num_threads());
                    }

                    if (config->dump_text_anno) {
                        dynamic_cast<annotate::BRWTCompressed<>*>(
                            annotation.get()
                        )->dump_columns(config->outfbase, false, get_num_threads());
                    }
                } else {
                    throw std::runtime_error("Dumping columns for this type not implemented");
                }

                if (config->verbose)
                    std::cout << timer.elapsed() << "sec" << std::endl;

                return 0;
            }

            /********************************************************/
            /***************** rename column labels *****************/
            /********************************************************/

            if (config->rename_instructions_file.size()) {
                std::unordered_map<std::string, std::string> dict;
                std::ifstream instream(config->rename_instructions_file);
                if (!instream.is_open()) {
                    std::cerr << "ERROR: Can't open file "
                              << config->rename_instructions_file << std::endl;
                    exit(1);
                }
                std::string old_name;
                std::string new_name;
                while (instream.good() && !(instream >> old_name).eof()) {
                    instream >> new_name;
                    if (instream.fail() || instream.eof()) {
                        std::cerr << "ERROR: wrong format of the rules for"
                                  << " renaming annotation columns passed in file "
                                  << config->rename_instructions_file << std::endl;
                        exit(1);
                    }
                    dict[old_name] = new_name;
                }

                auto annotation = initialize_annotation(files.at(0), *config);

                if (config->verbose)
                    std::cout << "Loading annotation..." << std::endl;

                // TODO: rename columns without loading the full annotation
                if (config->anno_type == Config::ColumnCompressed) {
                    if (!annotation->merge_load(files)) {
                        std::cerr << "ERROR: can't load annotations" << std::endl;
                        exit(1);
                    } else {
                        std::cout << annotation->num_objects() << " " << annotation->num_labels() << "\n";
                    }
                } else {
                    // Load annotation from disk
                    if (!annotation->load(files.at(0))) {
                        std::cerr << "ERROR: can't load annotation from file "
                              << files.at(0) << std::endl;
                        exit(1);
                    }
                }

                if (config->verbose) {
                    std::cout << "Annotation loaded in "
                              << timer.elapsed() << "sec" << std::endl;
                }

                if (config->verbose)
                    std::cout << "Renaming...\t" << std::flush;

                //TODO: could be made to work with streaming
                annotation->rename_labels(dict);

                annotation->serialize(config->outfbase);
                if (config->verbose)
                    std::cout << timer.elapsed() << "sec" << std::endl;

                return 0;
            }

            /********************************************************/
            /****************** convert annotation ******************/
            /********************************************************/

            const Config::AnnotationType input_anno_type
                = parse_annotation_type(files.at(0));

            if (config->anno_type == input_anno_type) {
                std::cerr << "Skipping conversion: same input and target type: "
                          << Config::annotype_to_string(config->anno_type)
                          << std::endl;
                exit(1);
            }

            if (input_anno_type == Config::ColumnCompressed && files.size() > 1) {
                std::cerr << "ERROR: conversion of multiple annotators only supported for ColumnCompressed" << std::endl;
                exit(1);
            }

            if (config->verbose) {
                std::cout << "Converting to " << Config::annotype_to_string(config->anno_type)
                          << " annotator..." << std::endl;
            }

            if (input_anno_type == Config::RowCompressed) {

                std::unique_ptr<const Annotator> target_annotator;

                switch (config->anno_type) {
                    case Config::RowFlat: {
                        auto annotator = annotate::convert<annotate::RowFlatAnnotator>(files.at(0));
                        target_annotator = std::move(annotator);
                        break;
                    }
                    case Config::RBFish: {
                        auto annotator = annotate::convert<annotate::RainbowfishAnnotator>(files.at(0));
                        target_annotator = std::move(annotator);
                        break;
                    }
                    case Config::BinRelWT_sdsl: {
                        auto annotator = annotate::convert<annotate::BinRelWT_sdslAnnotator>(files.at(0));
                        target_annotator = std::move(annotator);
                        break;
                    }
                    case Config::BinRelWT: {
                        auto annotator = annotate::convert<annotate::BinRelWTAnnotator>(files.at(0));
                        target_annotator = std::move(annotator);
                        break;
                    }
                    default:
                        std::cerr << "Error: Streaming conversion from RowCompressed annotation"
                                  << " is not implemented for the requested target type: "
                                  << Config::annotype_to_string(config->anno_type)
                                  << std::endl;
                        exit(1);
                }

                if (config->verbose) {
                    std::cout << "Annotation converted in "
                              << timer.elapsed() << "sec" << std::endl;
                }

                if (config->verbose) {
                    std::cout << "Serializing to " << config->outfbase
                              << "...\t" << std::flush;
                }

                target_annotator->serialize(config->outfbase);

                if (config->verbose) {
                    std::cout << timer.elapsed() << "sec" << std::endl;
                }

            } else if (input_anno_type == Config::ColumnCompressed) {
                auto annotation = initialize_annotation(files.at(0), *config);

                if (config->verbose)
                    std::cout << "Loading annotation..." << std::endl;

                // Load annotation from disk
                if (!annotation->merge_load(files)) {
                    std::cerr << "ERROR: can't load annotations" << std::endl;
                    exit(1);
                }

                if (config->verbose) {
                    std::cout << "Annotation loaded in "
                              << timer.elapsed() << "sec" << std::endl;
                }

                std::unique_ptr<annotate::ColumnCompressed<>> annotator {
                    dynamic_cast<annotate::ColumnCompressed<> *>(annotation.release())
                };
                assert(annotator);

                switch (config->anno_type) {
                    case Config::ColumnCompressed: {
                        assert(false);
                        break;
                    }
                    case Config::RowCompressed: {
                        if (config->fast) {
                            annotate::RowCompressed<> row_annotator(0);
                            annotator->convert_to_row_annotator(&row_annotator,
                                                                config->parallel);
                            annotator.reset();

                            if (config->verbose) {
                                std::cout << "Annotation converted in "
                                          << timer.elapsed() << "sec" << std::endl;
                            }

                            if (config->verbose) {
                                std::cout << "Serializing to " << config->outfbase
                                          << "...\t" << std::flush;
                            }

                            row_annotator.serialize(config->outfbase);

                            if (config->verbose) {
                                std::cout << timer.elapsed() << "sec" << std::endl;
                            }

                        } else {
                            annotator->convert_to_row_annotator(config->outfbase);
                            if (config->verbose) {
                                std::cout << "Annotation converted and serialized in "
                                          << timer.elapsed() << "sec" << std::endl;
                            }
                        }
                        break;
                    }
                    case Config::BRWT: {
                        auto brwt_annotator = config->greedy_brwt
                            ? annotate::convert_to_greedy_BRWT<annotate::BRWTCompressed<>>(
                                std::move(*annotator),
                                config->parallel_nodes,
                                config->parallel)
                            : annotate::convert_to_simple_BRWT<annotate::BRWTCompressed<>>(
                                std::move(*annotator),
                                config->arity_brwt,
                                config->parallel_nodes,
                                config->parallel);

                        annotator.reset();

                        if (config->verbose) {
                            std::cout << "Annotation converted in "
                                      << timer.elapsed() << "sec" << std::endl;
                        }

                        if (config->verbose) {
                            std::cout << "Serializing to " << config->outfbase
                                      << "...\t" << std::flush;
                        }

                        brwt_annotator->serialize(config->outfbase);

                        if (config->verbose) {
                            std::cout << timer.elapsed() << "sec" << std::endl;
                        }
                        break;
                    }
                    case Config::BinRelWT_sdsl: {
                        convert<annotate::BinRelWT_sdslAnnotator>(std::move(annotator), *config, timer);
                        break;
                    }
                    case Config::BinRelWT: {
                        convert<annotate::BinRelWTAnnotator>(std::move(annotator), *config, timer);
                        break;
                    }
                    case Config::RowFlat: {
                        convert<annotate::RowFlatAnnotator>(std::move(annotator), *config, timer);
                        break;
                    }
                    case Config::RBFish: {
                        convert<annotate::RainbowfishAnnotator>(std::move(annotator), *config, timer);
                        break;
                    }
                }

            } else {
                std::cerr << "Error: Conversion to other representations"
                          << " is not implemented for "
                          << Config::annotype_to_string(input_anno_type)
                          << " annotator." << std::endl;
                exit(1);
            }

            return 0;
        }
        case Config::TRANSFORM: {
            assert(files.size() == 1);
            assert(config->outfbase.size());

            if (config->initialize_bloom
                    && parse_graph_extension(files.at(0)) == Config::GraphType::SUCCINCT)
                std::filesystem::remove(
                    utils::remove_suffix(config->outfbase, ".bloom") + ".bloom"
                );

            Timer timer;
            if (config->verbose)
                std::cout << "Graph loading...\t" << std::flush;

            auto graph = load_critical_dbg(files.at(0));

            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(graph);

            if (!dbg_succ.get())
                throw std::runtime_error("Only implemented for DBGSuccinct");

            if (config->initialize_bloom) {
                assert(config->bloom_fpp > 0.0 && config->bloom_fpp <= 1.0);
                assert(config->bloom_bpk >= 0.0);
                assert(config->bloom_fpp < 1.0 || config->bloom_bpk > 0.0);

                if (config->verbose) {
                    std::cout << "Construct Bloom filter for nodes..." << std::endl;
                }

                timer.reset();

                if (config->bloom_fpp < 1.0) {
                    dbg_succ->initialize_bloom_filter_from_fpr(
                        config->bloom_fpp,
                        config->bloom_max_num_hash_functions
                    );
                } else {
                    dbg_succ->initialize_bloom_filter(
                        config->bloom_bpk,
                        config->bloom_max_num_hash_functions
                    );
                }

                if (config->verbose)
                    std::cout << timer.elapsed() << "sec" << std::endl;

                assert(dbg_succ->get_bloom_filter());

                auto prefix = utils::remove_suffix(config->outfbase, dbg_succ->bloom_filter_file_extension());
                std::ofstream bloom_outstream(
                    prefix + dbg_succ->bloom_filter_file_extension(), std::ios::binary
                );

                if (!bloom_outstream.good())
                    throw std::ios_base::failure("Can't write to file " + prefix + dbg_succ->bloom_filter_file_extension());

                dbg_succ->get_bloom_filter()->serialize(bloom_outstream);

                return 0;
            }

            if (config->clear_dummy) {
                if (config->verbose) {
                    std::cout << "Traverse source dummy edges and remove redundant ones..." << std::endl;
                }
                timer.reset();

                // remove redundant dummy edges and mark all other dummy edges
                dbg_succ->mask_dummy_kmers(config->parallel, true);

                if (config->verbose)
                    std::cout << "Done in " << timer.elapsed() << "sec" << std::endl;

                timer.reset();
            }

            if (config->to_adj_list) {
                if (config->verbose)
                    std::cout << "Converting graph to adjacency list...\t" << std::flush;

                auto *boss = &dbg_succ->get_boss();
                timer.reset();

                std::ofstream outstream(config->outfbase + ".adjlist");
                boss->print_adj_list(outstream);

                if (config->verbose)
                    std::cout << timer.elapsed() << "sec" << std::endl;

                return 0;
            }

            if (config->verbose) {
                std::cout << "Converting graph to state "
                          << Config::state_to_string(config->state)
                          << "...\t" << std::flush;
                timer.reset();
            }

            dbg_succ->switch_state(config->state);

            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            if (config->verbose) {
                std::cout << "Serializing transformed graph...\t" << std::flush;
                timer.reset();
            }
            dbg_succ->serialize(config->outfbase);
            if (config->verbose) {
                std::cout << timer.elapsed() << "sec" << std::endl;
            }

            return 0;
        }
        case Config::ASSEMBLE: {
            assert(files.size() == 1);
            assert(config->outfbase.size());

            Timer timer;
            if (config->verbose)
                std::cout << "Graph loading...\t" << std::flush;

            auto graph = load_critical_dbg(files.at(0));

            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            std::unique_ptr<AnnotatedDBG> anno_graph;
            if (config->infbase_annotators.size()) {
                anno_graph = initialize_annotated_dbg(graph, *config);

                if (config->verbose) {
                    std::cout << "Masking graph...\t" << std::flush;
                }

                graph = mask_graph(*anno_graph, config.get());

                if (config->verbose) {
                    std::cout << timer.elapsed() << "sec" << std::endl;
                }
            }

            if (config->verbose)
                std::cout << "Extracting sequences from graph...\t" << std::flush;

            timer.reset();

            if (config->to_gfa) {
                if (!config->unitigs) {
                    std::cerr << "'--unitigs' must be set for GFA output" << std::endl;
                    exit(1);
                }

                if (config->verbose)
                    std::cout << "Writing graph to GFA...\t" << std::flush;

                std::ofstream gfa_file(utils::remove_suffix(config->outfbase, ".gfa") + ".gfa");

                gfa_file << "H\tVN:Z:1.0" << std::endl;
                graph->call_unitigs(
                    [&](const auto &unitig, const auto &path) {
                        gfa_file << "S\t" << path.back() << "\t" << unitig << std::endl;
                        graph->adjacent_incoming_nodes(path.front(), [&](uint64_t node) {
                            gfa_file << "L\t" << node << "\t+\t" << path.back() << "\t+\t0M" << std::endl;
                        });
                    },
                    config->min_tip_size
                );
            }

            FastaWriter writer(utils::remove_suffix(config->outfbase, ".gz", ".fasta") + ".fasta.gz",
                               config->header, true);

            if (config->unitigs || config->min_tip_size > 1) {
                graph->call_unitigs([&](const auto &unitig, auto&&) { writer.write(unitig); },
                                    config->min_tip_size,
                                    config->kmers_in_single_form);
            } else {
                graph->call_sequences([&](const auto &contig, auto&&) { writer.write(contig); },
                                      config->kmers_in_single_form);
            }

            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            return 0;
        }
        case Config::RELAX_BRWT: {
            assert(files.size() == 1);
            assert(config->outfbase.size());

            Timer timer;

            auto annotator = std::make_unique<annotate::BRWTCompressed<>>();

            if (config->verbose)
                std::cout << "Loading annotator...\t" << std::flush;

            if (!annotator->load(files.at(0))) {
                std::cerr << "ERROR: can't load annotations from file "
                          << files.at(0) << std::endl;
                exit(1);
            }
            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            if (config->verbose)
                std::cout << "Relaxing BRWT tree...\t" << std::flush;

            annotate::relax_BRWT<annotate::BRWTCompressed<>>(annotator.get(),
                                                             config->relax_arity_brwt,
                                                             config->parallel);

            annotator->serialize(config->outfbase);
            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            return 0;
        }
        case Config::ALIGN: {
            assert(config->infbase.size());

            // initialize aligner
            auto graph = load_critical_dbg(config->infbase);
            auto dbg = std::dynamic_pointer_cast<DBGSuccinct>(graph);

            // This speeds up mapping, and allows for node suffix matching
            if (dbg)
                dbg->reset_mask();

            Timer timer;
            ThreadPool thread_pool(std::max(1u, config->parallel) - 1);
            std::mutex print_mutex;

            if (config->map_sequences) {
                if (!config->alignment_length) {
                    config->alignment_length = graph->get_k();
                } else if (config->alignment_length > graph->get_k()) {
                    std::cerr << "Warning: Mapping to k-mers"
                              << " longer than k is not supported." << std::endl;
                    config->alignment_length = graph->get_k();
                }

                if (!dbg && config->alignment_length != graph->get_k()) {
                    std::cerr << "Error: matching k-mers shorter than k only supported for DBGSuccinct"
                              << std::endl;
                    exit(1);
                }

                if (utils::get_verbose()) {
                    std::cout << "Map sequences against the de Bruijn graph with "
                              << "k = " << graph->get_k() << "\n"
                              << "Length of mapped k-mers: "
                              << config->alignment_length << std::endl;
                }

                for (const auto &file : files) {
                    if (utils::get_verbose())
                        std::cout << "Map sequences from file " << file << std::endl;

                    map_sequences_in_file(file,
                                          *graph,
                                          dbg,
                                          *config,
                                          timer,
                                          &thread_pool,
                                          &print_mutex);
                }

                thread_pool.join();

                return 0;
            }

            auto aligner = build_aligner(*graph, *config);

            for (const auto &file : files) {
                std::cout << "Align sequences from file " << file << std::endl;

                Timer data_reading_timer;

                std::ostream *outstream = config->outfbase.size()
                    ? new std::ofstream(config->outfbase)
                    : &std::cout;

                Json::StreamWriterBuilder builder;
                builder["indentation"] = "";

                read_fasta_file_critical(file, [&](kseq_t *read_stream) {
                    thread_pool.enqueue([&](std::string query,
                                            std::string header) {
                            auto paths = aligner->align(query);

                            std::ostringstream ostr;
                            if (!config->output_json) {
                                for (const auto &path : paths) {
                                    const auto& path_query = path.get_orientation()
                                        ? paths.get_query_reverse_complement()
                                        : paths.get_query();

                                    ostr << header << "\t"
                                         << path_query << "\t"
                                         << path
                                         << std::endl;
                                }

                                if (paths.empty())
                                    ostr << header << "\t"
                                         << query << "\t"
                                         << "*\t*\t"
                                         << config->alignment_min_path_score << "\t*\t*\t*"
                                         << std::endl;
                            } else {
                                bool secondary = false;
                                for (const auto &path : paths) {
                                    const auto& path_query = path.get_orientation()
                                        ? paths.get_query_reverse_complement()
                                        : paths.get_query();

                                    ostr << Json::writeString(
                                                builder,
                                                path.to_json(path_query,
                                                             *graph,
                                                             secondary,
                                                             header)
                                            )
                                         << std::endl;

                                    secondary = true;
                                }

                                if (paths.empty()) {
                                    ostr << Json::writeString(
                                                builder,
                                                DBGAligner<>::DBGAlignment().to_json(
                                                    query,
                                                    *graph,
                                                    secondary,
                                                    header)
                                            )
                                         << std::endl;
                                }
                            }

                            auto lock = std::lock_guard<std::mutex>(print_mutex);
                            *outstream << ostr.str();
                        },
                        std::string(read_stream->seq.s),
                        config->fasta_anno_comment_delim != Config::UNINITIALIZED_STR
                            && read_stream->comment.l
                                ? utils::join_strings(
                                    { read_stream->name.s, read_stream->comment.s },
                                    config->fasta_anno_comment_delim,
                                    true)
                                : std::string(read_stream->name.s)
                    );

                    if (config->verbose) {
                        std::cout << "File processed in "
                                  << data_reading_timer.elapsed()
                                  << "sec, current mem usage: "
                                  << (get_curr_RSS() >> 20) << " MiB"
                                  << ", total time: " << timer.elapsed()
                                  << "sec" << std::endl;
                    }
                });

                thread_pool.join();

                if (config->outfbase.size())
                    delete outstream;
            }

            return 0;
        }
        case Config::CALL_VARIANTS: {
            assert(config->infbase_annotators.size() == 1);

            std::unique_ptr<TaxIDMapper> taxid_mapper;
            if (config->taxonomy_map.length()) {
                taxid_mapper.reset(new TaxIDMapper());
                std::ifstream taxid_mapper_in(config->taxonomy_map, std::ios::binary);
                if (!taxid_mapper->load(taxid_mapper_in)) {
                    std::cerr << "ERROR: failed to read accession2taxid map" << std::endl;
                    exit(1);
                }
            }

            auto anno_graph = initialize_annotated_dbg(*config);
            auto masked_graph = mask_graph(*anno_graph, config.get());

            if (config->verbose) {
                std::cout << "Filter out:";
                for (const auto &out : config->label_filter) {
                    std::cout << " " << out;
                }
                std::cout << std::endl;
            }

            std::ostream *outstream = config->outfbase.size()
                ? new std::ofstream(config->outfbase)
                : &std::cout;

            std::unique_ptr<Json::StreamWriter> json_writer;
            if (config->output_json) {
                Json::StreamWriterBuilder builder;
                builder["indentation"] = "";
                json_writer.reset(builder.newStreamWriter());
            } else {
                *outstream << "Index"
                           << "\t" << "Ref"
                           << "\t" << "Var"
                           << "\t" << "Label";

                if (taxid_mapper.get())
                    *outstream << "\t" << "TaxID";

                *outstream << std::endl;
            }

            std::sort(config->label_filter.begin(), config->label_filter.end());

            ThreadPool thread_pool(std::max(1u, config->parallel) - 1);
            std::mutex print_label_mutex;
            std::atomic_uint64_t num_calls = 0;

            auto mask_in_labels = utils::join_strings(config->label_mask_in, ",");

            auto print_variant =
                [&](auto&& alignment, const std::string &query, auto&& vlabels) {
                    // filter out labels
                    std::sort(vlabels.begin(), vlabels.end());

                    auto it = config->label_filter.begin();
                    for (const auto &label : vlabels) {
                        while (it != config->label_filter.end() && *it < label)
                            ++it;

                        if (it == config->label_filter.end())
                            break;

                        if (*it == label)
                            return;
                    }

                    num_calls++;

                    auto label = utils::join_strings(vlabels, ",");

                    // map labels to Taxonomy IDs
                    if (taxid_mapper.get()) {
                        label += "\t" + std::accumulate(
                            vlabels.begin(),
                            vlabels.end(),
                            std::string(),
                            [&](std::string &taxids, const std::string &label) {
                                return std::move(taxids) + ","
                                    + std::to_string(taxid_mapper->gb_to_taxid(label));
                            }
                        );
                    }

                    // print labels
                    std::lock_guard<std::mutex> lock(print_label_mutex);
                    if (config->output_json) {
                        json_writer->write(
                            alignment.to_json(
                                query,
                                masked_graph->get_graph(),
                                false,
                                mask_in_labels + ":" + std::to_string(num_calls),
                                label
                            ),
                            outstream
                        );

                        *outstream << std::endl;
                    } else {
                        std::cout << alignment.front() << "\t"
                                  << query << "\t"
                                  << alignment.get_sequence() << "\t"
                                  << label << std::endl;
                    }
                };

            if (config->call_bubbles) {
                annotated_graph_algorithm::call_bubbles(
                    *masked_graph,
                    *anno_graph,
                    print_variant,
                    &thread_pool
                );
            } else if (config->call_breakpoints) {
                annotated_graph_algorithm::call_breakpoints(
                    *masked_graph,
                    *anno_graph,
                    print_variant,
                    &thread_pool
                );
            } else {
                std::cerr << "ERROR: no variant calling mode selected. Exiting" << std::endl;
                exit(1);
            }

            thread_pool.join();

            if (config->verbose) {
                std::cout << "# nodes checked: " << masked_graph->num_nodes()
                          << std::endl
                          << "# called: " << num_calls
                          << std::endl;
            }

            return 0;
        }
        case Config::PARSE_TAXONOMY: {
            TaxIDMapper taxid_mapper;
            if (config->accession2taxid.length()
                && !taxid_mapper.parse_accession2taxid(config->accession2taxid)) {
                std::cerr << "ERROR: failed to read accession2taxid file" << std::endl;
                exit(1);
            }

            if (config->taxonomy_nodes.length()
                && !taxid_mapper.parse_nodes(config->taxonomy_nodes)) {
                std::cerr << "ERROR: failed to read nodes.dmp file" << std::endl;
                exit(1);
            }

            std::ofstream out(config->outfbase + ".taxonomy.map", std::ios::binary);
            taxid_mapper.serialize(out);
            return 0;
        }
        case Config::NO_IDENTITY: {
            assert(false);
            break;
        }
    }

    return 0;
}
