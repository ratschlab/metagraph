#include <filesystem>
#include <typeinfo>

#include <fmt/format.h>
#include <ips4o.hpp>
#include <json/json.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <tsl/hopscotch_map.h>

#include "common/logger.hpp"

#include "common/unix_tools.hpp"
#include "cli/config/config.hpp"
#include "seq_io/sequence_io.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "graph/representation/succinct/boss_chunk.hpp"
#include "graph/representation/succinct/boss_merge.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/annotation_converters.hpp"
#include "common/serialization.hpp"
#include "common/algorithms.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "common/threads/threading.hpp"
#include "seq_io/kmc_parser.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/bitmap/dbg_bitmap_construct.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/alignment/aligner_methods.hpp"
#include "common/network/server.hpp"
#include "graph/graph_extensions/node_weights.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "common/taxid_mapper.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/utils/template_utils.hpp"
#include "common/vectors/int_vector_algorithm.hpp"
#include "seq_io/formats.hpp"
#include "cli/build.hpp"
#include "cli/augment.hpp"
#include "cli/clean.hpp"
#include "cli/load/load_graph.hpp"
#include "cli/load/load_annotation.hpp"

using mg::common::logger;
using utils::get_verbose;
using namespace mg::bitmap_graph;
using namespace mg::succinct;

typedef annotate::MultiLabelEncoded<std::string> Annotator;

static const size_t kRowBatchSize = 100'000;
const bool kPrefilterWithBloom = false;


void annotate_data(const std::vector<std::string> &files,
                   const std::string &ref_sequence_path,
                   AnnotatedDBG *anno_graph,
                   bool forward_and_reverse,
                   size_t min_count,
                   size_t max_count,
                   bool filename_anno,
                   bool annotate_sequence_headers,
                   const std::string &fasta_anno_comment_delim,
                   const std::string &fasta_header_delimiter,
                   const std::vector<std::string> &anno_labels) {
    size_t total_seqs = 0;

    Timer timer;

    ThreadPool thread_pool(get_num_threads() > 1 ? get_num_threads() : 0);

    // iterate over input files
    for (const auto &file : files) {
        Timer data_reading_timer;

        logger->trace("Parsing '{}'", file);

        std::vector<std::string> labels = anno_labels;
        if (filename_anno) {
            labels.push_back(file);
        }
        // remember the number of base labels to remove those unique to each sequence quickly
        const size_t num_base_labels = labels.size();

        if (file_format(file) == "VCF") {
            read_vcf_file_with_annotations_critical(
                file,
                ref_sequence_path,
                anno_graph->get_graph().get_k(),
                [&](auto&& seq, const auto &variant_labels) {
                    labels.insert(labels.end(),
                                  variant_labels.begin(), variant_labels.end());

                    thread_pool.enqueue(
                        [anno_graph](const std::string &sequence, const auto &labels) {
                            anno_graph->annotate_sequence(sequence, labels);
                        },
                        std::move(seq),
                        labels
                    );

                    total_seqs += 1;

                    if (logger->level() <= spdlog::level::level_enum::trace
                                                    && total_seqs % 10000 == 0) {
                        logger->trace(
                            "processed {} variants, last was annotated as <{}>, {} sec",
                            total_seqs, fmt::join(labels, "><"), timer.elapsed());
                    }

                    labels.resize(num_base_labels);
                },
                forward_and_reverse
            );
        } else if (file_format(file) == "KMC") {
            kmc::read_kmers(
                file,
                [&](std::string&& sequence) {
                    thread_pool.enqueue(
                        [anno_graph](const std::string &sequence, const auto &labels) {
                            anno_graph->annotate_sequence(sequence, labels);
                        },
                        std::move(sequence),
                        labels
                    );

                    total_seqs += 1;

                    if (logger->level() <= spdlog::level::level_enum::trace
                                                    && total_seqs % 10000 == 0) {
                        logger->trace(
                            "processed {} sequences, trying to annotate as <{}>, {} sec",
                            total_seqs, fmt::join(labels, "><"), timer.elapsed());
                    }
                },
                !anno_graph->get_graph().is_canonical_mode(),
                min_count,
                max_count
            );
        } else if (file_format(file) == "FASTA"
                    || file_format(file) == "FASTQ") {
            read_fasta_file_critical(
                file,
                [&](kseq_t *read_stream) {
                    // add sequence header to labels
                    if (annotate_sequence_headers) {
                        for (const auto &label
                                : utils::split_string(fasta_anno_comment_delim != Config::UNINITIALIZED_STR
                                                        ? utils::join_strings(
                                                            { read_stream->name.s, read_stream->comment.s },
                                                            fasta_anno_comment_delim,
                                                            true)
                                                        : read_stream->name.s,
                                                      fasta_header_delimiter)) {
                            labels.push_back(label);
                        }
                    }

                    thread_pool.enqueue(
                        [anno_graph](const std::string &sequence, const auto &labels) {
                            anno_graph->annotate_sequence(sequence, labels);
                        },
                        std::string(read_stream->seq.s),
                        labels
                    );

                    total_seqs += 1;

                    if (logger->level() <= spdlog::level::level_enum::trace
                                                    && total_seqs % 10000 == 0) {
                        logger->trace("processed {} sequences, last was {}, trying to annotate as <{}>, {} sec",
                                      total_seqs, read_stream->name.s, fmt::join(labels, "><"), timer.elapsed());
                    }

                    labels.resize(num_base_labels);
                },
                forward_and_reverse
            );
        } else {
            logger->error("Unknown filetype for file '{}'", file);
            exit(1);
        }

        logger->trace("File '{}' processed in {} sec, current mem usage: {} MiB, total time {} sec",
                      file, data_reading_timer.elapsed(), get_curr_RSS() >> 20, timer.elapsed());
    }

    thread_pool.join();
}


void annotate_coordinates(const std::vector<std::string> &files,
                          AnnotatedDBG *anno_graph,
                          bool forward_and_reverse,
                          size_t genome_bin_size) {
    size_t total_seqs = 0;

    Timer timer;

    const size_t k = anno_graph->get_graph().get_k();

    ThreadPool thread_pool(get_num_threads() > 1 ? get_num_threads() : 0);

    // iterate over input files
    for (const auto &file : files) {
        Timer data_reading_timer;

        logger->trace("Parsing '{}'", file);

        // open stream
        if (file_format(file) == "FASTA"
                    || file_format(file) == "FASTQ") {

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
                        thread_pool.enqueue(
                            [anno_graph](const std::string &sequence, const auto &labels) {
                                anno_graph->annotate_sequence(
                                    sequence, { utils::join_strings(labels, "\1"), }
                                );
                            },
                            forward_strand
                                ? sequence.substr(i, bin_size)
                                : sequence.substr(sequence.size() - i - bin_size, bin_size),
                            labels
                        );
                    }

                    total_seqs += 1;

                    if (logger->level() <= spdlog::level::level_enum::trace
                                                    && total_seqs % 10000 == 0) {
                        logger->trace("processed {} sequences, last was {}, trying to annotate as <{}>, {} sec",
                                      total_seqs, read_stream->name.s, fmt::join(labels, "><"), timer.elapsed());
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
            logger->error("The type of file '{}' is not supported", file);
            exit(1);
        }

        logger->trace("File '{}' processed in {} sec, current mem usage: {} MiB, total time {} sec",
                      file, data_reading_timer.elapsed(), get_curr_RSS() >> 20, timer.elapsed());
    }

    thread_pool.join();
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

        std::transform(alignments.begin(), alignments.end(),
                       std::back_inserter(sequences),
                       [](const auto &alignment) { return alignment.get_sequence(); });

        weights = alignments.get_alignment_weights(aligner->get_config());

        if (sequences.empty())
            aligner = nullptr;
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

std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(std::shared_ptr<DeBruijnGraph> graph,
                                                       const Config &config) {
    auto annotation_temp = config.infbase_annotators.size()
            ? initialize_annotation(parse_annotation_type(config.infbase_annotators.at(0)), config, 0)
            : initialize_annotation(config.anno_type, config, graph->max_index());

    if (config.infbase_annotators.size()
            && !annotation_temp->load(config.infbase_annotators.at(0))) {
        logger->error("Cannot load annotations for graph {}, file corrupted",
                      config.infbase);
        exit(1);
    }

    // load graph
    auto anno_graph = std::make_unique<AnnotatedDBG>(std::move(graph),
                                                     std::move(annotation_temp));

    if (!anno_graph->check_compatibility()) {
        logger->error("Graph and annotation are not compatible");
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
                           if (!exists)
                               logger->trace("Removing mask-in label {}", label);

                           return !exists;
                       }),
        config->label_mask_in.end()
    );

    config->label_mask_out.erase(
        std::remove_if(config->label_mask_out.begin(),
                       config->label_mask_out.end(),
                       [&](const auto &label) {
                           bool exists = anno_graph.label_exists(label);
                           if (!exists)
                               logger->trace("Removing mask-out label {}", label);

                           return !exists;
                       }),
        config->label_mask_out.end()
    );

    logger->trace("Masked in: {}", fmt::join(config->label_mask_in, " "));
    logger->trace("Masked out: {}", fmt::join(config->label_mask_out, " "));

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
    logger->trace("Converting annotation to {}...",
                  Config::annotype_to_string(config.anno_type));

    auto target_annotator = annotate::convert<AnnotatorTo>(std::move(*annotator));
    annotator.reset();
    logger->trace("Conversion done in {} sec", timer.elapsed());

    logger->trace("Serializing annotation to '{}'...", config.outfbase);
    target_annotator->serialize(config.outfbase);
    logger->trace("Serialization done in {} sec", timer.elapsed());
}


DBGAlignerConfig initialize_aligner_config(const DeBruijnGraph &graph, Config &config) {
    // fix seed length bounds
    if (!config.alignment_min_seed_length || config.alignment_seed_unimems)
        config.alignment_min_seed_length = graph.get_k();

    if (config.alignment_max_seed_length == std::numeric_limits<size_t>::max()
            && !config.alignment_seed_unimems)
        config.alignment_max_seed_length = graph.get_k();

    DBGAlignerConfig aligner_config;

    aligner_config.queue_size = config.alignment_queue_size;
    aligner_config.bandwidth = config.alignment_vertical_bandwidth;
    aligner_config.num_alternative_paths = config.alignment_num_alternative_paths;
    aligner_config.min_seed_length = config.alignment_min_seed_length;
    aligner_config.max_seed_length = config.alignment_max_seed_length;
    aligner_config.max_num_seeds_per_locus = config.alignment_max_num_seeds_per_locus;
    aligner_config.min_cell_score = config.alignment_min_cell_score;
    aligner_config.min_path_score = config.alignment_min_path_score;
    aligner_config.gap_opening_penalty = -config.alignment_gap_opening_penalty;
    aligner_config.gap_extension_penalty = -config.alignment_gap_extension_penalty;
    aligner_config.forward_and_reverse_complement = config.forward_and_reverse;
    aligner_config.alignment_edit_distance = config.alignment_edit_distance;
    aligner_config.alignment_match_score = config.alignment_match_score;
    aligner_config.alignment_mm_transition_score = config.alignment_mm_transition_score;
    aligner_config.alignment_mm_transversion_score = config.alignment_mm_transversion_score;

    logger->trace("Alignment settings:");
    logger->trace("\t Seeding: {}", (config.alignment_seed_unimems ? "unimems" : "nodes"));
    logger->trace("\t Alignments to report: {}", config.alignment_num_alternative_paths);
    logger->trace("\t Priority queue size: {}", config.alignment_queue_size);
    logger->trace("\t Min seed length: {}", config.alignment_min_seed_length);
    logger->trace("\t Max seed length: {}", config.alignment_max_seed_length);
    logger->trace("\t Max num seeds per locus: {}", config.alignment_max_num_seeds_per_locus);
    logger->trace("\t Gap opening penalty: {}",
                  int64_t(config.alignment_gap_opening_penalty));
    logger->trace("\t Gap extension penalty: {}",
                  int64_t(config.alignment_gap_extension_penalty));
    logger->trace("\t Min DP table cell score: {}", int64_t(config.alignment_min_cell_score));
    logger->trace("\t Min alignment score: {}", config.alignment_min_path_score);

    logger->trace("\t Scoring matrix: {}", config.alignment_edit_distance ? "unit costs" : "matrix");
    if (!config.alignment_edit_distance) {
        logger->trace("\t\t Match score: {}", int64_t(config.alignment_match_score));
        logger->trace("\t\t (DNA) Transition score: {}",
                      int64_t(config.alignment_mm_transition_score));
        logger->trace("\t\t (DNA) Transversion score: {}",
                      int64_t(config.alignment_mm_transversion_score));
    }

    aligner_config.set_scoring_matrix();

    return aligner_config;
}

std::unique_ptr<IDBGAligner> build_aligner(const DeBruijnGraph &graph, Config &config) {
    DBGAlignerConfig aligner_config = initialize_aligner_config(graph, config);

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
        return std::make_unique<DBGAligner<UniMEMSeeder<>>>(graph, aligner_config);

    } else if (config.alignment_min_seed_length < graph.get_k()) {
        if (!dynamic_cast<const DBGSuccinct*>(&graph)) {
            logger->error("SuffixSeeder can be used only with succinct graph representation");
            exit(1);
        }

        // Use the seeder that seeds to node suffixes
        return std::make_unique<DBGAligner<SuffixSeeder<>>>(graph, aligner_config);

    } else {
        return std::make_unique<DBGAligner<>>(graph, aligner_config);
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

    read_fasta_file_critical(file, [&](kseq_t *read_stream) {
        if (get_verbose())
            std::cout << "Sequence: " << read_stream->seq.s << "\n";

        if (config.query_presence
                && config.alignment_length == graph.get_k()) {

            bool found = graph.find(read_stream->seq.s,
                                    config.discovery_fraction);

            if (!config.filter_present) {
                std::cout << found << "\n";

            } else if (found) {
                std::cout << ">" << read_stream->name.s << "\n"
                                 << read_stream->seq.s << "\n";
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
                    std::string_view(read_stream->seq.s + i, config.alignment_length),
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
                std::cout << std::string_view(read_stream->seq.s + i, config.alignment_length)
                          << ": " << graphindices[i] << "\n";
            }
        } else {
            // map input subsequences to multiple nodes
            for (size_t i = 0; i + graph.get_k() <= read_stream->seq.l; ++i) {
                // TODO: make more efficient
                std::string_view subseq(read_stream->seq.s + i, config.alignment_length);

                dbg->call_nodes_with_suffix(subseq,
                                            [&](auto node, auto) {
                                                std::cout << subseq << ": "
                                                          << node
                                                          << "\n";
                                            },
                                            config.alignment_length);
            }
        }

    }, config.forward_and_reverse);

    logger->trace("File '{}' processed in {} sec, current mem usage: {} MiB, total time {} sec",
                  file, data_reading_timer.elapsed(), get_curr_RSS() >> 20, timer.elapsed());
}

typedef std::function<void(const std::string&)> SequenceCallback;

std::unique_ptr<AnnotatedDBG>
construct_query_graph(const AnnotatedDBG &anno_graph,
                      std::function<void(SequenceCallback)> call_sequences,
                      double discovery_fraction,
                      size_t num_threads) {
    const auto *full_dbg = &anno_graph.get_graph();
    if (!full_dbg)
        throw std::runtime_error("Error: batch queries are supported only for de Bruijn graphs");

    const auto &full_annotation = anno_graph.get_annotation();

    Timer timer;

    // construct graph storing all k-mers in query
    auto graph = std::make_shared<DBGHashOrdered>(full_dbg->get_k(), false);

    const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(full_dbg);
    if (kPrefilterWithBloom && dbg_succ) {
        if (dbg_succ->get_bloom_filter())
            logger->trace("[Query graph construction] Started indexing k-mers pre-filtered with Bloom filter");

        call_sequences([&graph,&dbg_succ](const std::string &sequence) {
            graph->add_sequence(sequence, get_missing_kmer_skipper(
                dbg_succ->get_bloom_filter(),
                sequence
            ));
        });
    } else {
        call_sequences([&graph](const std::string &sequence) {
            graph->add_sequence(sequence);
        });
    }

    logger->trace("[Query graph construction] k-mer indexing took {} sec", timer.elapsed());
    timer.reset();

    // pull contigs from query graph
    std::vector<std::pair<std::string, std::vector<DeBruijnGraph::node_index>>> contigs;
    graph->call_sequences(
        [&](const std::string &contig, const auto &path) { contigs.emplace_back(contig, path); },
        full_dbg->is_canonical_mode()
    );

    logger->trace("[Query graph construction] Contig extraction took {} sec", timer.elapsed());
    timer.reset();

    if (full_dbg->is_canonical_mode()) {
        // construct graph storing all distinct k-mers in query
        graph = std::make_shared<DBGHashOrdered>(full_dbg->get_k(), true);

        for (const auto &pair : contigs) {
            graph->add_sequence(pair.first);
        }

        logger->trace("[Query graph construction] k-mers reindexed in canonical mode in {} sec",
                      timer.elapsed());
        timer.reset();
    }

    // map contigs onto the full graph
    auto index_in_full_graph
        = std::make_shared<std::vector<uint64_t>>(graph->max_index() + 1, 0);

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

    logger->trace("[Query graph construction] Contigs mapped to graph in {} sec", timer.elapsed());
    timer.reset();

    contigs.clear();

    assert(!(*index_in_full_graph)[0]);

    if (discovery_fraction > 0) {
        sdsl::bit_vector mask(graph->max_index() + 1, false);

        call_sequences([&](const std::string &sequence) {
            if (sequence.length() < graph->get_k())
                return;

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

        logger->trace("[Query graph construction] Reduced k-mer dictionary in {} sec",
                      timer.elapsed());
        timer.reset();
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
                          utils::LessFirst(), num_threads);

    // initialize fast query annotation
    // copy annotations from the full graph to the query graph
    auto annotation = std::make_unique<annotate::RowCompressed<>>(
        graph->max_index(),
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

    logger->trace("[Query graph construction] Query annotation constructed in {} sec",
                  timer.elapsed());
    timer.reset();

    auto masked_graph = std::make_shared<MaskedDeBruijnGraph>(graph,
        [=](auto i) -> bool { return (*index_in_full_graph)[i]; }
    );

    // build annotated graph from the query graph and copied annotations
    return std::make_unique<AnnotatedDBG>(masked_graph, std::move(annotation));
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
                logger->error("Node weights are not compatible with graph");
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

        if (get_verbose()) {
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

template <class KmerHasher>
void print_bloom_filter_stats(const KmerBloomFilter<KmerHasher> *kmer_bloom) {
    if (!kmer_bloom)
        return;

    std::cout << "====================== BLOOM STATS =====================" << std::endl;
    std::cout << "Size (bits):\t" << kmer_bloom->size() << std::endl
              << "Num hashes:\t" << kmer_bloom->num_hash_functions() << std::endl;
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

    } else if (dynamic_cast<const annotate::MultiBRWTAnnotator *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::BRWT) << std::endl;
        const auto &brwt = dynamic_cast<const annotate::MultiBRWTAnnotator &>(annotation).get_matrix();
        std::cout << "=================== Multi-BRWT STATS ===================" << std::endl;
        std::cout << "num nodes: " << brwt.num_nodes() << std::endl;
        std::cout << "avg arity: " << brwt.avg_arity() << std::endl;
        std::cout << "shrinkage: " << brwt.shrinking_rate() << std::endl;
        if (get_verbose()) {
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
                logger->error("Bad json file:\n{}", errors);
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
            logger->error("No input sequences received from client");
            // TODO: no input sequences -> form an error message for the client
            throw std::domain_error("No input sequences");
        }

        return oss.str();

    } catch (const Json::LogicError &e) {
        logger->error("Bad json file: {}", e.what());
        //TODO: send errors in a json file
        throw;
    } catch (const std::exception &e) {
        logger->error("Processing request error: {}", e.what());
        //TODO: send errors in a json file
        throw;
    } catch (...) {
        logger->error("Processing request error");
        //TODO: send errors in a json file
        throw;
    }
}


int main(int argc, char *argv[]) {
    auto config = std::make_unique<Config>(argc, argv);

    logger->set_level(get_verbose() ? spdlog::level::trace : spdlog::level::info);
    //logger->set_pattern("%^date %x....%$  %v");
    //spdlog::set_pattern("[%H:%M:%S %z] [%n] [%^---%L---%$] [thread %t] %v");
    //console_sink->set_color(spdlog::level::trace, "\033[37m");
    spdlog::flush_every(std::chrono::seconds(1));

    logger->trace("Metagraph started");

    const auto &files = config->fnames;

    switch (config->identity) {
        case Config::EXPERIMENT: {
            break;
        }
        case Config::BUILD: {
            return build_graph(config.get());
        }
        case Config::EXTEND: {
            return augment_graph(config.get());
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
                              config->annotate_sequence_headers,
                              config->fasta_anno_comment_delim,
                              config->fasta_header_delimiter,
                              config->anno_labels);

                anno_graph->get_annotation().serialize(config->outfbase);

            } else {
                // |config->separately| is true

                size_t num_threads = 1;
                if (!config->files_sequentially) {
                    // annotate multiple files in parallel, each in a single thread
                    num_threads = get_num_threads();
                    set_num_threads(1);
                }

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
                                  config->annotate_sequence_headers,
                                  config->fasta_anno_comment_delim,
                                  config->fasta_header_delimiter,
                                  config->anno_labels);

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
                = std::make_unique<annotate::RowCompressed<>>(graph_temp->max_index());

            if (config->infbase_annotators.size()
                    && !annotation_temp->load(config->infbase_annotators.at(0))) {
                logger->error("Cannot load annotations from '{}'",
                              config->infbase_annotators.at(0));
                exit(1);
            }

            // load graph
            AnnotatedDBG anno_graph(graph_temp,
                                    std::move(annotation_temp),
                                    config->fast);

            if (!anno_graph.check_compatibility()) {
                logger->error("Graph and annotation are incompatible");
                exit(1);
            }

            annotate_coordinates(files,
                                 &anno_graph,
                                 config->forward_and_reverse,
                                 config->genome_binsize_anno);

            anno_graph.get_annotation().serialize(config->outfbase);

            return 0;
        }
        case Config::MERGE_ANNOTATIONS: {
            if (config->anno_type == Config::ColumnCompressed) {
                annotate::ColumnCompressed<> annotation(0, config->num_columns_cached, get_verbose());
                if (!annotation.merge_load(files)) {
                    logger->error("Cannot load annotations");
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
                        logger->error("Cannot load annotations from file '{}'", filename);
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
                annotate::merge<annotate::MultiBRWTAnnotator>(std::move(annotators), stream_files, config->outfbase);
            } else {
                logger->error("Merging of annotations to '{}' representation is not implemented",
                              config->annotype_to_string(config->anno_type));
                exit(1);
            }

            return 0;
        }
        case Config::QUERY: {
            assert(config->infbase_annotators.size() == 1);

            auto graph = load_critical_dbg(config->infbase);
            auto anno_graph = initialize_annotated_dbg(graph, *config);

            ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1);

            Timer timer;

            std::unique_ptr<IDBGAligner> aligner;
            // TODO: make aligner work with batch querying
            if (config->align_sequences && !config->fast)
                aligner.reset(build_aligner(*graph, *config).release());

            // iterate over input files
            for (const auto &file : files) {
                logger->trace("Parsing sequences from file '{}'", file);

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
                        get_num_threads()
                    );

                    graph_to_query = query_graph.get();

                    logger->trace("Query graph constructed for '{}' in {} sec",
                                  file, curr_timer.elapsed());
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

                logger->trace("File '{}' was processed in {} sec, total time: {}", file,
                              curr_timer.elapsed(), timer.elapsed());
            }

            return 0;
        }
        case Config::SERVER_QUERY: {
            assert(config->infbase_annotators.size() == 1);

            Timer timer;

            logger->info("Loading graph...");

            auto graph = load_critical_dbg(config->infbase);
            auto anno_graph = initialize_annotated_dbg(graph, *config);

            logger->info("Graph loaded in {} sec, current mem usage: {} MiB",
                         timer.elapsed(), get_curr_RSS() >> 20);

            std::unique_ptr<IDBGAligner> aligner;
            // TODO: make aligner work with batch querying
            if (config->align_sequences && !config->fast)
                aligner.reset(build_aligner(*graph, *config).release());

            const size_t num_threads = std::max(1u, get_num_threads());

            logger->info("Initializing a TCP service with {} threads"
                         ", listening on port {}", num_threads, config->port);

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
                for (size_t i = 0; i < std::max(1u, get_num_threads()); ++i) {
                    workers.emplace_back([&io_context]() { io_context.run(); });
                }
                for (auto &thread : workers) {
                    thread.join();
                }
            } catch (const std::exception &e) {
                logger->error("Exception: {}", e.what());
            } catch (...) {
                logger->error("Unknown exception");
            }

            return 0;
        }
        case Config::COMPARE: {
            assert(files.size());

            logger->info("Loading graph                '{}'", files.at(0));
            auto graph = load_critical_dbg(files.at(0));

            for (size_t f = 1; f < files.size(); ++f) {
                logger->info("Loading graph for comparison '{}'", files[f]);
                auto second = load_critical_dbg(files[f]);
                if (*graph == *second) {
                    logger->info("Graphs are identical");
                } else {
                    logger->info("Graphs are not identical");
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
                logger->error("No input files provided, nothing to concatenate");
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
                    auto p = BOSS::Chunk::build_boss_from_chunks(chunk_files, get_verbose());
                    auto dbg_succ = std::make_unique<DBGSuccinct>(p.first, p.second);

                    logger->trace("Chunks concatenated in {} sec", timer.elapsed());

                    if (config->clear_dummy) {
                        logger->trace("Traverse source dummy edges,"
                                      " remove redundant ones, and mark"
                                      " those that cannot be removed");
                        dbg_succ->mask_dummy_kmers(get_num_threads(), true);
                    }
                    graph = std::move(dbg_succ);
                    break;
                }
                case Config::GraphType::BITMAP: {
                    graph.reset(DBGBitmapConstructor::build_graph_from_chunks(
                        chunk_files, config->canonical, get_verbose()
                    ));
                    break;
                }
                default:
                    logger->error("Cannot concatenate chunks for this graph representation");
                    exit(1);
            }
            assert(graph);

            logger->trace("Graph was assembled in {} sec", timer.elapsed());

            if (logger->level() <= spdlog::level::level_enum::trace) {
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
                logger->info("Opening file '{}'", file);

                dbg_graphs.emplace_back(load_critical_graph_from_file<DBGSuccinct>(file));

                graphs.push_back(&dbg_graphs.back()->get_boss());

                if (get_verbose())
                    print_boss_stats(*graphs.back());

                config->canonical &= dbg_graphs.back()->is_canonical_mode();
            }

            logger->info("Graphs are loaded in {} sec", timer.elapsed());

            if (config->dynamic) {
                logger->info("Start merging traversal");
                timer.reset();

                graph = dbg_graphs.at(0)->release_boss();

                if (graph->get_state() != BOSS::State::DYN) {
                    logger->trace("Switching state of succinct graph to dynamic...");

                    graph->switch_state(BOSS::State::DYN);

                    logger->trace("Switching done in {} sec", timer.elapsed());
                }

                for (size_t i = 1; i < graphs.size(); ++i) {
                    graph->merge(dbg_graphs.at(i)->get_boss());

                    logger->info("Graph '{}' merged in {} sec", files[i], timer.elapsed());

                    dbg_graphs.at(i).reset();
                }
            } else if (get_num_threads() > 1 || config->parts_total > 1) {
                logger->info("Start merging blocks");
                timer.reset();

                auto *chunk = merge::merge_blocks_to_chunk(
                    graphs,
                    config->part_idx,
                    config->parts_total,
                    get_num_threads(),
                    config->num_bins_per_thread,
                    get_verbose()
                );
                if (!chunk) {
                    logger->error("ERROR when building chunk {}", config->part_idx);
                    exit(1);
                }
                logger->info("Blocks merged in {} sec", timer.elapsed());

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
                logger->info("Start merging graphs");
                timer.reset();

                graph = merge::merge(graphs, get_verbose());
            }
            dbg_graphs.clear();

            assert(graph);

            logger->info("Graphs merged in {} sec", timer.elapsed());

            // graph output
            DBGSuccinct(graph, config->canonical).serialize(config->outfbase);

            return 0;
        }
        case Config::CLEAN: {
            return clean_graph(config.get());
        }
        case Config::STATS: {
            for (const auto &file : files) {
                std::shared_ptr<DeBruijnGraph> graph;

                graph = load_critical_dbg(file);
                graph->load_extension<NodeWeights>(file);

                logger->info("Statistics for graph '{}'", file);

                print_stats(*graph);

                if (auto dbg_succ = dynamic_cast<DBGSuccinct*>(graph.get())) {
                    const auto &boss_graph = dbg_succ->get_boss();

                    print_boss_stats(boss_graph,
                                     config->count_dummy,
                                     get_num_threads(),
                                     get_verbose());

                    if (config->print_graph_internal_repr) {
                        logger->info("Printing internal representation");
                        boss_graph.print_internal_representation(std::cout);
                    }
                    print_bloom_filter_stats(dbg_succ->get_bloom_filter());
                }

                if (config->print_graph)
                    std::cout << *graph;
            }

            for (const auto &file : config->infbase_annotators) {
                auto annotation = initialize_annotation(file, *config);

                if (config->print_column_names) {
                    annotate::LabelEncoder<std::string> label_encoder;

                    logger->info("Scanning annotation '{}'", file);

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
                        logger->error("Cannot read label encoder from file '{}'", file);
                        exit(1);
                    }

                    std::cout << "Number of columns: " << label_encoder.size() << std::endl;
                    for (size_t c = 0; c < label_encoder.size(); ++c) {
                        std::cout << label_encoder.decode(c) << '\n';
                    }

                    continue;
                }

                if (!annotation->load(file)) {
                    logger->error("Cannot load annotations from file '{}'", file);
                    exit(1);
                }

                logger->info("Statistics for annotation '{}'", file);
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

            if (config->dump_text_anno) {
                const Config::AnnotationType input_anno_type
                    = parse_annotation_type(files.at(0));

                auto annotation = initialize_annotation(files.at(0), *config);

                logger->trace("Loading annotation...");

                if (config->anno_type == Config::ColumnCompressed) {
                    if (!annotation->merge_load(files)) {
                        logger->error("Cannot load annotations");
                        exit(1);
                    }
                } else {
                    // Load annotation from disk
                    if (!annotation->load(files.at(0))) {
                        logger->error("Cannot load annotations from file '{}'", files.at(0));
                        exit(1);
                    }
                }

                logger->trace("Annotation loaded in {} sec", timer.elapsed());
                logger->trace("Dumping annotators...\t");

                if (input_anno_type == Config::ColumnCompressed) {
                    assert(dynamic_cast<annotate::ColumnCompressed<>*>(annotation.get()));
                    dynamic_cast<annotate::ColumnCompressed<>*>(
                        annotation.get()
                    )->dump_columns(config->outfbase, get_num_threads());
                } else if (input_anno_type == Config::BRWT) {
                    assert(dynamic_cast<annotate::MultiBRWTAnnotator*>(annotation.get()));
                    dynamic_cast<annotate::MultiBRWTAnnotator*>(
                        annotation.get()
                    )->dump_columns(config->outfbase, get_num_threads());
                } else {
                    throw std::runtime_error("Dumping columns for this type not implemented");
                }

                logger->trace("Dumping done in {} sec", timer.elapsed());

                return 0;
            }

            /********************************************************/
            /***************** rename column labels *****************/
            /********************************************************/

            if (config->rename_instructions_file.size()) {
                tsl::hopscotch_map<std::string, std::string> dict;
                std::ifstream instream(config->rename_instructions_file);
                if (!instream.is_open()) {
                    logger->error("Cannot open file '{}'", config->rename_instructions_file);
                    exit(1);
                }
                std::string old_name;
                std::string new_name;
                while (instream.good() && !(instream >> old_name).eof()) {
                    instream >> new_name;
                    if (instream.fail() || instream.eof()) {
                        logger->error("Wrong format of the rules for renaming"
                                      " annotation columns passed in file '{}'",
                                      config->rename_instructions_file);
                        exit(1);
                    }
                    dict[old_name] = new_name;
                }

                auto annotation = initialize_annotation(files.at(0), *config);

                logger->trace("Loading annotation...");

                // TODO: rename columns without loading the full annotation
                if (config->anno_type == Config::ColumnCompressed) {
                    if (!annotation->merge_load(files)) {
                        logger->error("Cannot load annotations");
                        exit(1);
                    } else {
                        logger->info("Annotation #objects: {}\t#labels: {}",
                                     annotation->num_objects(), annotation->num_labels());
                    }
                } else {
                    // Load annotation from disk
                    if (!annotation->load(files.at(0))) {
                        logger->error("Cannot load annotations from file '{}'", files.at(0));
                        exit(1);
                    }
                }

                logger->trace("Annotation loaded in {} sec", timer.elapsed());
                logger->trace("Renaming...");

                //TODO: could be made to work with streaming
                annotation->rename_labels(dict);

                annotation->serialize(config->outfbase);
                logger->trace("Renaming done in {} sec", timer.elapsed());

                return 0;
            }

            /********************************************************/
            /****************** convert annotation ******************/
            /********************************************************/

            const Config::AnnotationType input_anno_type
                = parse_annotation_type(files.at(0));

            if (config->anno_type == input_anno_type) {
                logger->info("Skipping conversion: same input and target type: {}",
                              Config::annotype_to_string(config->anno_type));
                return 0;
            }

            if (input_anno_type == Config::ColumnCompressed && files.size() > 1) {
                logger->error("Conversion of multiple annotators only "
                              "supported for ColumnCompressed");
                exit(1);
            }

            logger->trace("Converting to {} annotator...",
                          Config::annotype_to_string(config->anno_type));

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
                        logger->error("Streaming conversion from RowCompressed "
                                      "annotation is not implemented for the requested "
                                      "target type: {}",
                                      Config::annotype_to_string(config->anno_type));
                        exit(1);
                }

                logger->trace("Annotation converted in {} sec", timer.elapsed());

                logger->trace("Serializing to '{}'...", config->outfbase);

                target_annotator->serialize(config->outfbase);

                logger->trace("Serialization done in {} sec", timer.elapsed());
            } else if (input_anno_type == Config::ColumnCompressed) {
                auto annotation = initialize_annotation(files.at(0), *config);

                logger->trace("Loading annotation...");

                // Load annotation from disk
                if (!annotation->merge_load(files)) {
                    logger->error("Cannot load annotations");
                    exit(1);
                }

                logger->trace("Annotation loaded in {} sec", timer.elapsed());

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
                            annotate::RowCompressed<> row_annotator(annotator->num_objects());
                            convert_to_row_annotator(*annotator,
                                                     &row_annotator,
                                                     get_num_threads());
                            annotator.reset();

                            logger->trace("Annotation converted in {} sec", timer.elapsed());
                            logger->trace("Serializing to '{}'...", config->outfbase);

                            row_annotator.serialize(config->outfbase);

                            logger->trace("Serialization done in {} sec", timer.elapsed());

                        } else {
                            convert_to_row_annotator(*annotator,
                                                     config->outfbase,
                                                     get_num_threads());
                            logger->trace("Annotation converted and serialized in {} sec",
                                          timer.elapsed());
                        }
                        break;
                    }
                    case Config::BRWT: {
                        auto brwt_annotator = config->greedy_brwt
                            ? annotate::convert_to_greedy_BRWT<annotate::MultiBRWTAnnotator>(
                                std::move(*annotator),
                                config->parallel_nodes,
                                get_num_threads())
                            : annotate::convert_to_simple_BRWT<annotate::MultiBRWTAnnotator>(
                                std::move(*annotator),
                                config->arity_brwt,
                                config->parallel_nodes,
                                get_num_threads());

                        annotator.reset();
                        logger->trace("Annotation converted in {} sec", timer.elapsed());

                        logger->trace("Serializing to '{}'", config->outfbase);

                        brwt_annotator->serialize(config->outfbase);

                        logger->trace("Serialization done in {} sec", timer.elapsed());
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
                logger->error("Conversion to other representations"
                              " is not implemented for {} annotator",
                              Config::annotype_to_string(input_anno_type));
                exit(1);
            }

            return 0;
        }
        case Config::TRANSFORM: {
            assert(files.size() == 1);
            assert(config->outfbase.size());

            if (config->initialize_bloom) {
                std::filesystem::remove(
                    utils::remove_suffix(config->outfbase, ".bloom") + ".bloom"
                );
            }

            Timer timer;
            logger->trace("Graph loading...");

            auto graph = load_critical_dbg(files.at(0));

            logger->trace("Graph loaded in {} sec", timer.elapsed());

            auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(graph);

            if (!dbg_succ.get())
                throw std::runtime_error("Only implemented for DBGSuccinct");

            if (config->initialize_bloom) {
                assert(config->bloom_fpp > 0.0 && config->bloom_fpp <= 1.0);
                assert(config->bloom_bpk >= 0.0);
                assert(config->bloom_fpp < 1.0 || config->bloom_bpk > 0.0);

                logger->trace("Construct Bloom filter for nodes...");

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

                logger->trace("Bloom filter constructed in {} sec", timer.elapsed());

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
                logger->trace("Traverse the tree of source dummy edges and remove redundant ones...");
                timer.reset();

                // remove redundant dummy edges and mark all other dummy edges
                dbg_succ->mask_dummy_kmers(get_num_threads(), true);

                logger->trace("The tree of source dummy edges traversed in {} sec", timer.elapsed());
                timer.reset();
            }

            if (config->to_adj_list) {
                logger->trace("Converting graph to adjacency list...");

                auto *boss = &dbg_succ->get_boss();
                timer.reset();

                std::ofstream outstream(config->outfbase + ".adjlist");
                boss->print_adj_list(outstream);

                logger->trace("Conversion done in {} sec", timer.elapsed());

                return 0;
            }

            logger->trace("Converting graph to state {}",
                          Config::state_to_string(config->state));
            timer.reset();

            dbg_succ->switch_state(config->state);

            logger->trace("Conversion done in {} sec", timer.elapsed());

            logger->trace("Serializing transformed graph...");
            dbg_succ->serialize(config->outfbase);
            logger->trace("Serialization done in {} sec", timer.elapsed());

            return 0;
        }
        case Config::ASSEMBLE: {
            assert(files.size() == 1);
            assert(config->outfbase.size());

            Timer timer;
            logger->trace("Graph loading...");

            auto graph = load_critical_dbg(files.at(0));

            logger->trace("Graph loaded in {} sec", timer.elapsed());

            std::unique_ptr<AnnotatedDBG> anno_graph;
            if (config->infbase_annotators.size()) {
                anno_graph = initialize_annotated_dbg(graph, *config);

                logger->trace("Masking graph...");

                graph = mask_graph(*anno_graph, config.get());

                logger->trace("Masked in {} sec", timer.elapsed());
            }

            logger->trace("Extracting sequences from graph...");

            timer.reset();

            if (config->to_gfa) {
                if (!config->unitigs) {
                    logger->error("Flag '--unitigs' must be set for GFA output");
                    exit(1);
                }

                logger->trace("Writing graph to GFA...");

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

            logger->trace("Sequences extracted in {} sec", timer.elapsed());

            return 0;
        }
        case Config::RELAX_BRWT: {
            assert(files.size() == 1);
            assert(config->outfbase.size());

            Timer timer;

            auto annotator = std::make_unique<annotate::MultiBRWTAnnotator>();

            logger->trace("Loading annotator...");

            if (!annotator->load(files.at(0))) {
                logger->error("Cannot load annotations from file '{}'", files.at(0));
                exit(1);
            }
            logger->trace("Annotator loaded in {} sec", timer.elapsed());

            logger->trace("Relaxing BRWT tree...");

            annotate::relax_BRWT<annotate::MultiBRWTAnnotator>(annotator.get(),
                                                               config->relax_arity_brwt,
                                                               get_num_threads());

            annotator->serialize(config->outfbase);
            logger->trace("BRWT relaxation done in {} sec", timer.elapsed());

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
            ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1);
            std::mutex print_mutex;

            if (config->map_sequences) {
                if (!config->alignment_length) {
                    config->alignment_length = graph->get_k();
                } else if (config->alignment_length > graph->get_k()) {
                    logger->warn("Mapping to k-mers longer than k is not supported");
                    config->alignment_length = graph->get_k();
                }

                if (!dbg && config->alignment_length != graph->get_k()) {
                    logger->error("Matching k-mers shorter than k only "
                                  "supported for DBGSuccinct");
                    exit(1);
                }

                logger->trace("Map sequences against the de Bruijn graph with k={}",
                              graph->get_k());
                logger->trace("Length of mapped k-mers: {}", config->alignment_length);

                for (const auto &file : files) {
                    logger->trace("Map sequences from file '{}'", file);

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
                logger->info("Align sequences from file '{}'", file);

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

                    logger->trace("File '{}' processed in {} sec, "
                                  "current mem usage: {} MiB, total time {} sec",
                                  file, data_reading_timer.elapsed(),
                                  get_curr_RSS() >> 20, timer.elapsed());
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
                    logger->error("Failed to read accession->taxid map");
                    exit(1);
                }
            }

            auto anno_graph = initialize_annotated_dbg(*config);
            auto masked_graph = mask_graph(*anno_graph, config.get());

            logger->trace("Filter out: {}", fmt::join(config->label_filter, " "));

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

            ThreadPool thread_pool(std::max(1u, get_num_threads()) - 1);
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
                        logger->info("{}\t{}\t{}\t{}", alignment.front(), query,
                                     alignment.get_sequence(), label);
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
                logger->error("No variant calling mode selected. Exiting");
                exit(1);
            }

            thread_pool.join();

            logger->trace("# nodes checked: {}", masked_graph->num_nodes());
            logger->trace("# called: {}", num_calls);

            return 0;
        }
        case Config::PARSE_TAXONOMY: {
            TaxIDMapper taxid_mapper;
            if (config->accession2taxid.length()
                && !taxid_mapper.parse_accession2taxid(config->accession2taxid)) {
                logger->error("Failed to read accession->taxid file");
                exit(1);
            }

            if (config->taxonomy_nodes.length()
                && !taxid_mapper.parse_nodes(config->taxonomy_nodes)) {
                logger->error("Failed to read nodes.dmp file");
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
