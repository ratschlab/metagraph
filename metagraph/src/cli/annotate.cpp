#include "annotate.hpp"

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "seq_io/formats.hpp"
#include "seq_io/sequence_io.hpp"
#include "seq_io/kmc_parser.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"


namespace mtg {
namespace cli {

using namespace mtg::seq_io;

using mtg::common::logger;


void annotate_data(const std::vector<std::string> &files,
                   const std::string &ref_sequence_path,
                   graph::AnnotatedDBG *anno_graph,
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
            read_kmers(
                file,
                [&](std::string_view sequence) {
                    thread_pool.enqueue(
                        [anno_graph](const std::string &sequence, const auto &labels) {
                            anno_graph->annotate_sequence(sequence, labels);
                        },
                        std::string(sequence),
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
                          graph::AnnotatedDBG *anno_graph,
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


int annotate_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase_annotators.size() <= 1);

    const auto graph = load_critical_dbg(config->infbase);

    if (graph->is_canonical_mode() || config->canonical)
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

int annotate_graph_with_genome_coordinates(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase_annotators.size() <= 1);

    auto graph_temp = load_critical_dbg(config->infbase);

    auto annotation_temp
        = std::make_unique<annot::RowCompressed<>>(graph_temp->max_index());

    if (config->infbase_annotators.size()
            && !annotation_temp->load(config->infbase_annotators.at(0))) {
        logger->error("Cannot load annotations from '{}'",
                      config->infbase_annotators.at(0));
        exit(1);
    }

    // load graph
    graph::AnnotatedDBG anno_graph(graph_temp,
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

} // namespace cli
} // namespace mtg
