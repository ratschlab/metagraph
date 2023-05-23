#include "annotate.hpp"

#include <filesystem>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/batch_accumulator.hpp"
#include "common/threads/threading.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "seq_io/formats.hpp"
#include "seq_io/sequence_io.hpp"
#include "seq_io/kmc_parser.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "graph/annotated_dbg.hpp"


namespace mtg {
namespace cli {

namespace fs = std::filesystem;

using namespace mtg::seq_io;

using mtg::common::logger;


template <class Callback>
void call_annotations(const std::string &file,
                      const std::string &ref_sequence_path,
                      const graph::DeBruijnGraph &graph,
                      bool forward_and_reverse,
                      size_t min_count,
                      size_t max_count,
                      bool filename_anno,
                      bool annotate_sequence_headers,
                      const std::string &fasta_anno_comment_delim,
                      const std::string &fasta_header_delimiter,
                      const std::vector<std::string> &anno_labels,
                      const Callback &callback) {
    size_t total_seqs = 0;

    Timer timer;

    logger->trace("Parsing '{}'", file);

    std::vector<std::string> labels = anno_labels;
    if (filename_anno) {
        labels.push_back(file);
    }
    // remember the number of base labels to remove those unique to each sequence quickly
    const size_t num_base_labels = labels.size();

    if (file_format(file) == "VCF") {
        read_vcf_file_with_annotations_critical(file, ref_sequence_path, graph.get_k(),
            [&](auto&& seq, const auto &variant_labels) {
                labels.insert(labels.end(),
                              variant_labels.begin(), variant_labels.end());

                callback(std::move(seq), labels);

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
                callback(std::string(sequence), labels);

                total_seqs += 1;

                if (logger->level() <= spdlog::level::level_enum::trace
                                                && total_seqs % 10000 == 0) {
                    logger->trace("processed {} sequences, annotated as <{}>, {} sec",
                                  total_seqs, fmt::join(labels, "><"), timer.elapsed());
                }
            },
            graph.get_mode() != graph::DeBruijnGraph::CANONICAL, // |call_both_from_canonical|
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

                callback(read_stream->seq.s, labels);

                total_seqs += 1;

                if (logger->level() <= spdlog::level::level_enum::trace
                                                && total_seqs % 10000 == 0) {
                    logger->trace("processed {} sequences, last was {}, annotated as <{}>, {} sec",
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

    logger->trace("File '{}' processed in {} sec, current mem usage: {} MiB",
                  file, timer.elapsed(), get_curr_RSS() >> 20);
}

template <class Callback>
void add_kmer_counts(const std::string &file,
                     const graph::DeBruijnGraph &graph,
                     bool forward_and_reverse,
                     bool filename_anno,
                     bool annotate_sequence_headers,
                     const std::string &fasta_anno_comment_delim,
                     const std::string &fasta_header_delimiter,
                     const std::vector<std::string> &anno_labels,
                     const Callback &callback) {
    size_t total_seqs = 0;

    Timer timer;

    std::vector<std::string> labels = anno_labels;
    if (filename_anno) {
        labels.push_back(file);
    }
    // remember the number of base labels to remove those unique to each sequence quickly
    const size_t num_base_labels = labels.size();

    const std::string &counts_fname
            = utils::remove_suffix(file, ".gz", ".fasta") + ".kmer_counts.gz";

    logger->trace("Parsing k-mer counts from '{}'", counts_fname);

    read_extended_fasta_file_critical<uint32_t>(file, "kmer_counts",
        [&](size_t k, const kseq_t *read_stream, const uint32_t *kmer_counts) {
            if (k != graph.get_k()) {
                logger->error("File '{}' contains counts for k-mers of "
                              "length {} but graph is constructed with k={}",
                              file, k, graph.get_k());
                exit(1);
            }
            assert(read_stream->seq.l >= k && "sequences can't be shorter than k-mers");

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

            callback(read_stream->seq.s, labels,
                     std::vector<uint64_t>(kmer_counts, kmer_counts + read_stream->seq.l - k + 1));

            total_seqs++;

            if (logger->level() <= spdlog::level::level_enum::trace
                                            && total_seqs % 10000 == 0) {
                logger->trace("processed {} sequences, last was {}, annotated as <{}>, {} sec",
                              total_seqs, read_stream->name.s, fmt::join(labels, "><"), timer.elapsed());
            }

            labels.resize(num_base_labels);
        },
        forward_and_reverse
    );
}

void annotate_data(std::shared_ptr<graph::DeBruijnGraph> graph,
                   const Config &config,
                   const std::vector<std::string> &files,
                   const std::string &annotator_filename,
                   size_t num_chunks_open = 2000) {
    auto anno_graph = initialize_annotated_dbg(graph, config, num_chunks_open);
    const size_t k = anno_graph->get_graph().get_k();

    bool forward_and_reverse = config.forward_and_reverse;
    if (anno_graph->get_graph().get_mode() == graph::DeBruijnGraph::CANONICAL) {
        logger->trace("Annotating canonical graph");
        forward_and_reverse = false;
    }

    Timer timer;

    ThreadPool thread_pool(get_num_threads() > 1 ? get_num_threads() : 0);

    // not too small, not too large
    const size_t batch_size = 1'000;
    const size_t batch_length = 100'000;

    if (config.coordinates) {
        for (const auto &file : files) {
            BatchAccumulator<std::tuple<std::string, std::vector<std::string>, uint64_t>> batcher(
                [&](auto&& data) {
                    thread_pool.enqueue([&](auto &data) {
                        anno_graph->annotate_kmer_coords(std::move(data));
                    }, std::move(data));
                },
                batch_size, batch_length, batch_size
            );

            logger->trace("Annotating k-mer coordinates for file {}", file);

            if (file_format(file) != "FASTA"
                    && file_format(file) != "FASTQ") {
                logger->error("Currently only FASTA or FASTQ format is supported"
                              " for annotating k-mer coordinates");
                exit(1);
            }

            uint64_t coord = 0;
            call_annotations(
                file,
                config.refpath,
                anno_graph->get_graph(),
                forward_and_reverse,
                config.min_count,
                config.max_count,
                config.filename_anno,
                config.annotate_sequence_headers,
                config.fasta_anno_comment_delim,
                config.fasta_header_delimiter,
                config.anno_labels,
                [&](std::string sequence, auto labels) {
                    if (config.num_kmers_in_seq
                            && config.num_kmers_in_seq + k - 1 != sequence.size()) {
                        logger->error("All input sequences must have the same"
                                      " length when flag --const-length is on");
                        exit(1);
                    }
                    if (sequence.size() >= k) {
                        uint64_t num_kmers = sequence.size() - k + 1;
                        batcher.push_and_pay(sequence.size(),
                                             std::move(sequence), std::move(labels), coord);
                        if (!config.annotate_sequence_headers)
                            coord += num_kmers;
                    }
                }
            );
        }

        thread_pool.join();

        anno_graph->get_annotator().serialize(annotator_filename);

        return;
    }

    // iterate over input files
    for (const auto &file : files) {
        BatchAccumulator<std::pair<std::string, std::vector<std::string>>> batcher(
            [&](auto&& data) {
                thread_pool.enqueue([&](auto &data) {
                    anno_graph->annotate_sequences(std::move(data));
                }, std::move(data));
            },
            batch_size, batch_length, batch_size
        );
        call_annotations(
            file,
            config.refpath,
            anno_graph->get_graph(),
            forward_and_reverse,
            config.min_count,
            config.max_count,
            config.filename_anno,
            config.annotate_sequence_headers,
            config.fasta_anno_comment_delim,
            config.fasta_header_delimiter,
            config.anno_labels,
            [&](std::string sequence, auto labels) {
                if (sequence.size() >= k) {
                    batcher.push_and_pay(sequence.size(),
                                         std::move(sequence), std::move(labels));
                }
            }
        );
    }

    thread_pool.join();

    if (config.count_kmers) {
        // add k-mer counts to existing binary annotations
        for (const auto &file : files) {
            logger->trace("Annotating k-mer counts for file {}", file);

            BatchAccumulator<std::tuple<std::string,
                                        std::vector<std::string>,
                                        std::vector<uint64_t>>> batcher(
                [&](auto&& data) {
                    using Batch = std::vector<std::tuple<std::string,
                                                         std::vector<std::string>,
                                                         std::vector<uint64_t>>>;
                    thread_pool.enqueue([&](Batch &data) {
                        for (auto &[seq, labels, kmer_counts] : data) {
                            anno_graph->add_kmer_counts(seq, labels, std::move(kmer_counts));
                        }
                    }, std::move(data));
                },
                batch_size, batch_length, batch_size
            );

            const std::string &counts_fname
                    = utils::remove_suffix(file, ".gz", ".fasta") + ".kmer_counts.gz";

            if (fs::exists(counts_fname)) {
                add_kmer_counts(
                    file,
                    anno_graph->get_graph(),
                    forward_and_reverse,
                    config.filename_anno,
                    config.annotate_sequence_headers,
                    config.fasta_anno_comment_delim,
                    config.fasta_header_delimiter,
                    config.anno_labels,
                    [&](std::string sequence,
                                std::vector<std::string> labels,
                                std::vector<uint64_t> kmer_counts) {
                        if (sequence.size() >= k) {
                            batcher.push_and_pay(sequence.size(),
                                                 std::move(sequence), std::move(labels),
                                                 std::move(kmer_counts));
                        }
                    }
                );
            } else {
                logger->warn("No k-mer counts found at '{}'. Every input k-mer"
                             " will have count 1.", counts_fname);
                call_annotations(
                    file,
                    config.refpath,
                    anno_graph->get_graph(),
                    forward_and_reverse,
                    config.min_count,
                    config.max_count,
                    config.filename_anno,
                    config.annotate_sequence_headers,
                    config.fasta_anno_comment_delim,
                    config.fasta_header_delimiter,
                    config.anno_labels,
                    [&](std::string sequence, auto labels) {
                        if (sequence.size() >= k) {
                            batcher.push_and_pay(sequence.size(),
                                                 std::move(sequence), std::move(labels),
                                                 std::vector<uint64_t>(sequence.size() - k + 1, 1));
                        }
                    }
                );
            }
        }
    }

    thread_pool.join();

    anno_graph->get_annotator().serialize(annotator_filename);
}


int annotate_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->infbase_annotators.size() <= 1);

    const auto graph = load_critical_dbg(config->infbase);

    if (!config->separately) {
        annotate_data(graph, *config, files, config->outfbase);

    } else {
        // |config->separately| is true

        // annotate multiple files in parallel, each with |parallel_each| threads
        size_t num_threads = get_num_threads();
        set_num_threads(std::max(1u, config->parallel_each));

        if (!config->outfbase.empty()) {
            try {
                fs::create_directory(config->outfbase);
            } catch (...) {
                logger->error("Failed to create directory {}", config->outfbase);
                throw;
            }
        }

        #pragma omp parallel for num_threads(num_threads) default(shared) schedule(dynamic, 1)
        for (size_t i = 0; i < files.size(); ++i) {
            annotate_data(graph, *config, { files[i] },
                config->outfbase.size()
                    ? config->outfbase + "/" + utils::split_string(files[i], "/").back()
                    : files[i],
                2000 / std::max((size_t)1, num_threads));
        }
    }

    return 0;
}

} // namespace cli
} // namespace mtg
