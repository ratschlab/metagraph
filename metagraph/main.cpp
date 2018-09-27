#include "annotated_dbg.hpp"
#include "dbg_succinct_chunk.hpp"
#include "dbg_succinct_construct.hpp"
#include "reads_filtering.hpp"
#include "config.hpp"
#include "utils.hpp"
#include "sequence_io.hpp"
#include "dbg_succinct_merge.hpp"
#include "annotate_column_compressed.hpp"
#include "annotate_column_compressed_fast.hpp"
#include "annotate_row_compressed.hpp"
#include "annotate_bloom_filter.hpp"
#include "unix_tools.hpp"
#include "serialization.hpp"

typedef annotate::MultiLabelAnnotation<uint64_t, std::string> Annotator;

const size_t kMaxNumParallelReadFiles = 5;

const size_t kNumCachedColumns = 10;


DBG_succ* load_critical_graph_from_file(const std::string &filename) {
    std::string filetype = ".dbg";
    std::string filebase =
        std::equal(filetype.rbegin(), filetype.rend(), filename.rbegin())
        ? filename.substr(0, filename.size() - filetype.size())
        : filename;

    DBG_succ *graph = new DBG_succ();
    if (!graph->load(filebase)) {
        std::cerr << "ERROR: input file "
                  << filename << " corrupted" << std::endl;
        delete graph;
        exit(1);
    }
    return graph;
}

std::string get_filter_filename(std::string filename,
                                size_t k,
                                size_t max_unreliable_abundance,
                                size_t unreliable_kmers_threshold,
                                bool critical = true) {
    filename = filename + ".filter_k" + std::to_string(k)
                        + "_s" + std::to_string(max_unreliable_abundance)
                        + (max_unreliable_abundance
                             ? std::string("_")
                                + std::to_string(unreliable_kmers_threshold)
                             : "");
    if (!critical || std::ifstream(filename).good()) {
        return filename;
    } else if (max_unreliable_abundance == 0) {
        return "";
    }

    std::cerr << "ERROR: read filter "
              << filename << " does not exist."
              << " Filter reads first." << std::endl;
    exit(1);
}


void annotate_data(const std::vector<std::string> &files,
                   const std::string &ref_sequence_path,
                   AnnotatedDBG *anno_graph,
                   bool reverse,
                   size_t filter_k,
                   size_t max_unreliable_abundance,
                   size_t unreliable_kmers_threshold,
                   bool filename_anno,
                   bool fasta_anno,
                   const std::string &fasta_header_delimiter,
                   const std::vector<std::string> &anno_labels,
                   bool verbose) {
    size_t total_seqs = 0;

    std::unique_ptr<Timer> timer;
    if (verbose) {
        timer.reset(new Timer());
    }

    // iterate over input files
    for (const auto &file : files) {
        if (verbose) {
            std::cout << std::endl << "Parsing " << file << std::endl;
        }
        // open stream
        if (utils::get_filetype(file) == "VCF") {
            std::vector<std::string> variant_labels;

            read_vcf_file_critical(file,
                                   ref_sequence_path,
                                   anno_graph->get_graph().get_k(),
                                   &variant_labels,
                [&](std::string &seq, std::vector<std::string> *variant_labels) {
                    assert(variant_labels);

                    if (filename_anno)
                        variant_labels->push_back(file);

                    for (const auto &label : anno_labels) {
                        variant_labels->push_back(label);
                    }

                    anno_graph->annotate_sequence(seq, *variant_labels);
                    if (reverse) {
                        reverse_complement(seq.begin(), seq.end());
                        anno_graph->annotate_sequence(seq, *variant_labels);
                    }
                    variant_labels->clear();
                }
            );
        } else if (utils::get_filetype(file) == "FASTA"
                    || utils::get_filetype(file) == "FASTQ") {
            read_fasta_file_critical(file,
                [&](kseq_t *read_stream) {
                    std::vector<std::string> labels;

                    if (fasta_anno) {
                        labels = utils::split_string(read_stream->name.s,
                                                     fasta_header_delimiter);
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
                        std::cout << ", " << timer->elapsed() << "sec" << std::endl;
                    }
                },
                reverse, timer.get(),
                get_filter_filename(
                    file, filter_k,
                    max_unreliable_abundance,
                    unreliable_kmers_threshold
                )
            );
        } else {
            std::cerr << "ERROR: Filetype unknown for file "
                      << file << std::endl;
            exit(1);
        }
    }

    // join threads if any were initialized
    anno_graph->join();
}


void annotate_coordinates(const std::vector<std::string> &files,
                          AnnotatedDBG *anno_graph,
                          bool reverse,
                          size_t filter_k,
                          size_t max_unreliable_abundance,
                          size_t unreliable_kmers_threshold,
                          size_t genome_bin_size,
                          bool verbose) {
    size_t total_seqs = 0;

    std::unique_ptr<Timer> timer;
    if (verbose)
        timer.reset(new Timer());

    // iterate over input files
    for (const auto &file : files) {
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
                        std::to_string(reverse && (total_seqs % 2)), // whether the read is reverse
                        "",
                    };

                    const std::string sequence(read_stream->seq.s);
                    for (size_t i = 0; i < sequence.size(); i += genome_bin_size) {
                        labels.back() = std::to_string(i);

                        // forward: |0 =>  |6 =>  |12=>  |18=>  |24=>  |30=>|
                        // reverse: |<=30|  <=24|  <=18|  <=12|  <= 6|  <= 0|
                        const size_t bin_size = std::min(
                            static_cast<size_t>(sequence.size() - i),
                            static_cast<size_t>(genome_bin_size + anno_graph->get_graph().get_k())
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
                        std::cout << ", " << timer->elapsed() << "sec" << std::endl;
                    }

                    // If we read both strands, the next sequence is
                    // either reverse (if the current one is forward)
                    // or new (if the current one is reverse), and therefore forward
                    if (reverse)
                        forward_strand = !forward_strand;
                },
                reverse, timer.get(),
                get_filter_filename(
                    file, filter_k,
                    max_unreliable_abundance,
                    unreliable_kmers_threshold
                )
            );
        } else {
            std::cerr << "ERROR: the type of file "
                      << file << " is not supported" << std::endl;
            exit(1);
        }
    }

    // join threads if any were initialized
    anno_graph->join();
}


/**
 * ACGT, ACG$, ..., $$$$ -- valid
 * AC$T, A$$T, ..., $AAA -- invalid
 */
bool valid_kmer_suffix(const std::string &suffix) {
    size_t last = suffix.rfind(DBG_succ::kSentinel);
    return last == std::string::npos
            || suffix.substr(0, last + 1)
                == std::string(last + 1, DBG_succ::kSentinel);
}


void execute_query(std::string seq_name,
                   std::string sequence,
                   bool count_labels,
                   bool suppress_unlabeled,
                   size_t num_top_labels,
                   double discovery_fraction,
                   std::string anno_labels_delimiter,
                   const AnnotatedDBG &anno_graph) {
    std::ostringstream oss;

    if (count_labels) {
        auto top_labels = anno_graph.get_top_labels(sequence, num_top_labels);

        if (!top_labels.size() && suppress_unlabeled)
            return;

        oss << seq_name << "\t";

        if (top_labels.size()) {
            oss << "<" << top_labels[0].first << ">:" << top_labels[0].second;
        }
        for (size_t i = 1; i < top_labels.size(); ++i) {
            oss << "\t<" << top_labels[i].first << ">:" << top_labels[i].second;
        }
        oss << "\n";
    } else {
        auto labels_discovered
                = anno_graph.get_labels(sequence, discovery_fraction);

        if (!labels_discovered.size() && suppress_unlabeled)
            return;

        oss << seq_name << "\t"
            << utils::join_strings(labels_discovered,
                                   anno_labels_delimiter) << "\n";
    }

    std::cout << oss.str();
}


int main(int argc, const char *argv[]) {
    // parse command line arguments and options
    std::unique_ptr<Config> config { new Config(argc, argv) };

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
            DBG_succ *graph = new DBG_succ(config->k);

            if (config->verbose)
                std::cout << "Build Succinct De Bruijn Graph with k-mer size k="
                          << graph->get_k() << std::endl;

            Timer timer;

            if (!config->dynamic) {
                if (config->verbose) {
                    std::cout << "Start reading data and extracting k-mers" << std::endl;
                }
                //enumerate all suffices
                assert(DBG_succ::alph_size > 1);
                assert(config->nsplits > 0);
                size_t suffix_len = std::min(
                    static_cast<size_t>(std::ceil(std::log2(config->nsplits)
                                                    / std::log2(DBG_succ::alph_size - 1))),
                    graph->get_k() - 1
                );
                std::deque<std::string> suffices;
                if (config->suffix.size()) {
                    suffices = { config->suffix };
                } else {
                    suffices = utils::generate_strings(
                        DBG_succ::alphabet.substr(0, DBG_succ::alph_size),
                        suffix_len
                    );
                }

                DBG_succ::Chunk graph_data(graph->get_k());

                //one pass per suffix
                for (const std::string &suffix : suffices) {
                    timer.reset();

                    if (suffix.size() > 0 || suffices.size() > 1) {
                        if (valid_kmer_suffix(suffix)) {
                            std::cout << "\nSuffix: " << suffix << std::endl;
                        } else {
                            std::cout << "\nSkipping suffix: " << suffix << std::endl;
                            continue;
                        }
                    }

                    std::unique_ptr<IChunkConstructor> constructor(
                        IChunkConstructor::initialize(
                            graph->get_k(),
                            suffix,
                            config->parallel,
                            static_cast<uint64_t>(config->memory_available) << 30,
                            config->verbose
                        )
                    );

                    // iterate over input files
                    for (unsigned int f = 0; f < files.size(); ++f) {
                        if (config->verbose) {
                            std::cout << std::endl << "Parsing " << files[f] << std::endl;
                        }

                        Timer data_reading_timer;

                        if (utils::get_filetype(files[f]) == "VCF") {
                            if (suffix.find('$') != std::string::npos)
                                continue;

                            Timer *timer_ptr = config->verbose ? &timer : NULL;
                            //assume VCF contains no noise
                            read_vcf_file_critical(
                                files[f], config->refpath, graph->get_k(), NULL,
                                [&](std::string &seq, auto *) {
                                    constructor->add_read(seq);
                                    if (config->reverse) {
                                        reverse_complement(seq.begin(), seq.end());
                                        constructor->add_read(seq);
                                    }
                                }, timer_ptr
                            );
                        } else if (utils::get_filetype(files[f]) == "FASTA"
                                    || utils::get_filetype(files[f]) == "FASTQ") {
                            Timer *timer_ptr = config->verbose ? &timer : NULL;

                            std::string filter_filename = get_filter_filename(
                                files[f], config->filter_k,
                                config->max_unreliable_abundance,
                                config->unreliable_kmers_threshold
                            );

                            if (files.size() >= config->parallel) {
                                auto reverse = config->reverse;
                                auto file = files[f];

                                // capture all required values by copying to be able
                                // to run task from other threads
                                constructor->add_reads([=](auto callback) {
                                    read_fasta_file_critical(file, [=](kseq_t *read_stream) {
                                        // add read to the graph constructor as a callback
                                        callback(read_stream->seq.s);
                                    }, reverse, timer_ptr, filter_filename);
                                });
                            } else {
                                read_fasta_file_critical(files[f], [&](kseq_t *read_stream) {
                                    // add read to the graph constructor as a callback
                                    constructor->add_read(read_stream->seq.s);
                                }, config->reverse, timer_ptr, filter_filename);
                            }
                        } else {
                            std::cerr << "ERROR: Filetype unknown for file "
                                      << files[f] << std::endl;
                            exit(1);
                        }
                        if (config->verbose) {
                            std::cout << "File processed in "
                                      << data_reading_timer.elapsed()
                                      << "sec, current mem usage: "
                                      << get_curr_mem2() / (1<<20) << " MB"
                                      << ", total time: " << timer.elapsed()
                                      << "sec" << std::endl;
                        }
                    }
                    std::cout << std::endl;

                    auto next_block = constructor->build_chunk();
                    std::cout << "Graph chunk with " << next_block->size()
                              << " k-mers was built in "
                              << timer.elapsed() << "sec" << std::endl;

                    if (config->outfbase.size() && config->suffix.size()) {
                        std::cout << "Serialize the graph chunk for suffix '"
                                  << suffix << "'...\t" << std::flush;
                        timer.reset();
                        next_block->serialize(config->outfbase + "." + suffix);
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    }

                    graph_data.extend(*next_block);
                    delete next_block;
                }

                if (!config->suffix.size())
                    graph_data.initialize_graph(graph);

            } else {
                //slower method
                for (unsigned int f = 0; f < files.size(); ++f) {
                    if (config->verbose) {
                        std::cout << std::endl << "Parsing " << files[f] << std::endl;
                    }
                    // open stream
                    if (utils::get_filetype(files[f]) == "VCF") {
                        read_vcf_file_critical(
                            files[f], config->refpath, graph->get_k(), NULL,
                            [&](std::string &seq, auto *) {
                                graph->add_sequence(seq);
                                if (config->reverse) {
                                    reverse_complement(seq.begin(), seq.end());
                                    graph->add_sequence(seq);
                                }
                            }
                        );
                    } else if (utils::get_filetype(files[f]) == "FASTA"
                                || utils::get_filetype(files[f]) == "FASTQ") {
                        read_fasta_file_critical(files[f],
                            [&](kseq_t *read_stream) {
                                graph->add_sequence(read_stream->seq.s);
                            },
                            config->reverse, NULL,
                            get_filter_filename(files[f], config->filter_k,
                                                config->max_unreliable_abundance,
                                                config->unreliable_kmers_threshold)
                        );
                    } else {
                        std::cerr << "ERROR: Filetype unknown for file "
                                  << files[f] << std::endl;
                        exit(1);
                    }
                }
            }

            graph->switch_state(config->state);

            // graph output
            if (config->print_graph_succ && !config->suffix.size())
                std::cout << *graph;
            if (!config->outfbase.empty() && !config->suffix.size())
                graph->serialize(config->outfbase);
            delete graph;

            return 0;
        }
        case Config::EXTEND: {
            Timer timer;

            // load graph
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(config->infbase)
            };

            if (config->verbose) {
                std::cout << "Succinct de Bruijn graph with k-mer size k="
                          << graph->get_k() << " has been loaded in "
                          << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();

            std::unique_ptr<bit_vector_dyn> inserted_edges;
            if (config->infbase_annotators.size())
                inserted_edges.reset(new bit_vector_dyn(graph->num_edges() + 1, 0));

            if (graph->get_state() != Config::DYN) {
                if (config->verbose)
                    std::cout << "Switching the graph state to dynamic..." << std::flush;
                graph->switch_state(Config::DYN);
                if (config->verbose)
                    std::cout << "\tdone in " << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();

            if (config->verbose)
                std::cout << "Start graph extension" << std::endl;

            for (const auto &file : files) {
                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << file << std::endl;
                }
                // open stream
                if (utils::get_filetype(file) == "VCF") {
                    read_vcf_file_critical(file, config->refpath, graph->get_k(), NULL,
                        [&](std::string &seq, auto *) {
                            graph->add_sequence(seq, true, inserted_edges.get());
                            if (config->reverse) {
                                reverse_complement(seq.begin(), seq.end());
                                graph->add_sequence(seq, true, inserted_edges.get());
                            }
                        }
                    );
                } else if (utils::get_filetype(file) == "FASTA"
                            || utils::get_filetype(file) == "FASTQ") {
                    read_fasta_file_critical(file,
                        [&](kseq_t *read_stream) {
                            graph->add_sequence(read_stream->seq.s, true, inserted_edges.get());
                        },
                        config->reverse, NULL,
                        get_filter_filename(file, config->filter_k,
                                            config->max_unreliable_abundance,
                                            config->unreliable_kmers_threshold)
                    );
                } else {
                    std::cerr << "ERROR: Filetype unknown for file "
                              << file << std::endl;
                    exit(1);
                }
            }

            if (config->verbose)
                std::cout << "Graph extension finished in "
                          << timer.elapsed() << "sec" << std::endl;

            // graph output
            if (config->print_graph_succ)
                std::cout << *graph;

            assert(config->outfbase.size());

            // serialize graph
            timer.reset();
            if (graph->get_state() != config->state) {
                if (config->verbose)
                    std::cout << "Switching state before dumping..." << std::flush;

                graph->switch_state(config->state);

                if (config->verbose)
                    std::cout << "\tdone in " << timer.elapsed() << "sec" << std::endl;
            }

            graph->serialize(config->outfbase);
            graph.reset();

            timer.reset();

            if (!config->infbase_annotators.size())
                return 0;

            std::unique_ptr<Annotator> annotation;

            if (config->use_row_annotator) {
                annotation.reset(
                    new annotate::RowCompressed<>(1, config->sparse)
                );
            } else {
                annotation.reset(
                    new annotate::ColumnCompressed<>(1, kNumCachedColumns,
                                                     config->verbose)
                );
            }

            if (!annotation->merge_load(config->infbase_annotators)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            } else if (config->verbose) {
                std::cout << "Annotation was loaded in "
                          << timer.elapsed() << "sec" << std::endl;
            }
            timer.reset();

            AnnotatedDBG anno_dbg(annotation.release());

            assert(inserted_edges.get());
            if (config->verbose)
                std::cout << "Insert empty rows to the annotation matrix..." << std::flush;

            anno_dbg.adjust_annotation(*inserted_edges);

            if (config->verbose)
                std::cout << "\tdone in " << timer.elapsed() << "sec" << std::endl;

            anno_dbg.get_annotation().serialize(config->outfbase);

            return 0;
        }
        case Config::FILTER: {
            if (config->verbose) {
                std::cout << "Filter out reads with rare k-mers" << std::endl;
                std::cout << "Start reading data and extracting k-mers" << std::endl;
            }

            utils::ThreadPool thread_pool_files(kMaxNumParallelReadFiles);
            utils::ThreadPool thread_pool(std::max(1u, config->parallel) - 1,
                                          std::max(1u, config->parallel));
            Timer timer;

            // iterate over input files
            for (const auto &file : files) {
                if (utils::get_filetype(file) != "FASTA"
                        && utils::get_filetype(file) != "FASTQ") {
                    std::cerr << "ERROR: Filetype unknown for file "
                              << file << std::endl;
                    exit(1);
                }

                Timer *timer_ptr = config->verbose ? &timer : NULL;

                if (config->generate_filtered_fasta || config->generate_filtered_fastq) {
                    thread_pool.enqueue([=](size_t k,
                                            size_t max_unreliable_abundance,
                                            size_t unreliable_kmers_threshold,
                                            bool out_fasta,
                                            bool out_fastq) {
                        assert(out_fasta || out_fastq);

                        auto filter_filename = get_filter_filename(
                            file, k, max_unreliable_abundance, unreliable_kmers_threshold
                        );

                        std::string filtered_reads_fasta = get_filter_filename(
                            utils::remove_suffix(file, ".gz", ".fasta", ".fastq"),
                            k, max_unreliable_abundance, unreliable_kmers_threshold, false
                        ) + ".fasta.gz";
                        std::string filtered_reads_fastq = get_filter_filename(
                            utils::remove_suffix(file, ".gz", ".fasta", ".fastq"),
                            k, max_unreliable_abundance, unreliable_kmers_threshold, false
                        ) + ".fastq.gz";

                        gzFile out_fasta_gz = Z_NULL;
                        gzFile out_fastq_gz = Z_NULL;

                        if (out_fasta) {
                            out_fasta_gz = gzopen(filtered_reads_fasta.c_str(), "w");
                            if (out_fasta_gz == Z_NULL) {
                                std::cerr << "ERROR: Can't write to "
                                          << filtered_reads_fasta << std::endl;
                                exit(1);
                            }
                        }
                        if (out_fastq) {
                            out_fastq_gz = gzopen(filtered_reads_fastq.c_str(), "w");
                            if (out_fastq_gz == Z_NULL) {
                                std::cerr << "ERROR: Can't write to "
                                          << filtered_reads_fastq << std::endl;
                                exit(1);
                            }
                        }

                        read_fasta_file_critical(file,
                            [&](kseq_t *read_stream) {
                                if (out_fasta_gz != Z_NULL
                                        && !write_fasta(out_fasta_gz, *read_stream)) {
                                    std::cerr << "ERROR: Can't write filtered reads to "
                                              << filtered_reads_fasta << std::endl;
                                    exit(1);
                                }
                                if (out_fastq_gz != Z_NULL
                                        && !write_fastq(out_fastq_gz, *read_stream)) {
                                    std::cerr << "ERROR: Can't write filtered reads to "
                                              << filtered_reads_fastq << std::endl;
                                    exit(1);
                                }
                            },
                            false, timer_ptr, filter_filename
                        );

                        if (out_fasta_gz != Z_NULL)
                            gzclose(out_fasta_gz);
                        if (out_fastq_gz != Z_NULL)
                            gzclose(out_fastq_gz);

                    },
                    config->filter_k,
                    config->max_unreliable_abundance,
                    config->unreliable_kmers_threshold,
                    config->generate_filtered_fasta,
                    config->generate_filtered_fastq);

                    continue;
                }

                auto *thread_pool_ptr = &thread_pool;
                // capture all required values by copying to be able
                // to run task from other threads
                thread_pool_files.enqueue([=](size_t k,
                                              size_t max_unreliable_abundance,
                                              size_t unreliable_kmers_threshold,
                                              bool verbose,
                                              bool reverse,
                                              bool use_kmc) {
                        // compute read filter, bit vector indicating filtered reads
                        std::vector<bool> filter = filter_reads([=](auto callback) {
                                read_fasta_file_critical(file, [=](kseq_t *read_stream) {
                                    // add read to the graph constructor as a callback
                                    callback(read_stream->seq.s);
                                }, reverse, timer_ptr);
                            },
                            k, max_unreliable_abundance, unreliable_kmers_threshold,
                            verbose, thread_pool_ptr,
                            (use_kmc && max_unreliable_abundance) ? file + ".kmc" : ""
                        );
                        if (reverse) {
                            assert(filter.size() % 2 == 0);

                            std::vector<bool> forward_filter(filter.size() / 2);
                            for (size_t i = 0; i < forward_filter.size(); ++i) {
                                forward_filter[i] = filter[i * 2] | filter[i * 2 + 1];
                            }
                            filter.swap(forward_filter);
                        }

                        // dump filter
                        std::ofstream outstream(get_filter_filename(
                            file, k, max_unreliable_abundance, unreliable_kmers_threshold, false
                        ));
                        serialize_number_vector(outstream, filter, 1);
                    },
                    config->filter_k,
                    config->max_unreliable_abundance,
                    config->unreliable_kmers_threshold,
                    config->verbose, config->reverse, config->use_kmc
                );
            }
            thread_pool_files.join();
            thread_pool.join();

            return 0;
        }
        case Config::FILTER_STATS: {
            if (!config->verbose) {
                std::cout << "File\tTotalReads\tRemainingReads\tRemainingRatio" << std::endl;
            }

            for (const auto &file : files) {
                auto filter_filename = get_filter_filename(
                    file,
                    config->filter_k,
                    config->max_unreliable_abundance,
                    config->unreliable_kmers_threshold,
                    true
                );

                std::ifstream instream(filter_filename);
                try {
                    auto filter = load_number_vector<bool>(instream);
                    size_t num_remaining = std::count(filter.begin(), filter.end(), true);

                    if (config->verbose) {
                        std::cout << "Statistics for filter file " << filter_filename << "\n"
                                  << "Total reads: " << filter.size() << "\n"
                                  << "Remaining reads: " << num_remaining << "\n"
                                  << "Remaining ratio: "
                                  << static_cast<double>(num_remaining) / filter.size()
                                  << std::endl;
                    } else {
                        std::cout << filter_filename << "\t"
                                  << filter.size() << "\t"
                                  << num_remaining << "\t"
                                  << static_cast<double>(num_remaining) / filter.size()
                                  << std::endl;
                    }
                } catch (...) {
                    std::cerr << "ERROR: Filter file " << filter_filename
                              << " is corrupted" << std::endl;
                    exit(1);
                }
            }

            return 0;
        }
        case Config::ANNOTATE: {
            // load graph
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(config->infbase)
            };

            // initialize empty annotation
            std::unique_ptr<Annotator> annotation;
            if (config->use_row_annotator) {
                annotation.reset(
                    new annotate::RowCompressed<>(graph->num_edges() + 1,
                                                  config->sparse)
                );
            } else {
                annotation.reset(
                    new annotate::ColumnCompressed<>(graph->num_edges() + 1,
                                                     kNumCachedColumns,
                                                     config->verbose)
                );
            }

            if (config->infbase_annotators.size()
                    && !annotation->merge_load(config->infbase_annotators)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }

            AnnotatedDBG annotated_dbg(graph.release(),
                                       annotation.release(),
                                       config->parallel);

            annotate_data(files,
                          config->refpath,
                          &annotated_dbg,
                          config->reverse,
                          config->filter_k,
                          config->max_unreliable_abundance,
                          config->unreliable_kmers_threshold,
                          config->filename_anno,
                          config->fasta_anno,
                          config->fasta_header_delimiter,
                          config->anno_labels,
                          config->verbose);

            annotated_dbg.get_annotation().serialize(config->outfbase);

            return 0;
        }
        case Config::ANNOTATE_COORDINATES: {
            // load graph
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(config->infbase)
            };

            // initialize empty annotation
            std::unique_ptr<Annotator> annotation;
            annotation.reset(
                new annotate::RowCompressed<>(graph->num_edges() + 1)
            );

            if (config->infbase_annotators.size()
                    && !annotation->merge_load(config->infbase_annotators)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }

            AnnotatedDBG annotated_dbg(graph.release(),
                                       annotation.release(),
                                       config->parallel);

            annotate_coordinates(files,
                                 &annotated_dbg,
                                 config->reverse,
                                 config->filter_k,
                                 config->max_unreliable_abundance,
                                 config->unreliable_kmers_threshold,
                                 config->genome_binsize_anno,
                                 config->verbose);

            annotated_dbg.get_annotation().serialize(config->outfbase);

            return 0;
        }
        case Config::ANNOTATE_BLOOM: {
            // load graph
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(config->infbase)
            };

            // initialize empty Bloom filter annotation
            std::unique_ptr<annotate::AnnotationCategoryBloom> annotation;

            if (config->bloom_fpp > -0.5) {
                // Expected FPP is set, optimize other parameters automatically
                annotation.reset(
                    new annotate::AnnotationCategoryBloom(
                        *graph,
                        config->bloom_fpp,
                        config->verbose
                    )
                );
            } else {
                assert(config->bloom_bits_per_edge >= 0);

                // Experiment mode, estimate FPP given other parameters,
                // optimize the number of hash functions if it's set to zero
                annotation.reset(
                    new annotate::AnnotationCategoryBloom(
                        *graph,
                        config->bloom_bits_per_edge,
                        config->bloom_num_hash_functions,
                        config->verbose
                    )
                );
            }

            AnnotatedDBG annotated_dbg(graph.release(),
                                       annotation.release(),
                                       config->parallel);

            annotate_data(files,
                          config->refpath,
                          &annotated_dbg,
                          config->reverse,
                          config->filter_k,
                          config->max_unreliable_abundance,
                          config->unreliable_kmers_threshold,
                          config->filename_anno,
                          config->fasta_anno,
                          config->fasta_header_delimiter,
                          config->anno_labels,
                          config->verbose);

            annotated_dbg.get_annotation().serialize(config->infbase);

            return 0;
        }
        case Config::MERGE_ANNOTATIONS: {
            std::unique_ptr<Annotator> annotation;
            if (config->use_row_annotator) {
                throw std::runtime_error("To be implemented");
                annotation.reset(new annotate::RowCompressed<>(0, config->sparse));
            } else {
                annotation.reset(
                    new annotate::ColumnCompressed<>(
                        0, kNumCachedColumns, config->verbose
                    )
                );
            }

            if (!annotation->merge_load(files)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }

            annotation->serialize(config->outfbase);

            return 0;
        }
        case Config::CLASSIFY: {
            // load graph
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(config->infbase)
            };

            std::unique_ptr<Annotator> annotation;
            if (config->use_row_annotator) {
                annotation.reset(new annotate::RowCompressed<>(0, config->sparse));
            } else if (config->fast) {
                annotation.reset(
                    new annotate::FastColumnCompressed<>(
                        0, kNumCachedColumns, config->verbose
                    )
                );
            } else {
                annotation.reset(
                    new annotate::ColumnCompressed<>(
                        0, kNumCachedColumns, config->verbose
                    )
                );
            }

            if (!annotation->merge_load(config->infbase_annotators)) {
                std::cerr << "ERROR: can't load annotations for graph "
                          << config->infbase + ".dbg"
                          << ", file corrupted" << std::endl;
                exit(1);
            }

            utils::ThreadPool thread_pool(std::max(1u, config->parallel) - 1);

            std::unique_ptr<Timer> timer;
            if (config->verbose) {
                timer.reset(new Timer());
            }

            AnnotatedDBG annotated_dbg(graph.release(),
                                       annotation.release(),
                                       config->parallel);

            // iterate over input files
            for (const auto &file : files) {
                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << file << std::endl;
                }

                size_t seq_count = 0;

                read_fasta_file_critical(file,
                    [&](kseq_t *read_stream) {
                        thread_pool.enqueue(execute_query,
                            std::to_string(seq_count++) + "\t"
                                + std::string(read_stream->name.s),
                            std::string(read_stream->seq.s),
                            config->count_labels,
                            config->suppress_unlabeled,
                            config->num_top_labels,
                            config->discovery_fraction,
                            config->anno_labels_delimiter,
                            std::ref(annotated_dbg)
                        );
                    },
                    config->reverse, timer.get(),
                    get_filter_filename(file, config->filter_k,
                                        config->max_unreliable_abundance,
                                        config->unreliable_kmers_threshold)
                );

                // wait while all threads finish processing the current file
                thread_pool.join();
            }

            return 0;
        }
        case Config::COMPARE: {
            assert(files.size());

            std::cout << "Opening file                " << files.at(0) << std::endl;
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(files.at(0))
            };

            for (size_t f = 1; f < files.size(); ++f) {
                std::cout << "Opening file for comparison " << files[f] << std::endl;
                DBG_succ *second = load_critical_graph_from_file(files[f]);
                if (config->internal
                        ? graph->equals_internally(*second, config->verbose)
                        : *graph == *second) {
                    std::cout << "Graphs are identical" << std::endl;
                } else {
                    std::cout << "Graphs are not identical" << std::endl;
                }
                delete second;
            }

            return 0;
        }
        case Config::CONCATENATE: {
            auto chunk_files = files;

            if (!files.size()) {
                assert(config->infbase.size());

                auto sorted_suffices = utils::generate_strings(
                    DBG_succ::alphabet.substr(0, DBG_succ::alph_size),
                    config->suffix_len
                );

                for (const std::string &suffix : sorted_suffices) {
                    assert(suffix.size() == config->suffix_len);

                    if (valid_kmer_suffix(suffix))
                        chunk_files.push_back(config->infbase + "." + suffix);
                }
            }

            for (auto &filename : chunk_files) {
                filename = utils::remove_suffix(filename, ".dbgchunk");
            }

            // collect results on an external merge or construction
            DBG_succ *graph = DBG_succ::Chunk::build_graph_from_chunks(
                chunk_files, config->verbose
            );
            assert(graph);

            if (config->verbose) {
                std::cout << "Succinct graph has been assembled" << std::endl;
                std::cout << "nodes: " << graph->num_nodes() << std::endl;
                std::cout << "edges: " << graph->num_edges() << std::endl;
                std::cout << "k: " << graph->get_k() << std::endl;
                std::cout << "state: " << Config::state_to_string(graph->get_state())
                          << std::endl;
            }

            // graph output
            if (graph && config->print_graph_succ)
                std::cout << *graph;
            if (graph && !config->outfbase.empty())
                graph->serialize(config->outfbase);

            return 0;
        }
        case Config::MERGE: {
            DBG_succ *graph = NULL;

            Timer timer;

            std::vector<const DBG_succ*> graphs;
            for (const auto &file : files) {
                std::cout << "Opening file " << file << std::endl;
                graphs.push_back(load_critical_graph_from_file(file));
                if (config->verbose) {
                    std::cout << "nodes: " << graphs.back()->num_nodes() << "\n";
                    std::cout << "edges: " << graphs.back()->num_edges() << "\n";
                    std::cout << "k: " << graphs.back()->get_k() << std::endl;
                }
            }
            std::cout << "Graphs are loaded\t" << timer.elapsed()
                                               << "sec" << std::endl;

            if (config->dynamic) {
                std::cout << "Start merging traversal" << std::endl;
                timer.reset();

                graph = const_cast<DBG_succ*>(graphs.at(0));
                graphs.erase(graphs.begin());

                for (size_t i = 0; i < graphs.size(); ++i) {
                    graph->merge(*graphs[i]);
                    std::cout << "traversal " << files[i + 1] << " done\t"
                              << timer.elapsed() << "sec" << std::endl;
                }
                std::cout << "Graphs merged\t" << timer.elapsed()
                                               << "sec" << std::endl;
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
                    graph = new DBG_succ(graphs[0]->get_k());
                    chunk->initialize_graph(graph);
                    std::cout << "Graphs merged\t" << timer.elapsed()
                              << "sec" << std::endl;
                }
                delete chunk;
            } else {
                std::cout << "Start merging graphs" << std::endl;
                timer.reset();

                graph = merge::merge(graphs, config->verbose);

                std::cout << "Graphs merged\t" << timer.elapsed()
                          << "sec" << std::endl;
            }

            for (auto *graph_ : graphs) {
                delete graph_;
            }

            std::cout << "... done merging." << std::endl;

            // graph output
            if (graph && config->print_graph_succ)
                std::cout << *graph;
            if (graph && !config->outfbase.empty())
                graph->serialize(config->outfbase);

            return 0;
        }
        case Config::STATS: {
            for (const auto &file : files) {
                std::unique_ptr<DBG_succ> graph {
                    load_critical_graph_from_file(file)
                };

                std::cout << "Statistics for graph " << file << std::endl;
                std::cout << "nodes: " << graph->num_nodes() << std::endl;
                std::cout << "edges: " << graph->num_edges() << std::endl;
                std::cout << "k: " << graph->get_k() << std::endl;
                std::cout << "state: " << Config::state_to_string(graph->get_state())
                          << std::endl;

                assert(graph->rank_W(graph->num_edges(), DBG_succ::alph_size) == 0);
                std::cout << "W stats: {'" << DBG_succ::decode(0) << "': "
                          << graph->rank_W(graph->num_edges(), 0);
                for (TAlphabet i = 1; i < DBG_succ::alph_size; ++i) {
                    std::cout << ", '" << DBG_succ::decode(i) << "': "
                              << graph->rank_W(graph->num_edges(), i)
                                    + graph->rank_W(graph->num_edges(), i + DBG_succ::alph_size);
                }
                std::cout << "}" << std::endl;

                assert(graph->get_F(0) == 0);
                std::cout << "F stats: {'";
                for (TAlphabet i = 1; i < DBG_succ::alph_size; ++i) {
                    std::cout << DBG_succ::decode(i - 1) << "': "
                              << graph->get_F(i) - graph->get_F(i - 1)
                              << ", '";
                }
                std::cout << DBG_succ::decode(DBG_succ::alph_size - 1) << "': "
                          << graph->num_edges() - graph->get_F(DBG_succ::alph_size - 1)
                          << "}" << std::endl;

                if (config->print_graph_internal_repr)
                    graph->print_internal_representation(std::cout);

                if (config->print_graph_succ)
                    std::cout << *graph;
            }

            for (const auto &file : config->infbase_annotators) {
                std::unique_ptr<Annotator> annotation;

                if (config->use_row_annotator) {
                    annotation.reset(new annotate::RowCompressed<>(0, config->sparse));
                } else {
                    annotation.reset(
                        new annotate::ColumnCompressed<>(
                            0, kNumCachedColumns, config->verbose
                        )
                    );
                }

                if (!annotation->load(file)) {
                    std::cerr << "ERROR: can't load annotation from file "
                              << file << std::endl;
                    exit(1);
                }

                std::cout << "Statistics for annotation " << file << std::endl;
                std::cout << "labels: " << annotation->num_labels() << std::endl;
                std::cout << "sparsity: " << std::scientific
                                          << annotation->sparsity() << std::endl;
            }

            return 0;
        }
        case Config::TRANSFORM_ANNOTATION: {
            Timer timer;

            auto annotator = std::make_unique<annotate::ColumnCompressed<>>(
                0, kNumCachedColumns, config->verbose
            );

            if (config->verbose)
                std::cout << "Loading annotator...\t" << std::flush;

            if (!annotator->merge_load(files)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }
            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            if (config->rename_instructions_file.size()) {
                if (config->verbose)
                    std::cout << "Renaming...\t" << std::flush;

                std::map<std::string, std::string> dict;
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
                annotator->rename_columns(dict);

                annotator->serialize(config->outfbase);
                if (config->verbose)
                    std::cout << timer.elapsed() << "sec" << std::endl;
            }

            if (config->to_row_annotator) {
                if (config->verbose)
                    std::cout << "Converting...\t" << std::flush;

                annotate::RowCompressed<> row_annotator(0);
                annotator->convert_to_row_annotator(&row_annotator,
                                                    config->parallel);
                annotator.reset();

                row_annotator.serialize(config->outfbase);
                if (config->verbose)
                    std::cout << timer.elapsed() << "sec" << std::endl;
            }
            return 0;
        }
        case Config::TRANSFORM: {
            Timer timer;
            if (config->verbose) {
                std::cout << "Graph loading...\t" << std::flush;
            }
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(files.at(0))
            };
            if (config->verbose) {
                std::cout << timer.elapsed() << "sec" << std::endl;
            }

            if (config->clear_dummy) {
                if (config->verbose) {
                    std::cout << "Traverse source dummy edges and remove redundant ones..." << std::endl;
                }
                timer.reset();
                graph->erase_redundant_dummy_edges(config->parallel, config->verbose);
                if (config->verbose) {
                    std::cout << "Done in " << timer.elapsed() << "sec" << std::endl;
                }
            }

            if (config->to_fasta) {
                if (config->verbose) {
                    std::cout << "Extracting sequences from graph...\t" << std::flush;
                }
                timer.reset();
                if (config->outfbase.size()) {

                    auto out_filename
                        = utils::remove_suffix(config->outfbase, ".gz", ".fasta")
                            + ".fasta.gz";

                    gzFile out_fasta_gz = gzopen(out_filename.c_str(), "w");

                    if (out_fasta_gz == Z_NULL) {
                        std::cerr << "ERROR: Can't write to " << out_filename << std::endl;
                        exit(1);
                    }

                    graph->call_sequences([&](const auto &sequence) {
                        if (!write_fasta(out_fasta_gz, "", sequence)) {
                            std::cerr << "ERROR: Can't write extracted sequences to "
                                      << out_filename << std::endl;
                            exit(1);
                        }
                    });

                    gzclose(out_fasta_gz);

                } else {
                    graph->call_sequences([&](const auto &sequence) {
                        std::cout << sequence << std::endl;
                    });
                }

                if (config->verbose) {
                    std::cout << timer.elapsed() << "sec" << std::endl;
                }
                return 0;
            }
            if (config->to_adj_list) {
                if (config->verbose) {
                    std::cout << "Converting graph to adjacency list...\t" << std::flush;
                }
                timer.reset();
                if (config->outfbase.size()) {
                    std::ofstream outstream(config->outfbase + ".adjlist");
                    graph->print_adj_list(outstream);
                } else {
                    graph->print_adj_list(std::cout);
                }
                if (config->verbose) {
                    std::cout << timer.elapsed() << "sec" << std::endl;
                }
                return 0;
            }

            if (config->verbose) {
                std::cout << "Converting graph to state "
                          << Config::state_to_string(config->state)
                          << "...\t" << std::flush;
                timer.reset();
            }

            graph->switch_state(config->state);

            if (config->verbose) {
                std::cout << timer.elapsed() << "sec" << std::endl;
            }

            if (config->outfbase.size()) {
                if (config->verbose) {
                    std::cout << "Serializing transformed graph...\t" << std::flush;
                    timer.reset();
                }
                graph->serialize(config->outfbase);
                if (config->verbose) {
                    std::cout << timer.elapsed() << "sec" << std::endl;
                }
            }

            return 0;
        }
        case Config::ALIGN: {
            assert(config->infbase.size());

            // load graph
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(config->infbase)
            };

            if (!config->alignment_length) {
                if (config->distance > 0) {
                    config->alignment_length = graph->get_k();
                } else {
                    config->alignment_length = graph->get_k() + 1;
                }
            }

            if (config->alignment_length > graph->get_k() + 1) {
                // std::cerr << "ERROR: Alignment patterns longer than"
                //           << " the graph edge size are not supported."
                //           << " Decrease the alignment length." << std::endl;
                // exit(1);
                std::cerr << "Warning: Aligning to k-mers"
                          << " longer than k+1 is not supported." << std::endl;
                config->alignment_length = graph->get_k() + 1;
            }
            if (config->distance > 0
                    && config->alignment_length != graph->get_k()) {
                std::cerr << "Warning: Aligning to k-mers longer or shorter than k"
                          << " is not supported for fuzzy alignment." << std::endl;
                config->alignment_length = graph->get_k();
            }

            std::unique_ptr<Timer> timer;
            if (config->verbose) {
                timer.reset(new Timer());
            }

            std::cout << "Align sequences against the de Bruijn graph with ";
            std::cout << "k=" << graph->get_k() << "\n"
                      << "Length of aligning k-mers: "
                      << config->alignment_length << std::endl;

            for (const auto &file : files) {
                std::cout << "Align sequences from file " << file << std::endl;

                read_fasta_file_critical(file, [&](kseq_t *read_stream) {
                    if (config->distance > 0) {
                        uint64_t aln_len = read_stream->seq.l;

                        // since we set aln_len = read_stream->seq.l, we only get a single hit vector
                        auto graphindices = graph->align_fuzzy(
                            std::string(read_stream->seq.s, read_stream->seq.l),
                            aln_len,
                            config->distance
                        );
                        //for (size_t i = 0; i < graphindices.size(); ++i) {
                        size_t i = 0;

                        int print_len = i + aln_len < read_stream->seq.l
                                        ? aln_len
                                        : (read_stream->seq.l - i);

                        printf("%.*s: ", print_len, read_stream->seq.s + i);

                        for (size_t l = 0; l < graphindices.at(i).size(); ++l) {
                            HitInfo curr_hit(graphindices.at(i).at(l));
                            for (size_t m = 0; m < curr_hit.path.size(); ++m) {
                                std::cout << curr_hit.path.at(m) << ':';
                            }
                            for (size_t m = curr_hit.rl; m <= curr_hit.ru; ++m) {
                                std::cout << m << " ";
                            }
                            std::cout << "[" << curr_hit.cigar << "] ";
                        }
                        //}
                        std::cout << std::endl;

                        return;
                    }

                    // Non-fuzzy mode

                    if (config->verbose) {
                        std::cout << "Sequence: " << read_stream->seq.s << "\n";
                    }

                    if (config->query_presence
                            && config->alignment_length == graph->get_k() + 1) {
                        if (config->filter_present) {
                            if (graph->find(read_stream->seq.s,
                                            config->discovery_fraction,
                                            config->kmer_mapping_mode))
                                std::cout << ">" << read_stream->name.s << "\n"
                                                 << read_stream->seq.s << "\n";
                        } else {
                            std::cout << graph->find(read_stream->seq.s,
                                                     config->discovery_fraction,
                                                     config->kmer_mapping_mode) << "\n";
                        }
                        return;
                    }

                    assert(config->alignment_length <= graph->get_k() + 1);

                    std::vector<uint64_t> graphindices;

                    if (config->alignment_length == graph->get_k() + 1) {
                        graphindices = graph->map_to_edges(read_stream->seq.s);
                    } else {
                        graphindices = graph->map_to_nodes(read_stream->seq.s,
                                                           config->alignment_length);
                    }

                    size_t num_discovered = std::count_if(
                        graphindices.begin(), graphindices.end(),
                        [](const auto &x) { return x > 0; }
                    );

                    const size_t num_kmers = graphindices.size();

                    if (config->query_presence) {
                        const size_t min_kmers_discovered =
                            num_kmers - num_kmers * (1 - config->discovery_fraction);
                        if (config->filter_present) {
                            if (num_discovered >= min_kmers_discovered)
                                std::cout << ">" << read_stream->name.s << "\n"
                                                 << read_stream->seq.s << "\n";
                        } else {
                            std::cout << (num_discovered >= min_kmers_discovered) << "\n";
                        }
                        return;
                    }

                    if (config->count_kmers_query) {
                        std::cout << num_discovered << "/" << num_kmers << "\n";
                        return;
                    }

                    for (size_t i = 0; i < graphindices.size(); ++i) {
                        assert(i + config->alignment_length <= read_stream->seq.l);
                        std::cout << std::string(read_stream->seq.s + i, config->alignment_length);
                        std::cout << ": " << graphindices[i] << "\n";
                    }
                }, config->reverse, timer.get(),
                    get_filter_filename(file, config->filter_k,
                                        config->max_unreliable_abundance,
                                        config->unreliable_kmers_threshold)
                );
            }

            return 0;
        }
        case Config::NO_IDENTITY: {
            assert(false);
            break;
        }
    }

    return 0;
}
