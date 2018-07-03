#include <zlib.h>

#include "dbg_succinct.hpp"
#include "dbg_succinct_chunk.hpp"
#include "dbg_succinct_construct.hpp"
#include "reads_filtering.hpp"
#include "config.hpp"
#include "helpers.hpp"
#include "utils.hpp"
#include "vcf_parser.hpp"
#include "traverse.hpp"
#include "dbg_succinct_merge.hpp"
#include "annotate_color_compressed.hpp"
#include "annotate_color_compressed_fast.hpp"
#include "annotate_row_compressed.hpp"
#include "annotate_bloom_filter.hpp"
#include "unix_tools.hpp"
#include "kmer.hpp"
#include "serialization.hpp"

typedef annotate::MultiColorAnnotation<uint64_t, std::string> Annotator;

const size_t kMaxNumParallelReadFiles = 5;

const size_t kNumCachedColors = 10;

const char kDefaultFastQualityChar = 'I';

KSEQ_INIT(gzFile, gzread);


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


bool write_fasta(gzFile gz_out, const kseq_t &kseq) {
    return gzputc(gz_out, '>') == '>'
        && gzwrite(gz_out, kseq.name.s, kseq.name.l)
                        == static_cast<int>(kseq.name.l)
        && (!kseq.comment.l
            || (gzputc(gz_out, ' ') == ' '
                && gzwrite(gz_out, kseq.comment.s, kseq.comment.l)
                                    == static_cast<int>(kseq.comment.l)))
        && gzputc(gz_out, '\n') == '\n'
        && gzwrite(gz_out, kseq.seq.s, kseq.seq.l)
                        == static_cast<int>(kseq.seq.l)
        && gzputc(gz_out, '\n') == '\n';
}

bool write_fastq(gzFile gz_out, const kseq_t &kseq) {
    std::string qual(kseq.qual.s, kseq.qual.l);

    if (!kseq.qual.l && kseq.seq.l)
        qual.assign(kseq.seq.l, kDefaultFastQualityChar);

    return gzputc(gz_out, '@') == '@'
        && gzwrite(gz_out, kseq.name.s, kseq.name.l)
                        == static_cast<int>(kseq.name.l)
        && (kseq.comment.l == 0
            || (gzputc(gz_out, ' ') == ' '
                && gzwrite(gz_out, kseq.comment.s, kseq.comment.l)
                                    == static_cast<int>(kseq.comment.l)))
        && gzputc(gz_out, '\n') == '\n'
        && gzwrite(gz_out, kseq.seq.s, kseq.seq.l)
                        == static_cast<int>(kseq.seq.l)
        && gzputc(gz_out, '\n') == '\n'
        && gzputc(gz_out, '+') == '+'
        && gzwrite(gz_out, kseq.name.s, kseq.name.l)
                        == static_cast<int>(kseq.name.l)
        && (kseq.comment.l == 0
            || (gzputc(gz_out, ' ') == ' '
                && gzwrite(gz_out, kseq.comment.s, kseq.comment.l)
                                    == static_cast<int>(kseq.comment.l)))
        && gzputc(gz_out, '\n') == '\n'
        && gzwrite(gz_out, qual.data(), qual.size())
                        == static_cast<int>(qual.size())
        && gzputc(gz_out, '\n') == '\n';
}


template <class Callback>
void read_fasta_file_critical(const std::string &filename,
                              Callback callback,
                              bool with_reverse = false,
                              Timer *timer = NULL,
                              const std::string &filter_filename = "") {
    std::vector<bool> filter;
    if (filter_filename.size()) {
        std::ifstream instream(filter_filename);
        try {
            filter = load_number_vector<bool>(instream);
        } catch (...) {
            std::cerr << "ERROR: Filter file " << filter_filename
                      << " is corrupted" << std::endl;
            exit(1);
        }
    }
    size_t seq_count = 0;

    gzFile input_p = gzopen(filename.c_str(), "r");
    if (input_p == Z_NULL) {
        std::cerr << "ERROR no such file " << filename << std::endl;
        exit(1);
    }
    //TODO: handle read_stream->qual
    kseq_t *read_stream = kseq_init(input_p);
    if (read_stream == NULL) {
        std::cerr << "ERROR while opening input file " << filename << std::endl;
        exit(1);
    }
    if (timer) {
        std::cout << "Start extracting sequences from file " << filename << std::endl;
    }
    while (kseq_read(read_stream) >= 0) {
        if (filter_filename.size() && filter.size() <= seq_count) {
            std::cerr << "ERROR: Filter file " << filter_filename
                      << " has fewer sequences" << std::endl;
            exit(1);
        }

        if (!filter_filename.size() || filter[seq_count]) {
            callback(read_stream);
            if (with_reverse) {
                reverse_complement(read_stream->seq);
                callback(read_stream);
            }
        }

        seq_count++;
    }
    if (filter_filename.size() && filter.size() != seq_count) {
        std::cerr << "ERROR: Filter file " << filter_filename
                  << " has more sequences" << std::endl;
        exit(1);
    }
    if (timer) {
        std::cout << "Finished extracting sequences from file " << filename
                  << " in " << timer->elapsed() << "sec"
                  << ", sequences extracted: " << seq_count << std::endl;
    }
    kseq_destroy(read_stream);
    gzclose(input_p);
}


template <class Callback>
void read_vcf_file_critical(const std::string &filename,
                            const std::string &ref_filename,
                            size_t k,
                            std::vector<std::string> *annotation,
                            Callback callback, Timer *timer = NULL) {

    //TODO: make this a configurable option
    //default list of tokens to extract as annotations
    //TODO: extract these guys directly from vcf parsed
    const std::vector<std::string> annots = {
      "AC_AFR", "AC_EAS", "AC_AMR", "AC_ASJ",
      "AC_FIN", "AC_NFE", "AC_SAS", "AC_OTH"
    };

    vcf_parser vcf;
    if (!vcf.init(ref_filename, filename, k)) {
        std::cerr << "ERROR reading VCF " << filename << std::endl;
        exit(1);
    }
    if (timer) {
        std::cout << "Extracting sequences from file " << filename << std::endl;
    }
    size_t seq_count = 0;
    while (vcf.get_seq(annots, annotation)) {
        callback(vcf.seq, annotation);
        seq_count++;
    }
    if (timer) {
        std::cout << "Finished extracting sequences from file " << filename
                  << " in " << timer->elapsed() << "sec"
                  << ", sequences extracted: " << seq_count << std::endl;
    }
}


template <typename... Args>
void annotate_sequence(const std::string &sequence,
                       const std::vector<std::string> &labels,
                       const DBG_succ &graph,
                       annotate::AnnotationCategoryBloom *annotator,
                       utils::ThreadPool *thread_pool,
                       std::mutex *annotation_mutex) {
    if (sequence.size() < graph.get_k())
        return;

    std::ignore = annotation_mutex;
    std::ignore = thread_pool;

    annotator->add_colors(sequence, labels, graph.num_edges());
}


void annotate_sequence_thread_safe(
        std::string sequence,
        std::vector<std::string> labels,
        const DBG_succ *graph,
        annotate::MultiColorAnnotation<uint64_t, std::string> *annotator,
        std::mutex *annotation_mutex) {

    std::vector<uint64_t> indices;

    graph->align(sequence, [&](uint64_t i) {
        if (i > 0)
            indices.push_back(i);
    });

    assert(annotation_mutex);

    std::lock_guard<std::mutex> lock(*annotation_mutex);
    annotator->add_colors(indices, labels);
}


void annotate_sequence(const std::string &sequence,
                       const std::vector<std::string> &labels,
                       const DBG_succ &graph,
                       annotate::MultiColorAnnotation<uint64_t, std::string> *annotator,
                       utils::ThreadPool *thread_pool,
                       std::mutex *annotation_mutex) {
    if (sequence.size() < graph.get_k())
        return;

    if (thread_pool && annotation_mutex) {
        thread_pool->enqueue(annotate_sequence_thread_safe,
                             sequence, labels,
                             &graph, annotator, annotation_mutex);
        return;
    }

    std::vector<uint64_t> indices;
    graph.align(sequence, [&](uint64_t i) {
        if (i > 0)
            indices.push_back(i);
    });
    annotator->add_colors(indices, labels);
}


std::string get_filter_filename(const std::string &filename,
                                size_t k,
                                size_t max_unreliable_abundance,
                                size_t unreliable_kmers_threshold,
                                bool critical = true) {
    std::string filter_filename = max_unreliable_abundance > 0
            ? filename + ".filter_k" + std::to_string(k)
                       + "_s" + std::to_string(max_unreliable_abundance)
                       + (max_unreliable_abundance
                            ? std::string("_") + std::to_string(unreliable_kmers_threshold)
                            : "")
            : "";

    if (critical
            && filter_filename.size()
            && !std::ifstream(filter_filename).good()) {
        std::cerr << "ERROR: read filter "
                  << filter_filename << " does not exist."
                  << " Filter reads first." << std::endl;
        exit(1);
    }

    return filter_filename;
}


template <class Annotator, typename... Args>
void annotate_data(const std::vector<std::string> &files,
                   const std::string &ref_sequence_path,
                   const DBG_succ &graph,
                   Annotator *annotator,
                   bool reverse,
                   size_t max_unreliable_abundance,
                   size_t unreliable_kmers_threshold,
                   bool filename_anno,
                   bool fasta_anno,
                   const std::string &fasta_header_delimiter,
                   const std::vector<std::string> &anno_labels,
                   size_t genome_bin_size,
                   bool verbose,
                   utils::ThreadPool *thread_pool = NULL,
                   std::mutex *annotation_mutex = NULL) {
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
                                   graph.get_k(),
                                   &variant_labels,
                [&](std::string &seq, std::vector<std::string> *variant_labels) {
                    assert(variant_labels);

                    if (filename_anno)
                        variant_labels->push_back(file);

                    for (const auto &label : anno_labels) {
                        variant_labels->push_back(label);
                    }

                    annotate_sequence(seq, *variant_labels, graph, annotator,
                                      thread_pool, annotation_mutex);
                    if (reverse) {
                        reverse_complement(seq.begin(), seq.end());
                        annotate_sequence(seq, *variant_labels, graph, annotator,
                                          thread_pool, annotation_mutex);
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

                    if (!genome_bin_size) {
                        annotate_sequence(read_stream->seq.s, labels, graph, annotator,
                                          thread_pool, annotation_mutex);
                    } else {
                        const std::string sequence(read_stream->seq.s);
                        labels.push_back("");

                        for (size_t i = 0; i < sequence.size(); i += genome_bin_size) {
                            labels.back() = std::string(read_stream->name.s)
                                                + "_" + std::to_string(i);
                            annotate_sequence(
                                sequence.substr(i, genome_bin_size + graph.get_k() - 1),
                                labels, graph, annotator, thread_pool, annotation_mutex
                            );
                        }
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
                },
                reverse, timer.get(),
                get_filter_filename(file, graph.get_k(),
                                    max_unreliable_abundance, unreliable_kmers_threshold)
            );
        } else {
            std::cerr << "ERROR: Filetype unknown for file "
                      << file << std::endl;
            exit(1);
        }
    }
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
                   const DBG_succ *graph,
                   const Annotator *annotation) {
    std::ostringstream oss;
    if (!suppress_unlabeled)
        oss << seq_name << "\t";

    if (count_labels) {
        auto top_labels = annotation->get_most_frequent_colors(
            graph->index(sequence),
            num_top_labels
        );

        if (top_labels.size()) {
            oss << "<" << top_labels[0].first << ">:" << top_labels[0].second;
        }
        for (size_t i = 1; i < top_labels.size(); ++i) {
            oss << "\t<" << top_labels[i].first << ">:" << top_labels[i].second;
        }
        oss << "\n";
    } else {
        std::vector<uint64_t> indices;
        size_t num_missing_kmers = 0;

        graph->align(sequence, [&](uint64_t i) {
            if (i > 0) {
                indices.push_back(i);
            } else {
                num_missing_kmers++;
            }
        });

        std::vector<std::string> labels_discovered;
        // discovery_fraction = discovered / all
        // new_discovered_fraction = discovered / (all - missing)
        if (indices.size() > 0
                && discovery_fraction * (indices.size() + num_missing_kmers)
                                                            <= indices.size()) {
            labels_discovered = annotation->aggregate_colors(
                indices,
                discovery_fraction
                    * (indices.size() + num_missing_kmers)
                    / indices.size()
            );
        }

        if (!suppress_unlabeled || labels_discovered.size()) {
            if (suppress_unlabeled)
                oss << seq_name << "\t";

            oss << utils::join_strings(labels_discovered,
                                       anno_labels_delimiter) << "\n";
        }
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
                    static_cast<uint64_t>(std::ceil(std::log2(config->nsplits)
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

                    if (valid_kmer_suffix(suffix)) {
                        std::cout << "\nSuffix: " << suffix << std::endl;
                    } else {
                        std::cout << "\nSkipping suffix: " << suffix << std::endl;
                        continue;
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
                                [&](std::string &seq, auto *variant_annotations) {
                                    constructor->add_read(seq);
                                    if (config->reverse) {
                                        reverse_complement(seq.begin(), seq.end());
                                        constructor->add_read(seq);
                                    }
                                    std::ignore = variant_annotations;
                                }, timer_ptr
                            );
                        } else if (utils::get_filetype(files[f]) == "FASTA"
                                    || utils::get_filetype(files[f]) == "FASTQ") {
                            Timer *timer_ptr = config->verbose ? &timer : NULL;

                            std::string filter_filename = get_filter_filename(
                                files[f], config->k,
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
                            [&](std::string &seq, auto *variant_annotations) {
                                graph->add_sequence(seq);
                                if (config->reverse) {
                                    reverse_complement(seq.begin(), seq.end());
                                    graph->add_sequence(seq);
                                }
                                std::ignore = variant_annotations;
                            }
                        );
                    } else if (utils::get_filetype(files[f]) == "FASTA"
                                || utils::get_filetype(files[f]) == "FASTQ") {
                        read_fasta_file_critical(files[f],
                            [&](kseq_t *read_stream) {
                                graph->add_sequence(read_stream->seq.s);
                            },
                            config->reverse, NULL,
                            get_filter_filename(files[f], config->k,
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
                graph->print_state();
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

            if (graph->state != Config::DYN) {
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
                        [&](std::string &seq, auto *variant_annotations) {
                            graph->add_sequence(seq, true, inserted_edges.get());
                            if (config->reverse) {
                                reverse_complement(seq.begin(), seq.end());
                                graph->add_sequence(seq, true, inserted_edges.get());
                            }
                            std::ignore = variant_annotations;
                        }
                    );
                } else if (utils::get_filetype(file) == "FASTA"
                            || utils::get_filetype(file) == "FASTQ") {
                    read_fasta_file_critical(file,
                        [&](kseq_t *read_stream) {
                            graph->add_sequence(read_stream->seq.s, true, inserted_edges.get());
                        },
                        config->reverse, NULL,
                        get_filter_filename(file, config->k,
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
                graph->print_state();

            assert(config->outfbase.size());

            // serialize graph
            timer.reset();
            if (graph->state != config->state) {
                if (config->verbose)
                    std::cout << "Switching state before dumping..." << std::flush;

                graph->switch_state(config->state);

                if (config->verbose)
                    std::cout << "\tdone in " << timer.elapsed() << "sec" << std::endl;
            }

            graph->serialize(config->outfbase);

            timer.reset();

            if (config->infbase_annotators.size()) {
                std::unique_ptr<Annotator> annotation;

                if (config->use_row_annotator) {
                    annotation.reset(
                        new annotate::RowCompressed<>(graph->num_edges() + 1,
                                                      config->sparse)
                    );
                } else {
                    annotation.reset(
                        new annotate::ColorCompressed<>(graph->num_edges() + 1,
                                                        kNumCachedColors,
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

                assert(inserted_edges.get());

                std::vector<uint64_t> inserted_edge_idx;
                for (uint64_t j = 1; j <= inserted_edges->get_num_set_bits(); ++j) {
                    inserted_edge_idx.push_back(inserted_edges->select1(j));
                }
                inserted_edges.reset();

                if (config->verbose)
                    std::cout << "Insert empty rows to the annotation matrix..." << std::flush;

                annotation->insert_rows(inserted_edge_idx);

                if (config->verbose)
                    std::cout << "\tdone in " << timer.elapsed() << "sec" << std::endl;

                annotation->serialize(config->outfbase);
            }

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
                    config->k,
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
                            use_kmc ? file + ".kmc" : ""
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
                    config->k,
                    config->max_unreliable_abundance,
                    config->unreliable_kmers_threshold,
                    config->verbose, config->reverse, config->use_kmc
                );
            }
            thread_pool_files.join();
            thread_pool.join();

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
                    new annotate::ColorCompressed<>(graph->num_edges() + 1,
                                                    kNumCachedColors,
                                                    config->verbose)
                );
            }

            if (config->infbase_annotators.size()
                    && !annotation->merge_load(config->infbase_annotators)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }

            std::unique_ptr<utils::ThreadPool> thread_pool;
            std::unique_ptr<std::mutex> annotation_mutex;

            if (config->parallel > 1) {
                thread_pool.reset(
                    new utils::ThreadPool(config->parallel - 1)
                );
                annotation_mutex.reset(new std::mutex());
            }

            annotate_data(files,
                          config->refpath,
                          *graph,
                          annotation.get(),
                          config->reverse,
                          config->max_unreliable_abundance,
                          config->unreliable_kmers_threshold,
                          config->filename_anno,
                          config->fasta_anno,
                          config->fasta_header_delimiter,
                          config->anno_labels,
                          config->genome_binsize_anno,
                          config->verbose,
                          thread_pool.get(),
                          annotation_mutex.get());

            // join threads if any were initialized
            thread_pool.reset();

            annotation->serialize(config->outfbase);

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

            annotate_data(files,
                          config->refpath,
                          *graph,
                          annotation.get(),
                          config->reverse,
                          config->max_unreliable_abundance,
                          config->unreliable_kmers_threshold,
                          config->filename_anno,
                          config->fasta_anno,
                          config->fasta_header_delimiter,
                          config->anno_labels,
                          config->genome_binsize_anno,
                          config->verbose);

            annotation->serialize(config->infbase);

            return 0;
        }
        case Config::MERGE_ANNOTATORS: {
            // load graph
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(config->infbase)
            };

            std::unique_ptr<Annotator> annotation;
            if (config->use_row_annotator) {
                throw std::runtime_error("To be implemented");
                annotation.reset(new annotate::RowCompressed<>(0, config->sparse));
            } else {
                annotation.reset(
                    new annotate::ColorCompressed<>(
                        0, kNumCachedColors, config->verbose
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
                    new annotate::FastColorCompressed<>(
                        0, kNumCachedColors, config->verbose
                    )
                );
            } else {
                annotation.reset(
                    new annotate::ColorCompressed<>(
                        0, kNumCachedColors, config->verbose
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
                            graph.get(),
                            annotation.get()
                        );
                    },
                    config->reverse, timer.get(),
                    get_filter_filename(file, graph->get_k(),
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
            std::cout << "Opening file " << files.at(0) << std::endl;
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(files.at(0))
            };

            for (size_t f = 1; f < files.size(); ++f) {
                std::cout << "Opening file for comparison ..." << files[f] << std::endl;
                DBG_succ *second = load_critical_graph_from_file(files[f]);
                if (*graph == *second) {
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
                std::cout << "state: " << graph->state << std::endl;
            }

            // graph output
            if (graph && config->print_graph_succ)
                graph->print_state();
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
                graph->print_state();
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
                std::cout << "state: " << graph->state << std::endl;

                if (config->print_graph_succ)
                    graph->print_state();
            }

            for (const auto &file : config->infbase_annotators) {
                std::unique_ptr<Annotator> annotation;

                if (config->use_row_annotator) {
                    annotation.reset(new annotate::RowCompressed<>(0, config->sparse));
                } else {
                    annotation.reset(
                        new annotate::ColorCompressed<>(
                            0, kNumCachedColors, config->verbose
                        )
                    );
                }

                if (!annotation->load(file)) {
                    std::cerr << "ERROR: can't load annotation from file "
                              << file << std::endl;
                    exit(1);
                }

                std::cout << "Statistics for annotation " << file << std::endl;
                std::cout << "colors: " << annotation->num_colors() << std::endl;
                std::cout << "sparsity: " << std::scientific
                                          << annotation->sparsity() << std::endl;
            }

            return 0;
        }
        case Config::TRANSFORM: {
            Timer timer;
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(files.at(0))
            };

            if (config->to_sequences) {
                if (config->verbose) {
                    std::cout << "Extracting sequences from graph...\t" << std::flush;
                }
                timer.reset();
                if (config->outfbase.size()) {
                    std::ofstream outstream(config->outfbase + ".adjlist");
                    graph->call_sequences([&](const auto &sequence) {
                        outstream << sequence << std::endl;
                    });
                } else {
                    graph->call_sequences([&](const auto &sequence) {
                        std::cout << sequence << std::endl;
                    });
                }
                if (config->verbose) {
                    std::cout << timer.elapsed() << "sec" << std::endl;
                }
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
            }
            if (config->sqlfbase.size()) {
                if (config->verbose) {
                    std::cout << "Converting graph to SQL DB...\t" << std::flush;
                }
                timer.reset();
                traverse::toSQL(graph.get(), files, config->sqlfbase);
                if (config->verbose) {
                    std::cout << timer.elapsed() << "sec" << std::endl;
                }
            }
            if (graph->state != config->state) {
                if (config->verbose) {
                    std::cout << "Converting graph to state " << config->state
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
            }

            return 0;
        }
        case Config::ALIGN: {
            assert(config->infbase.size());

            // load graph
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(config->infbase)
            };

            if (!config->alignment_length
                    || config->alignment_length > graph->get_k()) {
                config->alignment_length = graph->get_k();
            }

            std::unique_ptr<Timer> timer;
            if (config->verbose) {
                timer.reset(new Timer());
            }

            std::cout << "Align sequences against a de Bruijn graph with ";
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

                    if (config->verbose) {
                        std::cout << "Sequence: " << read_stream->seq.s << "\n";
                    }

                    if (config->query_presence
                            && config->alignment_length == graph->get_k()) {
                        if (config->filter_present) {
                            if (graph->find(read_stream->seq.s,
                                            config->discovery_fraction))
                                std::cout << ">" << read_stream->name.s << "\n"
                                                 << read_stream->seq.s << "\n";
                        } else {
                            std::cout << graph->find(read_stream->seq.s,
                                                     config->discovery_fraction) << "\n";
                        }
                        return;
                    }

                    auto graphindices = graph->index(read_stream->seq.s,
                                                     config->alignment_length);

                    size_t num_discovered = std::count_if(
                        graphindices.begin(), graphindices.end(),
                        [](const auto &x) { return x > 0; }
                    );

                    if (config->query_presence) {
                        const size_t num_kmers = read_stream->seq.l - graph->get_k() + 1;
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
                        std::cout << num_discovered << "/"
                                  << read_stream->seq.l
                                        - config->alignment_length + 1 << "\n";
                        return;
                    }

                    for (size_t i = 0; i < graphindices.size(); ++i) {
                        for (uint64_t j = 0; j < config->alignment_length; ++j) {
                            std::cout << read_stream->seq.s[i + j];
                        }
                        std::cout << ": " << graphindices[i] << "\n";
                    }
                }, config->reverse, timer.get(),
                   get_filter_filename(file, graph->get_k(),
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
