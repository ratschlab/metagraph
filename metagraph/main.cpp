#include <zlib.h>

#include "dbg_succinct.hpp"
#include "dbg_succinct_chunk.hpp"
#include "dbg_succinct_construct.hpp"
#include "config.hpp"
#include "helpers.hpp"
#include "utils.hpp"
#include "vcf_parser.hpp"
#include "traverse.hpp"
#include "dbg_succinct_merge.hpp"
#include "annotate_color_compressed.hpp"
#include "annotate_row_compressed.hpp"
#include "annotate_bloom_filter.hpp"
#include "unix_tools.hpp"
#include "kmer.hpp"

using libmaus2::util::NumberSerialisation;
typedef annotate::AnnotationCategory<std::set<std::string>> Annotator;

const size_t kMaxNumParallelReadFiles = 5;

const size_t kNumCachedColors = 10;

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


template <class Callback>
void read_fasta_file_critical(const std::string &filename,
                              Callback callback, Timer *timer = NULL,
                              const std::string &filter_filename = "") {
    bit_vector_stat filter;
    if (filter_filename.size()) {
        std::ifstream instream(filter_filename);
        if (!instream.good() || !filter.deserialise(instream)) {
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
        if (!filter_filename.size() || filter[seq_count])
            callback(read_stream);

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

template <class Annotator, typename... Args>
void annotate_data(const std::vector<std::string> &files,
                   const std::string &ref_sequence_path,
                   const DBG_succ &graph,
                   Annotator *annotator,
                   bool reverse,
                   bool fasta_anno,
                   const std::string &fasta_header_delimiter,
                   bool verbose,
                   Args... args) {
    size_t total_seqs = 0;

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

                    variant_labels->push_back(file);

                    annotator->add_labels(seq,
                        std::set<std::string>(variant_labels->begin(),
                                              variant_labels->end()),
                        args...
                    );
                    if (reverse) {
                        reverse_complement(seq.begin(), seq.end());
                        annotator->add_labels(seq,
                            std::set<std::string>(variant_labels->begin(),
                                                  variant_labels->end()),
                            args...
                        );
                    }
                    variant_labels->clear();
                }
            );
        } else if (utils::get_filetype(file) == "FASTA"
                    || utils::get_filetype(file) == "FASTQ") {
            read_fasta_file_critical(file, [&](kseq_t *read_stream) {
                std::set<std::string> labels = { file, };

                if (fasta_anno) {
                    auto header_labels = utils::split_string(read_stream->name.s,
                                                             fasta_header_delimiter);
                    labels.insert(header_labels.begin(), header_labels.end());
                }

                if (utils::get_filetype(file) == "FASTQ") {
                    annotator->add_labels(read_stream->seq.s, labels, args...);
                    if (reverse) {
                        reverse_complement(read_stream->seq);
                        annotator->add_labels(read_stream->seq.s, labels, args...);
                    }
                } else {
                    annotator->add_labels(read_stream->seq.s, labels);
                    if (reverse) {
                        reverse_complement(read_stream->seq);
                        annotator->add_labels(read_stream->seq.s, labels);
                    }
                }

                total_seqs += 1;
                if (verbose && total_seqs % 10000 == 0) {
                    std::cout << "added labels for " << total_seqs
                              << " sequences"
                              << ", last was " << read_stream->name.s
                              << ", annotated as ";
                    for (const auto &label : labels) {
                        std::cout << "<" << label << ">";
                    }
                    std::cout << std::endl;
                }
            });
        } else {
            std::cerr << "ERROR: Filetype unknown for file "
                      << file << std::endl;
            exit(1);
        }
    }
}


template <class Annotator>
std::vector<std::string> discover_labels(const DBG_succ &graph,
                                         const Annotator &annotator,
                                         const std::vector<std::string> labels,
                                         const std::string &sequence,
                                         double discovery_fraction) {
    const size_t max_kmers_missing =
        (sequence.size() - graph.get_k() + 1)
            * (1 - discovery_fraction);

    const size_t min_kmers_discovered =
        (sequence.size() - graph.get_k() + 1)
            - max_kmers_missing;

    std::map<std::string, size_t> labels_counter;
    for (const auto &label : labels) {
        labels_counter[label] = 0;
    }

    size_t kmers_checked = 0;
    std::vector<std::string> labels_discovered;

    graph.align(sequence,
        [&](uint64_t i) {
            kmers_checked++;
            for (auto it = labels_counter.begin(); it != labels_counter.end();) {
                if (i > 0 && annotator.has_label(i, { it->first }))
                    it->second++;

                if (it->second >= min_kmers_discovered) {
                    labels_discovered.push_back(it->first);
                    labels_counter.erase(it++);
                } else if (kmers_checked - it->second > max_kmers_missing) {
                    labels_counter.erase(it++);
                } else {
                    ++it;
                }
            }
        },
        [&]() { return labels_counter.size() == 0; }
    );

    return labels_discovered;
}


std::map<std::string, size_t> count_labels(const DBG_succ &graph,
                                           const Annotator &annotator,
                                           const std::string &sequence) {
    std::map<std::string, size_t> labels_counter;

    graph.align(sequence,
        [&](uint64_t i) {
            if (i) {
                const std::set<std::string> &coloring = annotator.get(i);
                for (const std::string &label : coloring) {
                    labels_counter[label]++;
                }
            }
        }
    );

    return labels_counter;
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

            if (config->fast) {
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
                std::deque<std::string> suffices = utils::generate_strings(
                    DBG_succ::alphabet.substr(0, DBG_succ::alph_size),
                    suffix_len
                );
                if (config->suffix.size())
                    suffices = { config->suffix };

                DBG_succ::Chunk graph_data;

                //one pass per suffix
                for (const std::string &suffix : suffices) {
                    timer.reset();

                    if (suffix.size()) {
                        std::cout << "\nSuffix: " << suffix << std::endl;
                    }

                    std::unique_ptr<KMerDBGSuccChunkConstructor> constructor(
                        new KMerDBGSuccChunkConstructor(
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
                            read_vcf_file_critical(files[f], config->refpath, graph->get_k(), NULL,
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

                            if (config->noise_kmer_frequency > 0
                                    || files.size() >= config->parallel) {
                                auto reverse = config->reverse;
                                auto file = files[f];

                                std::string filter_filename = config->noise_kmer_frequency > 0
                                    ? file + ".filter_k" + std::to_string(config->k)
                                           + "_s" + std::to_string(config->noise_kmer_frequency)
                                    : "";

                                if (filter_filename.size()
                                        && !std::ifstream(filter_filename).good()) {
                                    std::cerr << "ERROR: read filter "
                                              << filter_filename << " does not exist."
                                              << " Filter reads first." << std::endl;
                                    exit(1);
                                }

                                // capture all required values by copying to be able
                                // to run task from other threads
                                constructor->add_reads([=](auto callback) {
                                    read_fasta_file_critical(file, [=](kseq_t *read_stream) {
                                        // add read to the graph constructor as a callback
                                        callback(read_stream->seq.s);
                                        if (reverse) {
                                            reverse_complement(read_stream->seq);
                                            callback(read_stream->seq.s);
                                        }
                                    }, timer_ptr, filter_filename);
                                });
                            } else {
                                read_fasta_file_critical(files[f], [&](kseq_t *read_stream) {
                                    // add read to the graph constructor as a callback
                                    constructor->add_read(read_stream->seq.s);
                                    if (config->reverse) {
                                        reverse_complement(read_stream->seq);
                                        constructor->add_read(read_stream->seq.s);
                                    }
                                }, timer_ptr);
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
                    get_RAM();
                    std::cout << "Reading data finished\t" << timer.elapsed() << "sec" << std::endl;

                    std::cout << "Sorting kmers and appending succinct"
                              << " representation from current bin...\t" << std::flush;
                    timer.reset();
                    auto next_block = constructor->build_chunk();
                    std::cout << timer.elapsed() << "sec" << std::endl;

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
                        read_vcf_file_critical(files[f], config->refpath, graph->get_k(), NULL,
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
                        read_fasta_file_critical(files[f], [&](kseq_t *read_stream) {
                            graph->add_sequence(read_stream->seq.s);
                            if (config->reverse) {
                                reverse_complement(read_stream->seq);
                                graph->add_sequence(read_stream->seq.s);
                            }
                        });
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
                auto *thread_pool_ptr = &thread_pool;
                // capture all required values by copying to be able
                // to run task from other threads
                thread_pool_files.enqueue([=](size_t k,
                                              size_t noise_kmer_frequency,
                                              bool verbose,
                                              bool reverse) {
                        // compute read filter, bit vector indicating filtered reads
                        // TODO: fix for the case of reverse complement reads
                        bit_vector_stat filter(filter_reads([=](auto callback) {
                                read_fasta_file_critical(file, [=](kseq_t *read_stream) {
                                    // add read to the graph constructor as a callback
                                    callback(read_stream->seq.s);
                                    if (reverse) {
                                        reverse_complement(read_stream->seq);
                                        callback(read_stream->seq.s);
                                    }
                                }, timer_ptr);
                            },
                            k, noise_kmer_frequency, verbose, thread_pool_ptr
                        ));

                        // dump filter
                        std::ofstream outstream(
                            file + ".filter_k" + std::to_string(k)
                                 + "_s" + std::to_string(noise_kmer_frequency)
                        );
                        filter.serialise(outstream);
                    },
                    config->k, config->noise_kmer_frequency,
                    config->verbose, config->reverse
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
                annotation.reset(new annotate::RowCompressed(*graph));
            } else {
                annotation.reset(new annotate::ColorCompressed(*graph, kNumCachedColors));
            }

            annotate_data(files,
                          config->refpath,
                          *graph,
                          annotation.get(),
                          config->reverse,
                          config->fasta_anno,
                          config->fasta_header_delimiter,
                          config->verbose);

            annotation->serialize(config->infbase);

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
                          config->fasta_anno,
                          config->fasta_header_delimiter,
                          config->verbose,
                          graph->num_edges());

            annotation->serialize(config->infbase);

            return 0;
        }
        case Config::CLASSIFY: {
            // load graph
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(config->infbase)
            };

            std::unique_ptr<Annotator> annotation;
            if (config->use_row_annotator) {
                annotation.reset(new annotate::RowCompressed(*graph));
            } else {
                annotation.reset(new annotate::ColorCompressed(*graph, kNumCachedColors));
            }

            if (!annotation->load(config->infbase)) {
                std::cerr << "ERROR: can't load annotations for "
                          << config->infbase + ".dbg"
                          << ", file corrupted" << std::endl;
                exit(1);
            }

            // iterate over input files
            for (const auto &file : files) {
                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << file << std::endl;
                }

                read_fasta_file_critical(file, [&](kseq_t *read_stream) {
                    std::cout << read_stream->name.s << ":\t";

                    if (config->reverse)
                        reverse_complement(read_stream->seq);

                    if (config->count_labels) {
                        auto labels_counter = count_labels(*graph, *annotation,
                                                           read_stream->seq.s);

                        std::vector<std::pair<std::string, size_t>> counts(
                            labels_counter.begin(), labels_counter.end()
                        );
                        // sort in decreasing order
                        std::sort(counts.begin(), counts.end(),
                                  [](const auto &first, const auto &second) {
                                      return first.second > second.second;
                                  });

                        auto num_labels = std::min(
                            counts.size(),
                            static_cast<size_t>(config->num_top_labels)
                        );
                        if (num_labels) {
                            std::cout << "<" << counts[0].first << ">: "
                                      << counts[0].second;
                        }
                        for (size_t i = 1; i < num_labels; ++i) {
                            std::cout << ", <" << counts[i].first << ">: "
                                      << counts[i].second;
                        }
                        std::cout << "\n";
                    } else {
                        auto labels_discovered = discover_labels(
                            *graph, *annotation, annotation->get_label_names(),
                            read_stream->seq.s, config->discovery_fraction
                        );

                        std::cout << utils::join_strings(labels_discovered,
                                                         config->anno_labels_delimiter)
                                  << "\n";
                    }
                });
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
        case Config::MERGE: {
            DBG_succ *graph = NULL;

            // collect results on an external merge
            if (config->collect > 1) {
                graph = DBG_succ::Chunk::build_graph_from_chunks(config->k,
                                                                 files,
                                                                 config->verbose);
                assert(graph);

                if (config->verbose) {
                    std::cout << "Succinct graph has been assembled" << std::endl;
                    std::cout << "nodes: " << graph->num_nodes() << std::endl;
                    std::cout << "edges: " << graph->num_edges() << std::endl;
                    std::cout << "k: " << graph->get_k() << std::endl;
                    std::cout << "state: " << graph->state << std::endl;
                }
            } else {
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

                if (config->traversal_merge) {
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
            }

            // graph output
            if (graph && config->print_graph_succ)
                graph->print_state();
            if (graph && !config->outfbase.empty())
                graph->serialize(config->outfbase);

            return 0;
        }
        case Config::STATS: {
            std::ofstream outstream;
            if (!config->outfbase.empty()) {
                outstream.open(config->outfbase + ".dbgstats");
                outstream << "file\tnodes\tedges\tk" << std::endl;
            }
            for (const auto &file : files) {
                std::unique_ptr<DBG_succ> graph {
                    load_critical_graph_from_file(file)
                };

                if (!config->quiet) {
                    std::cout << "Statistics for graph " << file << std::endl;
                    std::cout << "nodes: " << graph->num_nodes() << std::endl;
                    std::cout << "edges: " << graph->num_edges() << std::endl;
                    std::cout << "k: " << graph->get_k() << std::endl;
                    std::cout << "state: " << graph->state << std::endl;
                }
                if (outstream.is_open()) {
                    outstream << file << "\t"
                              << graph->num_nodes() << "\t"
                              << graph->num_edges() << "\t"
                              << graph->get_k() << std::endl;
                }
                if (config->print_graph_succ)
                    graph->print_state();

                //TODO: fix this, options for different annotations
                /*
                std::ifstream instream(file + ".anno.dbg");
                if (instream.good()) {
                    size_t anno_size = NumberSerialisation::deserialiseNumber(instream);
                    std::cout << "annot: " << anno_size << std::endl;
                }
                */
            }

            return 0;
        }
        case Config::TRANSFORM: {
            Timer timer;
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(files.at(0))
            };

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

            std::cout << "Align sequences against a de Bruijn graph with ";
            std::cout << "k=" << graph->get_k() << "\n"
                      << "Length of aligning k-mers: " << config->alignment_length << std::endl;

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

                        for (size_t l = 0;  l < graphindices.at(i).size(); ++l) {
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

                        // try reverse
                        if (graphindices.at(i).size() == 0) {
                            reverse_complement(read_stream->seq);
                            graphindices = graph->align_fuzzy(
                                std::string(read_stream->seq.s, read_stream->seq.l),
                                aln_len, config->distance
                            );
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
                        }
                        std::cout << std::endl;
                        return;
                    }

                    if (config->verbose) {
                        std::cout << "Sequence: " << read_stream->seq.s << "\n";
                    }

                    if (config->query_presence
                            && config->alignment_length == graph->get_k()) {
                        std::cout << graph->find(read_stream->seq.s,
                                                 config->discovery_fraction) << "\n";
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

                        std::cout << (num_discovered >= min_kmers_discovered) << "\n";
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
                });
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
