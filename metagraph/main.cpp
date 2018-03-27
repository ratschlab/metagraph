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
#include "annotate.hpp"
#include "unix_tools.hpp"
#include "kmer.hpp"

using libmaus2::util::NumberSerialisation;

const size_t kMaxNumParallelReadFiles = 5;

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
            Timer timer;
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(files.at(0))
            };

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

                DBG_succ::VectorChunk graph_data;

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

            if (config->state == Config::DYN) {
                std::cout << "Converting static graph to dynamic...\t" << std::flush;
                timer.reset();
                graph->switch_state(Config::DYN);
                std::cout << timer.elapsed() << "sec" << std::endl;
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

            // initialize empty Bloom filter annotation
            hash_annotate::BloomAnnotator *bloom_annotation = NULL;
            annotate::DBGSuccAnnotWrapper *graph_wrapper = NULL;
            if (config->bloom_fpp > -0.5
                    || config->bloom_bits_per_edge > -0.5
                    || config->bloom_num_hash_functions > 0) {
                graph_wrapper = new annotate::DBGSuccAnnotWrapper(*graph);
                if (config->bloom_fpp > -0.5) {
                    // Expected FPP is set, optimize other parameters automatically
                    bloom_annotation = new hash_annotate::BloomAnnotator(
                        *graph_wrapper,
                        config->bloom_fpp,
                        config->verbose
                    );
                } else {
                    assert(config->bloom_bits_per_edge >= 0);
                    // Experiment mode, estimate FPP given other parameters,
                    // optimize the number of hash functions if it's set to zero
                    bloom_annotation = new hash_annotate::BloomAnnotator(
                        *graph_wrapper,
                        config->bloom_bits_per_edge,
                        config->bloom_num_hash_functions,
                        config->verbose
                    );
                }
            }

            // initialize empty annotation
            annotate::ColorCompressed annotation(1 + graph->num_edges());
            std::unordered_map<std::string, size_t> label_store;

            size_t total_seqs = 0;
            bool has_vcf = false;

            // iterate over input files
            for (const auto &file : files) {
                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << file << std::endl;
                }
                // open stream
                if (utils::get_filetype(file) == "VCF") {
                    has_vcf = true;
                    std::vector<std::string> variant_labels;
                    std::unordered_map<std::string, std::set<size_t>> label_map;
                    std::map<size_t, std::string> accum_variants;
                    read_vcf_file_critical(file, config->refpath, graph->get_k(), &variant_labels,
                        [&](std::string &seq, std::vector<std::string> *variant_labels) {
                            assert(variant_labels);
                            for (auto &label : *variant_labels) {
                                auto label_ins = label_store.insert(std::make_pair(label, label_store.size()));
                                if (bloom_annotation) {
                                    auto accum_variants_ins = accum_variants.insert(std::make_pair(
                                        label_ins.first->second,
                                        seq
                                    ));
                                    if (!accum_variants_ins.second) {
                                        accum_variants_ins.first->second += std::string("$") + seq;
                                    }
                                }
                            }
                            variant_labels->clear();
                            if (bloom_annotation && !config->bloom_test_num_kmers) {
                                return;
                            }
                            graph->align(seq, [&](uint64_t i) {
                                if (i > 0) {
                                    for (auto &label : *variant_labels) {
                                        label_map[label].insert(i);
                                    }
                                }
                            });
                            if (config->reverse) {
                                reverse_complement(seq.begin(), seq.end());
                                graph->align(seq, [&](uint64_t i) {
                                    if (i > 0) {
                                        for (auto &label : *variant_labels) {
                                            label_map[label].insert(i);
                                        }
                                    }
                                });
                            }
                        }
                    );
                    if (bloom_annotation) {
                        for (auto &variant : accum_variants) {
                            bloom_annotation->add_sequence(
                                variant.second,
                                variant.first,
                                (variant.second.length() - graph->get_k())
                                    * (static_cast<size_t>(config->reverse) + 1)
                            );
                            if (config->reverse) {
                                reverse_complement(variant.second.begin(), variant.second.end());
                                bloom_annotation->add_sequence(
                                    variant.second,
                                    variant.first,
                                    (variant.second.length() - graph->get_k())
                                        * (static_cast<size_t>(config->reverse) + 1)
                                );
                            }
                        }
                    }
                    if (!bloom_annotation || config->bloom_test_num_kmers) {
                        for (auto &label : label_map) {
                            for (auto &index : label.second) {
                                annotation.add_label(index, label.first);
                            }
                        }
                    }
                } else if (utils::get_filetype(file) == "FASTA"
                            || utils::get_filetype(file) == "FASTQ") {
                    std::string label;
                    std::pair<std::unordered_map<std::string, size_t>::iterator, bool> label_ins;
                    if (!config->fasta_anno) {
                        label = file;
                        label_ins = label_store.insert(std::make_pair(label, label_store.size()));
                    }
                    read_fasta_file_critical(file, [&](kseq_t *read_stream) {
                        if (config->fasta_anno) {
                            label = read_stream->name.s;
                            label_ins = label_store.insert(std::make_pair(label, label_store.size()));
                        }
                        size_t seq_length = (utils::get_filetype(file) == "FASTA") ?
                            (read_stream->seq.l - graph->get_k())
                                * (static_cast<size_t>(config->reverse) + 1) :
                            graph->num_edges();
                        if (bloom_annotation) {
                            bloom_annotation->add_sequence(
                                read_stream->seq.s,
                                label_ins.first->second,
                                seq_length
                            );
                        }
                        if (!bloom_annotation || config->bloom_test_num_kmers) {
                            graph->align(read_stream->seq.s, [&](uint64_t i) {
                                if (i > 0) annotation.add_label(i, label);
                            });
                        }
                        if (config->reverse) {
                            reverse_complement(read_stream->seq);
                            if (bloom_annotation) {
                                bloom_annotation->add_sequence(
                                    read_stream->seq.s,
                                    label_ins.first->second,
                                    seq_length
                                );
                            }
                            if (!bloom_annotation || config->bloom_test_num_kmers) {
                                graph->align(read_stream->seq.s, [&](uint64_t i) {
                                    if (i > 0) annotation.add_label(i, label);
                                });
                            }
                        }
                        total_seqs += 1;
                        if (config->verbose && total_seqs % 10000 == 0) {
                            std::cout << "added labels for " << total_seqs
                                      << " sequences"
                                      << ", last was " << read_stream->name.s
                                      << ", annotated as " << label << std::endl;
                        }
                    });
                } else {
                    std::cerr << "ERROR: Filetype unknown for file "
                              << file << std::endl;
                    exit(1);
                }
            }
            if (bloom_annotation && config->bloom_test_num_kmers) {
                std::cout << "Approximating FPP...\t" << std::flush;
                bloom_annotation->test_fp_all(
                    annotate::PreciseColorCompressedAnnotator(annotation),
                    config->bloom_test_num_kmers, has_vcf
                );
                std::cout << std::endl;
            }
            annotation.serialize(config->infbase + ".anno.dbg");
            if (config->dump_raw_anno) {
                annotation.serialize_uncompressed_rows(config->infbase + ".anno.rawrows.dbg");
            }
            if (bloom_annotation) {
                bloom_annotation->serialize(config->infbase + ".anno.bloom.dbg");
                assert(graph_wrapper);
                delete graph_wrapper;
                delete bloom_annotation;
                bloom_annotation = NULL;
                graph_wrapper = NULL;
            }

            return 0;
        }
        case Config::CLASSIFY: {
            // load graph
            std::unique_ptr<DBG_succ> graph {
                load_critical_graph_from_file(config->infbase)
            };

            // load annotatioun (if file does not exist, empty annotation is created)
            annotate::ColorCompressed annotation(1 + graph->num_edges());

            if (!annotation.load(config->infbase + ".anno.dbg")) {
                std::cerr << "ERROR: can't load annotations from "
                          << config->infbase + ".anno.dbg"
                          << " corrupted" << std::endl;
            }
            const auto labels = annotation.get_label_names();

            // iterate over input files
            for (unsigned int f = 0; f < files.size(); ++f) {
                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << files[f] << std::endl;
                }

                read_fasta_file_critical(files[f], [&](kseq_t *read_stream) {
                    std::cout << read_stream->name.s << ":\t";

                    if (config->reverse)
                        reverse_complement(read_stream->seq);

                    size_t max_kmers_missing =
                        (read_stream->seq.l - graph->get_k() + 1)
                            * (1 - config->discovery_fraction);

                    std::map<std::string, size_t> labels_counter;
                    for (const auto &label : labels) {
                        labels_counter[label] = 0;
                    }
                    graph->align(read_stream->seq.s,
                        [&](uint64_t i) {
                            for (auto it = labels_counter.begin(); it != labels_counter.end();) {
                                if ((i > 0 && annotation.has_label(i, it->first))
                                        || ++(it->second) <= max_kmers_missing) {
                                    ++it;
                                } else {
                                    labels_counter.erase(it++);
                                }
                            }
                        },
                        [&]() { return labels_counter.size() == 0; }
                    );
                    for (auto it = labels_counter.begin(); it != labels_counter.end(); ++it) {
                        std::cout << it->first << ",";
                    }
                    std::cout << std::endl;
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
                graph = merge::build_graph_from_chunks(
                    config->k,
                    std::vector<DBG_succ::Chunk*>(files.size(), NULL),
                    files
                );
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
                    std::cout << "Statistics for file " << file << std::endl;
                    std::cout << "nodes: " << graph->num_nodes() << std::endl;
                    std::cout << "edges: " << graph->num_edges() << std::endl;
                    std::cout << "k: " << graph->get_k() << std::endl;
                }
                if (outstream.is_open()) {
                    outstream << file << "\t"
                              << graph->num_nodes() << "\t"
                              << graph->num_edges() << "\t"
                              << graph->get_k() << std::endl;
                }
                if (config->print_graph_succ)
                    graph->print_state();

                std::ifstream instream(file + ".anno.dbg");
                if (instream.good()) {
                    size_t anno_size = NumberSerialisation::deserialiseNumber(instream);
                    std::cout << "annot: " << anno_size << std::endl;
                }
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
                    std::cout << timer.elapsed() << " sec" << std::endl;
                }
            }
            if (config->sqlfbase.size()) {
                if (config->verbose) {
                    std::cout << "Converting graph to SQL DB...\t" << std::flush;
                }
                timer.reset();
                traverse::toSQL(graph.get(), files, config->sqlfbase);
                if (config->verbose) {
                    std::cout << timer.elapsed() << " sec" << std::endl;
                }
            }

            if (config->state == Config::DYN) {
                if (config->verbose) {
                    std::cout << "Converting static graph to dynamic...\t" << std::flush;
                }
                timer.reset();
                graph->switch_state(Config::DYN);
                if (config->verbose) {
                    std::cout << timer.elapsed() << " sec" << std::endl;
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

            for (const auto &file : files) {
                std::cout << "Opening file for alignment ..." << file << std::endl;

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
                    } else {

                        auto graphindices = graph->index(read_stream->seq.s,
                                                         config->alignment_length);

                        if (config->query) {
                            size_t num_discovered = std::count_if(
                                graphindices.begin(), graphindices.end(),
                                [](const auto &x) { return x > 0; }
                            );
                            std::cout << num_discovered
                                      << " / "
                                      << read_stream->seq.l << std::endl;
                        } else {
                            for (size_t i = 0; i < graphindices.size(); ++i) {
                                for (uint64_t j = 0; j < graph->get_k(); ++j) {
                                    std::cout << read_stream->seq.s[i + j];
                                }
                                std::cout << ": " << graphindices[i] << std::endl;
                            }
                        }
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
