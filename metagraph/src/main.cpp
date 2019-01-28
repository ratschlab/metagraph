#include "unix_tools.hpp"
#include "config.hpp"
#include "sequence_io.hpp"
#include "reads_filtering.hpp"
#include "dbg_succinct_construct.hpp"
#include "dbg_succinct_chunk.hpp"
#include "dbg_succinct_merge.hpp"
#include "annotated_dbg.hpp"
#include "annotate_row_compressed.hpp"
#include "annotate_column_compressed.hpp"
#include "serialization.hpp"
#include "utils.hpp"
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"
#include "kmc_parser.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_bitmap.hpp"
#include "dbg_bitmap_construct.hpp"
#include "server.hpp"

typedef annotate::MultiLabelAnnotation<uint64_t, std::string> Annotator;

const size_t kMaxNumParallelReadFiles = 5;

const size_t kNumCachedColumns = 10;


Config::GraphType parse_graph_extension(const std::string &filename) {
    if (utils::ends_with(filename, ".dbg")) {
        return Config::GraphType::SUCCINCT;

    } else if (utils::ends_with(filename, ".orhashdbg")) {
        return Config::GraphType::HASH;

    } else if (utils::ends_with(filename, ".bitmapdbg")) {
        return Config::GraphType::BITMAP;

    } else {
        return Config::GraphType::INVALID;
    }
}

std::string remove_graph_extension(const std::string &filename) {
    return utils::remove_suffix(filename, ".dbg", ".orhashdbg", ".bitmapdbg");
}

template <class Graph = DBG_succ>
std::unique_ptr<Graph> load_critical_graph_from_file(const std::string &filename) {
    auto *graph = new Graph(2);
    if (!graph->load(filename)) {
        std::cerr << "ERROR: can't load graph from file " << filename << std::endl;
        delete graph;
        exit(1);
    }
    return std::unique_ptr<Graph> { graph };
}

template <class DefaultGraphType = DBGSuccinct>
std::unique_ptr<DeBruijnGraph> load_critical_dbg(const std::string &filename) {
    auto graph_type = parse_graph_extension(filename);
    switch (graph_type) {
        case Config::GraphType::SUCCINCT:
            return load_critical_graph_from_file<DBGSuccinct>(filename);

        case Config::GraphType::HASH:
            return load_critical_graph_from_file<DBGHashOrdered>(filename);

        case Config::GraphType::BITMAP:
            return load_critical_graph_from_file<DBGSD>(filename);

        case Config::GraphType::INVALID:
            return load_critical_graph_from_file<DefaultGraphType>(filename);
    }
    assert(false);
    exit(1);
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
                   bool use_kmc,
                   size_t filter_k,
                   size_t max_unreliable_abundance,
                   size_t unreliable_kmers_threshold,
                   bool filename_anno,
                   bool fasta_anno,
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
            std::vector<std::string> variant_labels;

            read_vcf_file_critical(
                file,
                ref_sequence_path,
                dynamic_cast<const DeBruijnGraph &>(anno_graph->get_graph()).get_k(),
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
        } else if (use_kmc) {
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
                max_unreliable_abundance + 1
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
                        std::cout << ", " << timer.elapsed() << "sec" << std::endl;
                    }
                },
                reverse,
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

        if (verbose) {
            std::cout << "File processed in "
                      << data_reading_timer.elapsed()
                      << "sec, current mem usage: "
                      << get_curr_mem2() / (1 << 20) << " MiB"
                      << ", total time: " << timer.elapsed()
                      << "sec" << std::endl;
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
                    if (reverse)
                        forward_strand = !forward_strand;
                },
                reverse,
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

        if (verbose) {
            std::cout << "File processed in "
                      << data_reading_timer.elapsed()
                      << "sec, current mem usage: "
                      << get_curr_mem2() / (1 << 20) << " MiB"
                      << ", total time: " << timer.elapsed()
                      << "sec" << std::endl;
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
                   const AnnotatedDBG &anno_graph,
                   std::ostream &output_stream) {
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

    output_stream << oss.str();
}


std::unique_ptr<Annotator> initialize_annotation(const std::string &filename,
                                                 const Config &config,
                                                 uint64_t num_rows = 0) {
    std::unique_ptr<Annotator> annotation;

    Config::AnnotationType anno_type = config.anno_type;

    if (utils::ends_with(filename, annotate::kColumnAnnotatorExtension)) {
        anno_type = Config::AnnotationType::ColumnCompressed;

    } else if (utils::ends_with(filename, annotate::kRowAnnotatorExtension)) {
        anno_type = Config::AnnotationType::RowCompressed;

    } else if (utils::ends_with(filename, annotate::kBRWTExtension)) {
        anno_type = Config::AnnotationType::BRWT;

    } else if (utils::ends_with(filename, annotate::kBinRelWT_sdslExtension)) {
        anno_type = Config::AnnotationType::BinRelWT_sdsl;

    } else if (utils::ends_with(filename, annotate::kBinRelWTExtension)) {
        anno_type = Config::AnnotationType::BinRelWT;

    } else if (utils::ends_with(filename, annotate::kRowPackedExtension)) {
        anno_type = Config::AnnotationType::RowFlat;

    } else if (utils::ends_with(filename, annotate::kRainbowfishExtension)) {
        anno_type = Config::AnnotationType::RBFish;
    }

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
            annotation.reset(new annotate::BRWTCompressed<>());
            break;
        }
        case Config::BinRelWT_sdsl: {
            annotation.reset(new annotate::BinRelWT_sdslAnnotator());
            break;
        }
        case Config::BinRelWT: {
            annotation.reset(new annotate::BinRelWTAnnotator());
            break;
        }
        case Config::RowFlat: {
            annotation.reset(new annotate::RowFlatAnnotator());
            break;
        }
        case Config::RBFish: {
            annotation.reset(new annotate::RainbowfishAnnotator());
            break;
        }
    }

    return annotation;
}


std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(const Config &config) {
    auto graph_temp = load_critical_dbg(config.infbase);

    auto annotation_temp = initialize_annotation(
        config.infbase_annotators.size()
            ? config.infbase_annotators.at(0)
            : "",
        config,
        graph_temp->num_nodes()
    );
    if (config.infbase_annotators.size()
            && !annotation_temp->load(config.infbase_annotators.at(0))) {
        std::cerr << "ERROR: can't load annotations for graph "
                  << config.infbase
                  << ", file corrupted" << std::endl;
        exit(1);
    }

    // load graph
    auto anno_graph = std::make_unique<AnnotatedDBG>(graph_temp.release(),
                                                     annotation_temp.release(),
                                                     config.parallel);

    if (!anno_graph->check_compatibility()) {
        std::cerr << "Error: graph and annotation are not compatible."
                  << std::endl;
        exit(1);
    }

    return anno_graph;
}


void print_boss_stats(const DBG_succ &boss_graph,
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
    std::cout << "========================================================" << std::endl;
}

void print_stats(const Annotator &annotation) {
    std::cout << "=================== ANNOTATION STATS ===================" << std::endl;
    std::cout << "labels:  " << annotation.num_labels() << std::endl;
    std::cout << "objects: " << annotation.num_objects() << std::endl;
    std::cout << "density: " << std::scientific
                              << static_cast<double>(annotation.num_relations())
                                    / annotation.num_objects()
                                    / annotation.num_labels() << std::endl;
    std::cout << "representation: ";

    if (dynamic_cast<const annotate::ColumnCompressed<std::string> *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::ColumnCompressed) << std::endl;

    } else if (dynamic_cast<const annotate::RowCompressed<std::string> *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::RowCompressed) << std::endl;

    } else if (dynamic_cast<const annotate::BRWTCompressed<std::string> *>(&annotation)) {
        std::cout << Config::annotype_to_string(Config::BRWT) << std::endl;

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


template <class Callback, class Loop>
void parse_sequences(const std::vector<std::string> &files,
                     const Config &config,
                     const Timer &timer,
                     Callback call_read,
                     Loop call_reads) {
    // iterate over input files
    for (unsigned int f = 0; f < files.size(); ++f) {
        if (config.verbose) {
            std::cout << std::endl << "Parsing " << files[f] << std::endl;
        }

        Timer data_reading_timer;

        if (utils::get_filetype(files[f]) == "VCF") {
            //assume VCF contains no noise
            read_vcf_file_critical(
                files[f], config.refpath, config.k, NULL,
                [&](std::string &seq, auto *) {
                    call_read(seq);
                    if (config.reverse) {
                        reverse_complement(seq.begin(), seq.end());
                        call_read(seq);
                    }
                }
            );
        } else if (config.use_kmc) {
            bool warning_different_k = false;
            kmc::read_kmers(
                files[f],
                [&](std::string&& sequence) {
                    if (!warning_different_k && sequence.size() != config.k) {
                        std::cerr << "Warning: k-mers parsed from KMC database "
                                  << files[f] << " have length " << sequence.size()
                                  << " but graph is constructed for k=" << config.k
                                  << std::endl;
                        warning_different_k = true;
                    }
                    call_read(sequence);
                },
                config.max_unreliable_abundance + 1
            );
        } else if (utils::get_filetype(files[f]) == "FASTA"
                    || utils::get_filetype(files[f]) == "FASTQ") {
            std::string filter_filename = get_filter_filename(
                files[f], config.filter_k,
                config.max_unreliable_abundance,
                config.unreliable_kmers_threshold
            );

            if (files.size() >= config.parallel) {
                auto reverse = config.reverse;
                auto file = files[f];

                // capture all required values by copying to be able
                // to run task from other threads
                call_reads([=](auto callback) {
                    read_fasta_file_critical(file, [=](kseq_t *read_stream) {
                        // add read to the graph constructor as a callback
                        callback(read_stream->seq.s);
                    }, reverse, filter_filename);
                });
            } else {
                read_fasta_file_critical(files[f], [&](kseq_t *read_stream) {
                    // add read to the graph constructor as a callback
                    call_read(read_stream->seq.s);
                }, config.reverse, filter_filename);
            }
        } else {
            std::cerr << "ERROR: Filetype unknown for file "
                      << files[f] << std::endl;
            exit(1);
        }

        if (config.verbose) {
            std::cout << "Finished extracting sequences from file " << files[f]
                      << " in " << timer.elapsed() << "sec" << std::endl;
        }
        if (config.verbose) {
            std::cout << "File processed in "
                      << data_reading_timer.elapsed()
                      << "sec, current mem usage: "
                      << get_curr_mem2() / (1 << 20) << " MiB"
                      << ", total time: " << timer.elapsed()
                      << "sec" << std::endl;
        }
    }

    if (config.verbose) {
        std::cout << std::endl;
    }
}


int main(int argc, const char *argv[]) {
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
            std::unique_ptr<DeBruijnGraph> graph;

            if (config->verbose)
                std::cout << "Build De Bruijn Graph with k-mer size k="
                          << config->k << std::endl;

            Timer timer;

            if (config->graph_type == Config::GraphType::SUCCINCT && !config->dynamic) {
                if (config->canonical)
                    config->reverse = true;

                std::unique_ptr<DBG_succ> boss_graph { new DBG_succ(config->k - 1) };

                if (config->verbose) {
                    std::cout << "Start reading data and extracting k-mers" << std::endl;
                }
                //enumerate all suffices
                assert(boss_graph->alph_size > 1);
                assert(config->nsplits > 0);
                size_t suffix_len = std::min(
                    static_cast<size_t>(std::ceil(std::log2(config->nsplits)
                                                    / std::log2(boss_graph->alph_size - 1))),
                    boss_graph->get_k() - 1
                );
                std::deque<std::string> suffices;
                if (config->suffix.size()) {
                    suffices = { config->suffix };
                } else {
                    suffices = utils::generate_strings(
                        boss_graph->alphabet,
                        suffix_len
                    );
                }

                DBG_succ::Chunk graph_data(boss_graph->get_k());

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

                    std::unique_ptr<IDBGBOSSChunkConstructor> constructor(
                        IDBGBOSSChunkConstructor::initialize(
                            boss_graph->get_k(),
                            suffix,
                            config->parallel,
                            static_cast<uint64_t>(config->memory_available) << 30,
                            config->verbose
                        )
                    );

                    parse_sequences(files, *config, timer,
                        [&](const auto &read) { constructor->add_sequence(read); },
                        [&](const auto &loop) { constructor->add_sequences(loop); }
                    );

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

                    if (config->suffix.size())
                        return 0;
                }

                graph_data.initialize_graph(boss_graph.get());
                graph.reset(new DBGSuccinct(boss_graph.release(), config->canonical));

            } else if (config->graph_type == Config::GraphType::BITMAP && !config->dynamic) {
                std::unique_ptr<DBGSD> sd_graph { new DBGSD(config->k) };

                if (config->verbose) {
                    std::cout << "Start reading data and extracting k-mers" << std::endl;
                }
                //enumerate all suffices
                assert(sd_graph->alphabet.size() > 1);
                assert(config->nsplits > 0);
                size_t suffix_len = std::min(
                    static_cast<size_t>(std::ceil(std::log2(config->nsplits)
                                                    / std::log2(sd_graph->alphabet.size() - 1))),
                    sd_graph->get_k() - 2
                );
                std::deque<std::string> suffices;
                if (config->suffix.size()) {
                    suffices = { config->suffix };
                } else {
                    suffices = utils::generate_strings(
                        sd_graph->alphabet,
                        suffix_len
                    );
                }

                std::unique_ptr<DBGSDConstructor> constructor;

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

                    constructor.reset(
                        new DBGSDConstructor(
                            sd_graph->get_k(),
                            config->canonical,
                            suffix,
                            config->parallel,
                            static_cast<uint64_t>(config->memory_available) << 30,
                            config->verbose
                        )
                    );

                    parse_sequences(files, *config, timer,
                        [&](const auto &read) { constructor->add_sequence(read); },
                        [&](const auto &loop) { constructor->add_sequences(loop); }
                    );

                    if (config->outfbase.size() && config->suffix.size()) {
                        const DBGSD::Chunk* next_block = constructor->build_chunk();
                        std::cout << "Graph chunk with " << next_block->num_set_bits()
                                  << " k-mers was built in "
                                  << timer.elapsed() << "sec" << std::endl;

                        std::cout << "Serialize the graph chunk for suffix '"
                                  << suffix << "'...\t" << std::flush;
                        timer.reset();
                        std::ofstream out(config->outfbase + "." + suffix + ".dbgsdchunk");
                        next_block->serialize(out);
                        std::cout << timer.elapsed() << "sec" << std::endl;
                        sd_graph.reset(new DBGSD(config->k));
                    }

                    if (config->suffix.size())
                        return 0;

                    constructor->build_graph(sd_graph.get());
                }

                graph.reset(sd_graph.release());

            } else {
                if (config->canonical)
                    config->reverse = false;

                //slower method

                switch (config->graph_type) {
                    case Config::GraphType::SUCCINCT:
                        graph.reset(new DBGSuccinct(config->k, config->canonical));
                        break;
                    case Config::GraphType::HASH:
                        graph.reset(new DBGHashOrdered(config->k, config->canonical));
                        break;
                    case Config::GraphType::BITMAP:
                        assert(false);
                    case Config::GraphType::INVALID:
                        assert(false);
                }
                assert(graph.get());

                parse_sequences(files, *config, timer,
                    [&graph](const auto &seq) { graph->add_sequence(seq); },
                    [&graph](const auto &loop) {
                        loop([&graph](const auto &seq) { graph->add_sequence(seq); });
                    }
                );
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
            }

            return 0;
        }
        case Config::EXTEND: {
            assert(config->infbase_annotators.size() <= 1);

            Timer timer;

            // load graph
            auto graph = load_critical_dbg(config->infbase);

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
            if (config->infbase_annotators.size())
                inserted_edges.reset(new bit_vector_dyn(graph->num_nodes() + 1, 0));

            timer.reset();

            if (config->verbose)
                std::cout << "Start graph extension" << std::endl;

            // Insert new k-mers
            parse_sequences(files, *config, timer,
                [&graph,&inserted_edges](const auto &seq) {
                    graph->add_sequence(seq, inserted_edges.get());
                },
                [&graph,&inserted_edges](const auto &loop) {
                    loop([&graph,&inserted_edges](const auto &seq) {
                        graph->add_sequence(seq, inserted_edges.get());
                    });
                }
            );
            // if (config->verbose) {
            //     std::cout << "Number of k-mers in graph: " << graph->num_nodes() << std::endl;
            // }

            if (config->verbose)
                std::cout << "Graph extension finished in "
                          << timer.elapsed() << "sec" << std::endl;

            assert(config->outfbase.size());

            // serialize graph
            timer.reset();

            graph->serialize(config->outfbase);
            graph.reset();

            if (config->verbose)
                std::cout << "Serialized in " << timer.elapsed() << "sec" << std::endl;

            timer.reset();

            auto annotation = initialize_annotation(config->infbase_annotators.size()
                                                        ? config->infbase_annotators.at(0)
                                                        : "",
                                                    *config);

            if (!annotation->load(config->infbase_annotators.at(0))) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            } else if (config->verbose) {
                std::cout << "Annotation was loaded in "
                          << timer.elapsed() << "sec" << std::endl;
            }

            timer.reset();

            assert(inserted_edges.get());

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

                Timer data_reading_timer;

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
                            false, filter_filename
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

                    if (config->verbose) {
                        std::cout << "File processed in "
                                  << data_reading_timer.elapsed()
                                  << "sec, current mem usage: "
                                  << get_curr_mem2() / (1 << 20) << " MiB"
                                  << ", total time: " << timer.elapsed()
                                  << "sec" << std::endl;
                    }

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
                                }, reverse);
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
                        std::ofstream outstream(
                            get_filter_filename(
                                file, k, max_unreliable_abundance, unreliable_kmers_threshold, false
                            ),
                            std::ios::binary
                        );
                        serialize_number_vector(outstream, filter, 1);
                    },
                    config->filter_k,
                    config->max_unreliable_abundance,
                    config->unreliable_kmers_threshold,
                    config->verbose, config->reverse, config->use_kmc
                );

                if (config->verbose) {
                    std::cout << "File processed in "
                              << data_reading_timer.elapsed()
                              << "sec, current mem usage: "
                              << get_curr_mem2() / (1 << 20) << " MiB"
                              << ", total time: " << timer.elapsed()
                              << "sec" << std::endl;
                }
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

                std::ifstream instream(filter_filename, std::ios::binary);
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
            assert(config->infbase_annotators.size() <= 1);

            auto anno_graph = initialize_annotated_dbg(*config);

            annotate_data(files,
                          config->refpath,
                          anno_graph.get(),
                          config->reverse,
                          config->use_kmc,
                          config->filter_k,
                          config->max_unreliable_abundance,
                          config->unreliable_kmers_threshold,
                          config->filename_anno,
                          config->fasta_anno,
                          config->fasta_header_delimiter,
                          config->anno_labels,
                          config->verbose);

            anno_graph->get_annotation().serialize(config->outfbase);

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
            AnnotatedDBG anno_graph(
                graph_temp.release(),
                annotation_temp.release(),
                config->parallel
            );

            if (!anno_graph.check_compatibility()) {
                std::cerr << "Error: graph and annotation are not compatible."
                          << std::endl;
                exit(1);
            }

            annotate_coordinates(files,
                                 &anno_graph,
                                 config->reverse,
                                 config->filter_k,
                                 config->max_unreliable_abundance,
                                 config->unreliable_kmers_threshold,
                                 config->genome_binsize_anno,
                                 config->verbose);

            anno_graph.get_annotation().serialize(config->outfbase);

            return 0;
        }
        case Config::MERGE_ANNOTATIONS: {
            std::unique_ptr<Annotator> annotation;
            if (config->anno_type == Config::ColumnCompressed) {
                annotation.reset(
                    new annotate::ColumnCompressed<>(
                        0, kNumCachedColumns, config->verbose
                    )
                );
            } else {
                throw std::runtime_error("To be implemented");
            }

            if (!annotation->merge_load(files)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }

            annotation->serialize(config->outfbase);

            return 0;
        }
        case Config::CLASSIFY: {
            assert(config->infbase_annotators.size() == 1);

            auto anno_graph = initialize_annotated_dbg(*config);

            utils::ThreadPool thread_pool(std::max(1u, config->parallel) - 1);

            Timer timer;

            // iterate over input files
            for (const auto &file : files) {
                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << file << std::endl;
                }

                Timer data_reading_timer;

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
                            std::ref(*anno_graph),
                            std::ref(std::cout)
                        );
                    },
                    config->reverse,
                    get_filter_filename(file, config->filter_k,
                                        config->max_unreliable_abundance,
                                        config->unreliable_kmers_threshold)
                );
                if (config->verbose) {
                    std::cout << "Finished extracting sequences from file " << file
                              << " in " << timer.elapsed() << "sec" << std::endl;
                }
                if (config->verbose) {
                    std::cout << "File processed in "
                              << data_reading_timer.elapsed()
                              << "sec, current mem usage: "
                              << get_curr_mem2() / (1 << 20) << " MiB"
                              << ", total time: " << timer.elapsed()
                              << "sec" << std::endl;
                }

                // wait while all threads finish processing the current file
                thread_pool.join();
            }

            return 0;
        }
        case Config::SERVER_CLASSIFY: {
            assert(config->infbase_annotators.size() == 1);

            Timer timer;

            std::cout << "Loading graph..." << std::endl;

            auto anno_graph = initialize_annotated_dbg(*config);

            std::cout << "Graph loaded in "
                      << timer.elapsed() << "sec, current mem usage: "
                      << get_curr_mem2() / (1 << 20) << " MiB" << std::endl;

            const size_t num_threads = std::max(1u, config->parallel);

            std::cout << "Initializing tcp service with "
                      << num_threads << " threads, listening port "
                      << config->port << std::endl;

            try {
                asio::io_context io_context(num_threads);

                asio::signal_set signals(io_context, SIGINT, SIGTERM);
                signals.async_wait([&](auto, auto) { io_context.stop(); });

                Server server(io_context, config->port,
                    [&](const std::string &input_sequence) {
                        std::ostringstream oss;
                        execute_query(input_sequence,
                                      input_sequence,
                                      config->count_labels,
                                      config->suppress_unlabeled,
                                      config->num_top_labels,
                                      config->discovery_fraction,
                                      config->anno_labels_delimiter,
                                      *anno_graph,
                                      oss);
                        return oss.str();
                    }
                );

                io_context.run();
            } catch (const std::exception &e) {
                std::cerr << "Exception: " << e.what() << std::endl;
            } catch (...) {
                std::cerr << "Error: Unknown exception" << std::endl;
            }

            return 0;
        }
        case Config::COMPARE: {
            assert(files.size());

            std::cout << "Opening file                " << files.at(0) << std::endl;
            auto graph = load_critical_graph_from_file(files.at(0));

            for (size_t f = 1; f < files.size(); ++f) {
                std::cout << "Opening file for comparison " << files[f] << std::endl;
                auto second = load_critical_graph_from_file(files[f]);
                if (config->internal
                        ? graph->equals_internally(*second, config->verbose)
                        : *graph == *second) {
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

            if (!files.size()) {
                assert(config->infbase.size());

                auto sorted_suffices = utils::generate_strings(
                    config->graph_type == Config::GraphType::SUCCINCT
                        ? KmerExtractor().alphabet
                        : KmerExtractor2Bit().alphabet,
                    config->suffix_len
                );

                for (const std::string &suffix : sorted_suffices) {
                    assert(suffix.size() == config->suffix_len);

                    if (valid_kmer_suffix(suffix))
                        chunk_files.push_back(config->infbase + "." + suffix);
                }
            }

            for (auto &filename : chunk_files) {
                filename = utils::remove_suffix(filename, ".dbgchunk", ".dbgsdchunk");
            }

            // collect results on an external merge or construction
            std::unique_ptr<DeBruijnGraph> graph;
            if (config->graph_type == Config::GraphType::SUCCINCT) {
                graph.reset(
                    new DBGSuccinct(DBG_succ::Chunk::build_graph_from_chunks(
                        chunk_files, config->verbose
                    ), config->canonical)
                );
            } else if (config->graph_type == Config::GraphType::BITMAP) {
                graph.reset(
                    DBGSDConstructor::build_graph_from_chunks(
                        chunk_files, config->canonical, config->verbose
                    )
                );
            }
            assert(graph.get());

            if (config->verbose) {
                std::cout << "Graph has been assembled" << std::endl;
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
            DBG_succ *graph = NULL;

            Timer timer;

            std::vector<std::unique_ptr<DBGSuccinct>> dbg_graphs;
            std::vector<const DBG_succ*> graphs;

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
                    graph = new DBG_succ(graphs[0]->get_k());
                    chunk->initialize_graph(graph);
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
        case Config::STATS: {
            for (const auto &file : files) {
                auto graph = load_critical_dbg(file);

                std::cout << "Statistics for graph " << file << std::endl;

                print_stats(*graph);

                if (!dynamic_cast<DBGSuccinct*>(graph.get()))
                    continue;

                const auto &boss_graph = dynamic_cast<DBGSuccinct&>(*graph).get_boss();

                print_boss_stats(boss_graph,
                                 config->count_dummy,
                                 config->parallel,
                                 config->verbose);

                if (config->print_graph_internal_repr)
                    boss_graph.print_internal_representation(std::cout);

                if (config->print_graph)
                    std::cout << boss_graph;
            }

            for (const auto &file : config->infbase_annotators) {
                auto annotation = initialize_annotation(file, *config);
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

            auto annotator = std::make_unique<annotate::ColumnCompressed<>>(
                0, kNumCachedColumns, config->verbose
            );

            if (config->verbose)
                std::cout << "Loading annotator...\t" << std::flush;

            if (!annotator->load(files.at(0))) {
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

            switch (config->anno_type) {
                case Config::ColumnCompressed:
                    break;
                case Config::RowCompressed: {
                    if (config->verbose)
                        std::cout << "Converting...\t" << std::flush;

                    annotate::RowCompressed<> row_annotator(0);
                    annotator->convert_to_row_annotator(&row_annotator,
                                                        config->parallel);
                    annotator.reset();

                    row_annotator.serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
                case Config::BRWT: {
                    if (config->verbose)
                        std::cout << "Converting...\t" << std::flush;

                    auto brwt_annotator = config->greedy_brwt
                        ? annotate::convert_to_greedy_BRWT<annotate::BRWTCompressed<>>(
                            std::move(*annotator),
                            config->parallel)
                        : annotate::convert_to_simple_BRWT<annotate::BRWTCompressed<>>(
                            std::move(*annotator),
                            config->arity_brwt,
                            config->parallel);

                    annotator.reset();

                    brwt_annotator->serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
                case Config::BinRelWT_sdsl: {
                    if (config->verbose)
                        std::cout << "Converting...\t" << std::flush;

                    auto binrelwt_sdsl_annotator
                            = annotate::convert<annotate::BinRelWT_sdslAnnotator>(
                        std::move(*annotator)
                    );
                    annotator.reset();
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;

                    if (config->verbose)
                        std::cout << "Serializing to " << config->outfbase
                                  << "...\t" << std::flush;
                    binrelwt_sdsl_annotator->serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
                case Config::BinRelWT: {
                    if (config->verbose)
                        std::cout << "Converting...\t" << std::flush;

                    auto binrelwt_annotator = annotate::convert<annotate::BinRelWTAnnotator>(
                        std::move(*annotator)
                    );
                    annotator.reset();
                    if (config->verbose)
                       std::cout << timer.elapsed() << "sec" << std::endl;

                    if (config->verbose)
                        std::cout << "Serializing to " << config->outfbase
                                  << "...\t" << std::flush;
                    binrelwt_annotator->serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
                case Config::RowFlat: {
                    if (config->verbose)
                        std::cout << "Converting to flat annotator...\t" << std::flush;

                    auto flat_annotator = annotate::convert<annotate::RowFlatAnnotator>(
                        std::move(*annotator)
                    );
                    annotator.reset();
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;

                    if (config->verbose)
                        std::cout << "Serializing to " << config->outfbase
                                  << "...\t" << std::flush;
                    flat_annotator->serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
                case Config::RBFish: {
                    if (config->verbose)
                        std::cout << "Converting to rainbowfish annotator...\t" << std::flush;

                    auto flat_annotator = annotate::convert<annotate::RainbowfishAnnotator>(
                        std::move(*annotator)
                    );
                    annotator.reset();

                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;

                    if (config->verbose)
                        std::cout << "Serializing to " << config->outfbase
                                  << "...\t" << std::flush;
                    flat_annotator->serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
            }
            return 0;
        }
        case Config::TRANSFORM: {
            assert(files.size() == 1);

            Timer timer;
            if (config->verbose) {
                std::cout << "Graph loading...\t" << std::flush;
            }
            auto dbg = load_critical_graph_from_file<DBGSuccinct>(files.at(0));
            auto *graph = &dbg->get_boss();

            if (config->verbose) {
                std::cout << timer.elapsed() << "sec" << std::endl;
            }

            if (config->clear_dummy) {
                if (config->verbose) {
                    std::cout << "Traverse source dummy edges and remove redundant ones..." << std::endl;
                }
                timer.reset();

                // remove redundant dummy edges and mark all other dummy edges
                dbg->mask_dummy_kmers(config->parallel, true);

                if (config->verbose) {
                    std::cout << "Done in " << timer.elapsed() << "sec" << std::endl;
                }
                timer.reset();
            }

            if (config->to_fasta) {
                if (config->verbose) {
                    std::cout << "Extracting sequences from graph...\t" << std::flush;
                }
                timer.reset();
                if (!config->outfbase.size()) {
                    std::cerr << "Error: no output file provided" << std::endl;
                    exit(1);
                }

                auto out_filename
                    = utils::remove_suffix(config->outfbase, ".gz", ".fasta")
                        + ".fasta.gz";

                gzFile out_fasta_gz = gzopen(out_filename.c_str(), "w");

                if (out_fasta_gz == Z_NULL) {
                    std::cerr << "ERROR: Can't write to " << out_filename << std::endl;
                    exit(1);
                }

                if (config->contigs || config->pruned_dead_end_size > 0) {
                    graph->call_contigs([&](const auto &sequence) {
                        if (!write_fasta(out_fasta_gz, "", sequence)) {
                            std::cerr << "ERROR: Can't write extracted sequences to "
                                      << out_filename << std::endl;
                            exit(1);
                        }
                    }, config->pruned_dead_end_size);
                } else {
                    graph->call_sequences([&](const auto &sequence) {
                        if (!write_fasta(out_fasta_gz, "", sequence)) {
                            std::cerr << "ERROR: Can't write extracted sequences to "
                                      << out_filename << std::endl;
                            exit(1);
                        }
                    });
                }

                gzclose(out_fasta_gz);

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

            dbg->switch_state(config->state);

            if (config->verbose) {
                std::cout << timer.elapsed() << "sec" << std::endl;
            }

            if (config->outfbase.size()) {
                if (config->verbose) {
                    std::cout << "Serializing transformed graph...\t" << std::flush;
                    timer.reset();
                }
                dbg->serialize(config->outfbase);
                if (config->verbose) {
                    std::cout << timer.elapsed() << "sec" << std::endl;
                }
            }

            return 0;
        }
        case Config::RELAX_BRWT: {
            assert(files.size() == 1);

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

            // load graph
            auto graph = load_critical_graph_from_file(config->infbase);

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

            Timer timer;

            std::cout << "Align sequences against the de Bruijn graph with ";
            std::cout << "k=" << graph->get_k() << "\n"
                      << "Length of aligning k-mers: "
                      << config->alignment_length << std::endl;

            for (const auto &file : files) {
                std::cout << "Align sequences from file " << file << std::endl;

                Timer data_reading_timer;

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
                }, config->reverse,
                    get_filter_filename(file, config->filter_k,
                                        config->max_unreliable_abundance,
                                        config->unreliable_kmers_threshold)
                );

                if (config->verbose) {
                    std::cout << "File processed in "
                              << data_reading_timer.elapsed()
                              << "sec, current mem usage: "
                              << get_curr_mem2() / (1 << 20) << " MiB"
                              << ", total time: " << timer.elapsed()
                              << "sec" << std::endl;
                }
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
