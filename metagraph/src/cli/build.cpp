#include "build.hpp"

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/file_utils.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/bitmap/dbg_bitmap_construct.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "graph/graph_extensions/node_weights.hpp"
#include "config/config.hpp"
#include "parse_sequences.hpp"
#include "stats.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;

using mtg::common::logger;
using mtg::common::get_verbose;

const uint64_t kBytesInGigabyte = 1'000'000'000;


template <class GraphConstructor>
void push_sequences(const std::vector<std::string> &files,
                    const Config &config,
                    const Timer &timer,
                    GraphConstructor *constructor) {
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic, 1)
    for (size_t i = 0; i < files.size(); ++i) {
        BatchAccumulator<std::pair<std::string, uint64_t>> batcher(
            [constructor](auto&& sequences) {
                constructor->add_sequences(std::move(sequences));
            },
            1'000'000 / sizeof(std::pair<std::string, uint64_t>),
            1'000'000
        );
        parse_sequences(files[i], config,
            [&](std::string_view seq) {
                batcher.push_and_pay(seq.size(), seq, 1);
            },
            [&](std::string_view seq, uint32_t count) {
                batcher.push_and_pay(seq.size(), seq, count);
            }
        );
        logger->trace("Extracted all sequences from file {} in {} sec",
                      files[i], timer.elapsed());
    }
}

int build_graph(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    std::unique_ptr<DeBruijnGraph> graph;

    logger->trace("Build De Bruijn Graph with k-mer size k={}", config->k);

    if (!config->outfbase.size()) {
        logger->error("No output file provided");
        exit(1);
    }

    Timer timer;

    if (config->forward_and_reverse) {
        switch (config->graph_mode) {
            case DeBruijnGraph::BASIC:
                logger->warn("Building a graph from both strands is probably undesired."
                              " To represent both strands, consider building in canonical"
                              " or primary mode.");
                break;
            case DeBruijnGraph::CANONICAL:
                config->forward_and_reverse = false;
                break;
            case DeBruijnGraph::PRIMARY:
                logger->error("Passed '--fwd-and-reverse', but primary graphs"
                              " by definition can't contain both strands.");
                exit(1);
        }
    }

    if (config->complete) {
        if (config->graph_type != Config::GraphType::BITMAP) {
            logger->error("Only bitmap-graph can be built in complete mode");
            exit(1);
        }

        graph.reset(new DBGBitmap(config->k, config->graph_mode));

    } else if (config->graph_type == Config::GraphType::SUCCINCT && !config->dynamic) {
        auto boss_graph = std::make_unique<boss::BOSS>(config->k - 1);

        logger->trace("Start reading data and extracting k-mers");
        //enumerate all suffixes
        assert(boss_graph->alph_size > 1);
        std::vector<std::string> suffixes;
        if (config->suffix.size()) {
            suffixes = { config->suffix };
        } else {
            suffixes = kmer::KmerExtractorBOSS::generate_suffixes(config->suffix_len);
        }

        boss::BOSS::Chunk graph_data;

        //one pass per suffix
        for (const std::string &suffix : suffixes) {
            timer.reset();

            if (suffix.size() > 0 || suffixes.size() > 1) {
                logger->info("k-mer suffix: '{}'", suffix);
            }

            auto constructor = boss::IBOSSChunkConstructor::initialize(
                boss_graph->get_k(),
                config->graph_mode == DeBruijnGraph::CANONICAL,
                config->count_width,
                suffix,
                get_num_threads(),
                config->memory_available * kBytesInGigabyte,
                config->tmp_dir.empty() ? kmer::ContainerType::VECTOR
                                        : kmer::ContainerType::VECTOR_DISK,
                config->tmp_dir.empty() ? std::filesystem::path(config->outfbase).remove_filename()
                                        : config->tmp_dir,
                config->disk_cap_bytes
            );

            push_sequences(files, *config, timer, constructor.get());

            boss::BOSS::Chunk next_chunk = constructor->build_chunk();
            logger->trace("Graph chunk with {} k-mers was built in {} sec",
                          next_chunk.size() - 1, timer.elapsed());

            if (config->suffix.size()) {
                logger->info("Serialize the graph chunk for suffix '{}'...", suffix);
                timer.reset();
                next_chunk.serialize(config->outfbase + "." + suffix);
                logger->info("Serialization done in {} sec", timer.elapsed());
                return 0;
            }

            if (graph_data.size()) {
                graph_data.extend(next_chunk);
            } else {
                graph_data = std::move(next_chunk);
            }
        }

        assert(graph_data.size());

        if (!config->mark_dummy_kmers && !config->node_suffix_length) {
            DBGSuccinct::serialize(std::move(graph_data), config->outfbase, config->graph_mode);
            logger->trace("Graph construction finished in {} sec", timer.elapsed());
            return 0;

        } else if (config->count_kmers) {
            sdsl::int_vector_buffer<> kmer_counts;
            graph_data.initialize_boss(boss_graph.get(), &kmer_counts);
            graph.reset(new DBGSuccinct(boss_graph.release(), config->graph_mode));
            NodeWeights::serialize(std::move(kmer_counts),
                    utils::remove_suffix(config->outfbase, DBGSuccinct::kExtension)
                                                            + DBGSuccinct::kExtension);
        } else {
            graph_data.initialize_boss(boss_graph.get());
            graph.reset(new DBGSuccinct(boss_graph.release(), config->graph_mode));
        }

    } else if (config->graph_type == Config::GraphType::BITMAP && !config->dynamic) {

        logger->trace("Start reading data and extracting k-mers");
        // enumerate all suffixes
        std::vector<std::string> suffixes;
        if (config->suffix.size()) {
            suffixes = { config->suffix };
        } else {
            suffixes = kmer::KmerExtractor2Bit().generate_suffixes(config->suffix_len);
        }

        std::unique_ptr<DBGBitmapConstructor> constructor;
        std::vector<std::string> chunk_filenames;

        //one pass per suffix
        for (const std::string &suffix : suffixes) {
            timer.reset();

            if ((suffix.size() > 0 || suffixes.size() > 1)) {
                logger->trace("k-mer suffix: '{}'", suffix);
            }

            constructor.reset(
                new DBGBitmapConstructor(
                    config->k,
                    config->graph_mode,
                    config->count_width,
                    suffix,
                    get_num_threads(),
                    config->memory_available * kBytesInGigabyte
                )
            );

            push_sequences(files, *config, timer, constructor.get());

            if (!suffix.size()) {
                assert(suffixes.size() == 1);

                auto *bitmap_graph = new DBGBitmap(config->k);
                constructor->build_graph(bitmap_graph);
                graph.reset(bitmap_graph);

            } else {
                DBGBitmap::Chunk chunk = constructor->build_chunk();
                logger->trace("Graph chunk with {} k-mers was built in {} sec",
                              chunk.num_set_bits(), timer.elapsed());

                logger->trace("Serialize the graph chunk for suffix '{}'...", suffix);

                chunk_filenames.push_back(
                        utils::join_strings({ config->outfbase, suffix }, ".")
                        + DBGBitmap::kChunkFileExtension
                );
                std::ofstream out(chunk_filenames.back(), std::ios::binary);
                chunk.serialize(out);
                logger->trace("Serialization done in {} sec", timer.elapsed());
            }

            // only one chunk had to be constructed
            if (config->suffix.size())
                return 0;
        }

        if (suffixes.size() > 1) {
            assert(chunk_filenames.size());
            timer.reset();
            graph.reset(constructor->build_graph_from_chunks(chunk_filenames,
                                                             config->graph_mode,
                                                             get_verbose()));
        }

    } else {
        //slower method
        switch (config->graph_type) {

            case Config::GraphType::SUCCINCT:
                graph.reset(new DBGSuccinct(config->k, config->graph_mode));
                break;

            case Config::GraphType::HASH:
                graph.reset(new DBGHashOrdered(config->k, config->graph_mode));
                break;

            case Config::GraphType::HASH_PACKED:
                graph.reset(new DBGHashOrdered(config->k, config->graph_mode, true));
                break;

            case Config::GraphType::HASH_FAST:
                graph.reset(new DBGHashFast(config->k, config->graph_mode, true));
                break;

            case Config::GraphType::HASH_STR:
                if (config->graph_mode != DeBruijnGraph::BASIC) {
                    logger->warn("String hash-based de Bruijn graph"
                                 " currently supports only basic mode."
                                 " Basic mode will be used.");
                }
                // TODO: implement other graph modes
                graph.reset(new DBGHashString(config->k/*, config->graph_mode*/));
                break;

            case Config::GraphType::BITMAP:
                logger->error("Bitmap-graph construction"
                              " in dynamic regime is not supported");
                exit(1);

            case Config::GraphType::INVALID:
                assert(false);
        }
        assert(graph);

        for (const auto &file : files) {
            parse_sequences(file, *config,
                [&graph](std::string_view seq) { graph->add_sequence(seq); },
                [&graph](std::string_view seq, uint32_t) { graph->add_sequence(seq); }
            );
            logger->trace("Extracted all sequences from file {} in {} sec",
                          file, timer.elapsed());
        }

        if (config->count_kmers) {
            graph->add_extension(std::make_shared<NodeWeights>(graph->max_index() + 1,
                                                               config->count_width));
            auto node_weights = graph->get_extension<NodeWeights>();
            assert(node_weights->is_compatible(*graph));

            // set counts for reverse complement k-mers as well
            if (graph->get_mode() == DeBruijnGraph::CANONICAL)
                config->forward_and_reverse = true;

            for (const auto &file : files) {
                parse_sequences(file, *config,
                    [&graph,&node_weights](std::string_view seq) {
                        graph->map_to_nodes_sequentially(seq,
                            [&](auto node) { node_weights->add_weight(node, 1); }
                        );
                    },
                    [&graph,&node_weights](std::string_view seq, uint32_t count) {
                        graph->map_to_nodes_sequentially(seq,
                            [&](auto node) { node_weights->add_weight(node, count); }
                        );
                    }
                );
                logger->trace("Extracted all sequences from file {} in {} sec",
                              file, timer.elapsed());
            }

            // set back to false
            config->forward_and_reverse = false;
        }
    }

    logger->trace("Graph construction finished in {} sec", timer.elapsed());

    if (auto *dbg_succ = dynamic_cast<DBGSuccinct*>(graph.get())) {
        if (config->mark_dummy_kmers) {
            logger->trace("Detecting all dummy k-mers...");

            timer.reset();
            dbg_succ->mask_dummy_kmers(get_num_threads(), false);

            logger->trace("Dummy k-mer detection done in {} sec", timer.elapsed());
        }

        size_t suffix_length = std::min((size_t)config->node_suffix_length,
                                        dbg_succ->get_boss().get_k());

        if (suffix_length * log2(dbg_succ->get_boss().alph_size - 1) > 63) {
            logger->warn("Node ranges for k-mer suffixes longer than {} cannot be indexed",
                         static_cast<int>(63 / log2(dbg_succ->get_boss().alph_size - 1)));

        } else if (suffix_length) {
            logger->trace("Index all node ranges for suffixes of length {} in {:.2f} MB",
                          suffix_length,
                          std::pow(dbg_succ->get_boss().alph_size - 1, suffix_length)
                                * 2. * sizeof(uint64_t) * 1e-6);
            timer.reset();
            dbg_succ->get_boss().index_suffix_ranges(suffix_length);

            logger->trace("Indexing of node ranges took {} sec", timer.elapsed());
        }
    }

    graph->serialize(config->outfbase);
    graph->serialize_extensions(config->outfbase);

    return 0;
}

int concatenate_graph_chunks(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(config->outfbase.size());

    auto chunk_files = files;

    Timer timer;

    if (!files.size()) {
        assert(config->infbase.size());

        const auto sorted_suffixes = config->graph_type == Config::GraphType::SUCCINCT
                ? kmer::KmerExtractorBOSS().generate_suffixes(config->suffix_len)
                : kmer::KmerExtractor2Bit().generate_suffixes(config->suffix_len);

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
                                        boss::BOSS::Chunk::kFileExtension,
                                        DBGBitmap::kChunkFileExtension);
    }

    // collect results on an external merge or construction
    std::unique_ptr<DeBruijnGraph> graph;
    switch (config->graph_type) {
        case Config::GraphType::SUCCINCT: {
            boss::BOSS *boss = boss::BOSS::Chunk::build_boss_from_chunks(chunk_files, get_verbose());
            auto dbg_succ = std::make_unique<DBGSuccinct>(boss, config->graph_mode);

            logger->trace("Chunks concatenated in {} sec", timer.elapsed());

            if (config->clear_dummy) {
                logger->trace("Traverse source dummy edges,"
                              " remove redundant ones, and mark"
                              " those that cannot be removed");
                dbg_succ->mask_dummy_kmers(get_num_threads(), true);
            }

            size_t suffix_length = std::min((size_t)config->node_suffix_length,
                                            dbg_succ->get_boss().get_k());

            if (suffix_length * log2(dbg_succ->get_boss().alph_size - 1) > 63) {
                logger->warn("Node ranges for k-mer suffixes longer than {} cannot be indexed",
                             static_cast<int>(63 / log2(dbg_succ->get_boss().alph_size - 1)));

            } else if (suffix_length) {
                logger->trace("Index all node ranges for suffixes of length {} in {:.2f} MB",
                              suffix_length,
                              std::pow(dbg_succ->get_boss().alph_size - 1, suffix_length)
                                    * 2. * sizeof(uint64_t) * 1e-6);
                timer.reset();
                dbg_succ->get_boss().index_suffix_ranges(suffix_length);

                logger->trace("Indexing of node ranges took {} sec", timer.elapsed());
            }

            graph = std::move(dbg_succ);
            break;
        }
        case Config::GraphType::BITMAP: {
            graph.reset(DBGBitmapConstructor::build_graph_from_chunks(
                chunk_files, config->graph_mode, get_verbose()
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
            print_boss_stats(dynamic_cast<DBGSuccinct&>(*graph).get_boss());
        }
    }

    // graph output
    graph->serialize(config->outfbase);

    return 0;
}

} // namespace cli
} // namespace mtg
