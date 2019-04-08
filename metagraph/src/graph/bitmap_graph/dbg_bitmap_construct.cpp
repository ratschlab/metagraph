#include "dbg_bitmap_construct.hpp"


template <class KmerExtractor>
inline std::vector<typename KmerExtractor::TAlphabet>
encode_filter_suffix(const std::string &filter_suffix) {
    KmerExtractor kmer_extractor;
    std::vector<typename KmerExtractor::TAlphabet> filter_suffix_encoded;
    // TODO: cleanup
    std::transform(
        filter_suffix.begin(), filter_suffix.end(),
        std::back_inserter(filter_suffix_encoded),
        [&kmer_extractor](char c) {
            return kmer_extractor.encode(c);
        }
    );
    return filter_suffix_encoded;
}

template <typename KMER>
SDChunkConstructor<KMER>
::SDChunkConstructor(size_t k,
                     bool canonical_mode,
                     const std::string &filter_suffix,
                     size_t num_threads,
                     double memory_preallocated,
                     bool verbose)
      : kmer_collector_(k,
                        canonical_mode,
                        encode_filter_suffix<KmerExtractor2Bit>(filter_suffix),
                        num_threads,
                        memory_preallocated,
                        verbose) {}

/**
 * Initialize graph chunk from a list of sorted kmers.
 */
template <typename KMER>
DBGSD::Chunk* SDChunkConstructor<KMER>
::build_chunk() {
    kmer_collector_.join();
    std::unique_ptr<DBGSD::Chunk> chunk {
        new DBGSD::Chunk(
            [&](const auto &index_callback) {
                kmer_collector_.call_kmers(
                    [&](const auto &kmer) {
                        index_callback(typename KMER::WordType(1u) + kmer.data());
                    }
                );
            },
            (1llu << (get_k() * KMER::kBitsPerChar)) + 1,
            kmer_collector_.size()
        )
    };
    assert(chunk.get());

    return chunk.release();
}

DBGSD* DBGSDConstructor
::build_graph_from_chunks(const std::vector<std::unique_ptr<DBGSD::Chunk>> &chunks,
                          bool canonical_mode,
                          bool verbose) {
    if (chunks.empty())
        return new DBGSD(2);

    uint64_t cumulative_size = 1;

    for (const auto &chunk : chunks) {
        cumulative_size += chunk->num_set_bits();
    }

    if (verbose)
        std::cout << "Cumulative size of chunks: "
                  << cumulative_size - 1 << std::endl;

    std::unique_ptr<DBGSD> graph{
        new DBGSD(2)
    };

    graph->kmers_ = DBGSD::Chunk(
        [&](const auto &index_callback) {
            index_callback(0);
            for (size_t i = 0; i < chunks.size(); ++i) {
                if (verbose) {
                    std::cout << "Chunk "
                              << i
                              << " loaded..." << std::flush;
                }

                chunks[i]->call_ones(index_callback);

                if (verbose) {
                    std::cout << " concatenated" << std::endl;
                }
            }
        },
        chunks[0]->size(), cumulative_size
    );

    graph->canonical_mode_ = canonical_mode;

    assert(!(sdsl::bits::hi(graph->kmers_.size()) % KmerExtractor2Bit::kLogSigma));

    graph->k_ = sdsl::bits::hi(graph->kmers_.size()) / KmerExtractor2Bit::kLogSigma;

    return graph.release();
}

DBGSD* DBGSDConstructor
::build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                          bool canonical_mode,
                          bool verbose) {
    if (chunk_filenames.empty())
        return new DBGSD(2);

    std::vector<std::unique_ptr<DBGSD::Chunk>> chunks;

    for (auto file : chunk_filenames) {
        file = utils::remove_suffix(file, ".dbgsdchunk") + ".dbgsdchunk";

        std::ifstream chunk_in(file, std::ios::binary);

        if (!chunk_in.good()) {
            std::cerr << "ERROR: input file " << file << " corrupted" << std::endl;
            exit(1);
        }

        chunks.emplace_back(new DBGSD::Chunk());
        chunks.back()->load(chunk_in);

        assert(chunks.empty() || chunks.back()->size() == chunks.front()->size());
    }

    return build_graph_from_chunks(chunks, canonical_mode, verbose);
}

ISDChunkConstructor*
ISDChunkConstructor
::initialize(size_t k,
             bool canonical_mode,
             const std::string &filter_suffix,
             size_t num_threads,
             double memory_preallocated,
             bool verbose) {
    using Extractor = KmerExtractor2Bit;

    if (k * Extractor::kLogSigma <= 64) {
        return new SDChunkConstructor<typename Extractor::Kmer64>(
            k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
        );
    } else if (k * Extractor::kLogSigma <= 128) {
        return new SDChunkConstructor<typename Extractor::Kmer128>(
            k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
        );
    } else {
        return new SDChunkConstructor<typename Extractor::Kmer256>(
            k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
        );
    }
}

void DBGSDConstructor::build_graph(DBGSD *graph) {
    auto chunk = constructor_->build_chunk();
    graph->k_ = constructor_->get_k();
    graph->canonical_mode_ = constructor_->is_canonical_mode();
    graph->kmers_ = decltype(graph->kmers_)(
        [&](const auto &index_callback) {
            index_callback(0);
            chunk->call_ones(index_callback);
        },
        chunk->size(), chunk->num_set_bits() + 1
    );
    delete chunk;
}

DBGSD::Chunk* DBGSDConstructor::build_chunk() {
    return constructor_->build_chunk();
}
