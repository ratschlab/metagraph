#include "dbg_bitmap_construct.hpp"

#include <progress_bar.hpp>

#include "common/logger.hpp"
#include "common/sorted_sets/sorted_set.hpp"
#include "common/sorted_sets/sorted_multiset.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/template_utils.hpp"
#include "kmer/kmer_collector.hpp"
#include "graph/graph_extensions/node_weights.hpp"


namespace mtg {
namespace graph {

using mtg::common::logger;
using mtg::kmer::KmerExtractor2Bit;


template <typename KmerCollector>
class BitmapChunkConstructor : public IBitmapChunkConstructor {
    friend IBitmapChunkConstructor;

    template <template <typename KMER> class Collector, typename... Args>
    friend IBitmapChunkConstructor*
    initialize_bitmap_chunk_constructor(size_t k, const Args& ...args);

  private:
    BitmapChunkConstructor(size_t k,
                           DeBruijnGraph::Mode mode,
                           const std::string &filter_suffix,
                           size_t num_threads,
                           double memory_preallocated);

    void add_sequence(std::string_view sequence, uint64_t count) {
        kmer_collector_.add_sequence(sequence, count);
    }

    void add_sequences(std::vector<std::string>&& sequences) {
        kmer_collector_.add_sequences(std::move(sequences));
    }

    void add_sequences(std::vector<std::pair<std::string, uint64_t>>&& sequences) {
        kmer_collector_.add_sequences(std::move(sequences));
    }

    size_t get_k() const { return kmer_collector_.get_k(); }

    DeBruijnGraph::Mode get_mode() const { return mode_; }

    DBGBitmap::Chunk build_chunk();

    sdsl::int_vector<> get_weights(uint8_t bits_per_count);

    DeBruijnGraph::Mode mode_;
    KmerCollector kmer_collector_;
};


template <class KmerExtractor>
inline std::vector<typename KmerExtractor::TAlphabet>
encode_filter_suffix(const std::string &filter_suffix) {
    KmerExtractor kmer_extractor;
    std::vector<typename KmerExtractor::TAlphabet> filter_suffix_encoded;
    for (char c : filter_suffix) {
        filter_suffix_encoded.push_back(kmer_extractor.encode(c));
    }
    return filter_suffix_encoded;
}

template <typename KmerCollector>
BitmapChunkConstructor<KmerCollector>
::BitmapChunkConstructor(size_t k,
                         DeBruijnGraph::Mode mode,
                         const std::string &filter_suffix,
                         size_t num_threads,
                         double memory_preallocated)
      : mode_(mode),
        kmer_collector_(k,
                        mode == DeBruijnGraph::CANONICAL
                            ? KmerCollector::BOTH
                            : KmerCollector::BASIC,
                        encode_filter_suffix<KmerExtractor2Bit>(filter_suffix),
                        num_threads,
                        memory_preallocated) {}

DBGBitmapConstructor::DBGBitmapConstructor(size_t k,
                                           DeBruijnGraph::Mode mode,
                                           uint8_t bits_per_count,
                                           const std::string &filter_suffix,
                                           size_t num_threads,
                                           double memory_preallocated) {
    constructor_.reset(
        IBitmapChunkConstructor::initialize(k,
                                            mode,
                                            bits_per_count,
                                            filter_suffix,
                                            num_threads,
                                            memory_preallocated)
    );
    bits_per_count_ = bits_per_count;
}

template <typename KmerCollector>
sdsl::int_vector<>
BitmapChunkConstructor<KmerCollector>::get_weights(uint8_t bits_per_count) {
    if constexpr(utils::is_pair_v<typename KmerCollector::Value>) {
        const auto &kmers = kmer_collector_.data();

        sdsl::int_vector<> weights(kmers.size() + 1, 0, bits_per_count);

        const uint64_t max_count = sdsl::bits::lo_set[weights.width()];

        for (size_t i = 0; i < kmers.size(); ++i) {
            weights[i + 1] = std::min(static_cast<uint64_t>(kmers[i].second), max_count);
        }

        return weights;

    } else {
        std::ignore = bits_per_count;
        throw std::runtime_error("Error: count k-mers to get weights");
    }
};

/**
 * Initialize graph chunk from a list of sorted kmers.
 */
template <typename KmerCollector>
DBGBitmap::Chunk BitmapChunkConstructor<KmerCollector>::build_chunk() {
    using KMER = typename KmerCollector::Kmer;

    const auto &kmers = kmer_collector_.data();

    return DBGBitmap::Chunk(
            [&](const auto &index_callback) {
                std::for_each(kmers.begin(), kmers.end(),
                    [&](const auto &kmer) { index_callback(utils::get_first(kmer) + 1); }
                );
            },
            (1llu << (get_k() * KMER::kBitsPerChar)) + 1,
            kmers.size()
    );
}

DBGBitmap* DBGBitmapConstructor
::build_graph_from_chunks(uint64_t size,
                          uint64_t num_kmers,
                          const std::function<DBGBitmap::Chunk(void)> &next_chunk,
                          DeBruijnGraph::Mode mode) {
    auto graph = std::make_unique<DBGBitmap>(2);


    graph->kmers_ = DBGBitmap::Chunk(
        [&](const auto &index_callback) {
            DBGBitmap::Chunk chunk;

            index_callback(0);

            while ((chunk = next_chunk()).size()) {
                chunk.call_ones(index_callback);
            }
        },
        size,
        num_kmers
    );

    graph->mode_ = mode;
    graph->complete_ = false;

    assert(!(sdsl::bits::hi(graph->kmers_.size()) % KmerExtractor2Bit::bits_per_char));

    graph->k_ = sdsl::bits::hi(graph->kmers_.size()) / KmerExtractor2Bit::bits_per_char;

    return graph.release();
}

DBGBitmap* DBGBitmapConstructor
::build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                          DeBruijnGraph::Mode mode,
                          bool verbose) {
    if (chunk_filenames.empty())
        return new DBGBitmap(2);

    uint64_t size = 0;
    uint64_t cumulative_size = 1;

    for (size_t i = 0; i < chunk_filenames.size(); ++i) {
        const auto &chunk_filename = chunk_filenames.at(i);
        DBGBitmap::Chunk chunk;

        std::unique_ptr<std::ifstream> chunk_in = utils::open_ifstream(chunk_filename);
        chunk.load(*chunk_in);

        if (!i) {
            size = chunk.size();
        } else if (size != chunk.size()) {
            logger->error("Inconsistent graph chunks");
            exit(1);
        }

        cumulative_size += chunk.num_set_bits();
    }
    assert(size > 0);

    if (verbose)
        logger->info("Cumulative size of chunks: {}", cumulative_size - 1);

    size_t chunk_idx = 0;

    ProgressBar progress_bar(chunk_filenames.size(), "Processing chunks", std::cout, !verbose);

    return build_graph_from_chunks(
        size,
        cumulative_size,
        [&]() {
            if (chunk_idx >= chunk_filenames.size())
                return DBGBitmap::Chunk();

            auto chunk_filename = utils::make_suffix(chunk_filenames.at(chunk_idx),
                                                     DBGBitmap::kChunkFileExtension);

            std::unique_ptr<std::ifstream> chunk_in = utils::open_ifstream(chunk_filename);

            if (!chunk_in->good()) {
                logger->error("Input file {} corrupted", chunk_filename);
                exit(1);
            }

            DBGBitmap::Chunk chunk;
            chunk.load(*chunk_in);

            chunk_idx++;
            ++progress_bar;
            return chunk;
        },
        mode
    );
}

template <template <typename KMER> class KmerCollector, typename... Args>
IBitmapChunkConstructor*
initialize_bitmap_chunk_constructor(size_t k, const Args& ...args) {
    if (k < 1 || k > 63 / KmerExtractor2Bit::bits_per_char) {
        logger->error("For bitmap graph, k must be between 1 and {}",
                      63 / KmerExtractor2Bit::bits_per_char);
        exit(1);
    }

    if (k * KmerExtractor2Bit::bits_per_char <= 64) {
        return new BitmapChunkConstructor<KmerCollector<KmerExtractor2Bit::Kmer64>>(k, args...);
    } else if (k * KmerExtractor2Bit::bits_per_char <= 128) {
        return new BitmapChunkConstructor<KmerCollector<KmerExtractor2Bit::Kmer128>>(k, args...);
    } else {
        return new BitmapChunkConstructor<KmerCollector<KmerExtractor2Bit::Kmer256>>(k, args...);
    }
}

template <typename KMER>
using KmerSet
        = kmer::KmerCollector<KMER, KmerExtractor2Bit, common::SortedSet<typename KMER::WordType>>;
template <typename KMER>
using KmerMultsetVector8 = kmer::KmerCollector<KMER,
                                               KmerExtractor2Bit,
                                               common::SortedMultiset<typename KMER::WordType, uint8_t>>;
template <typename KMER>
using KmerMultsetVector16 = kmer::KmerCollector<KMER,
                                                KmerExtractor2Bit,
                                                common::SortedMultiset<typename KMER::WordType, uint16_t>>;
template <typename KMER>
using KmerMultsetVector32 = kmer::KmerCollector<KMER,
                                                KmerExtractor2Bit,
                                                common::SortedMultiset<typename KMER::WordType, uint32_t>>;

IBitmapChunkConstructor*
IBitmapChunkConstructor::initialize(size_t k,
                                    DeBruijnGraph::Mode mode,
                                    uint8_t bits_per_count,
                                    const std::string &filter_suffix,
                                    size_t num_threads,
                                    double memory_preallocated) {
#define OTHER_ARGS k, mode, filter_suffix, num_threads, memory_preallocated

    if (!bits_per_count) {
        return initialize_bitmap_chunk_constructor<KmerSet>(OTHER_ARGS);
    } else if (bits_per_count <= 8) {
        return initialize_bitmap_chunk_constructor<KmerMultsetVector8>(OTHER_ARGS);
    } else if (bits_per_count <= 16) {
        return initialize_bitmap_chunk_constructor<KmerMultsetVector16>(OTHER_ARGS);
    } else if (bits_per_count <= 32) {
        return initialize_bitmap_chunk_constructor<KmerMultsetVector32>(OTHER_ARGS);
    } else {
        throw std::runtime_error("Error: trying to allocate too many bits per k-mer count");
    }
}

void DBGBitmapConstructor::build_graph(DBGBitmap *graph) {
    auto chunk = constructor_->build_chunk();
    graph->k_ = constructor_->get_k();
    graph->mode_ = constructor_->get_mode();
    graph->kmers_ = decltype(graph->kmers_)(
        [&](const auto &index_callback) {
            index_callback(0);
            chunk.call_ones(index_callback);
        },
        chunk.size(), chunk.num_set_bits() + 1
    );
    chunk = DBGBitmap::Chunk();
    graph->complete_ = false;

    if (bits_per_count_) {
        graph->add_extension(
            std::make_shared<NodeWeights>(constructor_->get_weights(bits_per_count_))
        );
    }
}

} // namespace graph
} // namespace mtg
