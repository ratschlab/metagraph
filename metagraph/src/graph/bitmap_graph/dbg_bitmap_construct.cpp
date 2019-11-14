#include "dbg_bitmap_construct.hpp"

#include <progress_bar.hpp>

#include "utils/template_utils.hpp"
#include "utils/algorithms.hpp"
#include "common/sorted_set.hpp"
#include "common/sorted_multiset.hpp"
#include "kmer/kmer_collector.hpp"
#include "node_weights.hpp"


template <typename KmerStorage>
class BitmapChunkConstructor : public IBitmapChunkConstructor {
    friend IBitmapChunkConstructor;

  private:
    BitmapChunkConstructor(size_t k,
                           bool canonical_mode = false,
                           const std::string &filter_suffix = "",
                           size_t num_threads = 1,
                           double memory_preallocated = 0,
                           bool verbose = false);

    void add_sequence(std::string&& sequence, uint64_t count) {
        kmer_collector_.add_sequence(std::move(sequence), count);
    }

    void add_sequences(std::function<void(CallString)> generate_sequences) {
        kmer_collector_.add_sequences(generate_sequences);
    }

    size_t get_k() const { return kmer_collector_.get_k(); }

    bool is_canonical_mode() const { return kmer_collector_.is_both_strands_mode(); }

    DBGBitmap::Chunk* build_chunk();

    sdsl::int_vector<> get_weights(uint8_t bits_per_count = 8);

    KmerStorage kmer_collector_;
};


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

template <typename KmerStorage>
BitmapChunkConstructor<KmerStorage>
::BitmapChunkConstructor(size_t k,
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

DBGBitmapConstructor::DBGBitmapConstructor(size_t k,
                                           bool canonical_mode,
                                           uint8_t bits_per_count,
                                           const std::string &filter_suffix,
                                           size_t num_threads,
                                           double memory_preallocated,
                                           bool verbose)
      : constructor_(IBitmapChunkConstructor::initialize(
            k,
            canonical_mode,
            bits_per_count > 0,
            filter_suffix,
            num_threads,
            memory_preallocated,
            verbose)
        ),
        bits_per_count_(bits_per_count) {}

template <typename KmerStorage>
sdsl::int_vector<>
BitmapChunkConstructor<KmerStorage>::get_weights(uint8_t bits_per_count) {
    if constexpr(utils::is_pair<typename KmerStorage::Value>::value) {
        const auto &kmers = kmer_collector_.data();

        sdsl::int_vector<> weights(kmers.size() + 1, 0, bits_per_count);

        const uint64_t max_count = utils::max_ull(weights.width());

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
template <typename KmerStorage>
DBGBitmap::Chunk* BitmapChunkConstructor<KmerStorage>
::build_chunk() {
    using KMER = typename KmerStorage::Key;

    const auto &kmers = kmer_collector_.data();
    std::unique_ptr<DBGBitmap::Chunk> chunk {
        new DBGBitmap::Chunk(
            [&](const auto &index_callback) {
                std::for_each(kmers.begin(), kmers.end(), [&](const typename KmerStorage::Value &kmer) {
                    if constexpr(utils::is_pair<typename KmerStorage::Value>::value) {
                        index_callback(typename KMER::WordType(1u) + kmer.first.data());
                    } else {
                        index_callback(typename KMER::WordType(1u) + kmer.data());
                    }
                });
            },
            (1llu << (get_k() * KMER::kBitsPerChar)) + 1,
            kmers.size()
        )
    };
    assert(chunk.get());

    return chunk.release();
}

DBGBitmap* DBGBitmapConstructor
::build_graph_from_chunks(uint64_t size,
                          uint64_t num_kmers,
                          const std::function<DBGBitmap::Chunk(void)> &next_chunk,
                          bool canonical_mode) {
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

    graph->canonical_mode_ = canonical_mode;
    graph->complete_ = false;

    assert(!(sdsl::bits::hi(graph->kmers_.size()) % KmerExtractor2Bit::bits_per_char));

    graph->k_ = sdsl::bits::hi(graph->kmers_.size()) / KmerExtractor2Bit::bits_per_char;

    return graph.release();
}

DBGBitmap* DBGBitmapConstructor
::build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                          bool canonical_mode,
                          bool verbose) {
    if (chunk_filenames.empty())
        return new DBGBitmap(2);

    uint64_t size = 0;
    uint64_t cumulative_size = 1;

    for (size_t i = 0; i < chunk_filenames.size(); ++i) {
        const auto &chunk_filename = chunk_filenames.at(i);
        std::unique_ptr<DBGBitmap::Chunk> chunk(new DBGBitmap::Chunk());

        std::ifstream chunk_in(chunk_filename, std::ios::binary);
        chunk->load(chunk_in);

        if (!i) {
            size = chunk->size();
        } else if (size != chunk->size()) {
            std::cerr << "ERROR: inconsistent graph chunks" << std::endl;
            exit(1);
        }

        cumulative_size += chunk->num_set_bits();
    }
    assert(size > 0);

    if (verbose)
        std::cout << "Cumulative size of chunks: "
                  << cumulative_size - 1 << std::endl;

    size_t chunk_idx = 0;

    //TODO configure stream for verbose output globally and refactor
    std::ofstream null_ofstream;
    ProgressBar progress_bar(chunk_filenames.size(), "Processing chunks", verbose ? std::cout : null_ofstream);

    return build_graph_from_chunks(
        size,
        cumulative_size,
        [&]() {
            if (chunk_idx >= chunk_filenames.size())
                return DBGBitmap::Chunk();

            auto chunk_filename = chunk_filenames.at(chunk_idx);
            chunk_filename = utils::remove_suffix(chunk_filename, DBGBitmap::kChunkFileExtension)
                                + DBGBitmap::kChunkFileExtension;

            std::ifstream chunk_in(chunk_filename, std::ios::binary);

            if (!chunk_in.good()) {
                std::cerr << "ERROR: input file " << chunk_filename << " corrupted" << std::endl;
                exit(1);
            }

            DBGBitmap::Chunk chunk;
            chunk.load(chunk_in);

            chunk_idx++;
            ++progress_bar;
            return chunk;
        },
        canonical_mode
    );
}

template <typename KMER>
using KmerCollector = KmerStorage<KMER,
                                  KmerExtractor2Bit,
                                  SortedSet<KMER, Vector<KMER>>>;

template <typename KMER,
          typename KmerCount = uint8_t>
using KmerCounter = KmerStorage<KMER,
                                KmerExtractor2Bit,
                                SortedMultiset<KMER, KmerCount, Vector<std::pair<KMER, KmerCount>>>>;

IBitmapChunkConstructor*
IBitmapChunkConstructor
::initialize(size_t k,
             bool canonical_mode,
             bool count_kmers,
             const std::string &filter_suffix,
             size_t num_threads,
             double memory_preallocated,
             bool verbose) {
    using Extractor = KmerExtractor2Bit;

    if (count_kmers) {
        if (k * Extractor::bits_per_char <= 64) {
            return new BitmapChunkConstructor<KmerCounter<typename Extractor::Kmer64, uint8_t>>(
                k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
            );
        } else if (k * Extractor::bits_per_char <= 128) {
            return new BitmapChunkConstructor<KmerCounter<typename Extractor::Kmer128, uint8_t>>(
                k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
            );
        } else {
            return new BitmapChunkConstructor<KmerCounter<typename Extractor::Kmer256, uint8_t>>(
                k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
            );
        }
    } else {
        if (k * Extractor::bits_per_char <= 64) {
            return new BitmapChunkConstructor<KmerCollector<typename Extractor::Kmer64>>(
                k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
            );
        } else if (k * Extractor::bits_per_char <= 128) {
            return new BitmapChunkConstructor<KmerCollector<typename Extractor::Kmer128>>(
                k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
            );
        } else {
            return new BitmapChunkConstructor<KmerCollector<typename Extractor::Kmer256>>(
                k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
            );
        }
    }
}

void DBGBitmapConstructor::build_graph(DBGBitmap *graph) {
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
    graph->complete_ = false;

    if (bits_per_count_) {
        graph->add_extension(
            std::make_shared<NodeWeights>(constructor_->get_weights(bits_per_count_))
        );
    }
}
