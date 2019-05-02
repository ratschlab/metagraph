#include "dbg_bitmap_construct.hpp"
#include <progress_bar.hpp>


template <typename KMER>
class BitmapChunkConstructor : public IBitmapChunkConstructor {
    friend IBitmapChunkConstructor;

  private:
    BitmapChunkConstructor(size_t k,
                           bool canonical_mode = false,
                           const std::string &filter_suffix = "",
                           size_t num_threads = 1,
                           double memory_preallocated = 0,
                           bool verbose = false);

    void add_sequence(const std::string &sequence) {
        kmer_collector_.add_sequence(sequence);
    }

    void add_sequences(std::function<void(CallString)> generate_sequences) {
        kmer_collector_.add_sequences(generate_sequences);
    }

    size_t get_k() const { return kmer_collector_.get_k(); }

    bool is_canonical_mode() const { return kmer_collector_.is_both_strands_mode(); }

    DBGBitmap::Chunk* build_chunk();

    KmerCollector<KMER, KmerExtractor2Bit> kmer_collector_;
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

template <typename KMER>
BitmapChunkConstructor<KMER>
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
                                           const std::string &filter_suffix,
                                           size_t num_threads,
                                           double memory_preallocated,
                                           bool verbose)
      : constructor_(IBitmapChunkConstructor::initialize(
            k,
            canonical_mode,
            filter_suffix,
            num_threads,
            memory_preallocated,
            verbose)
        ) {}

/**
 * Initialize graph chunk from a list of sorted kmers.
 */
template <typename KMER>
DBGBitmap::Chunk* BitmapChunkConstructor<KMER>
::build_chunk() {
    kmer_collector_.join();
    std::unique_ptr<DBGBitmap::Chunk> chunk {
        new DBGBitmap::Chunk(
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

DBGBitmap* DBGBitmapConstructor
::build_graph_from_chunks(const std::function<std::unique_ptr<DBGBitmap::Chunk>(void)> &next_chunk,
                          uint64_t size,
                          uint64_t cumulative_size,
                          bool canonical_mode) {
    std::unique_ptr<DBGBitmap> graph{
        new DBGBitmap(2)
    };

    graph->kmers_ = DBGBitmap::Chunk(
        [&](const auto &index_callback) {
            std::unique_ptr<DBGBitmap::Chunk> chunk;

            index_callback(0);

            while (chunk = next_chunk()) {
                chunk->call_ones(index_callback);
            }
        },
        size, cumulative_size
    );

    graph->canonical_mode_ = canonical_mode;

    graph->complete_ = false;

    assert(!(sdsl::bits::hi(graph->kmers_.size()) % KmerExtractor2Bit::kLogSigma));

    graph->k_ = sdsl::bits::hi(graph->kmers_.size()) / KmerExtractor2Bit::kLogSigma;

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

        if (i > 0)
            assert(chunk->size() == size);
        else
            size = chunk->size();

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

    return build_graph_from_chunks([&]() {
        if (chunk_idx < chunk_filenames.size()) {
            auto chunk_filename = chunk_filenames.at(chunk_idx);
            chunk_filename = utils::remove_suffix(chunk_filename, ".dbgsdchunk") + ".dbgsdchunk";

            std::ifstream chunk_in(chunk_filename, std::ios::binary);

            if (!chunk_in.good()) {
                std::cerr << "ERROR: input file " << chunk_filename << " corrupted" << std::endl;
                exit(1);
            }

            std::unique_ptr<DBGBitmap::Chunk> chunk(new DBGBitmap::Chunk());
            chunk->load(chunk_in);

            chunk_idx++;
            ++progress_bar;
            return std::move(chunk);
        } else {
            return std::unique_ptr<DBGBitmap::Chunk>(nullptr);
        }
    }, size, cumulative_size, canonical_mode);
}

IBitmapChunkConstructor*
IBitmapChunkConstructor
::initialize(size_t k,
             bool canonical_mode,
             const std::string &filter_suffix,
             size_t num_threads,
             double memory_preallocated,
             bool verbose) {
    using Extractor = KmerExtractor2Bit;

    if (k * Extractor::kLogSigma <= 64) {
        return new BitmapChunkConstructor<typename Extractor::Kmer64>(
            k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
        );
    } else if (k * Extractor::kLogSigma <= 128) {
        return new BitmapChunkConstructor<typename Extractor::Kmer128>(
            k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
        );
    } else {
        return new BitmapChunkConstructor<typename Extractor::Kmer256>(
            k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
        );
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
}
