#include "boss_chunk_construct.hpp"

#include <ips4o.hpp>

#include "unix_tools.hpp"
#include "boss_chunk.hpp"
#include "kmer_collector.hpp"


template <typename T>
void dont_erase_redundant_dummy_kmers(Vector<T> *) {}

template <typename KMER, typename COUNT>
inline KMER& get_kmer(Vector<std::pair<KMER, COUNT>> &kmers, uint64_t i) {
    return kmers[i].first;
}

template <typename KMER>
inline KMER& get_kmer(Vector<KMER> &kmers, uint64_t i) {
    return kmers[i];
}

template <typename T>
void erase_redundant_dummy_kmers(Vector<T> *kmers) {
    using KMER = std::remove_reference_t<decltype(get_kmer(*kmers, 0))>;

    assert(std::is_sorted(kmers->begin(), kmers->end(), utils::LessFirst<T>()));
    assert(std::unique(kmers->begin(), kmers->end(), utils::EqualFirst<T>()) == kmers->end());

    if (kmers->size() < 2)
        return;

    // last k-mer is never redundant. Start with the next one.
    uint64_t last = kmers->size() - 1;

    typename KMER::CharType edge_label, node_last_char;

    std::vector<uint64_t> last_kmer(1llu << KMER::kBitsPerChar, kmers->size());

    last_kmer[get_kmer(*kmers, last)[0]] = last;

    for (int64_t i = last - 1; i >= 0; --i) {
        const KMER &kmer = get_kmer(*kmers, i);
        node_last_char = kmer[1];
        edge_label = kmer[0];

        // assert((edge_label || node_last_char)
        //             && "dummy k-mer cannot be both source and sink dummy");

        if (!edge_label) {
            // sink dummy k-mer

            // skip if redundant
            if (node_last_char && KMER::compare_suffix(kmer, get_kmer(*kmers, last), 0))
                continue;

        } else if (!node_last_char) {
            // source dummy k-mer

            // skip if redundant
            if (last_kmer[edge_label] < kmers->size()
                    && KMER::compare_suffix(kmer, get_kmer(*kmers, last_kmer[edge_label]), 1))
                continue;
        }

        // the k-mer is either not dummy, or not redundant -> keep the k-mer
        kmers->at(--last) = kmers->at(i);
        last_kmer[edge_label] = last;
    }

    kmers->erase(kmers->begin(), kmers->begin() + last);
}


template <class T>
void sort_and_remove_duplicates(Vector<T> *array,
                                size_t num_threads,
                                size_t offset) {
    ips4o::parallel::sort(array->begin() + offset, array->end(),
                          utils::LessFirst<T>(),
                          num_threads);
    // remove duplicates
    auto unique_end = std::unique(array->begin() + offset, array->end(), utils::EqualFirst<T>());
    array->erase(unique_end, array->end());
}

template <typename T>
void shrink_kmers(Vector<T> *kmers,
                  size_t num_threads,
                  bool verbose,
                  size_t offset) {
    if (verbose) {
        std::cout << "Allocated capacity exceeded, filter out non-unique k-mers..."
                  << std::flush;
    }

    size_t prev_num_kmers = kmers->size();
    sort_and_remove_duplicates(kmers, num_threads, offset);

    if (verbose) {
        std::cout << " done. Number of kmers reduced from " << prev_num_kmers
                                                  << " to " << kmers->size() << ", "
                  << (kmers->size() * sizeof(T) >> 20) << "Mb" << std::endl;
    }
}

template <typename KMER, typename COUNT>
inline KMER& push_back(Vector<std::pair<KMER, COUNT>> &kmers, const KMER &kmer) {
    kmers.emplace_back(kmer, 0);
    return kmers.back().first;
}

template <typename KMER>
inline KMER& push_back(Vector<KMER> &kmers, const KMER &kmer) {
    kmers.push_back(kmer);
    return kmers.back();
}

// Although this function could be parallelized better,
// the experiments show it's already fast enough.
// k is node length
template <typename T>
void recover_source_dummy_nodes(size_t k,
                                Vector<T> *kmers,
                                size_t num_threads,
                                bool verbose) {
    using KMER = std::remove_reference_t<decltype(get_kmer(*kmers, 0))>;

    size_t dummy_begin = kmers->size();
    size_t num_dummy_parent_kmers = 0;

    for (size_t i = 0; i < dummy_begin; ++i) {
        const KMER &kmer = get_kmer(*kmers, i);
        // we never add reads shorter than k
        assert(kmer[1] != 0 || kmer[0] != 0 || kmer[k] == 0);

        auto node_last_char = kmer[1];
        auto edge_label = kmer[0];
        // check if it's not a source dummy kmer
        if (node_last_char || !edge_label)
            continue;

        num_dummy_parent_kmers++;

        if (kmers->size() + 1 > kmers->capacity())
            shrink_kmers(kmers, num_threads, verbose, dummy_begin);

        push_back(*kmers, kmer).to_prev(k + 1, BOSS::kSentinelCode);
    }
    if (verbose) {
        std::cout << "Number of dummy k-mers with dummy prefix of length 1: "
                  << num_dummy_parent_kmers << std::endl;
    }
    sort_and_remove_duplicates(kmers, num_threads, dummy_begin);

    if (verbose) {
        std::cout << "Number of dummy k-mers with dummy prefix of length 2: "
                  << kmers->size() - dummy_begin << std::endl;
    }

    for (size_t c = 3; c < k + 1; ++c) {
        size_t succ_dummy_begin = dummy_begin;
        dummy_begin = kmers->size();

        for (size_t i = succ_dummy_begin; i < dummy_begin; ++i) {
            if (kmers->size() + 1 > kmers->capacity())
                shrink_kmers(kmers, num_threads, verbose, dummy_begin);

            push_back(*kmers, get_kmer(*kmers, i)).to_prev(k + 1, BOSS::kSentinelCode);
        }
        sort_and_remove_duplicates(kmers, num_threads, dummy_begin);

        if (verbose) {
            std::cout << "Number of dummy k-mers with dummy prefix of length " << c
                      << ": " << kmers->size() - dummy_begin << std::endl;
        }
    }
    ips4o::parallel::sort(kmers->begin(), kmers->end(), utils::LessFirst<T>(), num_threads);
}

template <class KmerExtractor>
inline std::vector<typename KmerExtractor::TAlphabet>
encode_filter_suffix_boss(const std::string &filter_suffix) {
    KmerExtractor kmer_extractor;
    std::vector<typename KmerExtractor::TAlphabet> filter_suffix_encoded;
    std::transform(
        filter_suffix.begin(), filter_suffix.end(),
        std::back_inserter(filter_suffix_encoded),
        [&kmer_extractor](char c) {
            return c == BOSS::kSentinel
                            ? BOSS::kSentinelCode
                            : kmer_extractor.encode(c);
        }
    );
    return filter_suffix_encoded;
}

template <typename KmerStorage>
class BOSSChunkConstructor : public IBOSSChunkConstructor {
    friend IBOSSChunkConstructor;

  private:
    BOSSChunkConstructor(size_t k,
                         bool canonical_mode = false,
                         const std::string &filter_suffix = "",
                         size_t num_threads = 1,
                         double memory_preallocated = 0,
                         bool verbose = false)
          : kmer_storage_(k + 1,
                          canonical_mode,
                          encode_filter_suffix_boss<KmerExtractor>(filter_suffix),
                          num_threads,
                          memory_preallocated,
                          verbose,
                          filter_suffix.empty() ? erase_redundant_dummy_kmers<typename KmerStorage::Value>
                                                : dont_erase_redundant_dummy_kmers<typename KmerStorage::Value>) {
        if (filter_suffix == std::string(filter_suffix.size(), BOSS::kSentinel)) {
            kmer_storage_.insert_dummy(std::vector<KmerExtractor::TAlphabet>(k + 1, BOSS::kSentinelCode));
        }
    }

    void add_kmer(std::string&& kmer, uint32_t count) {
        assert(kmer.size() == get_k());
        kmer_storage_.add_kmer(std::move(kmer), count);
    }

    void add_sequence(std::string&& sequence) {
        kmer_storage_.add_sequence(std::move(sequence));
    }

    void add_sequences(std::function<void(CallString)> generate_sequences) {
        kmer_storage_.add_sequences(generate_sequences);
    }

    BOSS::Chunk* build_chunk() {
        auto &kmers = kmer_storage_.data();

        if (!kmer_storage_.suffix_length()) {
            if (kmer_storage_.verbose()) {
                std::cout << "Reconstructing all required dummy source k-mers..."
                          << std::endl;
            }
            Timer timer;

            // kmer_collector stores (BOSS::k_ + 1)-mers
            recover_source_dummy_nodes(kmer_storage_.get_k() - 1,
                                       &kmers,
                                       kmer_storage_.num_threads(),
                                       kmer_storage_.verbose());

            if (kmer_storage_.verbose())
                std::cout << "Dummy source k-mers were reconstructed in "
                          << timer.elapsed() << "sec" << std::endl;
        }

        // kmer_collector stores (BOSS::k_ + 1)-mers
        BOSS::Chunk *result = new BOSS::Chunk(kmer_storage_.alphabet_size(),
                                              kmer_storage_.get_k() - 1,
                                              kmers);
        kmer_storage_.clear();

        return result;
    }

    uint64_t get_k() const { return kmer_storage_.get_k() - 1; }

    KmerStorage kmer_storage_;
};

std::unique_ptr<IBOSSChunkConstructor>
IBOSSChunkConstructor
::initialize(size_t k,
             bool canonical_mode,
             bool count_kmers,
             const std::string &filter_suffix,
             size_t num_threads,
             double memory_preallocated,
             bool verbose) {
    using Extractor = KmerExtractor;

    if (count_kmers) {
        if ((k + 1) * Extractor::bits_per_char <= 64) {
            return std::unique_ptr<IBOSSChunkConstructor>(
                new BOSSChunkConstructor<KmerCounter<typename Extractor::Kmer64, Extractor, uint8_t>>(
                    k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
                )
            );
        } else if ((k + 1) * Extractor::bits_per_char <= 128) {
            return std::unique_ptr<IBOSSChunkConstructor>(
                new BOSSChunkConstructor<KmerCounter<typename Extractor::Kmer128, Extractor, uint8_t>>(
                    k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
                )
            );
        } else {
            return std::unique_ptr<IBOSSChunkConstructor>(
                new BOSSChunkConstructor<KmerCounter<typename Extractor::Kmer256, Extractor, uint8_t>>(
                    k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
                )
            );
        }
    } else {
        if ((k + 1) * Extractor::bits_per_char <= 64) {
            return std::unique_ptr<IBOSSChunkConstructor>(
                new BOSSChunkConstructor<KmerCollector<typename Extractor::Kmer64, Extractor>>(
                    k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
                )
            );
        } else if ((k + 1) * Extractor::bits_per_char <= 128) {
            return std::unique_ptr<IBOSSChunkConstructor>(
                new BOSSChunkConstructor<KmerCollector<typename Extractor::Kmer128, Extractor>>(
                    k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
                )
            );
        } else {
            return std::unique_ptr<IBOSSChunkConstructor>(
                new BOSSChunkConstructor<KmerCollector<typename Extractor::Kmer256, Extractor>>(
                    k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
                )
            );
        }
    }
}
