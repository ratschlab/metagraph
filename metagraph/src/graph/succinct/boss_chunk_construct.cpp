#include "boss_chunk_construct.hpp"

#include <ips4o.hpp>

#include "unix_tools.hpp"
#include "boss_chunk.hpp"
#include "kmer_collector.hpp"


template <typename KMER>
void dont_erase_redundant_dummy_kmers(Vector<KMER> *) {}

template <typename KMER>
void erase_redundant_dummy_kmers(Vector<KMER> *kmers) {

    assert(std::is_sorted(kmers->begin(), kmers->end()));
    assert(std::unique(kmers->begin(), kmers->end()) == kmers->end());

    if (kmers->size() < 2)
        return;

    // last k-mer is never redundant. Start with the next one.
    uint64_t last = kmers->size() - 1;

    KmerExtractor::TAlphabet edge_label, node_last_char;

    std::vector<uint64_t> last_kmer(1llu << KMER::kBitsPerChar, kmers->size());

    last_kmer[kmers->at(last)[0]] = last;

    for (int64_t i = last - 1; i >= 0; --i) {
        const KMER &kmer = kmers->at(i);
        node_last_char = kmer[1];
        edge_label = kmer[0];

        // assert((edge_label || node_last_char)
        //             && "dummy k-mer cannot be both source and sink dummy");

        if (!edge_label) {
            // sink dummy k-mer

            // skip if redundant
            if (node_last_char && KMER::compare_suffix(kmer, kmers->at(last), 0))
                continue;

        } else if (!node_last_char) {
            // source dummy k-mer

            // skip if redundant
            if (last_kmer[edge_label] < kmers->size()
                    && KMER::compare_suffix(kmer, kmers->at(last_kmer[edge_label]), 1))
                continue;
        }

        // the k-mer is either not dummy, or not redundant -> keep the k-mer
        kmers->at(--last) = kmer;
        last_kmer[edge_label] = last;
    }

    kmers->erase(kmers->begin(), kmers->begin() + last);
}


template <class V>
void sort_and_remove_duplicates(V *array,
                                size_t num_threads,
                                size_t offset) {
    ips4o::parallel::sort(array->begin() + offset, array->end(),
                          std::less<typename V::value_type>(),
                          num_threads);
    // remove duplicates
    auto unique_end = std::unique(array->begin() + offset, array->end());
    array->erase(unique_end, array->end());
}

template <typename KMER>
void shrink_kmers(Vector<KMER> *kmers,
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
                  << (kmers->size() * sizeof(KMER) >> 20) << "Mb" << std::endl;
    }
}

// Although this function could be parallelized better,
// the experiments show it's already fast enough.
// k is node length
template <typename KMER>
void recover_source_dummy_nodes(size_t k,
                                Vector<KMER> *kmers,
                                size_t num_threads,
                                bool verbose) {
    size_t dummy_begin = kmers->size();
    size_t num_dummy_parent_kmers = 0;

    for (size_t i = 0; i < dummy_begin; ++i) {
        const KMER &kmer = kmers->at(i);
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

        kmers->push_back(kmers->at(i));
        kmers->back().to_prev(k + 1, BOSS::kSentinelCode);
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

            kmers->push_back(kmers->at(i));
            kmers->back().to_prev(k + 1, BOSS::kSentinelCode);
        }
        sort_and_remove_duplicates(kmers, num_threads, dummy_begin);

        if (verbose) {
            std::cout << "Number of dummy k-mers with dummy prefix of length " << c
                      << ": " << kmers->size() - dummy_begin << std::endl;
        }
    }
    ips4o::parallel::sort(kmers->begin(), kmers->end(), std::less<KMER>(), num_threads);
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


template <typename KMER>
class BOSSChunkConstructor : public IBOSSChunkConstructor {
    friend IBOSSChunkConstructor;

  private:
    BOSSChunkConstructor(size_t k,
                         bool canonical_mode = false,
                         const std::string &filter_suffix = "",
                         size_t num_threads = 1,
                         double memory_preallocated = 0,
                         bool verbose = false);

    void add_sequence(std::string&& sequence) {
        kmer_collector_.add_sequence(std::move(sequence));
    }

    void add_sequences(std::function<void(CallString)> generate_sequences) {
        kmer_collector_.add_sequences(generate_sequences);
    }

    BOSS::Chunk* build_chunk();

    KmerCollector<KMER, KmerExtractor> kmer_collector_;
};

template <typename KMER>
BOSSChunkConstructor<KMER>
::BOSSChunkConstructor(size_t k,
                       bool canonical_mode,
                       const std::string &filter_suffix,
                       size_t num_threads,
                       double memory_preallocated,
                       bool verbose)
      : kmer_collector_(k + 1,
                        canonical_mode,
                        encode_filter_suffix_boss<KmerExtractor>(filter_suffix),
                        num_threads,
                        memory_preallocated,
                        verbose,
                        filter_suffix.empty() ? erase_redundant_dummy_kmers<KMER>
                                              : dont_erase_redundant_dummy_kmers<KMER>) {
    if (filter_suffix == std::string(filter_suffix.size(), BOSS::kSentinel)) {
        kmer_collector_.data().emplace_back(
            std::vector<KmerExtractor::TAlphabet>(k + 1, BOSS::kSentinelCode)
        );
    }
}

//TODO cleanup
// k is node length
template <typename KMER>
BOSS::Chunk* chunk_from_kmers(KmerExtractor::TAlphabet alph_size,
                              size_t k,
                              const Vector<KMER> &kmers) {
    assert(std::is_sorted(kmers.begin(), kmers.end()));

    // the array containing edge labels
    std::vector<KmerExtractor::TAlphabet> W(1 + kmers.size());
    W[0] = 0;
    // the bit array indicating last outgoing edges for nodes
    std::vector<bool> last(1 + kmers.size(), 1);
    last[0] = 0;
    // offsets
    std::vector<uint64_t> F(alph_size, 0);
    F.at(0) = 0;

    size_t curpos = 1;
    KmerExtractor::TAlphabet lastF = 0;

    for (size_t i = 0; i < kmers.size(); ++i) {
        KmerExtractor::TAlphabet curW = kmers[i][0];
        KmerExtractor::TAlphabet curF = kmers[i][k];

        assert(curW < alph_size);

        // check redundancy and set last
        if (i + 1 < kmers.size() && KMER::compare_suffix(kmers[i], kmers[i + 1])) {
            // skip redundant dummy edges
            if (curW == 0 && curF > 0)
                continue;

            last[curpos] = 0;
        }
        //set W
        if (i > 0) {
            for (size_t j = i - 1; KMER::compare_suffix(kmers[i], kmers[j], 1); --j) {
                if (curW > 0 && kmers[j][0] == curW) {
                    curW += alph_size;
                    break;
                }
                if (j == 0)
                    break;
            }
        }
        W[curpos] = curW;

        while (lastF + 1 < alph_size && curF != lastF) {
            F.at(++lastF) = curpos - 1;
        }
        curpos++;
    }
    while (++lastF < alph_size) {
        F.at(lastF) = curpos - 1;
    }

    W.resize(curpos);
    last.resize(curpos);

    return new BOSS::Chunk(k, std::move(W), std::move(last), std::move(F));
}

template <typename KMER>
BOSS::Chunk* BOSSChunkConstructor<KMER>::build_chunk() {
    auto &kmers = kmer_collector_.data();

    if (!kmer_collector_.suffix_length()) {
        if (kmer_collector_.verbose()) {
            std::cout << "Reconstructing all required dummy source k-mers..."
                      << std::endl;
        }
        Timer timer;

        // kmer_collector stores (BOSS::k_ + 1)-mers
        recover_source_dummy_nodes(kmer_collector_.get_k() - 1,
                                   &kmers,
                                   kmer_collector_.num_threads(),
                                   kmer_collector_.verbose());

        if (kmer_collector_.verbose())
            std::cout << "Dummy source k-mers were reconstructed in "
                      << timer.elapsed() << "sec" << std::endl;
    }

    // kmer_collector stores (BOSS::k_ + 1)-mers
    BOSS::Chunk *result = chunk_from_kmers(kmer_collector_.alphabet_size(),
                                           kmer_collector_.get_k() - 1,
                                           kmers);
    kmer_collector_.clear();

    return result;
}

IBOSSChunkConstructor*
IBOSSChunkConstructor
::initialize(size_t k,
             bool canonical_mode,
             const std::string &filter_suffix,
             size_t num_threads,
             double memory_preallocated,
             bool verbose) {
    using Extractor = KmerExtractor;

    if ((k + 1) * Extractor::kLogSigma <= 64) {
        return new BOSSChunkConstructor<typename Extractor::Kmer64>(
            k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
        );
    } else if ((k + 1) * Extractor::kLogSigma <= 128) {
        return new BOSSChunkConstructor<typename Extractor::Kmer128>(
            k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
        );
    } else {
        return new BOSSChunkConstructor<typename Extractor::Kmer256>(
            k, canonical_mode, filter_suffix, num_threads, memory_preallocated, verbose
        );
    }
}
