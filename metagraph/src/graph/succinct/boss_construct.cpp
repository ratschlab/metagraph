#include "boss_construct.hpp"

#include <ips4o.hpp>

#include "unix_tools.hpp"
#include "boss_chunk.hpp"
#include "kmer_collector.hpp"


template <typename KMER>
void erase_redundant_dummy_kmers(Vector<KMER> *kmers) {

    assert(std::is_sorted(kmers->begin(), kmers->end()));

    size_t cur_pos = 0;

    KmerExtractor::TAlphabet edge_label, node_last_char;

    for (size_t i = 0; i < kmers->size(); ++i) {
        const KMER &kmer = kmers->at(i);
        node_last_char = kmer[1];
        edge_label = kmer[0];

        // check if the k-mer isn't dummy
        if (edge_label && node_last_char) {
            kmers->at(cur_pos++) = kmer;
            continue;
        }

        // skip redundant sink dummy kmers
        if (node_last_char && !edge_label
                           && i + 1 < kmers->size()
                           && KMER::compare_suffix(kmer, kmers->at(i + 1), 0)) {
            continue;
        }

        // check if it isn't a source dummy kmer
        if (!edge_label || node_last_char) {
            kmers->at(cur_pos++) = kmer;
            continue;
        }

        bool redundant = false;
        for (size_t j = i + 1; j < kmers->size()
                                && KMER::compare_suffix(kmer, kmers->at(j), 1); ++j) {
            if (edge_label == kmers->at(j)[0]) {
                // This source dummy kmer is redundant and has to be erased
                redundant = true;
                break;
            }
        }

        // keep the dummy kmer in the list if it's not redundant
        if (!redundant)
            kmers->at(cur_pos++) = kmer;
    }

    kmers->resize(cur_pos);
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
    erase_redundant_dummy_kmers(kmers);

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
                        verbose) {
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
                              const KMER *kmers,
                              uint64_t num_kmers) {
    assert(std::is_sorted(kmers, kmers + num_kmers));

    // the array containing edge labels
    std::vector<KmerExtractor::TAlphabet> W(1 + num_kmers);
    W[0] = 0;
    // the bit array indicating last outgoing edges for nodes
    std::vector<bool> last(1 + num_kmers, 1);
    last[0] = 0;
    // offsets
    std::vector<uint64_t> F(alph_size, 0);
    F.at(0) = 0;

    size_t curpos = 1;
    KmerExtractor::TAlphabet lastF = 0;

    for (size_t i = 0; i < num_kmers; ++i) {
        KmerExtractor::TAlphabet curW = kmers[i][0];
        KmerExtractor::TAlphabet curF = kmers[i][k];

        assert(curW < alph_size);

        // check redundancy and set last
        if (i + 1 < num_kmers && KMER::compare_suffix(kmers[i], kmers[i + 1])) {
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
BOSS::Chunk* BOSSChunkConstructor<KMER>
::build_chunk() {
    auto &kmers = kmer_collector_.data();

    if (!kmer_collector_.suffix_length()) {
        if (kmer_collector_.verbose()) {
            std::cout << "Reconstructing all required dummy source k-mers...\t"
                      << std::flush;
        }
        Timer timer;

        // kmer_collector stores (BOSS::k_ + 1)-mers
        recover_source_dummy_nodes(kmer_collector_.get_k() - 1,
                                   &kmers,
                                   kmer_collector_.num_threads(),
                                   kmer_collector_.verbose());

        if (kmer_collector_.verbose())
            std::cout << timer.elapsed() << "sec" << std::endl;
    }

    // kmer_collector stores (BOSS::k_ + 1)-mers
    BOSS::Chunk *result = chunk_from_kmers(
        kmer_collector_.alphabet_size(),
        kmer_collector_.get_k() - 1,
        kmers.data(),
        kmers.size()
    );

    kmer_collector_.clear();

    return result;
}


BOSSConstructor::BOSSConstructor(size_t k,
                                 bool canonical_mode,
                                 const std::string &filter_suffix,
                                 size_t num_threads,
                                 double memory_preallocated,
                                 bool verbose)
      : constructor_(IBOSSChunkConstructor::initialize(
            k,
            canonical_mode,
            filter_suffix,
            num_threads,
            memory_preallocated,
            verbose)
      ) {}

void BOSSConstructor::build_graph(BOSS *graph) {
    auto chunk = constructor_->build_chunk();
    // initialize graph from the chunk built
    chunk->initialize_boss(graph);
    delete chunk;
}

BOSS* BOSSConstructor
::build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                          bool verbose) {
    // TODO: move from chunk to here?
    return BOSS::Chunk::build_boss_from_chunks(chunk_filenames, verbose);
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
