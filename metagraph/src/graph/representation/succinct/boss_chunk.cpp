#include "boss_chunk.hpp"

#include "common/threads/chunked_wait_queue.hpp"
#include "common/circular_buffer.hpp"
#include "common/algorithms.hpp"
#include "common/serialization.hpp"
#include "common/vector.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/utils/template_utils.hpp"

using namespace mg;


template <typename KMER, typename COUNT>
inline const KMER& get_kmer(const std::pair<KMER, COUNT> &pair) {
    return pair.first;
}

template <typename KMER>
inline const KMER& get_kmer(const KMER &kmer) {
    return kmer;
}

static_assert(utils::is_pair<std::pair<KmerExtractorBOSS::Kmer64,uint8_t>>::value);
static_assert(utils::is_pair<std::pair<KmerExtractorBOSS::Kmer128,uint8_t>>::value);
static_assert(utils::is_pair<std::pair<KmerExtractorBOSS::Kmer256,uint8_t>>::value);
static_assert(!utils::is_pair<KmerExtractorBOSS::Kmer64>::value);
static_assert(!utils::is_pair<KmerExtractorBOSS::Kmer128>::value);
static_assert(!utils::is_pair<KmerExtractorBOSS::Kmer256>::value);

//TODO cleanup
// k is node length
template <typename Iterator, typename TAlphabet>
void initialize_chunk(uint64_t alph_size,
                      Iterator begin, Iterator end,
                      size_t k,
                      std::vector<TAlphabet> *W,
                      std::vector<bool> *last,
                      std::vector<uint64_t> *F,
                      sdsl::int_vector<> *weights = nullptr) {
    using T = std::remove_const_t<std::remove_reference_t<decltype(*begin)>>;
    using KMER = std::remove_reference_t<decltype(get_kmer(*begin))>;

    static_assert(KMER::kBitsPerChar <= sizeof(TAlphabet) * 8);

    assert(2 * alph_size - 1 <= std::numeric_limits<TAlphabet>::max());
    assert(alph_size);
    assert(k);
    assert(W && last && F);
    assert(bool(weights) == utils::is_pair<T>::value);

    W->resize(end - begin + 1);
    last->assign(end - begin + 1, 1);
    F->assign(alph_size, 0);

    uint64_t max_count __attribute__((unused)) = 0;
    if (weights) {
        assert(utils::is_pair<T>::value);
        weights->resize(end - begin + 1);
        sdsl::util::set_to_value(*weights, 0);
        max_count = sdsl::bits::lo_set[weights->width()];
    }

    assert(std::is_sorted(begin, end, utils::LessFirst()));

    // the array containing edge labels
    W->at(0) = 0;
    // the bit array indicating last outgoing edges for nodes
    last->at(0) = 0;
    // offsets
    F->at(0) = 0;

    size_t curpos = 1;
    TAlphabet lastF = 0;

    for (Iterator it = begin; it != end; ++it) {
        const KMER &kmer = get_kmer(*it);
        TAlphabet curW = kmer[0];
        TAlphabet curF = kmer[k];

        assert(curW < alph_size);

        // check redundancy and set last
        if (it + 1 < end && KMER::compare_suffix(kmer, get_kmer(*(it + 1)))) {
            // skip redundant dummy sink edges
            if (curW == 0 && curF > 0)
                continue;

            (*last)[curpos] = 0;
        }
        // set W
        if (it != begin && curW > 0) {
            for (Iterator prev = it - 1; KMER::compare_suffix(kmer, get_kmer(*prev), 1);
                 --prev) {
                if (get_kmer(*prev)[0] == curW) {
                    curW += alph_size;
                    break;
                }
                if (prev == begin)
                    break;
            }
        }
        assert(curW < (1llu << (alph_size + 1)));
        (*W)[curpos] = curW;

        while (curF > lastF && lastF + 1 < static_cast<TAlphabet>(alph_size)) {
            F->at(++lastF) = curpos - 1;
        }

        if constexpr(utils::is_pair<T>::value) {
            // set weights for non-dummy k-mers
            if (weights && it->second && kmer[0] && kmer[1])
                (*weights)[curpos] = std::min(static_cast<uint64_t>(it->second), max_count);
        }

        curpos++;
    }
    while (++lastF < alph_size) {
        F->at(lastF) = curpos - 1;
    }

    W->resize(curpos);
    last->resize(curpos);

    if (weights)
        weights->resize(curpos);
}


/**
 * Wrapper class around a static initializer for a BOSS chunk in order to allow
 * partial template specialization for ChunkedWaitQueue, which exposes a more
 * restrictive interface than Vector
 */
template <class Container, class T, typename TAlphabet>
struct Init {
    static void initialize_chunk(uint64_t alph_size,
                                 const Container &container,
                                 size_t k,
                                 std::vector<TAlphabet> *W,
                                 std::vector<bool> *last,
                                 std::vector<uint64_t> *F,
                                 sdsl::int_vector<> *weights = nullptr) {
        return ::initialize_chunk(alph_size, container.begin(), container.end(),
                                  k, W, last, F, weights);
    }
};

// specialization of initialize_chunk for ChunkedWaitQueue
template <typename T, typename TAlphabet>
struct Init<typename common::ChunkedWaitQueue<T>, T, TAlphabet> {
    using Iterator = typename common::ChunkedWaitQueue<T>::iterator;

    template <typename KMER>
    static void set_weight(const size_t count,
                            const KMER &kmer,
                            size_t curpos,
                            uint64_t max_count,
                            sdsl::int_vector<> *weights) {
        if (weights->size() == curpos) {
            weights->resize(weights->size() * 1.5);
        }
        if (count && kmer[0] && kmer[1]) {
            (*weights)[curpos] = std::min(static_cast<uint64_t>(count), max_count);
        } else { // dummy k-mers have a weight of 0
            (*weights)[curpos] = 0;
        }
    }

    static void initialize_chunk(uint64_t alph_size,
                                 const common::ChunkedWaitQueue<T> &container,
                                 size_t k,
                                 std::vector<TAlphabet> *W,
                                 std::vector<bool> *last,
                                 std::vector<uint64_t> *F,
                                 sdsl::int_vector<> *weights = nullptr) {
        Iterator &begin = container.begin();
        Iterator &end = container.end();

        using KMER = std::decay_t<decltype(get_kmer(*begin))>;

        static_assert(KMER::kBitsPerChar <= sizeof(TAlphabet) * 8);

        assert(2 * alph_size - 1 <= std::numeric_limits<TAlphabet>::max());
        assert(alph_size);
        assert(k);
        assert(W && last && F);
        assert(weights == nullptr || utils::is_pair<T>::value);

        F->assign(alph_size, 0);

        uint64_t max_count __attribute__((unused)) = 0;
        if constexpr (utils::is_pair<T>::value) {
            weights->resize(100);
            max_count = sdsl::bits::lo_set[weights->width()];
        }

        W->push_back(0); // the array containing edge labels
        last->push_back(0); // the bit array indicating last outgoing edges for nodes
        F->at(0) = 0; // offsets for the first character

        size_t curpos = 1;
        TAlphabet lastF = 0;
        // last kmer for each label, so we can test multiple edges coming to same node
        std::vector<KMER> last_kmer(alph_size, typename KMER::WordType(0));
        for (Iterator &it = begin; it != end; ++it) {
            const KMER kmer = get_kmer(*it);
            TAlphabet curW = kmer[0];
            TAlphabet curF = kmer[k];

            assert(curW < alph_size);

            // peek at the next entry to check if this is a dummy sink (not source) edge
            // and to set #last
            ++it;
            if (it != end && KMER::compare_suffix(kmer, get_kmer(*it))) {
                // skip redundant dummy sink edges
                if (curW == 0 && curF > 0) {
                    --it;
                    continue;
                }

                last->push_back(false);
            } else {
                last->push_back(true);
            }
            --it;

            if (curW) {
                if (last_kmer[curW].data() && KMER::compare_suffix(kmer, last_kmer[curW], 1)) {
                    assert(last_kmer[curW][0] == curW);
                    // not the first incoming edge to the node, mark with -
                    curW += alph_size;
                } else {
                    last_kmer[curW] = kmer;
                }
            }
            W->push_back(curW);

            while (curF > lastF && lastF + 1 < static_cast<TAlphabet>(alph_size)) {
                F->at(++lastF) = curpos - 1;
            }

            if constexpr (utils::is_pair<T>::value) {
                set_weight((*it).second, kmer, curpos, max_count, weights);
            }

            curpos++;
        }

        while (++lastF < alph_size) {
            F->at(lastF) = curpos - 1;
        }

        W->resize(curpos);
        last->resize(curpos);

        if constexpr (utils::is_pair<T>::value)
            weights->resize(curpos);
    }
};


BOSS::Chunk::Chunk(uint64_t alph_size, size_t k, bool canonical)
      : alph_size_(alph_size), k_(k), canonical_(canonical),
        W_(1, 0), last_(1, 0), F_(alph_size_, 0) {

    assert(sizeof(TAlphabet) * 8 >= get_W_width());
    assert(alph_size_ * 2 <= 1llu << get_W_width());
}

template <typename Array>
BOSS::Chunk::Chunk(uint64_t alph_size, size_t k, bool canonical,
                   const Array &kmers)
      : alph_size_(alph_size), k_(k), canonical_(canonical) {

    assert(sizeof(TAlphabet) * 8 >= get_W_width());
    assert(alph_size_ * 2 <= 1llu << get_W_width());

    Init<Array, typename Array::value_type, TAlphabet>
    ::initialize_chunk(alph_size_, kmers, k_, &W_, &last_, &F_);
}

template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<KmerExtractorBOSS::Kmer64>&);
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<KmerExtractorBOSS::Kmer128>&);
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<KmerExtractorBOSS::Kmer256>&);

template <typename T>
using CWQ = mg::common::ChunkedWaitQueue<T>;
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const CWQ<KmerExtractorBOSS::Kmer64> &);
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const CWQ<KmerExtractorBOSS::Kmer128> &);
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const CWQ<KmerExtractorBOSS::Kmer256> &);

template <typename Array>
BOSS::Chunk::Chunk(uint64_t alph_size,
                   size_t k,
                   bool canonical,
                   const Array &kmers_with_counts,
                   uint8_t bits_per_count)
      : alph_size_(alph_size), k_(k), canonical_(canonical), weights_(0, 0, bits_per_count) {

    assert(sizeof(TAlphabet) * 8 >= get_W_width());
    assert(alph_size_ * 2 <= 1llu << get_W_width());

    Init<Array, typename Array::value_type, TAlphabet>
    ::initialize_chunk(alph_size_, kmers_with_counts, k_,
                       &W_, &last_, &F_, &weights_);
}

template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<std::pair<KmerExtractorBOSS::Kmer64, uint8_t>> &, uint8_t);
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<std::pair<KmerExtractorBOSS::Kmer128, uint8_t>> &, uint8_t);
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<std::pair<KmerExtractorBOSS::Kmer256, uint8_t>> &, uint8_t);

template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<std::pair<KmerExtractorBOSS::Kmer64, uint16_t>> &, uint8_t);
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<std::pair<KmerExtractorBOSS::Kmer128, uint16_t>> &, uint8_t);
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<std::pair<KmerExtractorBOSS::Kmer256, uint16_t>> &, uint8_t);

template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<std::pair<KmerExtractorBOSS::Kmer64, uint32_t>> &, uint8_t);
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<std::pair<KmerExtractorBOSS::Kmer128, uint32_t>> &, uint8_t);
template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<std::pair<KmerExtractorBOSS::Kmer256, uint32_t>> &, uint8_t);

#define BossWithPair(T, C) \
    BOSS::Chunk::Chunk(uint64_t, size_t, bool, const CWQ<std::pair<T, C>> &, uint8_t);

template BossWithPair(KmerExtractorBOSS::Kmer64, uint8_t);
template BossWithPair(KmerExtractorBOSS::Kmer128, uint8_t);
template BossWithPair(KmerExtractorBOSS::Kmer256, uint8_t);

template BossWithPair(KmerExtractorBOSS::Kmer64, uint16_t);
template BossWithPair(KmerExtractorBOSS::Kmer128, uint16_t);
template BossWithPair(KmerExtractorBOSS::Kmer256, uint16_t);

template BossWithPair(KmerExtractorBOSS::Kmer64, uint32_t);
template BossWithPair(KmerExtractorBOSS::Kmer128, uint32_t);
template BossWithPair(KmerExtractorBOSS::Kmer256, uint32_t);


void BOSS::Chunk::push_back(TAlphabet W, TAlphabet F, bool last) {
    assert(W < 2 * alph_size_);
    assert(F < alph_size_);
    assert(k_);

    W_.push_back(W);
    for (TAlphabet a = F + 1; a < F_.size(); ++a) {
        F_[a]++;
    }
    last_.push_back(last);

    assert(weights_.empty());
}

void BOSS::Chunk::extend(const BOSS::Chunk &other) {
    assert(!weights_.size() || weights_.size() == W_.size());
    assert(!other.weights_.size() || other.weights_.size() == other.W_.size());

    if (alph_size_ != other.alph_size_ || k_ != other.k_ || canonical_ != other.canonical_) {
        std::cerr << "ERROR: trying to concatenate incompatible graph chunks" << std::endl;
        exit(1);
    }

    if (!other.size())
        return;

    if (!size()) {
        *this = other;
        return;
    }

    assert(size() && other.size());

    if (weights_.empty() != other.weights_.empty()) {
        std::cerr << "ERROR: trying to concatenate weighted and unweighted blocks" << std::endl;
        exit(1);
    }

    W_.reserve(W_.size() + other.size());
    W_.insert(W_.end(), other.W_.begin() + 1, other.W_.end());

    last_.reserve(last_.size() + other.size());
    last_.insert(last_.end(), other.last_.begin() + 1, other.last_.end());

    assert(F_.size() == other.F_.size());
    for (size_t p = 0; p < other.F_.size(); ++p) {
        F_[p] += other.F_[p];
    }

    if (other.weights_.size()) {
        const auto start = weights_.size();
        weights_.resize(weights_.size() + other.size());
        std::copy(other.weights_.begin() + 1,
                  other.weights_.end(),
                  weights_.begin() + start);
    }

    assert(W_.size() == last_.size());
    assert(!weights_.size() || weights_.size() == W_.size());
}

void BOSS::Chunk::initialize_boss(BOSS *graph, sdsl::int_vector<> *weights) {
    assert(graph->W_);
    delete graph->W_;
    graph->W_ = new wavelet_tree_stat(get_W_width(), std::move(W_));
    W_ = decltype(W_)();

    assert(graph->last_);
    delete graph->last_;
    auto last_bv = to_sdsl(last_);
    last_ = decltype(last_)();
    graph->last_ = new bit_vector_stat(std::move(last_bv));

    graph->F_ = F_;

    graph->k_ = k_;

    graph->state = BOSS::State::STAT;

    if (weights)
        *weights = std::move(weights_);

    assert(graph->is_valid());
}

std::pair<BOSS*, bool>
BOSS::Chunk::build_boss_from_chunks(const std::vector<std::string> &chunk_filenames,
                                    bool verbose,
                                    sdsl::int_vector<> *weights) {
    assert(chunk_filenames.size());

    if (!chunk_filenames.size())
        return std::make_pair(nullptr, false);

    BOSS *graph = new BOSS();
    bool canonical = false;

    uint64_t cumulative_size = 1;

    for (auto file : chunk_filenames) {
        file = utils::remove_suffix(file, kFileExtension) + kFileExtension;

        std::ifstream chunk_in(file, std::ios::binary);

        if (!chunk_in.good()) {
            std::cerr << "ERROR: File corrupted. Cannot load graph chunk "
                      << file << std::endl;
            exit(1);
        }
        cumulative_size += get_number_vector_size(chunk_in) - 1;
    }

    if (verbose)
        std::cout << "Cumulative size of chunks: "
                  << cumulative_size << std::endl;

    sdsl::int_vector<> W;
    sdsl::bit_vector last;
    std::vector<uint64_t> F;

    uint64_t pos = 1;

    for (size_t i = 0; i < chunk_filenames.size(); ++i) {
        auto filename = utils::remove_suffix(chunk_filenames[i], kFileExtension)
                                                        + kFileExtension;
        BOSS::Chunk graph_chunk(1, 0, false);
        if (!graph_chunk.load(filename)) {
            std::cerr << "ERROR: File corrupted. Cannot load graph chunk "
                      << filename << std::endl;
            exit(1);

        } else if (!graph_chunk.k_
                    || !graph_chunk.alph_size_
                    || graph_chunk.last_.size() != graph_chunk.W_.size()
                    || (graph_chunk.weights_.size()
                            && graph_chunk.weights_.size() != graph_chunk.W_.size())) {
            std::cerr << "ERROR: trying to load invalid graph chunk from file "
                      << filename << std::endl;
            exit(1);

        } else if (weights && graph_chunk.weights_.empty()) {
            std::cerr << "ERROR: no weights in graph chunk "
                      << filename << std::endl;
            exit(1);

        } else if (i == 0) {
            W = sdsl::int_vector<>(cumulative_size, 0, graph_chunk.get_W_width());
            last = sdsl::bit_vector(cumulative_size, 0);
            F = std::vector<uint64_t>(graph_chunk.alph_size_, 0);

            graph->k_ = graph_chunk.k_;
            canonical = graph_chunk.canonical_;
            // TODO:
            // graph->alph_size = graph_chunk.alph_size_;

            if (weights) {
                (*weights) = graph_chunk.weights_;
                weights->resize(cumulative_size);
            }

        } else if (graph->k_ != graph_chunk.k_
                    || graph->alph_size != graph_chunk.alph_size_
                    || canonical != graph_chunk.canonical_) {
            std::cerr << "ERROR: trying to concatenate incompatible graph chunks"
                      << std::endl;
            exit(1);

        } else if (weights && weights->width() != graph_chunk.weights_.width()) {
            std::cerr << "ERROR: trying to concatenate chunks with inconsistent weights"
                      << std::endl;
            exit(1);
        }

        if (verbose) {
            std::cout << "Chunk " << filename << " loaded..." << std::flush;
        }

        std::copy(graph_chunk.W_.begin() + 1,
                  graph_chunk.W_.end(),
                  W.begin() + pos);

        std::copy(graph_chunk.last_.begin() + 1,
                  graph_chunk.last_.end(),
                  last.begin() + pos);

        if (weights && i) {
            std::copy(graph_chunk.weights_.begin() + 1,
                      graph_chunk.weights_.end(),
                      weights->begin() + pos);
        }

        pos += graph_chunk.size();

        assert(graph_chunk.F_.size() == F.size());
        for (size_t p = 0; p < F.size(); ++p) {
            F[p] += graph_chunk.F_[p];
        }

        if (verbose) {
            std::cout << " concatenated" << std::endl;
        }
    }

    assert(W.size());
    assert(last.size());
    assert(F.size());

    delete graph->W_;
    graph->W_ = new wavelet_tree_stat(W.width(), std::move(W));
    W = decltype(W)();

    delete graph->last_;
    graph->last_ = new bit_vector_stat(std::move(last));
    last = decltype(last)();

    graph->F_ = std::move(F);

    graph->state = BOSS::State::STAT;

    assert(graph->is_valid());

    return std::make_pair(graph, canonical);
}

bool BOSS::Chunk::load(const std::string &infbase) {
    try {
        std::ifstream instream(utils::remove_suffix(infbase, kFileExtension)
                                                                + kFileExtension,
                               std::ios::binary);

        if (!load_number_vector(instream, &W_)) {
            std::cerr << "ERROR: failed to load W vector" << std::endl;
            return false;
        }

        if (!load_number_vector(instream, &last_)) {
            std::cerr << "ERROR: failed to load L vector" << std::endl;
            return false;
        }

        if (!load_number_vector(instream, &F_)) {
            std::cerr << "ERROR: failed to load F vector" << std::endl;
            return false;
        }

        weights_.load(instream);

        alph_size_ = load_number(instream);
        k_ = load_number(instream);
        canonical_ = load_number(instream);

        return k_ && alph_size_ && W_.size() == last_.size()
                                && F_.size() == alph_size_
                                && (!weights_.size() || weights_.size() == W_.size());

    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load BOSS chunk from "
                  << infbase << "." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void BOSS::Chunk::serialize(const std::string &outbase) const {
    std::ofstream outstream(utils::remove_suffix(outbase, kFileExtension)
                                                                + kFileExtension,
                            std::ios::binary);

    serialize_number_vector(outstream, W_, get_W_width());
    serialize_number_vector(outstream, last_, 1);
    serialize_number_vector(outstream, F_);

    weights_.serialize(outstream);

    serialize_number(outstream, alph_size_);
    serialize_number(outstream, k_);
    serialize_number(outstream, canonical_);
}

uint8_t BOSS::Chunk::get_W_width() const {
    return alph_size_ ? sdsl::bits::hi(alph_size_ * 2 - 1) + 1 : 1;
}
