#include "boss_chunk.hpp"

#include <boost/multiprecision/integer.hpp>

#include "serialization.hpp"
#include "utils.hpp"


template <typename KMER, typename COUNT>
inline const KMER& get_kmer(const Vector<std::pair<KMER, COUNT>> &kmers, uint64_t i) {
    return kmers[i].first;
}

template <typename KMER>
inline const KMER& get_kmer(const Vector<KMER> &kmers, uint64_t i) {
    return kmers[i];
}

static_assert(utils::is_pair<std::pair<KmerExtractor::Kmer64,uint8_t>>::value);
static_assert(utils::is_pair<std::pair<KmerExtractor::Kmer128,uint8_t>>::value);
static_assert(utils::is_pair<std::pair<KmerExtractor::Kmer256,uint8_t>>::value);
static_assert(!utils::is_pair<KmerExtractor::Kmer64>::value);
static_assert(!utils::is_pair<KmerExtractor::Kmer128>::value);
static_assert(!utils::is_pair<KmerExtractor::Kmer256>::value);

//TODO cleanup
// k is node length
template <typename T, typename TAlphabet>
void initialize_chunk(uint64_t alph_size,
                      const Vector<T> &kmers,
                      size_t k,
                      std::vector<TAlphabet> *W,
                      std::vector<bool> *last,
                      std::vector<uint64_t> *F,
                      sdsl::int_vector<> *weights = nullptr) {
    using KMER = std::remove_reference_t<decltype(get_kmer(kmers, 0))>;

    static_assert(KMER::kBitsPerChar <= sizeof(TAlphabet) * 8);

    assert(2 * alph_size - 1 <= std::numeric_limits<TAlphabet>::max());
    assert(alph_size);
    assert(k);
    assert(W && last && F);

    W->resize(kmers.size() + 1);
    last->assign(kmers.size() + 1, 1);
    F->assign(alph_size, 0);

    uint64_t max_count __attribute__((unused)) = 0;
    if (weights) {
        assert(utils::is_pair<T>::value);
        weights->resize(kmers.size() + 1);
        max_count = ~uint64_t(0) >> (64 - weights->width());
    }

#ifndef NDEBUG
    utils::LessFirst<T> less;
    for (uint64_t i = 0; i + 1 < kmers.size(); ++i) {
        assert(less(kmers[i], kmers[i + 1]));
    }
#endif

    // the array containing edge labels
    W->at(0) = 0;
    // the bit array indicating last outgoing edges for nodes
    last->at(0) = 0;
    // offsets
    F->at(0) = 0;

    size_t curpos = 1;
    TAlphabet lastF = 0;

    for (size_t i = 0; i < kmers.size(); ++i) {
        const KMER &kmer = get_kmer(kmers, i);
        TAlphabet curW = kmer[0];
        TAlphabet curF = kmer[k];

        assert(curW < alph_size);

        // check redundancy and set last
        if (i + 1 < kmers.size() && KMER::compare_suffix(kmer, get_kmer(kmers, i + 1))) {
            // skip redundant dummy sink edges
            if (curW == 0 && curF > 0)
                continue;

            (*last)[curpos] = 0;
        }
        //set W
        if (i > 0) {
            for (size_t j = i - 1; KMER::compare_suffix(kmer, get_kmer(kmers, j), 1); --j) {
                if (curW > 0 && get_kmer(kmers, j)[0] == curW) {
                    curW += alph_size;
                    break;
                }
                if (j == 0)
                    break;
            }
        }
        assert(curW < (1llu << (alph_size + 1)));
        (*W)[curpos] = curW;

        while (curF > lastF && lastF + 1 < static_cast<TAlphabet>(alph_size)) {
            F->at(++lastF) = curpos - 1;
        }

        if constexpr(utils::is_pair<T>::value) {
            if (weights)
                (*weights)[curpos] = std::min(static_cast<uint64_t>(kmers[i].second), max_count);
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


BOSS::Chunk::Chunk(uint64_t alph_size, size_t k)
      : alph_size_(alph_size), k_(k), W_(1, 0), last_(1, 0), F_(alph_size_, 0) {

    assert(sizeof(TAlphabet) * 8 >= extended_alph_size());
    assert(alph_size_ * 2 <= 1llu << extended_alph_size());
}

template <typename KMER>
BOSS::Chunk::Chunk(uint64_t alph_size, size_t k, const Vector<KMER> &kmers)
      : alph_size_(alph_size), k_(k) {

    assert(sizeof(TAlphabet) * 8 >= extended_alph_size());
    assert(alph_size_ * 2 <= 1llu << extended_alph_size());

    initialize_chunk(alph_size_, kmers, k_, &W_, &last_, &F_);
}

template BOSS::Chunk::Chunk(uint64_t, size_t, const Vector<KmerExtractor::Kmer64>&);
template BOSS::Chunk::Chunk(uint64_t, size_t, const Vector<KmerExtractor::Kmer128>&);
template BOSS::Chunk::Chunk(uint64_t, size_t, const Vector<KmerExtractor::Kmer256>&);

template <typename KMER, typename COUNT>
BOSS::Chunk::Chunk(uint64_t alph_size,
                   size_t k,
                   const Vector<std::pair<KMER, COUNT>> &kmers,
                   uint8_t bits_per_count)
      : alph_size_(alph_size), k_(k), weights_(0, 0, bits_per_count) {

    assert(sizeof(TAlphabet) * 8 >= extended_alph_size());
    assert(alph_size_ * 2 <= 1llu << extended_alph_size());

    initialize_chunk(alph_size_, kmers, k_, &W_, &last_, &F_, &weights_);
}

template BOSS::Chunk::Chunk(uint64_t, size_t, const Vector<std::pair<KmerExtractor::Kmer64, uint8_t>> &, uint8_t);
template BOSS::Chunk::Chunk(uint64_t, size_t, const Vector<std::pair<KmerExtractor::Kmer128, uint8_t>> &, uint8_t);
template BOSS::Chunk::Chunk(uint64_t, size_t, const Vector<std::pair<KmerExtractor::Kmer256, uint8_t>> &, uint8_t);

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

    if (alph_size_ != other.alph_size_ || k_ != other.k_) {
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

    if (weights_.empty() ^ other.weights_.empty()) {
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
    graph->W_ = new wavelet_tree_stat(extended_alph_size(), std::move(W_));
    W_ = decltype(W_)();

    assert(graph->last_);
    delete graph->last_;
    auto last_bv = to_sdsl(last_);
    last_ = decltype(last_)();
    graph->last_ = new bit_vector_stat(std::move(last_bv));

    graph->F_ = F_;

    graph->k_ = k_;

    graph->state = Config::STAT;

    if (weights)
        *weights = std::move(weights_);

    assert(graph->is_valid());
}

BOSS*
BOSS::Chunk::build_boss_from_chunks(const std::vector<std::string> &chunk_filenames,
                                    bool verbose,
                                    sdsl::int_vector<> *weights) {
    assert(chunk_filenames.size());

    BOSS *graph = new BOSS();
    if (!chunk_filenames.size())
        return graph;

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
        BOSS::Chunk graph_chunk(1, 0);
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
            W = sdsl::int_vector<>(cumulative_size, 0, graph_chunk.extended_alph_size());
            last = sdsl::bit_vector(cumulative_size, 0);
            F = std::vector<uint64_t>(graph_chunk.alph_size_, 0);

            graph->k_ = graph_chunk.k_;
            // TODO:
            // graph->alph_size = graph_chunk.alph_size_;

            if (weights) {
                (*weights) = graph_chunk.weights_;
                weights->resize(cumulative_size);
            }

        } else if (graph->k_ != graph_chunk.k_
                    || graph->alph_size != graph_chunk.alph_size_) {
            std::cerr << "ERROR: trying to concatenate incompatible graph chunks"
                      << std::endl;
            exit(1);

        } else if (weights && weights->width() != graph_chunk.weights_.width()) {
            std::cerr << "ERROR: trying to concatenate chunks with incostistent weights"
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

    graph->state = Config::STAT;

    assert(graph->is_valid());

    return graph;
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

    serialize_number_vector(outstream, W_, extended_alph_size());
    serialize_number_vector(outstream, last_, 1);
    serialize_number_vector(outstream, F_);

    weights_.serialize(outstream);

    serialize_number(outstream, alph_size_);
    serialize_number(outstream, k_);
}

uint8_t BOSS::Chunk::extended_alph_size() const {
    return utils::code_length(alph_size_) + 1;
}
