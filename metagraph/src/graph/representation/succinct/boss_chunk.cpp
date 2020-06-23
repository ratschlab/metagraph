#include "boss_chunk.hpp"

#include "common/threads/chunked_wait_queue.hpp"
#include "common/serialization.hpp"
#include "common/vector.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/utils/template_utils.hpp"

using namespace mtg;
using utils::get_first;
using mtg::kmer::KmerExtractorBOSS;

static_assert(utils::is_pair_v<std::pair<KmerExtractorBOSS::Kmer64, uint8_t>>);
static_assert(utils::is_pair_v<std::pair<KmerExtractorBOSS::Kmer128, uint8_t>>);
static_assert(utils::is_pair_v<std::pair<KmerExtractorBOSS::Kmer256, uint8_t>>);
static_assert(!utils::is_pair_v<KmerExtractorBOSS::Kmer64>);
static_assert(!utils::is_pair_v<KmerExtractorBOSS::Kmer128>);
static_assert(!utils::is_pair_v<KmerExtractorBOSS::Kmer256>);

const double kGrowthFactor = 1.5;


// k is node length
template <typename Iterator>
void initialize_chunk(uint64_t alph_size,
                      Iterator *begin_ptr, Iterator *end_ptr,
                      size_t k,
                      sdsl::int_vector<> *W,
                      sdsl::bit_vector *last,
                      std::vector<uint64_t> *F,
                      sdsl::int_vector<> *weights = nullptr) {
    Iterator &it = *begin_ptr;
    Iterator &end = *end_ptr;

    using T = std::decay_t<decltype(*it)>;
    using KMER = utils::get_first_type_t<T>;
    using CharType = typename KMER::CharType;

    assert(KMER::kBitsPerChar <= W->width());
    assert(2 * alph_size - 1 <= sdsl::bits::lo_set[W->width()]);
    assert(alph_size);
    assert(k);
    assert(W && last && F);
    assert(bool(weights) == utils::is_pair_v<T>);

    uint64_t max_count __attribute__((unused)) = 0;

    if (W->size() < 1000)
        W->resize(1000);
    last->resize(W->size());
    if constexpr(utils::is_pair_v<T>) {
        weights->resize(W->size());
        (*weights)[0] = 0;
        max_count = sdsl::bits::lo_set[weights->width()];
    }
    (*W)[0] = 0; // the array containing edge labels
    (*last)[0] = 0; // the bit array indicating last outgoing edges for nodes
    F->assign(alph_size, 0); // the offsets for the last characters in the nodes of BOSS

    size_t curpos = 1;
    CharType lastF = 0;
    // last kmer for each label, so we can test multiple edges coming to same node
    std::vector<KMER> last_kmer(alph_size, typename KMER::WordType(0));

    while (it != end) {
        const T value = *it;
        const KMER &kmer = get_first(value);

        uint64_t curW = kmer[0];
        CharType curF = kmer[k];

        assert(curW < alph_size);

        assert(last->size() == W->size());
        assert(!utils::is_pair_v<T> || weights->size() == W->size());

        if (curpos == W->size()) {
            W->resize(curpos * kGrowthFactor);
            last->resize(curpos * kGrowthFactor);
            if constexpr(utils::is_pair_v<T>)
                weights->resize(curpos * kGrowthFactor);
        }

        // peek at the next entry to check if this is a dummy sink (not source) edge
        // and to set #last
        ++it;
        if (it != end && KMER::compare_suffix(kmer, get_first(*it))) {
            // skip redundant dummy sink edges
            if (curW == 0 && curF > 0)
                continue;

            (*last)[curpos] = false;
        } else {
            (*last)[curpos] = true;
        }

        // set W
        if (curW) {
            if (last_kmer[curW].data()
                    && KMER::compare_suffix(kmer, last_kmer[curW], 1)) {
                assert(last_kmer[curW][0] == curW);
                // not the first incoming edge to the node, mark with -
                curW += alph_size;
            } else {
                last_kmer[curW] = kmer;
            }
        }
        assert(curW <= sdsl::bits::lo_set[W->width()]);
        (*W)[curpos] = curW;

        while (curF > lastF && lastF + 1 < alph_size) {
            F->at(++lastF) = curpos - 1;
        }

        if constexpr(utils::is_pair_v<T>) {
            if (value.second && curW && kmer[1]) {
                (*weights)[curpos] = std::min(static_cast<uint64_t>(value.second),
                                              max_count);
            } else {
                (*weights)[curpos] = 0; // dummy k-mers have a weight of 0
            }
        }

        curpos++;
    }

    while (++lastF < alph_size) {
        F->at(lastF) = curpos - 1;
    }

    W->resize(curpos);
    last->resize(curpos);
    if constexpr(utils::is_pair_v<T>)
        weights->resize(curpos);
}


BOSS::Chunk::Chunk(uint64_t alph_size, size_t k, bool canonical)
      : alph_size_(alph_size), k_(k), canonical_(canonical),
        W_(1, 0, get_W_width()), last_(1, 0), F_(alph_size_, 0), size_(1) {}


template <typename Array>
BOSS::Chunk::Chunk(uint64_t alph_size,
                   size_t k,
                   bool canonical,
                   const Array &kmers_with_counts,
                   uint8_t bits_per_count)
      : Chunk(alph_size, k, canonical) {

    weights_.width(bits_per_count);

    if constexpr(utils::is_instance_v<Array, common::ChunkedWaitQueue>) {
        initialize_chunk(alph_size_,
                         &kmers_with_counts.begin(),
                         &kmers_with_counts.end(),
                         k_, &W_, &last_, &F_, bits_per_count ? &weights_ : NULL);
    } else {
        auto begin = kmers_with_counts.begin();
        auto end = kmers_with_counts.end();
        W_.resize(end - begin + 1);
        last_.resize(end - begin + 1);
        if (bits_per_count)
            weights_.resize(end - begin + 1);
        initialize_chunk(alph_size_, &begin, &end,
                         k_, &W_, &last_, &F_, bits_per_count ? &weights_ : NULL);
    }
    size_ = W_.size();
}

template <typename T>
using CWQ = common::ChunkedWaitQueue<T>;

#define INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(...) \
    template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const CWQ<__VA_ARGS__> &, uint8_t); \
    template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<__VA_ARGS__> &, uint8_t);

INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(KmerExtractorBOSS::Kmer64);
INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(KmerExtractorBOSS::Kmer128);
INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(KmerExtractorBOSS::Kmer256);

INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(std::pair<KmerExtractorBOSS::Kmer64, uint8_t>);
INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(std::pair<KmerExtractorBOSS::Kmer128, uint8_t>);
INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(std::pair<KmerExtractorBOSS::Kmer256, uint8_t>);

INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(std::pair<KmerExtractorBOSS::Kmer64, uint16_t>);
INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(std::pair<KmerExtractorBOSS::Kmer128, uint16_t>);
INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(std::pair<KmerExtractorBOSS::Kmer256, uint16_t>);

INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(std::pair<KmerExtractorBOSS::Kmer64, uint32_t>);
INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(std::pair<KmerExtractorBOSS::Kmer128, uint32_t>);
INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(std::pair<KmerExtractorBOSS::Kmer256, uint32_t>);


void BOSS::Chunk::push_back(TAlphabet W, TAlphabet F, bool last) {
    assert(W < 2 * alph_size_);
    assert(F < alph_size_);
    assert(k_);

    assert(last_.size() == W_.size());
    assert(weights_.empty());

    if (size_ == W_.size()) {
        W_.resize(std::max(size_ + 1, static_cast<uint64_t>(size_ * kGrowthFactor)));
        last_.resize(std::max(size_ + 1, static_cast<uint64_t>(size_ * kGrowthFactor)));
    }

    W_[size_] = W;
    last_[size_] = last;
    for (TAlphabet a = F + 1; a < F_.size(); ++a) {
        F_[a]++;
    }
    size_++;
}

void BOSS::Chunk::extend(const BOSS::Chunk &other) {
    assert(size_ && other.size_);
    assert(!weights_.size() || weights_.size() == W_.size());
    assert(!other.weights_.size() || other.weights_.size() == other.W_.size());

    if (alph_size_ != other.alph_size_
            || k_ != other.k_
            || canonical_ != other.canonical_
            || W_.width() != other.W_.width()) {
        std::cerr << "ERROR: trying to concatenate incompatible graph chunks" << std::endl;
        exit(1);
    }

    if (other.size_ == 1)
        return;

    if (size_ == 1) {
        *this = other;
        return;
    }

    if (weights_.empty() != other.weights_.empty()) {
        std::cerr << "ERROR: trying to concatenate weighted and unweighted blocks" << std::endl;
        exit(1);
    }

    W_.resize(size_ + other.size_ - 1);
    std::copy(other.W_.begin() + 1,
              other.W_.begin() + other.size_,
              W_.begin() + size_);

    last_.resize(size_ + other.size_ - 1);
    std::copy(other.last_.begin() + 1,
              other.last_.begin() + other.size_,
              last_.begin() + size_);

    assert(F_.size() == other.F_.size());
    for (size_t p = 0; p < other.F_.size(); ++p) {
        F_[p] += other.F_[p];
    }

    if (other.weights_.size()) {
        weights_.resize(size_ + other.size_ - 1);
        std::copy(other.weights_.begin() + 1,
                  other.weights_.begin() + other.size_,
                  weights_.begin() + size_);
    }

    size_ += other.size_ - 1;

    assert(W_.size() == last_.size());
    assert(!weights_.size() || weights_.size() == W_.size());
}

void BOSS::Chunk::initialize_boss(BOSS *graph, sdsl::int_vector<> *weights) {
    assert(size_ <= W_.size());
    assert(last_.size() == W_.size());
    assert(weights_.empty() || weights_.size() == W_.size());

    W_.resize(size_);
    last_.resize(size_);
    if (weights_.size())
        weights_.resize(size_);

    assert(graph->W_);
    delete graph->W_;
    graph->W_ = new wavelet_tree_small(get_W_width(), std::move(W_));
    W_ = decltype(W_)();

    assert(graph->last_);
    delete graph->last_;
    graph->last_ = new bit_vector_stat(std::move(last_));
    last_ = decltype(last_)();

    graph->F_ = F_;
    graph->recompute_NF();

    graph->k_ = k_;

    graph->state = BOSS::State::SMALL;

    if (weights) {
        weights->swap(weights_);
    }
    weights_ = decltype(weights_)();

    size_ = 0;

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
        cumulative_size += load_number(chunk_in) - 1;
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
                    || canonical != graph_chunk.canonical_
                    || W.width() != graph_chunk.W_.width()) {
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
                  graph_chunk.W_.begin() + graph_chunk.size_,
                  W.begin() + pos);

        std::copy(graph_chunk.last_.begin() + 1,
                  graph_chunk.last_.begin() + graph_chunk.size_,
                  last.begin() + pos);

        if (weights && i) {
            std::copy(graph_chunk.weights_.begin() + 1,
                      graph_chunk.weights_.begin() + graph_chunk.size_,
                      weights->begin() + pos);
        }

        pos += graph_chunk.size_ - 1;

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
    graph->W_ = new wavelet_tree_small(W.width(), std::move(W));
    W = decltype(W)();

    delete graph->last_;
    graph->last_ = new bit_vector_stat(std::move(last));
    last = decltype(last)();

    graph->F_ = std::move(F);
    graph->recompute_NF();

    graph->state = BOSS::State::SMALL;

    assert(graph->is_valid());

    return std::make_pair(graph, canonical);
}

bool BOSS::Chunk::load(const std::string &infbase) {
    try {
        std::ifstream instream(utils::remove_suffix(infbase, kFileExtension)
                                                                + kFileExtension,
                               std::ios::binary);
        size_ = load_number(instream);
        W_.load(instream);
        last_.load(instream);

        if (!load_number_vector(instream, &F_)) {
            std::cerr << "ERROR: failed to load F vector" << std::endl;
            return false;
        }

        weights_.load(instream);

        alph_size_ = load_number(instream);
        k_ = load_number(instream);
        canonical_ = load_number(instream);

        return k_ && alph_size_ && W_.size() == last_.size()
                                && size_ <= W_.size()
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
    serialize_number(outstream, size_);
    W_.serialize(outstream);
    last_.serialize(outstream);
    serialize_number_vector(outstream, F_);

    weights_.serialize(outstream);

    serialize_number(outstream, alph_size_);
    serialize_number(outstream, k_);
    serialize_number(outstream, canonical_);
}

uint8_t BOSS::Chunk::get_W_width() const {
    return alph_size_ ? sdsl::bits::hi(alph_size_ * 2 - 1) + 1 : 1;
}
