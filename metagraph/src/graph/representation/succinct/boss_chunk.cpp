#include "boss_chunk.hpp"

#include "common/threads/chunked_wait_queue.hpp"
#include "common/serialization.hpp"
#include "common/vector.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/utils/template_utils.hpp"
#include "common/utils/file_utils.hpp"

using namespace mtg;
using utils::get_first;
using mtg::kmer::KmerExtractorBOSS;

const uint64_t BUFFER_SIZE = 5 * 1024 * 1024;

static_assert(utils::is_pair_v<std::pair<KmerExtractorBOSS::Kmer64, uint8_t>>);
static_assert(utils::is_pair_v<std::pair<KmerExtractorBOSS::Kmer128, uint8_t>>);
static_assert(utils::is_pair_v<std::pair<KmerExtractorBOSS::Kmer256, uint8_t>>);
static_assert(!utils::is_pair_v<KmerExtractorBOSS::Kmer64>);
static_assert(!utils::is_pair_v<KmerExtractorBOSS::Kmer128>);
static_assert(!utils::is_pair_v<KmerExtractorBOSS::Kmer256>);


// k is node length
template <typename Iterator>
void initialize_chunk(uint64_t alph_size,
                      Iterator &it, const Iterator &end,
                      size_t k,
                      sdsl::int_vector_buffer<> *W,
                      sdsl::int_vector_buffer<1> *last,
                      std::vector<uint64_t> *F,
                      sdsl::int_vector_buffer<> *weights = nullptr) {
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

    W->reset();
    last->reset();
    if constexpr(utils::is_pair_v<T>) {
        weights->reset();
        weights->push_back(0);
        max_count = sdsl::bits::lo_set[weights->width()];
    }
    W->push_back(0); // the array containing edge labels
    last->push_back(0); // the bit array indicating last outgoing edges for nodes
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

        // peek at the next entry to check if this is a dummy sink (not source) edge
        // and to set #last
        ++it;
        if (it != end && KMER::compare_suffix(kmer, get_first(*it))) {
            // skip redundant dummy sink edges
            if (curW == 0 && curF > 0)
                continue;

            last->push_back(false);
        } else {
            last->push_back(true);
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
        W->push_back(curW);

        while (curF > lastF && lastF + 1 < alph_size) {
            F->at(++lastF) = curpos - 1;
        }

        if constexpr(utils::is_pair_v<T>) {
            if (value.second && curW && kmer[1]) {
                weights->push_back(std::min(static_cast<uint64_t>(value.second),
                                            max_count));
            } else {
                weights->push_back(0); // dummy k-mers have a weight of 0
            }
        }

        curpos++;
    }

    while (++lastF < alph_size) {
        F->at(lastF) = curpos - 1;
    }

    assert(W->size() == curpos);
    assert(last->size() == curpos);
    assert(!weights || weights->size() == curpos);
}

const std::string& create_int_vector(const std::string &filename) {
    sdsl::int_vector<> int_vector;
    std::ofstream out(filename, std::ios::binary);
    int_vector.serialize(out);
    return filename;
}

BOSS::Chunk::Chunk(uint64_t alph_size, size_t k, bool canonical,
                   const std::string &swap_dir)
      : alph_size_(alph_size), k_(k), canonical_(canonical),
        dir_(utils::create_temp_dir(swap_dir, "graph_chunk")),
        W_(create_int_vector(dir_ + "/W"), std::ios::in | std::ios::out, BUFFER_SIZE, get_W_width()),
        last_(create_int_vector(dir_ + "/last"), std::ios::in | std::ios::out, BUFFER_SIZE) {
    W_.push_back(0);
    last_.push_back(0);
    F_.assign(alph_size_, 0);
    size_ = 1;
}

BOSS::Chunk::~Chunk() {
    std::filesystem::remove_all(dir_);
}

template <typename Array>
BOSS::Chunk::Chunk(uint64_t alph_size,
                   size_t k,
                   bool canonical,
                   const Array &kmers_with_counts,
                   uint8_t bits_per_count,
                   const std::string &swap_dir)
      : Chunk(alph_size, k, canonical, swap_dir) {

    weights_ = sdsl::int_vector_buffer<>(create_int_vector(dir_ + "/weights"),
                                         std::ios::in | std::ios::out,
                                         BUFFER_SIZE,
                                         bits_per_count);

    if constexpr(utils::is_instance_v<Array, common::ChunkedWaitQueue>) {
        initialize_chunk(alph_size_,
                         kmers_with_counts.begin(),
                         kmers_with_counts.end(),
                         k_, &W_, &last_, &F_, bits_per_count ? &weights_ : NULL);
    } else {
        auto begin = kmers_with_counts.begin();
        auto end = kmers_with_counts.end();
        initialize_chunk(alph_size_, begin, end,
                         k_, &W_, &last_, &F_, bits_per_count ? &weights_ : NULL);
    }
    size_ = W_.size();
}

template <typename T>
using CWQ = common::ChunkedWaitQueue<T>;

#define INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(...) \
    template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const CWQ<__VA_ARGS__> &, \
                                uint8_t, const std::string &dir); \
    template BOSS::Chunk::Chunk(uint64_t, size_t, bool, const Vector<__VA_ARGS__> &, \
                                uint8_t, const std::string &dir);

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
    assert(weights_.size() == 0);

    W_.push_back(W);
    last_.push_back(last);
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

    if (bool(weights_.size()) != bool(other.weights_.size())) {
        std::cerr << "ERROR: trying to concatenate weighted and unweighted blocks" << std::endl;
        exit(1);
    }

    auto &other_W = const_cast<sdsl::int_vector_buffer<>&>(other.W_);
    for (auto it = other_W.begin() + 1; it != other_W.end(); ++it) {
        W_.push_back(*it);
    }

    auto &other_last = const_cast<sdsl::int_vector_buffer<1>&>(other.last_);
    for (auto it = other_last.begin() + 1; it != other_last.end(); ++it) {
        last_.push_back(*it);
    }

    assert(F_.size() == other.F_.size());
    for (size_t p = 0; p < other.F_.size(); ++p) {
        F_[p] += other.F_[p];
    }

    if (other.weights_.size()) {
        auto &other_weights = const_cast<sdsl::int_vector_buffer<>&>(other.weights_);
        for (auto it = other_weights.begin() + 1; it != other_weights.end(); ++it) {
            weights_.push_back(*it);
        }
    }

    size_ += other.size_ - 1;

    assert(W_.size() == last_.size());
    assert(!weights_.size() || weights_.size() == W_.size());
}

void BOSS::Chunk::initialize_boss(BOSS *graph, sdsl::int_vector<> *weights) {
    assert(size_ <= W_.size());
    assert(last_.size() == W_.size());
    assert(weights_.size() == 0 || weights_.size() == W_.size());

    assert(graph->W_);
    delete graph->W_;
    graph->W_ = new wavelet_tree_small(get_W_width(), W_);

    assert(graph->last_);
    delete graph->last_;
    sdsl::bit_vector last(last_.size());
    std::copy(last_.begin(), last_.end(), last.begin());
    graph->last_ = new bit_vector_stat(std::move(last));

    graph->F_ = F_;
    graph->recompute_NF();

    graph->k_ = k_;
    // TODO:
    // graph->alph_size = alph_size_;

    graph->state = BOSS::State::SMALL;

    assert(graph->is_valid());

    if (weights) {
        weights->resize(0);
        weights->width(weights_.width());
        weights->resize(weights_.size());
        std::copy(weights_.begin(), weights_.end(), weights->begin());
    }
}

std::pair<BOSS*, bool>
BOSS::Chunk::build_boss_from_chunks(const std::vector<std::string> &chunk_filenames,
                                    bool verbose,
                                    sdsl::int_vector<> *weights,
                                    const std::string &swap_dir) {
    assert(chunk_filenames.size());

    if (!chunk_filenames.size())
        return std::make_pair(nullptr, false);

    BOSS *graph = new BOSS();

    std::unique_ptr<BOSS::Chunk> full_chunk;

    for (size_t i = 0; i < chunk_filenames.size(); ++i) {
        auto filename = utils::remove_suffix(chunk_filenames[i], kFileExtension)
                                                        + kFileExtension;

        auto graph_chunk = std::make_unique<BOSS::Chunk>(1, 0, false, swap_dir);
        if (!graph_chunk->load(filename)) {
            std::cerr << "ERROR: File corrupted. Cannot load graph chunk "
                      << filename << std::endl;
            exit(1);
        }

        if (verbose) {
            std::cout << "Chunk " << filename << " loaded..." << std::flush;
        }

        if (i == 0) {
            full_chunk = std::move(graph_chunk);
        } else {
            full_chunk->extend(*graph_chunk);
        }

        if (verbose) {
            std::cout << " concatenated" << std::endl;
        }
    }

    full_chunk->initialize_boss(graph, weights);

    return std::make_pair(graph, full_chunk->canonical_);
}

bool BOSS::Chunk::load(const std::string &infbase) {
    try {
        std::ifstream instream(utils::remove_suffix(infbase, kFileExtension)
                                                                + kFileExtension,
                               std::ios::binary);
        size_ = load_number(instream);

        {
            W_ = sdsl::int_vector_buffer<>();
            sdsl::int_vector<> W_copy;
            W_copy.load(instream);
            std::ofstream outstream(dir_ + "/W", std::ios::binary);
            W_copy.serialize(outstream);
            W_ = sdsl::int_vector_buffer<>(dir_ + "/W",
                                           std::ios::in | std::ios::out,
                                           BUFFER_SIZE, W_copy.width());
        }

        {
            last_ = sdsl::int_vector_buffer<1>();
            sdsl::int_vector<1> last_copy;
            last_copy.load(instream);
            std::ofstream outstream(dir_ + "/last", std::ios::binary);
            last_copy.serialize(outstream);
            last_ = sdsl::int_vector_buffer<1>(dir_ + "/last",
                                               std::ios::in | std::ios::out,
                                               BUFFER_SIZE, last_copy.width());
        }

        if (!load_number_vector(instream, &F_)) {
            std::cerr << "ERROR: failed to load F vector" << std::endl;
            return false;
        }

        {
            weights_ = sdsl::int_vector_buffer<>();
            sdsl::int_vector<> weights_copy;
            weights_copy.load(instream);
            std::ofstream outstream(dir_ + "/weights", std::ios::binary);
            weights_copy.serialize(outstream);
            weights_ = sdsl::int_vector_buffer<>(dir_ + "/weights",
                                                 std::ios::in | std::ios::out,
                                                 BUFFER_SIZE, weights_copy.width());
        }

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

    sdsl::int_vector<> W_copy(W_.size(), 0, W_.width());
    auto &W = const_cast<sdsl::int_vector_buffer<>&>(W_);
    std::copy(W.begin(), W.end(), W_copy.begin());
    W_copy.serialize(outstream);
    W_copy = sdsl::int_vector<>();

    sdsl::bit_vector last_copy(last_.size(), 0, last_.width());
    auto &last = const_cast<sdsl::int_vector_buffer<1>&>(last_);
    std::copy(last.begin(), last.end(), last_copy.begin());
    last_copy.serialize(outstream);
    last_copy = sdsl::bit_vector();

    serialize_number_vector(outstream, F_);

    sdsl::int_vector<> weights_copy(weights_.size(), 0, weights_.width());
    auto &weights = const_cast<sdsl::int_vector_buffer<>&>(weights_);
    std::copy(weights.begin(), weights.end(), weights_copy.begin());
    weights_copy.serialize(outstream);
    weights_copy = sdsl::int_vector<>();

    serialize_number(outstream, alph_size_);
    serialize_number(outstream, k_);
    serialize_number(outstream, canonical_);
}

uint8_t BOSS::Chunk::get_W_width() const {
    return alph_size_ ? sdsl::bits::hi(alph_size_ * 2 - 1) + 1 : 1;
}
