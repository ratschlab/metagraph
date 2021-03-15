#include "boss_chunk.hpp"

#include <iostream>

#include "common/threads/chunked_wait_queue.hpp"
#include "common/serialization.hpp"
#include "common/vector.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/utils/template_utils.hpp"
#include "common/utils/file_utils.hpp"


namespace mtg {
namespace graph {
namespace boss {

using utils::get_first;
using mtg::kmer::KmerExtractorBOSS;
namespace fs = std::filesystem;

const uint64_t BUFFER_SIZE = 1024 * 1024; // 1 MiB

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

    W->flush();
    last->flush();
    if (weights)
        weights->flush();

    assert(W->size() == curpos);
    assert(last->size() == curpos);
    assert(!weights || weights->size() == curpos);
}

BOSS::Chunk::Chunk(uint64_t alph_size, size_t k, const std::string &swap_dir)
      : alph_size_(alph_size), k_(k),
        dir_(utils::create_temp_dir(swap_dir, "graph_chunk")),
        W_(dir_ + "/W", std::ios::out, BUFFER_SIZE, get_W_width()),
        last_(dir_ + "/last", std::ios::out, BUFFER_SIZE),
        weights_(dir_ + "/weights", std::ios::out, BUFFER_SIZE) {
    W_.push_back(0);
    last_.push_back(0);
    F_.assign(alph_size_, 0);
}

BOSS::Chunk::~Chunk() {
    try {
        W_.close(true);
        last_.close(true);
        weights_.close(true);

    } catch (const std::exception &e) {
        std::cerr << "ERROR: Failed to destruct BOSS::Chunk: "
                  << e.what() << std::endl;
    } catch (...) {
        std::cerr << "ERROR: Failed to destruct BOSS::Chunk";
    }

    // remove the temp directory, but only if it was initialized
    if (!dir_.empty())
        utils::remove_temp_dir(dir_);
}

template <typename Array>
BOSS::Chunk::Chunk(uint64_t alph_size,
                   size_t k,
                   const Array &kmers_with_counts,
                   uint8_t bits_per_count,
                   const std::string &swap_dir)
      : Chunk(alph_size, k, swap_dir) {

    weights_ = sdsl::int_vector_buffer<>(dir_ + "/weights",
                                         std::ios::out, BUFFER_SIZE,
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
}

template <typename T>
using CWQ = common::ChunkedWaitQueue<T>;

#define INSTANTIATE_BOSS_CHUNK_CONSTRUCTORS(...) \
    template BOSS::Chunk::Chunk(uint64_t, size_t, const CWQ<__VA_ARGS__> &, \
                                uint8_t, const std::string &dir); \
    template BOSS::Chunk::Chunk(uint64_t, size_t, const Vector<__VA_ARGS__> &, \
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
    assert(!weights_.size());

    W_.push_back(W);
    last_.push_back(last);
    for (TAlphabet a = F + 1; a < F_.size(); ++a) {
        F_[a]++;
    }
}

void BOSS::Chunk::extend(Chunk &other) {
    assert(!weights_.size() || weights_.size() == W_.size());
    assert(!other.weights_.size() || other.weights_.size() == other.W_.size());

    if (alph_size_ != other.alph_size_
            || k_ != other.k_
            || W_.width() != other.W_.width()) {
        std::cerr << "ERROR: trying to concatenate incompatible graph chunks" << std::endl;
        exit(1);
    }

    if (other.W_.size() == 1)
        return;

    if (bool(weights_.size()) != bool(other.weights_.size())) {
        std::cerr << "ERROR: trying to concatenate weighted and unweighted blocks" << std::endl;
        exit(1);
    }

    for (auto it = other.W_.begin() + 1; it != other.W_.end(); ++it) {
        W_.push_back(*it);
    }

    for (auto it = other.last_.begin() + 1; it != other.last_.end(); ++it) {
        last_.push_back(*it);
    }

    assert(F_.size() == other.F_.size());
    for (size_t p = 0; p < other.F_.size(); ++p) {
        F_[p] += other.F_[p];
    }

    if (other.weights_.size()) {
        for (auto it = other.weights_.begin() + 1; it != other.weights_.end(); ++it) {
            weights_.push_back(*it);
        }
    }

    assert(W_.size() == last_.size());
    assert(!weights_.size() || weights_.size() == W_.size());
}

void BOSS::Chunk::initialize_boss(BOSS *graph, sdsl::int_vector<> *weights) {
    assert(last_.size() == W_.size());
    assert(!weights_.size() || weights_.size() == W_.size());

    graph->initialize(this);

    assert(graph->is_valid());

    if (weights) {
        weights_.flush();
        std::ifstream in(weights_.filename(), std::ios::binary);
        weights->load(in);
    }
}

BOSS*
BOSS::Chunk::build_boss_from_chunks(const std::vector<std::string> &chunk_filenames,
                                    bool verbose,
                                    sdsl::int_vector<> *weights,
                                    const std::string &swap_dir) {
    assert(chunk_filenames.size());

    if (!chunk_filenames.size())
        return nullptr;

    BOSS *graph = new BOSS();

    Chunk full_chunk;

    for (size_t i = 0; i < chunk_filenames.size(); ++i) {
        auto filename = utils::remove_suffix(chunk_filenames[i], kFileExtension)
                                                        + kFileExtension;

        Chunk graph_chunk(1, 0, swap_dir);
        if (!graph_chunk.load(filename)) {
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
            full_chunk.extend(graph_chunk);
        }

        if (verbose) {
            std::cout << " concatenated" << std::endl;
        }
    }

    full_chunk.initialize_boss(graph, weights);

    return graph;
}

bool BOSS::Chunk::load(const std::string &infbase) {
    std::string fname
        = utils::remove_suffix(infbase, kFileExtension) + kFileExtension;

    try {
        W_.close(true);
        last_.close(true);
        weights_.close(true);

        fs::copy(fname + ".W", dir_ + "/W", fs::copy_options::overwrite_existing);
        fs::copy(fname + ".last", dir_ + "/last", fs::copy_options::overwrite_existing);
        fs::copy(fname + ".weights", dir_ + "/weights", fs::copy_options::overwrite_existing);

        W_ = sdsl::int_vector_buffer<>(dir_ + "/W",
                                       std::ios::in | std::ios::out, BUFFER_SIZE);
        last_ = sdsl::int_vector_buffer<1>(dir_ + "/last",
                                           std::ios::in | std::ios::out, BUFFER_SIZE);
        weights_ = sdsl::int_vector_buffer<>(dir_ + "/weights",
                                             std::ios::in | std::ios::out, BUFFER_SIZE);

        std::ifstream instream(fname, std::ios::binary);

        if (!load_number_vector(instream, &F_)) {
            std::cerr << "ERROR: failed to load F vector" << std::endl;
            return false;
        }

        alph_size_ = load_number(instream);
        k_ = load_number(instream);

        return k_ && alph_size_ && W_.size() == last_.size()
                                && F_.size() == alph_size_
                                && (!weights_.size() || weights_.size() == W_.size());

    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load BOSS chunk from "
                  << fname << "." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void BOSS::Chunk::serialize(const std::string &outbase) {
    std::string fname
        = utils::remove_suffix(outbase, kFileExtension) + kFileExtension;

    W_.flush();
    fs::copy(W_.filename(), fname + ".W", fs::copy_options::overwrite_existing);
    last_.flush();
    fs::copy(last_.filename(), fname + ".last", fs::copy_options::overwrite_existing);
    weights_.flush();
    fs::copy(weights_.filename(), fname + ".weights", fs::copy_options::overwrite_existing);

    std::ofstream outstream(fname, std::ios::binary);
    serialize_number_vector(outstream, F_);
    serialize_number(outstream, alph_size_);
    serialize_number(outstream, k_);
}

uint8_t BOSS::Chunk::get_W_width() const {
    return alph_size_ ? sdsl::bits::hi(alph_size_ * 2 - 1) + 1 : 1;
}

} // namespace boss
} // namespace graph
} // namespace mtg
