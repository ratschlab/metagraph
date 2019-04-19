#include "boss_chunk.hpp"

#include <boost/multiprecision/integer.hpp>

#include "serialization.hpp"
#include "utils.hpp"


BOSS::Chunk::Chunk(size_t k)
      : alph_size_(KmerExtractor::alphabet.size()),
        bits_per_char_W_(boost::multiprecision::msb(alph_size_ - 1) + 2),
        k_(k), W_(1, 0), last_(1, 0), F_(alph_size_, 0) {
    assert(sizeof(TAlphabet) * 8 >= bits_per_char_W_);
    assert(alph_size_ * 2 < 1llu << bits_per_char_W_);
}

BOSS::Chunk::Chunk(size_t k,
                   std::vector<TAlphabet>&& W,
                   std::vector<bool>&& last,
                   std::vector<uint64_t>&& F)
      : alph_size_(KmerExtractor::alphabet.size()),
        bits_per_char_W_(boost::multiprecision::msb(alph_size_ - 1) + 2),
        k_(k), W_(std::move(W)), last_(std::move(last)), F_(std::move(F)) {
    assert(sizeof(TAlphabet) * 8 >= bits_per_char_W_);
    assert(alph_size_ * 2 < 1llu << bits_per_char_W_);
}

void BOSS::Chunk::push_back(TAlphabet W, TAlphabet F, bool last) {
    W_.push_back(W);
    for (TAlphabet a = F + 1; a < F_.size(); ++a) {
        F_[a]++;
    }
    last_.push_back(last);
}

BOSS::Chunk::TAlphabet
BOSS::Chunk::get_W_back() const { return W_.back(); }

void BOSS::Chunk::alter_W_back(TAlphabet W) { W_.back() = W; }

void BOSS::Chunk::alter_last_back(bool last) { last_.back() = last; }

uint64_t BOSS::Chunk::size() const { return W_.size() - 1; }

void BOSS::Chunk::extend(const BOSS::Chunk &other) {
    assert(k_ == other.k_);

    W_.reserve(W_.size() + other.size());
    W_.insert(W_.end(), other.W_.begin() + 1, other.W_.end());

    last_.reserve(last_.size() + other.size());
    last_.insert(last_.end(), other.last_.begin() + 1, other.last_.end());

    assert(F_.size() == other.F_.size());
    for (size_t p = 0; p < other.F_.size(); ++p) {
        F_[p] += other.F_[p];
    }
}

void BOSS::Chunk::initialize_boss(BOSS *graph) const {
    delete graph->W_;
    graph->W_ = new wavelet_tree_stat(bits_per_char_W_, W_);

    delete graph->last_;
    graph->last_ = new bit_vector_stat(last_);

    graph->F_ = F_;

    graph->k_ = k_;

    graph->state = Config::STAT;

    assert(graph->is_valid());
}

BOSS*
BOSS::Chunk::build_boss_from_chunks(const std::vector<std::string> &chunk_filenames,
                                    bool verbose) {
    BOSS *graph = new BOSS();
    if (!chunk_filenames.size())
        return graph;

    uint64_t cumulative_size = 1;

    for (auto file : chunk_filenames) {
        file = utils::remove_suffix(file, ".dbgchunk") + ".dbgchunk";

        std::ifstream chunk_in(file, std::ios::binary);

        if (!chunk_in.good()) {
            std::cerr << "ERROR: input file " << file << " corrupted" << std::endl;
            exit(1);
        }
        cumulative_size += get_number_vector_size(chunk_in) - 1;
    }

    if (verbose)
        std::cout << "Cumulative size of chunks: "
                  << cumulative_size << std::endl;

    std::unique_ptr<sdsl::int_vector<>> W;
    std::unique_ptr<sdsl::bit_vector> last;
    std::unique_ptr<std::vector<uint64_t>> F;

    uint64_t pos = 1;

    if (verbose)
        std::cout << "Succinct arrays initialized" << std::endl;

    for (size_t i = 0; i < chunk_filenames.size(); ++i) {
        auto filename = utils::remove_suffix(chunk_filenames[i], ".dbgchunk")
                                                        + ".dbgchunk";
        BOSS::Chunk graph_chunk(0);
        if (!graph_chunk.load(filename)) {
            std::cerr << "ERROR: input file "
                      << filename << " corrupted" << std::endl;
            exit(1);
        } else if (!graph_chunk.k_) {
            std::cerr << "ERROR: trying to load invalid graph chunk from file "
                      << filename << std::endl;
            exit(1);
        } else if (i == 0) {
            W = std::make_unique<sdsl::int_vector<>>(cumulative_size, 0, graph_chunk.bits_per_char_W_);
            last = std::make_unique<sdsl::bit_vector>(cumulative_size, 0);
            F = std::make_unique<std::vector<uint64_t>>(graph_chunk.alph_size_, 0);

            graph->k_ = graph_chunk.k_;
        } else if (graph->k_ != graph_chunk.k_) {
            std::cerr << "ERROR: trying to build a graph with k=" << graph->k_
                      << " from chunk " << filename << " with k=" << graph_chunk.k_
                      << std::endl;
            exit(1);
        }

        if (verbose) {
            std::cout << "Chunk " << filename << " loaded..." << std::flush;
        }

        for (size_t i = 1; i < graph_chunk.W_.size(); ++i) {
            (*W)[pos] = graph_chunk.W_[i];
            (*last)[pos] = graph_chunk.last_[i];
            pos++;
        }

        assert(graph_chunk.F_.size() == F->size());
        for (size_t p = 0; p < F->size(); ++p) {
            (*F)[p] += graph_chunk.F_[p];
        }

        if (verbose) {
            std::cout << " concatenated" << std::endl;
        }
    }

    assert(W.get());
    assert(last.get());
    assert(F.get());

    delete graph->W_;
    graph->W_ = new wavelet_tree_stat(W->width(), std::move(*W));
    W.reset();

    delete graph->last_;
    graph->last_ = new bit_vector_stat(std::move(*last));
    last.reset();

    graph->F_ = std::move(*F);
    F.reset();

    graph->state = Config::STAT;

    assert(graph->is_valid());

    return graph;
}

bool BOSS::Chunk::load(const std::string &infbase) {
    try {
        std::ifstream instream(utils::remove_suffix(infbase, ".dbgchunk")
                                                                + ".dbgchunk",
                               std::ios::binary);
        W_ = load_number_vector<TAlphabet>(instream);
        last_ = load_number_vector<bool>(instream);
        F_ = load_number_vector<uint64_t>(instream);
        k_ = load_number(instream);

        instream.close();

        return F_.size() == alph_size_;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load BOSS chunk from "
                  << infbase << "." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void BOSS::Chunk::serialize(const std::string &outbase) const {
    std::ofstream outstream(utils::remove_suffix(outbase, ".dbgchunk")
                                                                + ".dbgchunk",
                            std::ios::binary);
    serialize_number_vector(outstream, W_, bits_per_char_W_);
    serialize_number_vector(outstream, last_, 1);
    serialize_number_vector(outstream, F_);
    serialize_number(outstream, k_);

    outstream.close();
}
