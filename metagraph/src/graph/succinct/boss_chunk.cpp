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
                   std::vector<uint64_t>&& F,
                   sdsl::int_vector<>&& weights)
      : alph_size_(KmerExtractor::alphabet.size()),
        bits_per_char_W_(boost::multiprecision::msb(alph_size_ - 1) + 2),
        k_(k), W_(std::move(W)), last_(std::move(last)), F_(std::move(F)),
        weights_(std::move(weights)) {
    assert(sizeof(TAlphabet) * 8 >= bits_per_char_W_);
    assert(alph_size_ * 2 < 1llu << bits_per_char_W_);
    assert(F_.back() < W_.size());
    assert(W_.size() == last_.size());
    assert(weights_.empty() || weights_.size() == W_.size());
}

void BOSS::Chunk::push_back(TAlphabet W, TAlphabet F, bool last) {
    W_.push_back(W);
    for (TAlphabet a = F + 1; a < F_.size(); ++a) {
        F_[a]++;
    }
    last_.push_back(last);

    assert(weights_.empty());
}

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

    if (other.weights_.size()) {
        const auto start = weights_.size();
        weights_.resize(weights_.size() + other.size());
        std::copy(other.weights_.begin() + 1,
                  other.weights_.end(),
                  weights_.begin() + start);
    }

    if (W_.size() != last_.size()
            || (weights_.size() && weights_.size() != W_.size())) {
        std::cerr << "ERROR: trying to concatenate incompatible graph chunks" << std::endl;
        exit(1);
    }
}

void BOSS::Chunk::initialize_boss(BOSS *graph, sdsl::int_vector<> *weights) {
    assert(graph->W_);
    delete graph->W_;
    graph->W_ = new wavelet_tree_stat(bits_per_char_W_, std::move(W_));
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
    BOSS *graph = new BOSS();
    if (!chunk_filenames.size())
        return graph;

    uint64_t cumulative_size = 1;

    for (auto file : chunk_filenames) {
        file = utils::remove_suffix(file, kFileExtension) + kFileExtension;

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
        auto filename = utils::remove_suffix(chunk_filenames[i], kFileExtension)
                                                        + kFileExtension;
        BOSS::Chunk graph_chunk(0);
        if (!graph_chunk.load(filename)) {
            std::cerr << "ERROR: input file "
                      << filename << " corrupted" << std::endl;
            exit(1);

        } else if (!graph_chunk.k_
                    || graph_chunk.last_.size() != graph_chunk.W_.size()
                    || (graph_chunk.weights_.size()
                            && graph_chunk.weights_.size() != graph_chunk.W_.size())) {
            std::cerr << "ERROR: trying to load invalid graph chunk from file "
                      << filename << std::endl;
            exit(1);

        } else if (weights && graph_chunk.weights_.empty()) {
            std::cerr << "ERROR: no weights stored in graph chunk "
                      << filename << std::endl;
            exit(1);

        } else if (i == 0) {
            W = std::make_unique<sdsl::int_vector<>>(cumulative_size, 0, graph_chunk.bits_per_char_W_);
            last = std::make_unique<sdsl::bit_vector>(cumulative_size, 0);
            F = std::make_unique<std::vector<uint64_t>>(graph_chunk.alph_size_, 0);

            graph->k_ = graph_chunk.k_;

            if (weights) {
                (*weights) = graph_chunk.weights_;
                weights->resize(cumulative_size);
            }

        } else if (graph->k_ != graph_chunk.k_) {
            std::cerr << "ERROR: trying to build a graph with k=" << graph->k_
                      << " from chunk " << filename << " with k=" << graph_chunk.k_
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
                  W->begin() + pos);

        std::copy(graph_chunk.last_.begin() + 1,
                  graph_chunk.last_.end(),
                  last->begin() + pos);

        if (weights && i) {
            std::copy(graph_chunk.weights_.begin() + 1,
                      graph_chunk.weights_.end(),
                      weights->begin() + pos);
        }

        pos += graph_chunk.size();

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
        std::ifstream instream(utils::remove_suffix(infbase, kFileExtension)
                                                                + kFileExtension,
                               std::ios::binary);
        W_ = load_number_vector<TAlphabet>(instream);
        last_ = load_number_vector<bool>(instream);
        F_ = load_number_vector<uint64_t>(instream);
        k_ = load_number(instream);

        weights_.load(instream);

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
    std::ofstream outstream(utils::remove_suffix(outbase, kFileExtension)
                                                                + kFileExtension,
                            std::ios::binary);
    serialize_number_vector(outstream, W_, bits_per_char_W_);
    serialize_number_vector(outstream, last_, 1);
    serialize_number_vector(outstream, F_);
    serialize_number(outstream, k_);

    weights_.serialize(outstream);

    outstream.close();
}
