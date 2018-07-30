#include "dbg_succinct_chunk.hpp"

#include "serialization.hpp"
#include "utils.hpp"


DBG_succ::Chunk::Chunk(size_t k)
      : k_(k), W_(1, 0), last_(1, 0), F_(DBG_succ::alph_size, 0) {}

DBG_succ::Chunk::Chunk(size_t k,
                       std::vector<TAlphabet>&& W,
                       std::vector<bool>&& last,
                       std::vector<uint64_t>&& F)
      : k_(k), W_(std::move(W)), last_(std::move(last)), F_(std::move(F)) {}

void DBG_succ::Chunk::push_back(TAlphabet W, TAlphabet F, bool last) {
    W_.push_back(W);
    for (TAlphabet a = F + 1; a < F_.size(); ++a) {
        F_[a]++;
    }
    last_.push_back(last);
}

TAlphabet DBG_succ::Chunk::get_W_back() const { return W_.back(); }

void DBG_succ::Chunk::alter_W_back(TAlphabet W) { W_.back() = W; }

void DBG_succ::Chunk::alter_last_back(bool last) { last_.back() = last; }

uint64_t DBG_succ::Chunk::size() const { return W_.size() - 1; }

void DBG_succ::Chunk::extend(const DBG_succ::Chunk &other) {
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

void DBG_succ::Chunk::initialize_graph(DBG_succ *graph) const {
    delete graph->W_;
    graph->W_ = new wavelet_tree_stat(kLogSigma, W_);

    delete graph->last_;
    graph->last_ = new bit_vector_stat(last_);

    graph->F_ = F_;

    graph->k_ = k_;

    graph->state = Config::STAT;

    assert(graph->is_valid());
}

DBG_succ*
DBG_succ::Chunk::build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                                         bool verbose) {
    DBG_succ *graph = new DBG_succ();
    if (!chunk_filenames.size())
        return graph;

    uint64_t cumulative_size = 1;

    for (auto file : chunk_filenames) {
        file = utils::remove_suffix(file, ".dbgchunk") + ".dbgchunk";

        std::ifstream chunk_in(file);

        if (!chunk_in.good()) {
            std::cerr << "ERROR: input file " << file << " corrupted" << std::endl;
            exit(1);
        }
        cumulative_size += load_number_vector_size(chunk_in) - 1;
    }

    if (verbose)
        std::cout << "Cumulative size of chunks: "
                  << cumulative_size << std::endl;

    sdsl::int_vector<> W(cumulative_size, 0, kLogSigma);
    sdsl::bit_vector last(cumulative_size, 0);
    std::vector<uint64_t> F(DBG_succ::alph_size, 0);
    uint64_t pos = 1;

    if (verbose)
        std::cout << "Succinct arrays initialized" << std::endl;

    for (size_t i = 0; i < chunk_filenames.size(); ++i) {
        auto filename = utils::remove_suffix(chunk_filenames[i], ".dbgchunk")
                                                        + ".dbgchunk";
        DBG_succ::Chunk graph_chunk(0);
        if (!graph_chunk.load(filename)) {
            std::cerr << "ERROR: input file "
                      << filename << " corrupted" << std::endl;
            exit(1);
        } else if (!graph_chunk.k_) {
            std::cerr << "ERROR: trying to load invalid graph chunk from file "
                      << filename << std::endl;
            exit(1);
        } else if (i == 0) {
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
            W[pos] = graph_chunk.W_[i];
            last[pos] = graph_chunk.last_[i];
            pos++;
        }

        assert(graph_chunk.F_.size() == F.size());
        for (size_t p = 0; p < F.size(); ++p) {
            F[p] += graph_chunk.F_[p];
        }

        if (verbose) {
            std::cout << " concatenated" << std::endl;
        }
    }

    delete graph->W_;
    graph->W_ = new wavelet_tree_stat(kLogSigma, std::move(W));

    delete graph->last_;
    graph->last_ = new bit_vector_stat(std::move(last));

    graph->F_ = std::move(F);

    graph->state = Config::STAT;

    assert(graph->is_valid());

    return graph;
}

bool DBG_succ::Chunk::load(const std::string &infbase) {
    try {
        std::ifstream instream(utils::remove_suffix(infbase, ".dbgchunk")
                                                                + ".dbgchunk");
        W_ = load_number_vector<TAlphabet>(instream);
        last_ = load_number_vector<bool>(instream);
        F_ = load_number_vector<uint64_t>(instream);
        k_ = load_number(instream);

        instream.close();

        return F_.size() == DBG_succ::alph_size;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load graph chunk from "
                  << infbase << "." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void DBG_succ::Chunk::serialize(const std::string &outbase) const {
    std::ofstream outstream(utils::remove_suffix(outbase, ".dbgchunk")
                                                                + ".dbgchunk");
    serialize_number_vector(outstream, W_, kLogSigma);
    serialize_number_vector(outstream, last_, 1);
    serialize_number_vector(outstream, F_);
    serialize_number(outstream, k_);

    outstream.close();
}
