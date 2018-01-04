#include "dbg_succinct_chunk.hpp"

/**
 * We use libmaus 2 for representing dynamic succint data structures
 * such as the dynamic bit array and the dynamic wavelet tree.
 */
#include "serialization.hpp"

using libmaus2::util::NumberSerialisation;
using libmaus2::util::NumberSerialisation;


void DBG_succ::DynamicChunk::push_back(TAlphabet W, TAlphabet F, bool last) {
    W_.insert(W_.size(), W);
    for (TAlphabet a = F + 1; a < F_.size(); ++a) {
        F_[a]++;
    }
    last_.insertBit(last_.size(), last);
}

TAlphabet DBG_succ::DynamicChunk::get_W_back() const { return W_[W_.size() - 1]; }

void DBG_succ::DynamicChunk::alter_W_back(TAlphabet W) { W_.set(W_.size() - 1, W); }

void DBG_succ::DynamicChunk::alter_last_back(bool last) { last_.set(last_.size() - 1, last); }

uint64_t DBG_succ::DynamicChunk::size() const { return W_.size(); }

void DBG_succ::DynamicChunk::extend(const Chunk &other) {
    extend(dynamic_cast<const DynamicChunk&>(other));
}

void DBG_succ::DynamicChunk::extend(const DynamicChunk &other) {
    for (uint64_t j = 0; j < other.W_.size(); ++j) {
        W_.insert(W_.size(), other.W_[j]);
        last_.insertBit(last_.size(), other.last_[j]);
    }

    assert(F_.size() == other.F_.size());
    for (size_t p = 0; p < other.F_.size(); ++p) {
        F_[p] += other.F_[p];
    }
}

void DBG_succ::DynamicChunk::initialize_graph(DBG_succ *graph) {
    delete graph->W;
    graph->W = new wavelet_tree_dyn(4, W_);

    delete graph->last;
    graph->last = new bit_vector_dyn(last_);

    graph->F = F_;
}

bool DBG_succ::DynamicChunk::load(const std::string &infbase) {
    F_.resize(0);

    try {
        std::ifstream instream(infbase + ".F.dbg");
        std::string cur_line;

        while (std::getline(instream, cur_line)) {
            F_.push_back(std::stoull(cur_line));
        }
        instream.close();

        if (F_.size() != DBG_succ::alph_size)
            return false;

        // load W and last arrays
        std::ifstream instream_W(infbase + ".W.dbg");
        std::ifstream instream_l(infbase + ".l.dbg");
        return W_.deserialise(instream_W) && last_.deserialise(instream_l);
    } catch (...) {
        return false;
    }
}

void DBG_succ::DynamicChunk::serialize(const std::string &outbase) const {
    std::ofstream outstream(outbase + ".W.dbg");
    W_.serialise(outstream);
    outstream.close();

    // write last array
    outstream.open(outbase + ".l.dbg");
    last_.serialise(outstream);
    outstream.close();

    // write F values and k
    outstream.open(outbase + ".F.dbg");
    for (size_t i = 0; i < F_.size(); ++i) {
        outstream << F_.at(i) << "\n";
    }
    outstream.close();
}


void DBG_succ::VectorChunk::push_back(TAlphabet W, TAlphabet F, bool last) {
    W_.push_back(W);
    for (TAlphabet a = F + 1; a < F_.size(); ++a) {
        F_[a]++;
    }
    last_.push_back(last);
}

TAlphabet DBG_succ::VectorChunk::get_W_back() const { return W_.back(); }

void DBG_succ::VectorChunk::alter_W_back(TAlphabet W) { W_.back() = W; }

void DBG_succ::VectorChunk::alter_last_back(bool last) { last_.back() = last; }

uint64_t DBG_succ::VectorChunk::size() const { return W_.size(); }

void DBG_succ::VectorChunk::extend(const DBG_succ::Chunk &other) {
    extend(dynamic_cast<const DBG_succ::VectorChunk&>(other));
}

void DBG_succ::VectorChunk::extend(const DBG_succ::VectorChunk &other) {
    W_.insert(W_.end(), other.W_.begin(), other.W_.end());
    last_.insert(last_.end(), other.last_.begin(), other.last_.end());

    assert(F_.size() == other.F_.size());
    for (size_t p = 0; p < other.F_.size(); ++p) {
        F_[p] += other.F_[p];
    }
}

void DBG_succ::VectorChunk::initialize_graph(DBG_succ *graph) {
    delete graph->W;
    graph->W = new wavelet_tree_dyn(4, W_);

    delete graph->last;
    graph->last = new bit_vector_dyn(last_);

    graph->F = F_;
}

bool DBG_succ::VectorChunk::load(const std::string &infbase) {
    try {
        std::ifstream instream_W(infbase + ".W.dbg");
        W_ = NumberSerialisation::deserialiseNumberVector<TAlphabet>(instream_W);
        instream_W.close();

        std::ifstream instream_l(infbase + ".l.dbg");
        last_ = NumberSerialisation::deserialiseNumberVector<bool>(instream_l);
        instream_l.close();

        std::ifstream instream_F(infbase + ".F.dbg");
        F_ = NumberSerialisation::deserialiseNumberVector<uint64_t>(instream_F);
        instream_F.close();

        return F_.size() == DBG_succ::alph_size;
    } catch (...) {
        return false;
    }
}

void DBG_succ::VectorChunk::serialize(const std::string &outbase) const {
    std::ofstream outstream(outbase + ".W.dbg");
    NumberSerialisation::serialiseNumberVector(outstream, W_);
    outstream.close();

    outstream.open(outbase + ".l.dbg");
    NumberSerialisation::serialiseNumberVector(outstream, last_);
    outstream.close();

    outstream.open(outbase + ".F.dbg");
    NumberSerialisation::serialiseNumberVector(outstream, F_);
    outstream.close();
}
