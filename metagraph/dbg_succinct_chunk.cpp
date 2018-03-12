#include "dbg_succinct_chunk.hpp"

/**
 * We use libmaus 2 for representing dynamic succint data structures
 * such as the dynamic bit array and the dynamic wavelet tree.
 */
#include "serialization.hpp"
#include "kmer.hpp"

using libmaus2::util::NumberSerialisation;


DBG_succ::VectorChunk::VectorChunk()
      : W_(1, 0), last_(1, 0), F_(DBG_succ::alph_size, 0) {}

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

uint64_t DBG_succ::VectorChunk::size() const { return W_.size() - 1; }

void DBG_succ::VectorChunk::extend(const DBG_succ::Chunk &other) {
    extend(dynamic_cast<const DBG_succ::VectorChunk&>(other));
}

void DBG_succ::VectorChunk::extend(const DBG_succ::VectorChunk &other) {
    W_.reserve(W_.size() + other.size());
    W_.insert(W_.end(), other.W_.begin() + 1, other.W_.end());

    last_.reserve(last_.size() + other.size());
    last_.insert(last_.end(), other.last_.begin() + 1, other.last_.end());

    assert(F_.size() == other.F_.size());
    for (size_t p = 0; p < other.F_.size(); ++p) {
        F_[p] += other.F_[p];
    }
}

void DBG_succ::VectorChunk::initialize_graph(DBG_succ *graph) {
    delete graph->W;
    graph->W = new wavelet_tree_stat(4, W_);

    delete graph->last;
    graph->last = new bit_vector_stat(last_);

    graph->F = F_;

    graph->state = Config::STAT;

    assert(graph->is_valid());
}

bool DBG_succ::VectorChunk::load(const std::string &infbase) {
    try {
        std::ifstream instream(infbase + ".dbgchunk");

        W_ = NumberSerialisation::deserialiseNumberVector<TAlphabet>(instream);
        last_ = NumberSerialisation::deserialiseNumberVector<bool>(instream);
        F_ = NumberSerialisation::deserialiseNumberVector<uint64_t>(instream);

        instream.close();

        return F_.size() == DBG_succ::alph_size;
    } catch (...) {
        return false;
    }
}

void DBG_succ::VectorChunk::serialize(const std::string &outbase) const {
    std::ofstream outstream(outbase + ".dbgchunk");
    NumberSerialisation::serialiseNumberVector(outstream, W_);
    NumberSerialisation::serialiseNumberVector(outstream, last_);
    NumberSerialisation::serialiseNumberVector(outstream, F_);
    outstream.close();
}

DBG_succ::VectorChunk* DBG_succ::VectorChunk::build_from_kmers(size_t k,
                                                               std::vector<KMer> *kmers) {
    assert(std::is_sorted(kmers->begin(), kmers->end()));

    DBG_succ::VectorChunk *result = new DBG_succ::VectorChunk();

    // the bit array indicating the last outgoing edge of a node (static container for full init)
    std::vector<uint8_t> last_stat_safe(1 + kmers->size(), 1);
    last_stat_safe[0] = 0;

    // the array containing the edge labels
    std::vector<TAlphabet> &W_stat = result->W_;
    W_stat.resize(1 + kmers->size());
    W_stat[0] = 0;

    result->F_.at(0) = 0;

    size_t curpos = 1;
    TAlphabet lastF = 0;

    for (size_t i = 0; i < kmers->size(); ++i) {
        TAlphabet curW = kmers->at(i)[0];
        TAlphabet curF = kmers->at(i)[k];

        assert(curW < DBG_succ::alph_size);

        // check redundancy and set last
        if (i + 1 < kmers->size()
                && KMer::compare_suffix(kmers->at(i), kmers->at(i + 1))) {
            // skip redundant dummy edges
            if (curW == 0 && curF > 0)
                continue;

            last_stat_safe[curpos] = 0;
        }
        //set W
        if (i > 0) {
            for (size_t j = i - 1; KMer::compare_suffix(kmers->at(i), kmers->at(j), 1); --j) {
                //TODO: recalculating W is probably faster than doing a pragma for ordered
                if (curW > 0 && kmers->at(j)[0] == curW) {
                    curW += DBG_succ::alph_size;
                    break;
                }
                if (j == 0)
                    break;
            }
        }
        W_stat[curpos] = curW;

        while (lastF + 1 < DBG_succ::alph_size && curF != lastF) {
            result->F_.at(++lastF) = curpos - 1;
        }
        curpos++;
    }
    while (++lastF < DBG_succ::alph_size) {
        result->F_.at(lastF) = curpos - 1;
    }

    kmers->clear();

    W_stat.resize(curpos);
    last_stat_safe.resize(curpos);
    result->last_.assign(last_stat_safe.begin(), last_stat_safe.end());

    return result;
}
