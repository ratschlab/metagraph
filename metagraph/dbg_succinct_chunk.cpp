#include "dbg_succinct_chunk.hpp"

#include <parallel/algorithm>
/**
 * We use libmaus 2 for representing dynamic succint data structures
 * such as the dynamic bit array and the dynamic wavelet tree.
 */
#include "serialization.hpp"
#include "kmer.hpp"

using libmaus2::util::NumberSerialisation;


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
    W_.reserve(W_.size() + other.W_.size());
    W_.insert(W_.end(), other.W_.begin(), other.W_.end());

    last_.reserve(last_.size() + other.last_.size());
    last_.insert(last_.end(), other.last_.begin(), other.last_.end());

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

DBG_succ::VectorChunk* DBG_succ::VectorChunk::build_from_kmers(size_t k,
                                                               std::vector<KMer> *kmers,
                                                               unsigned int parallel) {
    // parallel sort of all kmers
    omp_set_num_threads(std::max(static_cast<int>(parallel), 1));
    __gnu_parallel::sort(kmers->data(), kmers->data() + kmers->size());

    auto unique_end = std::unique(kmers->begin(), kmers->end());
    kmers->erase(unique_end, kmers->end());

    DBG_succ::VectorChunk *result = new DBG_succ::VectorChunk();

    // the bit array indicating the last outgoing edge of a node (static container for full init)
    std::vector<uint8_t> last_stat_safe { 0 };

    last_stat_safe.resize(last_stat_safe.size() + kmers->size(), 1);

    // the array containing the edge labels
    std::vector<TAlphabet> &W_stat = result->W_;
    W_stat.push_back(0);

    size_t curpos = W_stat.size();
    W_stat.resize(W_stat.size() + kmers->size());

    #pragma omp parallel num_threads(parallel)
    {
        #pragma omp for nowait
        for (size_t i = 0; i < kmers->size(); ++i) {
            //set last
            if (i + 1 < kmers->size()) {
                if (KMer::compare_kmer_suffix((*kmers)[i], (*kmers)[i + 1])) {
                    last_stat_safe[curpos + i] = 0;
                }
            }
            //set W
            uint8_t curW = (*kmers)[i][0];
            if (curW == 127) {
                std::cerr << "Failure decoding kmer " << i << "\n" << (*kmers)[i] << "\n"
                          << (*kmers)[i].to_string(DBG_succ::alphabet) << "\n";
                exit(1);
            }
            if (i) {
                for (size_t j = i - 1; KMer::compare_kmer_suffix((*kmers)[j], (*kmers)[i], 1); --j) {
                    //TODO: recalculating W is probably faster than doing a pragma for ordered
                    if ((*kmers)[j][0] == curW) {
                        curW += DBG_succ::alph_size;
                        break;
                    }
                    if (!j)
                        break;
                }
            }
            W_stat[curpos + i] = curW;
        }
    }
    result->last_.assign(last_stat_safe.begin(), last_stat_safe.end());

    size_t i;
    size_t lastlet = 0;

    result->F_.at(0) = 0;
    for (i = 0; i < kmers->size(); ++i) {
        while (DBG_succ::alphabet[lastlet] != DBG_succ::alphabet[(*kmers)[i][k]]
                    && lastlet + 1 < DBG_succ::alph_size) {
            result->F_.at(++lastlet) = curpos + i - 1;
        }
    }
    while (++lastlet < DBG_succ::alph_size) {
        result->F_.at(lastlet) = curpos + i - 1;
    }

    kmers->clear();

    return result;
}


std::vector<TAlphabet> prepare_sequence(const std::string &sequence,
                                        size_t k, bool add_bridge) {
    const std::string &seq = add_bridge
                                ? std::string(k, '$') + sequence + "$"
                                : sequence;

    std::vector<TAlphabet> sequence_encoded(seq.size());
    std::transform(seq.begin(), seq.end(),
                   sequence_encoded.begin(), DBG_succ::encode);

    return sequence_encoded;
}


void sequence_to_kmers(const std::vector<TAlphabet> &seq,
                       size_t k,
                       std::vector<KMer> *kmers,
                       const std::vector<TAlphabet> &suffix) {
    if (seq.size() < k + 1)
        return;

    for (size_t i = 0; i < seq.size() - k; ++i) {
        if (std::equal(suffix.begin(), suffix.end(),
                       seq.data() + i + k - suffix.size())) {
            kmers->emplace_back(seq.data() + i, k + 1);
        }
    }
}
