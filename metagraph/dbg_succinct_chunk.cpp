#include "dbg_succinct_chunk.hpp"

/**
 * We use libmaus 2 for representing dynamic succint data structures
 * such as the dynamic bit array and the dynamic wavelet tree.
 */
#include "serialization.hpp"
#include "kmer.hpp"

#ifdef _PROTEIN_GRAPH
const size_t kLogSigma = 6;
#else
const size_t kLogSigma = 4;
#endif

static_assert(sizeof(TAlphabet) * 8 >= kLogSigma,
              "Choose the TAlphabet type accordingly");


DBG_succ::Chunk::Chunk()
      : W_(1, 0), last_(1, 0), F_(DBG_succ::alph_size, 0) {}

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
    delete graph->W;
    graph->W = new wavelet_tree_stat(kLogSigma, W_);

    delete graph->last;
    graph->last = new bit_vector_stat(last_);

    graph->F = F_;

    graph->state = Config::STAT;

    assert(graph->is_valid());
}

DBG_succ* DBG_succ::Chunk::build_graph_from_chunks(size_t k,
                        const std::vector<std::string> &chunk_filenames,
                        bool verbose) {
    DBG_succ *graph = new DBG_succ(k);
    if (!chunk_filenames.size())
        return graph;

    uint64_t cumulative_size = 1;
    for (const auto &file : chunk_filenames) {
        std::ifstream chunk_in(file + ".dbgchunk");
        if (!chunk_in.good()) {
            std::cerr << "ERROR: input file "
                      << file + ".dbgchunk" << " corrupted" << std::endl;
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

    for (const auto &filename : chunk_filenames) {
        DBG_succ::Chunk graph_chunk;

        if (!graph_chunk.load(filename)) {
            std::cerr << "ERROR: input file "
                      << filename << " corrupted" << std::endl;
            exit(1);
        }
        if (verbose) {
            std::cout << "Chunk " << filename
                      << " loaded" << std::endl;
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
            std::cout << "Chunk " << filename
                      << " concatenated" << std::endl;
        }
    }

    delete graph->W;
    graph->W = new wavelet_tree_stat(std::move(W));

    delete graph->last;
    graph->last = new bit_vector_stat(std::move(last));

    graph->F = std::move(F);

    graph->state = Config::STAT;

    assert(graph->is_valid());

    return graph;
}

bool DBG_succ::Chunk::load(const std::string &infbase) {
    try {
        std::ifstream instream(infbase + ".dbgchunk");

        W_ = load_number_vector<TAlphabet>(instream);
        last_ = load_number_vector<bool>(instream);
        F_ = load_number_vector<uint64_t>(instream);

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
    std::ofstream outstream(outbase + ".dbgchunk");
    serialize_number_vector(outstream, W_, kLogSigma);
    serialize_number_vector(outstream, last_, 1);
    serialize_number_vector(outstream, F_);
    outstream.close();
}

DBG_succ::Chunk* DBG_succ::Chunk::build_from_kmers(size_t k,
                                                   std::vector<KMer> *kmers) {
    assert(std::is_sorted(kmers->begin(), kmers->end()));

    DBG_succ::Chunk *result = new DBG_succ::Chunk();

    // the bit array indicating the last outgoing edge of a node
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
            for (size_t j = i - 1; KMer::compare_suffix(kmers->at(i),
                                                        kmers->at(j), 1); --j) {
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
