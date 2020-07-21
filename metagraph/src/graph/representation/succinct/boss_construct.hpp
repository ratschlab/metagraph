#ifndef __BOSS_CONSTRUCT_HPP__
#define __BOSS_CONSTRUCT_HPP__

#include "graph/representation/base/dbg_construct.hpp"
#include "boss_chunk_construct.hpp"


namespace mtg {
namespace graph {
namespace boss {

class BOSSConstructor : public IGraphConstructor<BOSS> {
  public:
    // see input arguments in IBOSSChunkConstructor::initialize
    template <typename... Args>
    BOSSConstructor(const Args&... args)
      : constructor_(IBOSSChunkConstructor::initialize(args...)) {}

    void add_sequence(std::string_view sequence, uint64_t count = 1) {
        constructor_->add_sequence(sequence, count);
    }

    void add_sequences(const std::function<void(CallString)> &generate_sequences) {
        constructor_->add_sequences(generate_sequences);
    }

    void add_sequences(const std::function<void(CallStringCount)> &generate_sequences) {
        constructor_->add_sequences(generate_sequences);
    }

    void add_sequences(std::vector<std::string>&& sequences) {
        auto seqs = std::make_shared<std::vector<std::string>>(std::move(sequences));
        constructor_->add_sequences(
            [seqs](const CallString &callback) {
                std::for_each(seqs->begin(), seqs->end(), callback);
            }
        );
    }

    void build_graph(BOSS *graph) {
        auto chunk = constructor_->build_chunk();
        // initialize graph from the chunk built
        chunk->initialize_boss(graph);
        delete chunk;
    }

    void build_graph(BOSS *graph, sdsl::int_vector<> *weights) {
        auto chunk = constructor_->build_chunk();
        // initialize graph from the chunk built
        chunk->initialize_boss(graph, weights);
        delete chunk;
    }

    uint64_t get_k() const { return constructor_->get_k(); }

    static BOSS* build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                                         bool verbose = false,
                                         sdsl::int_vector<> *weights = nullptr) {
        return BOSS::Chunk::build_boss_from_chunks(chunk_filenames, verbose, weights).first;
    }

  private:
    std::unique_ptr<IBOSSChunkConstructor> constructor_;
};

} // namespace boss
} // namespace graph
} // namespace mtg

#endif // __BOSS_CONSTRUCT_HPP__
