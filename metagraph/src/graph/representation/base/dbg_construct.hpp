#ifndef __DBG_CONSTRUCT_HPP__
#define __DBG_CONSTRUCT_HPP__

#include <functional>
#include <string>

typedef std::function<void(const std::string&)> CallString;


template <class GraphChunk>
class IGraphChunkConstructor {
  public:
    virtual ~IGraphChunkConstructor() {}

    virtual void add_sequence(std::string_view sequence, uint64_t count = 1) = 0;
    virtual void add_sequences(std::function<void(CallString)> generate_sequences) = 0;

    virtual GraphChunk* build_chunk() = 0;
};


template <class Graph>
class IGraphConstructor {
  public:
    virtual ~IGraphConstructor() {}

    virtual void add_sequence(std::string_view sequence, uint64_t count = 1) = 0;
    virtual void add_sequences(std::function<void(CallString)> generate_sequences) = 0;

    virtual void build_graph(Graph *graph) = 0;
};

#endif // __DBG_CONSTRUCT_HPP__
