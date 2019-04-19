#ifndef __DBG_BITMAP_CONSTRUCT_HPP__
#define __DBG_BITMAP_CONSTRUCT_HPP__

#include "dbg_bitmap.hpp"
#include "dbg_construct.hpp"


class IBitmapChunkConstructor : public IGraphChunkConstructor<DBGBitmap::Chunk> {
  public:
    virtual ~IBitmapChunkConstructor() {}

    static IBitmapChunkConstructor* initialize(size_t k,
                                               bool canonical_mode = false,
                                               const std::string &filter_suffix = "",
                                               size_t num_threads = 1,
                                               double memory_preallocated = 0,
                                               bool verbose = false);

    virtual void add_sequence(const std::string &sequence) = 0;
    virtual void add_sequences(std::function<void(CallString)> generate_sequences) = 0;

    virtual DBGBitmap::Chunk* build_chunk() = 0;

    virtual size_t get_k() const = 0;
    virtual bool is_canonical_mode() const = 0;
};


class DBGBitmapConstructor : public IGraphConstructor<DBGBitmap> {
  public:
    DBGBitmapConstructor(size_t k,
                         bool canonical_mode = false,
                         const std::string &filter_suffix = "",
                         size_t num_threads = 1,
                         double memory_preallocated = 0,
                         bool verbose = false);

    void add_sequence(const std::string &sequence) {
        constructor_->add_sequence(sequence);
    }

    void add_sequences(const std::vector<std::string> &sequences) {
        constructor_->add_sequences([&sequences](const CallString &callback) {
            std::for_each(sequences.begin(), sequences.end(), callback);
        });
    }

    void add_sequences(std::function<void(CallString)> generate_sequences) {
        constructor_->add_sequences(generate_sequences);
    }

    void build_graph(DBGBitmap *graph);
    DBGBitmap::Chunk* build_chunk() { return constructor_->build_chunk(); }

    static DBGBitmap* build_graph_from_chunks(const std::vector<std::string> &chunk_filenames,
                                              bool canonical_mode = false,
                                              bool verbose = false);
  private:
    std::unique_ptr<IBitmapChunkConstructor> constructor_;
};

#endif // __DBG_BITMAP_CONSTRUCT_HPP__
