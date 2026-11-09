#ifndef __NODE_STREAM_HPP__
#define __NODE_STREAM_HPP__

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <span>
#include <vector>


namespace mtg {
namespace graph {

// Pull-based mapping of the k-mers of a sequence to graph nodes: one node per
// k-mer position, in input order, npos (0) for k-mers missing from the graph.
// The stream keeps a reference to the sequence; it must outlive the stream.
class NodeStream {
  public:
    typedef uint64_t node_index;

    virtual ~NodeStream() {}

    // Write the nodes of the next k-mers to |out| and return their number.
    // Returns less than out.size() only when the sequence is exhausted.
    virtual size_t next(std::span<node_index> out) = 0;
};

// A stream over nodes that were already computed.
class BufferedNodeStream final : public NodeStream {
  public:
    explicit BufferedNodeStream(std::vector<node_index>&& nodes) : nodes_(std::move(nodes)) {}

    size_t next(std::span<node_index> out) override {
        size_t n = std::min(out.size(), nodes_.size() - pos_);
        std::copy(nodes_.begin() + pos_, nodes_.begin() + pos_ + n, out.begin());
        pos_ += n;
        return n;
    }

  private:
    std::vector<node_index> nodes_;
    size_t pos_ = 0;
};

// Single-pass input range over a NodeStream. Reads the stream in small chunks
// to amortize the virtual call. Example:
//     for (node_index node : graph.node_range(sequence)) { ... }
class NodeRange {
  public:
    typedef NodeStream::node_index node_index;
    static constexpr size_t kChunkSize = 64;

    explicit NodeRange(std::unique_ptr<NodeStream> stream) : stream_(std::move(stream)) {}

    class iterator {
      public:
        using iterator_concept = std::input_iterator_tag;
        using value_type = node_index;
        using difference_type = std::ptrdiff_t;

        iterator() = default;
        explicit iterator(NodeRange *range) : range_(range) {}

        node_index operator*() const { return range_->buffer_[range_->pos_]; }
        iterator& operator++() { ++range_->pos_; return *this; }
        void operator++(int) { ++*this; }

        friend bool operator==(const iterator &it, std::default_sentinel_t) {
            return it.range_->pos_ == it.range_->size_ && !it.range_->refill();
        }

      private:
        NodeRange *range_ = nullptr;
    };

    iterator begin() { return iterator(this); }
    std::default_sentinel_t end() const { return {}; }

  private:
    std::unique_ptr<NodeStream> stream_;
    std::array<node_index, kChunkSize> buffer_;
    size_t pos_ = 0;
    size_t size_ = 0;

    bool refill() {
        pos_ = 0;
        size_ = stream_ ? stream_->next(buffer_) : 0;
        if (size_ < kChunkSize)
            stream_.reset();
        return size_;
    }
};

static_assert(std::input_iterator<NodeRange::iterator>);
static_assert(std::sentinel_for<std::default_sentinel_t, NodeRange::iterator>);

} // namespace graph
} // namespace mtg

#endif // __NODE_STREAM_HPP__
