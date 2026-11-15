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

// A stream of |size| nodes, where |get_node(i)| computes the i-th one.
// |get_node| is called once for each i, in increasing order.
template <class GetNode>
std::unique_ptr<NodeStream> make_node_stream(size_t size, GetNode get_node) {
    class Stream final : public NodeStream {
      public:
        Stream(size_t size, GetNode&& get_node) : size_(size), get_node_(std::move(get_node)) {}

        size_t next(std::span<node_index> out) override {
            size_t n = std::min(out.size(), size_ - i_);
            for (size_t j = 0; j < n; ++j) {
                out[j] = get_node_(i_++);
            }
            return n;
        }

      private:
        size_t size_;
        GetNode get_node_;
        size_t i_ = 0;
    };
    return std::make_unique<Stream>(size, std::move(get_node));
}

// A stream that applies |map| to every node of |stream|.
template <class Map>
std::unique_ptr<NodeStream> map_node_stream(std::unique_ptr<NodeStream> stream, Map map) {
    class Stream final : public NodeStream {
      public:
        Stream(std::unique_ptr<NodeStream> stream, Map&& map)
              : stream_(std::move(stream)), map_(std::move(map)) {}

        size_t next(std::span<node_index> out) override {
            size_t n = stream_->next(out);
            std::transform(out.begin(), out.begin() + n, out.begin(), map_);
            return n;
        }

      private:
        std::unique_ptr<NodeStream> stream_;
        Map map_;
    };
    return std::make_unique<Stream>(std::move(stream), std::move(map));
}

// Single-pass input range over a NodeStream. Reads the stream in chunks of
// doubling size (1 to kChunkSize) to amortize the virtual call, so a loop that
// stops early maps at most as many nodes in vain as it has consumed. Example:
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
            return it.at_end();
        }

      private:
        NodeRange *range_ = nullptr;

        bool at_end() const { return range_->pos_ == range_->size_ && !range_->refill(); }
    };

    iterator begin() {
        if (pos_ == size_)
            refill();
        return iterator(this);
    }
    std::default_sentinel_t end() const { return {}; }

  private:
    std::unique_ptr<NodeStream> stream_;
    std::array<node_index, kChunkSize> buffer_;
    size_t pos_ = 0;
    size_t size_ = 0;
    size_t chunk_size_ = 1;

    bool refill() {
        pos_ = 0;
        size_ = stream_ ? stream_->next({ buffer_.data(), chunk_size_ }) : 0;
        if (size_ < chunk_size_)
            stream_.reset();
        chunk_size_ = std::min(2 * chunk_size_, kChunkSize);
        return size_;
    }
};

static_assert(std::input_iterator<NodeRange::iterator>);
static_assert(std::sentinel_for<std::default_sentinel_t, NodeRange::iterator>);

} // namespace graph
} // namespace mtg

#endif // __NODE_STREAM_HPP__
