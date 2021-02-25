#ifndef __ALIGNER_DP_TABLE_HPP__
#define __ALIGNER_DP_TABLE_HPP__

#include <tsl/hopscotch_map.h>

#include "aligner_config.hpp"
#include "aligner_cigar.hpp"
#include "common/aligned_vector.hpp"
#include "common/utils/template_utils.hpp"


namespace mtg {
namespace graph {

class DeBruijnGraph;

namespace align {

template <typename NodeType>
class Alignment;

// dynamic programming table stores score columns and steps needed to reconstruct paths
template <typename NodeType = uint64_t>
class DPTable {
  public:
    typedef DBGAlignerConfig::score_t score_t;

    struct Column {
        Column() = default;

        Column(size_t size,
               score_t min_score,
               char start_char,
               size_t pos = 0,
               size_t priority_pos = 0,
               size_t start = 0,
               size_t end = std::numeric_limits<size_t>::max())
              : size_(size),
                min_score_(min_score),
                scores(std::min(end, size) - start + 8, min_score),
                gap_scores(scores.size(), min_score),
                ops(scores.size(), Cigar::CLIPPED),
                prev_nodes(scores.size(), 0),
                gap_prev_nodes(scores.size(), 0),
                gap_count(scores.size(), 0),
                last_char(start_char),
                best_pos(std::min(std::max(pos, start), start + scores.size() - (size_t)9)),
                last_priority_pos(std::min(std::max(priority_pos, start), start + scores.size() - (size_t)9)),
                start_index(start) {}

        size_t size_;
        score_t min_score_;
        AlignedVector<score_t> scores;
        AlignedVector<score_t> gap_scores;
        AlignedVector<Cigar::Operator> ops;
        AlignedVector<uint8_t> prev_nodes;
        AlignedVector<uint8_t> gap_prev_nodes;
        AlignedVector<int32_t> gap_count;
        mutable std::vector<NodeType> incoming_nodes;
        char last_char;
        size_t best_pos;
        size_t last_priority_pos;
        size_t start_index;

        void expand_to_cover(size_t begin, size_t end);

        const score_t& best_score() const { return scores.at(best_pos - start_index); }
        const score_t& last_priority_value() const { return scores.at(last_priority_pos - start_index); }
        const Cigar::Operator& best_op() const { return ops.at(best_pos - start_index); }

        bool operator<(const Column &other) const {
            return best_score() < other.best_score();
        }

        size_t cell_size() const {
            return sizeof(score_t) * 2 + sizeof(Cigar::Operator)
                 + sizeof(uint8_t) * 2 + sizeof(int32_t);
        }

        size_t bytes_taken() const {
            return sizeof(Column)
                + sizeof(score_t) * scores.capacity()
                + sizeof(score_t) * gap_scores.capacity()
                + sizeof(Cigar::Operator) * ops.capacity()
                + sizeof(uint8_t) * prev_nodes.capacity()
                + sizeof(uint8_t) * gap_prev_nodes.capacity()
                + sizeof(int32_t) * gap_count.capacity()
                + sizeof(NodeType) * incoming_nodes.capacity();
        }

        size_t size() const { return size_; }

        uint8_t rank_prev_node(NodeType node) const {
            assert(incoming_nodes.size() < std::numeric_limits<uint8_t>::max());
            for (uint8_t i = 0; i < incoming_nodes.size(); ++i) {
                if (incoming_nodes[i] == node)
                    return i + 1;
            }

            incoming_nodes.push_back(node);
            return incoming_nodes.size();
        }

        NodeType select_prev_node(uint8_t rank) const {
            assert(rank);
            return incoming_nodes.at(rank - 1);
        }
    };

    DPTable() {}

    bool add_seed(const Alignment<NodeType> &seed,
                  const DBGAlignerConfig &config,
                  size_t size,
                  size_t start_pos);

    typedef tsl::hopscotch_map<NodeType, Column> Storage;

    typedef NodeType key_type;
    typedef Column mapped_type;
    typedef typename Storage::value_type value_type;

    typedef typename Storage::iterator iterator;
    typedef typename Storage::const_iterator const_iterator;

    void expand_to_cover(iterator it, size_t begin, size_t end) {
        size_t old_size = it->second.bytes_taken();
        it.value().expand_to_cover(begin, end);
        num_bytes_ += it->second.bytes_taken() - old_size;
    }

    iterator begin() { return dp_table_.begin(); }
    iterator end() { return dp_table_.end(); }
    const_iterator begin() const { return dp_table_.begin(); }
    const_iterator end() const { return dp_table_.end(); }

    iterator find(const NodeType &node) { return dp_table_.find(node); }
    const_iterator find(const NodeType &node) const { return dp_table_.find(node); }

    size_t size() const { return dp_table_.size(); }

    size_t num_bytes() const { return num_bytes_; }

    void clear() {
        dp_table_.clear();
        num_bytes_ = 0;
    }

    template <typename... Args>
    std::pair<iterator, bool> emplace(Args&&... args) {
        auto pair = dp_table_.emplace(std::forward<Args>(args)...);

        if (pair.second)
            num_bytes_ += pair.first.value().bytes_taken();

        return pair;
    }

    void erase(NodeType key) { dp_table_.erase(key); }
    size_t count(NodeType key) const { return dp_table_.count(key); }

    void extract_alignments(const DeBruijnGraph &graph,
                            const DBGAlignerConfig &config,
                            std::string_view query_view,
                            std::function<void(Alignment<NodeType>&&, NodeType)> callback,
                            score_t min_path_score,
                            const Alignment<NodeType> &seed,
                            NodeType *node = nullptr);

    std::pair<NodeType, score_t> best_score() const {
        auto mx = std::max_element(begin(), end(), utils::LessSecond());
        return std::make_pair(mx->first, mx->second.best_score());
    }

    const Storage& data() const { return dp_table_; }
    NodeType get_start_node() const { return start_node_; }

  private:
    Storage dp_table_;
    NodeType start_node_;
    size_t num_bytes_ = 0;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_DP_TABLE_HPP__
