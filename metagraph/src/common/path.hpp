#ifndef __PATH_HPP__
#define __PATH_HPP__

#include <vector>

template <typename NodeType, typename VLabels>
class Path {

  public:
    Path(std::string::const_iterator query_it) : query_it_(query_it) {
        score_ = 0.0; }

    Path(const Path&) = default;
    Path(Path&&) = default;
    Path& operator= (const Path&) = default;
    Path& operator= (Path&&) = default;

    void seed(NodeType node, const VLabels &labels, std::string sequence, float score=0) {
        push_back(node, labels, score);
        sequence_ += sequence;
    }
    void extend(NodeType node, const VLabels &labels, char extention, float score=0) {
        push_back(node, labels, score);
        sequence_ += extention;
    }
    // Assume paths are appended on the correct order.
    void append_path(const Path& other, int64_t k) {
        if (nodes_.size() > 0)
            sequence_ += other.sequence_.substr(k - 2);
        else
            sequence_ = other.sequence_;
        for (auto node : other.nodes_)
            nodes_.push_back(node);
        score_ += other.score_;
        label_set_.insert(std::end(label_set_),
                          std::begin(other.label_set_),
                          std::end(other.label_set_));
        query_it_ = other.query_it_;
    }

    void update_total_score(float score) { score_ = score; }
    void set_query_it(std::string::const_iterator query_it) {
        query_it_ = query_it; }

    NodeType back() const { return nodes_.back(); }
    NodeType last_parent() const { return nodes_.at(nodes_.size() - 1); }
    size_t size() const { return nodes_.size(); }

    float get_total_score() const { return score_; }
    VLabels get_labels() const { return label_set_; }
    std::vector<NodeType> get_nodes() const { return nodes_; }

    std::string::const_iterator get_query_it() const { return query_it_; }
    std::string get_sequence() const { return sequence_; }

    // The paths are sorted in BoundedPriorityQueue in increasing order of score
    // per number of nodes. This gives paths with lower absolute score, but higher
    // score per node to appear at the top of the queue.
    bool operator< (const Path& other) const {
        return (this->score_/float(this->size()) < other.score_/float(other.size()));
    }

  private:
    float score_;
    std::vector<NodeType> nodes_;
    VLabels label_set_;
    std::string::const_iterator query_it_;
    std::string sequence_;

    void push_back(NodeType node, const VLabels &labels, float score=0) {
        nodes_.push_back(node);
        score_ += score;
        label_set_.insert(std::end(label_set_), std::begin(labels),
                          std::end(labels));
        ++ query_it_;
    }
};

#endif  // __PATH_HPP__
