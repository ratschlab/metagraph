#ifndef __PATH_HPP__
#define __PATH_HPP__

#include <vector>

template <typename NodeType, typename VLabels>
class Path {

  public:
    Path(std::string::const_iterator query_it) : query_it_(query_it) {
        loss_ = 0.0; }

    Path(const Path&) = default;
    Path(Path&&) = default;
    Path& operator= (const Path&) = default;
    Path& operator= (Path&&) = default;

    void push_back(NodeType node, const VLabels &labels, float error = 0) {
        nodes_.push_back(node);
        loss_ += error;
        label_set_.insert(std::end(label_set_), std::begin(labels),
                          std::end(labels));
        ++ query_it_;
    }
    void append_path(const Path& other) {
        for (auto node : other.nodes_)
            nodes_.push_back(node);
        loss_ += other.loss_;
        label_set_.insert(std::end(label_set_),
                          std::begin(other.label_set_),
                          std::end(other.label_set_));
        query_it_ = other.query_it_;
    }

    void update_total_loss(float loss) { loss_ = loss; }
    void set_sequence_it(std::string::const_iterator sequence_it) {
        query_it_ = sequence_it; }

    NodeType back() const { return nodes_.back(); }
    NodeType last_parent() const { return nodes_.at(nodes_.size() - 1); }
    size_t size() const { return nodes_.size(); }

    float get_total_loss() const { return loss_; }
    VLabels get_labels() const { return label_set_; }
    std::vector<NodeType> get_nodes() const { return nodes_; }

    std::string::const_iterator get_sequence_it() const { return query_it_; }

    bool operator< (const Path& other) const {
            return !(this->loss_/this->size() < other.loss_/other.size());
    }

  private:
    float loss_;
    std::vector<NodeType> nodes_;
    VLabels label_set_;
    std::string::const_iterator query_it_;
};

#endif  // __PATH_HPP__
