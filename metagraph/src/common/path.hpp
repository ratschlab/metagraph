#ifndef __PATH_HPP__
#define __PATH_HPP__

#include <vector>

template <typename NodeType, typename VLabels>
class Path {

  public:
    Path(std::string::const_iterator mapped_it) : mapped_it_(mapped_it) {
        loss_ = 0.0; }

    Path(const Path&) = default;
    Path(Path&&) = default;
    Path& operator= (const Path&) = default;
    Path& operator= (Path&&) = default;

    void push_back(NodeType node, float error, const VLabels &labels) {
        nodes_.push_back(node);
        loss_ += error;
        label_set_.insert(std::end(label_set_), std::begin(labels),
                          std::end(labels));
        ++ mapped_it_;
    }
    void set_sequence_it(std::string::const_iterator sequence_it) {
        mapped_it_ = sequence_it; }

    NodeType back() const { return nodes_.back(); }
    size_t size() const { return nodes_.size(); }
    float get_total_loss() const { return loss_; }
    VLabels get_labels() const { return label_set_; }
    std::vector<NodeType> get_nodes() const { return nodes_; }
    std::string::const_iterator get_sequence_it() const { return mapped_it_; }

    bool operator< (const Path& other) const {
            return (loss_ > other.loss_);
    }

    template <typename NType, typename VL> friend class PathCompare;

  private:
    float loss_;
    std::vector<NodeType> nodes_;
    VLabels label_set_;
    std::string::const_iterator mapped_it_;
};

#endif  // __PATH_HPP__
