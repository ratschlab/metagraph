#include <vector>


template <typename NodeType, typename VLabels>
class PathCompare;

template <typename NodeType, typename VLabels>
class Path {

  public:
    Path() { loss_ = 0.0; }

    Path(const Path&) = default;
    Path(Path&&) = default;
    Path& operator= (const Path&) = default;
    Path& operator= (Path&&) = default;

    void push_back(NodeType node, float error, VLabels labels) {
        nodes_.push_back(node);
        loss_ += error;
        label_set_.insert(std::end(label_set_), std::begin(labels),
                          std::end(labels));
    }
    NodeType back() { return nodes_.back(); }
    size_t size() { return nodes_.size(); }
    float get_total_loss() { return loss_; }
    VLabels get_labels() { return label_set_; }
    std::vector<NodeType> get_nodes() { return nodes_; }

    template <typename NType, typename VL> friend class PathCompare;

  private:
    float loss_;
    std::vector<NodeType> nodes_;
    VLabels label_set_;
};

template <typename NodeType, typename VLabels>
class PathCompare {
  public:
    bool operator() (const Path<NodeType, VLabels> &lhs, const Path<NodeType, VLabels> &rhs) const {
            return (lhs.loss_ > rhs.loss_);
    }
};
