#ifndef __PATH_HPP__
#define __PATH_HPP__

#include <vector>

template <typename NodeType, typename VLabels>
class Path {

  public:
    Path(std::string::const_iterator query_begin_it,
         std::string::const_iterator query_it) :
            query_begin_it_(query_begin_it),
            query_it_(query_it) {
        score_ = 0.0;
        is_similar_ = false; }

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

    void trim_with_score(uint64_t trim_length, float total_score) {
        nodes_.resize(nodes_.size() - trim_length);
        score_ = total_score;
        query_it_ -= trim_length;
        sequence_.resize(sequence_.size() - trim_length);
        // TODO: remove labels from trimmed nodes.
    }

    void update_total_score(float score) { score_ = score; }
    void set_query_begin_it(std::string::const_iterator query_begin_it) {
        query_begin_it_ = query_begin_it; }
    void set_query_it(std::string::const_iterator query_it) {
        query_it_ = query_it; }
    void set_cigar(const std::string& cigar) { cigar_ = cigar; }
    void set_similar() { is_similar_ = true; }

    NodeType back() const { return nodes_.back(); }
    NodeType last_parent() const { return nodes_.at(nodes_.size() - 1); }
    size_t size() const { return nodes_.size(); }

    float get_total_score() const { return score_; }
    VLabels get_labels() const { return label_set_; }
    std::vector<NodeType> get_nodes() const { return nodes_; }

    std::string::const_iterator get_query_it() const { return query_it_; }
    std::string::const_iterator get_query_begin_it() const { return query_begin_it_; }
    std::string get_sequence() const { return sequence_; }
    std::string get_cigar() const { return cigar_; }
    bool get_similar() const { return is_similar_; }

    // The paths are sorted in BoundedPriorityQueue in increasing order of score
    // per number of nodes. This gives paths with lower absolute score, but higher
    // score per node to appear at the top of the queue.
    bool operator< (const Path &other) const {
        return (this->score_/this->size() < other.score_/other.size());
        //return (this->score_ < other.score_);
    }

  private:
    float score_;
    std::vector<NodeType> nodes_;
    VLabels label_set_;
    std::string::const_iterator query_begin_it_;
    std::string::const_iterator query_it_;
    std::string sequence_;
    std::string cigar_;
    bool is_similar_;

    void push_back(NodeType node, const VLabels &labels, float score=0) {
        nodes_.push_back(node);
        score_ += score;
        label_set_.insert(std::end(label_set_), std::begin(labels),
                          std::end(labels));
        ++ query_it_;
    }
};

#endif  // __PATH_HPP__
