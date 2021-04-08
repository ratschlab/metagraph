#ifndef __ALIGNER_PREFIX_SUFFIX_HPP__
#define __ALIGNER_PREFIX_SUFFIX_HPP__

#include "aligner_alignment.hpp"

namespace mtg {
namespace graph {

class DeBruijnGraph;

namespace align {


template <typename NodeType>
class AlignmentSuffix {
  public:
    typedef typename Alignment<NodeType>::score_t score_t;

    AlignmentSuffix(const Alignment<NodeType> &alignment,
                    const DBGAlignerConfig &config,
                    const DeBruijnGraph &graph)
          : alignment_(&alignment),
            graph_(&graph),
            config_(&config),
            query_(alignment_->get_query()),
            ref_(alignment_->get_sequence()),
            cigar_begin_(alignment_->get_cigar(), alignment_->get_clipping()),
            cigar_end_(alignment_->get_cigar(),
                       alignment_->get_cigar().end()
                          - static_cast<bool>(alignment_->get_end_clipping())),
            cigar_it_(cigar_begin_),
            score_(alignment_->get_score()),
            trim_(0),
            added_offset_(0),
            node_it_(alignment_->begin()) {
        assert(cigar_begin_ <= cigar_end_);
        assert(cigar_it_ <= cigar_end_);
    }

    AlignmentSuffix& operator++();
    AlignmentSuffix& operator--();

    score_t get_score() const { return score_; }
    std::string_view get_query() const { return query_; }

    std::string_view get_sequence() const { return ref_; }

    Cigar::Operator get_front_op() const { return *cigar_it_; }
    const CigarOpIterator& get_op_it() const {
        assert(cigar_begin_ <= cigar_it_);
        assert(cigar_it_ <= cigar_end_);
        return cigar_it_;
    }
    typename std::vector<NodeType>::const_iterator get_node_begin_it() const { return node_it_; }
    size_t get_added_offset() const { return added_offset_; }

    bool eof() const { return cigar_it_ == cigar_end_; }
    bool reof() const { return cigar_it_ == cigar_begin_; }

    const Alignment<NodeType>& data() const { return *alignment_; }

    size_t get_trim() const { return trim_; }

    const DeBruijnGraph& get_graph() const { return *graph_; }

    void remove_query(size_t n);

    // step through all DELETION operations while pointing to the same query character
    void shift_path_ptr_right();

  private:
    const Alignment<NodeType> *alignment_;
    const DeBruijnGraph *graph_;
    const DBGAlignerConfig *config_;
    std::string_view query_;
    std::string_view ref_;
    CigarOpIterator cigar_begin_;
    CigarOpIterator cigar_end_;
    CigarOpIterator cigar_it_;
    score_t score_;
    size_t trim_;
    size_t added_offset_;
    typename std::vector<NodeType>::const_iterator node_it_;
};

template <typename NodeType>
std::ostream& operator<<(std::ostream &out, const AlignmentSuffix<NodeType> &suffix) {
    out << Alignment<NodeType>(suffix);
    return out;
}

template <typename NodeType>
class AlignmentPrefix {
  public:
    typedef typename Alignment<NodeType>::score_t score_t;

    AlignmentPrefix(const Alignment<NodeType> &alignment,
                    const DBGAlignerConfig &config,
                    const DeBruijnGraph &graph)
          : alignment_(&alignment),
            graph_(&graph),
            config_(&config),
            query_(alignment_->get_query()),
            ref_(alignment_->get_sequence()),
            cigar_rbegin_(std::make_reverse_iterator(
                CigarOpIterator(alignment_->get_cigar(),
                                alignment_->get_cigar().end()))),
            cigar_rend_(std::make_reverse_iterator(
                CigarOpIterator(alignment_->get_cigar()))),
            cigar_it_(cigar_rbegin_),
            score_(alignment_->get_score()),
            trim_(0), offset_(0), prefix_node_(0),
            node_it_(alignment_->rbegin()) {
        cigar_rbegin_ += alignment_->get_end_clipping();
        cigar_rend_ -= alignment_->get_clipping();
        cigar_it_ = cigar_rbegin_;

        assert(cigar_rbegin_ <= cigar_rend_);
        assert(cigar_it_ <= cigar_rend_);
    }

    AlignmentPrefix& operator++();
    AlignmentPrefix& operator--();

    score_t get_score() const { return score_; }
    std::string_view get_query() const { return query_; }

    std::string_view get_sequence() const { return ref_; }
    Cigar::Operator get_back_op() const { return *cigar_it_; }
    typename std::vector<NodeType>::const_reverse_iterator get_node_end_it() const {
        return node_it_;
    }

    bool eof() const {
        return cigar_it_ == cigar_rend_ || (offset_ && !prefix_node_);
    }
    bool reof() const { return cigar_it_ == cigar_rbegin_; }

    const Alignment<NodeType>& data() const { return *alignment_; }

    size_t get_trim() const { return trim_; }
    size_t get_offset() const { return offset_; }
    NodeType get_prefix_node() const { return prefix_node_; }

    const DeBruijnGraph& get_graph() const { return *graph_; }

    void extend_query(size_t n);

    // step through all DELETION operations while pointing to the same query character
    void shift_path_ptr_left();

  private:
    const Alignment<NodeType> *alignment_;
    const DeBruijnGraph *graph_;
    const DBGAlignerConfig *config_;
    std::string_view query_;
    std::string_view ref_;
    std::reverse_iterator<CigarOpIterator> cigar_rbegin_;
    std::reverse_iterator<CigarOpIterator> cigar_rend_;
    std::reverse_iterator<CigarOpIterator> cigar_it_;
    score_t score_;
    size_t trim_;
    size_t offset_;
    NodeType prefix_node_;
    typename std::vector<NodeType>::const_reverse_iterator node_it_;
};

template <typename NodeType>
std::ostream& operator<<(std::ostream &out, const AlignmentPrefix<NodeType> &prefix) {
    out << Alignment<NodeType>(prefix);
    return out;
}

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_PREFIX_SUFFIX_HPP__
