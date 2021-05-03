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

} // namespace align
} // namespace graph
} // namespace mtg

#endif // __ALIGNER_PREFIX_SUFFIX_HPP__
