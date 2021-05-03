#include "aligner_prefix_suffix.hpp"

#include "graph/representation/base/sequence_graph.hpp"

namespace mtg {
namespace graph {
namespace align {


template <typename NodeType>
void AlignmentSuffix<NodeType>::shift_path_ptr_right() {
    while (!eof() && *cigar_it_ == Cigar::DELETION) {
        operator++();
        assert(Alignment<NodeType>(*this).is_valid(*graph_, config_));
    }
}

template <typename NodeType>
void AlignmentSuffix<NodeType>::remove_query(size_t n) {
    for (size_t i = 0; i < n; ++i) {
        assert(!eof());
        shift_path_ptr_right();

        assert(!eof());
        operator++();
        assert(Alignment<NodeType>(*this).is_valid(*graph_, config_));

        shift_path_ptr_right();
    }
}


template <typename NodeType>
AlignmentSuffix<NodeType>& AlignmentSuffix<NodeType>::operator++() {
    assert(cigar_it_ != cigar_end_);

    switch (*cigar_it_) {
        case Cigar::MATCH:
        case Cigar::MISMATCH: {
            assert(query_.size());
            assert(ref_.size());
            score_ -= config_->get_row(query_[0])[ref_[0]];
            query_.remove_prefix(1);
            ref_.remove_prefix(1);
            if (node_it_ + 1 != alignment_->end()) {
                ++node_it_;
            } else {
                ++added_offset_;
            }
        } break;
        case Cigar::INSERTION: {
            assert(query_.size());
            if (cigar_it_.count_left() == 1) {
                score_ -= config_->gap_opening_penalty;
            } else {
                score_ -= config_->gap_extension_penalty;
            }

            query_.remove_prefix(1);
        } break;
        case Cigar::DELETION: {
            assert(ref_.size());
            if (cigar_it_.count_left() == 1) {
                score_ -= config_->gap_opening_penalty;
            } else {
                score_ -= config_->gap_extension_penalty;
            }

            ref_.remove_prefix(1);
            if (node_it_ + 1 != alignment_->end()) {
                ++node_it_;
            } else {
                ++added_offset_;
            }
        } break;
        case Cigar::CLIPPED: { assert(false); }
    }

    ++cigar_it_;
    ++trim_;

    assert(cigar_begin_ <= cigar_it_);
    assert(cigar_it_ <= cigar_end_);

    assert(Alignment<NodeType>(*this).is_valid(*graph_, config_));

    return *this;
}

template <typename NodeType>
AlignmentSuffix<NodeType>& AlignmentSuffix<NodeType>::operator--() {
    assert(cigar_it_ != cigar_begin_);

    --cigar_it_;

    switch (*cigar_it_) {
        case Cigar::MATCH:
        case Cigar::MISMATCH: {
            assert(query_.data() > alignment_->get_query().data());
            assert(ref_.data() > alignment_->get_sequence().data());
            query_ = { query_.data() - 1, query_.size() + 1 };
            ref_ = { ref_.data() - 1, ref_.size() + 1 };
            if (added_offset_) {
                --added_offset_;
            } else {
                assert(node_it_ != alignment_->begin());
                --node_it_;
            }
            score_ += config_->get_row(query_[0])[ref_[0]];
        } break;
        case Cigar::INSERTION: {
            assert(query_.data() > alignment_->get_query().data());
            if (cigar_it_.count_left() == 1) {
                score_ += config_->gap_opening_penalty;
            } else {
                score_ += config_->gap_extension_penalty;
            }

            query_ = { query_.data() - 1, query_.size() + 1 };
        } break;
        case Cigar::DELETION: {
            assert(cigar_it_.get_it() >= alignment_->get_cigar().begin());
            assert(ref_.data() > alignment_->get_sequence().data());
            if (cigar_it_.count_left() == 1) {
                score_ += config_->gap_opening_penalty;
            } else {
                score_ += config_->gap_extension_penalty;
            }

            ref_ = { ref_.data() - 1, ref_.size() + 1 };
            if (added_offset_) {
                --added_offset_;
            } else {
                assert(node_it_ != alignment_->begin());
                --node_it_;
            }
        } break;
        case Cigar::CLIPPED: { assert(false); }
    }

    --trim_;

    assert(cigar_begin_ <= cigar_it_);
    assert(cigar_it_ <= cigar_end_);

    assert(Alignment<NodeType>(*this).is_valid(*graph_, config_));

    return *this;
}

template class AlignmentSuffix<>;

} // namespace align
} // namespace graph
} // namespace mtg
