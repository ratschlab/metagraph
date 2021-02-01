#include "aligner_helper.hpp"

namespace mtg {
namespace graph {
namespace align {


template <typename NodeType>
AlignmentPrefix<NodeType>& AlignmentPrefix<NodeType>::operator++() {
    assert(cigar_it_ != cigar_rend_);

    switch (*cigar_it_) {
        case Cigar::MATCH:
        case Cigar::MISMATCH: {
            assert(end_it_ > begin_it_);
            assert(ref_end_it_ > ref_begin_it_);
            assert(node_it_ != alignment_->rend());
            --end_it_;
            --ref_end_it_;
            score_ -= config_->get_row(*end_it_)[*ref_end_it_];

            // if we're at the first node of the alignment, then traverse
            // backwards to reconstruct a path
            if (node_it_ + 1 == alignment_->rend()) {
                auto last_node = prefix_node_ ? prefix_node_ : alignment_->front();
                assert(last_node);
                prefix_node_ = DeBruijnGraph::npos;
                graph_->adjacent_incoming_nodes(last_node, [&](auto prev) {
                    // TODO: what if there are multiple?
                    if (!prefix_node_) {
                        prefix_node_ = prev;
                        assert(graph_->traverse(prev, *ref_end_it_) == last_node);
                    }
                });
                ++offset_;
            } else {
                ++node_it_;
            }
        } break;
        case Cigar::INSERTION: {
            assert(end_it_ > begin_it_);
            if ((--(cigar_it_.base())).offset()) {
                score_ -= config_->gap_extension_penalty;
            } else {
                score_ -= config_->gap_opening_penalty;
            }

            --end_it_;
        } break;
        case Cigar::DELETION: {
            assert((--(cigar_it_.base())).get_it() >= alignment_->get_cigar().begin());
            assert(ref_end_it_ > ref_begin_it_);
            assert(node_it_ != alignment_->rend());
            if ((--(cigar_it_.base())).offset()) {
                score_ -= config_->gap_extension_penalty;
            } else {
                score_ -= config_->gap_opening_penalty;
            }

            --ref_end_it_;

            // if we're at the first node of the alignment, then traverse
            // backwards to reconstruct a path
            if (node_it_ + 1 == alignment_->rend()) {
                auto last_node = prefix_node_ ? prefix_node_ : alignment_->front();
                assert(last_node);
                prefix_node_ = DeBruijnGraph::npos;
                graph_->adjacent_incoming_nodes(last_node, [&](auto prev) {
                    // TODO: what if there are multiple?
                    if (!prefix_node_) {
                        prefix_node_ = prev;
                        assert(graph_->traverse(prev, *ref_end_it_) == last_node);
                    }
                });
                ++offset_;

            } else {
                ++node_it_;
            }
        } break;
        case Cigar::MISSING: {
            assert(node_it_ != alignment_->rend());
            assert(!*node_it_);
            ++node_it_;
        } break;
        case Cigar::CLIPPED: { assert(false); }
    }

    ++cigar_it_;
    ++trim_;

    assert(cigar_rbegin_ <= cigar_it_);
    assert(cigar_it_ <= cigar_rend_);

    assert(Alignment<NodeType>(*this).is_valid(*graph_, config_));

    return *this;
}

template <typename NodeType>
AlignmentPrefix<NodeType>& AlignmentPrefix<NodeType>::operator--() {
    assert(cigar_it_ != cigar_rbegin_);

    --cigar_it_;

    switch (*cigar_it_) {
        case Cigar::MATCH:
        case Cigar::MISMATCH: {
            score_ += config_->get_row(*end_it_)[*ref_end_it_];
            ++end_it_;

            if (offset_) {
                assert(*node_it_ == alignment_->front());
                if (offset_ > 1) {
                    assert(*node_it_);
                    prefix_node_ = graph_->traverse(prefix_node_, *ref_end_it_);
                    assert(prefix_node_);
                }

                --offset_;
                assert(offset_ || graph_->traverse(prefix_node_, *ref_end_it_) == *node_it_);

                if (!offset_)
                    prefix_node_ = DeBruijnGraph::npos;
            } else {
                --node_it_;
            }

            ++ref_end_it_;
        } break;
        case Cigar::INSERTION: {
            ++cigar_it_;
            if (*cigar_it_ == Cigar::INSERTION) {
                score_ += config_->gap_extension_penalty;
            } else {
                score_ += config_->gap_opening_penalty;
            }
            --cigar_it_;

            ++end_it_;
        } break;
        case Cigar::DELETION: {
            ++cigar_it_;
            if (*cigar_it_ == Cigar::DELETION) {
                score_ += config_->gap_extension_penalty;
            } else {
                score_ += config_->gap_opening_penalty;
            }
            --cigar_it_;

            if (offset_) {
                assert(*node_it_ == alignment_->front());
                if (offset_ > 1) {
                    assert(*node_it_);
                    prefix_node_ = graph_->traverse(prefix_node_, *ref_end_it_);
                    assert(prefix_node_);
                }

                --offset_;
                assert(offset_ || graph_->traverse(prefix_node_, *ref_end_it_) == *node_it_);

                if (!offset_)
                    prefix_node_ = DeBruijnGraph::npos;
            } else {
                --node_it_;
            }
            ++ref_end_it_;
        } break;
        case Cigar::MISSING: {
            --end_it_;
        } break;
        case Cigar::CLIPPED: { assert(false); }
    }

    --trim_;

    assert(cigar_rbegin_ <= cigar_it_);
    assert(cigar_it_ <= cigar_rend_);

    assert(Alignment<NodeType>(*this).is_valid(*graph_, config_));

    return *this;
}

template <typename NodeType>
AlignmentSuffix<NodeType>& AlignmentSuffix<NodeType>::operator++() {
    assert(cigar_it_ != cigar_end_);

    switch (*cigar_it_) {
        case Cigar::MATCH:
        case Cigar::MISMATCH: {
            assert(begin_it_ < end_it_);
            assert(ref_begin_it_ < ref_end_it_);
            score_ -= config_->get_row(*begin_it_)[*ref_begin_it_];
            ++begin_it_;
            ++ref_begin_it_;
            if (node_it_ + 1 != alignment_->end()) {
                ++node_it_;
            } else {
                ++added_offset_;
            }
        } break;
        case Cigar::INSERTION: {
            assert(begin_it_ < end_it_);
            if (cigar_it_.count_left() == 1) {
                score_ -= config_->gap_opening_penalty;
            } else {
                score_ -= config_->gap_extension_penalty;
            }

            ++begin_it_;
        } break;
        case Cigar::DELETION: {
            assert(ref_begin_it_ < ref_end_it_);
            if (cigar_it_.count_left() == 1) {
                score_ -= config_->gap_opening_penalty;
            } else {
                score_ -= config_->gap_extension_penalty;
            }

            ++ref_begin_it_;
            if (node_it_ + 1 != alignment_->end()) {
                ++node_it_;
            } else {
                ++added_offset_;
            }
        } break;
        case Cigar::MISSING: {
            assert(!*node_it_);
            assert(node_it_ + 1 != alignment_->end());
            ++node_it_;
        } break;
        case Cigar::CLIPPED: { assert(false); }
    }

    ++cigar_it_;
    ++trim_;

    assert(cigar_begin_ <= cigar_it_);
    assert(cigar_it_ <= cigar_end_);

    return *this;
}

template <typename NodeType>
AlignmentSuffix<NodeType>& AlignmentSuffix<NodeType>::operator--() {
    assert(cigar_it_ != cigar_begin_);

    --cigar_it_;

    switch (*cigar_it_) {
        case Cigar::MATCH:
        case Cigar::MISMATCH: {
            assert(begin_it_ >= alignment_->get_query().data());
            assert(ref_begin_it_ >= alignment_->get_sequence().data());
            --begin_it_;
            --ref_begin_it_;
            if (added_offset_) {
                --added_offset_;
            } else {
                assert(node_it_ != alignment_->begin());
                --node_it_;
            }
            score_ += config_->get_row(*begin_it_)[*ref_begin_it_];
        } break;
        case Cigar::INSERTION: {
            assert(begin_it_ >= alignment_->get_query().data());
            if (cigar_it_.count_left() == 1) {
                score_ += config_->gap_opening_penalty;
            } else {
                score_ += config_->gap_extension_penalty;
            }

            --begin_it_;
        } break;
        case Cigar::DELETION: {
            assert(cigar_it_.get_it() >= alignment_->get_cigar().begin());
            assert(ref_begin_it_ >= alignment_->get_sequence().data());
            if (cigar_it_.count_left() == 1) {
                score_ += config_->gap_opening_penalty;
            } else {
                score_ += config_->gap_extension_penalty;
            }

            --ref_begin_it_;
            if (added_offset_) {
                --added_offset_;
            } else {
                assert(node_it_ != alignment_->begin());
                --node_it_;
            }
        } break;
        case Cigar::MISSING: {
            assert(!*node_it_);
            assert(node_it_ != alignment_->begin());
            --node_it_;
        } break;
        case Cigar::CLIPPED: { assert(false); }
    }

    --trim_;

    assert(cigar_begin_ <= cigar_it_);
    assert(cigar_it_ <= cigar_end_);

    return *this;
}

template class AlignmentPrefix<>;
template class AlignmentSuffix<>;

} // namespace align
} // namespace graph
} // namespace mtg
