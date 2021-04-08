#include "aligner_prefix_suffix.hpp"

#include "graph/representation/base/sequence_graph.hpp"

namespace mtg {
namespace graph {
namespace align {


template <typename NodeType>
AlignmentPrefix<NodeType>& AlignmentPrefix<NodeType>::operator++() {
    assert(cigar_it_ != cigar_rend_);

    switch (*cigar_it_) {
        case Cigar::MATCH:
        case Cigar::MISMATCH: {
            assert(query_.size());
            assert(ref_.size());
            assert(node_it_ != alignment_->rend());
            score_ -= config_->get_row(query_.back())[ref_.back()];

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
                        assert(graph_->traverse(prev, ref_.back()) == last_node);
                    }
                });
                ++offset_;
            } else {
                ++node_it_;
            }

            query_.remove_suffix(1);
            ref_.remove_suffix(1);
        } break;
        case Cigar::INSERTION: {
            assert(query_.size());
            if ((--(cigar_it_.base())).offset()) {
                score_ -= config_->gap_extension_penalty;
            } else {
                score_ -= config_->gap_opening_penalty;
            }

            query_.remove_suffix(1);
        } break;
        case Cigar::DELETION: {
            assert((--(cigar_it_.base())).get_it() >= alignment_->get_cigar().begin());
            assert(ref_.size());
            assert(node_it_ != alignment_->rend());
            if ((--(cigar_it_.base())).offset()) {
                score_ -= config_->gap_extension_penalty;
            } else {
                score_ -= config_->gap_opening_penalty;
            }

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
                        assert(graph_->traverse(prev, ref_.back()) == last_node);
                    }
                });
                ++offset_;

            } else {
                ++node_it_;
            }

            ref_.remove_suffix(1);

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
            ref_ = { ref_.data(), ref_.size() + 1 };
            query_ = { query_.data(), query_.size() + 1 };
            score_ += config_->get_row(query_.back())[ref_.back()];

            if (offset_) {
                assert(*node_it_ == alignment_->front());
                if (offset_ > 1) {
                    assert(*node_it_);
                    prefix_node_ = graph_->traverse(prefix_node_, ref_.back());
                    assert(prefix_node_);
                }

                --offset_;
                assert(offset_ || graph_->traverse(prefix_node_, ref_.back()) == *node_it_);

                if (!offset_)
                    prefix_node_ = DeBruijnGraph::npos;
            } else {
                --node_it_;
            }

        } break;
        case Cigar::INSERTION: {
            ++cigar_it_;
            if (*cigar_it_ == Cigar::INSERTION) {
                score_ += config_->gap_extension_penalty;
            } else {
                score_ += config_->gap_opening_penalty;
            }
            --cigar_it_;

            query_ = { query_.data(), query_.size() + 1 };
        } break;
        case Cigar::DELETION: {
            ++cigar_it_;
            if (*cigar_it_ == Cigar::DELETION) {
                score_ += config_->gap_extension_penalty;
            } else {
                score_ += config_->gap_opening_penalty;
            }
            --cigar_it_;

            ref_ = { ref_.data(), ref_.size() + 1 };

            if (offset_) {
                assert(*node_it_ == alignment_->front());
                if (offset_ > 1) {
                    assert(*node_it_);
                    prefix_node_ = graph_->traverse(prefix_node_, ref_.back());
                    assert(prefix_node_);
                }

                --offset_;
                assert(offset_ || graph_->traverse(prefix_node_, ref_.back()) == *node_it_);

                if (!offset_)
                    prefix_node_ = DeBruijnGraph::npos;
            } else {
                assert(node_it_ != alignment_->rbegin());
                --node_it_;
            }

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
void AlignmentPrefix<NodeType>::shift_path_ptr_left() {
    while (!reof() && *cigar_it_ == Cigar::DELETION) {
        operator--();
        assert(Alignment<NodeType>(*this).is_valid(*graph_, config_));
    }
}

template <typename NodeType>
void AlignmentPrefix<NodeType>::extend_query(size_t n) {
    for (size_t i = 0; i < n; ++i) {
        assert(!reof());
        shift_path_ptr_left();

        assert(!reof());
        operator--();
        assert(Alignment<NodeType>(*this).is_valid(*graph_, config_));

        shift_path_ptr_left();
    }
}

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

template class AlignmentPrefix<>;
template class AlignmentSuffix<>;

} // namespace align
} // namespace graph
} // namespace mtg
