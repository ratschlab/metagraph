#include "alignment.hpp"

#include "annotation_buffer.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "common/algorithms.hpp"
#include "common/logger.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/rc_dbg.hpp"


namespace mtg {
namespace graph {
namespace align {

using mtg::common::logger;

const Vector<Seed::Column> Seed::no_labels_ { std::numeric_limits<Seed::Column>::max() };

std::string Alignment::format_coords() const {
    if (!label_coordinates.size())
        return "";

    assert(label_coordinates.size() == get_columns(0).size());

    std::vector<std::string> decoded_labels = get_decoded_labels(0);
    for (size_t i = 0; i < decoded_labels.size(); ++i) {
        for (uint64_t coord : label_coordinates[i]) {
            // alignment coordinates are 1-based inclusive ranges
            decoded_labels[i]
                += fmt::format(":{}-{}", coord + 1, coord + sequence_.size());
        }
    }

    return fmt::format("{}", fmt::join(decoded_labels, ";"));
}

std::string Alignment::format_annotations() const {
    assert(has_annotation());
    std::string out = fmt::format("{}", fmt::join(get_decoded_labels(0), ";"));
    size_t count = 1;
    size_t last_cols = label_columns;
    for (size_t i = 0; i < label_column_diffs.size(); ++i) {
        if (label_column_diffs[i] == last_cols) {
            ++count;
        } else {
            out += fmt::format(":{}>{}", count, fmt::join(get_decoded_labels(i + 1), ";"));
            last_cols = label_column_diffs[i];
            count = 1;
        }
    }

    if (label_column_diffs.size())
        out += fmt::format(":{}", count);

    return out;
}

void Seed::set_columns(Vector<Column>&& columns) {
    if (columns.empty() || columns == no_labels_) {
        label_columns = 0;
        return;
    }

    assert(label_encoder);
    label_columns = label_encoder->cache_column_set(std::move(columns));
}

void Alignment::set_columns(Vector<Column>&& columns) {
    if (columns.empty() || columns == Seed::no_labels_) {
        label_columns = 0;
        return;
    }

    assert(label_encoder);
    label_columns = label_encoder->cache_column_set(std::move(columns));
}

auto Seed::get_columns() const -> const Vector<Column>& {
    if (!label_encoder)
        return no_labels_;

    return label_encoder->get_cached_column_set(label_columns);
}

auto Alignment::get_columns(size_t path_i) const -> const Vector<Column>& {
    if (!label_encoder)
        return Seed::no_labels_;

    assert(path_i < nodes_.size());
    assert(label_column_diffs.empty() || label_column_diffs.size() == nodes_.size() - 1);
    return label_encoder->get_cached_column_set(!path_i || label_column_diffs.empty()
        ? label_columns
        : label_column_diffs[path_i - 1]
    );
}

auto Alignment::get_column_union() const -> Vector<Column> {
    if (!label_encoder)
        return Seed::no_labels_;

    assert(label_column_diffs.empty() || label_column_diffs.size() == nodes_.size() - 1);
    Vector<Column> ret_val = label_encoder->get_cached_column_set(label_columns);
    for (size_t diff : label_column_diffs) {
        if (!diff)
            continue;

        Vector<Column> merge;
        const Vector<Column> &next = label_encoder->get_cached_column_set(diff);
        merge.reserve(ret_val.size() + next.size());
        std::set_union(ret_val.begin(), ret_val.end(), next.begin(), next.end(),
                       std::back_inserter(merge));
        std::swap(merge, ret_val);
    }
    return ret_val;
}

std::vector<std::string> Alignment::get_decoded_labels(size_t path_i) const {
    if (!label_encoder)
        return { "" };

    const auto &columns = get_columns(path_i);
    const auto &encoder = label_encoder->get_annotator().get_label_encoder();
    std::vector<std::string> result;
    result.reserve(columns.size());
    for (Column c : columns) {
        result.push_back(encoder.decode(c));
    }

    return result;
}

void Alignment::merge_annotations(const Alignment &other) {
    if (this == &other)
        return;

    assert(*this == other);
    assert(label_encoder);
    if (label_coordinates.size()) {
        assert(other.label_coordinates.size());
        assert(label_column_diffs.empty() && "label changes not supported");
        assert(extra_scores.empty());
        const auto &a_col = get_columns();
        const auto &b_col = other.get_columns();
        Vector<Column> col_union;
        CoordinateSet coord_union;
        auto add_col_coords = [&](Column c, const auto &coords) {
            col_union.push_back(c);
            coord_union.push_back(coords);
        };
        utils::match_indexed_values(
            a_col.begin(), a_col.end(), label_coordinates.begin(),
            b_col.begin(), b_col.end(), other.label_coordinates.begin(),
            [&](Column c, const auto &coords, const auto &other_coords) {
                col_union.push_back(c);
                Tuple merged_coords;
                std::set_union(coords.begin(), coords.end(),
                               other_coords.begin(), other_coords.end(),
                               std::back_inserter(merged_coords));
                coord_union.emplace_back(std::move(merged_coords));
            },
            add_col_coords, add_col_coords
        );
        std::swap(label_coordinates, coord_union);
        set_columns(std::move(col_union));
        return;
    }

    extra_scores.resize(std::max(extra_scores.size(), other.extra_scores.size()));

    if (other.label_column_diffs.size() && label_column_diffs.empty())
        label_column_diffs.resize(nodes_.size() - 1, label_columns);

    for (size_t i = 0; i < nodes_.size(); ++i) {
        if (!i || label_column_diffs.size()) {
            const auto &a_col = get_columns(i);
            const auto &b_col = other.get_columns(i);
            Vector<Column> col_union;
            std::set_union(a_col.begin(), a_col.end(), b_col.begin(), b_col.end(),
                           std::back_inserter(col_union));
            if (!i) {
                set_columns(std::move(col_union));
            } else {
                label_column_diffs[i - 1] = label_encoder->cache_column_set(std::move(col_union));
            }
        }
        if (i && i - 1 < extra_scores.size() && i - 1 < other.extra_scores.size())
            extra_scores[i - 1] += other.extra_scores[i - 1];
    }
    score_ += other.extra_score;
    extra_score += other.extra_score;
    assert(extra_scores.empty() || extra_scores.size() == nodes_.size() - 1);
}

bool Alignment::splice(Alignment&& other, score_t label_change_score) {
    if (empty()) {
        std::swap(*this, other);
        return has_annotation();
    }

    trim_end_clipping();
    other.trim_clipping();
    return append(std::move(other), label_change_score);
}

bool Alignment::append(Alignment&& other, score_t label_change_score) {
    assert(query_view_.data() + query_view_.size() + other.get_clipping()
            == other.query_view_.data());
    assert(orientation_ == other.orientation_);
    assert(nodes_.size());
    assert(other.nodes_.size());

    bool ret_val = false;

    if (label_coordinates.size() && other.label_coordinates.empty())
        label_coordinates.clear();

    if (has_annotation() && !other.has_annotation()) {
        label_columns = 0;
        label_column_diffs.clear();
        label_encoder = nullptr;
    }

    if (label_coordinates.size()) {
        assert(label_column_diffs.empty() && other.label_column_diffs.empty()
                && "label change not supported with coordinates");
        const auto &columns = get_columns(0);
        const auto &other_columns = other.get_columns(0);
        assert(columns.size() == label_coordinates.size());
        Vector<Column> merged_label_columns;
        CoordinateSet merged_label_coordinates;

        // if the alignments fit together without gaps, make sure that the
        // coordinates form a contiguous range
        utils::match_indexed_values(
            columns.begin(), columns.end(), label_coordinates.begin(),
            other_columns.begin(), other_columns.end(), other.label_coordinates.begin(),
            [&](auto col, const auto &coords, const auto &other_coords) {
                Tuple merged;
                utils::set_intersection(coords.begin(), coords.end(),
                                        other_coords.begin(), other_coords.end(),
                                        std::back_inserter(merged),
                                        sequence_.size());
                if (merged.size()) {
                    merged_label_columns.push_back(col);
                    merged_label_coordinates.push_back(std::move(merged));
                }
            }
        );

        if (merged_label_columns.empty()) {
            *this = Alignment();
            return true;
        }

        ret_val = merged_label_columns.size() < columns.size();

        if (!ret_val) {
            for (size_t i = 0; i < columns.size(); ++i) {
                if (merged_label_coordinates[i].size() < label_coordinates[i].size()) {
                    ret_val = true;
                    break;
                }
            }
        }

        label_columns = label_encoder->cache_column_set(std::move(merged_label_columns));
        std::swap(label_coordinates, merged_label_coordinates);

    } else if (has_annotation()) {
        if ((label_column_diffs.size() ? label_column_diffs.back() : label_columns)
                == other.label_columns) {
            label_change_score = 0;
        }

        if (label_change_score == DBGAlignerConfig::ninf) {
            DEBUG_LOG("Splice failed");
            *this = Alignment();
            return true;
        }

        if (other.label_column_diffs.empty()) {
            other.label_column_diffs.resize(other.nodes_.size(), other.label_columns);
        } else {
            other.label_column_diffs.insert(other.label_column_diffs.begin(), other.label_columns);
        }

        if (other.extra_scores.empty()) {
            other.extra_scores.resize(other.nodes_.size());
            other.extra_scores[0] = label_change_score;
        } else {
            assert(other.extra_scores.size() == other.get_nodes().size() - 1);
            other.extra_scores.insert(other.extra_scores.begin(), label_change_score);
        }
        other.extra_score += other.extra_scores[0];
        other.score_ += other.extra_scores[0];
    }

    if (other.extra_scores.size() && extra_scores.empty()) {
        assert(nodes_.size());
        extra_scores.resize(nodes_.size() - 1);
    }

    if (other.label_column_diffs.size() && label_column_diffs.empty()) {
        assert(nodes_.size());
        label_column_diffs.resize(nodes_.size() - 1, label_columns);
    }

    nodes_.insert(nodes_.end(), other.nodes_.begin(), other.nodes_.end());
    if (other.extra_scores.size())
        extra_scores.insert(extra_scores.end(), other.extra_scores.begin(), other.extra_scores.end());

    if (other.label_column_diffs.size())
        label_column_diffs.insert(label_column_diffs.end(), other.label_column_diffs.begin(), other.label_column_diffs.end());

    assert(extra_scores.empty() || extra_scores.size() == nodes_.size() - 1);
    assert(label_column_diffs.empty() || label_column_diffs.size() == nodes_.size() - 1);

    sequence_ += std::move(other.sequence_);
    score_ += other.score_;
    extra_score += other.extra_score;

    cigar_.append(std::move(other.cigar_));
    // expand the query window to cover both alignments
    query_view_ = std::string_view(query_view_.data(),
                                   other.query_view_.end() - query_view_.begin());
    return ret_val;
}

size_t Alignment::trim_offset(size_t num_nodes) {
    if (!offset_ || nodes_.size() <= 1)
        return 0;

    assert(extra_scores.empty() || extra_scores.size() == nodes_.size() - 1);
    assert(label_column_diffs.empty() || label_column_diffs.size() == nodes_.size() - 1);

    assert(nodes_.front());

    size_t first_dummy = (std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos)
        - nodes_.begin()) - 1;
    size_t trim = std::min({ num_nodes, offset_, nodes_.size() - 1, first_dummy });

    if (!trim)
        return trim;

    offset_ -= trim;
    nodes_.erase(nodes_.begin(), nodes_.begin() + trim);
    if (extra_scores.size()) {
        score_t removed_extra = std::accumulate(extra_scores.begin(),
                                                extra_scores.begin() + trim,
                                                score_t(0));
        extra_score -= removed_extra;
        score_ -= removed_extra;
        extra_scores.erase(extra_scores.begin(), extra_scores.begin() + trim);
    }

    if (label_column_diffs.size()) {
        std::swap(label_columns, label_column_diffs[trim - 1]);
        label_column_diffs.erase(label_column_diffs.begin(), label_column_diffs.begin() + trim);
    }

    assert(extra_scores.empty() || extra_scores.size() == nodes_.size() - 1);
    assert(label_column_diffs.empty() || label_column_diffs.size() == nodes_.size() - 1);
    assert(nodes_.front());
    return trim;
}

void Alignment::extend_offset(std::vector<node_index>&& path,
                             std::vector<size_t>&& columns,
                             std::vector<score_t>&& scores) {
    offset_ += path.size();
    if (columns.size()) {
        assert(columns.size() == path.size());
        if (label_column_diffs.empty())
            label_column_diffs.resize(nodes_.size() - 1, label_columns);

        std::rotate(columns.begin(), columns.begin() + 1, columns.end());
        std::swap(label_columns, columns.back());
        label_column_diffs.insert(label_column_diffs.begin(), columns.begin(), columns.end());
    } else if (label_column_diffs.size()) {
        label_column_diffs.insert(label_column_diffs.begin(), path.size(), label_columns);
    }

    if (scores.size()) {
        assert(scores.size() == path.size());
        if (extra_scores.empty())
            extra_scores.resize(nodes_.size() - 1);

        score_t added = std::accumulate(scores.begin(), scores.end(), score_t{0});
        extra_score += added;
        score_ += added;
        extra_scores.insert(extra_scores.begin(), scores.begin(), scores.end());
    } else if (extra_scores.size()) {
        extra_scores.insert(extra_scores.begin(), path.size(), 0);
    }

    nodes_.insert(nodes_.begin(), path.begin(), path.end());
    assert(extra_scores.empty() || extra_scores.size() == nodes_.size() - 1);
}

size_t Alignment::trim_query_prefix(size_t n,
                                    size_t node_overlap,
                                    const DBGAlignerConfig &config,
                                    bool trim_excess_deletions) {
    size_t clipping = get_clipping();
    const char *query_begin = query_view_.data() - clipping;
    auto it = cigar_.data().begin() + static_cast<bool>(clipping);
    size_t cigar_offset = 0;

    auto s_it = sequence_.begin();
    auto node_it = nodes_.begin();

    while (n || (trim_excess_deletions && it->first == Cigar::DELETION)) {
        if (it == cigar_.data().end()) {
            *this = Alignment();
            return 0;
        }

        switch (it->first) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                assert(s_it != sequence_.end());
                score_ -= config.score_matrix[query_view_[0]][*s_it];
                query_view_.remove_prefix(1);
                --n;
                assert(s_it != sequence_.end());
                ++s_it;
                if (offset_ < node_overlap) {
                    ++offset_;
                } else if (node_it + 1 < nodes_.end()) {
                    ++node_it;
                } else {
                    *this = Alignment();
                    return 0;
                }
            } break;
            case Cigar::INSERTION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                query_view_.remove_prefix(1);
                --n;
            } break;
            case Cigar::DELETION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                assert(s_it != sequence_.end());
                ++s_it;
                if (offset_ < node_overlap) {
                    ++offset_;
                } else if (node_it + 1 < nodes_.end()) {
                    ++node_it;
                } else {
                    *this = Alignment();
                    return 0;
                }
            } break;
            case Cigar::CLIPPED:
            case Cigar::NODE_INSERTION: {
                assert(false && "trimming chains not supported");
            } break;
        }

        ++cigar_offset;
        if (cigar_offset == it->second) {
            ++it;
            cigar_offset = 0;
        }
    }

    size_t seq_trim = s_it - sequence_.begin();
    for (auto &coords : label_coordinates) {
        for (auto &c : coords) {
            c += seq_trim;
        }
    }

    if (!clipping && (cigar_offset || it != cigar_.data().begin()))
        score_ -= config.left_end_bonus;

    nodes_.erase(nodes_.begin(), node_it);
    if (extra_scores.size() && node_it != nodes_.begin()) {
        score_t removed = std::accumulate(extra_scores.begin(),
                                          extra_scores.begin() + (node_it - nodes_.begin()),
                                          score_t(0));
        extra_score -= removed;
        score_ -= removed;
        extra_scores.erase(extra_scores.begin(), extra_scores.begin() + (node_it - nodes_.begin()));
    }

    if (label_column_diffs.size() && node_it != nodes_.begin()) {
        label_columns = label_column_diffs[node_it - nodes_.begin() - 1];
        label_column_diffs.erase(label_column_diffs.begin(), label_column_diffs.begin() + (node_it - nodes_.begin()));
    }

    assert(extra_scores.empty() || extra_scores.size() == nodes_.size() - 1);
    sequence_.erase(sequence_.begin(), s_it);
    it->second -= cigar_offset;
    cigar_.data().erase(cigar_.data().begin(), it);
    extend_query_begin(query_begin);

    return cigar_offset;
}

size_t Alignment::trim_query_suffix(size_t n,
                                    const DBGAlignerConfig &config,
                                    bool trim_excess_deletions) {
    size_t end_clipping = get_end_clipping();
    const char *query_end = query_view_.data() + query_view_.size() + end_clipping;

    trim_end_clipping();
    auto it = cigar_.data().rbegin();
    size_t cigar_offset = 0;

    auto s_it = sequence_.rbegin();
    auto node_it = nodes_.rbegin();

    auto consume_ref = [&]() {
        assert(s_it != sequence_.rend());
        ++s_it;
        if (node_it + 1 < nodes_.rend()) {
            ++node_it;
        } else {
            *this = Alignment();
        }
    };

    while (n || (trim_excess_deletions && it->first == Cigar::DELETION)) {
        if (it == cigar_.data().rend()) {
            *this = Alignment();
            return 0;
        }

        switch (it->first) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                assert(s_it != sequence_.rend());
                score_ -= config.score_matrix[query_view_.back()][*s_it];
                query_view_.remove_suffix(1);
                --n;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::INSERTION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                query_view_.remove_suffix(1);
                --n;
            } break;
            case Cigar::DELETION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::CLIPPED:
            case Cigar::NODE_INSERTION: {
                assert(false && "trimming chains not supported");
            } break;
        }

        ++cigar_offset;
        if (cigar_offset == it->second) {
            ++it;
            cigar_offset = 0;
        }
    }

    if (!end_clipping && (cigar_offset || it.base() != cigar_.data().end()))
        score_ -= config.right_end_bonus;

    nodes_.erase(node_it.base(), nodes_.end());
    if (extra_scores.size() >= nodes_.size()) {
        score_t removed = std::accumulate(extra_scores.begin() + nodes_.size() - 1,
                                          extra_scores.end(),
                                          score_t(0));
        extra_score -= removed;
        score_ -= removed;
        extra_scores.resize(nodes_.size() - 1);
    }

    if (label_column_diffs.size() >= nodes_.size())
        label_column_diffs.resize(nodes_.size() - 1);

    assert(extra_scores.empty() || extra_scores.size() == nodes_.size() - 1);
    sequence_.erase(s_it.base(), sequence_.end());
    it->second -= cigar_offset;
    cigar_.data().erase(it.base(), cigar_.data().end());

    extend_query_end(query_end);

    return cigar_offset;
}

size_t Alignment::trim_reference_prefix(size_t n,
                                        size_t node_overlap,
                                        const DBGAlignerConfig &config,
                                        bool trim_excess_insertions) {
    size_t clipping = get_clipping();
    const char *query_begin = query_view_.data() - clipping;

    auto it = cigar_.data().begin() + static_cast<bool>(clipping);
    size_t cigar_offset = 0;

    auto s_it = sequence_.begin();
    auto node_it = nodes_.begin();

    size_t seq_trim = 0;

    auto consume_ref = [&]() {
        assert(s_it != sequence_.end());
        assert(n);
        if (*s_it != '$') {
            ++seq_trim;
            --n;
        }

        ++s_it;
        if (offset_ < node_overlap) {
            ++offset_;
        } else if (node_it + 1 < nodes_.end()) {
            ++node_it;
        } else {
            *this = Alignment();
        }
    };

    while (n || (trim_excess_insertions && it->first == Cigar::INSERTION)) {
        if (it == cigar_.data().end()) {
            *this = Alignment();
            return 0;
        }

        switch (it->first) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                assert(s_it != sequence_.end());
                score_ -= config.score_matrix[query_view_[0]][*s_it];
                query_view_.remove_prefix(1);
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::INSERTION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                query_view_.remove_prefix(1);
            } break;
            case Cigar::DELETION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::NODE_INSERTION: {} break;
            case Cigar::CLIPPED: {
                assert(false && "this should not happen");
            } break;
        }

        ++cigar_offset;
        if (cigar_offset == it->second) {
            ++it;
            cigar_offset = 0;
        }
    }

    for (auto &coords : label_coordinates) {
        for (auto &c : coords) {
            c += seq_trim;
        }
    }

    if (!clipping && (cigar_offset || it != cigar_.data().begin()))
        score_ -= config.left_end_bonus;

    nodes_.erase(nodes_.begin(), node_it);
    if (extra_scores.size() && node_it != nodes_.begin()) {
        score_t removed = std::accumulate(extra_scores.begin(),
                                          extra_scores.begin() + (node_it - nodes_.begin()),
                                          score_t(0));
        extra_score -= removed;
        score_ -= removed;
        extra_scores.erase(extra_scores.begin(), extra_scores.begin() + (node_it - nodes_.begin()));
    }

    if (label_column_diffs.size() && node_it != nodes_.begin()) {
        label_columns = label_column_diffs[node_it - nodes_.begin() - 1];
        label_column_diffs.erase(label_column_diffs.begin(), label_column_diffs.begin() + (node_it - nodes_.begin()));
    }

    assert(extra_scores.empty() || extra_scores.size() == nodes_.size() - 1);
    sequence_.erase(sequence_.begin(), s_it);
    it->second -= cigar_offset;
    cigar_.data().erase(cigar_.data().begin(), it);

    extend_query_begin(query_begin);

    return cigar_offset;
}

size_t Alignment::trim_reference_suffix(size_t n,
                                        const DBGAlignerConfig &config,
                                        bool trim_excess_insertions) {
    size_t end_clipping = get_end_clipping();
    const char *query_end = query_view_.data() + query_view_.size() + end_clipping;

    trim_end_clipping();
    auto it = cigar_.data().rbegin();
    size_t cigar_offset = 0;

    auto s_it = sequence_.rbegin();
    auto node_it = nodes_.rbegin();

    auto consume_ref = [&]() {
        --n;
        assert(s_it != sequence_.rend());
        ++s_it;
        if (node_it + 1 < nodes_.rend()) {
            ++node_it;
        } else {
            *this = Alignment();
        }
    };

    while (n || (trim_excess_insertions && it->first == Cigar::INSERTION)) {
        if (it == cigar_.data().rend()) {
            *this = Alignment();
            return 0;
        }

        switch (it->first) {
            case Cigar::MATCH:
            case Cigar::MISMATCH: {
                assert(s_it != sequence_.rend());
                score_ -= config.score_matrix[query_view_.back()][*s_it];
                query_view_.remove_suffix(1);
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::INSERTION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                query_view_.remove_suffix(1);
            } break;
            case Cigar::DELETION: {
                score_ -= it->second - cigar_offset == 1
                    ? config.gap_opening_penalty
                    : config.gap_extension_penalty;
                consume_ref();
                if (empty())
                    return 0;
            } break;
            case Cigar::CLIPPED:
            case Cigar::NODE_INSERTION: {
                assert(false && "trimming chains not supported");
            } break;
        }

        ++cigar_offset;
        if (cigar_offset == it->second) {
            ++it;
            cigar_offset = 0;
        }
    }

    if (!end_clipping && (cigar_offset || it.base() != cigar_.data().end()))
        score_ -= config.right_end_bonus;

    nodes_.erase(node_it.base(), nodes_.end());
    if (extra_scores.size() >= nodes_.size()) {
        score_t removed = std::accumulate(extra_scores.begin() + nodes_.size() - 1,
                                          extra_scores.end(),
                                          score_t(0));
        extra_score -= removed;
        score_ -= removed;
        extra_scores.resize(nodes_.size() - 1);
    }

    if (label_column_diffs.size() >= nodes_.size())
        label_column_diffs.resize(nodes_.size() - 1);

    assert(extra_scores.empty() || extra_scores.size() == nodes_.size() - 1);
    sequence_.erase(s_it.base(), sequence_.end());
    it->second -= cigar_offset;
    cigar_.data().erase(it.base(), cigar_.data().end());

    extend_query_end(query_end);

    return cigar_offset;
}

void Alignment::reverse_complement(const DeBruijnGraph &graph,
                                   std::string_view query_rev_comp) {
    assert(query_view_.size() + get_end_clipping() == query_rev_comp.size() - get_clipping());

    trim_offset();
    assert(!offset_ || nodes_.size() == 1);

    if (label_column_diffs.size()) {
        // TODO: make more efficient
        std::reverse(label_column_diffs.begin(), label_column_diffs.end());
        label_column_diffs.push_back(label_columns);
        label_columns = label_column_diffs[0];
        label_column_diffs.erase(label_column_diffs.begin());
    }

    if (extra_scores.size())
        std::reverse(extra_scores.begin(), extra_scores.end());

    if (dynamic_cast<const RCDBG*>(&graph)) {
        if (offset_) {
            *this = Alignment();
        } else {
            std::reverse(cigar_.data().begin(), cigar_.data().end());
            std::reverse(nodes_.begin(), nodes_.end());
            ::reverse_complement(sequence_.begin(), sequence_.end());
            assert(query_rev_comp.size() >= get_clipping() + get_end_clipping());

            orientation_ = !orientation_;
            query_view_ = { query_rev_comp.data() + get_clipping(),
                            query_rev_comp.size() - get_clipping() - get_end_clipping() };
        }
        return;
    }

    if (!offset_) {
        reverse_complement_seq_path(graph, sequence_, nodes_);
    } else {
        assert(nodes_.size() == 1);
        // extract target sequence prefix
        sequence_ = graph.get_node_sequence(nodes_[0]).substr(0, offset_) + sequence_;

        if (sequence_[0] == boss::BOSS::kSentinel) {
            // If the alignment starts from a source k-mer, then traverse forwards
            // until a non-dummy k-mer is hit and check if its reverse complement exists.

            // An appropriate node may not exist if offset_ is greater than the
            // number of sentinel characters (i.e., if some non-sentinel characters
            // from the node prefix are not included).

            // TODO: this cascade of graph unwrapping is ugly, find a cleaner way to do it
            const DeBruijnGraph *base_graph = &graph;
            if (const auto *rc_dbg = dynamic_cast<const RCDBG*>(base_graph))
                base_graph = &rc_dbg->get_graph();

            const auto *canonical = dynamic_cast<const CanonicalDBG*>(base_graph);
            if (canonical)
                base_graph = &canonical->get_graph();

            const auto &dbg_succ = static_cast<const DBGSuccinct&>(*base_graph);

            size_t num_sentinels = sequence_.find_last_of(boss::BOSS::kSentinel) + 1;
            assert(offset_ >= num_sentinels);

            if (canonical && nodes_[0] != canonical->get_base_node(nodes_[0])) {
                // reverse complement of a sink dummy k-mer, no point in traversing
                assert(num_sentinels == 1);
                *this = Alignment();
                return;
            }

            size_t num_first_steps = canonical ? std::min(offset_, num_sentinels) : offset_;

            // the node is present in the underlying graph, so use
            // lower-level methods
            const auto &boss = dbg_succ.get_boss();
            boss::BOSS::edge_index edge = dbg_succ.kmer_to_boss_index(nodes_[0]);
            boss::BOSS::TAlphabet edge_label = boss.get_W(edge) % boss.alph_size;

            // TODO: This picks the node which is found by always traversing
            // the last outgoing edge. Is there a better way to pick a node?
            for (size_t i = 0; i < num_first_steps; ++i) {
                edge = boss.fwd(edge, edge_label);
                edge_label = boss.get_W(edge) % boss.alph_size;
                if (edge_label == boss::BOSS::kSentinelCode) {
                    // reverse complement not found
                    assert(offset_ > num_sentinels);
                    *this = Alignment();
                    return;
                }

                nodes_[0] = dbg_succ.boss_to_kmer_index(edge);
                assert(nodes_[0]);
                sequence_.push_back(boss.decode(edge_label));
                assert(graph.get_node_sequence(nodes_[0])
                    == sequence_.substr(sequence_.size() - graph.get_k()));
            }

            for (size_t i = num_first_steps; i < offset_; ++i) {
                node_index next_node = 0;
                char last_char;
                canonical->call_outgoing_kmers(nodes_[0], [&](node_index next, char c) {
                    if (c == boss::BOSS::kSentinel)
                        return;

                    next_node = next;
                    last_char = c;
                });

                if (!next_node) {
                    *this = Alignment();
                    return;
                } else {
                    nodes_[0] = next_node;
                    sequence_.push_back(last_char);
                    assert(graph.get_node_sequence(nodes_[0])
                        == sequence_.substr(sequence_.size() - graph.get_k()));
                }
            }

            assert(sequence_.size() == dbg_succ.get_k() + offset_);
            sequence_ = sequence_.substr(offset_);

            assert(nodes_ == map_to_nodes_sequentially(graph, sequence_));
            reverse_complement_seq_path(graph, sequence_, nodes_);

            assert(std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos)
                    == nodes_.end());

            sequence_.assign(sequence_.data() + offset_, graph.get_k() - offset_);

        } else {
            assert(nodes_.size() == 1);
            assert(nodes_ == map_to_nodes_sequentially(graph, sequence_));
            reverse_complement_seq_path(graph, sequence_, nodes_);

            assert(std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos)
                    == nodes_.end());

            // trim off ending from reverse complement (corresponding to the added prefix)
            for (size_t i = 0; i < offset_; ++i) {
                size_t indegree = 0;
                graph.adjacent_incoming_nodes(nodes_[0], [&](node_index prev) {
                    ++indegree;

                    // TODO: there are multiple possible reverse complements, which
                    // do we pick? Currently we pick the first one
                    if (indegree == 1)
                        nodes_[0] = prev;
                });

                if (!indegree) {
                    *this = Alignment();
                    return;
                }

                sequence_.pop_back();
                assert(graph.get_node_sequence(nodes_[0]).substr(i + 1)
                    == sequence_.substr(sequence_.size() - graph.get_k() + i + 1));
            }

            assert(sequence_.size() == graph.get_k() - offset_);
        }

        assert(graph.get_node_sequence(nodes_[0]).substr(offset_) == sequence_);
    }

    std::reverse(cigar_.data().begin(), cigar_.data().end());
    assert(query_rev_comp.size() >= get_clipping() + get_end_clipping());

    orientation_ = !orientation_;
    query_view_ = { query_rev_comp.data() + get_clipping(),
                    query_rev_comp.size() - get_clipping() - get_end_clipping() };
    assert(is_valid(graph));
}

// derived from:
// https://github.com/maickrau/GraphAligner/blob/236e1cf0514cfa9104e9a3333cdc1c43209c3c5a/src/vg.proto
Json::Value path_json(const std::vector<DeBruijnGraph::node_index> &nodes,
                      const Cigar &cigar,
                      size_t node_size,
                      std::string_view query_view,
                      size_t offset,
                      const std::string &label) {
    assert(nodes.size());

    Json::Value path;

    auto cigar_it = cigar.data().begin();
    if (cigar.size() && cigar_it->first == Cigar::CLIPPED) {
        cigar_it++;
    }

    size_t cigar_offset = 0;
    assert(cigar_it != cigar.data().end());

    int64_t rank = 1;
    const char *query_start = query_view.data();

#ifndef NDEBUG
    const char *query_end = query_view.data() + query_view.size();
#endif

    size_t cur_pos = rank == 1 ? offset : 0;

    Json::Value mapping;
    Json::Value position;
    position["node_id"] = nodes.front();

    if (cur_pos)
        position["offset"] = Json::Value::UInt64(cur_pos);

    // set to true if the node is the reverse complement of the query
    //position["is_reverse"] = false;

    mapping["position"] = position;

    // handle alignment to the first node
    while (cur_pos < node_size && cigar_it != cigar.data().end()) {
        assert(cigar_it->second > cigar_offset);
        size_t next_pos = std::min(node_size,
                                   cur_pos + (cigar_it->second - cigar_offset));
        size_t next_size = next_pos - cur_pos;
        assert(cigar_offset + next_size <= cigar_it->second);

        Json::Value edit;
        switch (cigar_it->first) {
            case Cigar::MISMATCH: {
                assert(query_start + next_size <= query_end);
                edit["from_length"] = Json::Value::UInt64(next_size);
                edit["to_length"] = Json::Value::UInt64(next_size);
                edit["sequence"] = std::string(query_start, next_size);
                query_start += next_size;
            } break;
            case Cigar::INSERTION: {
                assert(query_start + next_size <= query_end);
                // this assumes that INSERTIONs can't happen right after DELETIONs
                //edit["from_length"] = 0;
                edit["to_length"] = Json::Value::UInt64(next_size);
                edit["sequence"] = std::string(query_start, next_size);
                query_start += next_size;

                // the target is not consumed, so reset the position
                next_pos = cur_pos;
            } break;
            case Cigar::DELETION: {
                edit["from_length"] = Json::Value::UInt64(next_size);
                //edit["to_length"] = 0;
            } break;
            case Cigar::MATCH: {
                assert(query_start + next_size <= query_end);
                edit["from_length"] = Json::Value::UInt64(next_size);
                edit["to_length"] = Json::Value::UInt64(next_size);
                query_start += next_size;
            } break;
            case Cigar::CLIPPED: {
                ++cigar_it;
                cigar_offset = 0;
                assert(cigar_it == cigar.data().end());
                continue;
            } break;
            case Cigar::NODE_INSERTION: {
                assert(false && "this should not be reached");
            } break;
        }

        cigar_offset += next_size;
        cur_pos = next_pos;

        mapping["edit"].append(edit);

        if (cigar_offset == cigar_it->second) {
            ++cigar_it;
            cigar_offset = 0;
        }
    }

    mapping["rank"] = rank++;
    path["mapping"].append(mapping);

    // handle the rest of the alignment
    for (auto node_it = nodes.begin() + 1; node_it != nodes.end(); ++node_it) {
        assert(cigar_it != cigar.data().end());
        assert(cigar_it->second > cigar_offset);

        Json::Value mapping;
        Json::Value position;
        position["node_id"] = *node_it;
        position["offset"] = Json::Value::UInt64(node_size - 1);
        // set to true if the node is the reverse complement of the query
        //position["is_reverse"] = false;
        mapping["position"] = position;

        if (cigar_it->first == Cigar::INSERTION || cigar_it->first == Cigar::CLIPPED) {
            Json::Value edit;
            size_t length = cigar_it->second - cigar_offset;
            assert(query_start + length < query_end);
            // TODO: this assumes that INSERTIONs can't happen right after DELETIONs
            //edit["from_length"] = 0;
            edit["to_length"] = Json::Value::UInt64(length);
            edit["sequence"] = std::string(query_start, length);
            query_start += length;
            ++cigar_it;
            cigar_offset = 0;
            mapping["edit"].append(edit);
            assert(cigar_it != cigar.data().end());
        }

        Json::Value edit;
        switch (cigar_it->first) {
            case Cigar::MISMATCH: {
                assert(query_start < query_end);
                edit["from_length"] = 1;
                edit["to_length"] = 1;
                edit["sequence"] = std::string(query_start, 1);
                query_start++;
            } break;
            case Cigar::DELETION: {
                edit["from_length"] = 1;
                //edit["to_length"] = 0;
            } break;
            case Cigar::MATCH: {
                edit["from_length"] = 1;
                edit["to_length"] = 1;
                query_start++;
            } break;
            case Cigar::INSERTION:
            case Cigar::NODE_INSERTION:
            case Cigar::CLIPPED: assert(false && "this should not be reached"); break;
        }

        if (++cigar_offset == cigar_it->second) {
            cigar_offset = 0;
            ++cigar_it;
        }

        mapping["edit"].append(edit);
        mapping["rank"] = rank++;
        path["mapping"].append(mapping);
    }

    assert(query_start == query_end);
    assert(cigar_it == cigar.data().end()
            || (cigar_it + 1 == cigar.data().end() && cigar_it->first == Cigar::CLIPPED));

    path["length"] = Json::Value::UInt64(nodes.size());

    if (nodes.front() == nodes.back())
        path["is_circular"] = true;

    path["name"] = label;

    return path;
}

Json::Value Alignment::to_json(size_t node_size,
                               bool is_secondary,
                               const std::string &read_name,
                               const std::string &label) const {
    if (extra_score)
        throw std::runtime_error("Alignments from PSSMs not supported");

    if (sequence_.find("$") != std::string::npos
            || std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos) != nodes_.end()) {
        throw std::runtime_error("JSON output for chains not supported");
    }

    std::string_view full_query = get_full_query_view();

    // encode alignment
    Json::Value alignment;

    alignment["name"] = read_name;
    alignment["sequence"] = std::string(full_query);

    if (sequence_.size())
        alignment["annotation"]["ref_sequence"] = sequence_;

    if (query_view_.empty())
        return alignment;

    assert(query_view_.data() - cigar_.get_clipping() == full_query.data());
    assert(query_view_.size() + cigar_.get_clipping() + cigar_.get_end_clipping()
            == full_query.size());

    alignment["annotation"]["cigar"] = cigar_.to_string();

    // encode path
    if (nodes_.size())
        alignment["path"] = path_json(nodes_, cigar_, node_size, query_view_, offset_, label);

    alignment["score"] = static_cast<int32_t>(score_);

    if (cigar_.get_clipping()) {
        alignment["query_position"] = static_cast<int32_t>(cigar_.get_clipping());
        alignment["soft_clipped"] = static_cast<bool>(cigar_.get_clipping());
    }

    if (is_secondary)
        alignment["is_secondary"] = is_secondary;

    alignment["identity"] = query_view_.size()
        ? static_cast<double>(cigar_.get_num_matches()) / query_view_.size()
        : 0;

    alignment["read_mapped"] = static_cast<bool>(query_view_.size());

    if (orientation_)
        alignment["read_on_reverse_strand"] = orientation_;

    // Unused flags (for now)
    //alignment["quality"]; // bytes
    //alignment["mapping_quality"]; // int32
    //alignment["read_group"]; // string
    //alignment["fragment_prev"]; // Alignment
    //alignment["fragment_next"]; // Alignment
    //alignment["fragment"]; // Path
    //alignment["refpos"]; // Position
    //alignment["read_paired"]; // bool
    //alignment["mate_unmapped"]; // bool
    //alignment["mate_on_reverse_strand"]; // bool
    //alignment["discordant_insert_size"]; // bool
    //alignment["uniqueness"]; // double
    //alignment["correct"]; // double
    //alignment["secondary_score"]; // repeated int32
    //alignment["fragment_score"]; // double
    //alignment["mate_mapped_to_disjoint_subgraph"]; // bool
    //alignment["fragment_length_distribution"]; // string
    //alignment["haplotype_scored"]; // bool
    //alignment["haplotype_logprob"]; // double
    //alignment["time_used"]; // double
    //alignment["to_correct"]; // Position
    //alignment["correctly_mapped"]; // bool

    return alignment;
}

void Alignment::load_from_json(const Json::Value &alignment,
                               const DeBruijnGraph &graph,
                               std::string *query_sequence) {
    assert(query_sequence);

    cigar_ = Cigar();
    nodes_.clear();
    sequence_.clear();

    *query_sequence = alignment["sequence"].asString();
    orientation_ = alignment["read_on_reverse_strand"].asBool();
    score_ = alignment["score"].asInt();
    const Json::Value &mapping = alignment["path"]["mapping"];
    assert(mapping.size() == alignment["path"]["length"].asUInt64());

    cigar_.append(Cigar::CLIPPED, alignment["query_position"].asInt());

    const char *this_query_begin = query_sequence->c_str() + get_clipping();
    size_t this_query_size = 0;

    size_t path_steps = 0;
    Json::ArrayIndex i = 0;
    offset_ = mapping[i]["position"]["offset"].asUInt64();

    for ( ; i < mapping.size(); ++i) {
        nodes_.emplace_back(mapping[i]["position"]["node_id"].asUInt64());
        if (nodes_.size() == 1) {
            sequence_ = graph.get_node_sequence(nodes_.back()).substr(offset_);
        } else {
            graph.call_outgoing_kmers(*(nodes_.rbegin() + 1),
                                      [&](auto node, char c) {
                if (node == nodes_.back())
                    sequence_.push_back(c);
            });
        }
        const Json::Value &edits = mapping[i]["edit"];

        for (Json::ArrayIndex j = 0; j < edits.size(); ++j) {
            assert(edits[j]["from_length"].asUInt64()
                       || edits[j]["to_length"].asUInt64());

            if (edits[j]["from_length"] == edits[j]["to_length"]) {
                if (edits[j]["sequence"].asString().size()) {
                    cigar_.append(Cigar::MISMATCH, edits[j]["from_length"].asUInt64());
                } else {
                    cigar_.append(Cigar::MATCH, edits[j]["from_length"].asUInt64());
                }

                path_steps += edits[j]["from_length"].asUInt64();
                this_query_size += edits[j]["to_length"].asUInt64();

            } else if (edits[j]["from_length"].asUInt64()) {
                cigar_.append(Cigar::DELETION, edits[j]["from_length"].asUInt64());
                path_steps += edits[j]["from_length"].asUInt64();

            } else {
                assert(edits[j]["to_length"].asUInt64());
                cigar_.append(Cigar::INSERTION, edits[j]["to_length"].asUInt64());
                this_query_size += edits[j]["to_length"].asUInt64();
            }
        }
    }

    query_view_ = { this_query_begin, this_query_size };

    if (query_view_.size() + get_clipping() < query_sequence->size()) {
        cigar_.append(Cigar::CLIPPED,
                      query_sequence->size() - query_view_.size() - get_clipping());
    }

    if (alignment["annotation"]["cigar"]
            && cigar_ != Cigar(alignment["annotation"]["cigar"].asString())) {
        throw std::runtime_error("ERROR: CIGAR and mapping mismatch");
    }

    sequence_ = sequence_.substr(0, path_steps);
    assert(!alignment["annotation"]["ref_sequence"]
        || alignment["annotation"]["ref_sequence"].asString() == sequence_);

    if (!is_valid(graph))
        throw std::runtime_error("ERROR: JSON reconstructs invalid alignment");
}

void Alignment::splice_with_unknown(Alignment&& other,
                                    size_t num_unknown,
                                    size_t node_overlap,
                                    const DBGAlignerConfig &config) {
    assert(!empty());
    assert(num_unknown);
    assert(other.get_clipping() >= num_unknown);

    if (other.get_offset()) {
        logger->trace("Can't splice in sub-k alignment");
        *this = Alignment();
        return;
    }

    ssize_t overlap = static_cast<ssize_t>(get_clipping() + query_view_.size())
                            - other.get_clipping();

    other.trim_clipping();
    if (overlap <= 0) {
        // there is a normal gap between the alignments due to missing characters
        // in the graph (typically N)
        trim_end_clipping();
        size_t query_gap = -overlap;
        if (query_gap > num_unknown) {
            cigar_.append(Cigar::INSERTION, query_gap - num_unknown);
            score_ += static_cast<score_t>(config.gap_opening_penalty)
                            + static_cast<score_t>(query_gap - num_unknown - 1)
                                * static_cast<score_t>(config.gap_extension_penalty);
            query_view_ = std::string_view{ query_view_.data(),
                                            query_view_.size() + query_gap - num_unknown };
            query_gap = num_unknown;
        }

        const char *start = other.query_view_.data() - query_gap;
        if (query_gap) {
            other.cigar_.data().insert(other.cigar_.data().begin(),
                                       Cigar::value_type{ Cigar::MISMATCH, query_gap });
            other.score_ += config.score_sequences(std::string_view{ start, query_gap },
                                                   std::string(num_unknown, '$'));
        }

        if (num_unknown > query_gap) {
            other.cigar_.data().insert(other.cigar_.data().begin(),
                                       Cigar::value_type{ Cigar::DELETION,
                                                          num_unknown - query_gap });
            other.score_ += static_cast<score_t>(config.gap_opening_penalty)
                            + static_cast<score_t>(num_unknown - query_gap - 1)
                                * static_cast<score_t>(config.gap_extension_penalty);
        }

        other.cigar_.data().insert(other.cigar_.data().begin(),
                                   Cigar::value_type{ Cigar::NODE_INSERTION,
                                                      node_overlap + num_unknown - other.offset_ });
        other.score_ += config.node_insertion_penalty;
        other.query_view_ = std::string_view(start, other.query_view_.size() + query_gap);
        assert(query_view_.data() + query_view_.size() == other.query_view_.data());
    } else {
        // This can happen if there's a gap in the graph (due to N) at a point
        // where read has fewer copies of a repetitive sequence relative to the graph.
        // Correct for this by re-constructing the deletion.
        auto nodes = nodes_;
        auto seq = sequence_;
        trim_query_suffix(overlap, config, false);
        if (empty()) {
            // TODO: splice from both ends
            logger->trace("Not enough nodes in this alignment for splicing");
            *this = Alignment();
            return;
        }
        overlap = seq.size() - sequence_.size();

        trim_end_clipping();
        if (overlap) {
            cigar_.data().emplace_back(Cigar::DELETION, overlap);
            nodes_.insert(nodes_.end(), nodes.end() - overlap, nodes.end());
            if (extra_scores.size())
                extra_scores.resize(nodes_.size() - 1);

            if (label_column_diffs.size())
                label_column_diffs.resize(nodes_.size() - 1);

            sequence_ += std::string_view(seq.data() + seq.size() - overlap, overlap);
        }

        other.cigar_.data().insert(other.cigar_.data().begin(),
                                   Cigar::value_type{ Cigar::DELETION, num_unknown });
        other.score_ += static_cast<score_t>(config.node_insertion_penalty)
                            + static_cast<score_t>(config.gap_opening_penalty)
                                + static_cast<score_t>(num_unknown - 1)
                                    * static_cast<score_t>(config.gap_extension_penalty);
        other.cigar_.data().insert(other.cigar_.data().begin(),
                                   Cigar::value_type{ Cigar::NODE_INSERTION,
                                                      node_overlap + num_unknown });
    }

    other.sequence_ = std::string(num_unknown, '$') + other.sequence_;

    other.nodes_.insert(other.nodes_.begin(),
                        node_overlap + num_unknown - other.offset_,
                        DeBruijnGraph::npos);
    if (other.extra_scores.size()) {
        other.extra_scores.insert(other.extra_scores.begin(),
                                  node_overlap + num_unknown - other.offset_,
                                  0);
    }

    if (other.label_column_diffs.size()) {
        other.label_column_diffs.insert(other.label_column_diffs.begin(),
                                        node_overlap + num_unknown - other.offset_,
                                        0);
    }

    other.offset_ = node_overlap;
    for (auto &tuple : other.label_coordinates) {
        for (auto &c : tuple) {
            c -= num_unknown;
        }
    }

    append(std::move(other));
}

void Alignment::insert_gap_prefix(ssize_t gap_length,
                                  size_t node_overlap,
                                  const DBGAlignerConfig &config) {
    assert(size());
    size_t extra_nodes = node_overlap + 1;

    if (gap_length < 0) {
        // alignments overlap
        // extra_nodes = k - 1 - matching_overlap
        // e.g.,
        // k = 4
        // overlap = 4
        // matching overlap = 2
        // ATGCTATGCA
        //       ACCAACGACT

        trim_clipping();
        assert(extra_nodes + gap_length > 1);
        extra_nodes += gap_length - 1;

        if (offset_) {
            // if there are suffix-mapped nodes, only keep the ones that are
            // part of the overlap
            assert(static_cast<ssize_t>(offset_) >= -gap_length);
            assert(nodes_.size() > offset_ + gap_length);
            nodes_.erase(nodes_.begin(), nodes_.begin() + offset_ + gap_length);
            if (offset_ + gap_length) {
                if (extra_scores.size()) {
                    score_t removed = std::accumulate(extra_scores.begin(),
                                                      extra_scores.begin() + offset_ + gap_length,
                                                      score_t(0));
                    extra_score -= removed;
                    score_ -= removed;
                    extra_scores.erase(extra_scores.begin(), extra_scores.begin() + offset_ + gap_length);
                }

                if (label_column_diffs.size())
                    label_column_diffs.erase(label_column_diffs.begin(), label_column_diffs.begin() + offset_ + gap_length);
            }
        }

        if (extra_nodes) {
            // they can't be joined since the overlap is too small
            // ATGCTATGCA
            //           ACGACT
            //       TGCA
            //        GCAA - added
            //         CAAC
            //          AACG
            //           ACGA
            cigar_.data().insert(cigar_.data().begin(),
                                 Cigar::value_type{ Cigar::NODE_INSERTION, extra_nodes });
            score_ += config.node_insertion_penalty;
        }
    } else {
        // no overlap
        // extra_nodes = k
        // e.g.,
        // k = 4
        // gap = 2
        // ATGCTATGCA
        //             ACGTACGACT
        //       TGCA
        //        GCA$ - added
        //         CA$A - added
        //          A$AC - added
        //           $ACG - added
        //            ACGT

        trim_offset();
        if (offset_) {
            assert(false && "extra node addition to sub-k alignments not implemented");
        }
        assert(get_clipping() >= gap_length);
        trim_clipping();

        sequence_ = std::string(1, '$') + sequence_;
        cigar_.data().insert(cigar_.data().begin(), Cigar::value_type{ Cigar::DELETION, 1 });
        score_ += config.gap_opening_penalty;

        assert(extra_nodes >= 2);
        cigar_.data().insert(cigar_.data().begin(),
                             Cigar::value_type{ Cigar::NODE_INSERTION, extra_nodes - 1 });
        score_ += config.node_insertion_penalty;

        if (gap_length) {
            cigar_.data().insert(cigar_.data().begin(), Cigar::value_type{ Cigar::INSERTION, gap_length });
            score_ += config.gap_opening_penalty
                        + (gap_length - 1) * config.gap_extension_penalty;
            query_view_ = std::string_view(query_view_.data() - gap_length,
                                           query_view_.size() + gap_length);
        }
    }

    nodes_.insert(nodes_.begin(), extra_nodes, DeBruijnGraph::npos);

    if (extra_scores.size() && extra_nodes) {
        extra_scores.insert(extra_scores.begin(), extra_nodes, 0);
        assert(extra_scores.size() == nodes_.size() - 1);
    }

    if (extra_nodes && has_annotation()) {
        if (label_column_diffs.empty()) {
            label_column_diffs.resize(nodes_.size() - 1);
            std::fill(label_column_diffs.begin() + extra_nodes - 1, label_column_diffs.end(), label_columns);
            label_columns = 0;
        } else {
            label_column_diffs.insert(label_column_diffs.begin(), extra_nodes, 0);
            std::swap(label_column_diffs[extra_nodes - 1], label_columns);
        }
    }
    offset_ = node_overlap;

    assert(nodes_.size() == sequence_.size());
}

/**
 * Partition the alignment at the last k-mer. Return a pair containing the
 * alignment of all but the last k-mers, and the alignment of the last k-mer.
 */
std::pair<Alignment, Alignment> Alignment
::split_seed(size_t node_overlap, const DBGAlignerConfig &config) const {
    if (nodes_.size() <= 1
            || std::find(nodes_.begin(), nodes_.end(), DeBruijnGraph::npos) != nodes_.end()) {
        return std::make_pair(Alignment(), *this);
    }

    auto it = cigar_.data().rbegin() + static_cast<bool>(cigar_.data().back().first == Cigar::CLIPPED);
    if (it->first != Cigar::MATCH || it->second < 2)
        return std::make_pair(Alignment(), *this);

    size_t to_trim = std::min(static_cast<size_t>(it->second) - 1, nodes_.size() - 1);
    auto ret_val = std::make_pair(*this, *this);
    ret_val.second.trim_reference_prefix(sequence_.size() - to_trim, node_overlap, config);
    assert(ret_val.second.size());

    ret_val.first.trim_reference_suffix(to_trim, config, false);
    assert(ret_val.first.size());
    return ret_val;
}

// Return the string spelled by the path. This path may have disconnects (if it came)
// from a chain alignment), so this method handles that case. If there is an invalid
// edge, or if there is too long of a stretch of npos nodes, this throws a runtime error.
std::string spell_path(const DeBruijnGraph &graph,
                       const std::vector<DeBruijnGraph::node_index> &path,
                       size_t offset) {
    std::string seq;
    assert(offset < graph.get_k());

    if (path.empty())
        return seq;

    seq.reserve(path.size() + graph.get_k() - 1 - offset);

    size_t num_dummy = 0;
    size_t num_unknown = 0;
    if (path.front()) {
        seq += graph.get_node_sequence(path.front()).substr(offset);
    } else {
        num_unknown = graph.get_k() - offset;
        seq += std::string(num_unknown, '$');
        num_dummy = 1;
    }

    for (size_t i = 1; i < path.size(); ++i) {
        if (path[i]) {
            if (num_dummy) {
                seq += '$';
                ++num_unknown;
                std::string next_seq = graph.get_node_sequence(path[i]);
                std::string_view window(next_seq);
                if (next_seq.size() > seq.size())
                    window.remove_prefix(next_seq.size() - seq.size());

                std::transform(window.rbegin(), window.rend(), seq.rbegin(), [&](char c) {
                    if (c != '$')
                       --num_unknown;

                    return c;
                });
                num_dummy = 0;
            } else {
                char next = '\0';
                graph.call_outgoing_kmers(path[i - 1], [&](auto next_node, char c) {
                    if (next_node == path[i])
                        next = c;
                });

                if (!next) {
                    logger->error("Invalid edge {} {}\t{} {}",
                                  path[i - 1], path[i],
                                  graph.get_node_sequence(path[i - 1]),
                                  graph.get_node_sequence(path[i]));
                    throw std::runtime_error("");
                }

                seq += next;
                if (num_unknown) {
                    std::string next_seq = graph.get_node_sequence(path[i]);
                    auto it = seq.end() - next_seq.size();
                    for (char c : next_seq) {
                        if (*it == '$' && c != '$') {
                            --num_unknown;
                            *it = c;
                        }

                        ++it;
                    }
                }
            }
        } else {
            seq += '$';
            ++num_dummy;
            ++num_unknown;
        }
    }

    assert(seq.size() == path.size() + graph.get_k() - 1 - offset);
    return seq;
}

bool Alignment::is_valid(const DeBruijnGraph &graph, const DBGAlignerConfig *config) const {
    if (empty())
        return true;

    try {
        std::string spelling = spell_path(graph, nodes_, offset_);
        if (spelling != sequence_) {
            logger->error("Stored sequence is incorrect\n{}\t{}", spelling, *this);
            return false;
        }
    } catch (const std::runtime_error&) {
        logger->error("{}", *this);
        return false;
    }

    if (!cigar_.is_valid(sequence_, query_view_)) {
        logger->error("{}", *this);
        return false;
    }

    if (extra_scores.size() && extra_scores.size() != nodes_.size() - 1) {
        logger->error("Extra score array incorrect size: {} vs. {}\n{}",
                      extra_scores.size(), nodes_.size() - 1, *this);
        return false;
    }

    score_t change_score_sum = std::accumulate(extra_scores.begin(), extra_scores.end(),
                                               score_t(0));
    if (extra_score != change_score_sum) {
        logger->error("Mismatch between extra score array and extra score sum: {} {} vs. {}\n{}",
                      fmt::join(extra_scores, ","), change_score_sum, extra_score, *this);
        return false;
    }

    score_t cigar_score = config ? config->score_cigar(sequence_, query_view_, cigar_) : 0;
    cigar_score += extra_score;
    if (config && score_ != cigar_score) {
        logger->error("Mismatch between CIGAR and score\nCigar score: {} ({} from extra)\n{}\t{}",
                      cigar_score, extra_score, query_view_, *this);
        return false;
    }

    if (label_column_diffs.size() && label_column_diffs.size() != nodes_.size() - 1) {
        logger->error("Label storage array incorrect size: {} vs. {}\n{}",
                      label_column_diffs.size(), nodes_.size() - 1, *this);
        return false;
    }

    return true;
}


AlignmentResults::AlignmentResults(std::string_view query) {
    // Pad sequences for easier access in 64-bit blocks.
    // Take the max of the query size and sizeof(query_view_) to ensure that small-string
    // optimizations are disabled. Implementation taken from:
    // https://stackoverflow.com/questions/34788789/disable-stdstrings-sso
    query_.reserve(std::max(query.size(), sizeof(query_)) + 8);

    // TODO: use alphabet encoder
    // transform to upper and fix non-standard characters
    std::transform(query.begin(), query.end(), std::back_inserter(query_),
                   [](char c) { return c >= 0 ? toupper(c) : 127; });

    // fill padding with '\0'
    memset(query_.data() + query.size(), '\0', query_.capacity() - query.size());

    // set the reverse complement
    query_rc_ = query_;

    // fill padding just in case optimizations removed it
    query_rc_.reserve(query_.capacity());
    memset(query_rc_.data() + query.size(), '\0', query_rc_.capacity() - query.size());

    // reverse complement
    reverse_complement(query_rc_.begin(), query_rc_.end());
}

} // namespace align
} // namespace graph
} // namespace mtg
