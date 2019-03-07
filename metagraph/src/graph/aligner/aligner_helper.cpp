#include "aligner_helper.hpp"

#include "utils.hpp"


std::array<std::array<Cigar::Operator, 128>, 128> Cigar::initialize_opt_table() {
    std::array<std::array<Cigar::Operator, 128>, 128> char_to_op;

    for (auto& row : char_to_op) {
        row.fill(Operator::MISMATCH_TRANSVERSION);
    }

    char_to_op['a']['a'] = char_to_op['A']['a'] = char_to_op['a']['A'] = char_to_op['A']['A'] = Operator::MATCH;
    char_to_op['c']['c'] = char_to_op['C']['c'] = char_to_op['c']['C'] = char_to_op['C']['C'] = Operator::MATCH;
    char_to_op['g']['g'] = char_to_op['G']['g'] = char_to_op['g']['G'] = char_to_op['G']['G'] = Operator::MATCH;
    char_to_op['t']['t'] = char_to_op['T']['t'] = char_to_op['t']['T'] = char_to_op['T']['T'] = Operator::MATCH;
    char_to_op['a']['g'] = char_to_op['A']['g'] = char_to_op['a']['G'] = char_to_op['A']['G'] = Operator::MISMATCH_TRANSITION;
    char_to_op['g']['a'] = char_to_op['G']['a'] = char_to_op['g']['A'] = char_to_op['G']['A'] = Operator::MISMATCH_TRANSITION;
    char_to_op['c']['t'] = char_to_op['C']['t'] = char_to_op['c']['T'] = char_to_op['C']['T'] = Operator::MISMATCH_TRANSITION;
    char_to_op['t']['c'] = char_to_op['T']['c'] = char_to_op['t']['C'] = char_to_op['T']['C'] = Operator::MISMATCH_TRANSITION;
    char_to_op['a']['c'] = char_to_op['A']['c'] = char_to_op['a']['C'] = char_to_op['A']['C'] = Operator::MISMATCH_TRANSVERSION;
    char_to_op['a']['t'] = char_to_op['A']['t'] = char_to_op['a']['T'] = char_to_op['A']['T'] = Operator::MISMATCH_TRANSVERSION;
    char_to_op['g']['c'] = char_to_op['G']['c'] = char_to_op['g']['C'] = char_to_op['G']['C'] = Operator::MISMATCH_TRANSVERSION;
    char_to_op['g']['t'] = char_to_op['G']['t'] = char_to_op['g']['T'] = char_to_op['G']['T'] = Operator::MISMATCH_TRANSVERSION;
    char_to_op['c']['a'] = char_to_op['C']['a'] = char_to_op['c']['A'] = char_to_op['C']['A'] = Operator::MISMATCH_TRANSVERSION;
    char_to_op['c']['g'] = char_to_op['C']['g'] = char_to_op['c']['G'] = char_to_op['C']['G'] = Operator::MISMATCH_TRANSVERSION;
    char_to_op['t']['a'] = char_to_op['T']['a'] = char_to_op['t']['A'] = char_to_op['T']['A'] = Operator::MISMATCH_TRANSVERSION;
    char_to_op['t']['g'] = char_to_op['T']['g'] = char_to_op['t']['G'] = char_to_op['T']['G'] = Operator::MISMATCH_TRANSVERSION;

    return char_to_op;
}

const std::array<std::array<Cigar::Operator, 128>, 128> Cigar::char_to_op
    = Cigar::initialize_opt_table();

char opt_to_char(Cigar::Operator op) {
    switch (op) {
        case Cigar::Operator::MATCH: return '=';
        case Cigar::Operator::MISMATCH_TRANSITION: return 'X';
        case Cigar::Operator::MISMATCH_TRANSVERSION: return 'X';
        case Cigar::Operator::INSERTION: return 'I';
        case Cigar::Operator::DELETION: return 'D';
        case Cigar::Operator::CLIPPED: return 'S';
        default:
            assert(false);
            return '\0';
    }
}

std::string Cigar::to_string() const {
    if (cigar_.empty())
        return "";

    // fold down to CIGAR operations
    std::vector<std::pair<char, uint64_t>> converted;
    converted.reserve(cigar_.size());
    converted.emplace_back(opt_to_char(cigar_.front().first), cigar_.front().second);

    for (auto it = std::next(cigar_.begin()); it != cigar_.end(); ++it) {
        auto op_char = opt_to_char(it->first);
        if (converted.back().first != op_char) {
            converted.emplace_back(op_char, it->second);
        } else {
            converted.back().second += it->second;
        }
    }

    // generate string
    std::string cigar_string;
    for (const auto &pair : converted) {
        cigar_string += std::to_string(pair.second) + pair.first;
    }

    return cigar_string;
}

void Cigar::append(Operator op, LengthType num) {
    if (!num)
        return;

    if (cigar_.empty() || cigar_.back().first != op) {
        cigar_.emplace_back(op, num);
    } else {
        cigar_.back().second += num;
    }
}

void Cigar::append(Cigar&& other) {
    if (other.empty())
        return;

    append(other.cigar_.front().first, other.cigar_.front().second);
    cigar_.insert(cigar_.end(), std::next(other.cigar_.begin()), other.cigar_.end());
}


void DBGAlignerConfig::set_match_score(int8_t match_score) {
    match_score_ = match_score;
    score_matrix_['a']['a'] = score_matrix_['A']['a'] = score_matrix_['a']['A'] = score_matrix_['A']['A'] = match_score;
    score_matrix_['c']['c'] = score_matrix_['C']['c'] = score_matrix_['c']['C'] = score_matrix_['C']['C'] = match_score;
    score_matrix_['g']['g'] = score_matrix_['G']['g'] = score_matrix_['g']['G'] = score_matrix_['G']['G'] = match_score;
    score_matrix_['t']['t'] = score_matrix_['T']['t'] = score_matrix_['t']['T'] = score_matrix_['T']['T'] = match_score;
}

void DBGAlignerConfig::set_mismatch_transition_score(int8_t mm_transition) {
    mm_transition_ = mm_transition;
    score_matrix_['a']['g'] = score_matrix_['A']['g'] = score_matrix_['a']['G'] = score_matrix_['A']['G'] = mm_transition;
    score_matrix_['g']['a'] = score_matrix_['G']['a'] = score_matrix_['g']['A'] = score_matrix_['G']['A'] = mm_transition;
    score_matrix_['c']['t'] = score_matrix_['C']['t'] = score_matrix_['c']['T'] = score_matrix_['C']['T'] = mm_transition;
    score_matrix_['t']['c'] = score_matrix_['T']['c'] = score_matrix_['t']['C'] = score_matrix_['T']['C'] = mm_transition;
}

void DBGAlignerConfig::set_mismatch_transversion_score(int8_t mm_transversion) {
    mm_transversion_ = mm_transversion;

    // Set the default score to be the transversion penalty
    for (auto& row : score_matrix_) {
        row.fill(mm_transversion_);
    }

    // Uncomment these if the default value is anything other than the transversion penalty
    // score_matrix_['a']['c'] = score_matrix_['A']['c'] = score_matrix_['a']['C'] = score_matrix_['A']['C'] = mm_transversion;
    // score_matrix_['a']['t'] = score_matrix_['A']['t'] = score_matrix_['a']['T'] = score_matrix_['A']['T'] = mm_transversion;
    // score_matrix_['g']['c'] = score_matrix_['G']['c'] = score_matrix_['g']['C'] = score_matrix_['G']['C'] = mm_transversion;
    // score_matrix_['g']['t'] = score_matrix_['G']['t'] = score_matrix_['g']['T'] = score_matrix_['G']['T'] = mm_transversion;
    // score_matrix_['c']['a'] = score_matrix_['C']['a'] = score_matrix_['c']['A'] = score_matrix_['C']['A'] = mm_transversion;
    // score_matrix_['c']['g'] = score_matrix_['C']['g'] = score_matrix_['c']['G'] = score_matrix_['C']['G'] = mm_transversion;
    // score_matrix_['t']['a'] = score_matrix_['T']['a'] = score_matrix_['t']['A'] = score_matrix_['T']['A'] = mm_transversion;
    // score_matrix_['t']['g'] = score_matrix_['T']['g'] = score_matrix_['t']['G'] = score_matrix_['T']['G'] = mm_transversion;

    // Correct other scores
    set_match_score(match_score_);
    set_mismatch_transition_score(mm_transition_);
}







template <typename NodeType>
Alignment<NodeType>::Alignment(const DPTable &dp_table,
                               const typename DPTable::value_type *column,
                               size_t start_pos,
                               score_t score,
                               const char* path_end,
                               bool orientation)
          : query_begin_(NULL),
            query_end_(NULL),
            num_matches_(0),
            score_(score),
            orientation_(orientation) {
    assert(column != &*dp_table.end());

    auto [cigar_op, prev_node] = std::get<1>(column->second).at(start_pos);
    assert(cigar_op == Cigar::Operator::MATCH);

    auto i = start_pos;
    if (!i && prev_node == DeBruijnGraph::npos)
        return;

    std::vector<const typename DPTable::value_type*> out_columns;
    while (prev_node != DeBruijnGraph::npos) {
        cigar_.append(cigar_op);

        if (cigar_op != Cigar::Operator::INSERTION)
            out_columns.emplace_back(column);

        if (cigar_op != Cigar::Operator::DELETION)
            --i;

        if (cigar_op == Cigar::Operator::MATCH)
            num_matches_++;

        column = &*dp_table.find(prev_node);
        assert(column != &*dp_table.end());
        std::tie(cigar_op, prev_node) = std::get<1>(column->second).at(i);
    }

    if (UNLIKELY(i > std::numeric_limits<Cigar::LengthType>::max())) {
        throw std::runtime_error("Error: clipping length can't be stored in CIGAR");
    }

    cigar_.append(Cigar::Operator::CLIPPED, i);
    assert(cigar_.size());

    query_begin_ = path_end;
    query_end_ = path_end + start_pos;

    std::reverse(cigar_.begin(), cigar_.end());


    nodes_.resize(out_columns.size());
    std::transform(out_columns.rbegin(),
                   out_columns.rend(),
                   nodes_.begin(),
                   [](const auto &iter) { return iter->first; });

    sequence_.assign(out_columns.size(), 'N');
    std::transform(out_columns.rbegin(),
                   out_columns.rend(),
                   sequence_.begin(),
                   [](const auto &iter) { return std::get<2>(iter->second); });
}

template <typename NodeType>
void Alignment<NodeType>::append(Alignment&& other, size_t overlap, int8_t match_score) {
    assert(query_end_ - overlap == other.query_begin_);
    assert(orientation_ == other.orientation_);

    if (!overlap && other.cigar_.size()
            && other.cigar_.begin()->first == Cigar::Operator::CLIPPED) {
        // remove the seed if the extension has clipping
        assert(query_begin_ <= other.query_begin_);

        auto query_begin = query_begin_;
        *this = std::move(other);
        cigar_.begin()->second += query_begin_ - query_begin;
        query_begin_ = query_begin;
        return;
    }

    if (overlap) {
        assert(other.cigar_.size() && other.cigar_.begin()->first == Cigar::Operator::MATCH);
        assert(other.cigar_.size() && other.cigar_.begin()->second >= overlap);
        assert(other.sequence_.size() >= overlap);
        assert(match_score > 0);

        other.query_begin_ += overlap;
        if (other.cigar_.begin()->second == overlap) {
            other.cigar_.pop_front();
        } else {
            other.cigar_.begin()->second -= overlap;
        }

        other.sequence_ = other.sequence_.substr(overlap);
        other.score_ -= static_cast<score_t>(overlap) * match_score;
        other.num_matches_ -= overlap;
    }

    nodes_.insert(nodes_.end(),
                  other.nodes_.begin(),
                  other.nodes_.end());
    sequence_ += std::move(other.sequence_);
    score_ += other.score_;
    num_matches_ += other.num_matches_;

    cigar_.append(std::move(other.cigar_));
    query_end_ = other.query_end_;
}

template class Alignment<DeBruijnGraph::node_index>;
