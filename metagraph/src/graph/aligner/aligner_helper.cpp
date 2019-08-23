#include "aligner_helper.hpp"

#include "utils.hpp"
#include "alphabets.hpp"


std::pair<const char*, const char*>
get_alphabet_boundaries(const std::string &alphabet) {
    auto range = std::make_pair(&*alphabet.begin(), &*alphabet.end());

    // TODO: fix graph alphabets and these hacks
    if (*range.first == '$')
        ++range.first;

    // TODO: find a better solution for this hack
    // The last character in these alphabets is a wildcard ('N' or 'X')
    // make sure that wildcard-wildcard matches are counted as mismatches
    if (!strcmp(range.first, alphabets::kAlphabetProtein)
            || !strcmp(range.first, alphabets::kAlphabetDNA5))
        --range.second;

    return range;
}

Cigar::OperatorTable Cigar::char_to_op;

void Cigar::initialize_opt_table(const std::string &alphabet) {
    for (auto& row : char_to_op) {
        row.fill(Cigar::Operator::MISMATCH);
    }

    auto range = get_alphabet_boundaries(alphabet);

    for (const char *it = range.first; it != range.second; ++it) {
        char upper = toupper(*it);
        char lower = tolower(*it);

        char_to_op[upper][upper]
            = char_to_op[upper][lower]
            = char_to_op[lower][upper]
            = char_to_op[lower][lower] = Cigar::Operator::MATCH;
    }
}

char opt_to_char(Cigar::Operator op) {
    switch (op) {
        case Cigar::Operator::MATCH: return '=';
        case Cigar::Operator::MISMATCH: return 'X';
        case Cigar::Operator::INSERTION: return 'I';
        case Cigar::Operator::DELETION: return 'D';
        case Cigar::Operator::CLIPPED: return 'S';
    }

    assert(false);
    return '\0';
}

std::string Cigar::to_string() const {
    std::string cigar_string;

    for (const auto &pair : cigar_) {
        cigar_string += std::to_string(pair.second) + opt_to_char(pair.first);
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


DBGAlignerConfig::DBGAlignerConfig(const ScoreMatrix &score_matrix,
                                   int8_t gap_opening,
                                   int8_t gap_extension)
      : gap_opening_penalty(gap_opening),
        gap_extension_penalty(gap_extension),
        score_matrix_(score_matrix) { }

DBGAlignerConfig::DBGAlignerConfig(ScoreMatrix&& score_matrix,
                                   int8_t gap_opening,
                                   int8_t gap_extension)
      : gap_opening_penalty(gap_opening),
        gap_extension_penalty(gap_extension),
        score_matrix_(std::move(score_matrix)) { }

DBGAlignerConfig::DBGAlignerConfig(const Config &config, const DeBruijnGraph &graph)
      : queue_size(config.alignment_queue_size),
        num_alternative_paths(config.alignment_num_alternative_paths),
        min_seed_length(config.alignment_min_seed_length),
        max_seed_length(config.alignment_max_seed_length),
        max_num_seeds_per_locus(config.alignment_max_num_seeds_per_locus),
        min_cell_score(config.alignment_min_cell_score),
        gap_opening_penalty(-config.alignment_gap_opening_penalty),
        gap_extension_penalty(-config.alignment_gap_extension_penalty),
        score_matrix_(config.alignment_edit_distance
                          ? unit_scoring_matrix(1, graph.alphabet())
                          : scoring_matrix(config, graph)) { }

DBGAlignerConfig::ScoreMatrix DBGAlignerConfig
::scoring_matrix(const Config &config, const DeBruijnGraph &graph) {
    auto range = get_alphabet_boundaries(graph.alphabet());

    if (std::equal(range.first, range.second, alphabets::kAlphabetDNA)) {
        return dna_scoring_matrix(config.alignment_match_score,
                                  -config.alignment_mm_transition,
                                  -config.alignment_mm_transversion);
    } else if (std::equal(range.first, range.second, alphabets::kAlphabetProtein)) {
        return score_matrix_blosum62;
    }

    return unit_scoring_matrix(1, graph.alphabet());
}

DBGAlignerConfig::ScoreMatrix DBGAlignerConfig
::dna_scoring_matrix(int8_t match_score,
                     int8_t mm_transition,
                     int8_t mm_transversion) {
    ScoreMatrix score_matrix;
    for (auto& row : score_matrix) {
        row.fill(mm_transversion);
    }

    score_matrix['a']['g'] = score_matrix['A']['g'] = score_matrix['a']['G'] = score_matrix['A']['G'] = mm_transition;
    score_matrix['g']['a'] = score_matrix['G']['a'] = score_matrix['g']['A'] = score_matrix['G']['A'] = mm_transition;
    score_matrix['c']['t'] = score_matrix['C']['t'] = score_matrix['c']['T'] = score_matrix['C']['T'] = mm_transition;
    score_matrix['t']['c'] = score_matrix['T']['c'] = score_matrix['t']['C'] = score_matrix['T']['C'] = mm_transition;
    score_matrix['a']['a'] = score_matrix['A']['a'] = score_matrix['a']['A'] = score_matrix['A']['A'] = match_score;
    score_matrix['c']['c'] = score_matrix['C']['c'] = score_matrix['c']['C'] = score_matrix['C']['C'] = match_score;
    score_matrix['g']['g'] = score_matrix['G']['g'] = score_matrix['g']['G'] = score_matrix['G']['G'] = match_score;
    score_matrix['t']['t'] = score_matrix['T']['t'] = score_matrix['t']['T'] = score_matrix['T']['T'] = match_score;

    return score_matrix;
}

DBGAlignerConfig::ScoreMatrix DBGAlignerConfig
::unit_scoring_matrix(int8_t match_score, const std::string &alphabet) {
    ScoreMatrix score_matrix;
    for (auto& row : score_matrix) {
        row.fill(-match_score);
    }

    auto range = get_alphabet_boundaries(alphabet);

    for (const char *it = range.first; it != range.second; ++it) {
        char upper = toupper(*it);
        char lower = tolower(*it);

        score_matrix[upper][upper]
            = score_matrix[upper][lower]
            = score_matrix[lower][upper]
            = score_matrix[lower][lower] = match_score;
    }

    return score_matrix;
}

DBGAlignerConfig::ScoreMatrix blosum62_scoring_matrix() {
    std::string alphabet = "ARNDCQEGHILKMFPSTWYVBZX";

    std::vector<std::vector<int8_t>> scores = {
        {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0 },
        { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1 },
        { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1 },
        { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1 },
        {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2 },
        { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1 },
        { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1 },
        {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1 },
        { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1 },
        { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1 },
        { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1 },
        { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1 },
        { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1 },
        { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1 },
        { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2 },
        {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0 },
        {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0 },
        { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2 },
        { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1 },
        {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1 },
        { -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1 },
        { -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1 },
        {  0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1 }
    };

    DBGAlignerConfig::ScoreMatrix score_matrix;

    for (size_t i = 0; i < score_matrix.size(); ++i) {
        score_matrix[i].fill(-4);

        // meant to handle the letters J, O, U
        score_matrix[i][i] = 1;
    }

    for (size_t i = 0; i < alphabet.size(); ++i) {
        char a_upper = alphabet[i];
        char a_lower = tolower(alphabet[i]);

        for (size_t j = 0; j < alphabet.size(); ++j) {
            char b_upper = alphabet[j];
            char b_lower = tolower(alphabet[j]);

            score_matrix[a_lower][b_lower]
                = score_matrix[a_lower][b_upper]
                = score_matrix[a_upper][b_lower]
                = score_matrix[a_upper][b_upper] = scores[i][j];
        }
    }

    return score_matrix;
}

const DBGAlignerConfig::ScoreMatrix DBGAlignerConfig::score_matrix_blosum62
    = blosum62_scoring_matrix();

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

template <typename NodeType>
bool Alignment<NodeType>::is_valid(const DeBruijnGraph &graph) const {
    // TODO: check score

    if (query_begin_ > query_end_) {
        std::cerr << "ERROR: query begin after query end" << std::endl
                  << *this << std::endl;
        return false;
    }

    auto query_it = query_begin_;
    auto node_it = nodes_.begin();
    auto cigar_it = cigar_.begin() + static_cast<bool>(cigar_.get_clipping());

    size_t ref_counter = offset_;
    std::string path;
    if (nodes_.size())
        path = graph.get_node_sequence(nodes_.front()).substr(offset_);

    size_t path_steps = 0;

    for (; cigar_it != cigar_.end(); ++cigar_it) {
        if (query_it >= query_end_ && cigar_it->first != Cigar::Operator::DELETION) {
            std::cerr << "ERROR: end of query reached before end of CIGAR" << std::endl
                      << "Processed " << cigar_it - cigar_.begin()
                      << " of " << cigar_.size() << " operations" << std::endl
                      << *this << std::endl
                      << std::string(query_begin_, query_end_) << std::endl;
            return false;
        }

        if (node_it == nodes_.end() && cigar_it->first != Cigar::Operator::INSERTION) {
            std::cerr << "ERROR: end of nodes reached before end of CIGAR" << std::endl
                      << "Processed " << cigar_it - cigar_.begin()
                      << " of " << cigar_.size() << " operations" << std::endl
                      << *this << std::endl
                      << std::string(query_begin_, query_end_) << std::endl;
            return false;
        }

        switch (cigar_it->first) {
            case Cigar::Operator::MATCH:
            case Cigar::Operator::MISMATCH:
            case Cigar::Operator::DELETION: {
                auto cur_query_it = query_it;
                auto cur_path_steps = path_steps;

                if (cigar_it->first != Cigar::Operator::DELETION)
                    query_it += cigar_it->second;

                path_steps += cigar_it->second;
                size_t old_ref_counter = ref_counter;
                ref_counter += cigar_it->second;

                if (ref_counter >= graph.get_k()) {
                    size_t shift = 0;
                    if (old_ref_counter < graph.get_k()) {
                        assert(ref_counter - graph.get_k() + 1 <= cigar_it->second);
                        shift = ref_counter - graph.get_k() + 1;
                    } else {
                        shift = cigar_it->second;
                    }

                    for (size_t i = 0; i < shift; ++i) {
                        ++node_it;
                        if (node_it != nodes_.end()) {
                            graph.call_outgoing_kmers(
                                *(node_it - 1),
                                [&](auto node, char c) {
                                    if (node == *node_it)
                                        path += c;
                                }
                            );
                        }
                    }
                }

                if (cigar_it->first == Cigar::Operator::MATCH
                        && strncmp(cur_query_it,
                                   path.c_str() + cur_path_steps,
                                   path_steps - cur_path_steps)) {
                    std::cerr << "ERROR: mismatch found despite MATCH in CIGAR"
                              << *this << std::endl
                              << std::string(query_begin_, query_end_) << std::endl;
                    return false;
                } else if (cigar_it->first == Cigar::Operator::MISMATCH
                        && std::mismatch(cur_query_it,
                                         query_it,
                                         path.c_str() + cur_path_steps,
                                         path.c_str() + path_steps).first != cur_query_it) {
                    std::cerr << "ERROR: match found despite MISMATCH in CIGAR"
                              << *this << std::endl
                              << std::string(query_begin_, query_end_) << std::endl;
                    return false;
                }
            } break;
            case Cigar::Operator::INSERTION: query_it += cigar_it->second; break;
            case Cigar::Operator::CLIPPED: {
                std::cerr << "ERROR: internal clipping in CIGAR" << std::endl
                          << *this << std::endl
                          << std::string(query_begin_, query_end_) << std::endl;
                return false;
            }
        }
    }

    path = path.substr(0, path_steps);

    if (path_steps != sequence_.size()) {
        std::cerr << "ERROR: stored sequence is incorrect size" << std::endl
                  << *this << std::endl
                  << std::string(query_begin_, query_end_) << std::endl;
        return false;
    }

    if (path.size() != sequence_.size()) {
        std::cerr << "ERROR: stored sequence is incorrect" << std::endl
                  << "Reconstructed sequence: " << path << std::endl
                  << *this << std::endl
                  << std::string(query_begin_, query_end_) << std::endl;
        return false;
    }

    if (query_it != query_end_) {
        std::cerr << "ERROR: end of CIGAR reached before end of query" << std::endl
                  << "Processed " << query_it - query_begin_
                  << " of " << query_end_ - query_begin_ << " characters" << std::endl
                  << *this << std::endl
                  << std::string(query_begin_, query_end_) << std::endl;
        return false;
    }

    if (node_it != nodes_.end()) {
        std::cerr << "ERROR: end of CIGAR reached before end of path" << std::endl
                  << "Processed " << node_it - nodes_.begin()
                  << " of " << nodes_.size() << " nodes" << std::endl
                  << *this << std::endl
                  << std::string(query_begin_, query_end_) << std::endl;
        return false;
    }

    return true;
}

template class Alignment<DeBruijnGraph::node_index>;
