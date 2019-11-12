#include "aligner_helper.hpp"

#include "algorithms.hpp"
#include "alphabets.hpp"


Cigar::Cigar(const std::string &cigar_str) {
    std::string op_count;
    for (auto c : cigar_str) {
        switch (c) {
            case '=':
                cigar_.emplace_back(Cigar::Operator::MATCH, std::stol(op_count));
                op_count.clear();
                break;
            case 'X':
                cigar_.emplace_back(Cigar::Operator::MISMATCH, std::stol(op_count));
                op_count.clear();
                break;
            case 'I':
                cigar_.emplace_back(Cigar::Operator::INSERTION, std::stol(op_count));
                op_count.clear();
                break;
            case 'D':
                cigar_.emplace_back(Cigar::Operator::DELETION, std::stol(op_count));
                op_count.clear();
                break;
            case 'S':
                cigar_.emplace_back(Cigar::Operator::CLIPPED, std::stol(op_count));
                op_count.clear();
                break;
            default:
                op_count += c;
        }
    }
}

Cigar::OperatorTable Cigar::char_to_op;

void Cigar::initialize_opt_table(const std::string &alphabet, const uint8_t *encoding) {
    for (auto& row : char_to_op) {
        row.fill(Cigar::Operator::MISMATCH);
    }

    for (uint8_t c : alphabet) {
        if (encoding[c] == encoding[0])
            continue;

        char upper = toupper(c);
        char lower = tolower(c);

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

bool Cigar::is_valid(const char *reference_begin, const char *reference_end,
                     const char *query_begin, const char *query_end) const {
    const char *ref_it = reference_begin;
    const char *alt_it = query_begin;

    for (const auto &op : cigar_) {
        switch (op.first) {
            case Operator::CLIPPED: {
                if (ref_it != reference_begin || alt_it != query_begin) {
                    std::cerr << "Internal clipping found in CIGAR" << std::endl
                              << to_string() << std::endl
                              << std::string(reference_begin, reference_end) << std::endl
                              << std::string(query_begin, query_end) << std::endl;
                    return false;
                }
            } break;
            case Operator::MATCH: {
                if (ref_it > reference_end - op.second) {
                    std::cerr << "Reference too short" << std::endl
                              << to_string() << std::endl
                              << std::string(reference_begin, reference_end) << std::endl
                              << std::string(query_begin, query_end) << std::endl;
                    return false;
                }

                if (alt_it > query_end - op.second) {
                    std::cerr << "Query too short" << std::endl
                              << to_string() << std::endl
                              << std::string(reference_begin, reference_end) << std::endl
                              << std::string(query_begin, query_end) << std::endl;
                    return false;
                }

                if (strncmp(ref_it, alt_it, op.second)) {
                    std::cerr << "Mismatch despite MATCH in CIGAR" << std::endl
                              << to_string() << std::endl
                              << std::string(reference_begin, reference_end) << std::endl
                              << std::string(query_begin, query_end) << std::endl;
                    return false;
                }

                ref_it += op.second;
                alt_it += op.second;
            } break;
            case Operator::MISMATCH: {
                if (ref_it > reference_end - op.second) {
                    std::cerr << "Reference too short" << std::endl
                              << to_string() << std::endl
                              << std::string(reference_begin, reference_end) << std::endl
                              << std::string(query_begin, query_end) << std::endl;
                    return false;
                }

                if (alt_it > query_end - op.second) {
                    std::cerr << "Query too short" << std::endl
                              << to_string() << std::endl
                              << std::string(reference_begin, reference_end) << std::endl
                              << std::string(query_begin, query_end) << std::endl;
                    return false;
                }

                if (std::mismatch(ref_it, ref_it + op.second,
                                  alt_it, alt_it + op.second).first != ref_it) {
                    std::cerr << "Match despite MISMATCH in CIGAR" << std::endl
                              << to_string() << std::endl
                              << std::string(reference_begin, reference_end) << std::endl
                              << std::string(query_begin, query_end) << std::endl;
                    return false;
                }

                ref_it += op.second;
                alt_it += op.second;
            } break;
            case Operator::INSERTION: {
                if (alt_it > query_end - op.second) {
                    std::cerr << "Query too short" << std::endl
                              << to_string() << std::endl
                              << std::string(reference_begin, reference_end) << std::endl
                              << std::string(query_begin, query_end) << std::endl;
                    return false;
                }

                alt_it += op.second;
            } break;
            case Operator::DELETION: {
                if (ref_it > reference_end - op.second) {
                    std::cerr << "Reference too short" << std::endl
                              << to_string() << std::endl
                              << std::string(reference_begin, reference_end) << std::endl
                              << std::string(query_begin, query_end) << std::endl;
                    return false;
                }

                ref_it += op.second;
            } break;
        }
    }

    if (ref_it != reference_end) {
        std::cerr << "Reference end not reached" << std::endl
                  << to_string() << std::endl
                  << std::string(reference_begin, reference_end) << std::endl
                  << std::string(query_begin, query_end) << std::endl;
        return false;
    }

    if (alt_it != query_end) {
        std::cerr << "Query end not reached" << std::endl
                  << to_string() << std::endl
                  << std::string(reference_begin, reference_end) << std::endl
                  << std::string(query_begin, query_end) << std::endl;
        return false;
    }

    return true;
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

DBGAlignerConfig::DBGAlignerConfig(const Config &config)
      : queue_size(config.alignment_queue_size),
        bandwidth(config.alignment_vertical_bandwidth),
        num_alternative_paths(config.alignment_num_alternative_paths),
        min_seed_length(config.alignment_min_seed_length),
        max_seed_length(config.alignment_max_seed_length),
        max_num_seeds_per_locus(config.alignment_max_num_seeds_per_locus),
        min_cell_score(config.alignment_min_cell_score),
        min_path_score(config.alignment_min_path_score),
        gap_opening_penalty(-config.alignment_gap_opening_penalty),
        gap_extension_penalty(-config.alignment_gap_extension_penalty),
        forward_and_reverse_complement(config.forward_and_reverse),
        score_matrix_(scoring_matrix(config)) {}

DBGAlignerConfig::score_t DBGAlignerConfig
::score_cigar(const char *reference_begin, const char *reference_end,
              const char *query_begin, const char *query_end,
              const Cigar &cigar) const {
    score_t score = 0;

    assert(cigar.is_valid(reference_begin, reference_end, query_begin, query_end));
    std::ignore = reference_end;
    std::ignore = query_end;

    for (const auto &op : cigar) {
        switch (op.first) {
            case Cigar::Operator::CLIPPED:
                break;
            case Cigar::Operator::MATCH: {
                score += match_score(reference_begin, reference_begin + op.second);
                reference_begin += op.second;
                query_begin += op.second;
            } break;
            case Cigar::Operator::MISMATCH: {
                score += score_sequences(reference_begin,
                                         reference_begin + op.second,
                                         query_begin);
                reference_begin += op.second;
                query_begin += op.second;
            } break;
            case Cigar::Operator::INSERTION: {
                score += gap_opening_penalty + (op.second - 1) * gap_extension_penalty;
                query_begin += op.second;
            } break;
            case Cigar::Operator::DELETION: {
                score += gap_opening_penalty + (op.second - 1) * gap_extension_penalty;
                reference_begin += op.second;
            } break;
        }
    }

    return score;
}

DBGAlignerConfig::ScoreMatrix DBGAlignerConfig
::scoring_matrix(const Config &config) {
    if (config.alignment_edit_distance) {
        // TODO: REPLACE THIS
        #if _PROTEIN_GRAPH
            const auto *alphabet = alphabets::kAlphabetProtein;
            const auto *alphabet_encoding = alphabets::kCharToProtein;
        #elif _DNA_CASE_SENSITIVE_GRAPH
            const auto *alphabet = alphabets::kAlphabetDNA;
            const auto *alphabet_encoding = alphabets::kCharToDNA;
        #elif _DNA5_GRAPH
            const auto *alphabet = alphabets::kAlphabetDNA;
            const auto *alphabet_encoding = alphabets::kCharToDNA;
        #elif _DNA_GRAPH
            const auto *alphabet = alphabets::kAlphabetDNA;
            const auto *alphabet_encoding = alphabets::kCharToDNA;
        #else
            static_assert(false,
                "Define an alphabet: either "
                "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
            );
        #endif

        return unit_scoring_matrix(1, alphabet, alphabet_encoding);
    }

    #if _PROTEIN_GRAPH
        return score_matrix_blosum62;
    #elif _DNA_CASE_SENSITIVE_GRAPH
        return dna_scoring_matrix(config.alignment_match_score,
                                  -config.alignment_mm_transition,
                                  -config.alignment_mm_transversion);
    #elif _DNA5_GRAPH
        return dna_scoring_matrix(config.alignment_match_score,
                                  -config.alignment_mm_transition,
                                  -config.alignment_mm_transversion);
    #elif _DNA_GRAPH
        return dna_scoring_matrix(config.alignment_match_score,
                                  -config.alignment_mm_transition,
                                  -config.alignment_mm_transversion);
    #else
        static_assert(false,
            "Define an alphabet: either "
            "_DNA_GRAPH, _DNA5_GRAPH, _PROTEIN_GRAPH, or _DNA_CASE_SENSITIVE_GRAPH."
        );
    #endif
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
::unit_scoring_matrix(int8_t match_score,
                      const std::string &alphabet,
                      const uint8_t *encoding) {
    ScoreMatrix score_matrix;
    for (auto& row : score_matrix) {
        row.fill(-match_score);
    }

    for (uint8_t c : alphabet) {
        if (encoding[c] == encoding[0])
            continue;

        char upper = toupper(c);
        char lower = tolower(c);

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
Alignment<NodeType>::Alignment(const char* query_begin,
                               const char* query_end,
                               std::vector<NodeType>&& nodes,
                               std::string&& sequence,
                               score_t score,
                               size_t clipping,
                               bool orientation,
                               size_t offset)
      : query_begin_(query_begin),
        query_end_(query_end),
        nodes_(std::move(nodes)),
        sequence_(std::move(sequence)),
        num_matches_(0),
        score_(score),
        orientation_(orientation),
        offset_(offset) {
    assert(query_begin_);
    assert(query_end_);

    size_t query_size = query_end_ - query_begin_;
    size_t min_length = std::min(query_size, sequence_.size());

    cigar_ = std::inner_product(
        query_begin_,
        query_begin_ + min_length,
        sequence_.c_str(),
        Cigar(Cigar::Operator::CLIPPED, clipping),
        [&](Cigar &cigar, bool equal) -> Cigar& {
            num_matches_ += equal;
            cigar.append(equal
                  ? Cigar::Operator::MATCH
                  : Cigar::Operator::MISMATCH
              );
            return cigar;
        },
        std::equal_to<char>()
    );

    assert(!(query_size - min_length) || (sequence_.size() - min_length));
    cigar_.append(Cigar::Operator::INSERTION, query_size - min_length);
    cigar_.append(Cigar::Operator::DELETION, sequence_.size() - min_length);
}

template <typename NodeType>
Alignment<NodeType>::Alignment(const DPTable &dp_table,
                               const typename DPTable::value_type *column,
                               size_t start_pos,
                               score_t score,
                               const char* path_end,
                               bool orientation,
                               size_t offset)
      : query_begin_(NULL),
        query_end_(NULL),
        num_matches_(0),
        score_(score),
        orientation_(orientation),
        offset_(offset) {
    assert(column);

    auto i = start_pos;
    const auto* step = &column->second.steps.at(i);

    if (!i && step->prev_node == DeBruijnGraph::npos)
        return;

    std::vector<const typename DPTable::value_type*> out_columns;
    while (step->prev_node != DeBruijnGraph::npos) {
        cigar_.append(step->cigar_op);

        if (step->cigar_op != Cigar::Operator::INSERTION)
            out_columns.emplace_back(column);

        if (step->cigar_op != Cigar::Operator::DELETION)
            --i;

        if (step->cigar_op == Cigar::Operator::MATCH)
            num_matches_++;

        auto find = dp_table.find(step->prev_node);
        assert(find != dp_table.end());
        column = &*find;

        step = &column->second.steps.at(i);
    }

    if (UNLIKELY(i > std::numeric_limits<Cigar::LengthType>::max())) {
        throw std::runtime_error("Error: clipping length can't be stored in CIGAR");
    }

    cigar_.append(Cigar::Operator::CLIPPED, i);
    assert(cigar_.size());

    query_begin_ = path_end + i;
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
                   [](const auto &iter) { return iter->second.last_char; });
}

template <typename NodeType>
void Alignment<NodeType>::append(Alignment&& other, size_t overlap, int8_t match_score) {
    assert(query_end_ - overlap + other.get_clipping() == other.query_begin_);
    assert(orientation_ == other.orientation_);

    if (!overlap && other.cigar_.get_clipping()) {
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
void Alignment<NodeType>::recompute_score(const DBGAlignerConfig &config) {
    auto new_score = config.score_cigar(sequence_.c_str(),
                                        sequence_.c_str() + sequence_.size(),
                                        query_begin_,
                                        query_end_,
                                        cigar_);

    if (utils::get_verbose() && score_ != new_score)
        std::cerr << "WARNING: changing score from "
                  << score_ << " to " << new_score << std::endl;

    score_ = new_score;
}

// derived from:
// https://github.com/maickrau/GraphAligner/blob/236e1cf0514cfa9104e9a3333cdc1c43209c3c5a/src/vg.proto
template <typename NodeType>
Json::Value Alignment<NodeType>::path_json(size_t node_size,
                                           const std::string &label) const {
    assert(nodes_.size());

    Json::Value path;

    auto cigar_it = cigar_.begin();
    if (cigar_.size() && cigar_it->first == Cigar::Operator::CLIPPED) {
        cigar_it++;
    }

    size_t cigar_offset = 0;
    assert(cigar_it != cigar_.end());

    int64_t rank = 1;
    auto query_start = query_begin_;

    size_t cur_pos = rank == 1 ? offset_ : 0;

    Json::Value mapping;
    Json::Value position;
    position["node_id"] = nodes_.front();

    if (cur_pos)
        position["offset"] = Json::Value::UInt64(cur_pos);

    // set to true if the node is the reverse complement of the query
    //position["is_reverse"] = false;

    mapping["position"] = position;

    while (cur_pos < node_size && cigar_it != cigar_.end()) {
        assert(cigar_it->second > cigar_offset);
        size_t next_pos = std::min(node_size,
                                   cur_pos + (cigar_it->second - cigar_offset));
        size_t next_size = next_pos - cur_pos;
        assert(cigar_offset + next_size <= cigar_it->second);

        Json::Value edit;
        switch (cigar_it->first) {
            case Cigar::Operator::MISMATCH: {
                assert(query_start + next_size <= query_end_);
                edit["from_length"] = Json::Value::UInt64(next_size);
                edit["to_length"] = Json::Value::UInt64(next_size);
                edit["sequence"] = std::string(query_start, next_size);
                query_start += next_size;
            } break;
            case Cigar::Operator::INSERTION: {
                assert(query_start + next_size <= query_end_);
                // this assumes that INSERTIONS can't happen right after deletions
                //edit["from_length"] = 0;
                edit["to_length"] = Json::Value::UInt64(next_size);
                edit["sequence"] = std::string(query_start, next_size);
                query_start += next_size;
            } break;
            case Cigar::Operator::DELETION: {
                edit["from_length"] = Json::Value::UInt64(next_size);
                //edit["to_length"] = 0;
            } break;
            case Cigar::Operator::MATCH: {
                edit["from_length"] = Json::Value::UInt64(next_size);
                edit["to_length"] = Json::Value::UInt64(next_size);
                query_start += next_size;
            } break;
            case Cigar::Operator::CLIPPED: { assert(false); }
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

    for (auto node_it = nodes_.begin() + 1; node_it != nodes_.end(); ++node_it) {
        assert(cigar_it != cigar_.end());
        assert(cigar_it->second > cigar_offset);

        Json::Value mapping;
        Json::Value position;
        position["node_id"] = *node_it;
        position["offset"] = Json::Value::UInt64(node_size - 1);
        // set to true if the node is the reverse complement of the query
        //position["is_reverse"] = false;
        mapping["position"] = position;

        if (cigar_it->first == Cigar::Operator::INSERTION) {
            Json::Value edit;
            size_t length = cigar_it->second - cigar_offset;
            assert(query_start + length < query_end_);
            // this assumes that INSERTIONS can't happen right after deletions
            //edit["from_length"] = 0;
            edit["to_length"] = Json::Value::UInt64(length);
            edit["sequence"] = std::string(query_start, length);
            query_start += length;
            ++cigar_it;
            cigar_offset = 0;
            mapping["edit"].append(edit);
            assert(cigar_it != cigar_.end());
        }

        Json::Value edit;
        switch (cigar_it->first) {
            case Cigar::Operator::MISMATCH: {
                assert(query_start < query_end_);
                edit["from_length"] = 1;
                edit["to_length"] = 1;
                edit["sequence"] = std::string(query_start, 1);
                query_start++;
            } break;
            case Cigar::Operator::DELETION: {
                edit["from_length"] = 1;
                //edit["to_length"] = 0;
            } break;
            case Cigar::Operator::MATCH: {
                edit["from_length"] = 1;
                edit["to_length"] = 1;
                query_start++;
            } break;
            case Cigar::Operator::INSERTION:
            case Cigar::Operator::CLIPPED: { assert(false); }
        }

        if (++cigar_offset == cigar_it->second) {
            cigar_offset = 0;
            ++cigar_it;
        }

        mapping["edit"].append(edit);
        mapping["rank"] = rank++;
        path["mapping"].append(mapping);
    }

    assert(query_start == query_end_);
    assert(cigar_it == cigar_.end());

    path["length"] = Json::Value::UInt64(nodes_.size());
    //path["is_circular"]; // bool

    if (label.size())
        path["name"] = label;

    return path;
}

template <typename NodeType>
Json::Value Alignment<NodeType>::to_json(const std::string &query,
                                         const DeBruijnGraph &graph,
                                         bool is_secondary,
                                         const std::string &read_name,
                                         const std::string &label) const {
    assert(is_valid(graph));

    // encode alignment
    Json::Value alignment;

    alignment["name"] = read_name;
    alignment["sequence"] = query;

    if (sequence_.size())
        alignment["annotation"]["ref_sequence"] = sequence_;

    if (query_end_ == query_begin_)
        return alignment;

    auto query_start = query.c_str();
    assert(query_start + cigar_.get_clipping() == query_begin_);
    assert(query_end_ >= query_start);

    alignment["annotation"]["cigar"] = cigar_.to_string();

    // encode path
    if (nodes_.size())
        alignment["path"] = path_json(graph.get_k(), label);

    alignment["score"] = static_cast<int32_t>(score_);

    if (query_begin_ > query_start)
        alignment["query_position"] = static_cast<int32_t>(query_begin_ - query_start);

    if (is_secondary)
        alignment["is_secondary"] = is_secondary;

    alignment["identity"] = query_end_ != query_begin_
        ? static_cast<double>(num_matches_) / (query_end_ - query_begin_)
        : 0;

    alignment["read_mapped"] = (query_end_ != query_begin_);

    if (orientation_)
        alignment["read_on_reverse_strand"] = orientation_;


    if (cigar_.get_clipping())
        alignment["soft_clipped"] = static_cast<bool>(cigar_.get_clipping());

    // Unused flags (for now)
    //alignment["quality"]; // bytes
    //alignment["mapping_quality"]; // int32
    //alignment["sample_name"]; // string
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

template <typename NodeType>
std::shared_ptr<const std::string> Alignment<NodeType>
::load_from_json(const Json::Value &alignment,
                 const DeBruijnGraph &graph) {
    cigar_.clear();
    nodes_.clear();
    num_matches_ = 0;
    sequence_.clear();

    auto query_sequence = std::make_shared<const std::string>(
        alignment["sequence"].asString()
    );
    auto query_start = query_sequence->c_str();

    query_begin_ = query_start + alignment["query_position"].asInt();
    orientation_ = alignment["read_on_reverse_strand"].asBool();
    score_ = alignment["score"].asInt();

    const auto& mapping = alignment["path"]["mapping"];
    assert(mapping.size() == alignment["path"]["length"].asUInt64());

    Json::ArrayIndex i = 0;
    offset_ = mapping[i]["position"]["offset"].asUInt64();

    size_t path_steps = 0;

    query_end_ = query_begin_;

    cigar_.append(Cigar::Operator::CLIPPED, query_begin_ - query_start);
    for (; i < mapping.size(); ++i) {
        nodes_.emplace_back(mapping[i]["position"]["node_id"].asUInt64());
        if (nodes_.size() == 1) {
            sequence_ = graph.get_node_sequence(nodes_.back()).substr(offset_);
        } else {
            graph.call_outgoing_kmers(
                *(nodes_.rbegin() + 1),
                [&](auto node, char c) {
                    if (node == nodes_.back())
                        sequence_ += c;
                }
            );
        }
        const auto& edits = mapping[i]["edit"];

        for (Json::ArrayIndex j = 0; j < edits.size(); ++j) {
            assert(edits[j]["from_length"].asUInt64()
                       || edits[j]["to_length"].asUInt64());

            if (edits[j]["from_length"] == edits[j]["to_length"]) {
                if (edits[j]["sequence"].asString().size()) {
                    cigar_.append(Cigar::Operator::MISMATCH,
                                  edits[j]["from_length"].asUInt64());
                } else {
                    cigar_.append(Cigar::Operator::MATCH,
                                  edits[j]["from_length"].asUInt64());
                    num_matches_ += edits[j]["from_length"].asUInt64();
                }

                path_steps += edits[j]["from_length"].asUInt64();
                query_end_ += edits[j]["to_length"].asUInt64();
            } else if (edits[j]["from_length"].asUInt64()) {
                cigar_.append(Cigar::Operator::DELETION,
                              edits[j]["from_length"].asUInt64());

                path_steps += edits[j]["from_length"].asUInt64();
            } else {
                assert(edits[j]["to_length"].asUInt64());
                cigar_.append(Cigar::Operator::INSERTION,
                              edits[j]["to_length"].asUInt64());

                query_end_ += edits[j]["to_length"].asUInt64();
            }
        }
    }

    if (alignment["annotation"]["cigar"]) {
        assert(alignment["annotation"]["cigar"].asString() == cigar_.to_string());

        if (cigar_ != Cigar(alignment["annotation"]["cigar"].asString()))
            throw std::runtime_error("ERROR: CIGAR and mapping mismatch");
    }

    sequence_ = sequence_.substr(0, path_steps);
    assert(!alignment["annotation"]["ref_sequence"]
        || alignment["annotation"]["ref_sequence"].asString() == sequence_);

    if (!is_valid(graph))
        throw std::runtime_error("ERROR: JSON reconstructs invalid alignment");

    return query_sequence;
}

template <typename NodeType>
bool Alignment<NodeType>::is_valid(const DeBruijnGraph &graph,
                                   const DBGAlignerConfig *config) const {
    if (query_begin_ > query_end_) {
        std::cerr << "ERROR: query begin after query end" << std::endl
                  << *this << std::endl;
        return false;
    }

    if (!cigar_.is_valid(sequence_.c_str(), sequence_.c_str() + sequence_.size(),
                         query_begin_, query_end_)) {
        std::cerr << "ERROR: CIGAR invalid" << std::endl;
        return false;
    }

    if (config && score_ != config->score_cigar(sequence_.c_str(),
                                                sequence_.c_str() + sequence_.size(),
                                                query_begin_,
                                                query_end_,
                                                cigar_)) {
        std::cerr << "ERROR: mismatch between CIGAR and score" << std::endl
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
                      << *this << std::endl;
            return false;
        }

        if (node_it == nodes_.end() && cigar_it->first != Cigar::Operator::INSERTION) {
            std::cerr << "ERROR: end of nodes reached before end of CIGAR" << std::endl
                      << "Processed " << cigar_it - cigar_.begin()
                      << " of " << cigar_.size() << " operations" << std::endl
                      << *this << std::endl;
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
                            bool node_found = false;
                            graph.call_outgoing_kmers(
                                *(node_it - 1),
                                [&](auto node, char c) {
                                    if (node == *node_it) {
                                        path += c;
                                        node_found = true;
                                    }
                                }
                            );

                            if (!node_found) {
                                std::cerr << "ERROR: invalid node index "
                                          << *node_it << std::endl
                                          << "Processed " << cigar_it - cigar_.begin()
                                          << " of " << cigar_.size() << " operations" << std::endl
                                          << *this << std::endl;
                                return false;
                            }
                        }
                    }
                }

                if (cigar_it->first == Cigar::Operator::MATCH
                        && strncmp(cur_query_it,
                                   path.c_str() + cur_path_steps,
                                   path_steps - cur_path_steps)) {
                    std::cerr << "ERROR: mismatch found despite MATCH in CIGAR" << std::endl
                              << "Processed " << cigar_it - cigar_.begin()
                              << " of " << cigar_.size() << " operations" << std::endl
                              << *this << "\n\n"
                              << std::string(cur_query_it, path_steps - cur_path_steps) << "\n"
                              << std::string(path.c_str() + cur_path_steps, path_steps - cur_path_steps) << std::endl;
                    return false;
                } else if (cigar_it->first == Cigar::Operator::MISMATCH
                        && std::mismatch(cur_query_it,
                                         query_it,
                                         path.c_str() + cur_path_steps,
                                         path.c_str() + path_steps).first != cur_query_it) {
                    std::cerr << "ERROR: match found despite MISMATCH in CIGAR" << std::endl
                              << "Processed " << cigar_it - cigar_.begin()
                              << " of " << cigar_.size() << " operations" << std::endl
                              << *this << "\n\n"
                              << std::string(cur_query_it, path_steps - cur_path_steps) << "\n"
                              << std::string(path.c_str() + cur_path_steps, path_steps - cur_path_steps) << std::endl;
                    return false;
                }
            } break;
            case Cigar::Operator::INSERTION: query_it += cigar_it->second; break;
            case Cigar::Operator::CLIPPED: {
                std::cerr << "ERROR: internal clipping in CIGAR" << std::endl
                          << *this << std::endl;
                return false;
            }
        }
    }

    path = path.substr(0, path_steps);

    if (path_steps != sequence_.size()) {
        std::cerr << "ERROR: stored sequence is incorrect size" << std::endl
                  << *this << std::endl;
        return false;
    }

    if (path.size() != sequence_.size()) {
        std::cerr << "ERROR: stored sequence is incorrect" << std::endl
                  << "Reconstructed sequence: " << path << std::endl
                  << *this << std::endl;
        return false;
    }

    if (query_it != query_end_) {
        std::cerr << "ERROR: end of CIGAR reached before end of query" << std::endl
                  << "Processed " << query_it - query_begin_
                  << " of " << query_end_ - query_begin_ << " characters" << std::endl
                  << *this << std::endl;
        return false;
    }

    if (node_it != nodes_.end()) {
        std::cerr << "ERROR: end of CIGAR reached before end of path" << std::endl
                  << "Processed " << node_it - nodes_.begin()
                  << " of " << nodes_.size() << " nodes" << std::endl
                  << *this << std::endl;
        return false;
    }

    return true;
}

template class Alignment<>;
