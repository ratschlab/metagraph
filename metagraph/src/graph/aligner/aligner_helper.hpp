#ifndef __ALIGNER_HELPER_HPP__
#define __ALIGNER_HELPER_HPP__

#include <vector>
#include <numeric>
#include <string>
#include <cassert>
#include <array>
#include <ostream>
#include <unordered_map>

#include "config.hpp"
#include "sequence_graph.hpp"


class Cigar {
  public:
    enum Operator {
        CLIPPED,
        MATCH,
        MISMATCH_TRANSITION,
        MISMATCH_TRANSVERSION,
        INSERTION,
        DELETION
    };

    static const std::array<std::array<Operator, 128>, 128> char_to_op;
    typedef uint32_t LengthType;

    static const std::array<Operator, 128>& get_op_row(char a) { return char_to_op[a]; }

    Cigar(Operator op = Operator::CLIPPED, LengthType num = 0)
          : cigar_(num ? 1 : 0, std::make_pair(op, num)) { }

    size_t size() const { return cigar_.size(); }
    bool empty() const { return cigar_.empty(); }

    std::string to_string() const;

    void append(Operator op, LengthType num = 1);
    void append(Cigar&& other);

    void pop_front() {
        assert(cigar_.size());
        cigar_.erase(cigar_.begin(), cigar_.begin() + 1);
    }

    typedef typename std::vector<std::pair<Operator, LengthType>>::iterator iterator;
    typedef typename std::vector<std::pair<Operator, LengthType>>::const_iterator const_iterator;

    // This is essentially just a vector, so there's no reason not to have it editable
    iterator begin() { return cigar_.begin(); }
    iterator end() { return cigar_.end(); }
    const_iterator begin() const { return cigar_.cbegin(); }
    const_iterator end() const { return cigar_.cend(); }

    bool operator==(const Cigar &other) const {
        return std::equal(begin(), end(), other.begin(), other.end());
    }

    bool operator!=(const Cigar &other) const { return !(*this == other); }

  private:
    std::vector<std::pair<Operator, LengthType>> cigar_;

    static std::array<std::array<Operator, 128>, 128> initialize_opt_table();
};

typedef int32_t score_t;

class DBGAlignerConfig {
  public:
    typedef ::score_t score_t;

    DBGAlignerConfig() {
        set_match_score(match_score_);
        set_mismatch_transition_score(mm_transition_);
        set_mismatch_transversion_score(mm_transversion_);
    }

    explicit DBGAlignerConfig(const Config &config)
          : queue_size(config.alignment_queue_size),
            num_alternative_paths(config.alignment_num_alternative_paths),
            min_seed_length(config.alignment_min_seed_length),
            max_seed_length(config.alignment_max_seed_length),
            max_num_seeds_per_locus(config.alignment_max_num_seeds_per_locus),
            gap_opening_penalty(-config.alignment_gap_opening_penalty),
            gap_extension_penalty(-config.alignment_gap_extension_penalty),
            min_cell_score(config.alignment_min_cell_score),
            mm_transition_(-config.alignment_mm_transition),
            match_score_(config.alignment_match_score) {
        set_mismatch_transversion_score(-config.alignment_mm_transversion);
    }

    explicit DBGAlignerConfig(const std::array<std::array<int8_t, 128>, 128> &score_matrix,
                              size_t gap_opening,
                              size_t gap_extension)
          : gap_opening_penalty(gap_opening),
            gap_extension_penalty(gap_extension),
            score_matrix_(score_matrix) { }

    DBGAlignerConfig(std::array<std::array<int8_t, 128>, 128>&& score_matrix,
                     size_t gap_opening,
                     size_t gap_extension)
          : gap_opening_penalty(gap_opening),
            gap_extension_penalty(gap_extension),
            score_matrix_(std::move(score_matrix)) { }

    template <class StringIt>
    score_t score_sequences(StringIt a_begin, StringIt a_end,
                            StringIt b_begin, StringIt b_end) const {
        assert(a_end >= a_begin);
        assert(b_end >= b_begin);
        assert(a_end - a_begin == b_end - b_begin);

        std::ignore = b_end;

        return std::inner_product(
            a_begin, a_end,
            b_begin,
            score_t(0),
            std::plus<score_t>(),
            [&](const char &a, const char &b) -> score_t { return score_matrix_[a][b];}
        );
    }

    template <class StringIt>
    score_t match_score(StringIt begin, StringIt end) const {
        return score_sequences(begin, end, begin, end);
    }

    void set_match_score(int8_t match_score);

    void set_mismatch_transition_score(int8_t mm_transition);

    void set_mismatch_transversion_score(int8_t mm_transversion);

    int8_t get_match_score() const { return match_score_; }
    int8_t get_mismatch_transition_score() const { return mm_transition_; }
    int8_t get_mismatch_transversion_score() const { return mm_transversion_; }

    const std::array<int8_t, 128>& get_row(char char_in_query) const {
        return score_matrix_[char_in_query];
    }

    size_t queue_size = std::numeric_limits<size_t>::max();
    size_t num_alternative_paths = 1;
    size_t min_seed_length = 1;
    size_t max_seed_length = std::numeric_limits<size_t>::max();
    size_t max_num_seeds_per_locus = 1;

    // affine costs
    int8_t gap_opening_penalty = -3;
    int8_t gap_extension_penalty = -1;
    // unit costs
    // int8_t gap_opening_penalty = -1;
    // int8_t gap_extension_penalty = -1;

    score_t min_cell_score = 0;

  private:
    // affine costs
    int8_t mm_transition_ = -1;
    int8_t mm_transversion_ = -2;
    int8_t match_score_ = 2;
    // unit costs
    // int8_t mm_transition_ = -1;
    // int8_t mm_transversion_ = -1;
    // int8_t match_score_ = 1;

    std::array<std::array<int8_t, 128>, 128> score_matrix_;
};



template <typename NodeType>
class Alignment {
  public:
    typedef NodeType node_index;
    typedef ::score_t score_t;
    // dynamic programming table stores score columns and steps needed to reconstruct paths
    // operation, previous node
    typedef std::pair<Cigar::Operator, NodeType> Step;

    // score column, traversal column, last char of node, position of best score
    typedef std::tuple<std::vector<score_t>,
                       std::vector<Step>,
                       char,
                       size_t>
            Column;

    // storage for DP table
    typedef std::unordered_map<NodeType, Column> DPTable;

    // Used for constructing seeds
    Alignment(const char* query_begin,
              const char* query_end,
              std::vector<NodeType>&& nodes,
              score_t score,
              size_t clipping = 0,
              bool orientation = false)
          : Alignment(query_begin,
                      query_end,
                      std::move(nodes),
                      std::string(query_begin, query_end),
                      query_end - query_begin,
                      score,
                      Cigar(Cigar::Operator::MATCH, query_end - query_begin),
                      clipping,
                      orientation) {
        assert(clipping || is_exact_match());
    }

    // TODO: construct multiple alignments from the same starting point
    Alignment(const DPTable &dp_table,
              const typename DPTable::value_type *column,
              size_t start_pos,
              score_t score,
              const char* path_end,
              bool orientation);

    void append(Alignment&& other, size_t overlap = 0, int8_t match_score = 0);

    size_t size() const { return nodes_.size(); }
    const NodeType& back() const { return nodes_.back(); }
    const NodeType& front() const { return nodes_.front(); }
    const std::vector<NodeType>& get_nodes() const { return nodes_; }
    bool empty() const { return nodes_.empty(); }

    score_t get_score() const { return score_; }
    uint64_t get_num_matches() const { return num_matches_; }

    size_t query_size() const { return query_end_ - query_begin_; }
    const char* get_query_end() const { return query_end_; }
    const char* get_query_begin() const { return query_begin_; }

    const std::string& get_sequence() const { return sequence_; }

    const Cigar& get_cigar() const { return cigar_; }
    void set_cigar(Cigar&& cigar) { cigar_ = std::move(cigar); }

    bool get_orientation() const { return orientation_; }

    bool operator<(const Alignment &other) const { return score_ < other.score_; }
    bool operator>(const Alignment &other) const { return score_ > other.score_; }

    typedef typename std::vector<NodeType>::iterator iterator;
    typedef typename std::vector<NodeType>::const_iterator const_iterator;

    const_iterator begin() const { return nodes_.cbegin(); }
    const_iterator end() const { return nodes_.cend(); }

    const NodeType& operator[](size_t i) const { return nodes_[i]; }
    const NodeType& at(size_t i) const { return nodes_.at(i); }

    bool operator==(const Alignment &other) const {
        return orientation_ == other.orientation_
            && score_ == other.score_
            && num_matches_ == other.num_matches_
            && sequence_ == other.sequence_
            && std::equal(query_begin_, query_end_, other.query_begin_, other.query_end_)
            && cigar_ == other.cigar_;
    }

    bool operator!=(const Alignment &other) const { return !(*this == other); }

    bool is_exact_match() const {
        return cigar_.size() == 1
            && cigar_.begin()->first == Cigar::Operator::MATCH
            && query_begin_ + cigar_.begin()->second == query_end_;
    }

  private:
    Alignment(const char* query_begin,
              const char* query_end,
              std::vector<NodeType>&& nodes = {},
              std::string&& sequence = "",
              size_t num_matches = 0,
              score_t score = 0,
              Cigar&& cigar = Cigar(),
              size_t clipping = 0,
              bool orientation = false)
          : query_begin_(query_begin),
            query_end_(query_end),
            nodes_(std::move(nodes)),
            sequence_(std::move(sequence)),
            num_matches_(num_matches),
            score_(score),
            cigar_(Cigar::Operator::CLIPPED, clipping),
            orientation_(orientation) { cigar_.append(std::move(cigar)); }

    const char* query_begin_;
    const char* query_end_;
    std::vector<NodeType> nodes_;
    std::string sequence_;
    uint64_t num_matches_;
    score_t score_;
    Cigar cigar_;
    bool orientation_;
};

template <typename NodeType>
std::ostream& operator<<(std::ostream& out, const Alignment<NodeType> &alignment) {
    out << (alignment.get_orientation() ? "-" : "+") << "\t"
        << alignment.get_sequence() << "\t"
        << alignment.get_score() << "\t"
        << alignment.get_num_matches() << "\t"
        << alignment.get_cigar().to_string();

    return out;
}


#endif  // __ALIGNER_HELPER_HPP__
