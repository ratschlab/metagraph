#ifndef __ALIGNER_HELPER_HPP__
#define __ALIGNER_HELPER_HPP__

#include <vector>
#include <numeric>
#include <string>
#include <cassert>
#include <array>
#include <ostream>
#include <unordered_map>
#include <json/json.h>
#include <memory>

#include "config.hpp"
#include "sequence_graph.hpp"
#include "reverse_complement.hpp"


class Cigar {
  public:
    enum Operator {
        CLIPPED,
        MATCH,
        MISMATCH,
        INSERTION,
        DELETION
    };

    typedef uint32_t LengthType;
    typedef std::array<Operator, 128> OperatorTableRow;
    typedef std::array<OperatorTableRow, 128> OperatorTable;

    static OperatorTable char_to_op;
    static const OperatorTableRow& get_op_row(char a) { return char_to_op[a]; }
    static void initialize_opt_table(const std::string &alphabet);

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

    void clear() { cigar_.clear(); }

    LengthType get_clipping() const {
        return cigar_.size() && cigar_.front().first == Operator::CLIPPED
            ? cigar_.front().second
            : 0;
    }

  private:
    std::vector<std::pair<Operator, LengthType>> cigar_;
};

typedef int32_t score_t;

class DBGAlignerConfig {
  public:
    typedef ::score_t score_t;
    typedef std::array<int8_t, 128> ScoreMatrixRow;
    typedef std::array<ScoreMatrixRow, 128> ScoreMatrix;

    DBGAlignerConfig(const ScoreMatrix &score_matrix,
                     int8_t gap_opening = -3,
                     int8_t gap_extension = -1);

    explicit DBGAlignerConfig(const Config &config, const DeBruijnGraph &graph);

    DBGAlignerConfig(ScoreMatrix&& score_matrix,
                     int8_t gap_opening = -3,
                     int8_t gap_extension = -1);

    template <class StringItA, class StringItB>
    score_t score_sequences(StringItA a_begin, StringItA a_end,
                            StringItB b_begin) const {
        assert(a_end >= a_begin);

        return std::inner_product(
            a_begin,
            a_end,
            b_begin,
            score_t(0),
            std::plus<score_t>(),
            [&](const char &a, const char &b) -> score_t {
                return score_matrix_[a][b];
            }
        );
    }

    template <class StringIt>
    score_t match_score(StringIt begin, StringIt end) const {
        return score_sequences(begin, end, begin);
    }

    const ScoreMatrixRow& get_row(char char_in_query) const {
        return score_matrix_[char_in_query];
    }

    size_t queue_size = std::numeric_limits<size_t>::max();
    size_t bandwidth = std::numeric_limits<size_t>::max();
    size_t num_alternative_paths = 1;
    size_t min_seed_length = 1;
    size_t max_seed_length = std::numeric_limits<size_t>::max();
    size_t max_num_seeds_per_locus = 1;
    score_t min_cell_score = 0;

    int8_t gap_opening_penalty;
    int8_t gap_extension_penalty;

    bool forward_and_reverse_complement = false;

    static ScoreMatrix scoring_matrix(const Config &config,
                                      const DeBruijnGraph &graph);

    // Protein matrices
    static const ScoreMatrix score_matrix_blosum62;

    static ScoreMatrix dna_scoring_matrix(int8_t match_score,
                                          int8_t mm_transition,
                                          int8_t mm_transversion);

    static ScoreMatrix unit_scoring_matrix(int8_t match_score,
                                           const std::string &alphabet);

  private:
    ScoreMatrix score_matrix_;
};



template <typename NodeType>
class Alignment {
  public:
    typedef NodeType node_index;
    typedef ::score_t score_t;

    // dynamic programming table stores score columns and steps needed to reconstruct paths
    struct Step {
        Cigar::Operator cigar_op;
        NodeType prev_node;
    };

    struct Column {
        std::vector<score_t> scores;
        std::vector<Step> steps;
        char last_char;
        size_t best_pos;

        const score_t& best_score() const { return scores.at(best_pos); }
        const Step& best_step() const { return steps.at(best_pos); }

        bool operator<(const Column &other) const {
            return best_score() < other.best_score();
        }
    };

    // storage for DP table
    typedef std::unordered_map<NodeType, Column> DPTable;

    // Used for constructing seeds
    Alignment(const char* query_begin = nullptr,
              const char* query_end = nullptr,
              std::vector<NodeType>&& nodes = {},
              score_t score = 0,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : Alignment(query_begin,
                      query_end,
                      std::move(nodes),
                      std::string(query_begin, query_end),
                      query_end - query_begin,
                      score,
                      Cigar(Cigar::Operator::MATCH, query_end - query_begin),
                      clipping,
                      orientation,
                      offset) {
        assert(nodes.empty() || clipping || is_exact_match());
    }

    Alignment(const char* query_begin,
              const char* query_end,
              std::vector<NodeType>&& nodes,
              std::string&& sequence,
              score_t score,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0);

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

    template <class StringIt>
    void set_query_begin(const StringIt *begin) {
        assert((query_end_ - query_begin_) + begin >= begin);
        query_end_ = (query_end_ - query_begin_) + begin;
        query_begin_ = begin;
    }

    const std::string& get_sequence() const { return sequence_; }

    const Cigar& get_cigar() const { return cigar_; }
    void set_cigar(Cigar&& cigar) { cigar_ = std::move(cigar); }

    bool get_orientation() const { return orientation_; }
    size_t get_offset() const { return offset_; }
    Cigar::LengthType get_clipping() const { return cigar_.get_clipping(); }

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

    Json::Value to_json(const std::string &query,
                        const DeBruijnGraph &graph,
                        bool is_secondary = false,
                        const std::string &name = "",
                        const std::string &label = "") const;

    std::shared_ptr<const std::string>
    load_from_json(const Json::Value &alignment,
                   const DeBruijnGraph &graph);

    bool is_valid(const DeBruijnGraph &graph) const;

  private:
    Alignment(const char* query_begin,
              const char* query_end,
              std::vector<NodeType>&& nodes = {},
              std::string&& sequence = "",
              size_t num_matches = 0,
              score_t score = 0,
              Cigar&& cigar = Cigar(),
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : query_begin_(query_begin),
            query_end_(query_end),
            nodes_(std::move(nodes)),
            sequence_(std::move(sequence)),
            num_matches_(num_matches),
            score_(score),
            cigar_(Cigar::Operator::CLIPPED, clipping),
            orientation_(orientation),
            offset_(offset) { cigar_.append(std::move(cigar)); }

    Json::Value path_json(size_t node_size, const std::string &label = "") const;

    const char* query_begin_;
    const char* query_end_;
    std::vector<NodeType> nodes_;
    std::string sequence_;
    uint64_t num_matches_;
    score_t score_;
    Cigar cigar_;
    bool orientation_;
    size_t offset_;
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


template <typename NodeType>
class QueryAlignment {
  public:
    QueryAlignment(const std::string &query)
          : query_(query),
            query_rc_(query) {
        reverse_complement(const_cast<char*>(&*query_rc_.begin()),
                           const_cast<char*>(&*query_rc_.end()));
    }

    template <typename... Args>
    QueryAlignment(std::string&& query, Args&&... args) noexcept
          : query_(std::move(query)),
            query_rc_(query_),
            alignments_(std::forward<Args>(args)...) {
        reverse_complement(const_cast<char*>(&*query_rc_.begin()),
                           const_cast<char*>(&*query_rc_.end()));
    }

    size_t size() const { return alignments_.size(); }
    bool empty() const { return alignments_.empty(); }

    template <typename... Args>
    void emplace_back(Args&&... args) {
        alignments_.emplace_back(std::forward<Args>(args)...);

        assert(alignments_.back().get_orientation()
            || alignments_.back().get_query_begin() >= &*query_.begin());
        assert(alignments_.back().get_orientation()
            || alignments_.back().get_query_end() <= &*query_.end());
        assert(!alignments_.back().get_orientation()
            || alignments_.back().get_query_begin() >= &*query_rc_.begin());
        assert(!alignments_.back().get_orientation()
            || alignments_.back().get_query_end() <= &*query_rc_.end());
    }

    void pop_back() { alignments_.pop_back(); }
    void clear() { alignments_.clear(); }

    const std::string& get_query() const { return query_; }
    const std::string& get_query_reverse_complement() const { return query_rc_; }
    const Alignment<NodeType>& front() const { return alignments_.front(); }
    const Alignment<NodeType>& back() const { return alignments_.back(); }
    const Alignment<NodeType>& operator[](size_t i) const { return alignments_[i]; }

    typedef typename std::vector<Alignment<NodeType>>::iterator iterator;
    typedef typename std::vector<Alignment<NodeType>>::const_iterator const_iterator;

    iterator begin() { return alignments_.begin(); }
    iterator end() { return alignments_.end(); }

    const_iterator begin() const { return alignments_.cbegin(); }
    const_iterator end() const { return alignments_.cend(); }

    const_iterator cbegin() const { return alignments_.cbegin(); }
    const_iterator cend() const { return alignments_.cend(); }

    template <class Iterator>
    void erase(Iterator begin, Iterator end) { alignments_.erase(begin, end); }

    bool operator==(const QueryAlignment &other) const {
        return query_ == other.query_
            && std::equal(alignments_.begin(), alignments_.end(),
                          other.alignments_.begin(), other.alignments_.end());
    }

    bool operator!=(const QueryAlignment &other) const { return !(*this == other); }

  private:
    const std::string query_;
    const std::string query_rc_;
    std::vector<Alignment<NodeType>> alignments_;
};


#endif  // __ALIGNER_HELPER_HPP__
