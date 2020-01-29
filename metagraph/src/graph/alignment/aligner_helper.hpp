#ifndef __ALIGNER_HELPER_HPP__
#define __ALIGNER_HELPER_HPP__

#include <array>
#include <memory>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>
#include <cassert>

#include <json/json.h>
#include <tsl/ordered_map.h>

#include "common/seq_tools/reverse_complement.hpp"
#include "graph/representation/base/sequence_graph.hpp"


class Cigar {
  public:
    enum Operator : int32_t {
        CLIPPED,
        MATCH,
        MISMATCH,
        INSERTION,
        DELETION
    };

    typedef uint32_t LengthType;
    typedef std::array<Operator, 128> OperatorTableRow;
    typedef std::array<OperatorTableRow, 128> OperatorTable;
    typedef std::pair<Operator, LengthType> value_type;

    static OperatorTable char_to_op;
    static const OperatorTableRow& get_op_row(char a) { return char_to_op[a]; }
    static void initialize_opt_table(const std::string &alphabet, const uint8_t *encoding);

    Cigar(Operator op = Operator::CLIPPED, LengthType num = 0)
          : cigar_(num ? 1 : 0, std::make_pair(op, num)) { }

    // See section 1.4 in https://samtools.github.io/hts-specs/SAMv1.pdf for
    // a specification of the CIGAR string format.
    // e.g., 3=1X2I3D for 3 matches, 1 mismatch, 2 insertions, 3 deletions
    explicit Cigar(const std::string &cigar_str);

    size_t size() const { return cigar_.size(); }
    bool empty() const { return cigar_.empty(); }

    std::string to_string() const;

    void append(Operator op, LengthType num = 1);
    void append(Cigar&& other);

    void pop_front() {
        assert(cigar_.size());
        cigar_.erase(cigar_.begin(), cigar_.begin() + 1);
    }

    typedef typename std::vector<value_type>::iterator iterator;
    typedef typename std::vector<value_type>::const_iterator const_iterator;

    // This is essentially just a vector, so there's no reason not to have it editable
    iterator begin() { return cigar_.begin(); }
    iterator end() { return cigar_.end(); }
    const_iterator begin() const { return cigar_.cbegin(); }
    const_iterator end() const { return cigar_.cend(); }

    value_type& front() { return cigar_.front(); }
    value_type& back() { return cigar_.back(); }
    const value_type& front() const { return cigar_.front(); }
    const value_type& back() const { return cigar_.back(); }

    bool operator==(const Cigar &other) const { return cigar_ == other.cigar_; }

    bool operator!=(const Cigar &other) const { return !(*this == other); }

    void clear() { cigar_.clear(); }

    LengthType get_clipping() const {
        return cigar_.size() && cigar_.front().first == Operator::CLIPPED
            ? cigar_.front().second
            : 0;
    }

    LengthType get_end_clipping() const {
        return cigar_.size() && cigar_.back().first == Operator::CLIPPED
            ? cigar_.back().second
            : 0;
    }

    // Return true if the cigar is valid. reference_begin points to the first
    // character of the reference sequence after clipping is trimmed
    bool is_valid(const char *reference_begin, const char *reference_end,
                  const char *query_begin, const char *query_end) const;

  private:
    std::vector<value_type> cigar_;
};

typedef int32_t score_t;

class DBGAlignerConfig {
  public:
    typedef ::score_t score_t;
    typedef std::array<int8_t, 128> ScoreMatrixRow;
    typedef std::array<ScoreMatrixRow, 128> ScoreMatrix;

    // Set parameters manually and call `set_scoring_matrix()`
    DBGAlignerConfig() {}

    explicit DBGAlignerConfig(const ScoreMatrix &score_matrix,
                              int8_t gap_opening = -3,
                              int8_t gap_extension = -1);

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
            [&](char a, char b) -> score_t { return score_matrix_[a][b]; }
        );
    }

    template <class StringIt>
    score_t match_score(StringIt begin, StringIt end) const {
        return score_sequences(begin, end, begin);
    }

    score_t score_cigar(const char *reference_begin, const char *reference_end,
                        const char *query_begin, const char *query_end,
                        const Cigar &cigar) const;

    const ScoreMatrixRow& get_row(char char_in_query) const {
        return score_matrix_[char_in_query];
    }

    size_t queue_size = std::numeric_limits<size_t>::max();
    size_t bandwidth = std::numeric_limits<size_t>::max();
    size_t num_alternative_paths = 1;
    size_t min_seed_length = 1;
    size_t max_seed_length = std::numeric_limits<size_t>::max();
    size_t max_num_seeds_per_locus = 1;
    // thresholds for scores
    score_t min_cell_score = 0;
    score_t min_path_score = 0;

    int8_t gap_opening_penalty;
    int8_t gap_extension_penalty;

    bool forward_and_reverse_complement = false;

    bool alignment_edit_distance;
    int8_t alignment_match_score;
    int8_t alignment_mm_transition_score;
    int8_t alignment_mm_transversion_score;

    bool check_config_scores() const;

    void set_scoring_matrix();

    // Protein matrices
    static const ScoreMatrix score_matrix_blosum62;

    static ScoreMatrix dna_scoring_matrix(int8_t match_score,
                                          int8_t mm_transition_score,
                                          int8_t mm_transversion_score);

    static ScoreMatrix unit_scoring_matrix(int8_t match_score,
                                           const std::string &alphabet,
                                           const uint8_t *encoding);

  private:
    ScoreMatrix score_matrix_;
};


template <typename NodeType = SequenceGraph::node_index>
class DPTable;

// Note: this object stores pointers to the query sequence, so it is the user's
//       responsibility to ensure that the query sequence is not destroyed when
//       calling this class' methods
template <typename NodeType = SequenceGraph::node_index>
class Alignment {
  public:
    typedef NodeType node_index;
    typedef ::score_t score_t;
    typedef ::DPTable<NodeType> DPTable;

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
    // Since insertions into DPTable may invalidate iterators, use pointers instead
    Alignment(const DPTable &dp_table,
              const typename DPTable::value_type *column,
              size_t start_pos,
              score_t score,
              const char* path_end,
              bool orientation,
              size_t offset);

    void append(Alignment&& other, size_t overlap = 0, int8_t match_score = 0);

    size_t size() const { return nodes_.size(); }
    const NodeType& front() const { return nodes_.front(); }
    const NodeType& back() const { return nodes_.back(); }
    const std::vector<NodeType>& get_nodes() const { return nodes_; }
    bool empty() const { return nodes_.empty(); }

    score_t get_score() const { return score_; }
    uint64_t get_num_matches() const { return num_matches_; }

    void recompute_score(const DBGAlignerConfig &config);

    size_t query_size() const { return query_end_ - query_begin_; }
    const char* get_query_begin() const { return query_begin_; }
    const char* get_query_end() const { return query_end_; }

    void set_query_begin(const char *begin) {
        assert((query_end_ - query_begin_) + begin >= begin);
        query_end_ = (query_end_ - query_begin_) + begin;
        query_begin_ = begin;
    }

    void extend_query_end(const char *end) {
        assert(end >= query_end_);
        if (end == query_end_)
            return;

        cigar_.append(Cigar::Operator::CLIPPED, end - query_end_);
    }

    const std::string& get_sequence() const { return sequence_; }

    const Cigar& get_cigar() const { return cigar_; }
    void set_cigar(Cigar&& cigar) { cigar_ = std::move(cigar); }

    bool get_orientation() const { return orientation_; }
    size_t get_offset() const { return offset_; }
    Cigar::LengthType get_clipping() const { return cigar_.get_clipping(); }
    Cigar::LengthType get_end_clipping() const { return cigar_.get_end_clipping(); }

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
            && cigar_.front().first == Cigar::Operator::MATCH
            && query_begin_ + cigar_.front().second == query_end_;
    }

    Json::Value to_json(const std::string &query,
                        const DeBruijnGraph &graph,
                        bool is_secondary = false,
                        const std::string &name = "",
                        const std::string &label = "") const;

    std::shared_ptr<const std::string>
    load_from_json(const Json::Value &alignment,
                   const DeBruijnGraph &graph);

    bool is_valid(const DeBruijnGraph &graph, const DBGAlignerConfig *config = nullptr) const;

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
        << alignment.get_cigar().to_string() << "\t"
        << alignment.get_offset();

    return out;
}


template <typename NodeType = SequenceGraph::node_index>
class QueryAlignment {
  public:
    typedef Alignment<NodeType> value_type;

    explicit QueryAlignment(const std::string &query);

    QueryAlignment(const QueryAlignment &other);
    QueryAlignment(QueryAlignment&& other) noexcept;

    size_t size() const { return alignments_.size(); }
    bool empty() const { return alignments_.empty(); }

    template <typename... Args>
    void emplace_back(Args&&... args) {
        alignments_.emplace_back(std::forward<Args>(args)...);

        // sanity checks
        assert(alignments_.back().get_orientation()
            || alignments_.back().get_query_begin() >= query_.c_str());
        assert(alignments_.back().get_orientation()
            || alignments_.back().get_query_end() <= query_.c_str() + query_.size());
        assert(!alignments_.back().get_orientation()
            || alignments_.back().get_query_begin() >= query_rc_.c_str());
        assert(!alignments_.back().get_orientation()
            || alignments_.back().get_query_end() <= query_rc_.c_str() + query_rc_.size());
    }

    void push_back(const value_type &alignment) { emplace_back(alignment); }
    void push_back(value_type&& alignment) { emplace_back(std::move(alignment)); }

    void pop_back() { alignments_.pop_back(); }
    void clear() { alignments_.clear(); }

    const std::string& get_query() const { return query_; }
    const std::string& get_query_reverse_complement() const { return query_rc_; }
    const value_type& front() const { return alignments_.front(); }
    const value_type& back() const { return alignments_.back(); }
    const value_type& operator[](size_t i) const { return alignments_[i]; }

    typedef typename std::vector<value_type>::const_iterator const_iterator;

    const_iterator begin() const { return alignments_.cbegin(); }
    const_iterator end() const { return alignments_.cend(); }
    const_iterator cbegin() const { return alignments_.cbegin(); }
    const_iterator cend() const { return alignments_.cend(); }

    template <class Iterator>
    void erase(Iterator begin, Iterator end) { alignments_.erase(begin, end); }

    bool operator==(const QueryAlignment &other) const {
        return query_ == other.query_ && alignments_ == other.alignments_;
    }

    bool operator!=(const QueryAlignment &other) const { return !(*this == other); }

    std::vector<double> get_alignment_weights(const DBGAlignerConfig &config) const;

  private:
    // When a QueryAlignment is copied or moved, the pointers in the alignment
    // vector may be incorrect, so this corrects them
    void fix_pointers(const std::string &query, const std::string &query_rc);

    const std::string query_;
    const std::string query_rc_;
    std::vector<value_type> alignments_;
};


// dynamic programming table stores score columns and steps needed to reconstruct paths
template <typename NodeType>
class DPTable {
  public:
    typedef ::score_t score_t;

    struct Column {
        Column(size_t size,
               score_t min_score,
               char start_char,
               std::vector<NodeType>&& in_nodes,
               size_t pos = 0)
              : scores(size, min_score),
                ops(size),
                prev_nodes(size),
                last_char(start_char),
                incoming(std::move(in_nodes)),
                best_pos(pos) {}

        std::vector<score_t> scores;
        std::vector<Cigar::Operator> ops;
        std::vector<NodeType> prev_nodes;
        char last_char;
        const std::vector<NodeType> incoming;
        size_t best_pos;

        const score_t& best_score() const { return scores.at(best_pos); }
        const Cigar::Operator& best_op() const { return ops.at(best_pos); }
        const NodeType& best_prev_node() const { return prev_nodes.at(best_pos); }

        bool operator<(const Column &other) const {
            return best_score() < other.best_score();
        }
    };

    explicit DPTable(const SequenceGraph &graph,
                     NodeType start_node,
                     char start_char,
                     score_t initial_score,
                     score_t min_score,
                     size_t size,
                     int8_t gap_opening_penalty,
                     int8_t gap_extension_penalty);

    typedef tsl::ordered_map<NodeType, Column> Storage;

    typedef NodeType key_type;
    typedef Column mapped_type;
    typedef typename Storage::value_type value_type;

    typedef typename Storage::iterator iterator;
    typedef typename Storage::const_iterator const_iterator;

    iterator begin() { return dp_table_.begin(); }
    iterator end() { return dp_table_.end(); }
    const_iterator begin() const { return dp_table_.begin(); }
    const_iterator end() const { return dp_table_.end(); }

    iterator find(const NodeType &node) { return dp_table_.find(node); }
    const_iterator find(const NodeType &node) const { return dp_table_.find(node); }

    size_t size() const { return dp_table_.size(); }

    template <typename... Args>
    std::pair<iterator, bool> emplace(Args&&... args) {
        return dp_table_.emplace(std::forward<Args>(args)...);
    }

    std::vector<Alignment<NodeType>>
    extract_alignments(const DeBruijnGraph &graph,
                       const DBGAlignerConfig &config,
                       score_t start_score,
                       const char *align_start,
                       bool orientation,
                       score_t min_path_score);

    const Storage& data() const { return dp_table_; }

  private:
    Storage dp_table_;
};



#endif  // __ALIGNER_HELPER_HPP__
