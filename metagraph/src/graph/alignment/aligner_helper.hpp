#ifndef __ALIGNER_HELPER_HPP__
#define __ALIGNER_HELPER_HPP__

#include <cassert>
#include <array>
#include <memory>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

#include <json/json.h>
#include <tsl/hopscotch_map.h>

#include "common/aligned_vector.hpp"
#include "common/seq_tools/reverse_complement.hpp"
#include "common/utils/template_utils.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {
namespace align {

class Cigar {
  public:
    enum Operator : int8_t {
        CLIPPED,
        MISMATCH,
        MATCH,
        DELETION,
        INSERTION
    };

    typedef uint32_t LengthType;
    typedef std::array<Operator, 128> OperatorTableRow;
    typedef std::array<OperatorTableRow, 128> OperatorTable;
    typedef std::pair<Operator, LengthType> value_type;

    static OperatorTable char_to_op;
    static const OperatorTableRow& get_op_row(char a) { return char_to_op[a]; }

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

    void pop_back() {
        assert(cigar_.size());
        cigar_.pop_back();
    }

    typedef typename std::vector<value_type>::iterator iterator;
    typedef typename std::vector<value_type>::const_iterator const_iterator;

    // This is essentially just a vector, so there's no reason not to have it editable
    iterator begin() { return cigar_.begin(); }
    iterator end() { return cigar_.end(); }
    const_iterator begin() const { return cigar_.cbegin(); }
    const_iterator end() const { return cigar_.cend(); }

    template <typename... Args>
    void insert(iterator it, Args&&... args) {
        cigar_.insert(it, std::forward<Args>(args)...);
    }

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

    size_t get_num_matches() const {
        return std::accumulate(begin(), end(), 0, [&](size_t old, const value_type &op) {
            return old + (op.first == Operator::MATCH) * op.second;
        });
    }

    // Return true if the cigar is valid. reference_begin points to the first
    // character of the reference sequence after clipping is trimmed
    bool is_valid(const std::string_view reference, const std::string_view query) const;

    static char opt_to_char(Cigar::Operator op);

  private:
    std::vector<value_type> cigar_;

    static OperatorTable initialize_opt_table();
};


// TODO: move to a separate header file aligner_config.hpp
class DBGAlignerConfig {
  public:
    typedef int32_t score_t;
    typedef std::array<int8_t, 128> ScoreMatrixRow;
    typedef std::array<ScoreMatrixRow, 128> ScoreMatrix;

    // Set parameters manually and call `set_scoring_matrix()`
    DBGAlignerConfig() {}

    explicit DBGAlignerConfig(const ScoreMatrix &score_matrix,
                              int8_t gap_opening = -5,
                              int8_t gap_extension = -2);

    DBGAlignerConfig(ScoreMatrix&& score_matrix,
                     int8_t gap_opening = -5,
                     int8_t gap_extension = -2);

    score_t score_sequences(const std::string_view a, const std::string_view b) const {
        return std::inner_product(
            a.begin(), a.end(), b.begin(), score_t(0), std::plus<score_t>(),
            [&](char a, char b) -> score_t { return score_matrix_[a][b]; }
        );
    }

    score_t match_score(const std::string_view query) const {
        return score_sequences(query, query);
    }

    score_t score_cigar(const std::string_view reference,
                        const std::string_view query,
                        const Cigar &cigar) const;

    const ScoreMatrixRow& get_row(char char_in_query) const {
        return score_matrix_[char_in_query];
    }

    size_t queue_size = std::numeric_limits<size_t>::max();
    size_t bandwidth = std::numeric_limits<size_t>::max();
    size_t num_alternative_paths = 1;
    size_t min_seed_length = 1;
    size_t max_seed_length = std::numeric_limits<size_t>::max();
    // thresholds for scores
    score_t min_cell_score = 0;
    score_t min_path_score = 0;
    score_t xdrop = std::numeric_limits<score_t>::max();

    double exact_kmer_match_fraction = 0.0;
    double max_nodes_per_seq_char = std::numeric_limits<double>::max();
    double max_ram_per_alignment = std::numeric_limits<double>::max();

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


template <typename NodeType>
class DPTable;

// Note: this object stores pointers to the query sequence, so it is the user's
//       responsibility to ensure that the query sequence is not destroyed when
//       calling this class' methods
template <typename NodeType = SequenceGraph::node_index>
class Alignment {
  public:
    typedef NodeType node_index;
    typedef DBGAlignerConfig::score_t score_t;

    // Used for constructing seeds
    Alignment(const std::string_view query = {},
              std::vector<NodeType>&& nodes = {},
              score_t score = 0,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : Alignment(query,
                      std::move(nodes),
                      std::string(query),
                      score,
                      Cigar(Cigar::MATCH, query.size()),
                      clipping,
                      orientation,
                      offset) {
        assert(nodes.empty() || clipping || is_exact_match());
    }

    Alignment(const std::string_view query,
              std::vector<NodeType>&& nodes,
              std::string&& sequence,
              score_t score,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0);

    // TODO: construct multiple alignments from the same starting point
    Alignment(const DPTable<NodeType> &dp_table,
              const DBGAlignerConfig &config,
              const std::string_view query_view,
              typename DPTable<NodeType>::const_iterator column,
              size_t start_pos,
              size_t offset,
              NodeType *start_node,
              const Alignment &seed);

    void append(Alignment&& other);

    size_t size() const { return nodes_.size(); }
    const NodeType& front() const { return nodes_.front(); }
    const NodeType& back() const { return nodes_.back(); }
    const std::vector<NodeType>& get_nodes() const { return nodes_; }
    bool empty() const { return nodes_.empty(); }

    score_t get_score() const { return score_; }
    uint64_t get_num_matches() const { return cigar_.get_num_matches(); }

    const std::string_view get_query() const {
        return std::string_view(query_begin_, query_end_ - query_begin_);
    }
    const char* get_query_end() const { return query_end_; }

    void set_query_begin(const char *begin) {
        assert((query_end_ - query_begin_) + begin >= begin);
        query_end_ = (query_end_ - query_begin_) + begin;
        query_begin_ = begin;
    }

    void extend_query_begin(const char *begin) {
        size_t clipping = get_clipping();
        assert(begin <= query_begin_ - clipping);
        if (begin == query_begin_ - clipping)
            return;

        if (clipping) {
            cigar_.front().second += query_begin_ - clipping - begin;
            return;
        }

        cigar_.insert(
            cigar_.begin(),
            typename Cigar::value_type { Cigar::CLIPPED, query_begin_ - begin }
        );
    }

    void extend_query_end(const char *end) {
        size_t end_clipping = get_end_clipping();
        assert(end >= query_end_ + end_clipping);
        if (end == query_end_ + end_clipping)
            return;

        cigar_.append(Cigar::CLIPPED, end - query_end_ - end_clipping);
    }

    void trim_clipping() {
        if (get_clipping())
            cigar_.pop_front();
    }

    void trim_end_clipping() {
        if (get_end_clipping())
            cigar_.pop_back();
    }

    void trim_offset(const DeBruijnGraph *graph = nullptr);

    void reverse_complement(const DeBruijnGraph &graph,
                            const std::string_view query_rev_comp);

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
            && sequence_ == other.sequence_
            && std::equal(query_begin_, query_end_, other.query_begin_, other.query_end_)
            && cigar_ == other.cigar_;
    }

    bool operator!=(const Alignment &other) const { return !(*this == other); }

    bool is_exact_match() const {
        return cigar_.size() == 1
            && cigar_.front().first == Cigar::MATCH
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
    Alignment(const std::string_view query,
              std::vector<NodeType>&& nodes = {},
              std::string&& sequence = "",
              score_t score = 0,
              Cigar&& cigar = Cigar(),
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : query_begin_(query.data()),
            query_end_(query.data() + query.size()),
            nodes_(std::move(nodes)),
            sequence_(std::move(sequence)),
            score_(score),
            cigar_(Cigar::CLIPPED, clipping),
            orientation_(orientation),
            offset_(offset) { cigar_.append(std::move(cigar)); }

    Json::Value path_json(size_t node_size, const std::string &label = "") const;

    const char* query_begin_;
    const char* query_end_;
    std::vector<NodeType> nodes_;
    std::string sequence_;
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

    QueryAlignment(const std::string_view query);
    QueryAlignment(const QueryAlignment &other);
    QueryAlignment(QueryAlignment&& other) noexcept;

    size_t size() const { return alignments_.size(); }
    bool empty() const { return alignments_.empty(); }

    template <typename... Args>
    void emplace_back(Args&&... args) {
        alignments_.emplace_back(std::forward<Args>(args)...);

        // sanity checks
        assert(alignments_.back().get_orientation()
            || alignments_.back().get_query().data() >= query_.c_str());
        assert(alignments_.back().get_orientation()
            || alignments_.back().get_query_end() <= query_.c_str() + query_.size());
        assert(!alignments_.back().get_orientation()
            || alignments_.back().get_query().data() >= query_rc_.c_str());
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

  private:
    // When a QueryAlignment is copied or moved, the pointers in the alignment
    // vector may be incorrect, so this corrects them
    void fix_pointers(const std::string &query, const std::string &query_rc);

    const std::string query_;
    const std::string query_rc_;
    std::vector<value_type> alignments_;
};


// dynamic programming table stores score columns and steps needed to reconstruct paths
template <typename NodeType = SequenceGraph::node_index>
class DPTable {
  public:
    typedef DBGAlignerConfig::score_t score_t;

    struct Column {
        Column() = default;

        Column(size_t size,
               score_t min_score,
               char start_char,
               size_t pos = 0,
               size_t priority_pos = 0,
               size_t start = 0,
               size_t end = std::numeric_limits<size_t>::max())
              : size_(size),
                min_score_(min_score),
                scores(std::min(end, size) - start + 8, min_score),
                gap_scores(scores.size(), min_score),
                ops(scores.size(), Cigar::CLIPPED),
                prev_nodes(scores.size(), 0),
                gap_prev_nodes(scores.size(), 0),
                gap_count(scores.size(), 0),
                last_char(start_char),
                best_pos(std::min(std::max(pos, start), start + scores.size() - (size_t)9)),
                last_priority_pos(std::min(std::max(priority_pos, start), start + scores.size() - (size_t)9)),
                start_index(start) {}

        size_t size_;
        score_t min_score_;
        AlignedVector<score_t> scores;
        AlignedVector<score_t> gap_scores;
        AlignedVector<Cigar::Operator> ops;
        AlignedVector<uint8_t> prev_nodes;
        AlignedVector<uint8_t> gap_prev_nodes;
        AlignedVector<int32_t> gap_count;
        mutable std::vector<NodeType> incoming_nodes;
        char last_char;
        size_t best_pos;
        size_t last_priority_pos;
        size_t start_index;

        void expand_to_cover(size_t begin, size_t end) {
            assert(best_pos >= start_index);
            assert(best_pos - start_index < scores.size());
            assert(last_priority_pos >= start_index);
            assert(last_priority_pos - start_index < scores.size());

            if (begin >= start_index) {
                // the current range already covers [begin, end)
                if (end <= start_index + scores.size() - 8)
                    return;

                // extend the range to the right to reach end
                scores.resize(end + 8 - start_index, min_score_);
                gap_scores.resize(end + 8 - start_index, min_score_);
                ops.resize(end + 8 - start_index, Cigar::CLIPPED);
                prev_nodes.resize(end + 8 - start_index, 0);
                gap_prev_nodes.resize(end + 8 - start_index, 0);
                gap_count.resize(end + 8 - start_index, 0);
            } else if (end <= start_index + scores.size() - 8) {
                // extend the range to the left to reach begin
                size_t shift = start_index - begin;
                start_index = begin;
                scores.insert(scores.begin(), shift, min_score_);
                gap_scores.insert(gap_scores.begin(), shift, min_score_);
                ops.insert(ops.begin(), shift, Cigar::CLIPPED);
                prev_nodes.insert(prev_nodes.begin(), shift, 0);
                gap_prev_nodes.insert(gap_prev_nodes.begin(), shift, 0);
                gap_count.insert(gap_count.begin(), shift, 0);
            } else {
                // extend the range in both directions
                size_t shift = start_index - begin;
                start_index = begin;
                end += 8;

                size_t new_size = end - begin;

                scores.reserve(new_size);
                gap_scores.reserve(new_size);
                ops.reserve(new_size);
                prev_nodes.reserve(new_size);
                gap_prev_nodes.reserve(new_size);
                gap_count.reserve(new_size);

                scores.insert(scores.begin(), shift, min_score_);
                gap_scores.insert(gap_scores.begin(), shift, min_score_);
                ops.insert(ops.begin(), shift, Cigar::CLIPPED);
                prev_nodes.insert(prev_nodes.begin(), shift, 0);
                gap_prev_nodes.insert(gap_prev_nodes.begin(), shift, 0);
                gap_count.insert(gap_count.begin(), shift, 0);

                scores.resize(new_size, min_score_);
                gap_scores.resize(new_size, min_score_);
                ops.resize(new_size, Cigar::CLIPPED);
                prev_nodes.resize(new_size, 0);
                gap_prev_nodes.resize(new_size, 0);
                gap_count.resize(new_size, 0);
            }

            assert(best_pos >= start_index);
            assert(best_pos - start_index < scores.size());
            assert(last_priority_pos >= start_index);
            assert(last_priority_pos - start_index < scores.size());
        }

        const score_t& best_score() const { return scores.at(best_pos - start_index); }
        const score_t& last_priority_value() const { return scores.at(last_priority_pos - start_index); }
        const Cigar::Operator& best_op() const { return ops.at(best_pos - start_index); }

        bool operator<(const Column &other) const {
            return best_score() < other.best_score();
        }

        size_t cell_size() const {
            return sizeof(score_t) * 2 + sizeof(Cigar::Operator)
                 + sizeof(uint8_t) * 2 + sizeof(int32_t);
        }

        size_t bytes_taken() const {
            return sizeof(Column)
                + sizeof(score_t) * scores.capacity()
                + sizeof(score_t) * gap_scores.capacity()
                + sizeof(Cigar::Operator) * ops.capacity()
                + sizeof(uint8_t) * prev_nodes.capacity()
                + sizeof(uint8_t) * gap_prev_nodes.capacity()
                + sizeof(int32_t) * gap_count.capacity()
                + sizeof(NodeType) * incoming_nodes.capacity();
        }

        size_t size() const { return size_; }

        uint8_t rank_prev_node(NodeType node) const {
            assert(incoming_nodes.size() < std::numeric_limits<uint8_t>::max());
            for (uint8_t i = 0; i < incoming_nodes.size(); ++i) {
                if (incoming_nodes[i] == node)
                    return i + 1;
            }

            incoming_nodes.push_back(node);
            return incoming_nodes.size();
        }

        NodeType select_prev_node(uint8_t rank) const {
            assert(rank);
            return incoming_nodes.at(rank - 1);
        }
    };

    DPTable() {}

    bool add_seed(const Alignment<NodeType> &seed,
                  const DBGAlignerConfig &config,
                  size_t size,
                  size_t start_pos,
                  size_t query_offset = 0);

    typedef tsl::hopscotch_map<NodeType, Column> Storage;

    typedef NodeType key_type;
    typedef Column mapped_type;
    typedef typename Storage::value_type value_type;

    typedef typename Storage::iterator iterator;
    typedef typename Storage::const_iterator const_iterator;

    void expand_to_cover(iterator it, size_t begin, size_t end) {
        size_t old_size = it->second.bytes_taken();
        it.value().expand_to_cover(begin, end);
        num_bytes_ += it->second.bytes_taken() - old_size;
    }

    iterator begin() { return dp_table_.begin(); }
    iterator end() { return dp_table_.end(); }
    const_iterator begin() const { return dp_table_.begin(); }
    const_iterator end() const { return dp_table_.end(); }

    iterator find(const NodeType &node) { return dp_table_.find(node); }
    const_iterator find(const NodeType &node) const { return dp_table_.find(node); }

    size_t size() const { return dp_table_.size(); }

    size_t num_bytes() const { return num_bytes_; }

    void clear() {
        dp_table_.clear();
        num_bytes_ = 0;
    }

    template <typename... Args>
    std::pair<iterator, bool> emplace(Args&&... args) {
        auto pair = dp_table_.emplace(std::forward<Args>(args)...);

        if (pair.second)
            num_bytes_ += pair.first.value().bytes_taken();

        return pair;
    }

    void erase(NodeType key) { dp_table_.erase(key); }
    size_t count(NodeType key) const { return dp_table_.count(key); }

    void extract_alignments(const DeBruijnGraph &graph,
                            const DBGAlignerConfig &config,
                            const std::string_view query_view,
                            std::function<void(Alignment<NodeType>&&, NodeType)> callback,
                            score_t min_path_score,
                            const Alignment<NodeType> &seed,
                            NodeType *node = nullptr);

    std::pair<NodeType, score_t> best_score() const {
        auto mx = std::max_element(begin(), end(), utils::LessSecond());
        return std::make_pair(mx->first, mx->second.best_score());
    }

    const Storage& data() const { return dp_table_; }
    size_t get_query_offset() const { return query_offset_; }
    NodeType get_start_node() const { return start_node_; }

  private:
    Storage dp_table_;
    NodeType start_node_;
    size_t query_offset_ = 0;
    size_t num_bytes_ = 0;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif  // __ALIGNER_HELPER_HPP__
