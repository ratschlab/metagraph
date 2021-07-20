#ifndef __ALIGNER_ALIGNMENT_HPP__
#define __ALIGNER_ALIGNMENT_HPP__

#include <cassert>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <json/json.h>

#include "aligner_cigar.hpp"
#include "aligner_config.hpp"
#include "graph/representation/base/sequence_graph.hpp"


namespace mtg {
namespace graph {
namespace align {

// Note: this object stores pointers to the query sequence, so it is the user's
//       responsibility to ensure that the query sequence is not destroyed when
//       calling this class' methods
class Alignment {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef DBGAlignerConfig::score_t score_t;

    Alignment() {}

    Alignment(std::string_view query,
              std::vector<node_index>&& nodes,
              std::string&& sequence,
              score_t score,
              Cigar&& cigar,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : query_(query), nodes_(std::move(nodes)), sequence_(std::move(sequence)),
            score_(score), cigar_(Cigar::CLIPPED, clipping), orientation_(orientation),
            offset_(offset) { cigar_.append(std::move(cigar)); }

    // Used for constructing gapless Alignments
    Alignment(std::string_view query,
              std::vector<node_index>&& nodes,
              std::string&& sequence,
              score_t score,
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0);

    void append(Alignment&& other);

    size_t size() const { return nodes_.size(); }
    bool empty() const { return nodes_.empty(); }
    const std::vector<node_index>& get_nodes() const { return nodes_; }

    score_t get_score() const { return score_; }

    std::string_view get_query() const { return query_; }

    void extend_query_begin(const char *begin) {
        size_t clipping = get_clipping();
        const char *full_query_begin = query_.data() - clipping;
        assert(begin <= full_query_begin);
        if (begin == full_query_begin)
            return;

        if (clipping) {
            cigar_.front().second += full_query_begin - begin;
        } else {
            cigar_.insert(cigar_.begin(),
                          std::make_pair(Cigar::CLIPPED, full_query_begin - begin));
        }
    }

    void extend_query_end(const char *end) {
        const char *full_query_end = query_.data() + query_.size() + get_end_clipping();
        assert(end >= full_query_end);
        if (end > full_query_end)
            cigar_.append(Cigar::CLIPPED, end - full_query_end);
    }

    void trim_clipping() {
        if (get_clipping())
            cigar_.pop_front();
    }

    void trim_end_clipping() {
        if (get_end_clipping())
            cigar_.pop_back();
    }

    size_t trim_offset();

    void reverse_complement(const DeBruijnGraph &graph, std::string_view query_rev_comp);

    const std::string& get_sequence() const { return sequence_; }
    const Cigar& get_cigar() const { return cigar_; }
    bool get_orientation() const { return orientation_; }
    size_t get_offset() const { return offset_; }
    Cigar::LengthType get_clipping() const { return cigar_.get_clipping(); }
    Cigar::LengthType get_end_clipping() const { return cigar_.get_end_clipping(); }

    bool operator==(const Alignment &other) const {
        return orientation_ == other.orientation_
            && offset_ == other.offset_
            && score_ == other.score_
            && query_ == other.query_
            && sequence_ == other.sequence_
            && cigar_ == other.cigar_
            && nodes_ == other.nodes_;
    }

    bool operator!=(const Alignment &other) const { return !(*this == other); }

    bool is_exact_match() const {
        return cigar_.size() == 1
            && cigar_.front() == Cigar::value_type(Cigar::MATCH, query_.size());
    }

    Json::Value to_json(std::string_view query,
                        const DeBruijnGraph &graph,
                        bool is_secondary = false,
                        std::string_view name = {},
                        std::string_view label = {}) const;

    // returns a shared_ptr of the query string which is referenced in this object
    std::shared_ptr<const std::string> load_from_json(const Json::Value &alignment,
                                                      const DeBruijnGraph &graph);

    bool is_valid(const DeBruijnGraph &graph, const DBGAlignerConfig *config = nullptr) const;

  private:
    Json::Value path_json(size_t node_size, std::string_view label = {}) const;

    std::string_view query_;
    std::vector<node_index> nodes_;
    std::string sequence_;
    score_t score_;
    Cigar cigar_;
    bool orientation_;
    size_t offset_;
};

std::ostream& operator<<(std::ostream& out, const Alignment &alignment);

struct LocalAlignmentLess {
    bool operator()(const Alignment &a, const Alignment &b) const {
        // 1) score is less, or
        // 2) more of the query is covered, or
        // 3) if it is in the reverse orientation, or
        // 4) if the starting point is later in the query
        return std::make_tuple(b.get_score(), a.get_query().size(),
                               a.get_orientation(), a.get_clipping())
            > std::make_tuple(a.get_score(), b.get_query().size(),
                              b.get_orientation(), b.get_clipping());
    }
};

struct LocalAlignmentGreater {
    bool operator()(const Alignment &a, const Alignment &b) const {
        // 1) score is higher, or
        // 2) less of the query is covered, or
        // 3) if it is in the forward orientation, or
        // 4) if the starting point is earlier in the query
        return std::make_tuple(a.get_score(), b.get_query().size(),
                               b.get_orientation(), b.get_clipping())
            > std::make_tuple(b.get_score(), a.get_query().size(),
                              a.get_orientation(), a.get_clipping());
    }
};

// A container holding many alignments to a shared query sequence. Each alignment
// only holds a string_view to the query, so this class ensures that the query sequence
// is always accessible.
// TODO: rename to AlignmentResults
class QueryAlignment {
  public:
    typedef std::vector<Alignment>::iterator iterator;
    typedef std::vector<Alignment>::const_iterator const_iterator;

    explicit QueryAlignment(std::string_view query, bool is_reverse_complement = false);

    explicit QueryAlignment(std::shared_ptr<const std::string> query,
                            std::shared_ptr<const std::string> query_rc)
          : query_(query), query_rc_(query_rc) {}

    template <typename... Args>
    void emplace_back(Args&&... args) {
        alignments_.emplace_back(std::forward<Args>(args)...);

        assert(alignments_.back().get_query().data()
            >= get_query(alignments_.back().get_orientation()).c_str());
        assert(alignments_.back().get_query().data() + alignments_.back().get_query().size()
            <= get_query(alignments_.back().get_orientation()).c_str()
                + get_query(alignments_.back().get_orientation()).size());
    }

    std::shared_ptr<const std::string> get_query_ptr(bool reverse_complement = false) const {
        return !reverse_complement ? query_ : query_rc_;
    }

    const std::string& get_query(bool reverse_complement = false) const {
        return *get_query_ptr(reverse_complement);
    }

    size_t size() const { return alignments_.size(); }
    bool empty() const { return alignments_.empty(); }
    const Alignment& operator[](size_t i) const { return alignments_[i]; }
    iterator begin() { return alignments_.begin(); }
    iterator end() { return alignments_.end(); }
    const_iterator begin() const { return alignments_.begin(); }
    const_iterator end() const { return alignments_.end(); }

    std::vector<Alignment>& data() { return alignments_; }
    const std::vector<Alignment>& data() const { return alignments_; }

  private:
    std::shared_ptr<const std::string> query_;
    std::shared_ptr<const std::string> query_rc_;
    std::vector<Alignment> alignments_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif  // __ALIGNER_ALIGNMENT_HPP__
