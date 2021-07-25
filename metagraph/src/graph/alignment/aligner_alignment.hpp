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

class DeBruijnGraph;

namespace align {

// Note: this object stores pointers to the query sequence, so it is the user's
//       responsibility to ensure that the query sequence is not destroyed when
//       calling this class' methods
class Alignment {
  public:
    typedef DeBruijnGraph::node_index node_index;
    typedef DBGAlignerConfig::score_t score_t;

    Alignment(std::string_view query,
              std::vector<node_index>&& nodes = {},
              std::string&& sequence = "",
              score_t score = 0,
              Cigar&& cigar = Cigar(),
              size_t clipping = 0,
              bool orientation = false,
              size_t offset = 0)
          : query_(query),
            nodes_(std::move(nodes)),
            sequence_(std::move(sequence)),
            score_(score),
            cigar_(Cigar::CLIPPED, clipping),
            orientation_(orientation),
            offset_(offset) { cigar_.append(std::move(cigar)); }

    // Used for constructing seeds
    Alignment(std::string_view query = {},
              std::vector<node_index>&& nodes = {},
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
        assert(nodes.empty() || clipping || cigar_.is_exact_match(query_.size()));
    }

    // Used for constructing exact match seeds
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
    const node_index& operator[](size_t i) const { return nodes_[i]; }
    const node_index& front() const { return nodes_.front(); }
    const node_index& back() const { return nodes_.back(); }

    score_t get_score() const { return score_; }
    uint64_t get_num_matches() const { return cigar_.get_num_matches(); }

    std::string_view get_query() const { return query_; }

    void set_query_begin(const char *begin) { query_ = { begin, query_.size() }; }

    void extend_query_begin(const char *begin) {
        const char *full_query_begin = query_.data() - get_clipping();
        assert(full_query_begin >= begin);
        if (full_query_begin > begin)
            cigar_.extend_clipping(full_query_begin - begin);
    }

    void extend_query_end(const char *end) {
        const char *full_query_end = query_.data() + query_.size() + get_end_clipping();
        assert(full_query_end <= end);
        if (full_query_end < end)
            cigar_.append(Cigar::CLIPPED, end - full_query_end);
    }

    inline size_t trim_clipping() { return cigar_.trim_clipping(); }
    inline size_t trim_end_clipping() { return cigar_.trim_end_clipping(); }

    size_t trim_offset();

    void reverse_complement(const DeBruijnGraph &graph,
                            std::string_view query_rev_comp);

    const std::string& get_sequence() const { return sequence_; }
    const Cigar& get_cigar() const { return cigar_; }
    bool get_orientation() const { return orientation_; }
    size_t get_offset() const { return offset_; }
    Cigar::LengthType get_clipping() const { return cigar_.get_clipping(); }
    Cigar::LengthType get_end_clipping() const { return cigar_.get_end_clipping(); }

    typedef typename std::vector<node_index>::iterator iterator;
    typedef typename std::vector<node_index>::const_iterator const_iterator;

    const_iterator begin() const { return nodes_.cbegin(); }
    const_iterator end() const { return nodes_.cend(); }

    bool operator==(const Alignment &other) const {
        return orientation_ == other.orientation_
            && score_ == other.score_
            && sequence_ == other.sequence_
            && query_ == other.query_
            && cigar_ == other.cigar_;
    }

    bool operator!=(const Alignment &other) const { return !(*this == other); }

    Json::Value to_json(std::string_view query,
                        const DeBruijnGraph &graph,
                        bool is_secondary = false,
                        std::string_view name = {},
                        std::string_view label = {}) const;

    std::shared_ptr<const std::string>
    load_from_json(const Json::Value &alignment,
                   const DeBruijnGraph &graph);

    bool is_valid(const DeBruijnGraph &graph, const DBGAlignerConfig *config = nullptr) const;

  private:
    Json::Value path_json(size_t node_size, std::string_view label = {}) const;

    // TODO: rename to query_view_
    std::string_view query_;
    std::vector<node_index> nodes_;
    std::string sequence_;
    score_t score_;
    Cigar cigar_;
    bool orientation_;
    size_t offset_;
};

inline std::ostream& operator<<(std::ostream& out, const Alignment &alignment) {
    out << (alignment.get_orientation() ? "-" : "+") << "\t"
        << alignment.get_sequence() << "\t"
        << alignment.get_score() << "\t"
        << alignment.get_num_matches() << "\t"
        << alignment.get_cigar().to_string() << "\t"
        << alignment.get_offset();

    return out;
}

bool spell_path(const DeBruijnGraph &graph,
                const std::vector<uint64_t> &path,
                std::string &seq,
                size_t offset = 0);

struct LocalAlignmentLess {
    bool operator()(const Alignment &a, const Alignment &b) {
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
    bool operator()(const Alignment &a, const Alignment &b) {
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


class QueryAlignment {
  public:
    typedef typename std::vector<Alignment>::const_iterator const_iterator;

    QueryAlignment(std::string_view query, bool is_reverse_complement = false);

    size_t size() const { return alignments_.size(); }
    bool empty() const { return alignments_.empty(); }

    template <typename... Args>
    void emplace_back(Args&&... args) {
        alignments_.emplace_back(std::forward<Args>(args)...);

#ifndef NDEBUG
        const auto &added = alignments_.back();
        const std::string &this_query = get_query(added.get_orientation());
        assert(added.get_query().data() >= this_query.c_str());
        assert(added.get_query().data() + added.get_query().size()
                    <= this_query.c_str() + this_query.size());
#endif
    }

    void pop_back() { alignments_.pop_back(); }
    void clear() { alignments_.clear(); }

    const std::string& get_query(bool reverse_complement = false) const {
        return !reverse_complement ? *query_ : *query_rc_;
    }

    const Alignment& operator[](size_t i) const { return alignments_[i]; }
    const_iterator begin() const { return alignments_.cbegin(); }
    const_iterator end() const { return alignments_.cend(); }
    const_iterator cbegin() const { return alignments_.cbegin(); }
    const_iterator cend() const { return alignments_.cend(); }

  private:
    std::shared_ptr<std::string> query_;
    std::shared_ptr<std::string> query_rc_;
    std::vector<Alignment> alignments_;
};

} // namespace align
} // namespace graph
} // namespace mtg

#endif  // __ALIGNER_ALIGNMENT_HPP__
