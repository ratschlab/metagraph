#ifndef __CIGAR_HPP__
#define __CIGAR_HPP__

#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <sdsl/int_vector.hpp>

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
        INSERTION,
        NODE_INSERTION
    };

    typedef uint32_t LengthType;
    typedef std::pair<Operator, LengthType> value_type;

    Cigar(Operator op = CLIPPED, LengthType num = 0)
          : cigar_(num ? 1 : 0, std::make_pair(op, num)) { }

    // See section 1.4 in https://samtools.github.io/hts-specs/SAMv1.pdf for
    // a specification of the CIGAR string format.
    // e.g., 3=1X2I3D for 3 matches, 1 mismatch, 2 insertions, 3 deletions
    // The symbol 'G' is introduced to indicate the insertion of a graph node.
    Cigar(std::string_view cigar_str);

    size_t size() const { return cigar_.size(); }
    bool empty() const { return cigar_.empty(); }

    std::string to_string() const;
    std::string to_md_string(std::string_view reference) const;

    void append(Operator op, LengthType num = 1);
    void append(Cigar&& other);

    LengthType trim_clipping() {
        if (cigar_.size() && cigar_.front().first == CLIPPED) {
            LengthType ret_val = cigar_.front().second;
            cigar_.erase(cigar_.begin(), cigar_.begin() + 1);
            return ret_val;
        } else {
            return 0;
        }
    }

    LengthType trim_end_clipping() {
        if (cigar_.size() && cigar_.back().first == CLIPPED) {
            LengthType ret_val = cigar_.back().second;
            cigar_.pop_back();
            return ret_val;
        } else {
            return 0;
        }
    }

    LengthType get_clipping() const {
        return cigar_.size() && cigar_.front().first == CLIPPED ? cigar_.front().second : 0;
    }

    LengthType get_end_clipping() const {
        return cigar_.size() && cigar_.back().first == CLIPPED ? cigar_.back().second : 0;
    }

    size_t get_coverage() const;

    void extend_clipping(LengthType n) {
        assert(cigar_.size());
        if (cigar_.front().first != CLIPPED) {
            cigar_.insert(cigar_.begin(), value_type(CLIPPED, n));
        } else {
            cigar_.front().second += n;
        }
    }

    std::vector<value_type>& data() { return cigar_; }
    const std::vector<value_type>& data() const { return cigar_; }

    bool operator==(const Cigar &other) const { return cigar_ == other.cigar_; }
    bool operator!=(const Cigar &other) const { return !(*this == other); }

    size_t get_num_matches() const {
        return std::accumulate(cigar_.begin(), cigar_.end(), 0,
                               [&](size_t old, const value_type &op) {
            return old + (op.first == MATCH) * op.second;
        });
    }

    // Return true if the cigar is valid
    bool is_valid(std::string_view reference, std::string_view query) const;

    bool is_exact_match(size_t query_size) const {
        return cigar_.size() == 1 && cigar_.front() == value_type{ MATCH, query_size };
    }

    void mark_exact_matches(sdsl::bit_vector &mask) const;

    static constexpr char opt_to_char(Cigar::Operator op) { return op_str_[op]; }

  private:
    static constexpr char op_str_[] = "SX=DIG";
    std::vector<value_type> cigar_;
};


typedef std::array<std::array<Cigar::Operator, 128>, 128> OperatorTable;
OperatorTable initialize_opt_table();
static const OperatorTable kCharToOp = initialize_opt_table();


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __CIGAR_HPP__
