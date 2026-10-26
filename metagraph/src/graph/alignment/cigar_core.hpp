#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <numeric>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mtg::graph::detail {

// Common representation and operations; each aligner retains its query contract.
class CigarCore {
  public:
    enum Operator : int8_t {
        CLIPPED,
        MISMATCH,
        MATCH,
        DELETION,
        INSERTION,
        NODE_INSERTION
    };

    using LengthType = uint32_t;
    using value_type = std::pair<Operator, LengthType>;

    size_t size() const { return cigar_.size(); }
    bool empty() const { return cigar_.empty(); }

    std::string to_string() const;
    void append(Operator op, LengthType num = 1);

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

    size_t get_num_matches() const {
        return std::accumulate(cigar_.begin(), cigar_.end(), size_t{0},
                               [&](size_t old, const value_type &op) {
            return old + (op.first == MATCH) * op.second;
        });
    }

    bool is_exact_match(size_t query_size) const {
        return cigar_.size() == 1 && cigar_.front() == value_type{ MATCH, query_size };
    }

    static constexpr char opt_to_char(Operator op) { return op_str_[op]; }

  protected:
    CigarCore(Operator op, LengthType num)
          : cigar_(num ? 1 : 0, std::make_pair(op, num)) {}

    explicit CigarCore(std::vector<value_type> operations)
          : cigar_(std::move(operations)) {}

    // SAM operations plus G for graph-node insertion; lengths must be positive.
    explicit CigarCore(std::string_view cigar_str);

    void append_cigar(CigarCore&& other);

    std::vector<value_type> cigar_;

  private:
    static constexpr char op_str_[] = "SX=DIG";
};

} // namespace mtg::graph::detail
