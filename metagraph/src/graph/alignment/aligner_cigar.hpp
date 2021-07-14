#ifndef __CIGAR_HPP__
#define __CIGAR_HPP__

#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

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

    Cigar(Operator op = Operator::CLIPPED, LengthType num = 0)
          : cigar_(num ? 1 : 0, std::make_pair(op, num)) { }

    // See section 1.4 in https://samtools.github.io/hts-specs/SAMv1.pdf for
    // a specification of the CIGAR string format.
    // e.g., 3=1X2I3D for 3 matches, 1 mismatch, 2 insertions, 3 deletions
    Cigar(std::string_view cigar_str);

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

    template <typename... Args>
    void erase(Args&&... args) { cigar_.erase(std::forward<Args>(args)...); }

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
    bool is_valid(std::string_view reference, std::string_view query) const;

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
