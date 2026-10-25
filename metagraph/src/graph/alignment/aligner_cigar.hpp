#ifndef __CIGAR_HPP__
#define __CIGAR_HPP__

#include <array>
#include <iostream>

#include "cigar_core.hpp"

namespace mtg {
namespace graph {
namespace align {

class Cigar : public detail::CigarCore {
  public:
    Cigar(Operator op = CLIPPED, LengthType num = 0) : CigarCore(op, num) {}
    Cigar(std::string_view cigar_str) : CigarCore(cigar_str) {}
    explicit Cigar(std::vector<value_type> operations) : CigarCore(std::move(operations)) {}

    using CigarCore::append;
    void append(Cigar&& other) { append_cigar(std::move(other)); }

    bool operator==(const Cigar &other) const { return cigar_ == other.cigar_; }
    bool operator!=(const Cigar &other) const { return !(*this == other); }

    // Query is the aligned slice: terminal clipping is not consumed.
    bool is_valid(std::string_view reference, std::string_view query) const;
};


typedef std::array<std::array<Cigar::Operator, 128>, 128> OperatorTable;
OperatorTable initialize_opt_table();
static const OperatorTable kCharToOp = initialize_opt_table();


} // namespace align
} // namespace graph
} // namespace mtg

#endif // __CIGAR_HPP__
