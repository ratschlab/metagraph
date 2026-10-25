#pragma once

#include <iostream>

#include "graph/alignment/cigar_core.hpp"

namespace mtg::graph::align_redone {

class Cigar : public detail::CigarCore {
  public:
    Cigar(Operator op = CLIPPED, LengthType num = 0) : CigarCore(op, num) {}
    Cigar(std::string_view cigar_str) : CigarCore(cigar_str) {}
    explicit Cigar(std::vector<value_type> operations) : CigarCore(std::move(operations)) {}

    using CigarCore::append;
    void append(Cigar&& other) { append_cigar(std::move(other)); }

    void extend_clipping(LengthType n) {
        if (n)
            CigarCore::extend_clipping(n);
    }

    bool operator==(const Cigar &other) const { return cigar_ == other.cigar_; }
    bool operator!=(const Cigar &other) const { return !(*this == other); }

    size_t get_num_query() const {
        return std::accumulate(cigar_.begin(), cigar_.end(), size_t{0},
                               [&](size_t old, const value_type &op) {
            return old + (op.first == MATCH || op.first == MISMATCH || op.first == INSERTION || op.first == CLIPPED) * op.second;
        });
    }

    size_t get_num_ref() const {
        return std::accumulate(cigar_.begin(), cigar_.end(), size_t{0},
                               [&](size_t old, const value_type &op) {
            return old + (op.first == MATCH || op.first == MISMATCH || op.first == DELETION) * op.second;
        });
    }

    // Query includes every clipped base, including permitted graph-jump clips.
    bool is_valid(std::string_view reference, std::string_view query) const;
};


} // namespace mtg::graph::align_redone
