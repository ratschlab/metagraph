#ifndef __DATATYPES_HPP__
#define __DATATYPES_HPP__

#include <vector>
#include <string>


struct HitInfo {
    uint64_t rl;
    uint64_t ru;
    uint64_t str_pos;
    uint64_t graph_pos;
    uint64_t distance;
    std::string cigar;
    std::vector<uint64_t> path;
};


class HitInfoCompare {
  public:
    explicit HitInfoCompare(bool is_reverse) : is_reverse_(is_reverse) {}
    HitInfoCompare() : HitInfoCompare::HitInfoCompare(false) {}

    bool operator()(const HitInfo &lhs, const HitInfo &rhs) const {
        return is_reverse_ ? lhs.distance < rhs.distance
                           : lhs.distance > rhs.distance;
    }

  private:
    bool is_reverse_;
};


#endif // __DATATYPES_HPP__
