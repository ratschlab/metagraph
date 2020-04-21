#ifndef __RANGE_PARTITION_HPP__
#define __RANGE_PARTITION_HPP__

#include <vector>
#include <cassert>
#include <cstdint>
#include <iostream>


// Partitions a range of numbers in [0,n) into groups.
// The groups must not be empty.
class RangePartition {
  public:
    typedef uint32_t T;
    typedef uint32_t G;
    typedef uint32_t R;

    RangePartition() {}
    RangePartition(const std::vector<uint64_t> &arrangement,
                   const std::vector<size_t> &group_sizes);
    explicit RangePartition(std::vector<std::vector<uint64_t>>&& partition);

    // get group that contains value
    inline G group(T value) const;

    // get index of value in its group
    inline R rank(T value) const;

    // get value given its group and the rank
    inline T get(G group, R rank) const;

    inline size_t num_groups() const;
    inline size_t group_size(G group) const;
    inline size_t size() const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

  private:
    // Based on |partition_|, initializes groups and ranks.
    // Returns false if partition is invalid.
    bool initialize_groups_and_ranks();

    std::vector<std::vector<T>> partition_;
    std::vector<G> groups_;
    std::vector<R> ranks_;
};

RangePartition::G RangePartition::group(T value) const {
    assert(value < groups_.size());
    return groups_[value];
}

RangePartition::R RangePartition::rank(T value) const {
    assert(value < ranks_.size());
    return ranks_[value];
}

RangePartition::T RangePartition::get(G group, R rank) const {
    assert(group < partition_.size());
    assert(rank < partition_[group].size());
    return partition_[group][rank];
}

size_t RangePartition::num_groups() const {
    return partition_.size();
}

size_t RangePartition::group_size(G group) const {
    assert(group < partition_.size());
    return partition_[group].size();
}

size_t RangePartition::size() const {
    assert(groups_.size() == ranks_.size());
    return ranks_.size();
}

#endif // __RANGE_PARTITION_HPP__
