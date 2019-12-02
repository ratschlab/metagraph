#ifndef __RANGE_PARTITION_HPP__
#define __RANGE_PARTITION_HPP__

#include <vector>
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

    explicit RangePartition(const RangePartition &) = default;
    RangePartition& operator=(const RangePartition &) = default;
    RangePartition(RangePartition&&) = default;
    RangePartition& operator=(RangePartition&&) = default;

    // get group that contains value
    G group(T value) const;

    // get index of value in its group
    R rank(T value) const;

    // get value given its group and the rank
    T get(G group, R rank) const;

    uint64_t num_groups() const;
    uint64_t size() const;

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

#endif // __RANGE_PARTITION_HPP__
