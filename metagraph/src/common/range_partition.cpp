#include "range_partition.hpp"

#include <cassert>

#include "common/serialization.hpp"


RangePartition::RangePartition(const std::vector<uint64_t> &arrangement,
                               const std::vector<size_t> &group_sizes) {
    size_t offset = 0;
    for (size_t group_size : group_sizes) {
        assert(group_size && "partition blocks must not be empty");
        partition_.emplace_back(arrangement.begin() + offset,
                                arrangement.begin() + offset + group_size);
        offset += group_size;
    }

    assert(initialize_groups_and_ranks());
    initialize_groups_and_ranks();
}

RangePartition::RangePartition(std::vector<std::vector<uint64_t>>&& partition) {
    partition_.reserve(partition.size());
    for (auto &group : partition) {
        assert(group.size() && "partition blocks must not be empty");
        partition_.emplace_back(group.begin(), group.end());
        group.clear();
    }
    partition.clear();

    assert(initialize_groups_and_ranks());
    initialize_groups_and_ranks();
}

bool RangePartition::initialize_groups_and_ranks() {
    uint64_t set_size = 0;
    for (const auto &group : partition_) {
        set_size += group.size();
        if (!group.size())
            return false;
    }

    groups_.assign(set_size, -1);
    ranks_.assign(set_size, -1);

    for (size_t g = 0; g < partition_.size(); ++g) {
        const auto &group = partition_[g];

        for (size_t i = 0; i < group.size(); ++i) {
            auto value = group[i];
            if (value >= set_size
                    || groups_.at(value) != static_cast<T>(-1)
                    || ranks_.at(value) != static_cast<T>(-1))
                return false;

            groups_[value] = g;
            ranks_[value] = i;
        }
    }
    return true;
}

bool RangePartition::load(std::istream &in) {
    if (!in.good())
        return false;

    partition_.assign(load_number(in), {});
    for (auto &group : partition_) {
        if (!load_number_vector(in, &group))
            return false;
    }
    return initialize_groups_and_ranks();
}

void RangePartition::serialize(std::ostream &out) const {
    if (!out.good())
        throw std::ofstream::failure("Bad stream");

    serialize_number(out, partition_.size());
    for (const auto &group : partition_) {
        serialize_number_vector(out, group);
    }
}
