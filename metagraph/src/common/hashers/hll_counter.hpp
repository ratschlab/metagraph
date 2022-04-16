#ifndef __HLL_COUNTER_HPP__
#define __HLL_COUNTER_HPP__

#include <functional>
#include <iostream>
#include <vector>

#include <sdsl/int_vector.hpp>

#include <count/hll.h>


class HLLCounter {
  public:
    HLLCounter() {}
    HLLCounter(double error);
    HLLCounter(const HLLCounter &other) { *this = other; }

    HLLCounter& operator=(const HLLCounter &other);

    void insert(uint64_t hash);
    void insert(const uint64_t *hashes_begin, const uint64_t *hashes_end);

    bool check(uint64_t hash) const;
    void check(const uint64_t *hashes_begin,
               const uint64_t *hashes_end,
               const std::function<void(size_t)> &present_index_callback) const;

    inline sdsl::bit_vector check(const uint64_t *hashes_begin,
                                  const uint64_t *hashes_end) const {
        assert(hashes_end >= hashes_begin);

        sdsl::bit_vector presence(hashes_end - hashes_begin, false);
        check(hashes_begin, hashes_end, [&](size_t i) { presence[i] = true; });
        return presence;
    }

    void merge(const HLLCounter &other);

    void serialize(std::ostream &out) const;
    bool load(std::istream &in);

    double estimate_cardinality() const;
    double estimate_union_cardinality(const HLLCounter &other) const;
    double estimate_intersection_cardinality(const HLLCounter &other) const;
    double estimate_jaccard(const HLLCounter &other) const;

    void reset();

    const libcount::HLL& data() const { return *counter_; }
    libcount::HLL& data() { return *counter_; }

    // bool operator==(const HLLCounter &other) const {
    //     return std::equal(counter_->data(),
    //                       counter_->data() + counter_->get_register_count(),
    //                       other.counter_->data(),
    //                       other.counter_->data() + other.counter_->get_register_count());
    // }

    // bool operator!=(const HLLCounter &other) const {
    //     return !operator==(other);
    // }

  private:
    int precision_;
    std::unique_ptr<libcount::HLL> counter_;

    explicit HLLCounter(int precision);
};

#endif // __HLL_COUNTER_HPP__
