#include "hll_counter.hpp"

#include <numeric>
#include <cmath>

#if __AVX2__
#include <immintrin.h>
#endif

#include "common/serialization.hpp"
#include "common/utils/simd_utils.hpp"


inline int compute_precision(double error) {
    // want 1.04 / std::sqrt(m) < error, where m = num_buckets
    // 1.04 / error < std::sqrt(m)
    // 1.04 * 1.04 / error / error < m = 2^precision
    // std::log2(1.04 * 1.04 / error / error) < precision

    if (error < 0.00203125)
        std::cerr << "WARNING: precision changed to 0.00203125" << std::endl;

    return std::min(18, std::max(4, static_cast<int>(
        std::ceil(std::log2(1.04 * 1.04 / error / error))
    )));
}

HLLCounter::HLLCounter(double error) : HLLCounter(compute_precision(error)) {}
HLLCounter::HLLCounter(int precision) : precision_(precision) {
    if (precision_ < 4 || precision_ > 18)
        throw std::runtime_error("HLL precision not set correctly: " + std::to_string(precision_));

    assert(precision_ >= 4);
    assert(precision_ <= 18);

    counter_ = libcount::HLL::Create(precision_);
    assert(counter_);
}

HLLCounter& HLLCounter::operator=(const HLLCounter &other) {
    precision_ = other.precision_;
    counter_.reset();

    if (!other.counter_)
        return *this;

    counter_ = libcount::HLL::Create(precision_);

    if (!counter_)
        return *this;

    auto *registers = counter_->data();
    const auto *other_registers = other.counter_->data();

    memcpy(registers, other_registers, counter_->get_register_count());

    return *this;
}

void HLLCounter::insert(uint64_t hash) {
    assert(counter_);
    assert(precision_ >= 4);
    assert(precision_ <= 18);

    counter_->Update(hash);
}

void HLLCounter::insert(const uint64_t *hashes_begin, const uint64_t *hashes_end) {
    assert(counter_);
    assert(precision_ >= 4);
    assert(precision_ <= 18);

    counter_->UpdateMany(hashes_begin, hashes_end - hashes_begin);
}

bool HLLCounter::check(uint64_t hash) const {
    // TODO: better implementation
    assert(counter_);
    assert(precision_ >= 4);
    assert(precision_ <= 18);

    HLLCounter query(precision_);
    query.insert(hash);
    return std::round(estimate_intersection_cardinality(query)) >= 1.0;
}

void HLLCounter::check(const uint64_t *hashes_begin,
                       const uint64_t *hashes_end,
                       const std::function<void(size_t)> &present_index_callback) const {
    // TODO: AVX2
    assert(hashes_end >= hashes_begin);

    size_t count = hashes_end - hashes_begin;
    for (size_t i = 0; i < count; ++i) {
        if (check(hashes_begin[i]))
            present_index_callback(i);
    }
}

void HLLCounter::merge(const HLLCounter &other) {
    assert(counter_);
    assert(other.counter_);
    assert(precision_ >= 4);
    assert(precision_ <= 18);
    assert(other.precision_ >= 4);
    assert(other.precision_ <= 18);

    counter_->Merge(other.counter_.get());
}

double HLLCounter::estimate_cardinality() const {
    return counter_->Estimate();
}

double HLLCounter::estimate_union_cardinality(const HLLCounter &other) const {
    auto merged = libcount::HLL::Create(precision_);
    merged->Merge(counter_.get());
    merged->Merge(other.counter_.get());
    return merged->Estimate();
}

double HLLCounter::estimate_intersection_cardinality(const HLLCounter &other) const {
    // size(union(a,b)) = size(a) + size(b) - size(intersection(a,b))
    // size(intersection(a,b)) = size(a) + size(b) - size(union(a,b))

    return estimate_cardinality() + other.estimate_cardinality()
        - estimate_union_cardinality(other);
}

double HLLCounter::estimate_jaccard(const HLLCounter &other) const {
    // size(intersection(a,b)) / size(union(a,b))
    // size(size(a) + size(b) - size of union(a,b)) / size(union(a,b))
    double union_cardinality = estimate_union_cardinality(other);
    return (estimate_cardinality() + other.estimate_cardinality() - union_cardinality)
        / union_cardinality;
}

void HLLCounter::serialize(std::ostream &out) const {
    assert(counter_);
    assert(precision_ >= 4);
    assert(precision_ <= 18);

    serialize_number(out, counter_->get_precision());
    serialize_number(out, counter_->get_register_count());
    out.write(reinterpret_cast<const char*>(counter_->data()),
              counter_->get_register_count());
}

bool HLLCounter::load(std::istream &in) {
    try {
        precision_ = load_number(in);
        int register_count = load_number(in);
        std::ignore = register_count;

        if (precision_ < 4 || precision_ > 18) {
            throw std::runtime_error(
                "Invalid HLL precision read: " + std::to_string(precision_)
            );
        }

        counter_ = libcount::HLL::Create(precision_);

        if (!counter_)
            return false;

        assert(counter_->get_precision() == precision_);
        assert(counter_->get_register_count() == register_count);
        in.read(reinterpret_cast<char*>(counter_->data()), counter_->get_register_count());

        return true;
    } catch (...) {
        return false;
    }
}

void HLLCounter::reset() {
    assert(counter_);
    std::fill(counter_->data(), counter_->data() + counter_->get_register_count(), 0);
}
