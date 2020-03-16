#ifndef __BIT_VECTOR_STAT_HPP__
#define __BIT_VECTOR_STAT_HPP__

#include <cstdint>

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/select_support_mcl.hpp>

#include "common/serialization.hpp"
#include "vector_algorithm.hpp"
#include "bit_vector.hpp"


class bit_vector_stat : public bit_vector {
    friend bit_vector;

    static constexpr size_t MAX_ITER_BIT_VECTOR_STAT = 1000;
    static constexpr double BITS_PER_BIT_IF_SPARSE = 1.07;
    // FYI: rank support takes 0.0625 bits per bit, the other ~0.31 are taken by select
    static constexpr double BITS_PER_BIT_IF_DENSE = 1.37;

  public:
    inline explicit bit_vector_stat(uint64_t size = 0, bool value = 0);
    inline explicit bit_vector_stat(const sdsl::bit_vector &vector);
    inline explicit bit_vector_stat(const bit_vector_stat &other);
    inline bit_vector_stat(sdsl::bit_vector&& vector);
    inline bit_vector_stat(bit_vector_stat&& other) noexcept;
    inline bit_vector_stat(std::initializer_list<bool> init);

    inline bit_vector_stat& operator=(const bit_vector_stat &other);
    inline bit_vector_stat& operator=(bit_vector_stat&& other) noexcept;

    inline std::unique_ptr<bit_vector> copy() const override;

    inline uint64_t rank1(uint64_t id) const override;
    inline uint64_t select1(uint64_t id) const override;
    inline uint64_t select0(uint64_t id) const override;

    inline uint64_t next1(uint64_t id) const override;
    inline uint64_t prev1(uint64_t id) const override;

    inline bool operator[](uint64_t id) const override;
    inline uint64_t get_int(uint64_t id, uint32_t width) const override;

    inline bool load(std::istream &in) override;
    inline void serialize(std::ostream &out) const override;

    inline uint64_t size() const override { return vector_.size(); }
    inline uint64_t num_set_bits() const override { return num_set_bits_; }

    inline void call_ones_in_range(uint64_t begin, uint64_t end,
                                   const VoidCall<uint64_t> &callback) const override;

    inline void add_to(sdsl::bit_vector *other) const override;

    inline sdsl::bit_vector to_vector() const override { return vector_; }

    inline const sdsl::bit_vector& data() const { return vector_; }

    /**
     * Predict space taken by the vector with given its parameters in bits.
     */
    static uint64_t predict_size(uint64_t size, uint64_t num_set_bits) {
        return BITS_PER_BIT_IF_SPARSE * size
                + (BITS_PER_BIT_IF_DENSE - BITS_PER_BIT_IF_SPARSE) * num_set_bits;
    }

  private:
    sdsl::bit_vector vector_;
    uint64_t num_set_bits_ = 0;

    // maintain rank/select operations
    sdsl::rank_support_v5<1> rk_;
    sdsl::select_support_mcl<1> slct_;
    sdsl::select_support_scan<0> slct_zero_;
};


bit_vector_stat::bit_vector_stat(uint64_t size, bool value)
      : bit_vector_stat(sdsl::bit_vector(size, value)) {}

bit_vector_stat::bit_vector_stat(const sdsl::bit_vector &vector)
      : bit_vector_stat(sdsl::bit_vector(vector)) {}

bit_vector_stat::bit_vector_stat(const bit_vector_stat &other) {
    *this = other;
}

bit_vector_stat::bit_vector_stat(sdsl::bit_vector&& vector)
      : vector_(std::move(vector)) {
    rk_ = sdsl::rank_support_v5<1>(&vector_);
    slct_ = sdsl::select_support_mcl<1>(&vector_);
    slct_zero_ = sdsl::select_support_scan<0>(&vector_);
    num_set_bits_ = rk_(vector_.size());
    assert(num_set_bits_ == num_set_bits());
}

bit_vector_stat::bit_vector_stat(bit_vector_stat&& other) noexcept {
    *this = std::move(other);
}

bit_vector_stat::bit_vector_stat(std::initializer_list<bool> init)
      : bit_vector_stat(sdsl::bit_vector(init)) {}

bit_vector_stat& bit_vector_stat::operator=(const bit_vector_stat &other) {
    vector_ = other.vector_;
    num_set_bits_ = other.num_set_bits_;

    rk_ = other.rk_;
    rk_.set_vector(&vector_);
    slct_ = other.slct_;
    slct_.set_vector(&vector_);
    slct_zero_ = other.slct_zero_;
    slct_zero_.set_vector(&vector_);

    return *this;
}

bit_vector_stat& bit_vector_stat::operator=(bit_vector_stat&& other) noexcept {
    vector_ = std::move(other.vector_);
    num_set_bits_ = other.num_set_bits_;

    rk_ = std::move(other.rk_);
    rk_.set_vector(&vector_);
    slct_ = std::move(other.slct_);
    slct_.set_vector(&vector_);
    slct_zero_ = std::move(other.slct_zero_);
    slct_zero_.set_vector(&vector_);

    return *this;
}

std::unique_ptr<bit_vector> bit_vector_stat::copy() const {
    return std::make_unique<bit_vector_stat>(*this);
}

uint64_t bit_vector_stat::rank1(uint64_t id) const {
    //the rank method in SDSL does not include id in the count
    return rk_(id >= this->size() ? this->size() : id + 1);
}

uint64_t bit_vector_stat::select1(uint64_t id) const {
    assert(id > 0 && id <= num_set_bits_);
    return slct_(id);
}

uint64_t bit_vector_stat::select0(uint64_t id) const {
    assert(id > 0 && id + num_set_bits_ <= size());
    return slct_zero_(id);
}

uint64_t bit_vector_stat::next1(uint64_t pos) const {
    assert(pos < size());

    auto next = next_bit(vector_, pos, MAX_ITER_BIT_VECTOR_STAT);
    if (next < vector_.size())
        return next;

    if (vector_.size() - pos <= MAX_ITER_BIT_VECTOR_STAT)
        return size();

    uint64_t rk = rank1(pos) + 1;
    return rk <= num_set_bits() ? select1(rk) : size();
}

uint64_t bit_vector_stat::prev1(uint64_t pos) const {
    assert(pos < size());

    auto prev = prev_bit(vector_, pos, MAX_ITER_BIT_VECTOR_STAT);
    if (prev <= pos)
        return prev;

    if (pos < MAX_ITER_BIT_VECTOR_STAT)
        return size();

    uint64_t rk = rank1(pos);
    return rk ? select1(rk) : size();
}

bool bit_vector_stat::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
}

uint64_t bit_vector_stat::get_int(uint64_t id, uint32_t width) const {
    return vector_.get_int(id, width);
}

void bit_vector_stat::serialize(std::ostream &out) const {
    vector_.serialize(out);

    serialize_number(out, num_set_bits_);
    rk_.serialize(out);
    slct_.serialize(out);
    slct_zero_.serialize(out);
}

bool bit_vector_stat::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);

        num_set_bits_ = load_number(in);
        rk_.load(in, &vector_);
        slct_.load(in, &vector_);
        slct_zero_.load(in, &vector_);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_stat." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void bit_vector_stat::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());
    *other |= vector_;
}

void bit_vector_stat::call_ones_in_range(uint64_t begin, uint64_t end,
                                         const VoidCall<uint64_t> &callback) const {
    assert(begin <= end);
    assert(end <= size());
    ::call_ones(vector_, begin, end, callback);
}

#endif // __BIT_VECTOR_STAT_HPP__
