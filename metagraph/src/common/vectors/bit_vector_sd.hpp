#ifndef __BIT_VECTOR_SD_HPP__
#define __BIT_VECTOR_SD_HPP__

#include <cstdint>

#include <sdsl/sd_vector.hpp>

#include "vector_algorithm.hpp"
#include "bit_vector.hpp"


class bit_vector_sd : public bit_vector {
    static constexpr size_t MAX_ITER_BIT_VECTOR_SD = 10;

  public:
    inline explicit bit_vector_sd(uint64_t size = 0, bool value = false);
    inline explicit bit_vector_sd(const sdsl::bit_vector &vector);
    inline explicit bit_vector_sd(const sdsl::bit_vector &vector, uint64_t num_set_bits);
    inline explicit bit_vector_sd(const bit_vector_sd &other);

    inline bit_vector_sd(bit_vector_sd&& other);
    inline bit_vector_sd(std::initializer_list<bool> init);
    inline bit_vector_sd(const std::function<void(const VoidCall<uint64_t>&)> &call_ones,
                         uint64_t size,
                         uint64_t num_set_bits);

    inline bit_vector_sd& operator=(const bit_vector_sd &other);
    inline bit_vector_sd& operator=(bit_vector_sd&& other) noexcept;

    inline std::unique_ptr<bit_vector> copy() const override;

    inline uint64_t rank1(uint64_t id) const override;
    inline uint64_t conditional_rank1(uint64_t id) const override;
    inline uint64_t select0(uint64_t id) const override;
    inline uint64_t select1(uint64_t id) const override;

    inline uint64_t next1(uint64_t id) const override;
    inline uint64_t prev1(uint64_t id) const override;

    inline bool operator[](uint64_t id) const override;
    inline uint64_t get_int(uint64_t id, uint32_t width) const override;

    inline bool load(std::istream &in) override;
    inline void serialize(std::ostream &out) const override;

    inline uint64_t size() const override { return vector_.size(); }

    inline bool is_inverted() const { return inverted_; }

    inline void call_ones_in_range(uint64_t begin, uint64_t end,
                                   const VoidCall<uint64_t> &callback) const override;

    inline void add_to(sdsl::bit_vector *other) const override;

    inline sdsl::bit_vector to_vector() const override;

    inline const sdsl::sd_vector<>& data() const { return vector_; }

    /**
     * Predict space taken by the vector with given its parameters in bits.
     */
    static uint64_t predict_size(uint64_t size, uint64_t num_set_bits) {
        num_set_bits = std::min(num_set_bits, size - num_set_bits);
        return sizeof(bit_vector_sd) * 8 + footprint_sd_vector(size, num_set_bits);
    }

  private:
    bool inverted_;
    sdsl::sd_vector<> vector_;
    sdsl::sd_vector<>::rank_1_type rk1_;
    sdsl::sd_vector<>::select_1_type slct1_;
    sdsl::sd_vector<>::select_0_type slct0_;
};


bit_vector_sd::bit_vector_sd(uint64_t size, bool value)
      : inverted_(value && size) {
    sdsl::sd_vector_builder builder(size, 0);
    vector_ = decltype(vector_)(builder);
    rk1_ = decltype(rk1_)(&vector_);
    slct1_ = decltype(slct1_)(&vector_);
    slct0_ = decltype(slct0_)(&vector_);
}

bit_vector_sd::bit_vector_sd(const sdsl::bit_vector &vector)
      : bit_vector_sd(vector, sdsl::util::cnt_one_bits(vector)) {}

bit_vector_sd::bit_vector_sd(const sdsl::bit_vector &vector, uint64_t num_set_bits) {
    assert(num_set_bits == sdsl::util::cnt_one_bits(vector));

    // check if it needs to be inverted
    if (num_set_bits <= vector.size() / 2) {
        // vector is sparse, no need to invert
        vector_ = sdsl::sd_vector<>(vector);
        inverted_ = false;
    } else {
        // invert
        sdsl::sd_vector_builder builder(vector.size(), vector.size() - num_set_bits);
        ::call_zeros(vector, [&builder](uint64_t i) { builder.set(i); });
        vector_ = sdsl::sd_vector<>(builder);
        inverted_ = true;
    }
    slct0_ = decltype(slct0_)(&vector_);
    slct1_ = decltype(slct1_)(&vector_);
    rk1_ = decltype(rk1_)(&vector_);
}

bit_vector_sd::bit_vector_sd(const bit_vector_sd &other) {
    *this = other;
}

bit_vector_sd::bit_vector_sd(bit_vector_sd&& other) {
    *this = std::move(other);
}

bit_vector_sd
::bit_vector_sd(const std::function<void(const VoidCall<uint64_t>&)> &call_ones,
                uint64_t size,
                uint64_t num_set_bits)
      : inverted_(num_set_bits > size / 2) {
    sdsl::sd_vector_builder builder(size, !inverted_
                                            ? num_set_bits
                                            : size - num_set_bits);
    if (inverted_) {
        uint64_t last_pos = 0;
        call_ones([&](uint64_t pos) {
            while (last_pos < pos) {
                builder.set(last_pos++);
            }
            ++last_pos;
        });
        while (last_pos < size) {
            builder.set(last_pos++);
        }
    } else {
        call_ones([&](uint64_t pos) { builder.set(pos); });
    }

    assert(builder.items() == builder.capacity());

    vector_ = sdsl::sd_vector<>(builder);
    slct0_ = decltype(slct0_)(&vector_);
    slct1_ = decltype(slct1_)(&vector_);
    rk1_ = decltype(rk1_)(&vector_);
}

bit_vector_sd::bit_vector_sd(std::initializer_list<bool> init)
      : bit_vector_sd(sdsl::bit_vector(init)) {}

bit_vector_sd& bit_vector_sd::operator=(const bit_vector_sd &other) {
    inverted_ = other.inverted_;
    vector_ = other.vector_;
    rk1_ = other.rk1_;
    rk1_.set_vector(&vector_);
    slct1_ = other.slct1_;
    slct1_.set_vector(&vector_);
    slct0_ = other.slct0_;
    slct0_.set_vector(&vector_);
    return *this;
}

bit_vector_sd& bit_vector_sd::operator=(bit_vector_sd&& other) noexcept {
    inverted_ = other.inverted_;
    vector_ = std::move(other.vector_);
    rk1_ = std::move(other.rk1_);
    rk1_.set_vector(&vector_);
    slct1_ = std::move(other.slct1_);
    slct1_.set_vector(&vector_);
    slct0_ = std::move(other.slct0_);
    slct0_.set_vector(&vector_);
    return *this;
}

std::unique_ptr<bit_vector> bit_vector_sd::copy() const {
    return std::make_unique<bit_vector_sd>(*this);
}

uint64_t bit_vector_sd::rank1(uint64_t id) const {
    //the rank method in SDSL does not include id in the count
    id = id >= this->size() ? this->size() : id + 1;
    return !inverted_ ? rk1_(id)
                      : id - rk1_(id);
}

uint64_t bit_vector_sd::conditional_rank1(uint64_t i) const {
    assert(i <= vector_.size());
    uint64_t high_val = (i >> (vector_.wl));
    uint64_t sel_high = vector_.high_0_select(high_val + 1);
    uint64_t rank_low = sel_high - high_val; //
    if (!rank_low)
        return !inverted_ ? 0 : i + 1;
    uint64_t val_low = i & sdsl::bits::lo_set[vector_.wl];
    // now since rank_low > 0 => sel_high > 0
    do {
        if (!sel_high)
            return !inverted_ ? 0 : i + 1;
        --sel_high;
        --rank_low;
    } while (vector_.high[sel_high] && vector_.low[rank_low] > val_low);

    if (vector_.high[sel_high] && vector_.low[rank_low] == val_low) {
        return !inverted_ ? rank_low + 1 : 0;
    } else {
        return !inverted_ ? 0 : i - rank_low;
    }
}

uint64_t bit_vector_sd::select0(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= rank0(size() - 1));
    assert(num_set_bits() == rank1(size() - 1));

    return !inverted_ ? slct0_(id) : slct1_(id);
}

uint64_t bit_vector_sd::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));

    return !inverted_ ? slct1_(id) : slct0_(id);
}

uint64_t bit_vector_sd::next1(uint64_t pos) const {
    assert(pos < size());

    return ::next1(*this, pos, !is_inverted() ? 0 : MAX_ITER_BIT_VECTOR_SD);
}

uint64_t bit_vector_sd::prev1(uint64_t pos) const {
    assert(pos < size());

    return ::prev1(*this, pos, !is_inverted() ? 0 : MAX_ITER_BIT_VECTOR_SD);
}

bool bit_vector_sd::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id] != inverted_;
}

uint64_t bit_vector_sd::get_int(uint64_t id, uint32_t width) const {
    if (inverted_)
        return ~vector_.get_int(id, width) & sdsl::bits::lo_set[width];

    return vector_.get_int(id, width);
}

bool bit_vector_sd::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);
        inverted_ = in.get();
        if (!in.good())
            return false;
        rk1_ = decltype(rk1_)(&vector_);
        slct1_ = decltype(slct1_)(&vector_);
        slct0_ = decltype(slct0_)(&vector_);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_sd." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void bit_vector_sd::serialize(std::ostream &out) const {
    vector_.serialize(out);
    if (!out.put(inverted_).good())
        throw std::ofstream::failure("Error when dumping bit_vector_sd");
}

sdsl::bit_vector bit_vector_sd::to_vector() const {
    sdsl::bit_vector vector(size(), inverted_);
    uint64_t max_rank = rk1_(size());
    for (uint64_t i = 1; i <= max_rank; ++i) {
        vector[slct1_(i)] = !inverted_;
    }
    return vector;
}

void bit_vector_sd::call_ones_in_range(uint64_t begin, uint64_t end,
                                       const VoidCall<uint64_t> &callback) const {
    assert(begin <= end);
    assert(end <= size());

    bit_vector::call_ones_adaptive(begin, end, callback, 2);
}

void bit_vector_sd::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());

    bit_vector::add_to_adaptive(other, 2);
}

#endif // __BIT_VECTOR_SD_HPP__
