#ifndef __BIT_VECTOR_DYN_HPP__
#define __BIT_VECTOR_DYN_HPP__

#include <cstdint>

#include <dynamic/dynamic.hpp>

#include "vector_algorithm.hpp"
#include "bit_vector.hpp"


class bit_vector_dyn : public bit_vector {
    // sequential access for bit_vector_dyn is 18 times faster than select
    static constexpr size_t SEQ_ACCESS_VS_SELECT_FACTOR_DYN = 18;
    static constexpr size_t MAX_ITER_BIT_VECTOR_DYN = 10;

  public:
    inline explicit bit_vector_dyn(uint64_t size = 0, bool value = 0);
    inline explicit bit_vector_dyn(const sdsl::bit_vector &vector);
    inline bit_vector_dyn(std::initializer_list<bool> init);

    inline std::unique_ptr<bit_vector> copy() const override;

    inline uint64_t rank1(uint64_t id) const override;
    inline uint64_t select1(uint64_t id) const override;
    inline uint64_t select0(uint64_t id) const override;

    inline uint64_t next1(uint64_t id) const override;
    inline uint64_t prev1(uint64_t id) const override;

    inline bool operator[](uint64_t id) const override;
    inline uint64_t get_int(uint64_t id, uint32_t width) const override;

    inline void insert_bit(uint64_t id, bool val);
    inline void delete_bit(uint64_t id);
    inline void set(uint64_t id, bool val);

    inline bool load(std::istream &in) override;
    inline void serialize(std::ostream &out) const override;

    inline uint64_t size() const override { return vector_.size(); }

    inline void call_ones_in_range(uint64_t begin, uint64_t end,
                                   const VoidCall<uint64_t> &callback) const override;

    inline void add_to(sdsl::bit_vector *other) const override;

    inline sdsl::bit_vector to_vector() const override;

    /**
     * Predict space taken by the vector with given its parameters in bits.
     */
    static uint64_t predict_size(uint64_t size, uint64_t /*num_set_bits*/) {
        return size * 1.084;
    }

  private:
    dyn::suc_bv vector_;
};


bit_vector_dyn::bit_vector_dyn(uint64_t size, bool value) {
    uint64_t i;
    for (i = 0; i + 64 <= size; i += 64) {
        vector_.push_word(value ? ~uint64_t(0) : 0, 64);
    }
    if (i < size)
        vector_.push_word(value ? (uint64_t(1) << (size - i)) - 1 : 0, size - i);
}

bit_vector_dyn::bit_vector_dyn(const sdsl::bit_vector &v) {
    uint64_t i;
    for (i = 0; i + 64 <= v.size(); i += 64) {
        vector_.push_word(v.get_int(i, 64), 64);
    }
    if (i < v.size())
        vector_.push_word(v.get_int(i, v.size() - i), v.size() - i);
}

bit_vector_dyn::bit_vector_dyn(std::initializer_list<bool> init)
      : bit_vector_dyn(sdsl::bit_vector(init)) {}

std::unique_ptr<bit_vector> bit_vector_dyn::copy() const {
    return std::make_unique<bit_vector_dyn>(*this);
}

uint64_t bit_vector_dyn::rank1(uint64_t id) const {
    return vector_.rank1(std::min(id + 1, size()));
}

uint64_t bit_vector_dyn::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= rank1(size() - 1));
    return vector_.select1(id - 1);
}

uint64_t bit_vector_dyn::select0(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= rank0(size() - 1));
    return vector_.select0(id - 1);
}

uint64_t bit_vector_dyn::next1(uint64_t pos) const {
    assert(pos < size());

    return ::next1(*this, pos, num_set_bits() < size() / 3
                                ? 0 : MAX_ITER_BIT_VECTOR_DYN);
}

uint64_t bit_vector_dyn::prev1(uint64_t pos) const {
    assert(pos < size());

    return ::prev1(*this, pos, num_set_bits() < size() / 3
                                ? 0 : MAX_ITER_BIT_VECTOR_DYN);
}

void bit_vector_dyn::set(uint64_t id, bool val) {
    vector_.set(id, val);
}

bool bit_vector_dyn::operator[](uint64_t id) const {
    assert(id < size());
    return vector_.at(id);
}

uint64_t bit_vector_dyn::get_int(uint64_t id, uint32_t width) const {
    assert(id + width <= size());
    assert(width);
    // TODO: implement get_int in dyn::suc_bv and use it here instead of
    // querying one bit at a time
    uint64_t word = 0;
    for (int64_t pos = id + width - 1; pos >= static_cast<int64_t>(id); --pos) {
        word = (word << 1) + vector_.at(pos);
    }
    return word;
}

void bit_vector_dyn::insert_bit(uint64_t id, bool val) {
    assert(id <= size());
    vector_.insert(id, val);
}

void bit_vector_dyn::delete_bit(uint64_t id) {
    assert(id < size());
    vector_.remove(id);
}

bool bit_vector_dyn::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        //TODO: catch reading errors
        vector_.load(in);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_sd." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void bit_vector_dyn::serialize(std::ostream &out) const {
    vector_.serialize(out);
}

void bit_vector_dyn::call_ones_in_range(uint64_t begin, uint64_t end,
                                        const VoidCall<uint64_t> &callback) const {
    assert(begin <= end);
    assert(end <= size());

    bit_vector::call_ones_adaptive(begin, end, callback, SEQ_ACCESS_VS_SELECT_FACTOR_DYN);
}

sdsl::bit_vector bit_vector_dyn::to_vector() const {
    return bit_vector::to_vector_adaptive(SEQ_ACCESS_VS_SELECT_FACTOR_DYN);
}

void bit_vector_dyn::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());

    bit_vector::add_to_adaptive(other, SEQ_ACCESS_VS_SELECT_FACTOR_DYN);
}

#endif // __BIT_VECTOR_DYN_HPP__
