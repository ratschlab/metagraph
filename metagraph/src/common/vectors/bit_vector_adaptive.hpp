#ifndef __BIT_VECTOR_ADAPTIVE_HPP__
#define __BIT_VECTOR_ADAPTIVE_HPP__

#include <cstdint>

#include "bit_vector_sdsl.hpp"
#include "bit_vector_sd.hpp"


class bit_vector_adaptive : public bit_vector {
    friend bit_vector;

  public:
    virtual ~bit_vector_adaptive() {}

    virtual uint64_t rank1(uint64_t id) const override final { return vector_->rank1(id); }
    virtual uint64_t select1(uint64_t id) const override final { return vector_->select1(id); }
    virtual uint64_t select0(uint64_t id) const override final { return vector_->select0(id); }
    virtual std::pair<bool, uint64_t> inverse_select(uint64_t id) const override final {
        return vector_->inverse_select(id);
    }
    virtual uint64_t conditional_rank1(uint64_t id) const override final {
        return vector_->conditional_rank1(id);
    }

    virtual uint64_t next1(uint64_t id) const override final { return vector_->next1(id); }
    virtual uint64_t prev1(uint64_t id) const override final { return vector_->prev1(id); }

    virtual bool operator[](uint64_t id) const override final { return (*vector_)[id]; }
    virtual uint64_t get_int(uint64_t id, uint32_t width) const override final {
        return vector_->get_int(id, width);
    }

    inline bool load(std::istream &in) override final;
    inline void serialize(std::ostream &out) const override final;

    virtual uint64_t size() const override final { return vector_->size(); }

    virtual sdsl::bit_vector to_vector() const override final { return vector_->to_vector(); }

    virtual void add_to(sdsl::bit_vector *other) const override final { vector_->add_to(other); }

    virtual void call_ones_in_range(uint64_t begin, uint64_t end,
                                    const VoidCall<uint64_t> &callback) const override final {
        vector_->call_ones_in_range(begin, end, callback);
    }

    enum VectorCode {
    // FUI: don't change order of the variables!
    // Add new ones to the end, if any.
    // Otherwise, serialized vectors will be not loadable
        RRR_VECTOR = 0,
        SD_VECTOR,
        STAT_VECTOR,
        IL4096_VECTOR
    };

    typedef VectorCode (*DefineRepresentation)(uint64_t /* size */,
                                               uint64_t /* num_set_bits */);

    inline VectorCode representation_tag() const;

    const bit_vector& data() const { return *vector_; }

  protected:
    bit_vector_adaptive() {}

    inline explicit bit_vector_adaptive(const bit_vector_adaptive &other);
    bit_vector_adaptive(bit_vector_adaptive&& other) = default;
    inline bit_vector_adaptive& operator=(const bit_vector_adaptive &other);
    bit_vector_adaptive& operator=(bit_vector_adaptive&& other) = default;

    std::unique_ptr<bit_vector> vector_;
};


bit_vector_adaptive::bit_vector_adaptive(const bit_vector_adaptive &other)
      : vector_(other.vector_->copy()) {}

bit_vector_adaptive&
bit_vector_adaptive::operator=(const bit_vector_adaptive &other) {
    vector_ = other.vector_->copy();
    return *this;
}

bit_vector_adaptive::VectorCode
bit_vector_adaptive::representation_tag() const {
    if (dynamic_cast<const bit_vector_sd*>(vector_.get())) {
        return VectorCode::SD_VECTOR;

    } else if (dynamic_cast<const bit_vector_rrr<>*>(vector_.get())) {
        return VectorCode::RRR_VECTOR;

    } else if (dynamic_cast<const bit_vector_stat*>(vector_.get())) {
        return VectorCode::STAT_VECTOR;

    } else if (dynamic_cast<const bit_vector_il<4096>*>(vector_.get())) {
        return VectorCode::IL4096_VECTOR;

    } else {
        throw std::runtime_error("Unsupported type");
    }
}

bool bit_vector_adaptive::load(std::istream &in) {
    switch (load_number(in)) {
        case VectorCode::SD_VECTOR:
            vector_.reset(new bit_vector_sd());
            break;
        case VectorCode::RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>());
            break;
        case VectorCode::STAT_VECTOR:
            vector_.reset(new bit_vector_stat());
            break;
        case VectorCode::IL4096_VECTOR:
            vector_.reset(new bit_vector_il<4096>());
            break;
        default:
            return false;
    }
    return vector_->load(in);
}

void bit_vector_adaptive::serialize(std::ostream &out) const {
    serialize_number(out, representation_tag());
    vector_->serialize(out);
}


/**
 * static hybrid vector:
 *      the internal representation is defined in constructor.
 */
template <bit_vector_adaptive::DefineRepresentation optimal_representation>
class bit_vector_adaptive_stat : public bit_vector_adaptive {
    friend bit_vector;

  public:
    inline bit_vector_adaptive_stat(uint64_t size = 0, bool value = false);

    inline explicit bit_vector_adaptive_stat(const bit_vector &vector);
    inline explicit bit_vector_adaptive_stat(const sdsl::bit_vector &vector);

    inline bit_vector_adaptive_stat(bit_vector&& vector);
    inline bit_vector_adaptive_stat(sdsl::bit_vector&& vector);

    inline bit_vector_adaptive_stat(const VoidCall<const VoidCall<uint64_t>&> &call_ones,
                                    uint64_t size,
                                    uint64_t num_set_bits);

    bit_vector_adaptive_stat(std::initializer_list<bool> init)
      : bit_vector_adaptive_stat(sdsl::bit_vector(init)) {}

    inline std::unique_ptr<bit_vector> copy() const override final;

    /**
     * Predict space taken by the vector with given its parameters in bits.
     */
    inline static uint64_t predict_size(uint64_t size, uint64_t num_set_bits);
};


template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(uint64_t size, bool value) {
    switch (optimal_representation(size, size * static_cast<uint64_t>(value))) {
        case SD_VECTOR:
            vector_.reset(new bit_vector_sd(size, value));
            break;
        case RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>(size, value));
            break;
        case STAT_VECTOR:
            vector_.reset(new bit_vector_stat(size, value));
            break;
        case IL4096_VECTOR:
            vector_.reset(new bit_vector_il<4096>(size, value));
            break;
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(const bit_vector &vector) {
    switch (optimal_representation(vector.size(), vector.num_set_bits())) {
        case SD_VECTOR:
            vector_.reset(new bit_vector_sd(vector.copy_to<bit_vector_sd>()));
            break;
        case RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>(vector.copy_to<bit_vector_rrr<>>()));
            break;
        case STAT_VECTOR:
            vector_.reset(new bit_vector_stat(vector.copy_to<bit_vector_stat>()));
            break;
        case IL4096_VECTOR:
            vector_.reset(new bit_vector_il<4096>(vector.copy_to<bit_vector_il<4096>>()));
            break;
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(const sdsl::bit_vector &vector) {
    uint64_t num_set_bits = sdsl::util::cnt_one_bits(vector);

    switch (optimal_representation(vector.size(), num_set_bits)) {
        case SD_VECTOR:
            vector_.reset(new bit_vector_sd(vector, num_set_bits));
            break;
        case RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>(vector));
            break;
        case STAT_VECTOR:
            vector_.reset(new bit_vector_stat(vector));
            break;
        case IL4096_VECTOR:
            vector_.reset(new bit_vector_il<4096>(vector));
            break;
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(bit_vector&& vector) {
    switch (optimal_representation(vector.size(), vector.num_set_bits())) {
        case SD_VECTOR:
            vector_.reset(new bit_vector_sd(vector.convert_to<bit_vector_sd>()));
            break;
        case RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>(vector.convert_to<bit_vector_rrr<>>()));
            break;
        case STAT_VECTOR:
            vector_.reset(new bit_vector_stat(vector.convert_to<bit_vector_stat>()));
            break;
        case IL4096_VECTOR:
            vector_.reset(new bit_vector_il<4096>(vector.convert_to<bit_vector_il<4096>>()));
            break;
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(sdsl::bit_vector&& vector) {
    uint64_t num_set_bits = sdsl::util::cnt_one_bits(vector);

    switch (optimal_representation(vector.size(), num_set_bits)) {
        case SD_VECTOR:
            vector_.reset(new bit_vector_sd(std::move(vector), num_set_bits));
            break;
        case RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>(std::move(vector)));
            break;
        case STAT_VECTOR:
            vector_.reset(new bit_vector_stat(std::move(vector)));
            break;
        case IL4096_VECTOR:
            vector_.reset(new bit_vector_il<4096>(std::move(vector)));
            break;
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(const VoidCall<const VoidCall<uint64_t>&> &call_ones,
                           uint64_t size,
                           uint64_t num_set_bits) {
    switch (optimal_representation(size, num_set_bits)) {
        case SD_VECTOR: {
            vector_.reset(new bit_vector_sd(call_ones, size, num_set_bits));
            break;
        }
        case RRR_VECTOR: {
            sdsl::bit_vector vector(size, false);
            call_ones([&](uint64_t i) { vector[i] = true; });
            vector_.reset(new bit_vector_rrr<>(std::move(vector)));
            break;
        }
        case STAT_VECTOR: {
            sdsl::bit_vector vector(size, false);
            call_ones([&](uint64_t i) { vector[i] = true; });
            vector_.reset(new bit_vector_stat(std::move(vector)));
            break;
        }
        case IL4096_VECTOR: {
            sdsl::bit_vector vector(size, false);
            call_ones([&](uint64_t i) { vector[i] = true; });
            vector_.reset(new bit_vector_il<4096>(std::move(vector)));
            break;
        }
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
std::unique_ptr<bit_vector>
bit_vector_adaptive_stat<optimal_representation>::copy() const {
    auto copy = std::make_unique<bit_vector_adaptive_stat>();
    copy->vector_ = vector_->copy();
    return copy;
}


/**
 * hybrid vector: the smallest representation
 * combines:
 *    - bit_vector_sd
 *    - bit_vector_rrr
 */
inline bit_vector_adaptive::VectorCode
smallest_representation(uint64_t size, uint64_t num_set_bits) {
    assert(num_set_bits <= size);

    return bit_vector_sd::predict_size(size, num_set_bits)
            < bit_vector_rrr<>::predict_size(size, num_set_bits)
                  ? bit_vector_adaptive::VectorCode::SD_VECTOR
                  : bit_vector_adaptive::VectorCode::RRR_VECTOR;
}

typedef bit_vector_adaptive_stat<smallest_representation> bit_vector_small;

/**
 * hybrid vector: a good tradeoff between the speed and size
 * combines:
 *    - bit_vector_sd
 *    - bit_vector_stat
 */
inline bit_vector_adaptive::VectorCode
smart_representation(uint64_t size, uint64_t num_set_bits) {
    assert(num_set_bits <= size);

    return bit_vector_sd::predict_size(size, num_set_bits)
            < bit_vector_stat::predict_size(size, num_set_bits)
                  ? bit_vector_adaptive::VectorCode::SD_VECTOR
                  : bit_vector_adaptive::VectorCode::STAT_VECTOR;
}

typedef bit_vector_adaptive_stat<smart_representation> bit_vector_smart;

/**
 * hybrid vector: the smallest bit vector with fast rank support
 * combines:
 *    - bit_vector_sd
 *    - bit_vector_rrr<63>
 *    - bit_vector_il<4096>
 */
inline bit_vector_adaptive::VectorCode
smallrank_representation(uint64_t size, uint64_t num_set_bits) {
    assert(num_set_bits <= size);

    // use bit_vector_il<4096> if density is between (0.3, 0.7)
    if (std::min(num_set_bits, size - num_set_bits) > size * 0.3)
        return bit_vector_adaptive::VectorCode::IL4096_VECTOR;

    return bit_vector_sd::predict_size(size, num_set_bits)
            < bit_vector_rrr<>::predict_size(size, num_set_bits)
                  ? bit_vector_adaptive::VectorCode::SD_VECTOR
                  : bit_vector_adaptive::VectorCode::RRR_VECTOR;
}

typedef bit_vector_adaptive_stat<smallrank_representation> bit_vector_smallrank;


template <bit_vector_adaptive::DefineRepresentation optimal_representation>
uint64_t
bit_vector_adaptive_stat<optimal_representation>
::predict_size(uint64_t size, uint64_t num_set_bits) {
    switch (optimal_representation(size, num_set_bits)) {
        case SD_VECTOR:
            return bit_vector_sd::predict_size(size, num_set_bits);
        case RRR_VECTOR:
            return bit_vector_rrr<>::predict_size(size, num_set_bits);
        case STAT_VECTOR:
            return bit_vector_stat::predict_size(size, num_set_bits);
        case IL4096_VECTOR:
            return bit_vector_il<4096>::predict_size(size, num_set_bits);
        default:
            assert(false);
            return 0;
    }
}

#endif // __BIT_VECTOR_ADAPTIVE_HPP__
