#include "bit_vector.hpp"

#include <type_traits>
#include <cassert>

#include "bit_vector_adaptive.hpp"
#include "bit_vector_sdsl.hpp"
#include "bit_vector_dyn.hpp"
#include "bit_vector_sd.hpp"
#include "vector_algorithm.hpp"


std::ostream& operator<<(std::ostream &os, const bit_vector &bv) {
    return os << bv.to_vector();
}

template <class Vector>
Vector bit_vector::convert_to() {
    static_assert(!std::is_same<Vector, bit_vector_smart>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_small>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_smallrank>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_adaptive>::value, "");

    if (dynamic_cast<Vector*>(this)) {
        // to the same type, no conversion
        return dynamic_cast<Vector&&>(*this);

    } else if (dynamic_cast<bit_vector_adaptive*>(this)) {
        // adaptive(x) -> anything
        return dynamic_cast<bit_vector_adaptive*>(this)->vector_->convert_to<Vector>();

    } else {
        uint64_t n_set_bits = num_set_bits();
        sdsl::bit_vector bv;
        if (auto *bv_stat = dynamic_cast<bit_vector_stat*>(this)) {
            // stat -> anything else
            bv = std::move(bv_stat->vector_);
        } else {
            // anything -> anything (slower: with full reconstruction)
            bv = to_vector();
        }

        if constexpr(std::is_same_v<Vector, bit_vector_sd>) {
            return Vector(std::move(bv), n_set_bits);
        } else {
            std::ignore = n_set_bits;
            return Vector(std::move(bv));
        }
    }
}
template bit_vector_dyn bit_vector::convert_to<bit_vector_dyn>();
template bit_vector_stat bit_vector::convert_to<bit_vector_stat>();
template bit_vector_sd bit_vector::convert_to<bit_vector_sd>();
template bit_vector_il<> bit_vector::convert_to<bit_vector_il<>>();
template bit_vector_hyb<> bit_vector::convert_to<bit_vector_hyb<>>();
template bit_vector_rrr<3> bit_vector::convert_to<bit_vector_rrr<3>>();
template bit_vector_rrr<8> bit_vector::convert_to<bit_vector_rrr<8>>();
template bit_vector_rrr<15> bit_vector::convert_to<bit_vector_rrr<15>>();
template bit_vector_rrr<31> bit_vector::convert_to<bit_vector_rrr<31>>();
template bit_vector_rrr<63> bit_vector::convert_to<bit_vector_rrr<63>>();
template bit_vector_rrr<127> bit_vector::convert_to<bit_vector_rrr<127>>();
template bit_vector_rrr<255> bit_vector::convert_to<bit_vector_rrr<255>>();
template sdsl::bit_vector bit_vector::convert_to<sdsl::bit_vector>();
template<> bit_vector_small bit_vector::convert_to() {
    return bit_vector_small(std::move(*this));
}
template<> bit_vector_smart bit_vector::convert_to() {
    return bit_vector_smart(std::move(*this));
}
template<> bit_vector_smallrank bit_vector::convert_to() {
    return bit_vector_smallrank(std::move(*this));
}

template <class Vector>
Vector bit_vector::copy_to() const {
    static_assert(!std::is_same<Vector, bit_vector_smart>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_small>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_smallrank>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_adaptive>::value, "");

    if (dynamic_cast<const Vector*>(this)) {
        // copy to the same type, no conversion
        return Vector(dynamic_cast<const Vector&>(*this));

    } else if (dynamic_cast<const bit_vector_adaptive*>(this)) {
        // copy adaptive(x) -> anything
        return dynamic_cast<const bit_vector_adaptive*>(this)->vector_->copy_to<Vector>();

    } else {
        sdsl::bit_vector bv;
        if (auto *bv_stat = dynamic_cast<const bit_vector_stat*>(this)) {
            // copy stat -> anything else
            bv = bv_stat->vector_;
        } else {
            // anything -> anything (slower: with full reconstruction)
            bv = to_vector();
        }

        if constexpr(std::is_same_v<Vector, bit_vector_sd>) {
            return Vector(std::move(bv), num_set_bits());
        } else {
            return Vector(std::move(bv));
        }
    }
}
template bit_vector_dyn bit_vector::copy_to<bit_vector_dyn>() const;
template bit_vector_stat bit_vector::copy_to<bit_vector_stat>() const;
template bit_vector_sd bit_vector::copy_to<bit_vector_sd>() const;
template bit_vector_il<> bit_vector::copy_to<bit_vector_il<>>() const;
template bit_vector_hyb<> bit_vector::copy_to<bit_vector_hyb<>>() const;
template bit_vector_rrr<3> bit_vector::copy_to<bit_vector_rrr<3>>() const;
template bit_vector_rrr<8> bit_vector::copy_to<bit_vector_rrr<8>>() const;
template bit_vector_rrr<15> bit_vector::copy_to<bit_vector_rrr<15>>() const;
template bit_vector_rrr<31> bit_vector::copy_to<bit_vector_rrr<31>>() const;
template bit_vector_rrr<63> bit_vector::copy_to<bit_vector_rrr<63>>() const;
template bit_vector_rrr<127> bit_vector::copy_to<bit_vector_rrr<127>>() const;
template bit_vector_rrr<255> bit_vector::copy_to<bit_vector_rrr<255>>() const;
template sdsl::bit_vector bit_vector::copy_to<sdsl::bit_vector>() const;
template<> bit_vector_small bit_vector::copy_to() const {
    return bit_vector_small(*this);
}
template<> bit_vector_smart bit_vector::copy_to() const {
    return bit_vector_smart(*this);
}
template<> bit_vector_smallrank bit_vector::copy_to() const {
    return bit_vector_smallrank(*this);
}

void bit_vector::call_ones_adaptive(uint64_t begin, uint64_t end,
                                    const VoidCall<uint64_t> &callback,
                                    double WORD_ACCESS_VS_SELECT_FACTOR) const {
    assert(begin <= end);
    assert(end <= size());

    // TODO: store num_set_bits() to avoid this call
    const uint64_t m = num_set_bits();
    if (m <= size() / WORD_ACCESS_VS_SELECT_FACTOR) {
        // sparse
        uint64_t num_ones = end ? rank1(end - 1) : 0;
        for (uint64_t r = begin ? rank1(begin - 1) + 1 : 1; r <= num_ones; ++r) {
            callback(select1(r));
        }
    } else if (size() - m <= size() / WORD_ACCESS_VS_SELECT_FACTOR) {
        // dense
        uint64_t one_pos = begin;
        uint64_t zero_pos = 0;
        uint64_t num_zeros = end ? rank0(end - 1) : 0;
        for (uint64_t r = begin ? rank0(begin - 1) + 1 : 1; r <= num_zeros; ++r) {
            zero_pos = select0(r);
            while (one_pos < zero_pos) {
                callback(one_pos++);
            }
            one_pos++;
        }
        while (one_pos < end) {
            callback(one_pos++);
        }
    } else {
        // moderate density
        ::call_ones(*this, begin, end, callback);
    }
}

sdsl::bit_vector
bit_vector::to_vector_adaptive(double WORD_ACCESS_VS_SELECT_FACTOR) const {
    sdsl::bit_vector result;

    if (num_set_bits() <= size() / WORD_ACCESS_VS_SELECT_FACTOR) {
        // sparse
        result = sdsl::bit_vector(size(), false);

        uint64_t num_ones = num_set_bits();
        for (uint64_t r = 1; r <= num_ones; ++r) {
            result[select1(r)] = true;
        }

    } else if ((size() - num_set_bits())
                <= size() / WORD_ACCESS_VS_SELECT_FACTOR) {
        // dense
        result = sdsl::bit_vector(size(), true);

        uint64_t num_zeros = size() - num_set_bits();
        for (uint64_t r = 1; r <= num_zeros; ++r) {
            result[select0(r)] = false;
        }

    } else {
        // moderate density
        result.resize(size());

        uint64_t i;
        const uint64_t end = size();
        uint64_t *data = result.data();
        for (i = 0; i + 64 <= end; i += 64, ++data) {
            *data = get_int(i, 64);
        }
        if (i < size()) {
            *data = get_int(i, size() - i);
        }
    }

    return result;
}

void bit_vector::add_to_adaptive(sdsl::bit_vector *other,
                                 double WORD_ACCESS_VS_SELECT_FACTOR) const {
    assert(other);
    assert(other->size() == size());

    if (std::min(num_set_bits(), size() - num_set_bits())
            <= size() / WORD_ACCESS_VS_SELECT_FACTOR) {
        // for very sparse or very dense vectors
        call_ones_adaptive(0, size(),
                           [other](auto i) { (*other)[i] = true; },
                           WORD_ACCESS_VS_SELECT_FACTOR);
    } else {
        // moderate density
        uint64_t i;
        const uint64_t end = size();
        uint64_t *data = other->data();
        for (i = 0; i + 64 <= end; i += 64, ++data) {
            *data |= get_int(i, 64);
        }
        if (i < size()) {
            *data |= get_int(i, size() - i);
        }
    }
}
