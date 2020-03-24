#include "bit_vector.hpp"

#include <type_traits>
#include <cassert>

#include "bit_vector_adaptive.hpp"
#include "bit_vector_stat.hpp"
#include "bit_vector_sdsl.hpp"
#include "bit_vector_dyn.hpp"
#include "bit_vector_sd.hpp"


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

    } else if (dynamic_cast<bit_vector_stat*>(this)) {
        // stat -> anything else
        return Vector(std::move(dynamic_cast<bit_vector_stat*>(this)->vector_));

    } else if (dynamic_cast<bit_vector_adaptive*>(this)) {
        // adaptive(x) -> anything
        return dynamic_cast<bit_vector_adaptive*>(this)->vector_->convert_to<Vector>();

    } else {
        // anything -> anything (slower: with full reconstruction)
        return Vector(to_vector());
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

    } else if (dynamic_cast<const bit_vector_stat*>(this)) {
        // copy stat -> anything else
        auto bv = dynamic_cast<const bit_vector_stat*>(this)->vector_;
        return Vector(std::move(bv));

    } else if (dynamic_cast<const bit_vector_adaptive*>(this)) {
        // copy adaptive(x) -> anything
        return dynamic_cast<const bit_vector_adaptive*>(this)->vector_->copy_to<Vector>();

    } else {
        // anything -> anything (slower: with full reconstruction)
        return Vector(to_vector());
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
