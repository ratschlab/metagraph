#include "elias_fano.hpp"

namespace mg{
namespace common {
template <class T, class Enable = void>
struct Unaligned;

/**
 * Representation of an unaligned value of a POD type.
 */
template <class T>
struct Unaligned<T, typename std::enable_if<std::is_pod<T>::value>::type> {
    Unaligned() = default; // uninitialized
    /* implicit */ Unaligned(T v) : value(v) {}
    T value;
} __attribute__((__packed__));

/**
 * Read an unaligned value of type T and return it.
 */
template <class T>
inline T load_unaligned(const void *p) {
    static_assert(sizeof(Unaligned<T>) == sizeof(T), "Invalid unaligned size");
    static_assert(alignof(Unaligned<T>) == 1, "Invalid alignment");
    return static_cast<const Unaligned<T> *>(p)->value;
}

/**
 * Inform the compiler that the argument can be assumed true. It is
 * undefined behavior if the argument is not actually true, so use
 * with care.
 *
 * Implemented as a function instead of a macro because
 * __builtin_assume does not evaluate its argument at runtime, so it
 * cannot be used with expressions that have side-effects.
 */
inline __attribute__((__always_inline__)) void assume(bool cond) {
#if defined(__clang__) // Must go first because Clang also defines __GNUC__.
    __builtin_assume(cond);
#elif defined(__GNUC__)
    if (!cond) {
        __builtin_unreachable();
    }
#endif
}

/**
 * Write an unaligned value of type T.
 */
template <class T>
inline void store_unaligned(void *p, T value) {
    static_assert(sizeof(Unaligned<T>) == sizeof(T), "Invalid unaligned size");
    static_assert(alignof(Unaligned<T>) == 1, "Invalid alignment");
    // Prior to C++14, the spec says that a placement new like this
    // is required to check that p is not nullptr, and to do nothing
    // if p is a nullptr. By assuming it's not a nullptr, we get a
    // nice loud segfault in optimized builds if p is nullptr, rather
    // than just silently doing nothing.
    assume(p != nullptr);
    new (p) Unaligned<T>(value);
}



}
}
