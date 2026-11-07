#ifndef METAGRAPH_COMMON_UTILS_ATOMIC_HPP
#define METAGRAPH_COMMON_UTILS_ATOMIC_HPP

#include <atomic>
#include <cstdint>

namespace utils {

static_assert(std::atomic_ref<uint64_t>::required_alignment <= alignof(uint64_t),
              "Atomic bit-vector operations require naturally aligned words");

inline uint64_t atomic_load(const uint64_t &word, std::memory_order order) {
    // Some libc++ versions cannot instantiate atomic_ref<const T>::load.
    // This reference is only used for loading; the object is never modified.
    return std::atomic_ref<uint64_t>(const_cast<uint64_t&>(word)).load(order);
}

} // namespace utils

#endif // METAGRAPH_COMMON_UTILS_ATOMIC_HPP
