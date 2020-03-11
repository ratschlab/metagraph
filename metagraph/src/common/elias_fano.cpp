#include "elias_fano.hpp"
namespace mg {
namespace common {

void write_bits(uint8_t *data, size_t pos, uint8_t __attribute__((unused)) len, uint64_t value) {
    assert(uint32_t(len) < 56);
    assert(0 == (value & ~((uint64_t(1) << len) - 1)));
    unsigned char *const ptr = data + (pos / 8);
    uint64_t ptrv = folly::loadUnaligned<uint64_t>(ptr);
    ptrv |= value << (pos % 8);
    folly::storeUnaligned<uint64_t>(ptr, ptrv);
}

uint8_t get_num_lower_bits(size_t max_value, size_t size) {
    if (size == 0 || max_value < size) {
        return 0;
    }
    // Result that should be returned is "floor(log(upperBound / size))".
    // In order to avoid expensive division, we rely on
    // "floor(a) - floor(b) - 1 <= floor(a - b) <= floor(a) - floor(b)".
    // Assuming "candidate = floor(log(upperBound)) - floor(log(upperBound))",
    // then result is either "candidate - 1" or "candidate".
    auto candidate = folly::findLastSet(max_value) - folly::findLastSet(size);
    // NOTE: As size != 0, "candidate" is always < 64.
    return (size > (max_value >> candidate)) ? candidate - 1 : candidate;
}

uint64_t clear_high_bits(uint64_t value, uint32_t index) {
    assert(index < 64);
    return value & ((uint64_t(1) << index) - 1);
}

} // namespace common
} // namespace mg
