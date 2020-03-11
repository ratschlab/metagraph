#include "elias_fano.hpp"
namespace mg {
namespace common {

/** Writes #value (with len up to 56 bits) to #data starting at the #pos-th bit. */
void write_bits(uint8_t *data, size_t pos, uint8_t len, uint64_t value) {
    assert(uint32_t(len) < 56);
    assert(0 == (value & ~((uint64_t(1) << len) - 1)));
    unsigned char *const ptr = data + (pos / 8);
    uint64_t ptrv = folly::loadUnaligned<uint64_t>(ptr);
    ptrv |= value << (pos % 8);
    folly::storeUnaligned<uint64_t>(ptr, ptrv);
    len = 0; // TODO: avoid unused warnings in a smarter way
}

uint64_t bzhi(uint64_t value, uint32_t index) {
    assert(index < 64);
    return value & ((uint64_t(1) << index) - 1);
}

} // namespace common
} // namespace mg
