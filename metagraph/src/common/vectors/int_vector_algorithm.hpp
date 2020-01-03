#ifndef __INT_VECTOR_ALGORITHM_HPP__
#define __INT_VECTOR_ALGORITHM_HPP__

#include <functional>
#include <cassert>

#include <sdsl/int_vector.hpp>


sdsl::bit_vector to_sdsl(const std::vector<bool> &vector);

template <class Callback>
void call_ones(const sdsl::bit_vector &vector,
               uint64_t begin, uint64_t end,
               Callback callback) {
    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for (; i < end && i & 0x3F; ++i) {
        if (vector[i])
            callback(i);
    }
    uint64_t word;
    for (uint64_t j = i + 64; j <= end; j += 64) {
        word = vector.get_int(i);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for (; i < j; ++i) {
            if (vector[i])
                callback(i);
        }
    }
    for (; i < end; ++i) {
        if (vector[i])
            callback(i);
    }
}

template <class Callback>
void call_ones(const sdsl::bit_vector &vector, Callback callback) {
    call_ones(vector, 0, vector.size(), callback);
}

template <class Callback>
void call_zeros(const sdsl::bit_vector &vector,
                uint64_t begin, uint64_t end,
                Callback callback) {
    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for (; i < end && i & 0x3F; ++i) {
        if (!vector[i])
            callback(i);
    }
    uint64_t word;
    for (uint64_t j = i + 64; j <= end; j += 64) {
        word = ~vector.get_int(i);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for (; i < j; ++i) {
            if (!vector[i])
                callback(i);
        }
    }
    for (; i < end; ++i) {
        if (!vector[i])
            callback(i);
    }
}

template <class Callback>
void call_zeros(const sdsl::bit_vector &vector, Callback callback) {
    call_zeros(vector, 0, vector.size(), callback);
}

uint64_t count_ones(const sdsl::bit_vector &vector,
                    uint64_t begin, uint64_t end);

uint64_t inner_prod(const sdsl::bit_vector &first,
                    const sdsl::bit_vector &second);


// Call (uint64_t index, uint64_t value) for each non-zero value in [begin, end)
template <class Callback>
void call_nonzeros(const sdsl::int_vector<> &vector,
                   uint64_t begin, uint64_t end,
                   Callback callback) {
    if (begin >= end)
        return;

    uint64_t i = begin;
    switch (vector.width()) {
        case 64:
            std::for_each(reinterpret_cast<const uint64_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint64_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        case 32:
            std::for_each(reinterpret_cast<const uint32_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint32_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        case 16:
            std::for_each(reinterpret_cast<const uint16_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint16_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        case 8:
            std::for_each(reinterpret_cast<const uint8_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint8_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        default:
            // vector.width() is not a power of two
            assert(vector.width() < 64);

            size_t begin_word = begin * vector.width() / 64;
            size_t end_word = (end * vector.width() + 63) / 64;

            for (uint64_t w = begin_word, it = begin; w < end_word; ++w) {
                if (!vector.data()[w])
                    continue;

                it = std::max(it, w * 64 / vector.width());

                auto it_end = std::min(end,
                    ((w + 1) * 64 + vector.width() - 1) / vector.width());

                for (; it < it_end; ++it) {
                    if (vector[it])
                        callback(it, vector[it]);
                }
            }
    }
}

// Call (uint64_t index, uint64_t value) for each non-zero value in |vector|.
template <class Callback>
void call_nonzeros(const sdsl::int_vector<> &vector, Callback callback) {
    return call_nonzeros(vector, 0, vector.size(), callback);
}

#endif // __INT_VECTOR_ALGORITHM_HPP__
