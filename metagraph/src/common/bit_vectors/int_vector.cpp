#include "int_vector.hpp"

template <typename Int, class Callback>
void call_range(const uint64_t* data,
                uint64_t begin, uint64_t end,
                Callback callback) {
    auto end_it = reinterpret_cast<const Int*>(data) + end;

    for (auto it = reinterpret_cast<const Int*>(data) + begin; it < end_it; ++it, ++begin) {
        if (*it)
            callback(begin, *it);
    }
}

void call_nonzeros(const sdsl::int_vector<> &vector,
                   uint64_t begin, uint64_t end,
                   std::function<void(uint64_t, uint64_t)> callback) {
    if (begin >= end)
        return;

    size_t begin_word = (begin * vector.width()) >> 6;
    size_t end_word = ((end * vector.width() + 63) >> 6);
    assert((end_word << 6) >= end * vector.width());
    assert(begin_word < end_word);

    switch (vector.width()) {
        case 64: { call_range<uint64_t>(vector.data(), begin, end, callback); } break;
        case 32: { call_range<uint32_t>(vector.data(), begin, end, callback); } break;
        case 16: { call_range<uint16_t>(vector.data(), begin, end, callback); } break;
        case 8: { call_range<uint8_t>(vector.data(), begin, end, callback); } break;
        default: {
            // vector.width() is not a power of two
            auto it = begin;
            size_t lcm = std::lcm(64, vector.width());
            for (size_t i = begin_word; i < end_word; ++i) {
                if (vector.data()[i]) {
                    it = std::max(it, (i << 6) / vector.width());

                    // set end iterator to next common multiple of vector.width() and 64
                    auto end_it = ((it + 1) * vector.width() + lcm - 1)
                        / lcm * lcm / vector.width();
                    end_it = std::min(end, end_it);

                    for (; it < end_it; ++it) {
                        if (vector[it] != 0)
                            callback(it, vector[it]);
                    }

                    if (it == end)
                        break;

                    i = ((it * vector.width()) >> 6) - 1;
                }
            }
        }
    }
}

void call_nonzeros(const sdsl::int_vector<> &vector,
                   std::function<void(uint64_t, uint64_t)> callback) {
    return call_nonzeros(vector, 0, vector.size(), callback);
}
