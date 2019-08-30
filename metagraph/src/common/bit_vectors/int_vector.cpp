#include "int_vector.hpp"


void call_nonzeros(const sdsl::int_vector<> &vector,
                   uint64_t begin, uint64_t end,
                   std::function<void(uint64_t /* index */,
                                      uint64_t /* value */)> callback) {
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

void call_nonzeros(const sdsl::int_vector<> &vector,
                   std::function<void(uint64_t /* index */,
                                      uint64_t /* value */)> callback) {
    return call_nonzeros(vector, 0, vector.size(), callback);
}

template <class Vector>
void insert_new_indexes(Vector &vector, bit_vector_dyn *new_indexes) {
    size_t curpos = vector.size() - 1;
    assert(new_indexes->size() - vector.size() == new_indexes->num_set_bits());
    vector.resize(new_indexes->size());
    size_t i = vector.size() - 1;

    while (true) {
        if ((*new_indexes)[i]) {
            vector[i] = 0;
        } else {
            assert(curpos < vector.size());
            vector[i] = vector[curpos];
            curpos--;
        }
        if (0 == i)
            break;
        i--;
    }
}


template void insert_new_indexes<sdsl::int_vector<>>(sdsl::int_vector<> &vector, bit_vector_dyn *new_indexes);
