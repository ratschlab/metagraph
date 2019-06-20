#include "int_vector.hpp"

void call_nonzeros(const sdsl::int_vector<> &vector,
                   uint64_t begin, uint64_t end,
                   std::function<void(uint64_t, uint64_t)> callback) {
    if (begin >= end)
        return;

    size_t begin_word = (begin * vector.width()) >> 6;
    size_t end_word = ((end * vector.width() + 63) >> 6);
    assert((end_word << 6) >= end * vector.width());
    assert(begin_word < end_word);

    size_t lcm = std::lcm(64, vector.width());

    for (size_t i = begin_word; i < end_word; ++i) {
        // TODO: if the last word of vector.data() is not fully initialized
        //       this may incorrecly evaluate to true. Though this shouldn't be
        //       a problem since the rest of the code computes the correct bounds
        if (vector.data()[i]) {
            // iterate until end, or until the next common multiple of word
            // size and vector width
            auto it = std::max(begin, (i << 6) / vector.width());
            auto end_it = (((it + 1) * vector.width() + lcm - 1) / lcm) * lcm / vector.width();
            i = ((end_it * vector.width()) >> 6) - 1;

            end_it = std::min(end_it, end);
            for (; it < end_it; ++it) {
                if (vector[it] != 0)
                    callback(it, vector[it]);
            }
        }
    }
}

void call_nonzeros(const sdsl::int_vector<> &vector,
                   std::function<void(uint64_t, uint64_t)> callback) {
    return call_nonzeros(vector, 0, vector.size(), callback);
}
