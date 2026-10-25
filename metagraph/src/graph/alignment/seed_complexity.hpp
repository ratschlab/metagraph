#pragma once

#include <cstdlib>
#include <memory>
#include <string_view>

#include <sdust.h>

namespace mtg::graph::align {

inline bool is_low_complexity(std::string_view sequence, int threshold = 20, int window = 64) {
#if ! _PROTEIN_GRAPH
    int count;
    std::unique_ptr<uint64_t, decltype(std::free)*> regions {
        sdust(nullptr, reinterpret_cast<const uint8_t*>(sequence.data()), sequence.size(),
              threshold, window, &count),
        std::free
    };
    return count > 0;
#else
    // DUST is a DNA filter; protein seed filtering is not implemented.
    return false;
#endif
}

} // namespace mtg::graph::align
