#include "algorithms.hpp"

#include <algorithm>


namespace utils {

std::vector<uint64_t> sample_indexes(uint64_t universe_size,
                                     uint64_t sample_size,
                                     std::mt19937 &gen) {
    if (!universe_size)
        return {};

    sample_size = std::min(universe_size, sample_size);

    std::vector<uint64_t> indexes;
    indexes.reserve(3 * sample_size);

    if (sample_size * 10 < universe_size) {
        std::uniform_int_distribution<uint64_t> dis(0, universe_size - 1);

        while (indexes.size() < sample_size) {
            indexes.clear();
            for (size_t i = 0; i < 1.5 * sample_size; ++i) {
                indexes.push_back(dis(gen));
            }
            std::sort(indexes.begin(), indexes.end());
            indexes.erase(std::unique(indexes.begin(), indexes.end()),
                          indexes.end());
        }
    } else {
        std::bernoulli_distribution dis(2.0 * sample_size / universe_size);

        while (indexes.size() < sample_size) {
            indexes.clear();
            for (size_t i = 0; i < universe_size; ++i) {
                if (dis(gen)) {
                    indexes.push_back(i);
                }
            }
        }
    }

    std::shuffle(indexes.begin(), indexes.end(), gen);

    return std::vector<uint64_t>(indexes.begin(), indexes.begin() + sample_size);
}

} // namespace utils
