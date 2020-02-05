#include "string_utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>


namespace utils {

static std::random_device device;
static std::mt19937 generator(device());

bool ends_with(const std::string &str, const std::string &suffix) {
    auto actual_suffix = str.substr(
        std::max(0, static_cast<int>(str.size())
                    - static_cast<int>(suffix.size()))
    );
    return actual_suffix == suffix;
}

std::string remove_suffix(const std::string &str, const std::string &suffix) {
    return ends_with(str, suffix)
            ? str.substr(0, str.size() - suffix.size())
            : str;
}

std::string join_strings(const std::vector<std::string> &strings,
                         const std::string &delimiter,
                         bool discard_empty_strings) {
    auto it = std::find_if(strings.begin(), strings.end(),
        [&](const auto &str) { return !discard_empty_strings || !str.empty(); }
    );

    std::string result;

    for (; it != strings.end(); ++it) {
        if (it->size() || !discard_empty_strings) {
            result += *it;
            result += delimiter;
        }
    }
    // remove last appended delimiter
    if (result.size())
        result.resize(result.size() - delimiter.size());

    return result;
}

std::vector<std::string> split_string(const std::string &string,
                                      const std::string &delimiter) {
    if (!string.size())
        return {};

    if (!delimiter.size())
        return { string, };

    std::vector<std::string> result;

    size_t current_pos = 0;
    size_t delimiter_pos;

    while ((delimiter_pos = string.find(delimiter, current_pos))
                                             != std::string::npos) {
        if (delimiter_pos > current_pos)
            result.push_back(string.substr(current_pos, delimiter_pos - current_pos));
        current_pos = delimiter_pos + delimiter.size();
    }
    if (current_pos < string.size()) {
        result.push_back(string.substr(current_pos));
    }

    assert(result.size());
    return result;
}

/**
 * Given a minimum number of splits,
 * generate a list of suffixes from the alphabet.
 */
std::deque<std::string> generate_strings(const std::string &alphabet,
                                         size_t length) {

    std::deque<std::string> suffixes = { "" };
    while (suffixes[0].length() < length) {
        for (const char c : alphabet) {
            suffixes.push_back(c + suffixes[0]);
        }
        suffixes.pop_front();
    }
    assert(suffixes.size() == std::pow(alphabet.size(), length));
    return suffixes;
}

std::string random_string(size_t length) {
    auto get_random_char = []() -> char {
        static constexpr char charset[]
                = "0123456789"
                  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                  "abcdefghijklmnopqrstuvwxyz";
        static constexpr size_t max_index = (sizeof(charset) - 2);
        static std::uniform_int_distribution<int> distribution(0, max_index);
        return charset[distribution(generator)];
    };
    std::string str(length, 0);
    std::generate_n(str.begin(), length, get_random_char);
    return str;
}

} // namespace utils
