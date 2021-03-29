#include "string_utils.hpp"

#include <cmath>
#include <cassert>
#include <algorithm>


namespace utils {

bool ends_with(const std::string &str, const std::string &suffix) {
    auto actual_suffix = str.substr(
        std::max(0, static_cast<int>(str.size())
                    - static_cast<int>(suffix.size()))
    );
    return actual_suffix == suffix;
}

bool starts_with(const std::string &str, const std::string &prefix) {
    if (prefix.size() > str.size()) {
        return false;
    }
    return prefix == str.substr(0, static_cast<int>(prefix.size()));
}

std::string remove_suffix(const std::string &str, const std::string &suffix) {
    return ends_with(str, suffix)
            ? str.substr(0, str.size() - suffix.size())
            : str;
}

std::string make_suffix(const std::string &str, const std::string &suffix) {
    return remove_suffix(str, suffix) + suffix;
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
                                      const std::string &delimiter,
                                      bool skip_empty_parts) {
    if (!string.size())
        return {};

    if (!delimiter.size())
        return { string, };

    std::vector<std::string> result;

    size_t current_pos = 0;
    size_t delimiter_pos;

    while ((delimiter_pos = string.find(delimiter, current_pos))
                                             != std::string::npos) {
        if (delimiter_pos > current_pos || !skip_empty_parts)
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

} // namespace utils
