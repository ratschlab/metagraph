#include "string_utils.hpp"

#include <cmath>
#include <regex>
#include <cassert>
#include <algorithm>


namespace utils {

bool starts_with(const std::string &str, const std::string &prefix) {
    if (prefix.size() > str.size()) {
        return false;
    }
    return prefix == std::string_view(str).substr(0, prefix.size());
}

bool ends_with(const std::string &str, const std::string &suffix) {
    auto actual_suffix = str.substr(
        std::max(0, static_cast<int>(str.size())
                    - static_cast<int>(suffix.size()))
    );
    return actual_suffix == suffix;
}

bool valid_utf8(std::string_view s) {
    size_t i = 0;
    while (i < s.size()) {
        const unsigned char c = s[i];
        size_t n;
        uint32_t cp;
        if (c < 0x80) { i++; continue; }
        if ((c & 0xE0) == 0xC0) { n = 1; cp = c & 0x1F; }
        else if ((c & 0xF0) == 0xE0) { n = 2; cp = c & 0x0F; }
        else if ((c & 0xF8) == 0xF0) { n = 3; cp = c & 0x07; }
        else return false;
        if (i + n >= s.size())
            return false;
        for (size_t j = 1; j <= n; ++j) {
            const unsigned char d = s[i + j];
            if ((d & 0xC0) != 0x80)
                return false;
            cp = (cp << 6) | (d & 0x3F);
        }
        if ((n == 1 && cp < 0x80) || (n == 2 && cp < 0x800) || (n == 3 && cp < 0x10000)
                || cp > 0x10FFFF || (cp >= 0xD800 && cp <= 0xDFFF))
            return false;
        i += n + 1;
    }
    return true;
}

// Try to parse ka:f:[abundance] or km:f:[abundance] from header
// https://github.com/IndexThePlanet/Logan/blob/main/Unitigs.md
std::optional<uint64_t> parse_abundance(const std::string &comment) {
    std::smatch match;
    std::regex abundance_regex(R"((ka|km):f:([0-9.eE+-]+))");
    if (std::regex_search(comment, match, abundance_regex) && match.size() > 2) {
        return std::max<uint64_t>(1, std::llround(std::stod(match[2].str())));
    } else {
        return std::nullopt;
    }
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
