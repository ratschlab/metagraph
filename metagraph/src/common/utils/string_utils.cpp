#include "string_utils.hpp"

#include <cmath>
#include <regex>
#include <cassert>
#include <algorithm>


namespace utils {

size_t utf8_sequence_at(std::string_view s, size_t i) {
    const unsigned char c = s[i];
    if (c < 0x80)
        return 1;
    size_t n;
    uint32_t cp;
    if ((c & 0xE0) == 0xC0) { n = 1; cp = c & 0x1F; }
    else if ((c & 0xF0) == 0xE0) { n = 2; cp = c & 0x0F; }
    else if ((c & 0xF8) == 0xF0) { n = 3; cp = c & 0x07; }
    else return 0;
    if (i + n >= s.size())
        return 0;
    for (size_t j = 1; j <= n; ++j) {
        const unsigned char d = s[i + j];
        if ((d & 0xC0) != 0x80)
            return 0;
        cp = (cp << 6) | (d & 0x3F);
    }
    if ((n == 1 && cp < 0x80) || (n == 2 && cp < 0x800) || (n == 3 && cp < 0x10000)
            || cp > 0x10FFFF || (cp >= 0xD800 && cp <= 0xDFFF))
        return 0;
    return n + 1;
}

bool valid_utf8(std::string_view s) {
    for (size_t i = 0; i < s.size(); ) {
        const size_t n = utf8_sequence_at(s, i);
        if (!n)
            return false;
        i += n;
    }
    return true;
}

std::string to_hex(std::string_view s) {
    static constexpr char digits[] = "0123456789abcdef";
    std::string out;
    out.reserve(s.size() * 2);
    for (unsigned char c : s) {
        out.push_back(digits[c >> 4]);
        out.push_back(digits[c & 0xF]);
    }
    return out;
}

std::string clip_for_message(std::string_view text, size_t max_bytes) {
    if (text.size() <= max_bytes)
        return std::string(text);
    // never inside a UTF-8 sequence: back to its lead byte, at most 3 bytes back (a name that is
    // no UTF-8 is cut where it falls), so that a clipped name keeps at least max_bytes - 3
    size_t cut = max_bytes;
    const size_t lowest = max_bytes > 3 ? max_bytes - 3 : 0;
    while (cut > lowest && (static_cast<unsigned char>(text[cut]) & 0xC0) == 0x80) {
        --cut;
    }
    return std::string(text.substr(0, cut)) + "... (" + std::to_string(text.size()) + " bytes)";
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
    return str.ends_with(suffix)
            ? str.substr(0, str.size() - suffix.size())
            : str;
}

std::string make_suffix(const std::string &str, const std::string &suffix) {
    return str.ends_with(suffix) ? str : str + suffix;
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
