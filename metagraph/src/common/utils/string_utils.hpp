#ifndef __STRING_UTILS_HPP__
#define __STRING_UTILS_HPP__

#include <string>
#include <string_view>
#include <deque>
#include <vector>
#include <cstdint>
#include <optional>


namespace utils {

bool starts_with(const std::string &str, const std::string &prefix);

// well-formed UTF-8: no overlong forms, no surrogates, at most U+10FFFF. What a JSON or MGT
// string may carry; a label name (a FASTA header, a file name) need not be one
bool valid_utf8(std::string_view s);

bool ends_with(const std::string &str, const std::string &suffix);

std::optional<uint64_t> parse_abundance(const std::string &comment);

std::string remove_suffix(const std::string &str, const std::string &suffix);

template <typename... String>
std::string remove_suffix(const std::string &str, const std::string &suffix,
                                                  const String&... other_suffixes) {
    return remove_suffix(remove_suffix(str, suffix), other_suffixes...);
}

std::string make_suffix(const std::string &str, const std::string &suffix);

std::string join_strings(const std::vector<std::string> &strings,
                         const std::string &delimiter,
                         bool discard_empty_strings = false);

std::vector<std::string> split_string(const std::string &string,
                                      const std::string &delimiter,
                                      bool skip_empty_parts = true);

/**
 * Given a minimum number of splits,
 * generate a list of suffixes from the alphabet.
 */
std::deque<std::string> generate_strings(const std::string &alphabet,
                                         size_t length);

} // namespace utils

#endif // __STRING_UTILS_HPP__
