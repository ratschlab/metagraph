#ifndef __STRING_UTILS_HPP__
#define __STRING_UTILS_HPP__

#include <string>
#include <deque>
#include <vector>


namespace utils {

bool ends_with(const std::string &str, const std::string &suffix);
bool starts_with(const std::string &str, const std::string &prefix);

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
