#ifndef __STRING_UTILS_HPP__
#define __STRING_UTILS_HPP__

#include <string>
#include <string_view>
#include <deque>
#include <vector>
#include <cstdint>
#include <optional>


namespace utils {

// well-formed UTF-8: no overlong forms, no surrogates, at most U+10FFFF. What a JSON or MGT
// string may carry; a label name (a FASTA header, a file name) need not be one
bool valid_utf8(std::string_view s);

// The length in bytes of the well-formed UTF-8 sequence that begins at s[i] (1 for an ASCII
// byte, 2 to 4 for a multi-byte one), 0 when none does: a byte that is no lead byte (80-BF, C0,
// C1, F5-FF), a missing or wrong continuation byte, an overlong form, a surrogate (U+D800 to
// U+DFFF), a code point above U+10FFFF. valid_utf8's rule, one sequence at a time; i < s.size()
size_t utf8_sequence_at(std::string_view s, size_t i);

// the bytes of |s| as lowercase hex, two digits each: how a name that is not valid UTF-8 is
// reported, since its raw bytes could not be shown as text
std::string to_hex(std::string_view s);

// The most bytes of a client's string an error text quotes (clip_for_message)
constexpr size_t kMessageQuoteBytes = 1024;

// |text| as an error text quotes it: whole when at most |max_bytes| long, else its first
// |max_bytes| bytes (never cut inside a UTF-8 sequence) followed by "... (N bytes)", N its
// length. A message quoting a client's string (a label, a member name, a graph name) then stays
// small whatever the string: the server writes it into its answer, which escapes it again, and
// into its log
std::string clip_for_message(std::string_view text, size_t max_bytes = kMessageQuoteBytes);

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
