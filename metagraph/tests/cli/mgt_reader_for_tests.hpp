#ifndef __TESTS_CLI_MGT_READER_FOR_TESTS_HPP__
#define __TESTS_CLI_MGT_READER_FOR_TESTS_HPP__

/**
 * The reader of MGT v1 tokens (DESIGN-traverse-graphlet.md §2) and the reference encoders of
 * RANGES and SETEXPR it checks against, for the codec's conformance tests: the golden vectors
 * (codec_vectors.tsv) and the test Reader. The server and the CLI only write MGT, through
 * encode_float, pct_escape, front_code, encode_kvalue and the graphlet writer (cli/traverse.hpp);
 * the Python library parses it.
 */

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "cli/traverse.hpp"


namespace mtg {
namespace cli {
namespace mgt {

class FormatError : public std::invalid_argument {
  public:
    using std::invalid_argument::invalid_argument;
};

namespace reader_detail {

[[noreturn]] inline void fail(const std::string &what, std::string_view token) {
    throw FormatError(what + ": '" + std::string(token) + "'");
}

inline bool is_digit(char c) { return c >= '0' && c <= '9'; }

// a canonical unsigned decimal: no sign, no leading zero, fits 64 bits
inline bool parse_uint(std::string_view s, uint64_t *x) {
    if (s.empty() || (s.size() > 1 && s[0] == '0'))
        return false;
    uint64_t v = 0;
    for (char c : s) {
        if (!is_digit(c))
            return false;
        const uint64_t d = c - '0';
        if (v > (std::numeric_limits<uint64_t>::max() - d) / 10)
            return false;
        v = v * 10 + d;
    }
    *x = v;
    return true;
}

inline void append_uint(std::string &out, uint64_t x) {
    char buf[24];
    auto res = std::to_chars(buf, buf + sizeof(buf), x);
    out.append(buf, res.ptr);
}

inline bool is_continuation(unsigned char c) { return (c & 0xC0) == 0x80; }

inline int hex_value(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;   // lowercase hex is not canonical
}

} // namespace reader_detail

// the float a canonical token denotes (encode_float's spelling, or "inf")
inline double decode_float(std::string_view token) {
    if (token == "inf")
        return std::numeric_limits<double>::infinity();
    // (0|[1-9][0-9]*)(\.[0-9]*[1-9])?
    size_t i = 0;
    if (token.empty() || !reader_detail::is_digit(token[0]))
        reader_detail::fail("not a canonical float", token);
    if (token[0] == '0') {
        i = 1;
    } else {
        while (i < token.size() && reader_detail::is_digit(token[i])) i++;
    }
    if (i < token.size()) {
        if (token[i] != '.' || i + 1 == token.size() || token.back() == '0')
            reader_detail::fail("not a canonical float", token);
        for (size_t j = i + 1; j < token.size(); ++j) {
            if (!reader_detail::is_digit(token[j]))
                reader_detail::fail("not a canonical float", token);
        }
    }
    // the C locale's strtod is correctly rounded; whether the token is THE spelling of
    // the double it denotes is decided by re-encoding it (rejects the exact expansion of
    // 0.1, 1e23 written by to_chars(fixed), an overflow to inf and an underflow to 0)
    const std::string s(token);
    const double x = std::strtod(s.c_str(), nullptr);
    if (encode_float(x) != token)
        reader_detail::fail("not the shortest spelling of its float", token);
    return x;
}


// strictly ascending ids; every maximal run of >= 2 written a-b; "." = empty
inline std::string encode_ranges(const std::vector<uint64_t> &ids) {
    std::string out;
    if (ids.empty())
        return ".";
    constexpr uint64_t kMax = std::numeric_limits<uint64_t>::max();
    for (size_t i = 0; i < ids.size(); ) {
        if (i && ids[i] <= ids[i - 1])
            throw std::logic_error("RANGES: ids must be strictly ascending");
        size_t j = i;
        // adjacency without overflow: nothing follows UINT64_MAX (MAX + 1 wraps to 0, which
        // made {MAX, 0} one descending run "MAX-0")
        while (j + 1 < ids.size() && ids[j] != kMax && ids[j + 1] == ids[j] + 1) j++;
        if (!out.empty())
            out.push_back(',');
        reader_detail::append_uint(out, ids[i]);
        if (j > i) {
            out.push_back('-');
            reader_detail::append_uint(out, ids[j]);
        }
        i = j + 1;
    }
    return out;
}

// |bound|: every id must be below it, checked BEFORE a run is expanded, so that a few
// bytes of a corrupt or hostile token ("0-3000000000") cannot allocate in proportion to
// the range (a reader passes the number of L records)
inline std::vector<uint64_t> decode_ranges(std::string_view token,
                                           std::optional<uint64_t> bound = std::nullopt) {
    std::vector<uint64_t> out;
    if (token == ".")
        return out;
    if (token.empty())
        reader_detail::fail("empty RANGES", token);
    size_t pos = 0;
    while (true) {
        size_t comma = token.find(',', pos);
        std::string_view item = token.substr(pos, comma == std::string_view::npos
                                                      ? std::string_view::npos : comma - pos);
        const size_t dash = item.find('-');
        uint64_t a, b;
        if (dash == std::string_view::npos) {
            if (!reader_detail::parse_uint(item, &a))
                reader_detail::fail("malformed RANGES", token);
            b = a;
        } else if (!reader_detail::parse_uint(item.substr(0, dash), &a)
                       || !reader_detail::parse_uint(item.substr(dash + 1), &b) || b <= a) {
            reader_detail::fail("malformed RANGES", token);
        }
        // ascending, and every maximal run collapsed: an item may not continue the last.
        // Without overflow: nothing follows UINT64_MAX (out.back() + 1 wrapped to 0, so
        // "18446744073709551615,1" was accepted as the descending {MAX, 1})
        if (!out.empty() && (out.back() == std::numeric_limits<uint64_t>::max()
                             || a <= out.back() + 1))
            reader_detail::fail("RANGES not ascending or not collapsed", token);
        if (bound && b >= *bound)
            reader_detail::fail("RANGES id beyond the id space", token);
        for (uint64_t x = a; ; ++x) {
            out.push_back(x);
            if (x == b) break;
        }
        if (comma == std::string_view::npos)
            break;
        pos = comma + 1;
    }
    return out;
}

// against |base| (null: none, explicit only); the shorter of explicit and delta, a tie
// explicit
inline std::string encode_setexpr(const std::vector<uint64_t> &ids,
                                  const std::vector<uint64_t> *base) {
    std::string explicit_form = encode_ranges(ids);
    if (!base)
        return explicit_form;
    std::vector<uint64_t> removed, added;
    std::set_difference(base->begin(), base->end(), ids.begin(), ids.end(),
                        std::back_inserter(removed));
    std::set_difference(ids.begin(), ids.end(), base->begin(), base->end(),
                        std::back_inserter(added));
    std::string delta = "!";
    if (!removed.empty())
        delta += encode_ranges(removed);
    if (!added.empty())
        delta += "+" + encode_ranges(added);
    return delta.size() < explicit_form.size() ? delta : explicit_form;
}

inline std::vector<uint64_t> decode_setexpr(std::string_view token,
                                            const std::vector<uint64_t> *base,
                                            std::optional<uint64_t> bound = std::nullopt) {
    if (token.empty() || token[0] != '!')
        return decode_ranges(token, bound);
    if (!base)
        reader_detail::fail("a delta SETEXPR needs a base", token);
    std::string_view rest = token.substr(1);
    const size_t plus = rest.find('+');
    std::string_view r = rest.substr(0, plus);
    std::string_view a = plus == std::string_view::npos ? std::string_view() : rest.substr(plus + 1);
    // an empty part is omitted, never written '.'
    if (r == "." || (plus != std::string_view::npos && (a.empty() || a == ".")))
        reader_detail::fail("malformed SETEXPR", token);
    std::vector<uint64_t> removed = r.empty() ? std::vector<uint64_t>() : decode_ranges(r, bound);
    std::vector<uint64_t> added = a.empty() ? std::vector<uint64_t>() : decode_ranges(a, bound);
    // a mismatch with the base means the reader derived another base than the writer
    if (!std::includes(base->begin(), base->end(), removed.begin(), removed.end()))
        reader_detail::fail("SETEXPR removes ids its base does not hold", token);
    std::vector<uint64_t> kept, out;
    std::set_difference(base->begin(), base->end(), removed.begin(), removed.end(),
                        std::back_inserter(kept));
    for (uint64_t x : added) {
        if (std::binary_search(base->begin(), base->end(), x))
            reader_detail::fail("SETEXPR adds ids its base holds", token);
    }
    std::set_union(kept.begin(), kept.end(), added.begin(), added.end(), std::back_inserter(out));
    return out;
}


// the inverse of pct_escape: %25 %0A %0D, a raw line break refused
inline std::string pct_unescape(std::string_view token) {
    std::string out;
    out.reserve(token.size());
    for (size_t i = 0; i < token.size(); ++i) {
        const char c = token[i];
        if (c == '\n' || c == '\r')
            reader_detail::fail("raw line break in a free-text field", token);
        if (c != '%') {
            out.push_back(c);
            continue;
        }
        std::string_view esc = token.substr(i + 1, 2);
        if (esc == "25") out.push_back('%');
        else if (esc == "0A") out.push_back('\n');
        else if (esc == "0D") out.push_back('\r');
        else reader_detail::fail("not one of %25 %0A %0D", token);
        i += 2;
    }
    return out;
}


// the full name of a front-coded L record after |previous| (front_code)
inline std::string front_decode(std::string_view previous, std::string_view token) {
    const size_t space = token.find(' ');
    uint64_t p;
    if (space == std::string_view::npos
            || !reader_detail::parse_uint(token.substr(0, space), &p))
        reader_detail::fail("malformed front-coded name", token);
    if (p > previous.size()
            || (p < previous.size() && reader_detail::is_continuation(previous[p])))
        reader_detail::fail("prefix beyond the previous name or inside a character", token);
    return std::string(previous.substr(0, p)) + pct_unescape(token.substr(space + 1));
}


// a typed K value (encode_kvalue)
inline KValue decode_kvalue(std::string_view token) {
    KValue v;
    if (token == "u") {
        v.type = KValue::UNLIMITED;
        return v;
    }
    if (token.size() < 2 || token[1] != ':')
        reader_detail::fail("malformed K value", token);
    std::string_view body = token.substr(2);
    switch (token[0]) {
        case 'i':
            v.type = KValue::INTEGER;
            if (!reader_detail::parse_uint(body, &v.integer))
                reader_detail::fail("malformed K integer", token);
            return v;
        case 'f':
            v.type = KValue::FLOAT;
            v.real = decode_float(body);
            return v;
        case 's':
            v.type = KValue::STRING;
            for (size_t i = 0; i < body.size(); ++i) {
                const unsigned char c = body[i];
                if (c != '%') {
                    if (c < 0x21 || c > 0x7E || c == ',')
                        reader_detail::fail("a K string byte that must be escaped", token);
                    v.string.push_back(c);
                    continue;
                }
                const bool whole = i + 2 < body.size();
                const int hi = whole ? reader_detail::hex_value(body[i + 1]) : -1;
                const int lo = whole ? reader_detail::hex_value(body[i + 2]) : -1;
                if (hi < 0 || lo < 0)
                    reader_detail::fail("malformed escape in a K string", token);
                const unsigned char b = hi * 16 + lo;
                // an escape is canonical only for a byte that cannot be written raw
                if (b >= 0x21 && b <= 0x7E && b != '%' && b != ',')
                    reader_detail::fail("escape of a byte that is written raw", token);
                v.string.push_back(b);
                i += 2;
            }
            if (!valid_utf8(v.string))
                reader_detail::fail("a K string that is not UTF-8", token);
            return v;
    }
    reader_detail::fail("unknown K value type", token);
}

} // namespace mgt
} // namespace cli
} // namespace mtg

#endif // __TESTS_CLI_MGT_READER_FOR_TESTS_HPP__
