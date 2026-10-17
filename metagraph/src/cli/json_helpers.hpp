#ifndef __METAGRAPH_CLI_JSON_HELPERS_HPP__
#define __METAGRAPH_CLI_JSON_HELPERS_HPP__

/**
 * JSON helpers shared by the routes' request parsers and answer builders (/traverse, /resolve,
 * /pattern and the predicate language): the strict reader of a request object, small value
 * builders, and the digit and escape counts the text-size estimates use. Header-only.
 */

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <iomanip>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <json/json.h>


namespace mtg {
namespace cli {

/**
 * |v| as a non-negative integer, nullopt when it is not one: an integer, or a number without a
 * fraction below 2^64, that is not negative. The sign is read with asInt64() only where
 * isInt64() holds: on an integral value above INT64_MAX asInt64() throws Json::LogicError, which
 * is no route's refusal (it was answered as a 400 with jsoncpp's text and no field path), and
 * asUInt64() on a negative one throws too. Each caller keeps its own message and range.
 */
inline std::optional<uint64_t> non_negative_integer(const Json::Value &v) {
    if (!v.isIntegral() || (v.isInt64() && v.asInt64() < 0))
        return std::nullopt;
    return v.asUInt64();
}

/**
 * Strict access to one JSON object of a request: every field read is remembered, and finish()
 * refuses the first field nothing read ("unknown field"), after the known ones were checked, so
 * that a field the server does not know is never ignored. Paths in messages are |path| and
 * path(key) = |path|.key.
 *
 * |Refusal| makes the exception a refusal throws: Refusal()(message) returns it (each route
 * answers a malformed request with its own error type and body).
 */
template <class Refusal>
class StrictObject {
  public:
    StrictObject(const Json::Value &value, std::string path)
          : v_(value), path_(std::move(path)) {
        if (!v_.isObject())
            refuse(path_ + ": expected an object");
    }

    bool has(const std::string &key) {
        seen_.insert(key);
        return v_.isMember(key);
    }
    const Json::Value& raw(const std::string &key) {
        seen_.insert(key);
        return v_[key];
    }
    std::string path(const std::string &key) const { return path_ + "." + key; }

    void finish() const {
        for (const std::string &name : v_.getMemberNames()) {
            if (!seen_.count(name))
                refuse(path_ + ": unknown field '" + name + "'");
        }
    }

    // Typed fields: |def| when the field is omitted, a refusal when it has another type or
    // lies out of range

    std::string str(const std::string &key, const std::string &def) {
        if (!has(key))
            return def;
        if (!v_[key].isString())
            refuse(path(key) + ": expected a string");
        return v_[key].asString();
    }

    bool boolean(const std::string &key, bool def) {
        if (!has(key))
            return def;
        if (!v_[key].isBool())
            refuse(path(key) + ": expected a boolean");
        return v_[key].asBool();
    }

    uint64_t uint(const std::string &key, uint64_t def, uint64_t min = 0,
                  uint64_t max = std::numeric_limits<uint64_t>::max()) {
        if (!has(key))
            return def;
        const std::optional<uint64_t> value = non_negative_integer(v_[key]);
        if (!value)
            refuse(path(key) + ": expected a non-negative integer");
        const uint64_t x = *value;
        if (x < min || x > max)
            refuse(path(key) + ": out of range [" + std::to_string(min) + ", "
                       + std::to_string(max) + "]");
        return x;
    }

    double number(const std::string &key, double def, double min = 0,
                  double max = std::numeric_limits<double>::infinity()) {
        if (!has(key))
            return def;
        if (!v_[key].isNumeric())
            refuse(path(key) + ": expected a number");
        const double x = v_[key].asDouble();
        if (!(x >= min && x <= max))
            refuse(path(key) + ": out of range");
        return x;
    }

    std::vector<std::string> strings(const std::string &key) {
        std::vector<std::string> out;
        if (!has(key))
            return out;
        if (!v_[key].isArray())
            refuse(path(key) + ": expected an array of strings");
        for (const Json::Value &x : v_[key]) {
            if (!x.isString())
                refuse(path(key) + ": expected an array of strings");
            out.push_back(x.asString());
        }
        return out;
    }

    // one of |values| by name; the refusal lists the names in their order, '|'-separated
    template <class E>
    E enumeration(const std::string &key, E def,
                  const std::vector<std::pair<std::string, E>> &values) {
        if (!has(key))
            return def;
        const std::string s = str(key, "");
        std::string allowed;
        for (const auto &[name, value] : values) {
            if (name == s)
                return value;
            allowed += (allowed.empty() ? "" : "|") + name;
        }
        refuse(path(key) + ": expected one of " + allowed);
    }

  protected:
    [[noreturn]] static void refuse(const std::string &message) { throw Refusal()(message); }

  private:
    const Json::Value &v_;
    std::string path_;
    std::set<std::string> seen_;
};

inline Json::Value uint_json(uint64_t x) { return Json::Value(static_cast<Json::UInt64>(x)); }

// A number of milliseconds or bits as JSON: an integer when it is one (the flags are integers
// and a client compares them as written), else the double
inline Json::Value number_json(double x) {
    if (x >= 0 && x == std::floor(x) && x <= 9007199254740991.0)
        return uint_json(static_cast<uint64_t>(x));
    return Json::Value(x);
}

// A number of milliseconds in a message: 250.5, not std::to_string's 250.500000 nor a cast's 250
inline std::string ms_text(double x) {
    std::ostringstream out;
    out << std::setprecision(15) << x;
    return out.str();
}

// |s| as a JSON string, null when it is empty
inline Json::Value string_or_null(const std::string &s) {
    return s.empty() ? Json::Value() : Json::Value(s);
}

// {"reason": reason}
inline Json::Value reason_json(const std::string &reason) {
    Json::Value r;
    r["reason"] = reason;
    return r;
}

inline Json::Value strings_json(std::initializer_list<const char*> values) {
    Json::Value a(Json::arrayValue);
    for (const char *v : values) {
        a.append(v);
    }
    return a;
}

// Appends {field, requested, effective} to |clamped|, the list of request values lowered to a
// server cap (limits.clamped). The values keep the type of their field: an integer field as an
// integer, a time budget as the number it was given as.
inline void note_clamped(Json::Value *clamped, const char *field, Json::Value requested,
                         Json::Value effective) {
    Json::Value c;
    c["field"] = field;
    c["requested"] = std::move(requested);
    c["effective"] = std::move(effective);
    clamped->append(std::move(c));
}

// The decimal digits of |x| (1 for 0)
inline uint64_t decimal_digits(uint64_t x) {
    uint64_t digits = 1;
    while (x >= 10) {
        x /= 10;
        ++digits;
    }
    return digits;
}

// The length jsoncpp's writer gives |s| inside a JSON string, quotes excluded (its
// escaping with emitUTF8 off, which every writer here uses): '"' '\\' and \b \f \n \r \t
// take two bytes, any other control character six (\u00XX), a code point above U+007F six
// (twelve as a surrogate pair), and a malformed sequence the six of U+FFFD for the bytes
// the writer consumes with it. Exact for valid UTF-8 and for most bytes that are not, from
// above for the rest: a 4-byte sequence that decodes below U+10000 (lead 0xF0, the next byte's
// low six bits below 0x10) is one U+FFFD (6), counted as a pair (12). The one model of that
// escaping: /traverse's delivery costs and compact_json_size, /pattern's answer text
// (string_text_bytes) and its predicate echo (quoted_bytes) price names by it.
// GraphletCodec.JsonEscapedSizeIsTheWriters checks it against the writer itself.
inline uint64_t json_escaped_size(std::string_view s) {
    uint64_t n = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        if (c == '"' || c == '\\' || c == '\b' || c == '\f' || c == '\n' || c == '\r' || c == '\t') {
            n += 2;
        } else if (c < 0x20) {
            n += 6;
        } else if (c < 0x80) {
            n += 1;
        } else if (c >= 0xF8) {
            n += 6;     // no lead byte: one replacement character
        } else {
            // a lead byte claims up to three more bytes, as jsoncpp's utf8ToCodepoint does
            // (without checking that they are continuation bytes); a sequence cut by the
            // end of the string is one replacement character for its lead byte
            const size_t width = c < 0xE0 ? 2 : c < 0xF0 ? 3 : 4;
            if (s.size() - i < width) {
                n += 6;
                continue;
            }
            n += width == 4 ? 12 : 6;
            i += width - 1;
        }
    }
    return n;
}

} // namespace cli
} // namespace mtg

#endif // __METAGRAPH_CLI_JSON_HELPERS_HPP__
