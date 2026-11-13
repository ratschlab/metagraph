#ifndef __METAGRAPH_CLI_JSON_HELPERS_HPP__
#define __METAGRAPH_CLI_JSON_HELPERS_HPP__

/**
 * JSON helpers shared by the routes' request parsers and answer builders (/traverse, /resolve,
 * /pattern and the predicate language): the strict reader of a request body and of a request
 * object, small value builders, and the digit and escape counts the text-size estimates use.
 * Header-only.
 */

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <json/json.h>

#include "common/threads/thread_usage.hpp"
#include "common/utils/string_utils.hpp"


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

// The deepest nesting of arrays and objects a request body may have (jsoncpp's stackLimit, its
// default): parse_strict_json refuses a deeper one
constexpr int kMaxBodyNesting = 1000;

/**
 * A body with more JSON values than the server reads (--max-request-json-values; GET
 * /capabilities states it as request_intake.max_request_json_values): thrown by
 * parse_strict_json before the text is parsed, so that jsoncpp never builds the tree of such a
 * body (about 100 bytes and more per value, also for a body refused afterwards). Each route
 * answers it alike (400 too_many_json_values, server_checks.hpp); the CLI, which sets no cap,
 * never throws it.
 */
class TooManyJsonValues : public std::invalid_argument {
  public:
    explicit TooManyJsonValues(uint64_t max_values)
          : std::invalid_argument("the request body holds more than "
                                  + std::to_string(max_values)
                                  + " JSON values (--max-request-json-values): it was not "
                                    "parsed and nothing of it ran; send fewer values (split "
                                    "the request)"),
            max_values_(max_values) {}
    uint64_t max_values() const { return max_values_; }

  private:
    uint64_t max_values_;
};

/**
 * jsoncpp's message of a text it refused, its quotes of the client's text clipped as every
 * error text clips a client's string (clip_for_message: whole up to 1,024 bytes of its own,
 * else its first 1,024 bytes and "... (N bytes)"). jsoncpp 1.9.3 quotes in two messages only:
 * "Duplicate key: '<name>'" (the name decoded, its closing quote the message's last) and
 * "'<token>' is not a number." (a raw number token, which holds no quote); it stops at the
 * first fault, so a message holds one of them at most
 */
inline std::string clip_reader_quotes(std::string errors) {
    auto clip = [&](size_t from, size_t to) {
        if (to > from && to - from > ::utils::kMessageQuoteBytes) {
            errors = errors.substr(0, from)
                    + ::utils::clip_for_message(std::string_view(errors).substr(from, to - from))
                    + errors.substr(to);
        }
    };
    static const std::string kDuplicate = "Duplicate key: '";
    if (const size_t at = errors.find(kDuplicate); at != std::string::npos) {
        const size_t to = errors.rfind('\'');
        if (to != std::string::npos)
            clip(at + kDuplicate.size(), to);
    }
    static const std::string kNotNumber = "' is not a number.";
    if (const size_t to = errors.rfind(kNotNumber); to != std::string::npos && to > 0) {
        const size_t open = errors.rfind('\'', to - 1);
        if (open != std::string::npos)
            clip(open + 1, to);
    }
    return errors;
}

/**
 * A request body as one RFC 8259 JSON text whose member names are unique: no comment, no
 * trailing comma, nothing after the value, no duplicated member name (jsoncpp would keep the
 * last value without a word: an id or a budget given twice would be replaced silently), at most
 * kMaxBodyNesting nested arrays and objects. Throws Refusal()(|prefix| + what is wrong)
 * otherwise. Not CharReaderBuilder::strictMode(): its strictRoot would refuse a root that is
 * not an object, which each route refuses in its own order. The reader of /pattern, /resolve,
 * /traverse and /traverse/cancel and of the pattern and traverse CLIs; /search and /align keep
 * theirs (parse_json_string), /search with its values counted (count_lenient_json_values).
 *
 * |max_values| (0: unlimited): the most JSON values the text may hold — the root, each array
 * element and each member's value (a member's name is none; an empty array or object holds
 * none) —, counted in the scan below before anything is parsed: one more throws
 * TooManyJsonValues at the value that passes it. Exact for every JSON text; a text that is
 * not JSON may be refused by either check. The scan refuses at the first fault it meets: a
 * comment, a raw control character inside a string (RFC 8259 asks it escaped), a string that
 * is not well-formed UTF-8 (::utils::utf8_sequence_at), a \u escape of a surrogate that is not
 * half of a pair, or the value past the count. A backslash that begins no JSON escape ends the
 * scan: jsoncpp refuses the text at that escape, or at a fault of its own before it.
 */
template <class Refusal>
Json::Value parse_strict_json(std::string_view content, const std::string &prefix,
                              uint64_t max_values = 0) {
    // jsoncpp skips a comment before an object's member name whatever allowComments says: a
    // '/' outside a string is no JSON at all, so the text is refused before it is parsed.
    // The same scan counts the values: one per ',' outside a string (the next element or
    // member) and one per '[' or '{' that does not close at once (its first), the root
    // counted before
    uint64_t values = 1;
    auto count_value = [&]() {
        if (++values > max_values && max_values)
            throw TooManyJsonValues(max_values);
    };
    bool in_string = false;
    // an array or object was just opened: its first value is counted at the next character
    // that is not whitespace, unless that character closes it
    bool opened = false;
    static constexpr char kHex[] = "0123456789abcdef";
    auto hex_of = [](unsigned char c) { return std::string { kHex[c >> 4], kHex[c & 0xF] }; };
    // the code unit of the four hex digits at |at|, -1 when there are not four
    auto code_unit_at = [&](size_t at) {
        if (at + 4 > content.size())
            return -1;
        int unit = 0;
        for (size_t j = at; j < at + 4; ++j) {
            const char d = content[j];
            const int v = d >= '0' && d <= '9' ? d - '0'
                        : d >= 'a' && d <= 'f' ? d - 'a' + 10
                        : d >= 'A' && d <= 'F' ? d - 'A' + 10 : -1;
            if (v < 0)
                return -1;
            unit = unit * 16 + v;
        }
        return unit;
    };
    for (size_t i = 0; i < content.size(); ++i) {
        const char c = content[i];
        if (in_string) {
            if (c == '\\') {
                // a \u escape of a surrogate is half of a pair: a high one (D800-DBFF) followed
                // at once by a low one (DC00-DFFF). jsoncpp refuses a lone high one, but decodes
                // a lone low one to bytes that are no UTF-8 and joins a high one with whatever
                // escape follows ("\ud800\u0041" becomes U+10041)
                const char escaped = i + 1 < content.size() ? content[i + 1] : '\0';
                const int unit = escaped == 'u' ? code_unit_at(i + 2) : -1;
                if (unit >= 0xD800 && unit <= 0xDFFF) {
                    const int low = unit <= 0xDBFF && i + 7 < content.size()
                                        && content[i + 6] == '\\' && content[i + 7] == 'u'
                            ? code_unit_at(i + 8) : -1;
                    if (low < 0xDC00 || low > 0xDFFF) {
                        throw Refusal()(prefix + "an unpaired surrogate escape (\\u"
                                        + std::string(content.substr(i + 2, 4))
                                        + ") inside a string at byte " + std::to_string(i)
                                        + ": a code point above U+FFFF is written as a high "
                                          "surrogate escape (\\ud800-\\udbff) followed by a "
                                          "low one (\\udc00-\\udfff)");
                    }
                    i += 11;    // both escapes
                } else if (unit >= 0) {
                    i += 5;     // the escape and its four digits
                } else if (escaped == '"' || escaped == '\\' || escaped == '/' || escaped == 'b'
                            || escaped == 'f' || escaped == 'n' || escaped == 'r'
                            || escaped == 't') {
                    ++i;        // the escaped character
                } else {
                    // no JSON escape (jsoncpp reads these and \u with four hex digits, nothing
                    // else): the text is no JSON from here on. The scan ends and leaves it to
                    // jsoncpp, which refuses it at this escape ("Bad escape sequence in string",
                    // "Bad unicode escape sequence in string: ...") or at a fault of its own
                    // before it, so that no later fault's text, or a byte the escape was meant
                    // to take, is named instead. jsoncpp builds no value past the escape: every
                    // value it builds was counted above
                    break;
                }
            } else if (c == '"') {
                in_string = false;
            } else if (static_cast<unsigned char>(c) < 0x20) {
                // RFC 8259 §7: a control character inside a string is escaped. jsoncpp keeps a
                // raw one, and an answer quoting the string writes it as \u00XX: 6 bytes for 1
                const std::string hex = hex_of(static_cast<unsigned char>(c));
                throw Refusal()(prefix + "a control character (byte 0x" + hex
                                + ") inside a string at byte " + std::to_string(i)
                                + ": RFC 8259 requires it escaped (\\u00" + hex + ")");
            } else if (static_cast<unsigned char>(c) >= 0x80) {
                // RFC 8259 §8.1: the text is UTF-8. jsoncpp keeps any byte, and its writer turns
                // one that begins no well-formed sequence into U+FFFD (6 bytes for 1 at 0xF8 to
                // 0xFF) or swallows up to three bytes after it, so that distinct names echo alike
                const size_t n = ::utils::utf8_sequence_at(content, i);
                if (!n) {
                    throw Refusal()(prefix + "invalid UTF-8 (byte 0x"
                                    + hex_of(static_cast<unsigned char>(c))
                                    + ") inside a string at byte " + std::to_string(i)
                                    + ": RFC 8259 requires UTF-8");
                }
                i += n - 1;
            }
            continue;
        }
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
            continue;
        if (opened && c != ']' && c != '}')
            count_value();
        opened = false;
        if (c == '"') {
            in_string = true;
        } else if (c == '/') {
            throw Refusal()(prefix + "a comment (or a '/' outside a string) at byte "
                            + std::to_string(i));
        } else if (c == '[' || c == '{') {
            opened = true;
        } else if (c == ',') {
            count_value();
        }
    }
    Json::CharReaderBuilder builder;
    builder["allowComments"] = false;
    builder["allowTrailingCommas"] = false;
    builder["failIfExtra"] = true;
    builder["rejectDupKeys"] = true;
    builder["stackLimit"] = kMaxBodyNesting;
    std::unique_ptr<Json::CharReader> reader { builder.newCharReader() };
    Json::Value json;
    std::string errors;
    bool parsed = false;
    try {
        parsed = reader->parse(content.data(), content.data() + content.size(), &json, &errors);
    } catch (const Json::Exception &e) {
        // jsoncpp throws rather than returns past its nesting limit: the client's to fix, as
        // every other fault of the text
        throw Refusal()(prefix + e.what() + " (more than " + std::to_string(kMaxBodyNesting)
                        + " nested arrays and objects)");
    }
    if (!parsed)
        throw Refusal()(prefix + clip_reader_quotes(std::move(errors)));
    return json;
}

/**
 * The JSON values that jsoncpp's default reader (CharReaderBuilder's defaults: the reader of
 * /search, parse_json_string) builds of |content| (the root, each array element and each
 * member's value, as parse_strict_json counts them), counted before the text is parsed: one
 * more than |max_values| (0: unlimited) throws TooManyJsonValues at the value that passes it.
 * Nothing else is checked: the text is the lenient reader's to accept or refuse afterwards, as
 * before. The count follows jsoncpp 1.9.3's OurReader token by token (readToken, readValue,
 * readArray, readObject with collectComments, allowComments, allowTrailingCommas, no single
 * quotes, failIfExtra false, stackLimit 1000), so that a text it accepts is counted at the
 * number of values it builds, its quirks included:
 * - a comment opens at a '/' or a '\'' (a '\'' is read as a '/') followed by a '*' or a '/';
 * - a trailing comma closes an array or an object (the bracket right after it);
 * - in an object, a member's value followed by a comment takes the next token, whatever it is,
 *   for the separator;
 * - nothing after the root value is read.
 * A member named twice counts each time (the reader builds both values, the tree keeps the
 * last). Where the reader stops at a fault of the text's structure (a token where none may
 * stand, nesting past the stack limit) the count stops too; a fault it finds only in a token's
 * content (a number, an escape, a comment it cannot attach) is not looked for, so that such a
 * text, which it refuses, may be counted higher than what it built before the fault, never
 * lower.
 */
class LenientJsonValues {
  public:
    static uint64_t count(std::string_view content, uint64_t max_values) {
        LenientJsonValues reader(content, max_values);
        // the root
        reader.add();
        reader.value(1);
        return reader.values_;
    }

  private:
    enum class Token { kObjectBegin, kObjectEnd, kArrayBegin, kArrayEnd, kString, kNumber,
                       kLiteral, kComment, kArraySeparator, kMemberSeparator, kEnd, kError };

    // jsoncpp's stackLimit (its default, as parse_json_string's): it throws past it
    static constexpr size_t kStackLimit = 1000;

    LenientJsonValues(std::string_view content, uint64_t max_values)
          : s_(content), max_(max_values) {}

    void add() {
        if (++values_ > max_ && max_)
            throw TooManyJsonValues(max_);
    }

    void skip_spaces() {
        while (at_ < s_.size()
                && (s_[at_] == ' ' || s_[at_] == '\t' || s_[at_] == '\r' || s_[at_] == '\n')) {
            ++at_;
        }
    }

    // the next character, '\0' at the end (not advancing there), as jsoncpp's getNextChar
    char next_char() { return at_ < s_.size() ? s_[at_++] : '\0'; }

    bool match(std::string_view pattern) {
        if (s_.substr(at_, pattern.size()) != pattern)
            return false;
        at_ += pattern.size();
        return true;
    }

    // readNumber: digits, a fraction and an exponent, each optional; a '-' or a '+' before
    // an 'I' is no number (special floats are off)
    Token number(bool sign) {
        if (sign && at_ < s_.size() && s_[at_] == 'I') {
            ++at_;
            return Token::kError;
        }
        size_t p = at_;
        auto next = [&]() {
            at_ = p;
            return p < s_.size() ? s_[p++] : '\0';
        };
        auto digit = [](char c) { return c >= '0' && c <= '9'; };
        char c = '0';
        while (digit(c)) {
            c = next();
        }
        if (c == '.') {
            c = next();
            while (digit(c)) {
                c = next();
            }
        }
        if (c == 'e' || c == 'E') {
            c = next();
            if (c == '+' || c == '-')
                c = next();
            while (digit(c)) {
                c = next();
            }
        }
        return Token::kNumber;
    }

    // readComment, after its opening '/' or '\''
    Token comment() {
        const char c = next_char();
        if (c == '*') {
            while (at_ + 1 < s_.size()) {
                const char d = s_[at_++];
                if (d == '*' && s_[at_] == '/')
                    break;
            }
            return next_char() == '/' ? Token::kComment : Token::kError;
        }
        if (c == '/') {
            while (at_ < s_.size()) {
                const char d = s_[at_++];
                if (d == '\n')
                    break;
                if (d == '\r') {
                    if (at_ < s_.size() && s_[at_] == '\n')
                        ++at_;
                    break;
                }
            }
            return Token::kComment;
        }
        return Token::kError;
    }

    Token token() {
        skip_spaces();
        const char c = next_char();
        switch (c) {
            case '{': return Token::kObjectBegin;
            case '}': return Token::kObjectEnd;
            case '[': return Token::kArrayBegin;
            case ']': return Token::kArrayEnd;
            case '"': {
                char d = 0;
                while (at_ < s_.size()) {
                    d = s_[at_++];
                    if (d == '\\')
                        next_char();
                    else if (d == '"')
                        break;
                }
                return d == '"' ? Token::kString : Token::kError;
            }
            case '\'':
            case '/': return comment();
            case '0': case '1': case '2': case '3': case '4':
            case '5': case '6': case '7': case '8': case '9': return number(false);
            case '-':
            case '+': return number(true);
            case 't': return match("rue") ? Token::kLiteral : Token::kError;
            case 'f': return match("alse") ? Token::kLiteral : Token::kError;
            case 'n': return match("ull") ? Token::kLiteral : Token::kError;
            case ',': return Token::kArraySeparator;
            case ':': return Token::kMemberSeparator;
            case '\0': return Token::kEnd;
            default: return Token::kError;
        }
    }

    // readValue at nesting |depth| (the root 1), its value counted by the caller: false at a
    // fault (the reader builds nothing more)
    bool value(size_t depth) {
        if (depth > kStackLimit)
            return false;
        Token t = token();
        while (t == Token::kComment) {
            t = token();
        }
        switch (t) {
            case Token::kObjectBegin: return object(depth);
            case Token::kArrayBegin: return array(depth);
            case Token::kString:
            case Token::kNumber:
            case Token::kLiteral: return true;
            default: return false;
        }
    }

    bool array(size_t depth) {
        for (;;) {
            skip_spaces();
            if (at_ < s_.size() && s_[at_] == ']') {
                ++at_;
                return true;
            }
            // the element exists before its value is read
            add();
            if (!value(depth + 1))
                return false;
            Token t = token();
            while (t == Token::kComment) {
                t = token();
            }
            if (t == Token::kArrayEnd)
                return true;
            if (t != Token::kArraySeparator)
                return false;
        }
    }

    bool object(size_t depth) {
        for (;;) {
            Token name = token();
            while (name == Token::kComment) {
                name = token();
            }
            if (name == Token::kObjectEnd)
                return true;
            if (name != Token::kString || token() != Token::kMemberSeparator)
                return false;
            add();
            if (!value(depth + 1))
                return false;
            Token separator = token();
            if (separator != Token::kObjectEnd && separator != Token::kArraySeparator
                    && separator != Token::kComment) {
                return false;
            }
            // after a comment, any next token is taken for the separator (jsoncpp's readObject)
            while (separator == Token::kComment) {
                separator = token();
            }
            if (separator == Token::kObjectEnd)
                return true;
        }
    }

    std::string_view s_;
    uint64_t max_;
    size_t at_ = 0;
    uint64_t values_ = 0;
};

inline uint64_t count_lenient_json_values(std::string_view content, uint64_t max_values) {
    return LenientJsonValues::count(content, max_values);
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
                refuse(path_ + ": unknown field '" + ::utils::clip_for_message(name) + "'");
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

/**
 * The volatile figures of a span of a thread's work (common::ThreadUsage), into a `timing`
 * block of /traverse, /resolve or /pattern: `cpu_ms`, the thread's CPU time (user and system),
 * and `page_faults` {minor, major} where the platform reads them per thread (Linux). A figure
 * the platform does not give is absent, never 0. Both vary between runs, as everything in
 * `timing` does.
 */
inline void put_thread_usage(Json::Value *timing, const common::ThreadUsage &d) {
    if (d.has_cpu)
        (*timing)["cpu_ms"] = d.cpu_ms();
    if (d.has_faults) {
        Json::Value faults;
        faults["minor"] = uint_json(d.minor_faults);
        faults["major"] = uint_json(d.major_faults);
        (*timing)["page_faults"] = std::move(faults);
    }
}

// The figures put_thread_usage wrote into |timing| (a part of a request measured on another
// thread, read back from its answer); a figure absent there is absent here
inline common::ThreadUsage thread_usage_of(const Json::Value &timing) {
    common::ThreadUsage u;
    if (!timing.isObject())
        return u;
    if (timing.isMember("cpu_ms") && timing["cpu_ms"].isNumeric()) {
        u.has_cpu = true;
        u.cpu_ns = static_cast<uint64_t>(std::llround(timing["cpu_ms"].asDouble() * 1e6));
    }
    const Json::Value &faults = timing.isMember("page_faults") ? timing["page_faults"]
                                                               : Json::Value::nullSingleton();
    if (faults.isObject() && faults["minor"].isIntegral() && faults["major"].isIntegral()) {
        u.has_faults = true;
        u.minor_faults = faults["minor"].asUInt64();
        u.major_faults = faults["major"].asUInt64();
    }
    return u;
}

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
