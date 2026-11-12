#ifndef __STRICT_JSON_HPP__
#define __STRICT_JSON_HPP__

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <json/json.h>


namespace mtg {
namespace cli {

// A request that is not valid under the schema of its route (an unknown field included)
class InvalidRequest : public std::invalid_argument {
  public:
    using std::invalid_argument::invalid_argument;
};

/**
 * A JSON object read field by field: every field read is known, and done() refuses any
 * other, so that a misspelled or unsupported field is an error rather than silently ignored.
 */
class Strict {
  public:
    Strict(const Json::Value &value, std::string path) : value_(value), path_(std::move(path)) {
        if (!value_.isObject())
            throw InvalidRequest(path_ + ": expected an object");
    }

    std::string child(const char *key) const { return path_ + "." + key; }

    bool has(const char *key) {
        known_.insert(key);
        return value_.isMember(key);
    }

    const Json::Value& raw(const char *key) {
        known_.insert(key);
        return value_[key];
    }

    std::string str(const char *key, const std::string &def) {
        if (!has(key))
            return def;
        if (!value_[key].isString())
            throw InvalidRequest(child(key) + ": expected a string");
        return value_[key].asString();
    }

    uint64_t uint(const char *key, uint64_t def, uint64_t min = 0,
                  uint64_t max = std::numeric_limits<uint64_t>::max()) {
        if (!has(key))
            return def;
        const Json::Value &v = value_[key];
        if (!v.isUInt64() || (v.isDouble() && v.asDouble() != std::floor(v.asDouble())))
            throw InvalidRequest(child(key) + ": expected a non-negative integer");
        const uint64_t x = v.asUInt64();
        if (x < min || x > max) {
            throw InvalidRequest(child(key) + ": out of range [" + std::to_string(min) + ", "
                                 + std::to_string(max) + "]");
        }
        return x;
    }

    double number(const char *key, double def, double min, double max) {
        if (!has(key))
            return def;
        const Json::Value &v = value_[key];
        if (!v.isNumeric())
            throw InvalidRequest(child(key) + ": expected a number");
        const double x = v.asDouble();
        if (!(x >= min && x <= max))
            throw InvalidRequest(child(key) + ": out of range");
        return x;
    }

    std::vector<std::string> strings(const char *key) {
        std::vector<std::string> out;
        if (!has(key))
            return out;
        const Json::Value &v = value_[key];
        if (!v.isArray())
            throw InvalidRequest(child(key) + ": expected an array of strings");
        for (const Json::Value &x : v) {
            if (!x.isString())
                throw InvalidRequest(child(key) + ": expected an array of strings");
            out.push_back(x.asString());
        }
        return out;
    }

    template <class T>
    T enumeration(const char *key, T def,
                  std::initializer_list<std::pair<const char*, T>> values) {
        if (!has(key))
            return def;
        std::string names;
        if (value_[key].isString()) {
            for (const auto &[name, x] : values) {
                if (value_[key].asString() == name)
                    return x;
            }
        }
        for (const auto &[name, x] : values) {
            names += std::string(names.empty() ? "" : ", ") + "\"" + name + "\"";
        }
        throw InvalidRequest(child(key) + ": expected one of " + names);
    }

    // refuse the fields that were not read
    void done() const {
        for (const std::string &name : value_.getMemberNames()) {
            if (!known_.count(name)) {
                throw InvalidRequest(child(name.c_str()) + ": unknown field");
            }
        }
    }

  private:
    const Json::Value &value_;
    std::string path_;
    std::set<std::string> known_;
};

} // namespace cli
} // namespace mtg

#endif // __STRICT_JSON_HPP__
