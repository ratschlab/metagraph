#ifndef __SERVER_REQUEST_HPP__
#define __SERVER_REQUEST_HPP__

#include <json/json.h>

#include "config/config.hpp"


namespace mtg {
namespace cli {

// Resolve request overrides once for both query execution and cache identity.
struct SearchRequest {
    std::string fasta;
    Config config;
    bool align;
};

SearchRequest parse_search_request(const Json::Value &json, const Config &server_config);

} // namespace cli
} // namespace mtg

#endif // __SERVER_REQUEST_HPP__
