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
    // a unit (a FASTA record on a pair) past the answer budget is refused alone and reported in
    // its cut entry, not the whole request (SPEC-multigraph-envelope.md section 3, lane T7-P1)
    bool partial_refusal = false;
};

SearchRequest parse_search_request(const Json::Value &json, const Config &server_config);

} // namespace cli
} // namespace mtg

#endif // __SERVER_REQUEST_HPP__
