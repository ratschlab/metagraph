#include "server_request.hpp"

#include <stdexcept>


namespace mtg {
namespace cli {

SearchRequest parse_search_request(const Json::Value &json, const Config &server_config) {
    if (json["FASTA"].isNull())
        throw std::invalid_argument("No input sequences received from client");

    // `partial_refusal` changes what a 200 means (a unit may be answered as refused): only a
    // JSON boolean is read, anything else given is refused (lane T7-P1, review m4)
    const Json::Value &partial_refusal = json["partial_refusal"];
    if (json.isMember("partial_refusal") && !partial_refusal.isBool())
        throw std::invalid_argument("partial_refusal should be a boolean (true or false)");

    SearchRequest request{json["FASTA"].asString(), server_config,
                          json.get("align", false).asBool(),
                          partial_refusal.isBool() && partial_refusal.asBool()};
    auto &config = request.config;
    config.discovery_fraction
            = json.get("discovery_fraction", config.discovery_fraction).asDouble();
    config.alignment_min_exact_match
            = json.get("min_exact_match", config.alignment_min_exact_match).asDouble();
    config.alignment_max_nodes_per_seq_char = json.get(
            "max_num_nodes_per_seq_char", config.alignment_max_nodes_per_seq_char).asDouble();

    if (config.discovery_fraction < 0.0 || config.discovery_fraction > 1.0) {
        throw std::invalid_argument(
                "Discovery fraction should be within [0, 1.0]. Instead got "
                + std::to_string(config.discovery_fraction));
    }
    if (config.alignment_min_exact_match < 0.0
            || config.alignment_min_exact_match > 1.0) {
        throw std::invalid_argument(
                "Minimum exact match should be within [0, 1.0]. Instead got "
                + std::to_string(config.alignment_min_exact_match));
    }

    config.num_top_labels = json.get("top_labels", config.num_top_labels).asInt();
    if (json.get("query_coords", false).asBool()) {
        config.query_mode = COORDS;
    } else if (json.get("query_counts", false).asBool()) {
        config.query_mode = COUNTS;
    } else if (json.get("with_signature", false).asBool()) {
        config.query_mode = SIGNATURE;
    } else if (json.get("abundance_sum", false).asBool()) {
        config.query_mode = COUNTS_SUM;
    } else {
        config.query_mode = MATCHES;
    }
    return request;
}

} // namespace cli
} // namespace mtg
