#ifndef __CONFIGURED_PAIR_ID_HPP__
#define __CONFIGURED_PAIR_ID_HPP__

#include <boost/uuid/name_generator.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <optional>
#include <string>
#include <string_view>

namespace mtg::common {

// Configured-pair identity, not a file-content fingerprint. Keep the namespace
// and UTF-8 byte encoding in sync with the API's Python UUIDv5 implementation.
inline std::string configured_pair_id(std::optional<std::string_view> graph,
                                      std::optional<std::string_view> graph_path,
                                      std::optional<std::string_view> annotation_path) {
    static const auto ns = boost::uuids::string_generator()(
        "d3565c73-e7e9-5bca-b6a9-64cb5c4c9f93");
    std::string name;
    for (const auto &field : { graph, graph_path, annotation_path }) {
        if (field) {
            name += std::to_string(field->size()) + ":";
            name.append(field->data(), field->size());
        } else {
            name += "-1:";
        }
    }
    boost::uuids::name_generator generate(ns);
    return "n1-" + boost::uuids::to_string(generate(name.data(), name.size()));
}

} // namespace mtg::common

#endif // __CONFIGURED_PAIR_ID_HPP__
