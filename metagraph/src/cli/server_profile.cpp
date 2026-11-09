#include "server_profile.hpp"

#include <cstdio>
#include <stdexcept>
#include <utility>

#include "graph/alignment/pattern_search.hpp"
#include "graph/representation/bitmap/dbg_bitmap.hpp"
#include "graph/representation/hash/dbg_hash_fast.hpp"
#include "graph/representation/hash/dbg_hash_ordered.hpp"
#include "graph/representation/hash/dbg_hash_string.hpp"
#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "annotation/coord_to_header.hpp"
#include "graph/traversal/label_oracle.hpp"
#include "graph/traversal/resolve.hpp"
#include "cli/config/config.hpp"
#include "cli/json_helpers.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "cli/load/load_annotation.hpp"
#include "pattern.hpp"
#include "traverse.hpp"


namespace mtg {
namespace cli {

std::string annotation_form(const graph::AnnotatedDBG::Annotator &annotation) {
    // parse_annotation_type decides by the file's extension, as the loader did; every type
    // round-trips (GraphletServer.InventoryTableMatchesTheLoadersTypes)
    return Config::annotype_to_string(parse_annotation_type("x" + annotation.file_extension()));
}

std::string graph_representation(const graph::DeBruijnGraph &graph) {
    using namespace graph;
    // the loaders' extensions (load_critical_dbg), each with its `--graph` name
    static const std::pair<const char*, const char*> kRepresentations[] = {
        { DBGSuccinct::kExtension, "succinct" },
        { DBGHashOrdered::kExtension, "hash" },
        { DBGHashString::kExtension, "hashstr" },
        { DBGHashFast::kExtension, "hashfast" },
        { DBGBitmap::kExtension, "bitmap" },
        { DBGSSHash::kExtension, "sshash" },
    };
    // a CanonicalDBG states the extension of the graph it wraps (DBGWrapper)
    const std::string extension = graph.file_extension();
    for (const auto &[ext, name] : kRepresentations) {
        if (extension == ext)
            return name;
    }
    throw std::logic_error("graph_representation: no representation has the extension "
                           + extension);
}

const std::vector<std::string>& support_fields() {
    // graph representation, regime, alphabet, k, annotation form, coordinates, .seqs and what
    // its load verified of its column order (seqs_order: what a header answer certifies), trace,
    // decode charged, counts and the pattern fields (graph_mode, available,
    // unavailable_reason, mask, counting); not the identity (graph, graph_path,
    // annotation_path, index_ns, index_fp, index_meta_fp) nor the size (traversal.num_labels)
    static const std::vector<std::string> kFields = {
        "alphabet",
        "annotation_form",
        "available",
        "counting",
        "decode_charged",
        "graph_mode",
        "graph_representation",
        "k",
        "mask",
        "seqs_order",
        "traversal.counts",
        "traversal.has_coord_to_header",
        "traversal.has_coordinates",
        "traversal.regime",
        "traversal.supports_trace",
        "unavailable_reason",
    };
    return kFields;
}

std::string support_text(const Json::Value &pair) {
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    std::string text = std::string(kSupportVersion) + "\n";
    for (const std::string &name : support_fields()) {
        const size_t dot = name.find('.');
        const Json::Value &value = dot == std::string::npos
                ? pair.get(name, Json::Value())
                : pair.get(name.substr(0, dot), Json::Value())
                      .get(name.substr(dot + 1), Json::Value());
        if (value.isArray() || value.isObject())
            throw std::logic_error("support_fp: the field " + name + " is not a scalar");
        text += name + "=" + Json::writeString(builder, value) + "\n";
    }
    return text;
}

std::string support_fingerprint(const Json::Value &pair) {
    const uint64_t h = graph::traversal::fnv1a64(support_text(pair));
    char buf[17];
    std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(h));
    return buf;
}

Json::Value pair_summary_json(const graph::AnnotatedDBG &index, const IndexIdentity &identity,
                              const Json::Value &graph, const Json::Value &graph_path,
                              const Json::Value &annotation_path) {
    const graph::DeBruijnGraph &dbg = index.get_graph();
    Json::Value p;
    p["graph"] = graph;
    p["graph_path"] = graph_path;
    p["annotation_path"] = annotation_path;
    p["index_ns"] = string_or_null(identity.name);
    p["index_fp"] = string_or_null(identity.fp);
    p["index_meta_fp"] = identity.meta_fp;
    p["k"] = uint_json(dbg.get_k());
    p["graph_representation"] = graph_representation(dbg);
    p["annotation_form"] = annotation_form(index.get_annotator());
    // the pattern search, as the pair's block states it (pattern_capabilities_json)
    const graph::pattern::GraphSupport support = route_support(dbg);
    const bool recognised = support.supported || support.reason == "alphabet_unsupported"
                                || support.reason == "alphabet_untested"
                                || support.reason == "mask_invalid";
    p["available"] = support.supported;
    p["unavailable_reason"] = support.supported ? Json::Value() : Json::Value(support.reason);
    p["graph_mode"] = recognised ? Json::Value(to_string(support.mode)) : Json::Value();
    p["mask"] = !recognised ? Json::Value()
              : !support.mask_present ? Json::Value("absent")
              : mask_built_at_load(dbg) ? Json::Value("built_at_load")
                                        : Json::Value("file");
    p["counting"] = !support.supported ? Json::Value()
                  : support.mask_present ? Json::Value("exact")
                                         : Json::Value("upper_bound");
    // what a traversal reads, as the capabilities of every /traverse response state it
    const graph::traversal::LabelOracle oracle(index);
    p["alphabet"] = oracle.graph().alphabet();
    // whether the annotation reads are charged under a memory budget (the oracle's predicate)
    p["decode_charged"] = oracle.decode_charged();
    // what the load verified of the record mapping's column order (check_record_mapping): its
    // stated column labels matched ("verified"), or it states none and its column count alone
    // was checked ("count_only"); null without a mapping
    const annot::CoordToHeader *mapping = oracle.coord_to_header();
    p["seqs_order"] = !mapping ? Json::Value()
                    : mapping->column_order_verified() ? Json::Value("verified")
                                                       : Json::Value("count_only");
    Json::Value t;
    t["regime"] = to_string(oracle.regime());
    t["num_labels"] = uint_json(oracle.num_columns());
    t["has_coordinates"] = oracle.has_coordinates();
    t["has_coord_to_header"] = oracle.coord_to_header() != nullptr;
    t["supports_trace"] = oracle.has_coordinates()
            && oracle.regime() == graph::traversal::Regime::BASIC;
    // a count annotation is read as presence (the probe states counts only there; null here)
    t["counts"] = oracle.has_counts() ? Json::Value("presence_only") : Json::Value();
    p["traversal"] = std::move(t);
    p["support_fp"] = support_fingerprint(p);
    return p;
}

} // namespace cli
} // namespace mtg
