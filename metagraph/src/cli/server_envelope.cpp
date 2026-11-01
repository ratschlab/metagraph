#include "server_envelope.hpp"

#include <algorithm>
#include <initializer_list>
#include <unordered_set>
#include <utility>

#include "graph/alignment/pattern_search.hpp"
#include "traverse.hpp"
#include "traverse_attempts.hpp"


namespace mtg {
namespace cli {

namespace {

// the tags of a pair's entry, set on |entry|
void tag_pair(Json::Value *entry, const SelectedPair &pair, const std::string &index_fp) {
    (*entry)["graph"] = pair.name;
    (*entry)["graph_path"] = pair.pair.first;
    (*entry)["annotation_path"] = pair.pair.second;
    (*entry)["index_fp"] = index_fp.empty() ? Json::Value() : Json::Value(index_fp);
}

// the pair an entry is tagged with, as not_reached and unfinished list it
Json::Value ref_of_entry(const Json::Value &entry) {
    Json::Value ref(Json::objectValue);
    for (const char *key : { "graph", "graph_path", "annotation_path" }) {
        ref[key] = entry[key];
    }
    return ref;
}

// The entries' outcomes, counted, and the pairs of the refusals section 4 lists apart
struct Outcomes {
    Json::UInt answered = 0;
    Json::UInt refused = 0;
    Json::Value not_reached = Json::Value(Json::arrayValue);
    Json::Value unfinished = Json::Value(Json::arrayValue);
    uint64_t memory_bytes = 0;
};

Outcomes outcomes_of(const std::vector<Json::Value> &entries) {
    Outcomes out;
    for (const Json::Value &entry : entries) {
        if (entry["outcome"].asString() != "refused") {
            out.answered++;
            // the answered pairs' accounts together: the pairs run at once on the pool, so at
            // most what they held together
            const Json::Value &bytes = entry["work"]["memory_bytes"];
            if (bytes.isUInt64())
                out.memory_bytes += bytes.asUInt64();
            continue;
        }
        out.refused++;
        const Json::Value &code = entry["refusal"]["code"];
        if (code == "not_reached") {
            out.not_reached.append(ref_of_entry(entry));
        } else if (code == "unfinished") {
            out.unfinished.append(ref_of_entry(entry));
        }
    }
    return out;
}

Json::Value iso_or_null(const std::optional<std::chrono::system_clock::time_point> &t) {
    return t ? Json::Value(iso_utc(*t)) : Json::Value();
}

Json::Value ms_or_null(const std::optional<uint64_t> &ms) {
    return ms ? Json::Value(Json::UInt64(*ms)) : Json::Value();
}

} // namespace

Json::Value envelope_pair_ref(const SelectedPair &pair) {
    return pair_ref_json(pair);
}

Json::Value envelope_entry_answered(const SelectedPair &pair, const std::string &index_fp,
                                    Json::Value answer) {
    tag_pair(&answer, pair, index_fp);
    answer["outcome"] = "answered";
    return answer;
}

Json::Value envelope_entry_refused(const SelectedPair &pair, const std::string &index_fp,
                                   int http_status, const Json::Value &refusal_body) {
    Json::Value entry(Json::objectValue);
    tag_pair(&entry, pair, index_fp);
    entry["outcome"] = "refused";
    Json::Value refusal = refusal_body;
    refusal["http_status"] = http_status;
    entry["refusal"] = std::move(refusal);
    return entry;
}

Json::Value deadline_json(std::optional<uint64_t> not_after_ms,
                          std::optional<uint64_t> finish_by_ms,
                          std::optional<std::chrono::system_clock::time_point> effective,
                          std::optional<std::chrono::system_clock::time_point> cutoff) {
    Json::Value deadline(Json::objectValue);
    deadline["not_after_ms"] = ms_or_null(not_after_ms);
    deadline["finish_by_ms"] = ms_or_null(finish_by_ms);
    deadline["effective_deadline_at"] = iso_or_null(effective);
    deadline["cutoff_at"] = iso_or_null(cutoff);
    return deadline;
}

Json::Value post_envelope(const std::string &route,
                          const std::vector<std::string> &names,
                          std::vector<Json::Value> entries,
                          const EnvelopeContext &context,
                          double elapsed_ms) {
    Outcomes outcomes = outcomes_of(entries);
    Json::Value out(Json::objectValue);
    out["envelope_version"] = kEnvelopeVersion;
    out["route"] = route;
    if (route == "/pattern")
        out["pattern_contract_version"] = graph::pattern::kPatternContractVersion;
    Json::Value graphs(Json::arrayValue);
    for (const std::string &name : names) {
        graphs.append(name);
    }
    out["graphs"] = std::move(graphs);
    out["answered"] = outcomes.answered;
    out["refused"] = outcomes.refused;
    out["not_reached"] = std::move(outcomes.not_reached);
    out["unfinished"] = std::move(outcomes.unfinished);
    out["request_id"] = context.request_id ? Json::Value(*context.request_id) : Json::Value();
    out["server_instance"] = context.server_instance;
    out["deadline"] = context.deadline.isNull()
            ? deadline_json(std::nullopt, std::nullopt, std::nullopt, std::nullopt)
            : context.deadline;
    if (!context.usage.isNull())
        out["usage"] = context.usage;
    Json::Value answers(Json::arrayValue);
    for (Json::Value &entry : entries) {
        answers.append(std::move(entry));
    }
    out["answers"] = std::move(answers);
    out["timing"]["elapsed_ms"] = elapsed_ms;
    out["work"]["memory_bytes"] = Json::Value(Json::UInt64(outcomes.memory_bytes));
    return out;
}

Json::Value get_envelope(const std::string &route,
                         std::vector<Json::Value> entries,
                         const std::string &server_instance) {
    const Outcomes outcomes = outcomes_of(entries);
    Json::Value out(Json::objectValue);
    out["envelope_version"] = kEnvelopeVersion;
    out["route"] = route;
    // the entries' names: in byte order, as the entries come (ordered_pairs)
    Json::Value graphs(Json::arrayValue);
    std::unordered_set<std::string> seen;
    for (const Json::Value &entry : entries) {
        std::string name = entry["graph"].asString();
        if (seen.insert(name).second)
            graphs.append(std::move(name));
    }
    out["graphs"] = std::move(graphs);
    out["answered"] = outcomes.answered;
    out["refused"] = outcomes.refused;
    out["server_instance"] = server_instance;
    Json::Value answers(Json::arrayValue);
    for (Json::Value &entry : entries) {
        answers.append(std::move(entry));
    }
    out["answers"] = std::move(answers);
    return out;
}

Json::Value composite_fp_json(std::vector<std::string> index_fps) {
    // one pair without a manifest leaves the set of indexes unidentified: no composite
    if (index_fps.empty()
            || std::any_of(index_fps.begin(), index_fps.end(),
                           [](const std::string &fp) { return fp.empty(); })) {
        return Json::Value();
    }
    // the order of the list's lines and names does not enter it: the same indexes served in
    // another order, or under other names, state the same composite (the names are the list
    // beside it, `graphs` and graph_summary)
    std::sort(index_fps.begin(), index_fps.end());
    std::string joined;
    for (size_t i = 0; i < index_fps.size(); ++i) {
        if (i)
            joined += '\n';
        joined += index_fps[i];
    }
    return Json::Value(sha256_hex(joined));
}

Json::Value multigraph_capabilities_json() {
    auto strings = [](std::initializer_list<const char *> values) {
        Json::Value out(Json::arrayValue);
        for (const char *value : values) {
            out.append(value);
        }
        return out;
    };
    Json::Value m(Json::objectValue);
    m["envelope_version"] = kEnvelopeVersion;
    m["envelope_routes"] = strings({ "/resolve", "/traverse", "/pattern",
                                     "/traverse/capabilities", "/pattern/capabilities" });
    // the selectors of section 6, as this release serves them: /search has no re-ask, so no
    // pairs; a /traverse request runs on one pair however it is named; a capabilities probe
    // names its one pair by its parameters
    Json::Value selectors(Json::objectValue);
    selectors["/search"] = strings({ "graphs" });
    selectors["/resolve"] = strings({ "graphs", "pairs" });
    selectors["/pattern"] = strings({ "graphs", "pairs" });
    selectors["/traverse"] = strings({ "graph", "graphs", "pairs" });
    selectors["/traverse/capabilities"] = strings({ "graph", "graph_path", "annotation_path" });
    selectors["/pattern/capabilities"] = strings({ "graph", "graph_path", "annotation_path" });
    m["selectors"] = std::move(selectors);
    m["traverse_pairs_per_request"] = 1;
    // the deadlines (section 4) and the request handle (section 5) by route: /traverse's
    // handle is its attempt_id
    m["not_after_ms"] = strings({ "/traverse", "/resolve", "/pattern" });
    m["finish_by_ms"] = strings({ "/traverse", "/resolve", "/pattern" });
    m["request_id"] = strings({ "/resolve", "/pattern" });
    m["requests_route"] = kRequestsRoute;
    m["unknown_labels"] = true;
    m["search_cuts"] = true;
    m["resolve_boundary_ties"] = true;
    // a discovery under support "trace" ranks by traced k-mers (section 7), as every
    // /traverse capabilities block states it
    m["discovery_trace_scale"] = "traced";
    return m;
}

} // namespace cli
} // namespace mtg
