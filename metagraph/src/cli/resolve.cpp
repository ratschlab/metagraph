#include "resolve.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <limits>

#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/traversal/resolve.hpp"

namespace mtg::cli {
namespace {
using namespace graph::traversal;

void members(const Json::Value &value, std::initializer_list<std::string_view> names) {
    if (!value.isObject()) throw std::invalid_argument("Expected a JSON object");
    for (const auto &name : value.getMemberNames())
        if (std::find(names.begin(), names.end(), name) == names.end())
            throw std::invalid_argument("Unsupported resolve field: " + name);
}
uint64_t number(const Json::Value &v, const char *name, uint64_t fallback,
                uint64_t minimum = 0) {
    if (!v.isMember(name)) return fallback;
    if (!v[name].isUInt64() || v[name].asUInt64() < minimum)
        throw std::invalid_argument(std::string("Invalid nonnegative integer: ") + name);
    return v[name].asUInt64();
}
std::string string(const Json::Value &v, const char *name, const char *fallback) {
    if (!v.isMember(name)) return fallback;
    if (!v[name].isString()) throw std::invalid_argument(std::string("Expected string: ") + name);
    return v[name].asString();
}
std::vector<std::string> strings(const Json::Value &v, const char *name) {
    std::vector<std::string> out;
    if (!v.isMember(name)) return out;
    if (!v[name].isArray()) throw std::invalid_argument(std::string("Expected array: ") + name);
    std::set<std::string> seen;
    for (const auto &item : v[name]) {
        if (!item.isString() || !seen.insert(item.asString()).second)
            throw std::invalid_argument(std::string("Expected distinct strings: ") + name);
        out.push_back(item.asString());
    }
    return out;
}
Json::Value uint_json(uint64_t n) { return Json::UInt64(n); }
Json::Value interval_json(const KmerInterval &iv) {
    Json::Value out(Json::arrayValue); out.append(uint_json(iv.begin)); out.append(uint_json(iv.end));
    return out;
}

// Profile/selection serialization extracted from the prototype's resolve route.
Json::Value profile_to_json(const SupportProfile &p, const SeedSelection *sel, const std::string &run_format,
                            const char *regime) {
    Json::Value j;
    j["k"] = uint_json(p.k);
    j["regime"] = regime;
    j["num_kmers"] = uint_json(p.num_kmers);
    j["support"] = "kmer";
    auto runs_json = [&](const std::vector<KmerInterval> &runs) -> Json::Value {
        if (run_format == "rle")
            return Json::Value(encode_runs(runs, p.num_kmers));
        Json::Value arr(Json::arrayValue);
        for (const auto &r : runs) arr.append(interval_json(r));
        return arr;
    };
    j["graph_runs"] = runs_json(p.graph_runs);
    Json::Value labels(Json::arrayValue);
    for (const auto &l : p.labels) {
        Json::Value lj;
        lj["label"] = l.label.name;
        lj["kind"] = "column";
        lj["kmers_supported"] = uint_json(l.kmers_supported);
        lj["runs"] = runs_json(l.runs);

        labels.append(std::move(lj));
    }
    j["labels"] = std::move(labels);
    if (p.labels_truncated) {
        Json::Value t;
        t["kept"] = uint_json(p.labels_truncated->kept);
        t["total"] = uint_json(p.labels_truncated->total);
        t["min_kept_kmers"] = uint_json(p.labels_truncated->min_kept_kmers);
        t["max_dropped_kmers"] = uint_json(p.labels_truncated->max_dropped_kmers);
        t["dropped_full_length"] = uint_json(p.labels_truncated->dropped_full_length);
        j["labels_truncated"] = std::move(t);
    } else {
        j["labels_truncated"] = Json::Value::null;
    }
    Json::Value cands(Json::arrayValue);
    for (const auto &c : p.candidates) {
        Json::Value cj;
        cj["kmer_interval"] = interval_json(c.kmers);
        auto [a, b] = p.bp_interval(c.kmers);
        cj["bp_interval"] = interval_json({ a, b });
        Json::Value ls(Json::arrayValue);
        for (LabelId l : c.labels) ls.append(p.labels[l].label.name);
        cj["labels"] = std::move(ls);
        cands.append(std::move(cj));
    }
    j["candidates"] = std::move(cands);
    if (sel) {
        Json::Value s;
        Json::Value pol;
        pol["policy"] = sel->policy.policy == SelectionPolicy::LONGEST_FIRST ? "longest_first"
                      : sel->policy.policy == SelectionPolicy::MAX_SUPPORT ? "max_support" : "explicit";
        pol["max_seeds"] = uint_json(sel->policy.max_seeds);
        pol["min_block_bp"] = uint_json(sel->policy.min_block_bp);
        pol["max_labels_per_seed"] = uint_json(sel->policy.max_labels_per_seed);
        pol["label_order"] = sel->policy.label_order == SelectionPolicy::HASH ? "hash"
                           : sel->policy.label_order == SelectionPolicy::COLUMN_ID ? "column_id" : "kmers_supported";
        pol["sample_seed"] = uint_json(sel->policy.sample_seed);
        pol["merge_overlapping"] = sel->policy.merge_overlapping;
        s["policy"] = std::move(pol);
        Json::Value counts;
        counts["candidates"] = uint_json(sel->num_candidates);
        counts["eligible"] = uint_json(sel->num_eligible);
        counts["selected"] = uint_json(sel->seeds.size());
        s["counts"] = std::move(counts);
        Json::Value seeds(Json::arrayValue);
        for (const auto &seed : sel->seeds) {
            Json::Value sj;
            sj["seed_id"] = seed.seed_id;
            sj["sequence"] = seed.sequence;
            sj["kmer_interval"] = interval_json(seed.kmers);
            Json::Value ls(Json::arrayValue);
            for (const auto &l : seed.labels) ls.append(l);
            sj["labels"] = std::move(ls);
            Json::Value pop;
            pop["supporting_total"] = uint_json(seed.population.supporting_total);
            pop["included"] = uint_json(seed.population.included);
            pop["dropped_count"] = uint_json(seed.population.dropped_count);
            if (seed.population.dropped_count) {
                pop["dropped_digest"] = seed.population.dropped_digest;
                Json::Value dl(Json::arrayValue);
                for (const auto &d : seed.population.dropped) dl.append(d);
                pop["dropped"] = std::move(dl);
            }
            sj["label_population"] = std::move(pop);
            Json::Value ov(Json::arrayValue);
            for (size_t o : seed.overlaps_with) ov.append(uint_json(o));
            sj["overlaps_with"] = std::move(ov);
            if (!seed.labels_not_covering.empty()) {
                Json::Value nc(Json::arrayValue);
                for (const auto &l : seed.labels_not_covering) nc.append(l);
                sj["labels_not_covering"] = std::move(nc);
            }
            seeds.append(std::move(sj));
        }
        s["seeds"] = std::move(seeds);
        j["selection"] = std::move(s);
    }
    return j;
}
} // namespace

Json::Value process_resolve_request(const Json::Value &request, const graph::AnnotatedDBG &graph) {
    members(request, {"sequence", "labels", "discover", "support", "min_block_kmers",
                      "run_format", "select", "bounds"});
    auto sequence = string(request, "sequence", "");
    if (sequence.empty()) throw std::invalid_argument("sequence is required");
    uint64_t max_query_bp = 1'000'000;
    if (request.isMember("bounds")) {
        members(request["bounds"], {"max_query_bp"});
        max_query_bp = number(request["bounds"], "max_query_bp", max_query_bp, 1);
    }
    if (sequence.size() > max_query_bp)
        throw std::invalid_argument("sequence exceeds bounds.max_query_bp");
    if (string(request, "support", "kmer") != "kmer")
        throw std::invalid_argument("This resolve increment supports kmer membership only");
    auto format = string(request, "run_format", "intervals");
    if (format != "intervals" && format != "rle") throw std::invalid_argument("Unknown run_format");
    ResolveOptions options;
    options.labels = strings(request, "labels");
    options.min_block_kmers = number(request, "min_block_kmers", 1, 1);
    if (request.isMember("discover")) {
        members(request["discover"], {"max_labels", "kind"});
        options.discover = true;
        options.discover_max_labels = number(request["discover"], "max_labels", 1000, 1);
        if (string(request["discover"], "kind", "column") != "column")
            throw std::invalid_argument("This resolve increment supports column labels only");
    }
    SelectionPolicy policy;
    if (request.isMember("select")) {
        const auto &p = request["select"];
        members(p, {"policy", "max_seeds", "min_block_bp", "max_labels_per_seed",
                    "label_order", "sample_seed", "merge_overlapping", "seeds"});
        auto choice = string(p, "policy", "max_support");
        if (choice == "longest_first") policy.policy = SelectionPolicy::LONGEST_FIRST;
        else if (choice == "max_support") policy.policy = SelectionPolicy::MAX_SUPPORT;
        else if (choice == "explicit") policy.policy = SelectionPolicy::EXPLICIT;
        else throw std::invalid_argument("Unknown selection policy");
        policy.max_seeds = number(p, "max_seeds", 10, 1);
        policy.min_block_bp = number(p, "min_block_bp", 0);
        policy.max_labels_per_seed = number(p, "max_labels_per_seed", 1000, 1);
        policy.sample_seed = number(p, "sample_seed", 0);
        auto order = string(p, "label_order", "hash");
        if (order == "hash") policy.label_order = SelectionPolicy::HASH;
        else if (order == "column_id") policy.label_order = SelectionPolicy::COLUMN_ID;
        else if (order == "kmers_supported") policy.label_order = SelectionPolicy::KMERS_SUPPORTED;
        else throw std::invalid_argument("Unknown label_order");
        if (p.isMember("merge_overlapping")) {
            if (!p["merge_overlapping"].isBool()) throw std::invalid_argument("Expected boolean: merge_overlapping");
            policy.merge_overlapping = p["merge_overlapping"].asBool();
        }
        if (p.isMember("seeds")) {
            if (policy.policy != SelectionPolicy::EXPLICIT)
                throw std::invalid_argument("seeds requires explicit selection");
            if (!p["seeds"].isArray()) throw std::invalid_argument("Expected array: seeds");
            for (const auto &value : p["seeds"]) {
                members(value, {"kmer_interval", "labels"});
                const auto &iv = value["kmer_interval"];
                if (!iv.isArray() || iv.size() != 2 || !iv[0].isUInt64() || !iv[1].isUInt64())
                    throw std::invalid_argument("Expected nonnegative [begin,end] kmer_interval");
                policy.explicit_seeds.push_back({{iv[0].asUInt64(), iv[1].asUInt64()}, strings(value, "labels")});
            }
        }
        if (policy.policy == SelectionPolicy::EXPLICIT && policy.explicit_seeds.empty())
            throw std::invalid_argument("Explicit selection requires seeds");
    }
    auto profile = resolve_support(graph, sequence, options);
    auto mode = graph.get_graph().get_mode();
    bool primary = dynamic_cast<const graph::CanonicalDBG *>(&graph.get_graph());
    const char *regime = primary || mode == graph::DeBruijnGraph::PRIMARY ? "primary"
                       : mode == graph::DeBruijnGraph::CANONICAL ? "canonical" : "basic";
    std::optional<SeedSelection> selection;
    if (request.isMember("select"))
        selection = select_seeds(profile, sequence, policy, mode != graph::DeBruijnGraph::BASIC);
    return profile_to_json(profile, selection ? &*selection : nullptr, format, regime);
}

int resolve_graph(Config *config) {
    auto graph = initialize_annotated_dbg(*config);
    Json::StreamWriterBuilder writer; writer["indentation"] = "";
    int status = 0;
    for (const auto &path : config->fnames) {
        Json::Value response;
        try {
            std::ifstream input(path);
            if (!input) throw std::invalid_argument("Cannot open request file: " + path);
            Json::CharReaderBuilder reader;
            Json::CharReaderBuilder::strictMode(&reader.settings_);
            Json::Value request; std::string errors;
            if (!Json::parseFromStream(reader, input, &request, &errors))
                throw std::invalid_argument(errors);
            response = process_resolve_request(request, *graph);
        } catch (const std::exception &error) {
            response["error"] = error.what(); status = 1;
        }
        std::cout << Json::writeString(writer, response) << '\n';
    }
    return status;
}
} // namespace mtg::cli
