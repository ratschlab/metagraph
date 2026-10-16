// Experimental JSON-lines driver for the coordinate-backed context cursor.
// Uses the existing query loader (including mmap and row-diff annotation setup).
#include "metagraph.hpp"
#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "graph/traversal/context_cursor.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <json/json.h>

using mtg::graph::traversal::ContextCursor;

namespace {
void check_members(const Json::Value &object, std::initializer_list<std::string_view> allowed) {
    if (!object.isObject())
        throw std::invalid_argument("Expected an object");
    for (const auto &name : object.getMemberNames())
        if (std::ranges::find(allowed, name) == allowed.end())
            throw std::invalid_argument("Unknown field: " + name);
}

size_t bounded_size(const Json::Value &request, const char *name, size_t fallback, size_t limit) {
    if (!request.isMember(name))
        return fallback;
    if (!request[name].isUInt64() || request[name].asUInt64() > limit)
        throw std::invalid_argument(std::string(name) + " is outside the example's limit");
    return request[name].asUInt64();
}

const char *stop_name(ContextCursor::StopReason reason) {
    using enum ContextCursor::StopReason;
    switch (reason) {
        case NO_OCCURRENCES: return "no_occurrences";
        case RADIUS_LIMIT: return "radius_limit";
        case RECORD_END: return "record_end";
        case UNINDEXED: return "unindexed";
        case BRANCH: return "branch";
    }
    throw std::logic_error("Invalid closure stop reason");
}

Json::Value closure_json(const ContextCursor::Closure &closure) {
    Json::Value out;
    out["left_added"] = Json::UInt64(closure.left);
    out["right_added"] = Json::UInt64(closure.right);
    out["left_stop"] = stop_name(closure.left_stop);
    out["right_stop"] = stop_name(closure.right_stop);
    return out;
}

Json::Value snapshot(ContextCursor &cursor, const mtg::graph::AnnotatedDBG &graph, size_t seed_offset) {
    Json::Value out;
    out["sequence"] = cursor.sequence();
    out["seed_offset"] = Json::UInt64(seed_offset);
    out["occurrences"] = Json::Value(Json::arrayValue);
    out["choices"] = Json::Value(Json::arrayValue);
    for (const auto &hit : cursor.occurrences()) {
        Json::Value value;
        value["column"] = graph.get_annotator().get_label_encoder().decode(hit.column);
        value["record"] = Json::UInt64(hit.record);
        value["header"] = graph.get_coord_to_header()->get_headers(hit.column)[hit.record];
        value["offset"] = Json::UInt64(hit.offset);
        out["occurrences"].append(value);
    }
    for (const auto &choice : cursor.choices()) {
        Json::Value value;
        value["left"] = choice.left ? std::string(1, choice.left) : "";
        value["right"] = choice.right ? std::string(1, choice.right) : "";
        value["count"] = Json::UInt64(choice.count);
        out["choices"].append(value);
    }
    return out;
}

Json::Value query(const Json::Value &request, const mtg::graph::AnnotatedDBG &graph) {
    check_members(request, {"seed", "columns", "max_occurrences", "radius", "steps"});
    if (!request["seed"].isString() || request["seed"].asString().size() > 10000)
        throw std::invalid_argument("Expected a seed string of at most 10000 bases");
    const auto radius = bounded_size(request, "radius", 1000, 10000);
    const auto max_occurrences = bounded_size(request, "max_occurrences", 100000, 100000);
    if (!max_occurrences)
        throw std::invalid_argument("max_occurrences must be positive");
    std::vector<std::string> columns;
    if (request.isMember("columns")) {
        if (!request["columns"].isArray())
            throw std::invalid_argument("columns must be an array");
        for (const auto &column : request["columns"]) {
            if (!column.isString())
                throw std::invalid_argument("columns must contain strings");
            columns.push_back(column.asString());
        }
    }
    struct Step { ContextCursor::Side side; char base; size_t radius; };
    std::vector<Step> steps;
    if (request.isMember("steps")) {
        if (!request["steps"].isArray() || request["steps"].size() > 100)
            throw std::invalid_argument("steps must be an array of at most 100 choices");
        for (const auto &step : request["steps"]) {
            check_members(step, {"side", "base", "radius"});
            if (!step["base"].isString() || step["base"].asString().size() != 1
                    || (step["side"] != "left" && step["side"] != "right"))
                throw std::invalid_argument("Each step needs side (left/right) and one base");
            char base = step["base"].asString()[0];
            if (std::string_view("ACGTacgt").find(base) == std::string_view::npos)
                throw std::invalid_argument("Step base must be A, C, G or T");
            steps.push_back({step["side"] == "left" ? ContextCursor::Side::LEFT : ContextCursor::Side::RIGHT,
                             base, bounded_size(step, "radius", radius, 10000)});
        }
    }

    // Validate the complete request before decoding its seed or executing steps.
    ContextCursor cursor(graph, request["seed"].asString(), columns, max_occurrences);
    auto closure = cursor.close(radius);
    size_t seed_offset = closure.left;
    auto result = snapshot(cursor, graph, seed_offset);
    result["closure"] = closure_json(closure);
    result["steps"] = Json::Value(Json::arrayValue);
    for (const auto &step : steps) {
        bool extended = cursor.extend(step.side, step.base);
        if (extended) {
            seed_offset += step.side == ContextCursor::Side::LEFT;
            closure = cursor.close(step.radius);
            seed_offset += closure.left;
        }
        auto value = snapshot(cursor, graph, seed_offset);
        value["extended"] = extended;
        value["closure"] = extended ? closure_json(closure) : Json::Value();
        result["steps"].append(value);
    }
    return result;
}
} // namespace

int main(int argc, char **argv) {
    if (argc < 2 || std::string(argv[1]) != "query") {
        std::cerr << "Usage: context_cursor query [--mmap] -i graph.dbg -a annotation.annodbg requests.jsonl\n";
        return 1;
    }
    mtg::cli::Config config(argc, argv);
    auto graph = mtg::cli::initialize_annotated_dbg(config);
    Json::StreamWriterBuilder writer;
    writer["indentation"] = "";
    bool failed = false;
    for (const auto &filename : config.fnames) {
        std::ifstream input(filename);
        if (!input) {
            std::cerr << "Cannot open " << filename << '\n';
            return 1;
        }
        std::string line;
        while (std::getline(input, line)) {
            if (line.empty())
                continue;
            Json::Value request, result;
            try {
                Json::CharReaderBuilder builder;
                Json::CharReaderBuilder::strictMode(&builder.settings_);
                std::unique_ptr<Json::CharReader> reader(builder.newCharReader());
                std::string error;
                if (!reader->parse(line.data(), line.data() + line.size(), &request, &error))
                    throw std::invalid_argument(error);
                result = query(request, *graph);
            } catch (const std::exception &error) {
                result["error"] = error.what();
                failed = true;
            }
            std::cout << Json::writeString(writer, result) << '\n';
        }
    }
    return failed ? 1 : 0;
}
