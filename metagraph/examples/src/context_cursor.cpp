// Experimental JSON-lines driver for the coordinate-backed context cursor.
// Uses the existing query loader (including mmap and row-diff annotation setup).
#include "metagraph.hpp"
#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "graph/traversal/context_cursor.hpp"

#include <fstream>
#include <iostream>
#include <json/json.h>

using mtg::graph::traversal::ContextCursor;

namespace {
size_t bounded_size(const Json::Value &request, const char *name, size_t fallback, size_t limit) {
    if (!request.isMember(name))
        return fallback;
    if (!request[name].isUInt64() || request[name].asUInt64() > limit)
        throw std::invalid_argument(std::string(name) + " is outside the example's limit");
    return request[name].asUInt64();
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
    if (!request.isObject() || !request["seed"].isString() || request["seed"].asString().size() > 10000)
        throw std::invalid_argument("Expected a seed string of at most 10000 bases");
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
    ContextCursor cursor(graph, request["seed"].asString(), columns,
                         bounded_size(request, "max_occurrences", 100000, 100000));
    size_t seed_offset = cursor.close(bounded_size(request, "radius", 1000, 10000)).first;
    auto result = snapshot(cursor, graph, seed_offset);
    result["steps"] = Json::Value(Json::arrayValue);
    if (request.isMember("steps")) {
        if (!request["steps"].isArray() || request["steps"].size() > 100)
            throw std::invalid_argument("steps must be an array of at most 100 choices");
        for (const auto &step : request["steps"]) {
            if (!step.isObject() || !step["base"].isString() || step["base"].asString().size() != 1
                    || (step["side"] != "left" && step["side"] != "right"))
                throw std::invalid_argument("Each step needs side (left/right) and one base");
            bool extended = cursor.extend(step["side"] == "left" ? ContextCursor::Side::LEFT
                                                                  : ContextCursor::Side::RIGHT,
                                          step["base"].asString()[0]);
            if (extended) {
                seed_offset += step["side"] == "left";
                seed_offset += cursor.close(bounded_size(step, "radius", 1000, 10000)).first;
            }
            auto value = snapshot(cursor, graph, seed_offset);
            value["extended"] = extended;
            result["steps"].append(value);
        }
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
