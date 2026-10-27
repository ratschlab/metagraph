#include "traverse_server.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <filesystem>
#include <memory>
#include <set>

#include "config/config.hpp"
#include "json_helpers.hpp"
#include "traverse_attempts.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/traversal/label_oracle.hpp"

namespace mtg::cli {

GraphListEntry parse_graph_list_line(const std::string &text, size_t line) {
    std::vector<std::string> columns;
    size_t from = 0;
    while (true) {
        const size_t comma = text.find(',', from);
        columns.push_back(text.substr(from, comma == std::string::npos ? std::string::npos
                                                                       : comma - from));
        if (comma == std::string::npos)
            break;
        from = comma + 1;
    }
    auto bad = [&](const std::string &what) {
        return std::invalid_argument("line " + std::to_string(line) + " of the graph list ('"
                                     + text + "'): " + what);
    };
    if (columns.size() < 3)
        throw bad("expected name,graph_path,annotation_path[,manifest_path[,index_ns]]");
    if (columns.size() > 5) {
        throw bad(std::to_string(columns.size()) + " columns, at most five are read "
                  "(name,graph_path,annotation_path,manifest_path,index_ns): a column the "
                  "server would not read is refused rather than ignored");
    }
    GraphListEntry e;
    e.line = line;
    e.name = columns[0];
    e.graph_path = columns[1];
    e.annotation_path = columns[2];
    if (columns.size() > 3)
        e.manifest_path = columns[3];
    if (columns.size() > 4) {
        e.index_ns = columns[4];
        // the name is an MGT token and goes into file names (as --index-name)
        bool ok = !e.index_ns.empty();
        for (unsigned char c : e.index_ns) {
            ok &= std::isalnum(c) || c == '.' || c == '_' || c == '-';
        }
        if (!e.index_ns.empty() && !ok)
            throw bad("index_ns '" + e.index_ns + "' does not match [A-Za-z0-9._-]+");
    }
    return e;
}

std::map<std::pair<std::string, std::string>, std::pair<std::string, std::string>>
graph_list_identities(const std::vector<GraphListEntry> &entries,
                      const std::function<std::string(const GraphListEntry &)> &fingerprint,
                      const std::function<std::vector<std::string>(const GraphListEntry &)> &inventory) {
    using Pair = std::pair<std::string, std::string>;
    struct Stated {
        std::string value;
        size_t line = 0;
    };
    std::map<Pair, std::pair<Stated, Stated>> stated;     // (ns, fp) and the lines stating them
    std::map<std::pair<Pair, std::string>, std::string> digests;   // (pair, manifest) -> fp
    auto agree = [](Stated *have, const std::string &value, size_t line, const char *what,
                    const Pair &pair) {
        if (value.empty())
            return;
        if (have->value.empty()) {
            *have = { value, line };
            return;
        }
        if (have->value != value) {
            throw std::invalid_argument(
                    std::string("lines ") + std::to_string(have->line) + " and "
                    + std::to_string(line) + " of the graph list state different " + what
                    + " for the same index (" + pair.first + ", " + pair.second + "): '"
                    + have->value + "' and '" + value + "'; one index has one identity");
        }
    };
    // One index is one set of loaded files, however its paths are spelled: the lines are
    // grouped by the real paths of the pair's COMPLETE loader inventory (index_load_inventory:
    // the main files and every identity sidecar, derived from the spelling as the loaders derive
    // them; the graph's derived data, its mask and Bloom filter, is not part of an identity),
    // so that two spellings of one index agree, and two pairs whose main files are symlinks to
    // the same files but whose sidecars differ are two indexes, each validated against its
    // manifest (grouped by the main files' real paths, the second pair's own .seqs would never
    // be checked and both would state one index_fp)
    auto real = [](const std::string &path) {
        std::error_code ec;
        const auto absolute = std::filesystem::absolute(path, ec);
        if (ec)
            return path;
        const auto canonical = std::filesystem::weakly_canonical(absolute, ec);
        return ec ? absolute.lexically_normal().string() : canonical.string();
    };
    // the bundle of a listed pair: its files' real paths in the inventory's order, and the
    // real main files (for messages)
    using Bundle = std::vector<std::string>;
    std::map<Pair, Bundle> bundle_of;                      // the pair as listed -> its bundle
    std::map<Bundle, Pair> pair_of;                        // a bundle -> the first pair listing it
    for (const GraphListEntry &e : entries) {
        const Pair listed { e.graph_path, e.annotation_path };
        auto [b, fresh_pair] = bundle_of.try_emplace(listed);
        if (fresh_pair) {
            const std::vector<std::string> files
                = inventory ? inventory(e) : index_bundle_files(e.graph_path, e.annotation_path);
            for (const std::string &f : files) {
                b->second.push_back(real(f));
            }
            pair_of.try_emplace(b->second, listed);
        }
        const Bundle &bundle = b->second;
        // keyed by the bundle's first listing: the messages name a pair as listed
        const Pair pair = pair_of.at(bundle);
        auto &[ns, fp] = stated[pair];
        agree(&ns, e.index_ns, e.line, "index_ns", listed);
        if (e.manifest_path.empty())
            continue;
        // every distinct bundle is validated against each manifest a line gives for it
        auto [it, fresh] = digests.try_emplace({ pair, e.manifest_path });
        if (fresh)
            it->second = fingerprint(e);
        agree(&fp, it->second, e.line, "index_fp (manifest digests)", listed);
    }
    // Two different indexes never state one fingerprint: the server checks sizes, not
    // contents, so one manifest whose files have the base names and sizes of another pair's
    // (two annotations of one graph with swapped memberships, written by one tool run) would
    // let a client take one index for the other. Copies of one index under
    // two paths are refused too (telling them apart would need hashing): list one path
    std::map<std::string, std::pair<Pair, size_t>> by_fp;
    for (const auto &[pair, s] : stated) {
        if (s.second.value.empty())
            continue;
        auto [it, fresh] = by_fp.try_emplace(s.second.value, pair, s.second.line);
        if (!fresh) {
            throw std::invalid_argument(
                    "lines " + std::to_string(it->second.second) + " and "
                    + std::to_string(s.second.line) + " of the graph list state one index_fp "
                    + s.second.value + " for two different indexes (" + it->second.first.first
                    + ", " + it->second.first.second + ") and (" + pair.first + ", "
                    + pair.second + "): a manifest states the identity of one graph with one "
                    "annotation; give each pair its own manifest, or list copies of one "
                    "index by one path");
        }
    }
    std::map<Pair, std::pair<std::string, std::string>> out;
    for (const auto &[listed, bundle] : bundle_of) {
        const auto &s = stated.at(pair_of.at(bundle));
        out[listed] = { s.first.value, s.second.value };
    }
    return out;
}



GraphSelection traverse_graph_selection(const Json::Value &request, bool multi_graph) {
    GraphSelection selection;
    if (!request.isObject())
        return selection;
    if (!multi_graph) {
        if (request.isMember("graph") || request.isMember("graph_path")) {
            throw std::invalid_argument("Bad request: this server hosts a single graph; "
                                        "remove the 'graph' / 'graph_path' field");
        }
        if (request.isMember("graphs")) {
            throw std::invalid_argument("Bad request: this server hosts a single graph; "
                                        "remove the 'graphs' field");
        }
        return selection;
    }
    if (request.isMember("graphs")) {
        // /search's field, one name: a traversal reads one graph
        const Json::Value &graphs = request["graphs"];
        if (request.isMember("graph")) {
            throw std::invalid_argument("Bad request: give the graph as 'graph' or as "
                                        "'graphs', not both");
        }
        if (!graphs.isArray() || graphs.size() != 1 || !graphs[0].isString()) {
            throw std::invalid_argument("Bad request: 'graphs' names the one graph a traversal "
                                        "reads: expected [name]");
        }
        selection.name = graphs[0].asString();
    } else {
        if (!request.isMember("graph") || !request["graph"].isString())
            throw std::invalid_argument("Bad request: 'graph' (index name) is required in "
                                        "multi-graph mode");
        selection.name = request["graph"].asString();
    }
    if (request.isMember("graph_path")) {
        if (!request["graph_path"].isString())
            throw InvalidRequest("request.graph_path: expected a string");
        selection.graph_path = request["graph_path"].asString();
    }
    return selection;
}


const GraphPair& select_traverse_pair(const std::string &name,
                                      const std::optional<std::string> &graph_path,
                                      std::span<const GraphPair> pairs) {
    const GraphPair *found = nullptr;
    for (const auto &pair : pairs) {
        if (graph_path && pair.first != *graph_path)
            continue;
        if (found && *found != pair)
            throw InvalidRequest("Index '" + name + "' is ambiguous; select one graph_path "
                                 "and list different annotations under separate names");
        found = &pair;
    }
    if (!found)
        throw InvalidRequest("No graph matching the requested path in index '" + name + "'");
    return *found;
}

void add_traversal_routes(HttpServer &server, const Config &config,
                          TraversalLookup lookup, std::function<Json::Value()> catalog,
                          std::atomic<size_t> &request_ids,
                          const std::atomic<bool> &stopping) {
    AttemptSettings settings;
    settings.allowance_ms = config.traverse_attempt_allowance_ms;
    settings.retention_s = config.traverse_attempt_retention_s;
    settings.retention_count = config.traverse_attempt_retention;
    settings.tombstone_max_s = config.traverse_attempt_tombstone_max_s;
    settings.clock_skew_ms = config.traverse_clock_skew_ms;
    settings.observed_window_s = config.traverse_observed_window_s;
    settings.delivery_compress_mbps = config.traverse_delivery_compress_mbps;
    settings.delivery_build_mbps = config.traverse_delivery_build_mbps;
    settings.delivery_stop_ms = config.traverse_chunk_target_ms + 950.0;
    auto attempts = std::make_shared<AttemptRegistry>(settings);
    const ResolveTimeLimits resolve_time {
        config.traverse_max_time_ms > 0
            ? std::min(config.traverse_max_time_ms, 899'000.0) : 899'000.0,
        kResolveFinalizeMs
    };
    auto gone_of = [&stopping](const std::shared_ptr<HttpServer::Request> &request) {
        return [weak = std::weak_ptr<HttpServer::Request>(request), &stopping]() {
            auto r = weak.lock();
            return stopping.load() || !r || client_gone(*r);
        };
    };
    auto parse = [](const std::string &content) {
        return parse_strict_json<RefuseRequest>(content, "Bad json received: ");
    };
    auto read_index = [lookup](const Json::Value &json) {
        if (json.isMember("in_ram") && (!json["in_ram"].isBool() || json["in_ram"].asBool()))
            throw InvalidRequest("Traversal uses the server's loaded index; in_ram must be false "
                                 "or omitted. Start the server without --mmap to load into RAM.");
        return lookup(json);
    };
    auto put_contract = [=, &config](Json::Value &caps) {
        caps["algorithm_version"] = kTraverseAlgorithmVersion;
        caps["feature_level"] = kTraverseFeatureLevel;
        caps["attempts"] = attempts->capabilities_json();
        caps["resolve"] = resolve_capabilities_json(resolve_time);
        caps["content_encodings"] = Json::arrayValue;
        caps["content_encodings"].append("gzip");
        caps["content_encodings"].append("deflate");
        caps["compression_level"] = config.traverse_compression_level;
        auto &deadline = caps["deadline_check"];
        deadline["chunk_target_ms"] = Json::UInt64(config.traverse_chunk_target_ms);
        deadline["observed_window_s"] = Json::UInt64(config.traverse_observed_window_s);
        deadline["observed_max_uninterruptible_ms"]
            = Json::UInt64(attempts->observed_max_uninterruptible_ms());
        deadline["poll_stride"] = Json::UInt64(attempts->settings().poll_stride);
        caps["in_ram"]["available"] = false;
        caps["in_ram"]["reason"] = "Traversal uses the index already loaded by the server";
    };
    server.resource["^/capabilities$"]["GET"] = [=, &request_ids, &config](auto response, auto request) {
        process_request(response, request, request_ids++, [&](const std::string &) {
            auto caps = catalog();
            put_contract(caps);
            caps["features"] = Json::arrayValue;
            for (const auto *feature : { "search", "resolve", "traverse", "attempts" })
                caps["features"].append(feature);
            auto &routes = caps["routes"];
            routes["capabilities"] = "GET /capabilities";
            routes["search"] = "POST /search";
            routes["resolve"] = "POST /resolve";
            routes["traverse"] = "POST /traverse";
            routes["traverse_capabilities"] = config.fnames.empty()
                ? "GET /traverse/capabilities"
                : "GET /traverse/capabilities?graph={name}[&graph_path={path}]";
            routes["cancel"] = "POST /traverse/cancel";
            routes["attempt"] = "GET /traverse/attempt/{attempt_id}";
            routes["column_labels"] = "GET /column_labels";
            routes["stats"] = "GET /stats";
            if (config.fnames.empty()) {
                caps["features"].append("align");
                routes["align"] = "POST /align";
            }
            return caps;
        }, true);
    };
    // Keep registry and configuration alive through the server-owned route closures.
    server.resource["^/resolve$"]["POST"] = [=, &request_ids, &config](auto response, auto request) {
        const auto gone = gone_of(request);
        ResolveWriting writing(gone, std::chrono::milliseconds(attempts->settings().client_check_ms),
                               config.traverse_compression_level);
        process_request(response, request, request_ids++, [&](const std::string &content) {
            const auto json = parse(content);
            const auto index = read_index(json);
            try {
                auto out = process_resolve_request(json, index.graph, config.index_release,
                                                   config.resolve_max_query_bp, &index.identity,
                                                   gone, resolve_time, &writing.delivery);
                out["timing"]["prepare_ms"] = index.prepare_ms;
                return out;
            } catch (const graph::traversal::AttemptAborted &e) {
                throw ClientGone(e.what());
            } catch (const ResolveDeadline &e) {
                throw HttpError(503, resolve_deadline_body(e));
            }
        }, true, &writing.control);
    };

    server.resource["^/traverse$"]["POST"] = [=, &request_ids, &config](auto response, auto request) {
        const size_t request_id = request_ids++;
        auto attempt = std::make_shared<Attempt>(request_id, request->header_read_time,
                attempts->settings(), attempts->server_instance(), true, gone_of(request));
        attempt->set_measured(attempts->measured());
        bool registered = false;
        auto with_usage = [&](int status, const std::string &message, const std::string &reason) {
            Json::Value body;
            body["error"] = message;
            body["usage"] = attempt->usage_json(reason);
            return HttpError(status, std::move(body));
        };
        ResultTexts texts;
        ResponseControl control;
        control.compression_level = config.traverse_compression_level;
        control.write = [&](const Json::Value &envelope, const std::function<void()> &check) {
            double gap = 0;
            auto text = texts.active
                ? assemble_traverse_response(envelope, std::move(texts.texts), check, &gap)
                : json_text(envelope, true, check, &gap);
            attempt->note_delivery_gap_ms(gap);
            return text;
        };
        control.check = [&]() {
            try {
                attempt->check_delivery();
            } catch (const graph::traversal::AttemptAborted &e) {
                throw ClientGone(e.what());
            } catch (const AttemptAtBound &e) {
                throw with_usage(503, e.what(), "deadline");
            }
        };
        control.on_compressed = [&](size_t bytes, double seconds) {
            if (bytes >= kMeasuredTextBytes && seconds > 0)
                attempts->note_compress_rate(bytes / seconds / 1e6);
        };
        control.on_written = [&](int status, std::optional<size_t> bytes) {
            attempts->note_uninterruptible(attempt->max_uninterruptible_ms());
            attempts->note_uninterruptible(attempt->max_delivery_gap_ms());
            attempts->note_build_rate(attempt->own_build_mbps());
            attempts->note_account_per_text_byte(attempt->delivery_detail(),
                                                 attempt->own_account_per_text_byte());
            attempts->note_stop_latency(attempt->own_stop_ms());
            if (registered) {
                attempts->finish(attempt, !status ? "client_gone" : status == 503 ? "deadline"
                    : status >= 400 ? "error" : attempt->walk_reason(), status, bytes);
            }
        };
        process_request(response, request, request_id, [&](const std::string &content) {
            const auto json = parse(content);
            attempt->set_ids(attempt_ids(json));
            if (attempt->managed()) {
                if (auto refused = attempts->start(attempt)) {
                    if (refused->instance_mismatch || refused->expired)
                        throw HttpError(409, std::move(refused->body));
                    Json::Value body;
                    body["error"] = "attempt_id is running, retained or cancelled; request not started";
                    body["attempt"] = std::move(refused->body);
                    throw HttpError(409, std::move(body));
                }
                registered = true;
            } else if (not_after_passed(attempt->ids(), attempts->now_ms())) {
                throw HttpError(409, expired_json(attempt->ids(), attempts->now_ms(),
                                                  attempts->server_instance()));
            }
            try {
                const auto index = read_index(json);
                TraverseLimits limits;
                limits.max_time_ms = config.traverse_max_time_ms;
                limits.max_seeds = config.traverse_max_seeds;
                limits.max_seed_bp = config.traverse_max_seed_bp;
                limits.max_seed_labels = config.traverse_max_seed_labels;
                limits.max_memory_mb = config.traverse_max_memory_mb;
                limits.max_work_units = config.traverse_max_work_units;
                limits.chunk_target_ms = config.traverse_chunk_target_ms;
                limits.path_cache_bytes = config.traverse_path_cache_mb << 20;
                limits.not_in_graph_per_seed = !config.fnames.empty();
                auto out = process_traverse_request(json, index.graph, config.index_release, limits,
                                                    &index.identity, attempt.get(), &texts);
                out["timing"]["prepare_ms"] = index.prepare_ms;
                return out;
            } catch (const graph::traversal::AttemptAborted &e) {
                throw ClientGone(e.what());
            } catch (const AttemptAtBound &e) {
                throw with_usage(503, e.what(), "deadline");
            } catch (const CurrentlyInitializingError &) {
                if (registered)
                    throw with_usage(503, "Server is currently initializing", "error");
                throw;
            } catch (const std::exception &e) {
                if (registered)
                    throw with_usage(400, e.what(), "error");
                throw;
            }
        }, true, &control);
    };

    server.resource["^/traverse/cancel$"]["POST"] = [=, &request_ids](auto response, auto request) {
        process_request(response, request, request_ids++, [&](const std::string &content) {
            const auto cancel = parse_cancel_request(parse(content));
            auto [status, body] = attempts->cancel(cancel.attempt_id, cancel.wait_ms, cancel.not_after_ms);
            if (status != 200)
                throw HttpError(status, std::move(body));
            return body;
        }, true);
    };
    server.resource["^/traverse/attempt/([^/]+)$"]["GET"] = [=, &request_ids](auto response, auto request) {
        process_request(response, request, request_ids++, [&](const std::string &) {
            const std::string id = request->path_match[1].str();
            if (!valid_attempt_id(id))
                throw InvalidRequest("Invalid attempt_id");
            auto [status, body] = attempts->state(id);
            if (status != 200)
                throw HttpError(status, std::move(body));
            return body;
        }, true);
    };
    server.resource["^/traverse/capabilities$"]["GET"] = [=, &request_ids, &config](auto response, auto request) {
        process_request(response, request, request_ids++, [&](const std::string &) {
            Json::Value query(Json::objectValue);
            for (const auto &[key, value] : request->parse_query_string()) {
                if ((key != "graph" && key != "graph_path") || query.isMember(key))
                    throw InvalidRequest("Unknown or repeated capability parameter: " + key);
                query[key] = value;
            }
            const auto index = read_index(query);
            graph::traversal::LabelOracle oracle(index.graph);
            auto caps = capabilities_to_json(oracle, config.index_release, &index.identity);
            put_contract(caps);
            caps["coordinates"] = coordinates_capabilities_json(oracle);
            caps["max_time_ms"] = config.traverse_max_time_ms;
            caps["max_seeds"] = Json::UInt64(config.traverse_max_seeds);
            caps["max_seed_bp"] = Json::UInt64(config.traverse_max_seed_bp);
            caps["max_seed_labels"] = Json::UInt64(config.traverse_max_seed_labels);
            caps["max_query_bp"] = Json::UInt64(config.resolve_max_query_bp);
            caps["max_memory_mb"] = Json::UInt64(config.traverse_max_memory_mb);
            caps["max_work_units"] = Json::UInt64(config.traverse_max_work_units);
            caps["memory_bound"] = "soft";
            caps["work_check_interval"] = Json::UInt64(graph::traversal::kWorkCheckInterval);
            caps["decode_cache"]["path_cache_mb"] = Json::UInt64(config.traverse_path_cache_mb);
            if (query.isMember("graph"))
                caps["graph"] = query["graph"];
            if (query.isMember("graph_path"))
                caps["graph_path"] = query["graph_path"];
            return caps;
        }, true);
    };
}

} // namespace mtg::cli
