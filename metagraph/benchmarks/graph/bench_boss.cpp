#include <random>

#include <benchmark/benchmark.h>

#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/dbg_succinct_cached.hpp"
#include "graph/representation/canonical_dbg.hpp"


namespace {

using namespace mtg::graph;

using mtg::graph::boss::BOSS;

const std::string filename = "./benchmark_graph.dbg";

constexpr uint64_t NUM_DISTINCT_INDEXES = 1 << 21;
constexpr uint64_t PATH_SIZE = 1 << 7;
constexpr uint64_t CACHE_SIZE = 1 << 10;

template <typename Graph = DBGSuccinct, typename OutGraph = DBGSuccinct>
std::shared_ptr<OutGraph> load_graph(benchmark::State &state) {
    auto graph = std::make_shared<DBGSuccinct>(2);
    if (!std::getenv("GRAPH")) {
        state.SkipWithError("Set environment variable GRAPH");
    } else if (!graph->load(std::getenv("GRAPH"))) {
        state.SkipWithError((std::string("Can't load the graph from ")
                                + std::getenv("GRAPH")).c_str());
    }

    if constexpr(std::is_same_v<OutGraph, DBGSuccinct>) {
        return graph;
    } else {
        static_assert(std::is_same_v<OutGraph, DeBruijnGraph>);
        if constexpr(std::is_same_v<Graph, CanonicalDBG>) {
            return std::make_shared<CanonicalDBG>(graph);
        } else if constexpr(std::is_same_v<Graph, DBGSuccinctCached>) {
            std::shared_ptr<DeBruijnGraph> wrapped_graph
                = make_cached_dbgsuccinct(graph, CACHE_SIZE);
            if (graph->get_mode() == DeBruijnGraph::PRIMARY)
                wrapped_graph = std::make_shared<CanonicalDBG>(wrapped_graph);

            return wrapped_graph;
        }
    }
}

// generate a deterministic sequence of pseudo-random numbers
std::vector<uint64_t> random_numbers(size_t size, uint64_t min, uint64_t max) {
    std::mt19937 gen(32);
    std::uniform_int_distribution<uint64_t> dis(min, max);

    std::vector<uint64_t> numbers(size);
    for (uint64_t &number : numbers) {
        number = dis(gen);
    }
    return numbers;
}

std::vector<uint64_t> random_traversal_numbers(const DeBruijnGraph &graph,
                                               size_t size,
                                               size_t path_size) {
    std::mt19937 gen(32);
    std::uniform_int_distribution<uint64_t> dis(1, graph.num_nodes());

    std::vector<uint64_t> numbers;
    numbers.reserve(size);

    while (numbers.size() < size) {
        numbers.push_back(dis(gen));
        for (size_t j = 1; j < path_size && numbers.size() < size; ++j) {
            bool found = false;
            graph.adjacent_outgoing_nodes(numbers.back(), [&](auto next) {
                found = true;
                numbers.push_back(next);
            });
            if (!found)
                break;
        }
    }

    if (numbers.size() > size)
        throw std::runtime_error("Number generation failed");

    return numbers;
}


#define DEFINE_BOSS_BENCHMARK(NAME, OPERATION, VECTOR, SIZE, ...) \
static void BM_BOSS_##NAME(benchmark::State& state) { \
    auto graph = load_graph(state); \
    const BOSS &boss = graph->get_boss(); \
 \
    auto indexes = random_numbers(NUM_DISTINCT_INDEXES, 1, boss.VECTOR().SIZE() - 1); \
    size_t i = 0; \
    for (auto _ : state) { \
        uint64_t edge = indexes[i++ % NUM_DISTINCT_INDEXES]; \
        benchmark::DoNotOptimize(boss.OPERATION(edge, ##__VA_ARGS__)); \
    } \
} \
BENCHMARK(BM_BOSS_##NAME) -> Unit(benchmark::kMicrosecond); \

DEFINE_BOSS_BENCHMARK(get_node_last_value, get_node_last_value, get_W,    size);
DEFINE_BOSS_BENCHMARK(get_W,               get_W,               get_W,    size);
DEFINE_BOSS_BENCHMARK(rank_W_$,            rank_W,              get_W,    size, 0);
DEFINE_BOSS_BENCHMARK(rank_W_A,            rank_W,              get_W,    size, 1);
DEFINE_BOSS_BENCHMARK(rank_W_C,            rank_W,              get_W,    size, 2);
DEFINE_BOSS_BENCHMARK(rank_W_G,            rank_W,              get_W,    size, 3);
DEFINE_BOSS_BENCHMARK(rank_W_T,            rank_W,              get_W,    size, 4);
DEFINE_BOSS_BENCHMARK(get_last,            get_last,            get_last, size);
DEFINE_BOSS_BENCHMARK(rank_last,           rank_last,           get_last, size);
DEFINE_BOSS_BENCHMARK(select_last,         select_last,         get_last, num_set_bits);
DEFINE_BOSS_BENCHMARK(pred_last,           pred_last,           get_last, size);
DEFINE_BOSS_BENCHMARK(succ_last,           succ_last,           get_last, size);
DEFINE_BOSS_BENCHMARK(bwd,                 bwd,                 get_W,    size);


#define DEFINE_BOSS_CACHED_CYCLE_BENCHMARK(NAME, OPERATION, GRAPH_TYPE, OUT_TYPE, ...) \
static void BM_BOSS_##NAME(benchmark::State& state) { \
    auto graph = load_graph<GRAPH_TYPE, OUT_TYPE>(state); \
 \
    auto indexes = random_numbers(NUM_DISTINCT_INDEXES, 1, graph->num_nodes()); \
    size_t i = 0; \
    for (auto _ : state) { \
        DeBruijnGraph::node_index node = indexes[i++ % NUM_DISTINCT_INDEXES]; \
        benchmark::DoNotOptimize(graph->OPERATION(node, ##__VA_ARGS__)); \
    } \
} \
BENCHMARK(BM_BOSS_##NAME) -> Unit(benchmark::kMicrosecond); \

#define DEFINE_BOSS_CACHED_PATH_BENCHMARK(NAME, OPERATION, GRAPH_TYPE, OUT_TYPE) \
static void BM_BOSS_##NAME(benchmark::State& state) { \
    auto graph = load_graph<GRAPH_TYPE, OUT_TYPE>(state); \
 \
    size_t size = NUM_DISTINCT_INDEXES >> 2; \
    auto indexes = random_traversal_numbers(*graph, size, PATH_SIZE); \
    size_t i = 0; \
    for (auto _ : state) { \
        DeBruijnGraph::node_index node = indexes[i++ % size]; \
        graph->OPERATION(node, [](auto node, char c) { \
            benchmark::DoNotOptimize(node); \
            benchmark::DoNotOptimize(c); \
        }); \
    } \
} \
BENCHMARK(BM_BOSS_##NAME) -> Unit(benchmark::kMicrosecond); \

DEFINE_BOSS_CACHED_CYCLE_BENCHMARK(get_node_sequence_uncached_distinct, get_node_sequence, CanonicalDBG, DeBruijnGraph);
DEFINE_BOSS_CACHED_CYCLE_BENCHMARK(get_node_sequence_cached_distinct, get_node_sequence, DBGSuccinctCached, DeBruijnGraph);
DEFINE_BOSS_CACHED_PATH_BENCHMARK(call_outgoing_kmers_uncached_path, call_outgoing_kmers, CanonicalDBG, DeBruijnGraph);
DEFINE_BOSS_CACHED_PATH_BENCHMARK(call_outgoing_kmers_cached_path, call_outgoing_kmers, DBGSuccinctCached, DeBruijnGraph);


static void BM_BOSS_get_W_and_fwd(benchmark::State &state) {
    auto graph = load_graph(state);
    const BOSS &boss = graph->get_boss();

    std::mt19937 gen(32);
    std::uniform_int_distribution<uint64_t> dis(1, boss.get_W().size() - 1);

    std::vector<uint64_t> edges;
    edges.reserve(NUM_DISTINCT_INDEXES);
    while (edges.size() < NUM_DISTINCT_INDEXES) {
        uint64_t edge = dis(gen);
        if (boss.get_W().size() <= 2 || boss.get_W(edge))
            edges.push_back(edge);
    }

    size_t i = 0;
    for (auto _ : state) {
        uint64_t edge = edges[i++ % NUM_DISTINCT_INDEXES];
        benchmark::DoNotOptimize(boss.fwd(edge, boss.get_W(edge) % boss.alph_size));
    }
}
BENCHMARK(BM_BOSS_get_W_and_fwd) -> Unit(benchmark::kMicrosecond);


static void BM_BOSS_fwd(benchmark::State &state) {
    auto graph = load_graph(state);
    const BOSS &boss = graph->get_boss();

    std::mt19937 gen(32);
    std::uniform_int_distribution<uint64_t> dis(1, boss.get_W().size() - 1);

    std::vector<std::pair<uint64_t, BOSS::TAlphabet>> edges;
    edges.reserve(NUM_DISTINCT_INDEXES);
    while (edges.size() < NUM_DISTINCT_INDEXES) {
        uint64_t edge = dis(gen);
        if (boss.get_W().size() <= 2 || boss.get_W(edge))
            edges.emplace_back(edge, boss.get_W(edge) % boss.alph_size);
    }

    size_t i = 0;
    for (auto _ : state) {
        auto edge = edges[i++ % NUM_DISTINCT_INDEXES];
        benchmark::DoNotOptimize(boss.fwd(edge.first, edge.second));
    }
}
BENCHMARK(BM_BOSS_fwd) -> Unit(benchmark::kMicrosecond);


template <BOSS::TAlphabet edge_label>
static void BM_BOSS_pick_edge(benchmark::State &state) {
    auto graph = load_graph(state);
    const BOSS &boss = graph->get_boss();

    std::mt19937 gen(32);
    std::uniform_int_distribution<uint64_t> dis(1, boss.get_W().size() - 1);

    std::vector<uint64_t> edges;
    edges.reserve(NUM_DISTINCT_INDEXES);
    while (edges.size() < NUM_DISTINCT_INDEXES) {
        uint64_t edge = dis(gen);
        if (boss.get_last(edge))
            edges.push_back(edge);
    }

    size_t i = 0;
    for (auto _ : state) {
        uint64_t edge = edges[i++ % NUM_DISTINCT_INDEXES];
        benchmark::DoNotOptimize(boss.pick_edge(edge, edge_label));
    }
}
BENCHMARK_TEMPLATE(BM_BOSS_pick_edge, 0) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_BOSS_pick_edge, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_BOSS_pick_edge, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_BOSS_pick_edge, 3) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_BOSS_pick_edge, 4) -> Unit(benchmark::kMicrosecond);


template <BOSS::TAlphabet edge_label>
static void BM_BOSS_fwd_and_pick_edge(benchmark::State &state) {
    auto graph = load_graph(state);
    const BOSS &boss = graph->get_boss();

    std::mt19937 gen(32);
    std::uniform_int_distribution<uint64_t> dis(1, boss.get_W().size() - 1);

    std::vector<std::pair<uint64_t, BOSS::TAlphabet>> edges;
    edges.reserve(NUM_DISTINCT_INDEXES);
    while (edges.size() < NUM_DISTINCT_INDEXES) {
        uint64_t edge = dis(gen);
        if (boss.get_W().size() <= 2 || boss.get_W(edge))
            edges.emplace_back(edge, boss.get_W(edge) % boss.alph_size);
    }

    size_t i = 0;
    for (auto _ : state) {
        auto edge = edges[i++ % NUM_DISTINCT_INDEXES];
        auto next_edge = boss.fwd(edge.first, edge.second);
        benchmark::DoNotOptimize(boss.pick_edge(next_edge, edge_label));
    }
}
BENCHMARK_TEMPLATE(BM_BOSS_fwd_and_pick_edge, 0) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_BOSS_fwd_and_pick_edge, 1) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_BOSS_fwd_and_pick_edge, 2) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_BOSS_fwd_and_pick_edge, 3) -> Unit(benchmark::kMicrosecond);
BENCHMARK_TEMPLATE(BM_BOSS_fwd_and_pick_edge, 4) -> Unit(benchmark::kMicrosecond);


DEFINE_BOSS_BENCHMARK(pred_W_$,            pred_W,              get_W,    size, 0);
DEFINE_BOSS_BENCHMARK(pred_W_A,            pred_W,              get_W,    size, 1);
DEFINE_BOSS_BENCHMARK(pred_W_C,            pred_W,              get_W,    size, 2);
DEFINE_BOSS_BENCHMARK(pred_W_G,            pred_W,              get_W,    size, 3);
DEFINE_BOSS_BENCHMARK(pred_W_T,            pred_W,              get_W,    size, 4);
DEFINE_BOSS_BENCHMARK(pred_W_AAm,          pred_W,              get_W,    size, 1, 6);
DEFINE_BOSS_BENCHMARK(pred_W_CCm,          pred_W,              get_W,    size, 2, 7);
DEFINE_BOSS_BENCHMARK(pred_W_GGm,          pred_W,              get_W,    size, 3, 8);
DEFINE_BOSS_BENCHMARK(pred_W_TTm,          pred_W,              get_W,    size, 4, 9);

DEFINE_BOSS_BENCHMARK(succ_W_$,            succ_W,              get_W,    size, 0);
DEFINE_BOSS_BENCHMARK(succ_W_A,            succ_W,              get_W,    size, 1);
DEFINE_BOSS_BENCHMARK(succ_W_C,            succ_W,              get_W,    size, 2);
DEFINE_BOSS_BENCHMARK(succ_W_G,            succ_W,              get_W,    size, 3);
DEFINE_BOSS_BENCHMARK(succ_W_T,            succ_W,              get_W,    size, 4);
DEFINE_BOSS_BENCHMARK(succ_W_AAm,          succ_W,              get_W,    size, 1, 6);
DEFINE_BOSS_BENCHMARK(succ_W_CCm,          succ_W,              get_W,    size, 2, 7);
DEFINE_BOSS_BENCHMARK(succ_W_GGm,          succ_W,              get_W,    size, 3, 8);
DEFINE_BOSS_BENCHMARK(succ_W_TTm,          succ_W,              get_W,    size, 4, 9);


static void BM_BOSS_num_incoming_to_target(benchmark::State &state) {
    auto graph = load_graph(state);
    const BOSS &boss = graph->get_boss();

    std::mt19937 gen(32);
    std::uniform_int_distribution<uint64_t> dis(1, boss.get_W().size() - 1);

    std::vector<std::pair<uint64_t, BOSS::TAlphabet>> edges;
    edges.reserve(NUM_DISTINCT_INDEXES);
    while (edges.size() < NUM_DISTINCT_INDEXES) {
        uint64_t edge = boss.bwd(dis(gen));
        edges.emplace_back(edge, boss.get_W(edge));
    }

    size_t i = 0;
    for (auto _ : state) {
        auto edge = edges[i++ % NUM_DISTINCT_INDEXES];
        benchmark::DoNotOptimize(boss.num_incoming_to_target(edge.first, edge.second));
    }
}
BENCHMARK(BM_BOSS_num_incoming_to_target) -> Unit(benchmark::kMicrosecond);


DEFINE_BOSS_BENCHMARK(get_node_seq,  get_node_seq,              get_W,    size);

} // namespace
