#include <random>

#include <benchmark/benchmark.h>

#include "graph/representation/succinct/boss.hpp"


const std::string filename = "./benchmark_graph.dbg";

constexpr uint64_t NUM_DISTINCT_INDEXES = 1 << 24;


namespace mg {
namespace bm {

void load_graph(benchmark::State &state, BOSS *boss) {
    if (!std::getenv("GRAPH")) {
        state.SkipWithError("Set environment variable GRAPH");
    } else if (!boss->load(std::getenv("GRAPH"))) {
        state.SkipWithError((std::string("Can't load the graph from ")
                                + std::getenv("GRAPH")).c_str());
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


#define DEFINE_BOSS_BENCHMARK(NAME, OPERATION, VECTOR, SIZE, ...) \
static void BM_BOSS_##NAME(benchmark::State& state) { \
    BOSS boss; \
    load_graph(state, &boss); \
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
DEFINE_BOSS_BENCHMARK(rank_last,           rank_last,           get_last, size);
DEFINE_BOSS_BENCHMARK(select_last,         select_last,         get_last, num_set_bits);
DEFINE_BOSS_BENCHMARK(bwd,                 bwd,                 get_W,    size);


static void BM_BOSS_get_W_and_fwd(benchmark::State &state) {
    BOSS boss;
    load_graph(state, &boss);

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
    BOSS boss;
    load_graph(state, &boss);

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
    BOSS boss;
    load_graph(state, &boss);

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
    BOSS boss;
    load_graph(state, &boss);

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


} // namespace bm
} // namespace mg
