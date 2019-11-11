#include "partitionings.hpp"

#include <ips4o.hpp>
#include <progress_bar.hpp>

#include "algorithms.hpp"

const uint64_t kNumRowsSampled = 1'000'000;


std::vector<sdsl::bit_vector>
get_submatrix(const BRWTBottomUpBuilder::VectorsPtr &columns,
              const std::vector<uint64_t> &row_indexes,
              size_t num_threads) {
    assert(std::is_sorted(row_indexes.begin(), row_indexes.end()));

    if (!columns.size())
        return {};

    assert(row_indexes.size() <= columns[0]->size());

    std::vector<sdsl::bit_vector> submatrix(columns.size());

    ProgressBar progress_bar(columns.size(), "Subsampling",
                             std::cerr, !utils::get_verbose());

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < columns.size(); ++i) {
        submatrix[i] = subvector(*columns[i], row_indexes);

#ifndef NDEBUG
        for (size_t j = 0; j < row_indexes.size(); ++j) {
            assert(submatrix[i][j] == (*columns[i])[row_indexes[j]]);
        }
#endif
        ++progress_bar;
    }

    return submatrix;
}

// returns shrinked columns
std::vector<sdsl::bit_vector>
random_submatrix(const BRWTBottomUpBuilder::VectorsPtr &columns,
                 uint64_t num_rows_sampled,
                 int seed,
                 size_t num_threads) {
    if (!columns.size())
        return {};

    std::mt19937 gen;

    if (seed)
        gen.seed(seed);

    auto indexes = utils::sample_indexes(columns[0]->size(),
                                         num_rows_sampled,
                                         gen);
    // sort indexes
    std::sort(indexes.begin(), indexes.end());
    // check if indexes are sampled without replacement
    assert(std::unique(indexes.begin(), indexes.end()) == indexes.end());

    return get_submatrix(columns, indexes, num_threads);
}


// Partitionings for BRWT

// input: columns
// output: partition, for instance -- a set of column pairs
std::vector<BRWTBottomUpBuilder::Column>
inverted_arrangement(const BRWTBottomUpBuilder::VectorsPtr &vectors) {
    auto init_arrangement
            = utils::arange<BRWTBottomUpBuilder::Column>(0, vectors.size());

    return { init_arrangement.rbegin(), init_arrangement.rend() };
}

std::vector<std::vector<double>>
correlation_similarity(const std::vector<sdsl::bit_vector> &cols,
                       size_t num_threads) {
    std::vector<std::vector<double>> similarities(cols.size());

    for (size_t j = 1; j < cols.size(); ++j) {
        similarities[j] = std::vector<double>(j, 0.);
    }

    ProgressBar progress_bar(cols.size() * (cols.size() - 1) / 2, "Computing correlations",
                             std::cerr, !utils::get_verbose());

    #pragma omp parallel for num_threads(num_threads) collapse(2) schedule(static, 5)
    for (size_t j = 1; j < cols.size(); ++j) {
        for (size_t k = 0; k < cols.size(); ++k) {
            if (k >= j)
                continue;

            similarities[j][k] = inner_prod(cols[j], cols[k]);
            ++progress_bar;
        }
    }

    return similarities;
}

std::vector<std::vector<double>>
jaccard_similarity(const std::vector<sdsl::bit_vector> &cols,
                   size_t num_threads) {
    std::vector<uint64_t> num_set_bits(cols.size(), 0);

    #pragma omp parallel for num_threads(num_threads)
    for (size_t j = 0; j < cols.size(); ++j) {
        num_set_bits[j] = sdsl::util::cnt_one_bits(cols[j]);
    }

    auto similarities = correlation_similarity(cols, num_threads);

    #pragma omp parallel for num_threads(num_threads) collapse(2) schedule(static, 5)
    for (size_t j = 0; j < cols.size(); ++j) {
        for (size_t k = 0; k < cols.size(); ++k) {
            if (k >= j)
                continue;

            similarities[j][k] /= (num_set_bits[j]
                                    + num_set_bits[k]
                                    - similarities[j][k]);
            std::cout << similarities[j][k] << std::endl;
        }
    }

    return similarities;
}

// For each vector j return similarities with vectors 0, ..., j-1
std::vector<std::vector<double>>
estimate_similarities(const BRWTBottomUpBuilder::VectorsPtr &vectors,
                      size_t num_threads) {
    if (!vectors.size())
        return {};

    uint64_t num_sampled_rows = std::min(kNumRowsSampled, vectors[0]->size());

    return correlation_similarity(
        random_submatrix(vectors, num_sampled_rows, 1, num_threads),
        num_threads
    );
}

template <typename T>
inline T dist(T first, T second) {
    return first > second
                ? first - second
                : second - first;
}

// input: columns
// output: partition, for instance -- a set of column pairs
BRWTBottomUpBuilder::Partition
parallel_binary_grouping_greedy(const BRWTBottomUpBuilder::VectorsPtr &columns,
                                size_t num_threads) {
    if (!columns.size())
        return {};

    auto similarities = estimate_similarities(columns, num_threads);

    ProgressBar progress_bar(columns.size() * (columns.size() - 1), "Clustering",
                             std::cerr, !utils::get_verbose());

    std::vector<std::tuple<size_t, size_t, uint64_t>> candidates;
    candidates.reserve(columns.size() * (columns.size() - 1) / 2);

    for (size_t j = 1; j < similarities.size(); ++j) {
        for (size_t k = 0; k < j; ++k) {
            candidates.emplace_back(j, k, similarities[j][k]);
            ++progress_bar;
        }
    }

    // pick either a pair of the most similar columns,
    // or pair closest in the initial arrangement
    ips4o::parallel::sort(candidates.begin(), candidates.end(),
        [](const auto &first_pair, const auto &second_pair) {
              return std::get<2>(first_pair) > std::get<2>(second_pair)
                || (std::get<2>(first_pair) == std::get<2>(second_pair)
                        && dist(std::get<0>(first_pair), std::get<1>(first_pair))
                            < dist(std::get<0>(second_pair), std::get<1>(second_pair)));
        },
        num_threads
    );

    BRWTBottomUpBuilder::Partition partition;
    partition.reserve((columns.size() + 1) / 2);

    std::vector<bool> matched(columns.size(), false);

    for (const auto &next_candidate : candidates) {
        auto i = std::get<0>(next_candidate);
        auto j = std::get<1>(next_candidate);
        if (!matched[i] && !matched[j]) {
            matched[i] = matched[j] = true;
            partition.push_back({ i, j });
        }
        ++progress_bar;
    }

    for (size_t i = 0; i < columns.size(); ++i) {
        if (!matched[i])
            partition.push_back({ i });
    }

    return partition;
}

BRWTBottomUpBuilder::Partition
binary_grouping_greedy(const BRWTBottomUpBuilder::VectorsPtr &columns) {
    return parallel_binary_grouping_greedy(columns, 1);
}

std::function<BRWTBottomUpBuilder::Partition(const BRWTBottomUpBuilder::VectorsPtr &)>
get_parallel_binary_grouping_greedy(size_t num_threads) {
    return [num_threads](const BRWTBottomUpBuilder::VectorsPtr &columns) {
        return parallel_binary_grouping_greedy(columns, num_threads);
    };
}
