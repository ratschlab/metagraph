#include "partitionings.hpp"

#include <ips4o.hpp>
#include <progress_bar.hpp>

#include "common/algorithms.hpp"
#include "common/vectors/vector_algorithm.hpp"

typedef std::vector<std::vector<uint64_t>> Partition;
typedef std::vector<const bit_vector *> VectorPtrs;


std::vector<sdsl::bit_vector>
get_submatrix(const VectorPtrs &columns,
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
random_submatrix(const VectorPtrs &columns,
                 uint64_t num_rows_sampled, int seed, size_t num_threads) {
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


// Partitionings for Multi-BRWT

// input: columns
// output: partition, for instance -- a set of column pairs
std::vector<uint64_t> inverted_arrangement(const VectorPtrs &vectors) {
    auto init_arrangement
            = utils::arange<uint64_t>(0, vectors.size());

    return { init_arrangement.rbegin(), init_arrangement.rend() };
}

std::vector<std::vector<uint64_t>>
correlation_similarity(const std::vector<sdsl::bit_vector> &cols,
                       size_t num_threads) {
    std::vector<std::vector<uint64_t>> similarities(cols.size());

    for (size_t j = 1; j < cols.size(); ++j) {
        similarities[j].assign(j, 0);
    }

    ProgressBar progress_bar(cols.size() * (cols.size() - 1) / 2, "Correlations",
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
jaccard_similarity(const std::vector<sdsl::bit_vector> &cols, size_t num_threads) {
    std::vector<std::vector<double>> similarities(cols.size());

    for (size_t j = 1; j < cols.size(); ++j) {
        similarities[j].assign(j, 0);
    }

    std::vector<uint64_t> num_set_bits(cols.size(), 0);

    #pragma omp parallel for num_threads(num_threads)
    for (size_t j = 0; j < cols.size(); ++j) {
        num_set_bits[j] = sdsl::util::cnt_one_bits(cols[j]);
    }

    ProgressBar progress_bar(cols.size() * (cols.size() - 1) / 2, "Jaccard",
                             std::cerr, !utils::get_verbose());

    #pragma omp parallel for num_threads(num_threads) collapse(2) schedule(static, 5)
    for (size_t j = 0; j < cols.size(); ++j) {
        for (size_t k = 0; k < cols.size(); ++k) {
            if (k >= j)
                continue;

            uint64_t intersect = inner_prod(cols[j], cols[k]);
            similarities[j][k]
                = intersect / (num_set_bits[j] + num_set_bits[k] - intersect);
            ++progress_bar;
        }
    }

    return similarities;
}

// For each vector j return similarities with vectors 0, ..., j-1
std::vector<std::vector<uint64_t>>
estimate_similarities(const VectorPtrs &vectors,
                      size_t num_threads,
                      uint64_t num_rows_subsampled) {
    if (!vectors.size())
        return {};

    uint64_t num_sampled_rows = std::min(num_rows_subsampled, vectors[0]->size());

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

template <typename P>
inline bool first_closest(P first_pair, P second_pair) {
    auto first_dist = dist(std::get<0>(first_pair), std::get<1>(first_pair));
    auto second_dist = dist(std::get<0>(second_pair), std::get<1>(second_pair));
    return first_dist < second_dist
        || (first_dist == second_dist
                && std::min(std::get<0>(first_pair), std::get<1>(first_pair))
                    < std::min(std::get<0>(second_pair), std::get<1>(second_pair)));
}

// input: columns
// output: partition, for instance -- a set of column pairs
Partition greedy_matching(const VectorPtrs &columns,
                          size_t num_threads,
                          uint64_t num_rows_subsampled) {
    if (!columns.size())
        return {};

    if (columns.size() > std::numeric_limits<uint32_t>::max()) {
        std::cerr << "ERROR: too many columns" << std::endl;
        exit(1);
    }

    auto similarities = estimate_similarities(columns, num_threads, num_rows_subsampled);

    ProgressBar progress_bar(columns.size() * (columns.size() - 1), "Clustering",
                             std::cerr, !utils::get_verbose());

    std::vector<std::tuple<uint32_t, uint32_t, uint64_t>> candidates;
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
                        && first_closest(first_pair, second_pair));
        },
        num_threads
    );

    Partition partition;
    partition.reserve((columns.size() + 1) / 2);

    std::vector<uint_fast8_t> matched(columns.size(), false);

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
