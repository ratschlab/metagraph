#include "clustering.hpp"

#include <ips4o.hpp>
#include <progress_bar.hpp>

#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "common/vectors/vector_algorithm.hpp"


namespace mtg {
namespace annot {
namespace binmat {

using mtg::common::logger;

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
                             std::cerr, !common::get_verbose());

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < columns.size(); ++i) {
        const bit_vector &col = *columns[i];
        sdsl::bit_vector &subvector = submatrix[i];

        assert(row_indexes.size() <= col.size());

        subvector = sdsl::bit_vector(row_indexes.size(), false);

        for (size_t j = 0; j < row_indexes.size(); ++j) {
            if (col[row_indexes[j]])
                subvector[j] = true;
        }

        ++progress_bar;
    }

    return submatrix;
}

// returns shrunk columns
std::vector<uint64_t>
sample_row_indexes(uint64_t num_rows, uint64_t num_samples, int seed) {
    std::mt19937 gen(seed);

    num_samples = std::min(num_samples, num_rows);

    auto indexes = utils::sample_indexes(num_rows, num_samples, gen);
    // sort indexes
    std::sort(indexes.begin(), indexes.end());
    // check if indexes are sampled without replacement
    assert(std::unique(indexes.begin(), indexes.end()) == indexes.end());

    return indexes;
}

// returns shrinked columns
std::vector<sdsl::bit_vector>
random_submatrix(const VectorPtrs &columns,
                 uint64_t num_rows_sampled, size_t num_threads, int seed) {
    if (!columns.size())
        return {};

    auto indexes = sample_row_indexes(columns[0]->size(), num_rows_sampled, seed);

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

std::vector<std::tuple<uint32_t, uint32_t, float>>
correlation_similarity(const std::vector<sdsl::bit_vector> &cols,
                       size_t num_threads) {
    if (cols.size() > std::numeric_limits<uint32_t>::max()) {
        std::cerr << "ERROR: too many columns" << std::endl;
        exit(1);
    }

    if (!cols.size())
        return {};

    std::vector<std::tuple<uint32_t, uint32_t, float>>
            similarities(cols.size() * (cols.size() - 1) / 2);

    ProgressBar progress_bar(similarities.size(), "Correlations",
                             std::cerr, !common::get_verbose());

    #pragma omp parallel for num_threads(num_threads) collapse(2) schedule(static, 5)
    for (uint64_t j = 1; j < cols.size(); ++j) {
        for (uint64_t i = 0; i < cols.size(); ++i) {
            if (i < j) {
                float sim = inner_prod(cols[i], cols[j]);
                similarities[(j - 1) * j / 2 + i] = std::tie(i, j, sim);
                ++progress_bar;
            }
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
                             std::cerr, !common::get_verbose());

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

template <typename T>
inline T dist(T first, T second) {
    return first > second
                ? first - second
                : second - first;
}

template <typename P>
inline bool first_closest(const P &first, const P &second) {
    auto first_dist = dist(std::get<0>(first), std::get<1>(first));
    auto second_dist = dist(std::get<0>(second), std::get<1>(second));
    return first_dist < second_dist
        || (first_dist == second_dist
                && std::min(std::get<0>(first), std::get<1>(first))
                    < std::min(std::get<0>(second), std::get<1>(second)));
}

// input: columns
// output: partition, for instance -- a set of column pairs
Partition greedy_matching(const std::vector<sdsl::bit_vector> &columns,
                          size_t num_threads) {
    if (!columns.size())
        return {};

    if (columns.size() > std::numeric_limits<uint32_t>::max()) {
        std::cerr << "ERROR: too many columns" << std::endl;
        exit(1);
    }

    auto similarities = correlation_similarity(columns, num_threads);

    ProgressBar progress_bar(similarities.size(), "Matching",
                             std::cerr, !common::get_verbose());

    // pick either a pair of the most similar columns,
    // or pair closest in the initial arrangement
    ips4o::parallel::sort(similarities.begin(), similarities.end(),
        [](const auto &first, const auto &second) {
              return std::get<2>(first) > std::get<2>(second)
                || (std::get<2>(first) == std::get<2>(second)
                        && first_closest(first, second));
        },
        num_threads
    );

    Partition partition;
    partition.reserve((columns.size() + 1) / 2);

    std::vector<uint_fast8_t> matched(columns.size(), false);

    for (const auto &[i, j, sim] : similarities) {
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

LinkageMatrix
agglomerative_greedy_linkage(std::vector<sdsl::bit_vector>&& columns,
                             size_t num_threads) {
    if (columns.empty())
        return LinkageMatrix(0, 4);

    LinkageMatrix linkage_matrix(columns.size() - 1, 4);
    size_t i = 0;

    uint64_t num_clusters = columns.size();
    std::vector<uint64_t> column_ids
            = utils::arange<uint64_t>(0, columns.size());

    for (size_t level = 1; columns.size() > 1; ++level) {
        logger->trace("Clustering: level {}", level);

        Partition groups = greedy_matching(columns, num_threads);

        assert(groups.size() > 0);
        assert(groups.size() < columns.size());

        std::vector<sdsl::bit_vector> cluster_centers(groups.size());
        std::vector<uint64_t> cluster_ids(groups.size());

        ProgressBar progress_bar(groups.size(), "Merging clusters",
                                 std::cerr, !common::get_verbose());

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t g = 0; g < groups.size(); ++g) {
            // merge into new clusters
            cluster_centers[g] = sdsl::bit_vector(columns[0].size(), 0);
            for (size_t i : groups[g]) {
                cluster_centers[g] |= columns[i];
            }

            uint64_t num_set_bits = sdsl::util::cnt_one_bits(cluster_centers[g]);

            #pragma omp critical
            {
                if (groups[g].size() > 1) {
                    assert(groups[g].size() == 2);
                    cluster_ids[g] = num_clusters;
                    linkage_matrix(i, 0) = column_ids[groups[g][0]];
                    linkage_matrix(i, 1) = column_ids[groups[g][1]];
                    linkage_matrix(i, 2) = num_set_bits;
                    linkage_matrix(i, 3) = cluster_ids[g];
                    num_clusters++;
                    i++;
                } else {
                    assert(groups[g].size() == 1);
                    cluster_ids[g] = column_ids[groups[g][0]];
                }
            }

            ++progress_bar;
        }

        columns.swap(cluster_centers);
        column_ids.swap(cluster_ids);
    }

    assert(i == static_cast<size_t>(linkage_matrix.rows()));

    return linkage_matrix;
}

LinkageMatrix agglomerative_linkage_trivial(size_t num_columns) {
    if (!num_columns)
        return LinkageMatrix(0, 4);

    LinkageMatrix linkage_matrix(num_columns - 1, 4);
    size_t i = 0;

    uint64_t num_clusters = num_columns;
    std::vector<uint64_t> column_ids
            = utils::arange<uint64_t>(0, num_columns);

    for (size_t level = 1; column_ids.size() > 1; ++level) {
        logger->trace("Clustering: level {}", level);

        Partition groups((column_ids.size() - 1) / 2 + 1);
        for (size_t j = 0; j < column_ids.size(); ++j) {
            groups[j / 2].push_back(j);
        }

        assert(groups.size() > 0);
        assert(groups.size() < column_ids.size());

        std::vector<uint64_t> cluster_ids(groups.size());

        for (size_t g = 0; g < groups.size(); ++g) {
            // merge into new clusters
            if (groups[g].size() > 1) {
                assert(groups[g].size() == 2);
                cluster_ids[g] = num_clusters;
                linkage_matrix(i, 0) = column_ids[groups[g][0]];
                linkage_matrix(i, 1) = column_ids[groups[g][1]];
                linkage_matrix(i, 2) = 0;
                linkage_matrix(i, 3) = cluster_ids[g];
                num_clusters++;
                i++;
            } else {
                assert(groups[g].size() == 1);
                cluster_ids[g] = column_ids[groups[g][0]];
            }
        }

        column_ids.swap(cluster_ids);
    }

    assert(i == static_cast<size_t>(linkage_matrix.rows()));

    return linkage_matrix;
}

} // namespace binmat
} // namespace annot
} // namespace mtg
