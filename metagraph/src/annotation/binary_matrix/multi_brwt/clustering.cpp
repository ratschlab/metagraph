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

double inner_prod(const sdsl::bit_vector &first, const sdsl::bit_vector &second) {
    assert(first.size() == second.size());
    return static_cast<double>(::inner_prod(first, second)) / first.size();
}

double inner_prod(const SparseColumn &first, const SparseColumn &second) {
    const auto &[size_1, col_1] = first;
    const auto &[size_2, col_2] = second;
    auto it_1 = col_1.begin();
    auto it_2 = col_2.begin();

    if (!size_1 || !size_2)
        throw std::runtime_error("Vector size must be non-zero");

    uint64_t prod = 0;

    while (it_1 != col_1.end() && it_2 != col_2.end()) {
        assert(*it_1 < size_1 && *it_2 < size_2);
        if (*it_1 < *it_2) {
            ++it_1;
        } else if (*it_1 > *it_2) {
            ++it_2;
        } else {
            prod++;
            ++it_1;
            ++it_2;
        }
    }

    return (double)prod / std::min(size_1, size_2);
}

template <class T>
std::vector<std::tuple<uint32_t, uint32_t, float>>
correlation_similarity(const std::vector<T> &cols, size_t num_threads) {
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

// input: columns, where each column `T` is either `sdsl::bit_vector` or
// `std::pair<uint32_t, std::vector<uint32>>` storing column size and positions
// of its set bits
// output: partition -- a set of column pairs greedily matched
template <class T>
Partition greedy_matching(const std::vector<T> &columns, size_t num_threads) {
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

void merge(const sdsl::bit_vector &first, sdsl::bit_vector *second) {
    assert(second);
    assert(first.size() == second->size());
    *second |= first;
}

void merge(const SparseColumn &first, SparseColumn *second) {
    assert(second);

    const auto &[size_first, col_first] = first;
    const auto &[size_second, col_second] = *second;

    SparseColumn merged;
    auto &[size_merged, col_merged] = merged;

    size_merged = std::min(size_first, size_second);

    col_merged.reserve(col_first.size() + col_second.size());

    auto it_first = col_first.begin();
    auto it_second = col_second.begin();

    while (it_first != col_first.end() && it_second != col_second.end()) {
        if (*it_first < *it_second) {
            col_merged.push_back(*it_first);
            ++it_first;
        } else if (*it_first > *it_second) {
            col_merged.push_back(*it_second);
            ++it_second;
        } else {
            col_merged.push_back(*it_first);
            ++it_first;
            ++it_second;
        }
    }
    while (it_first != col_first.end() && *it_first < size_merged) {
        col_merged.push_back(*it_first);
        ++it_first;
    }
    while (it_second != col_second.end() && *it_second < size_merged) {
        col_merged.push_back(*it_second);
        ++it_second;
    }

    second->swap(merged);
}

uint64_t count_set_bits(const sdsl::bit_vector &v) {
    return sdsl::util::cnt_one_bits(v);
}

uint64_t count_set_bits(const SparseColumn &v) {
    return v.second.size();
}

template <class T>
LinkageMatrix agglomerative_greedy_linkage(std::vector<T>&& columns, size_t num_threads) {
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

        std::vector<T> cluster_centers(groups.size());
        std::vector<uint64_t> cluster_ids(groups.size());

        ProgressBar progress_bar(groups.size(), "Merging clusters",
                                 std::cerr, !common::get_verbose());

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t g = 0; g < groups.size(); ++g) {
            // merge into new clusters
            cluster_centers[g] = std::move(columns[groups[g].at(0)]);
            for (size_t i = 1; i < groups[g].size(); ++i) {
                merge(columns[groups[g][i]], &cluster_centers[g]);
                columns[groups[g][i]] = T();
            }

            uint64_t num_set_bits = count_set_bits(cluster_centers[g]);

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

template
LinkageMatrix agglomerative_greedy_linkage(std::vector<sdsl::bit_vector>&&, size_t);

template
LinkageMatrix agglomerative_greedy_linkage(std::vector<SparseColumn>&&, size_t);


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
