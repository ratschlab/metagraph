#include <tsl/hopscotch_set.h>
#include <algorithm>

#include "common/data_generation.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "common/vectors/transpose.hpp"

std::vector<uint64_t>
DataGenerator::generate_random_ints(uint64_t n, uint64_t begin, uint64_t end) {
    std::vector<uint64_t> ints;
    ints.reserve(n);

    std::uniform_int_distribution<> dis(begin, end - 1);
    while (n--) {
        ints.emplace_back(dis(gen));
    }

    return ints;
}

sdsl::bit_vector
generate_random_column_uncompressed(std::mt19937 &gen, uint64_t n, double d) {
    assert(d >= 0. && d <= 1.);

    bool invert = false;
    if (d > 0.5) {
        d = 1 - d;
        invert = true;
    }

    sdsl::bit_vector builder(n, false);
    std::bernoulli_distribution dis(d);
    for (size_t i = 0; i < builder.size(); ++i) {
        if (dis(gen))
            builder[i] = true;
    }

    if (invert)
        builder.flip();

    return builder;
}

void
generate_random_noise(sdsl::bit_vector *vector, std::mt19937 &gen, double d) {
    assert(vector);

    std::bernoulli_distribution dis(d);
    for (size_t i = 0; i < vector->size(); ++i) {
        if (dis(gen))
            (*vector)[i] = true;
    }
}

sdsl::bit_vector DataGenerator::generate_random_column(uint64_t n, double d) {
    return generate_random_column_uncompressed(gen, n, d);
}

std::unique_ptr<bit_vector>
DataGenerator::generate_random_column_fixed_size(uint64_t n, uint64_t n_set_bits) {
    std::uniform_int_distribution<> dis(0, n - 1);

    sdsl::bit_vector bv(n, false);

    while (n_set_bits) {
        uint64_t pos = dis(gen);
        if (!bv[pos]) {
            bv[pos] = true;
            --n_set_bits;
        }
    }

    return std::make_unique<bit_vector_stat>(std::move(bv));
}

std::vector<std::unique_ptr<bit_vector>>
DataGenerator
::replicate_shuffle(const std::vector<std::unique_ptr<bit_vector>> &vectors,
                    const std::vector<uint32_t> &frequencies) {
    assert(vectors.size());
    assert(vectors.size() == frequencies.size());

    std::vector<std::unique_ptr<bit_vector>> vectors_dup;
    vectors_dup.reserve(std::accumulate(frequencies.begin(),
                                        frequencies.end(), 0));

    for (uint64_t i = 0; i < vectors.size(); ++i) {
        auto vect = vectors[i]->copy_to<bit_vector_stat>();

        std::generate_n(std::back_inserter(vectors_dup), frequencies[i],
                        [&vect]() { return std::make_unique<bit_vector_stat>(vect); });
    }

    shuffle(vectors_dup.begin(), vectors_dup.end());

    return vectors_dup;
}

std::vector<std::unique_ptr<bit_vector>>
DataGenerator
::generate_random_columns(uint64_t n, uint64_t m,
                          const std::vector<double> &column_densities) {
    assert(column_densities.size() == m);

    std::vector<std::unique_ptr<bit_vector>> columns;
    columns.reserve(m);
    for (double density : column_densities) {
        columns.emplace_back(new bit_vector_stat(generate_random_column(n, density)));
    }
    return columns;
}

std::vector<std::unique_ptr<bit_vector>>
DataGenerator
::generate_random_columns(uint64_t n, uint64_t m, double d) {
    std::vector<std::unique_ptr<bit_vector>> columns(m);
    for (size_t i = 0; i < m; ++i) {
        columns[i] = std::make_unique<bit_vector_stat>(generate_random_column(n, d));
    }
    return columns;
}

std::vector<std::unique_ptr<bit_vector>>
DataGenerator
::generate_random_rows(uint64_t n_distinct,
                       uint64_t m,
                       double d,
                       const std::vector<uint32_t> &row_frequencies) {
    return generate_random_rows(n_distinct, m,
                                std::vector<double>(m, d),
                                row_frequencies);
}

std::vector<std::unique_ptr<bit_vector>>
DataGenerator
::generate_random_rows(uint64_t n_distinct,
                       uint64_t m,
                       const std::vector<double> &column_densities,
                       const std::vector<uint32_t> &row_frequencies) {
    assert(column_densities.size() == m);
    assert(row_frequencies.size() == n_distinct);

    auto replicated_rows = replicate_shuffle(
        utils::transpose<bit_vector_stat>(
            generate_random_columns(n_distinct, m, column_densities)
        ),
        row_frequencies
    );
    return utils::transpose<bit_vector_stat>(replicated_rows);
}

std::vector<std::unique_ptr<bit_vector>>
DataGenerator
::generate_random_columns(uint64_t n, uint64_t m_generating_columns,
                          double d,
                          const std::vector<uint32_t> &column_frequencies) {
    return generate_random_columns(n, m_generating_columns,
                                   std::vector<double>(m_generating_columns, d),
                                   column_frequencies);
}

std::vector<std::unique_ptr<bit_vector>>
DataGenerator
::generate_random_columns(uint64_t n, uint64_t m_generating_columns,
                          const std::vector<double> &column_densities,
                          const std::vector<uint32_t> &column_frequencies) {
    assert(column_densities.size() == m_generating_columns);
    assert(column_frequencies.size() == m_generating_columns);

    auto replicated_columns = replicate_shuffle(
        generate_random_columns(n, m_generating_columns, column_densities),
        column_frequencies
    );

    return replicated_columns;
}

std::vector<uint32_t>
DataGenerator::generate_row_counts_index_inverse(uint64_t size,
                                                 double starting_count,
                                                 uint64_t n) {
    assert(size);

    double factor = 0;
    for (uint64_t i = 1; i <= size; ++i) {
        factor += 1.0 / i;
    }
    factor = static_cast<double>(n) / factor / starting_count;

    std::vector<uint32_t> frequencies;
    for (uint64_t i = 1; i <= size; ++i) {
        frequencies.emplace_back(starting_count * factor / i);
        if (!frequencies.back()) {
            frequencies.back() = 1;
            frequencies.insert(frequencies.end(), size - frequencies.size(), 1);
            break;
        }
    }

    uint64_t final_num_rows = std::accumulate(frequencies.begin(),
                                              frequencies.end(),
                                              0u);

    if (final_num_rows < n) {
        std::cerr << "Expected: " << n << ", computed: " << final_num_rows << std::endl;
        throw std::runtime_error("ERROR: initial count too low");
    }

    assert(frequencies.size() == size);
    return frequencies;
}

std::vector<uint32_t>
DataGenerator
::generate_row_counts_logistic(uint64_t size,
                               double supremum,
                               double decay,
                               double initial) {
    assert(size);
    double logsup = std::log10(supremum);
    std::vector<double> logfrequencies = { logsup * initial };
    while (logfrequencies.size() < size) {
        double cur = logfrequencies.back();
        assert(cur <= logsup);
        assert(cur * decay * (logsup - cur)/logsup >= 0);
        logfrequencies.emplace_back(
            cur - cur * decay * (logsup - cur)/logsup
        );
        assert(logfrequencies.back() < logsup);
        assert(logfrequencies.back() >= 0);
    }

    std::vector<uint32_t> frequencies;
    frequencies.reserve(size);
    std::transform(logfrequencies.begin(), logfrequencies.end(),
                   std::back_inserter(frequencies),
                   [](double logfreq) -> uint64_t {
                       return std::pow(10.0, logfreq);
                   });

    assert(std::all_of(frequencies.begin(), frequencies.end(),
        [&supremum](uint64_t freq) -> bool {
            return freq >= 1 && freq < supremum;
        }));
    return frequencies;
}
