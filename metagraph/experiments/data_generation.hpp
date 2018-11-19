#ifndef __DATA_GENERATION_HPP__
#define __DATA_GENERATION_HPP__

#include "bit_vector.hpp"


class DataGenerator {
  public:
    DataGenerator() : gen(rd()) { gen.seed(0); }

    void set_seed(int seed) { gen.seed(seed); }

    std::unique_ptr<bit_vector>
    generate_random_column(uint64_t n, double d);

    std::vector<uint64_t>
    generate_random_ints(uint64_t n, uint64_t begin, uint64_t end);

    std::unique_ptr<bit_vector>
    generate_random_column_fixed_size(uint64_t n, uint64_t n_set_bits);

    std::vector<std::unique_ptr<bit_vector>>
    generate_random_columns(uint64_t n, uint64_t m, double d);

    std::vector<std::unique_ptr<bit_vector>>
    generate_random_columns(uint64_t n, uint64_t m,
                            const std::vector<double> &column_densities);

    // Duplicated rows
    std::vector<std::unique_ptr<bit_vector>>
    generate_random_rows(uint64_t num_generating_rows,
                         uint64_t m,
                         double d,
                         const std::vector<uint32_t> &row_frequencies);

    std::vector<std::unique_ptr<bit_vector>>
    generate_random_rows(uint64_t num_generating_rows,
                         uint64_t m,
                         const std::vector<double> &column_densities,
                         const std::vector<uint32_t> &row_frequencies);

    // Duplicated columns
    std::vector<std::unique_ptr<bit_vector>>
    generate_random_columns(uint64_t n,
                            uint64_t m_generating_columns,
                            double d,
                            const std::vector<uint32_t> &column_frequencies);

    std::vector<std::unique_ptr<bit_vector>>
    generate_random_columns(uint64_t n,
                            uint64_t m_generating_columns,
                            const std::vector<double> &column_densities,
                            const std::vector<uint32_t> &column_frequencies);


    std::vector<uint32_t>
    static generate_row_counts_index_inverse(uint64_t n_distinct,
                                             double starting_count,
                                             uint64_t n);

    std::vector<uint32_t>
    static generate_row_counts_logistic(uint64_t n_distinct,
                                        double supremum,
                                        double decay,
                                        double initial = 0.975);

    template <typename Iterator>
    void shuffle(Iterator begin, Iterator end) {
        std::shuffle(begin, end, gen);
    }

  private:
    std::vector<std::unique_ptr<bit_vector>>
    replicate_shuffle(const std::vector<std::unique_ptr<bit_vector>> &vectors,
                      const std::vector<uint32_t> &frequencies);

    std::random_device rd;
    std::mt19937 gen;
};

#endif // __DATA_GENERATION_HPP__
