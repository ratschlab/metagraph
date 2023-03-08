#include <cmath>
#include <cstdint>
#include <vector>
#include <tuple>

#include "differential_tests_notes.hpp" // TODO rename file when running cmake

class LogFactorialTable
{
public:
    LogFactorialTable(size_t size);

    double operator[](size_t i)
    {
        if (i < m_size) return m_table[i];
        return log_factorial(i);
    }

private:
    double log_factorial(size_t k);

private:
    std::vector<double> m_table;
    size_t m_size;
};

class DifferentialTest
{
private:
    size_t preload_table_size;
    Chi2PLookup Chi2PLookupTable; // generated with https://github.com/MoseleyBioinformaticsLab/chi2plookup
    LogFactorialTable lf_table{preload_table_size};
    double family_wise_error_rate;
    size_t total_hypotheses;
    size_t out_total_kmers;
    size_t in_total_kmers;
    double gamma = std::sqrt(/* pi */ std::atan(1)*4);

public:
    DifferentialTest(double family_wise_error_rate, size_t total_hypotheses, size_t preload_table_size,
                     size_t in_total_kmers, size_t out_total_kmers);

    double chi_sqrd_1_cdf(double x);

    static double lower_incomplete_gamma_function(double x);

    double poisson_prob(int k, double lambda);

    bool bonferroni_correction(double& pvalue);

    std::tuple<double, bool> likelihood_ratio_test(double in_sum, double out_sum);
};
