#include <cmath>
#include <cstdint>
#include <vector>
#include <tuple>

namespace mtg {
//class LogFactorialTable {
//  public:
//    LogFactorialTable(size_t size);
//
//    double operator[](size_t k) {
//        if (k < m_size) return m_table[k];
//        return approximate_log_factorial(k);
//    }
//
//  private:
//    double log2e_v = std::log(std::exp(1.0));
//    double log_factorial(size_t k);
//    double approximate_log_factorial(size_t k);
//
//  private:
//    std::vector<double> m_table;
//    size_t m_size;
//};

class DifferentialTest {
  private:
    // size_t preload_table_size;
    //LogFactorialTable lf_table { preload_table_size };  // TODO : turn of, this is temporary to do distribution tests.
    double family_wise_error_rate;
    size_t total_hypotheses;
    size_t out_total_kmers;
    size_t in_total_kmers;
    double gamma = std::sqrt(/* pi */ std::atan(1) * 4);
    double likelihood_ratio_threshold = DifferentialTest::lrt_threshold();

  public:
    DifferentialTest(double family_wise_error_rate,
                     size_t total_hypotheses,
                     size_t preload_table_size,
                     size_t in_total_kmers,
                     size_t out_total_kmers);

    double poisson_prob(int k, double lambda);

    bool bonferroni_correction(double &pvalue);

    double lrt_threshold();

    std::tuple<bool, double> likelihood_ratio_test(double in_sum, double out_sum);
};
} // namespace mtg