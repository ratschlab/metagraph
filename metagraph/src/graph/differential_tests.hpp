#include <cmath>
#include <cstdint>
#include <vector>
#include <tuple>

namespace mtg {

class DifferentialTest {
  private:
    double family_wise_error_rate;
    size_t total_hypotheses;
    size_t in_total_kmers;
    size_t out_total_kmers;
    size_t in_unique_kmers;
    size_t out_unique_kmers;
    const double pi = std::atan(1) * 4;
    double gamma = std::sqrt(pi);
    double likelihood_ratio_threshold = DifferentialTest::lrt_threshold();

  public:
    DifferentialTest(double family_wise_error_rate,
                     size_t total_hypotheses,
                     size_t preload_table_size,
                     size_t in_total_kmers,
                     size_t out_total_kmers);

    double poisson_prob(int k, double lambda);

    std::tuple<std::vector<int>, std::vector<double>> sorted_pvalues();

    bool bonferroni_correction(double &pvalue);

    double lrt_threshold();

    std::tuple<bool, double> likelihood_ratio_test(double in_sum, double out_sum);
};
} // namespace mtg