#include <cmath>
#include <cstdint>
#include <vector>
#include <tuple>
#include <numeric>
#include "look_up_table_t_dist.hpp"

namespace mtg {

class DifferentialTest {
  private:
    double family_wise_error_rate;
    size_t total_hypotheses;
    size_t in_total_kmers;
    size_t out_total_kmers;
    size_t in_unique_kmers;
    size_t out_unique_kmers;
    TDistributionTable t_table;
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
    std::tuple<bool, double> brunner_munzel_test(std::vector<int> in_counts, std::vector<int> out_counts);

    int benjamini_yekutieli(std::vector<std::pair<double, int>> likelihood_ratios, double lrt_threshold);

    bool bonferroni_correction(double &pvalue);

    double lrt_threshold();

    double get_t_test_alpha(int df, double alpha);

    std::vector<double> get_midranks(std::vector<int> in_counts, int size_in_counts);

    std::tuple<bool, double> likelihood_ratio_test(double in_sum, double out_sum);
};
} // namespace mtg