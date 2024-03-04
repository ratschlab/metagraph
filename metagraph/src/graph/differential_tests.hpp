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
    double alpha_precalc = -1;
    int df_precalc = -1;
    const double pi = std::atan(1) * 4;
    double gamma = std::sqrt(pi);
    double likelihood_ratio_threshold = DifferentialTest::lrt_threshold();

  public:
    int64_t var_0 = 0;
    DifferentialTest(double family_wise_error_rate,
                     size_t total_hypotheses,
                     size_t preload_table_size,
                     size_t in_total_kmers,
                     size_t out_total_kmers);

    double poisson_prob(int k, double lambda);
    std::tuple<bool, double> brunner_munzel_test(std::vector<double> in_counts, std::vector<double> out_counts);

    int benjamini_yekutieli(std::vector<std::pair<double, int>> likelihood_ratios, double lrt_threshold);

    bool bonferroni_correction(double &pvalue);

    double lrt_threshold();

    double get_t_test_alpha(int df, double alpha);

    std::tuple<std::vector<double>, std::vector<double>> get_overall_midrank(std::vector<double> in_counts, std::vector<double> out_counts, int N);

    double get_var_midranks(std::vector<double> mid_ranks_within, std::vector<double> mid_ranks_overall, int m);

    std::vector<double> get_midranks(std::vector<double> in_counts, int size_in_counts);

    std::tuple<bool, double> likelihood_ratio_test(double in_sum, double out_sum);

    int get_df_approx(std::vector<double> in_counts, std::vector<double> out_counts);

    int get_df_conservative(std::vector<double> in_counts, std::vector<double> out_counts);

    double get_var(std::vector<double> counts, int n);
};
} // namespace mtg