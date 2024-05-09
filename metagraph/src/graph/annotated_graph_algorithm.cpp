#include "annotated_graph_algorithm.hpp"

#include <typeinfo>

#include <sdust.h>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/tools/roots.hpp>

#include "common/utils/string_utils.hpp"
#include "common/logger.hpp"
#include "common/vectors/bitmap.hpp"
#include "common/vector_map.hpp"
#include "common/vectors/transpose.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/graph_cleaning.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "differential_tests.hpp"

namespace mtg {
namespace graph {

using mtg::common::logger;

typedef AnnotatedDBG::node_index node_index;
typedef AnnotatedDBG::Annotator Annotator;
typedef AnnotatedDBG::Annotator::Label Label;
using Column = annot::matrix::BinaryMatrix::Column;
using PairContainer = std::vector<std::pair<uint64_t, uint64_t>>;
using ValuesContainer = sdsl::int_vector_buffer<>;

double CDF_CUTOFF = 0.95;
uint64_t N_BUCKETS_FOR_ESTIMATION = 3;

template <typename ForwardIterator>
inline auto quartiles(ForwardIterator first, ForwardIterator last)
{
    auto m = std::distance(first,last);
    BOOST_MATH_ASSERT_MSG(m >= 3, "At least 3 samples are required to compute the interquartile range.");
    auto k = m/4;
    auto j = m - (4*k);
    // m = 4k+j.
    // If j = 0 or j = 1, then there are an even number of samples below the median, and an even number above the median.
    //    Then we must average adjacent elements to get the quartiles.
    // If j = 2 or j = 3, there are an odd number of samples above and below the median, these elements may be directly extracted to get the quartiles.

    double Q1;
    double Q3;
    if (j==2 || j==3)
    {
        auto q1 = first + k;
        auto q3 = first + 3*k + j - 1;
        std::nth_element(first, q1, last);
        Q1 = *q1;
        std::nth_element(q1, q3, last);
        Q3 = *q3;
    } else {
        // j == 0 or j==1:
        auto q1 = first + k - 1;
        auto q3 = first + 3*k - 1 + j;
        std::nth_element(first, q1, last);
        double a = *q1;
        std::nth_element(q1, q1 + 1, last);
        double b = *(q1 + 1);
        Q1 = (a+b)/2;
        std::nth_element(q1, q3, last);
        a = *q3;
        std::nth_element(q3, q3 + 1, last);
        b = *(q3 + 1);
        Q3 = (a+b)/2;
    }

    // return std::make_tuple(Q1, boost::math::statistics::median(first, last), Q3);
    return std::make_pair(Q1, Q3);
}

inline bool is_low_complexity(std::string_view s, int T = 20, int W = 64) {
    int n;
    std::unique_ptr<uint64_t, decltype(std::free)*> r {
        sdust(0, (const uint8_t*)s.data(), s.size(), T, W, &n),
        std::free
    };
    return n > 0;
}


std::pair<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>>
mask_nodes_by_label_dual(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                         const std::vector<std::string> &files,
                         const tsl::hopscotch_set<Label> &labels_in,
                         const tsl::hopscotch_set<Label> &labels_out,
                         const DifferentialAssemblyConfig &config,
                         size_t num_threads,
                         size_t num_parallel_files) {
    bool is_primary = graph_ptr->get_mode() == DeBruijnGraph::PRIMARY;
    common::logger->trace("Graph mode: {}", is_primary ? "PRIMARY" : "other");
    common::logger->trace("Labels in: {}", fmt::join(labels_in, ","));
    common::logger->trace("Labels out: {}", fmt::join(labels_out, ","));

    std::vector<bool> groups(labels_in.size() + labels_out.size());
    std::vector<std::unique_ptr<bit_vector>> columns_all(groups.size());
    std::vector<std::unique_ptr<ValuesContainer>> column_values_all(groups.size());

    uint8_t max_width = 0;

    sdsl::bit_vector ignored(AnnotatedDBG::graph_to_anno_index(graph_ptr->max_index() + 1), false);

    annot::ColumnCompressed<>::load_columns_and_values(
        files,
        [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> &&column, ValuesContainer&& column_values) {
            groups[offset] = labels_out.count(label);
            max_width = std::max(max_width, column_values.width());
            std::swap(columns_all[offset], column);
            column_values_all[offset] = std::make_unique<ValuesContainer>(std::move(column_values));
        },
        num_parallel_files
    );


    std::vector<uint64_t> min_counts(groups.size(), 0);
    std::vector<uint64_t> check_cutoff(groups.size(), std::numeric_limits<uint64_t>::max());

    auto adjust_count = [&](uint64_t raw_count, size_t j, uint64_t row_i) {
        std::ignore = row_i;
        return raw_count >= min_counts[j] && raw_count <= check_cutoff[j]
            ? raw_count
            : 0;
    };

    if (config.clean && !config.min_count) {
        common::logger->trace("Cleaning count columns");
        for (size_t j = 0; j < groups.size(); ++j) {
            const auto &column_values = *column_values_all[j];

            // set cutoff for lower end of distribution
            auto [mean_est, nzeros_est] = estimate_ztp_mean(
                [&](const auto &callback) {
                    for (size_t i = 0; i < column_values.size(); ++i) {
                        callback(column_values[i]);
                    }
                },
                max_width,
                0,
                N_BUCKETS_FOR_ESTIMATION
            );

            min_counts[j] = boost::math::quantile(
                boost::math::poisson(mean_est), exp(-mean_est) - CDF_CUTOFF * expm1(-mean_est)
            );

            // quartile calculation modifies the underlying vector, so we have to make a copy
            std::vector<uint64_t> vals;
            vals.reserve(column_values.size());
            for (size_t i = 0; i < column_values.size(); ++i) {
                if (column_values[i] >= min_counts[j])
                    vals.emplace_back(column_values[i]);
            }

            // set cutoff for upper end of distribution
            auto [q1, q3] = quartiles(vals.begin(), vals.end());
            double outlier_cutoff = q3 + 1.5 * (q3 - q1);
            double repeats_mu = 0;
            double repeats_cutoff = 0;
            size_t num_large = 0;
            for (uint64_t c : vals) {
                if (c >= outlier_cutoff) {
                    repeats_mu += c;
                    ++num_large;
                }
            }
            if (num_large) {
                repeats_mu /= num_large;

                repeats_cutoff = boost::math::quantile(boost::math::complement(
                    boost::math::poisson(repeats_mu), CDF_CUTOFF
                ));

                if (repeats_cutoff > outlier_cutoff)
                    check_cutoff[j] = repeats_cutoff;
            }
        }
    }

    common::logger->trace("Marking discarded k-mers");
    if (groups.size()) {
        utils::call_rows<std::unique_ptr<bit_vector>,
                         std::unique_ptr<ValuesContainer>,
                         PairContainer>(columns_all, column_values_all,
                                        [&](uint64_t row_i, const auto &row) {
            if (row.empty()) {
                ignored[row_i] = true;
            } else if (config.clean) {
                uint64_t sum = 0;
                for (const auto &[j, raw_c] : row) {
                    if (uint64_t c = adjust_count(raw_c, j, row_i))
                        sum += c;
                }

                ignored[row_i] = sum < config.min_count;
            }
        });
    }

    common::logger->trace("Computing aggregate statistics");
    double in_kmers = 0;
    double out_kmers = 0;
    std::vector<double> medians;
    std::vector<uint64_t> max_obs_vals;
    std::vector<size_t> n_kmers;
    std::vector<uint64_t> sums;
    std::vector<uint64_t> sums_of_squares;
    for (size_t j = 0; j < groups.size(); ++j) {
        const auto &column = *columns_all[j];
        const auto &column_values = *column_values_all[j];
        uint64_t sum = 0;
        uint64_t sum_of_squares = 0;
        size_t i = 0;
        size_t n = 0;

        std::vector<uint64_t> vals;
        vals.reserve(column_values.size());
        column.call_ones([&](uint64_t row_i) {
            if (uint64_t c = adjust_count(column_values[i], j, row_i)) {
                vals.emplace_back(c);
                sum += c;
                sum_of_squares += c * c;
                ++n;
            }
            ++i;
        });

        max_obs_vals.emplace_back(*std::max_element(vals.begin(), vals.end()));
        medians.push_back(boost::math::statistics::median(vals));

        sums.push_back(sum);
        if (groups[j]) {
            out_kmers += sum;
        } else {
            in_kmers += sum;
        }
        sums_of_squares.push_back(sum_of_squares);
        n_kmers.push_back(n);

        common::logger->trace("{}: sum: {}\tmedian: {}\tmax_obs: {}\tmin_cutoff: {}\tmax_cutof: {}",
                              j, sum, medians[j], max_obs_vals[j], min_counts[j], check_cutoff[j]);
    }

    double total_kmers = in_kmers + out_kmers;
    double nelem = ignored.size() - sdsl::util::cnt_one_bits(ignored);

    common::logger->trace("Number of kept unique k-mers: {}\tNumber of kept k-mers: {}",
                          nelem, total_kmers);

    common::logger->trace("Running differential tests");
    sdsl::bit_vector indicator_in(graph_ptr->max_index() + 1, false);
    sdsl::bit_vector indicator_out(graph_ptr->max_index() + 1, false);

    std::vector<double> pvals;
    pvals.reserve(graph_ptr->max_index() + 1);
    pvals.emplace_back(1.1);

    VectorMap<size_t, std::vector<uint64_t>> m;
    using RowStats = std::tuple<double, double, double>;
    std::function<RowStats(const PairContainer&)> compute_pval;

    if (config.test_type == "likelihoodratio" || config.test_type == "likelihoodratio_unitig"
            || config.test_type == "cmh"
            || config.test_type == "fisher") {
        boost::math::chi_squared dist(1);
        compute_pval = [&,dist](const auto &row) {
            if (row.empty())
                return std::make_tuple(1.1, 0.0, 0.0);

            double in_sum = 0;
            double out_sum = 0;

            for (const auto &[j, c] : row) {
                if (groups[j]) {
                    out_sum += c;
                } else {
                    in_sum += c;
                }
            }

            double chi_stat = 0;
            double pval = -1;
            if (config.test_type == "cmh") {
                if (in_sum + out_sum == total_kmers || total_kmers <= 1 || in_sum + out_sum == 0) {
                    common::logger->error("Invalid counts: in_sum: {}\tout_sum: {}\ttotal_kmers: {}",
                                          in_sum, out_sum, total_kmers);
                }

                double xi = in_sum * out_kmers - out_sum * in_kmers;
                xi *= xi * (total_kmers - 1) / in_kmers / out_kmers / (in_sum + out_sum) / (total_kmers - in_sum - out_sum);
                chi_stat = xi;
            } else if (config.test_type == "fisher") {
                double lp_base = lgamma(in_kmers + 1) - lgamma(in_kmers - in_sum + 1) - lgamma(in_sum + 1);
                lp_base += lgamma(out_kmers + 1) - lgamma(out_kmers - out_sum + 1) - lgamma(out_sum + 1);

                double denom_in = lgamma(total_kmers + 1) - lgamma(total_kmers - in_sum - out_sum + 1) - lgamma(in_sum + out_sum + 1);

                pval = exp(lp_base - denom_in);
                chi_stat = boost::math::quantile(boost::math::complement(dist, pval));
            } else {
                double log_likelihood_in = in_sum > 0 ? -in_sum + in_sum * log(in_sum) : 0.0;
                double log_likelihood_out = out_sum > 0 ? -out_sum + out_sum * log(out_sum) : 0.0;

                double theta = (in_sum + out_sum) / total_kmers;
                double mean_denom_in = theta * in_kmers;
                double mean_denom_out = theta * out_kmers;
                double log_likelihood_null = mean_denom_in > 0 && in_sum >= 0 ? -mean_denom_in + in_sum * log(mean_denom_in) : 0.0;
                log_likelihood_null += mean_denom_out > 0 && out_sum >= 0 ? -mean_denom_out + out_sum * log(mean_denom_out) : 0.0;

                chi_stat = (log_likelihood_in + log_likelihood_out - log_likelihood_null) * 2;
            }

            if (chi_stat == 0)
                return std::make_tuple(1.0, in_sum, out_sum);

            if (chi_stat < 0) {
                common::logger->error("Detected likelihood ratio {} < 0", chi_stat);
                throw std::runtime_error("Test failed");
            }

            if (pval == -1)
                pval = boost::math::cdf(boost::math::complement(dist, chi_stat));

            return std::make_tuple(pval, in_sum, out_sum / out_kmers * in_kmers);
        };
    } else if (config.test_type == "nbinom") {
        // use moments to fit the r parameter of a beta negative binomial
        uint64_t sum_of_cubes = 0;
        for (size_t j = 0; j < groups.size(); ++j) {
            const auto &column = *columns_all[j];
            const auto &column_values = *column_values_all[j];
            size_t i = 0;
            column.call_ones([&](uint64_t row_i) {
                if (!ignored[row_i]) {
                    if (int64_t c = adjust_count(column_values[i], j, row_i))
                        sum_of_cubes += c * (c - 1) * (c - 2);
                }
                ++i;
            });
        }
        double mu = static_cast<double>(total_kmers) / groups.size() / nelem;
        double mu2 = static_cast<double>(std::accumulate(sums_of_squares.begin(),
                                                         sums_of_squares.end(), uint64_t(0))) / groups.size() / nelem - mu;
        double mu3 = static_cast<double>(sum_of_cubes) / groups.size() / nelem;
        double a = mu3 * mu + mu * mu * mu2 + 2 * mu2 * mu2;
        double b = 2 * mu3 * mu * mu + mu3 * mu - mu2 * mu3 - mu * mu2 * mu2 + 3 * mu * mu + mu2 - 4 * mu2 * mu2;
        double c = -2 * mu * (mu2 * mu2 - mu * mu2 - mu3 * mu);

        auto [r_low, r_high] = boost::math::tools::quadratic_roots(a, b, c);

        common::logger->trace("mu: {}\tmu2: {}\tmu3: {}\tr_low: {}\tr_high: {}",
                              mu, mu2, mu3, r_low, r_high);

        // check for NaNs
        // https://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
        if (r_low != r_low && r_high != r_high) {
            common::logger->error("{},{},{}", a, b, c);
            throw std::runtime_error("Failed to fit, no real roots");
        }

        double r_low_var_sqerr = 0;
        double r_high_var_sqerr = 0;
        for (size_t j = 0; j < groups.size(); ++j) {
            double mu = static_cast<double>(sums[j]) / nelem;
            double var = static_cast<double>(sums_of_squares[j]) / nelem - mu * mu;
            double p_low = r_low / (r_low + mu);
            double p_high = r_high / (r_high + mu);
            double t_var_low = mu / p_low;
            double t_var_high = mu / p_high;
            r_low_var_sqerr += pow(t_var_low - var, 2.0);
            r_high_var_sqerr += pow(t_var_high - var, 2.0);
            common::logger->trace("label: {}\tmean: {}\tvar: {}\tp_low: {}\tvar_low: {}\tp_high: {}\tvar_high: {}",
                                  j, mu, var, p_low, t_var_low, p_high, t_var_high);
        }

        common::logger->trace("r low: {} (sqerr: {})\t r high: {} (sqerr: {})",
                              r_low, r_low_var_sqerr, r_high, r_high_var_sqerr);

        double r = r_low_var_sqerr < r_high_var_sqerr ? r_low : r_high;
        common::logger->trace("r: {}", r);

        boost::math::chi_squared dist(1);

        compute_pval = [&,r,dist](const auto &row) {
            if (row.empty())
                return std::make_tuple(1.1, 0.0, 0.0);

            double in_sum = 0;
            double out_sum = 0;
            sdsl::bit_vector found(groups.size());
            for (const auto &[j, c] : row) {
                found[j] = true;
                if (groups[j]) {
                    out_sum += c;
                } else {
                    in_sum += c;
                }
            }

            std::vector<size_t> unfound;
            for (size_t i = 0; i < found.size(); ++i) {
                if (!found[i])
                    unfound.emplace_back(i);
            }

            double total_sum = in_sum + out_sum;
            size_t ndigits = 30;

            auto ll_dx_ddx_pick = [&](const auto &picker) {
                return [&](double lambda) {
                    double val = 0;
                    double dval = 0;
                    for (const auto &[j, c] : row) {
                        if (picker(j)) {
                            double prop = static_cast<double>(sums[j]) * lambda;
                            double shift_prop = r + prop;
                            val += (c - prop) / shift_prop;
                            dval -= (r + c) / shift_prop / shift_prop * sums[j];
                        }
                    }

                    for (size_t j : unfound) {
                        if (picker(j)) {
                            double prop = static_cast<double>(sums[j]) * lambda;
                            double shift_prop = r + prop;
                            val -= prop / shift_prop;
                            dval -= r / shift_prop / shift_prop * sums[j];
                        }
                    }

                    return std::make_pair(val, dval);
                };
            };

            double lambda_null = boost::math::tools::newton_raphson_iterate(
                ll_dx_ddx_pick([&](size_t) { return true; }),
                total_sum / total_kmers, -0.01, 1.01, ndigits
            );

            double chi_stat = 0;
            double loglikelihood_null = 0;
            // if (false) {
            //     // score test
            //     auto [val_in, dval_in] = ll_dx_ddx_in(lambda_null);
            //     auto [val_out, dval_out] = ll_dx_ddx_out(lambda_null);

            //     val_in *= r / lambda_null;
            //     val_out *= r / lambda_null;

            //     dval_in = (val_in / lambda_null - dval_in * r / lambda_null) / labels_in.size();
            //     dval_out = (val_out / lambda_null - dval_out * r / lambda_null) / labels_out.size();
            //     chi_stat = val_in * val_in / dval_in + val_out * val_out / dval_out;
            // } else {
                // log likelihood test
                double lambda_in = boost::math::tools::newton_raphson_iterate(
                    ll_dx_ddx_pick([&](size_t j) { return !groups[j]; }),
                    in_sum / in_kmers, -0.01, 1.01, ndigits
                );

                double lambda_out = boost::math::tools::newton_raphson_iterate(
                    ll_dx_ddx_pick([&](size_t j) { return groups[j]; }),
                    out_sum / out_kmers, -0.01, 1.01, ndigits
                );

                double loglikelihood_alt = 0;

                for (const auto &[j, c] : row) {
                    double lambda = groups[j] ? lambda_out : lambda_in;
                    double prop = lambda * sums[j];
                    loglikelihood_alt += log(lambda) * c - log(r + prop) * (r + c);

                    double prop_null = lambda_null * sums[j];
                    loglikelihood_null += log(lambda_null) * c - log(r + prop_null) * (r + c);
                }

                for (size_t j : unfound) {
                    double lambda = groups[j] ? lambda_out : lambda_in;
                    double prop = lambda * sums[j];
                    loglikelihood_alt -= log(r + prop) * r;

                    double prop_null = lambda_null * sums[j];
                    loglikelihood_null -= log(r + prop_null) * r;
                }

                chi_stat = (loglikelihood_alt - loglikelihood_null) * 2;

                if (chi_stat < 0) {
                    common::logger->error("Detected likelihood ratio {} < 0", chi_stat);
                    common::logger->error("Theta null: {}\ttheta in: {}\ttheta out: {}",
                                          lambda_null, lambda_in, lambda_out);
                    throw std::runtime_error("Test failed");
                }
            // }

            double pval = chi_stat > 0 ? boost::math::cdf(boost::math::complement(dist, chi_stat)) : 1.0;

            return std::make_tuple(pval, in_sum, out_sum / out_kmers * in_kmers);
        };
    } else if (config.test_type == "bp") {
        double sum_of_squares = 0;
        std::vector<double> counts;
        counts.reserve(nelem);
        utils::call_rows<std::unique_ptr<bit_vector>,
                         std::unique_ptr<ValuesContainer>,
                         PairContainer>(columns_all, column_values_all,
                                        [&](uint64_t row_i, const auto &row) {
            if (ignored[row_i])
                return;

            double sum_in = 0;
            double sum_out = 0;
            for (const auto &[j, raw_c] : row) {
                if (uint64_t c = adjust_count(raw_c, j, row_i)) {
                    if (groups[j]) {
                        sum_in += c;
                    } else {
                        sum_out += c;
                    }
                }
            }

            if (sum_in + sum_out > 0) {
                counts.emplace_back(sum_in + sum_out);
                sum_of_squares += pow(counts.back(), 2.0);
            }
        });

        auto get_a = [&](double k, double k2) {
            double mu = k / nelem;
            // model with a negative binomial, similar to
            // https://github.com/JBHilton/beta-poisson-epidemics/blob/master/functions.py#L140
            auto get_dll = [&](double r) {
                double lp = log(r) - log(r + mu);
                double dll = nelem * lp - nelem * boost::math::digamma(r);
                for (double c : counts) {
                    dll += boost::math::digamma(r + c);
                }

                return dll;
            };

            double var = k2 / nelem - mu * mu;
            double r_guess = mu * mu / (var - mu);

            long unsigned int max_iter = 100;
            auto [r_min,r_max] = boost::math::tools::bracket_and_solve_root(
                get_dll, r_guess, 2.0, false, [&](double r_min, double r_max) {
                    return abs(r_min - r_max) / r_max < 1e-5;
                },
                max_iter
            );
            double a = (r_min + r_max) / 2.0 * nelem / (nelem - 1);

            double b = (nelem - 1) * a;

            common::logger->trace("mean: {}\tvar: {}\ta: {}\tb: {}\texp: {}\tvar: {}",
                                  mu, k2 / nelem - mu * mu,
                                  a, b,
                                  k * a / (a + b),
                                  mu + mu * mu * (nelem - 1) / (nelem * a + 1));

            return std::make_pair(a, b);
        };

        auto [a, b] = get_a(total_kmers, sum_of_squares);

        if (a <= 1) {
            common::logger->error("a <= 1");
            throw std::runtime_error("Fit failed");
        }

        auto get_theta = [&,a,b](double total, double k_sum) {
            double qb = 2.0 - total - k_sum - a - b;
            double qc = k_sum + a - 1.0;

            auto [theta_low, theta_high] = boost::math::tools::quadratic_roots(total, qb, qc);

            // Note: comparisons of thetas to themselves ensure that they are not NaN
            // https://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
            bool low_valid = theta_low == theta_low && theta_low > 0 && theta_low < 1.0;
            bool high_valid = theta_high == theta_high && theta_low != theta_high && theta_high > 0 && theta_high < 1.0;

            double stat_low = low_valid
                ? qc * log(theta_low) + (b - 1) * log1p(-theta_low) - total * theta_low
                : -std::numeric_limits<double>::infinity();
            double stat_high = high_valid
                ? qc * log(theta_high) + (b - 1) * log1p(-theta_high) - total * theta_high
                : -std::numeric_limits<double>::infinity();

            double stat_max = std::max({ stat_low, stat_high });

            bool stat_is_low = low_valid && stat_max == stat_low;
            bool stat_is_high = high_valid && stat_max == stat_high;

            if (stat_is_low && stat_is_high) {
                common::logger->error("Multiple valid thetas\t{}", stat_max);
                common::logger->error("{}\t{}", theta_low, theta_high);
                common::logger->error("{}\t{}", stat_low, stat_high);
                common::logger->error("{}\t{}", stat_is_low, stat_is_high);
                common::logger->error("Ratio: {} / {}", k_sum, total);
                common::logger->error("{},{}\t{},{},{}", a,b, total, qb, qc);
                throw std::runtime_error("Fail");
            }

            if (stat_is_low || stat_is_high)
                return std::make_tuple(theta_low, theta_high, stat_is_low, stat_is_high);

            common::logger->error("No valid thetas\t{}", stat_max);
            common::logger->error("{}\t{}", theta_low, theta_high);
            common::logger->error("{}\t{}", stat_low, stat_high);
            common::logger->error("{}\t{}", stat_is_low, stat_is_high);
            common::logger->error("Ratio: {} / {}", k_sum, total);
            common::logger->error("{},{}\t{},{},{}", a,b,total, qb, qc);
            throw std::runtime_error("Fail");
        };

        compute_pval = [&,a,b,nelem,get_theta](const auto &row) {
            if (row.empty())
                return std::make_tuple(1.1, 0.0, 0.0);

            double in_sum = 0;
            double out_sum = 0;
            for (const auto &[j, c] : row) {
                if (groups[j]) {
                    out_sum += c;
                } else {
                    in_sum += c;
                }
            }

            double null_sum = in_sum + out_sum;

            auto [theta_null_low, theta_null_high, null_is_low, null_is_high]
                = get_theta(total_kmers, null_sum);
            auto [theta_in_low, theta_in_high, in_is_low, in_is_high]
                = get_theta(in_kmers, in_sum);
            auto [theta_out_low, theta_out_high, out_is_low, out_is_high]
                = get_theta(out_kmers, out_sum);

            double theta_null = null_is_low ? theta_null_low : theta_null_high;
            double theta_in = in_is_low ? theta_in_low : theta_in_high;
            double theta_out = out_is_low ? theta_out_low : theta_out_high;

            auto get_stat_log = [&,a,b](double theta, double k, double total) -> double {
                return (k + a - 1) * log(theta) + (b - 1) * log1p(-theta) - total * theta;
            };

            double map_pval_in = exp(get_stat_log(theta_null, in_sum, in_kmers)
                                    - get_stat_log(theta_in, in_sum, in_kmers));

            double map_pval_out = exp(get_stat_log(theta_null, out_sum, out_kmers)
                                    - get_stat_log(theta_out, out_sum, out_kmers));

            if (map_pval_in > 1.0 || map_pval_in < 0.0
                    || map_pval_out > 1.0 || map_pval_out < 0.0) {
                return std::make_tuple(1.0, in_sum,
                                       out_sum / out_kmers * in_kmers);
                common::logger->error("Failed to compute MAP of thetas");
                common::logger->error("Theta in: {}\tTheta out: {}\tTheta null: {}",
                                      theta_in, theta_out, theta_null);
                common::logger->error("Theta in low: {}\tTheta out low: {}\tTheta null low: {}",
                                      theta_in_low, theta_out_low, theta_null_low);
                common::logger->error("Theta in high: {}\tTheta out high: {}\tTheta null high: {}",
                                      theta_in_high, theta_out_high, theta_null_high);
                common::logger->error("In sum: {}\tOut sum: {}", in_sum, out_sum);
                common::logger->error("In p-val: {}\tOut p-val: {}", map_pval_in, map_pval_out);
                throw std::runtime_error("Fail");
            }

            return std::make_tuple(std::min(std::min(map_pval_in, map_pval_out) * 2.0, 1.0),
                                   in_sum,
                                   out_sum / out_kmers * in_kmers);
        };
    } else if (config.test_type == "nonparametric" || config.test_type == "nonparametric_u" || config.test_type == "quantile") {
        compute_pval = [&](const auto &row) {
            if (row.empty())
                return std::make_tuple(1.1, 0.0, 0.0);

            auto norm = boost::math::normal_distribution();
            auto rankdata = [&](const std::vector<double> &c) {
                auto c_p = c;
                std::sort(c_p.begin(), c_p.end());
                tsl::hopscotch_map<double, std::pair<size_t, size_t>> counts;
                for (size_t i = 0; i < c_p.size(); ++i) {
                    auto &[rank_sum, nranks] = counts[c_p[i]];
                    rank_sum += i + 1;
                    ++nranks;
                }

                std::vector<double> ranks;
                ranks.reserve(c.size());
                for (double v : c) {
                    ranks.emplace_back(static_cast<double>(counts[v].first)/counts[v].second);
                }

                return std::make_pair(std::move(ranks), std::move(counts));
            };

            std::vector<double> in_counts;
            in_counts.reserve(labels_in.size());
            std::vector<double> out_counts;
            out_counts.reserve(labels_out.size());

            for (const auto &[j, c] : row) {
                auto &bucket = groups[j] ? out_counts : in_counts;
                bucket.emplace_back(static_cast<double>(c) / medians[j]);
            }

            in_counts.resize(labels_in.size());
            out_counts.resize(labels_out.size());

            double p = 0.0;
            size_t nx = in_counts.size();
            size_t ny = out_counts.size();

            std::vector<double> c;
            c.reserve(in_counts.size() + out_counts.size());
            std::copy(in_counts.begin(), in_counts.end(), std::back_inserter(c));
            std::copy(out_counts.begin(), out_counts.end(), std::back_inserter(c));
            auto [rankc, countsc] = rankdata(c);
            double rankcx_mean = boost::math::statistics::mean(rankc.begin(), rankc.begin() + nx);
            double rankcy_mean = boost::math::statistics::mean(rankc.begin() + nx, rankc.end());

            if (config.test_type == "nonparametric") {
                if (nx > 1 && ny > 1) {
                    auto [rankx, countsx] = rankdata(in_counts);
                    auto [ranky, countsy] = rankdata(out_counts);
                    double rankx_mean = boost::math::statistics::mean(rankx.begin(), rankx.end());
                    double ranky_mean = boost::math::statistics::mean(ranky.begin(), ranky.end());

                    double Sx = 0;
                    for (size_t i = 0; i < nx; ++i) {
                        Sx += std::pow(rankc[i] - rankx[i] - rankcx_mean + rankx_mean, 2.0);
                    }
                    Sx /= nx - 1;
                    double Sy = 0;
                    for (size_t i = 0; i < ny; ++i) {
                        Sy += std::pow(rankc[nx + i] - ranky[i] - rankcy_mean + ranky_mean, 2.0);
                    }
                    Sy /= ny - 1;

                    double wbfn = (rankcy_mean - rankcx_mean) * nx * ny;
                    wbfn /= std::sqrt(Sx * nx + Sy * ny) * (nx + ny);

                    double df_numer = std::pow(Sx * nx + Sy * ny, 2.0);
                    double df_denom = std::pow(Sx * nx, 2.0) / (nx - 1) + std::pow(Sy * ny, 2.0) / (ny - 1);

                    if (df_denom != 0) {
                        wbfn = abs(wbfn);
                        double df = df_numer / df_denom;
                        boost::math::students_t dist(df);
                        p = 2 * boost::math::cdf(boost::math::complement(dist, wbfn));
                    }
                }
            } else {
                // https://datatab.net/tutorial/mann-whitney-u-test
                double rankcx_sum = std::accumulate(rankc.begin(), rankc.begin() + nx, double(0.0));
                double rankcy_sum = std::accumulate(rankc.begin() + nx, rankc.end(), double(0.0));
                double u1 = nx * ny + static_cast<double>(nx * (nx + 1))/2 - rankcx_sum;
                double u2 = nx * ny + static_cast<double>(ny * (ny + 1))/2 - rankcy_sum;
                double u = std::min(u1, u2);

                double n = nx + ny;
                // normal approximation
                double mu = static_cast<double>(nx * ny) / 2;
                double corr = 0;
                for (const auto &[val,rks] : countsc) {
                    const auto &[rk, c] = rks;
                    if (c > 1)
                        corr += c * c * c - c;
                }
                double sd = sqrt(static_cast<double>(nx * ny) / n / (n - 1)) * sqrt((n * n * n - n)/12 - corr / 12);

                double z = abs((u - mu) / sd);
                p = boost::math::cdf(boost::math::complement(norm, z)) * 2;
            }

            return std::make_tuple(p, rankcx_mean, rankcy_mean);
        };
    } else {
        throw std::runtime_error("Test not implemented");
    }

    utils::call_rows<std::unique_ptr<bit_vector>,
                     std::unique_ptr<ValuesContainer>,
                     PairContainer>(columns_all, column_values_all,
                                    [&](uint64_t row_i, const auto &raw_row) {
        PairContainer row;
        if (!ignored[row_i]) {
            for (const auto &[j, raw_c] : raw_row) {
                if (uint64_t c = adjust_count(raw_c, j, row_i))
                    row.emplace_back(j, c);
            }
        }

        auto [pval, in_stat, out_stat] = compute_pval(row);
        pvals.emplace_back(pval);

        if (pval > 1.1 || in_stat == out_stat)
            return;

        bool in_kmer = in_stat > out_stat;
        bool out_kmer = in_stat < out_stat;
        uint64_t k = std::numeric_limits<uint64_t>::max();

        if (config.test_type != "nonparametric") {
            double pval_min = 0;
            if (config.test_type == "nonparametric_u") {
                pval_min = 2.0 / exp(lgamma(groups.size() + 1) - lgamma(labels_in.size() + 1) - lgamma(labels_out.size() + 1));
            } else {
                PairContainer best;
                for (const auto &[j, c] : row) {
                    if ((in_kmer && !groups[j]) || (out_kmer && groups[j])) {
                        best.emplace_back(j, max_obs_vals[j]);
                    } else {
                        best.emplace_back(j, 1);
                    }
                }

                auto [pm, in_sum_m, out_sum_m] = compute_pval(best);
                pval_min = pm;
            }

            if (pval_min > pval) {
                common::logger->error("Min p-val estimate too high: min {} > cur {}", pval_min, pval);
                throw std::runtime_error("Test failed");
            }

            if (pval_min >= 0.05)
                return;

            if (pval_min > 0) {
                double lkd = log2(0.05) - log2(pval_min);
                if (lkd <= 64)
                    k = pow(2.0, lkd);

                if (k == 0) {
                    common::logger->error("k: {}\tlog k: {}\tpval_min: {}\tpval: {}", k, lkd, pval_min, pval);
                    throw std::runtime_error("Min failed");
                }
            }
        }

        node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
        m[k].emplace_back(node);
        auto &indicator = in_kmer ? indicator_in : indicator_out;
        if (pval < 0.05)
            indicator[node] = true;
    });

    for (auto &col : columns_all) {
        col.reset();
    }

    for (auto &col_vals : column_values_all) {
        col_vals.reset();
    }

    uint64_t total_sig = sdsl::util::cnt_one_bits(indicator_in) + sdsl::util::cnt_one_bits(indicator_out);
    auto m_data = const_cast<std::vector<std::pair<size_t, std::vector<uint64_t>>>&&>(m.values_container());

    std::sort(m_data.begin(), m_data.end(), utils::LessFirst());
    size_t acc = std::accumulate(m_data.begin(), m_data.end(), size_t(0),
                                 [](size_t sum, const auto &a) { return sum + a.second.size(); });
    size_t total_tests = acc;
    common::logger->trace("Correcting {}/{} significant p-values", total_sig, acc);

    if (total_sig) {
        auto begin = m_data.begin();
        size_t k = 0;
        size_t last_k = 0;
        for ( ; begin != m_data.end(); ++begin) {
            const auto &[cur_k, bucket] = *begin;
            common::logger->trace("k: {}\tm(k): {}", cur_k, acc);
            if (acc <= cur_k) {
                k = std::max(acc, last_k + 1);
                break;
            }

            last_k = cur_k;
            acc -= bucket.size();
        }

        if (k == 0) {
            common::logger->trace("No good k value found, defaulting to Bonferroni cutoff");
            k = total_tests;
            common::logger->trace("k: {}\talpha_corr: {}", k, 0.05 / k);
            for (auto jt = m_data.begin(); jt != m_data.end(); ++jt) {
                const auto &[cur_k, bucket] = *jt;
                for (uint64_t i : bucket) {
                    if (pvals[i] < 0.05 && pvals[i] * k >= 0.05) {
                        indicator_in[i] = false;
                        indicator_out[i] = false;
                    }
                }
            }
        } else {
            common::logger->trace("k: {}\talpha_corr: {}", k, 0.05 / k);
            for (auto jt = m_data.begin(); jt != m_data.end(); ++jt) {
                const auto &[cur_k, bucket] = *jt;
                for (uint64_t i : bucket) {
                    if (jt < begin || pvals[i] * k >= 0.05) {
                        indicator_in[i] = false;
                        indicator_out[i] = false;
                    }
                }
            }
        }
    }

    auto masked_graph_in = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr, std::make_unique<bitmap_vector>(std::move(indicator_in)), true,
        is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
    );
    masked_graph_in->likelihood_ratios = pvals;

    auto masked_graph_out = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr, std::make_unique<bitmap_vector>(std::move(indicator_out)), true,
        is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
    );
    masked_graph_out->likelihood_ratios = std::move(pvals);

    return std::make_pair(masked_graph_in, masked_graph_out);
}



typedef std::function<size_t()> LabelCountCallback;

constexpr std::memory_order MO_RELAXED = std::memory_order_relaxed;

uint64_t atomic_fetch(const sdsl::int_vector<> &vector,
                      uint64_t i,
                      std::mutex &backup_mutex,
                      int mo) {
    const uint8_t width = vector.width();
    if (width == 64)
        return __atomic_load_n(vector.data() + i, mo);

    if (width == 32)
        return __atomic_load_n((uint32_t*)vector.data() + i, mo);

    if (width == 16)
        return __atomic_load_n((uint16_t*)vector.data() + i, mo);

    if (width == 8)
        return __atomic_load_n((uint8_t*)vector.data() + i, mo);

    std::lock_guard<std::mutex> lock(backup_mutex);
    return vector[i];
}

uint64_t atomic_fetch_and_add(sdsl::int_vector<> &vector,
                              uint64_t i,
                              uint64_t val,
                              std::mutex &backup_mutex,
                              int mo) {
    const uint8_t width = vector.width();
    if (width == 64)
        return __atomic_fetch_add(vector.data() + i, val, mo);

    if (width == 32)
        return __atomic_fetch_add((uint32_t*)vector.data() + i, val, mo);

    if (width == 16)
        return __atomic_fetch_add((uint16_t*)vector.data() + i, val, mo);

    if (width == 8)
        return __atomic_fetch_add((uint8_t*)vector.data() + i, val, mo);

    std::lock_guard<std::mutex> lock(backup_mutex);
    uint64_t old_val = vector[i];
    vector[i] += val;
    return old_val;
}


/**
 * Return an int_vector<>, bit_vector of lengths anno_graph.get_graph().max_index() * 2
 * and anno_graph.get_graph().max_index(), respectively, and a bit_vector of length
 * equal to the total number of labels.
 * For an index i, the int_vector at indices 2*i and 2*i + 1 represent the
 * number of labels in labels_in and labels_out which the k-mer of index i is
 * annotated with, respectively. The width of the int_vector<> is computed to be
 * wide enough to contain counts up to num_labels.
 * The returned bit_vector is a k-mer mask indicating those k-mers annotated
 * with at least one in-label. If add_out_labels_to_mask is true, then it also indicates
 * which those k-mers with at least one out-label.
 * The second bit vector indicates which of the labels correspond to an "other"
 * label (i.e., not in- or out-).
 */


using ValueCallback = std::function<void(uint64_t /* row index */, uint64_t /* row value */)>;
using ValueGenerator = std::function<void(const ValueCallback&)>;
using ColumnCallback = std::function<void(uint64_t /* offset */, const Label& /* annotation column */, const ValueGenerator&)>;
using ColumnGenerator = std::function<void(const ColumnCallback&)>;

template <typename T>
std::string outstring(std::vector<T> out_array);
template <typename T>
std::string outstring(std::vector<T> out_array){
    std::string out = "";
    for (auto &item: out_array) { std::string x = std::to_string(item); out = out + x + ", "; }
    return out;
}

void kmer_distribution_table(const ColumnGenerator &generate_columns,
                             const tsl::hopscotch_set<Label> &labels_in,
                             const tsl::hopscotch_set<Label> &labels_out,
                             size_t max_index,
                             MaskedDeBruijnGraph &masked_graph, // std::shared_ptr<const DeBruijnGraph> graph_ptr,
                             const DifferentialAssemblyConfig &config,
                             size_t num_threads);  // For kmer distribution.

std::tuple<sdsl::int_vector<>, sdsl::bit_vector, sdsl::bit_vector, std::tuple<size_t, size_t>>
construct_diff_label_count_vector(const ColumnGenerator &generate_columns,
                                  const std::vector<std::string> &files,
                                  const tsl::hopscotch_set<Label> &labels_in,
                                  const tsl::hopscotch_set<Label> &labels_out,
                                  const DifferentialAssemblyConfig &config,
                                  size_t max_index,
                                  size_t num_labels,
                                  size_t num_threads,
                                  bool add_out_labels_to_mask);

// Regions of a graph mask which should be kept (i.e., masked in)
typedef std::vector<std::pair<size_t, size_t>> Intervals;

// Returns a vector of kept regions given a unitig and its corresponding path
typedef std::function<Intervals(const std::string&, const std::vector<node_index>&)> GetKeptIntervals;

// Compare the first element of two pairs (used for sorting)
bool sortByPValue(const std::pair<int, double>& a, const std::pair<int, double>& b) {
    return a.first < b.first; // Sort in ascending order of p-values
}

// Function for calculating
// the median
template <class Container>
double findMedian(const Container &a,
                    int threshold) {
    logger->trace("non zero kmers: {}", a.size());
    std::vector<double> a_filtered;
    for (size_t i = 0; i < a.size(); i++){
        int value = a[i];
        if (value > threshold){
            a_filtered.push_back(value);
        }
    }
    // auto it = std::copy_if(a.begin(), a.end(), std::back_inserter(a_filtered), [&](int i) {
    //     return (i > threshold);
    // });
    logger->trace("kmers above threshold: {}", a_filtered.size());
    int n = a_filtered.size();
    // If size of the arr[] is even
    if (n % 2 == 0) {

        // Applying nth_element
        // on n/2th index
        std::nth_element(a_filtered.begin(),
                    a_filtered.begin() + n / 2,
                    a_filtered.end());

        // Applying nth_element
        // on (n-1)/2 th index
        std::nth_element(a_filtered.begin(),
                    a_filtered.begin() + (n - 1) / 2,
                    a_filtered.end());

        // Find the average of value at
        // index N/2 and (N-1)/2
        return (double)(a_filtered[(n - 1) / 2]
                        + a_filtered[n / 2])
               / 2.0;
    }

    // If size of the arr[] is odd
    else {

        // Applying nth_element
        // on n/2
        std::nth_element(a_filtered.begin(),
                    a_filtered.begin() + n / 2,
                    a_filtered.end());

        // Value at index (N/2)th
        // is the median
        return (double)a_filtered[n / 2];
    }
}

bool filtering_row(std::vector<double> in_counts_non_zero,
                    std::vector<double> out_counts_non_zero,
                    int total_samples,
                    std::shared_ptr<const DeBruijnGraph> graph_ptr,
                    int64_t row_id);

std::shared_ptr<MaskedDeBruijnGraph> brunner_munzel_test(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                                                        std::shared_ptr<MaskedDeBruijnGraph> masked_graph,
                                                        const std::vector<std::string> &files,
                                                        const tsl::hopscotch_set<Label> &labels_in,
                                                        const tsl::hopscotch_set<Label> &labels_out,
                                                        const DifferentialAssemblyConfig &config,
                                                        size_t total_unitigs,
                                                        std::vector<double> medians) {

    const auto &in_mask = static_cast<const bitmap_vector &>(masked_graph->get_mask()).data();
    sdsl::bit_vector mask = in_mask;
    // brunner munzel test
    // call kmer matrix per row
    std::vector<std::unique_ptr<bit_vector>> columns_all;
    std::vector<std::unique_ptr<sdsl::int_vector_buffer<>>> column_values_all;
    std::vector<bool> column_label;
    annot::ColumnCompressed<>::load_columns_and_values(files,
        [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> && column, sdsl::int_vector_buffer<> && column_values) {
            columns_all.push_back(std::move(column));
            column_values_all.push_back(std::make_unique<sdsl::int_vector_buffer<>>(std::move(column_values)));
            if (labels_in.find(label) != labels_in.end()) // if the label is in the in labels set column_label to 1, else 0
                column_label.push_back(1);
            else if (labels_out.find(label) != labels_out.end())
                column_label.push_back(0);
            else throw std::runtime_error("Label not found in labels_in or labels_out");
        }
    );
    int64_t num_tests = 0;
    int64_t kept_nodes = 0;
    // calculate the p-value for the brunner munzel test
    auto statistical_model = DifferentialTest(config.family_wise_error_rate, total_unitigs, 0, 0, 0);
    utils::call_rows<std::unique_ptr<bit_vector>,
                    std::unique_ptr<sdsl::int_vector_buffer<>>,
                    std::vector<std::pair<uint64_t, uint8_t>>>(columns_all, column_values_all,
                                                               [&](uint64_t row_id, const auto &row) {
                        // auto kmer_string = graph_ptr->get_node_sequence(AnnotatedDBG::anno_to_graph_index(row_id));
                        if (row.size() == 0) { // if the kmer is not present in any sample
                            mask[AnnotatedDBG::anno_to_graph_index(row_id)] = false;
                            return;
                        }
                        else{
                            // populate the vector for the brunner munzel test with the kmer counts for each sample
                            std::vector<double> in_counts;
                            std::vector<double> out_counts;
                            for (size_t i = 0; i < row.size(); i++){
                                double median = medians[row[i].first];
                                if (median == 0) median = 1; // if the median is 0, set it to 1 (to avoid division by 0)
                                if (column_label[row[i].first])
                                    in_counts.push_back(row[i].second/median);
                                else
                                    out_counts.push_back(row[i].second/median);
                            }
                            // add 0 to in_counts and out_counts if the kmer is not present in the sample
                            if (filtering_row(in_counts, out_counts, files.size(), graph_ptr, row_id)) {
                                num_tests++;
                                if (in_counts.size() < labels_in.size())
                                    in_counts.resize(labels_in.size(), 0);
                                if (out_counts.size() < labels_out.size())
                                    out_counts.resize(labels_out.size(), 0);
                                // run test and adjust masked graph accordingly
                                auto [keep, pvalue] = statistical_model.brunner_munzel_test(in_counts, out_counts);
                                if (keep == false)
                                    mask[AnnotatedDBG::anno_to_graph_index(row_id)] = false;
                                else{
                                    kept_nodes++;
                                }
                            }
                            else{
                                mask[AnnotatedDBG::anno_to_graph_index(row_id)] = false;
                            }

                        }
    });
    logger->trace("number of tests with 0 var: {}", statistical_model.var_0);
    logger->trace("number of tests: {}", num_tests);
    masked_graph->set_mask(new bitmap_vector(std::move(mask)));
    logger->trace("number of kept nodes: {}", kept_nodes);
    return masked_graph;
}
// returns true if the kmer passes the filtering criteria
bool filtering_row(std::vector<double> in_counts_non_zero,
                    std::vector<double> out_counts_non_zero,
                    int total_samples,
                    std::shared_ptr<const DeBruijnGraph> graph_ptr,
                    int64_t row_id){
    // minimum present in 20% of samples
    bool is_present_in_20_percent = (in_counts_non_zero.size() + out_counts_non_zero.size()) > 0.2 * double(total_samples);
    if (is_present_in_20_percent == false)
        return false;
    // present more in in than out
    bool is_present_more_in_foreground = in_counts_non_zero.size() > out_counts_non_zero.size();
    // low complexity filter
    if (is_present_more_in_foreground == false)
        return false;
    std::string kmer_string = graph_ptr->get_node_sequence(AnnotatedDBG::anno_to_graph_index(row_id));
    bool high_complexity =  !is_low_complexity(kmer_string);
    return (high_complexity);
}
std::vector<double> filtering_column(const std::vector<std::string> &files,
               const tsl::hopscotch_set<Label> &labels_in,
               const tsl::hopscotch_set<Label> &labels_out,
               const DifferentialAssemblyConfig &config,
               const ColumnCallback &column_callback
               ) {
    // median normalisation
    // threshold of 2
    double threshold = config.min_count; // TODO get from config
    assert(files.size() > 0);
    // generate a vector of the total number of kmers per sample for each group (in and out)
    std::vector<double> medians;
    logger->trace("filtering_column");
    annot::ColumnCompressed<>::load_columns_and_values(files,
        [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> && column, sdsl::int_vector_buffer<> && column_values) { // goes through all files in a group
            uint64_t max_col_width = std::pow(2, 16); // max value for int16_t
            column_callback(offset, label, [&](const ValueCallback &value_callback) {
                assert(std::accumulate(column_values.begin(), column_values.end(), 0) > 0); // at least one kmer in sample
                assert(column->num_set_bits() > 0); // at least one kmer present in sample
                assert(column->num_set_bits() == column_values.size()); // same number of kmers in column as in column_values
                call_ones(*column, [&](uint64_t i) {
                    uint64_t column_value = column_values[column->rank1(i)-1];
                    if (column_value > threshold){
                        column_value = std::min(column_value, max_col_width); // normalize value
                        value_callback(i, column_value);
                    }
                });
                double median = findMedian(std::move(column_values), threshold);
                medians.push_back(median);
                logger->trace("median: {}", median);
            });
        }
    );
    return medians;
}

// Assemble unitigs from the masked graph, then use get_kept_intervals to generate
// regions which should be masked in. Update the graph mask accordingly.
void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const GetKeptIntervals &get_kept_intervals,
                                   size_t num_threads);

// Given an initial mask and counts, generate a masked graph. If add_complement
// is true, then add the reverse complements of all nodes to the graph as well.
std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                          sdsl::int_vector<> &counts,
                          sdsl::bit_vector&& mask,
                          bool add_complement,
                          size_t num_threads);

// Helper function for calls in interface
std::shared_ptr<MaskedDeBruijnGraph>
mask_nodes_by_label(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                    const AnnotatedDBG *anno_graph,
                    sdsl::int_vector<>&& counts,
                    sdsl::bit_vector&& mask,
                    sdsl::bit_vector&& other_mask,
                    size_t num_labels,
                    const std::vector<std::string> &files,
                    const tsl::hopscotch_set<Label> &labels_in,
                    const tsl::hopscotch_set<Label> &labels_out,
                    const tsl::hopscotch_set<Label> &labels_in_round2,
                    const tsl::hopscotch_set<Label> &labels_out_round2,
                    std::tuple<size_t, size_t> total_kmers,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads,
                    std::vector<double> medians);

std::shared_ptr<MaskedDeBruijnGraph> // This version of mask_nodes_by_label is for when the columns are already loaded.
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const tsl::hopscotch_set<Label> &labels_in,
                    const tsl::hopscotch_set<Label> &labels_out,
                    const tsl::hopscotch_set<Label> &labels_in_round2,
                    const tsl::hopscotch_set<Label> &labels_out_round2,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads,
                    size_t num_parallel_files) {
    num_parallel_files = std::min(num_threads, num_parallel_files);

    auto graph_ptr = std::static_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );

    bool check_other = config.label_mask_other_unitig_fraction != 1.0;
    bool unitig_mode = check_other || labels_in_round2.size() || labels_out_round2.size()
            || config.label_mask_in_unitig_fraction != 0.0
            || config.label_mask_out_unitig_fraction != 1.0
            || config.label_mask_other_unitig_fraction != 1.0;

    bool add_complement = graph_ptr->get_mode() == DeBruijnGraph::CANONICAL
        && (config.add_complement || unitig_mode);
    logger->trace("add_complement: {}", add_complement);

    logger->trace("Generating initial mask");

    auto count_vector = construct_diff_label_count_vector(
        [&](const ColumnCallback &column_callback) {
            if (config.count_kmers) {
                throw std::runtime_error("AnnotatedDBG with counts not supported yet");
            } else {
                size_t offset = 0;
                for (const auto &label : labels_in) {
                    column_callback(offset, label, [&](const ValueCallback &value_callback) {
                        anno_graph.call_annotated_nodes(label, [&](node_index i) {
                            i = AnnotatedDBG::graph_to_anno_index(i);
                            value_callback(i, 1);
                        });
                    });
                    ++offset;
                }
                for (const auto &label : labels_out) {
                    column_callback(offset, label, [&](const ValueCallback &value_callback) {
                        anno_graph.call_annotated_nodes(label, [&](node_index i) {
                            i = AnnotatedDBG::graph_to_anno_index(i);
                            value_callback(i, 1);
                        });
                    });
                    ++offset;
                }
            };
        },
        { }, labels_in, labels_out, config, graph_ptr->max_index(),
        std::max(labels_in.size(), labels_out.size()), num_parallel_files,
        add_complement
    );

    auto &[counts, init_mask, other_labels, total_kmers] = count_vector; //ERROR  only 3 names provided for structured binding, std::vector<long unsigned int, std::allocator<long unsigned int> > >’ decomposes into 4 elements
    sdsl::bit_vector other_mask(init_mask.size() * check_other, false);
    return mask_nodes_by_label(graph_ptr, &anno_graph,
                               std::move(counts), std::move(init_mask),
                               std::move(other_mask),
                               anno_graph.get_annotator().num_labels(),
                               {},
                               labels_in, labels_out,
                               labels_in_round2, labels_out_round2, total_kmers,
                               config, num_threads, {});
}

std::shared_ptr<DeBruijnGraph>
mask_nodes_by_label(std::shared_ptr<const DeBruijnGraph> graph_ptr, // Myrthe: this version i don't need to change at first
                    const std::vector<std::string> &files,
                    const tsl::hopscotch_set<Label> &labels_in,
                    const tsl::hopscotch_set<Label> &labels_out,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads,
                    size_t num_parallel_files) {
    bool check_other = config.label_mask_other_unitig_fraction != 1.0;
    bool unitig_mode = check_other
            || config.label_mask_in_unitig_fraction != 0.0
            || config.label_mask_out_unitig_fraction != 1.0
            || config.label_mask_other_unitig_fraction != 1.0;

    bool add_complement = graph_ptr->get_mode() == DeBruijnGraph::CANONICAL
        && (config.add_complement || unitig_mode);
    logger->trace("add_complement: {}", add_complement);
    logger->trace("Generating initial mask");
    bool filter = config.filter;
    std::vector<double> medians;
    auto count_vector = construct_diff_label_count_vector( // called for each group (in and out)
            [&](const ColumnCallback &column_callback) {
            if (config.count_kmers) {
                if (filter){
                    medians = filtering_column(files, labels_in, labels_out, config, column_callback);
                }
                else{
                    assert(files.size() > 0);
                    // generate a vector of the total number of kmers per sample for each group (in and out)
                    annot::ColumnCompressed<>::load_columns_and_values(files,
                        [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> && column, sdsl::int_vector_buffer<> && column_values) { // goes through all files in a group
                            int max_col_width = std::pow(2, 16); // max value for int16_t
                            column_callback(offset, label, [&](const ValueCallback &value_callback) {
                                assert(std::accumulate(column_values.begin(), column_values.end(), 0) > 0); // at least one kmer in sample
                                assert(column->num_set_bits() > 0); // at least one kmer present in sample
                                assert(column->num_set_bits() == column_values.size()); // same number of kmers in column as in column_values
                                call_ones(*column, [&](uint64_t i) {
                                    uint64_t column_value = column_values[column->rank1(i)-1];
                                    column_value = std::min((int) column_value, max_col_width); // normalize value
                                    if (column_value > config.min_count)
                                        value_callback(i, column_value);
                                });
                            });
                        }
                    );
                }
            } else {
                if (!config.clean) {
                    annot::ColumnCompressed<>::merge_load(files,
                        [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> && column) {
                            column_callback(offset, label, [&](const ValueCallback &value_callback) {
                                call_ones(*column, [&](uint64_t i) {
                                    value_callback(i, 1);
                                });
                            });
                        }
                    );
                } else {
                    annot::ColumnCompressed<>::load_columns_and_values(files,
                        [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> && column, sdsl::int_vector<>&& column_values) {
                            column_callback(offset, label, [&](const ValueCallback &value_callback) {
                                call_ones(*column, [&](uint64_t i) {
                                    auto val = column_values[column->rank1(i) - 1];
                                    if (val > config.min_count)
                                        value_callback(i, 1);
                                });
                            });
                            // size_t min_count = config.min_count;

                            // auto mask = std::make_unique<bitmap_lazy>([&](node_index node) {
                            //     if (!node)
                            //         return false;

                            //     return (*column)[AnnotatedDBG::graph_to_anno_index(node)];
                            // }, column->size() + 1);
                            // bool is_primary = (graph_ptr->get_mode() == DeBruijnGraph::PRIMARY);
                            // std::shared_ptr<DeBruijnGraph> masked
                            //     = std::make_shared<MaskedDeBruijnGraph>(graph_ptr, std::move(mask), true,
                            //                            is_primary
                            //                             ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC);

                            // std::shared_ptr<CanonicalDBG> wrapped;
                            // if (is_primary) {
                            //     wrapped = std::make_shared<CanonicalDBG>(masked);
                            //     masked = wrapped;
                            // }
                            // column_callback(offset, label, [&](const ValueCallback &value_callback) {
                            //     assert(column->num_set_bits() == column_values.size());
                            //     // call_ones(*column, [&](uint64_t i) {
                            //     //     auto val = column_values[column->rank1(i) - 1];
                            //     //     if (val >= min_count && val < max_count)
                            //     //         value_callback(i, 1);
                            //     // });
                            //     masked->call_unitigs([&](const std::string&, const auto &path) {
                            //         std::vector<uint64_t> counts;
                            //         std::vector<uint64_t> rows;
                            //         counts.reserve(path.size());
                            //         rows.reserve(path.size());
                            //         for (node_index node : path) {
                            //             auto row = AnnotatedDBG::graph_to_anno_index(
                            //                 is_primary ? wrapped->get_base_node(node) : node
                            //             );
                            //             rows.emplace_back(row);
                            //             counts.emplace_back(column_values[column->rank1(row) - 1]);
                            //         }
                            //         std::sort(counts.begin(), counts.end());
                            //         double median = counts[counts.size() / 2];

                            //         if (counts.size() % 2 == 1)
                            //             median = (median + counts[counts.size() / 2 - 1])/2;

                            //         if (median >= min_count && median < max_count) {
                            //             for (auto row : rows) {
                            //                 value_callback(row, 1);
                            //             }
                            //         }
                            //     });
                            // });
                        }
                    );
                }
            }
        },
        files, labels_in, labels_out, config, graph_ptr->max_index(),
        std::max(labels_in.size(), labels_out.size()), num_parallel_files,
        config.test_by_unitig
    );

    auto &[counts, init_mask, other_labels, total_kmers] = count_vector;
    // total_kmers: tuple of in_total_kmers and out_total_kmers (number of kmers in foreground and background)
    // init_mask: bit_vector of length graph_ptr->max_index() indicating which nodes are in the foreground (size = number of unique kmers total)
    // counts: int_vector of length graph_ptr->max_index() * 2 indicating the number of in and out labels for each node (size = number of unique kmers total * 2)
    // other_labels: bit_vector of length num_labels indicating which labels are not in or out labels
/*
    if (config.evaluate_assembly){     // Myrthe temporary: evaluate what k-mers overlap

        logger->trace("Evaluate confusion table");
        assert(counts.size() > 1);
        assert( std::accumulate(counts.begin(), counts.end(), 0) > 0);
        int true_positive = 0; int true_negative = 0; int false_positive = 0; int false_negative =0;
        for (size_t i = 0; i < counts.size() ; i+=2){
            auto &[counts, init_mask, other_labels, total_kmers] = count_vector;
            uint64_t true_count = counts[i];
            uint64_t found_count = counts[i + 1];
            if (true_count == found_count){
                true_positive += true_count; // if 0 nothing will be added.
            }else if (true_count > found_count){
                false_negative += (true_count - found_count);
                true_positive += found_count; // add to true_positive
            }else{false_positive += (found_count - true_count);
                true_positive += true_count;} // TODO :Why are there no false positives?
        };
        // true_negatives: total number of possible k-mers? Total number of k-mers in the original foreground samples from which the result had been calculated ?
        auto precision = true_positive / ( true_positive - true_negative); // accuracy/precision = Ratio of true positives to total predicted positives.
        auto recall = true_positive / (true_positive + false_negative);  //recall / sensitivity:  ratio of true positives to total (actual) positives in the data.
        //auto specificity = true_negative / (true_negative + false_positive); // Specificity is the Ratio of true negatives to total negatives in the data.
        std::cout<< "precision" << precision << std::endl << std::flush;
        std::cout<< "recall" << recall << std::endl << std::flush;
        std::cout << std::accumulate(counts.begin(), counts.end(), 0) << std::flush;
        assert(true_positive + false_negative + false_positive == std::accumulate(counts.begin(), counts.end(), 0));
    } // in mask probably only has sequences
*/

    sdsl::bit_vector other_mask(init_mask.size() * check_other, false);
    if (check_other && sdsl::util::cnt_one_bits(other_labels)) { // currently not used
        bool parallel = num_parallel_files > 1;
        size_t j = 0;
        std::atomic_thread_fence(std::memory_order_release);
        for (size_t i = 0; i < files.size(); ++i) {
            size_t num_labels_per_file = annot::ColumnCompressed<>::read_label_encoder(files[i]).size();
            if (count_ones(other_labels, j, j + num_labels_per_file)) {
                annot::ColumnCompressed<>::merge_load({ files[i] },
                    [&](size_t offset, const auto &, auto&& column) {
                        auto &[counts, init_mask, other_labels, total_kmers] = count_vector;
                        if (other_labels[j + offset]) {
                            call_ones(init_mask, [&](size_t i) { // check how many set bits are also in other labels.
                                if ((*column)[AnnotatedDBG::graph_to_anno_index(i)])
                                    set_bit(other_mask.data(), i, parallel, MO_RELAXED);
                            });
                        }
                    },
                    num_parallel_files
                );
            }

            j += num_labels_per_file;
        }
        std::atomic_thread_fence(std::memory_order_acquire);
    }

/*
    if (config.family_wise_error_rate == 0){  // Temporary code to make k-mer distribution
        logger->trace("K-mer distribution");

        auto masked_graph = make_initial_masked_graph(graph_ptr, counts, std::move(init_mask),
                                                      add_complement, num_threads);
        kmer_distribution_table(
                [&](const ColumnCallback &column_callback) {
                    annot::ColumnCompressed<>::load_columns_and_values(files,
                       [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> && column, auto&& column_values) {
                           column_callback(offset, label, [&](const ValueCallback &value_callback) {

                               if (config.test_by_unitig == false) {
                                   call_ones(*column, [&](uint64_t i) {
                                       value_callback(i, column_values[column->rank1(i)- 1]);
                                   });
                               }
                               else {
                                   for(size_t i = 0; i < column->size(); i++){
                                       double value = (*column)[i];
                                       if (value) value = column_values[column->rank1(i)- 1];
                                       value_callback(i, value);
                                   }
                               }
                           });
                       }
                    );
                },
                labels_in, labels_out, graph_ptr->max_index(), *masked_graph, config, num_threads
        );
        std::exit(EXIT_SUCCESS);
    }
*/
    auto masked_graph = mask_nodes_by_label(graph_ptr, nullptr,
                               std::move(counts), std::move(init_mask),
                               std::move(other_mask),
                               files.size(),
                               files,
                               labels_in, labels_out, {}, {},
                               total_kmers,
                               config, num_threads, medians);

    if (masked_graph->get_mode() != DeBruijnGraph::PRIMARY)
        return masked_graph;

    return std::make_shared<CanonicalDBG>(masked_graph);
}

// does differential testing
std::shared_ptr<MaskedDeBruijnGraph> // TODO add all columns and columns label as param
mask_nodes_by_label(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                    const AnnotatedDBG *anno_graph,
                    sdsl::int_vector<>&& counts,
                    sdsl::bit_vector&& init_mask,
                    sdsl::bit_vector&& other_mask,
                    size_t num_labels,
                    const std::vector<std::string> &files,
                    const tsl::hopscotch_set<Label> &labels_in,
                    const tsl::hopscotch_set<Label> &labels_out,
                    const tsl::hopscotch_set<Label> &labels_in_round2,
                    const tsl::hopscotch_set<Label> &labels_out_round2,
                    std::tuple<size_t, size_t> total_kmers,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads,
                    std::vector<double> medians) {
    // in and out counts are stored interleaved in the counts vector
    logger->trace("aaaaaaaaaaaaaaaaaaaaaA counts {}", counts.size());
    assert(counts.size() == init_mask.size() * 2);

    bool parallel = num_threads > 1;
    bool check_other = config.label_mask_other_unitig_fraction != 1.0;
    bool unitig_mode = check_other || labels_in_round2.size() || labels_out_round2.size()
            || config.label_mask_in_unitig_fraction != 0.0
            || config.label_mask_out_unitig_fraction != 1.0
            || config.label_mask_other_unitig_fraction != 1.0;

    size_t num_in_labels = labels_in.size() + labels_in_round2.size(); // number of in labels before filtering
    size_t num_out_labels = labels_out.size() + labels_out_round2.size(); // number of out labels before filtering

    bool add_complement = graph_ptr->get_mode() == DeBruijnGraph::CANONICAL
            && (config.add_complement || unitig_mode);
    logger->trace("add_complement: {}", add_complement);
    auto masked_graph = make_initial_masked_graph(graph_ptr, counts, std::move(init_mask),
                                                  add_complement, num_threads);

    auto mask_or = [&](sdsl::bit_vector &a, const sdsl::bit_vector &b,
                       const std::vector<node_index> &id_map) {
        call_ones(b, [&](size_t i) {
            if (id_map[i])
                set_bit(a.data(), id_map[i], parallel, MO_RELAXED);
        });
    };

    // check all other labels and round 2 labels
    if (anno_graph && (check_other || labels_in_round2.size() || labels_out_round2.size())) { // currently not used? only in binary?

        sdsl::bit_vector union_mask
                = static_cast<const bitmap_vector &>(masked_graph->get_mask()).data();
        std::mutex vector_backup_mutex;
        std::atomic_thread_fence(std::memory_order_release);

        auto count_merge
                = [&](sdsl::bit_vector
                              &a, // Myrthe ANSWER: what does count_merge do? Is a hack, I do not have to touch this now. Count_add would be better name. If you have a dense column, and include it in foreground or background. Wate time, Instead first create mask.  You only access dense columns, for instance a reference.
                      const sdsl::bit_vector &b, const std::vector<node_index> &id_map,
                      size_t offset = 0) {
                      call_ones(b, [&](size_t i) {
                          if (id_map[i]) {
                              set_bit(a.data(), id_map[i], parallel, MO_RELAXED);
                              atomic_fetch_and_add(counts, id_map[i] * 2 + offset, 1,
                                                   vector_backup_mutex, MO_RELAXED);
                          }
                      });
                  };
        logger->trace("Checking shared and other labels");
        masked_graph->call_sequences(
                [&](const std::string &contig, const std::vector<node_index> &path) {
                    for (const auto &[label, sig] :
                         anno_graph->get_top_label_signatures(contig, num_labels)) {
                        bool found_in = labels_in.count(label);
                        bool found_out = labels_out.count(label);
                        bool found_in_round2 = labels_in_round2.count(label);
                        bool found_out_round2 = labels_out_round2.count(label);
                        if (!found_in && !found_out && !found_in_round2
                            && !found_out_round2 && check_other) {
                            mask_or(other_mask, sig, path);
                        }

                        if (found_in_round2)
                            count_merge(union_mask, sig, path);

                        if (found_out_round2)
                            count_merge(union_mask, sig, path, 1);
                    }
                },
                num_threads);

        std::atomic_thread_fence(std::memory_order_acquire);

        masked_graph->set_mask(new bitmap_vector(std::move(union_mask)));
    }


    // Filter unitigs from masked graph based on filtration criteria
    logger->trace("Filtering out background");
    if (config.count_kmers) {
        auto &[in_total_kmers, out_total_kmers] = total_kmers;
        logger->trace("I'm here");
        if (config.test_by_unitig == false) { // Statistical testing part when k-mer counts are included.
            std::atomic<uint64_t> total_unitigs(0);
            if (config.test_type == "brunner_munzel"){
                logger->trace("Brunner Munzel test");
                masked_graph->likelihood_ratios.resize(counts.size()/2+1);
                masked_graph->call_unitigs([&](const std::string &, const std::vector<node_index> &) {
                    ++total_unitigs;
                    // total_unitigs.fetch_add(1, MO_RELAXED);
                }, get_num_threads());
                masked_graph = brunner_munzel_test(graph_ptr, masked_graph,files,labels_in,labels_out,config,total_unitigs, medians);
                masked_graph->likelihood_ratios.resize(counts.size()/2+1);
            }
            else if (config.test_type == "likelihoodratio"){
                const auto &in_mask = static_cast<const bitmap_vector &>(masked_graph->get_mask()).data();
                total_unitigs = static_cast<const bitmap_vector &>(masked_graph->get_mask()).num_set_bits();
                logger->trace("Likelihood ratio test");
                logger->trace("ooooo total unitigs {}", total_unitigs); // the total number of hypotheses tested: number of unitigs, to accounts for dependence
                auto statistical_model = DifferentialTest(config.family_wise_error_rate, total_unitigs,
                                                        std::min((int) std::distance(counts.begin(), std::max_element(counts.begin(), counts.end())), (int) 100000),
                                                        in_total_kmers, out_total_kmers);
                logger->trace("test_by_unitig == false");
                sdsl::bit_vector mask = in_mask;
                size_t total_nodes = masked_graph->num_nodes();
                size_t kept_nodes = 0;
                //std::vector<int> nodes;
                std::vector<std::pair<double,int>> likelihood_ratios_nodes;
                call_ones(in_mask, [&](node_index node) { // only run if counts_case > counts_control
                    uint64_t in_sum = counts[node * 2];
                    uint64_t out_sum = counts[node * 2 + 1];
                    // filter kmers in samples less noise
                    // if (in_sum + out_sum >= int(total_nodes/20)) { // check for minimum occurance of k-mer
                        double out_sum_normalized = (double) out_sum * in_total_kmers / out_total_kmers;
                        if (in_sum > out_sum_normalized) { // check that kmer occurs more in foreground than in background
                            // check for high complexity
                            // auto kmer_string = graph_ptr->get_node_sequence(node);
                            // auto complexity_low = is_low_complexity(kmer_string);
                            // logger->trace("kmer_string: {}", kmer_string);
                            // logger->trace("complexity_low: {}", complexity_low);
                            auto [keep, likelihood_ratio] = statistical_model.likelihood_ratio_test(in_sum, out_sum);
                            // likelihood_ratios_nodes.push_back(std::make_pair(likelihood_ratio, node));
                            if (keep){
                                mask[node] = true;
                                kept_nodes++;
                            }
                            else
                                mask[node] = false;
                        }
                        else
                            mask[node] = false;
                    // }
                    // else
                    //     mask[node] = false;
                });
                // // get index for benjamini-yekutieli correction
                // std::sort(likelihood_ratios_nodes.begin(), likelihood_ratios_nodes.end(), sortByPValue);
                // // double alpha = statistical_model.get_t_test_alpha(statistical_model.get_df_approx(), 0.05) // 0.05 is the alpha value, TODO: add params.
                // int max_i = statistical_model.benjamini_yekutieli(likelihood_ratios_nodes, statistical_model.lrt_threshold());
                // // reject all hypotheses with index <= max_i
                // for (int i = 0; i <= max_i; i++){
                //     mask[likelihood_ratios_nodes[i].second] = false;
                // }
                // logger->trace("max_i {}", max_i);
                // logger->trace("likelihood_ratios__nodes_elem {}, {}, {}", likelihood_ratios_nodes[0].first, likelihood_ratios_nodes[1].first, likelihood_ratios_nodes[2].first);
                // logger->trace("likelihood_ratios_nodes_size {}", likelihood_ratios_nodes.size());
                masked_graph->set_mask(new bitmap_vector(std::move(mask)));
                masked_graph->likelihood_ratios.resize(counts.size()/2+1);
                logger->trace("Kept {} out of {} nodes", kept_nodes, total_nodes);
            }

        } else {
            logger->trace("test_by_unitig == true");
            masked_graph->likelihood_ratios.resize(counts.size()/2+1);
            std::fill(masked_graph->likelihood_ratios.begin(), masked_graph->likelihood_ratios.end(), 0);
            std::atomic<uint64_t> total_unitigs(0);
            masked_graph->call_unitigs([&](const std::string &, const std::vector<node_index> &) {
                total_unitigs.fetch_add(1, MO_RELAXED);
            });
            auto statistical_model = DifferentialTest(config.family_wise_error_rate, total_unitigs, // total number of hypotheses equals the number of unitigs
                                               std::min((int) std::distance(counts.begin(),
                                               std::max_element(counts.begin(), counts.end())), (int) 100000), // Myrthe: the factorial should increase for the unitig mode.
                                               in_total_kmers, out_total_kmers); // Myrthe: I don't know what to do with the last two arguments.
            update_masked_graph_by_unitig(*masked_graph, // intitial mask
                [&](const std::string &, const std::vector<node_index> &path) -> Intervals { // return a set of intervals to keep in the graph
                    int in_kmer_count_unitig = 0;
                    int out_kmer_count_unitig = 0;
                    for (size_t i = 0; i < path.size(); ++i) {
                        in_kmer_count_unitig += counts[path[i] * 2];
                        out_kmer_count_unitig += counts[path[i] * 2 + 1];
                    }
                    auto [significance, likelihood_ratio] = statistical_model.likelihood_ratio_test(in_kmer_count_unitig, out_kmer_count_unitig);

                    if (significance){
                        for (size_t i = 0; i < path.size(); ++i) { // if significant add the likelihood ratio or the effect (in_kmer_count_unitig - out_kmer_count_unitig) to the likelihood_ratios vector. This is a temporary solution
                            masked_graph->likelihood_ratios[path[i]] = in_kmer_count_unitig - out_kmer_count_unitig; //likelihood_ratio;
                        }
                        return { std::make_pair(0, path.size()) };
                    }
                    else
                        return {};
                },
            num_threads);
        }
            //std::fill(masked_graph->likelihood_ratios.begin(), masked_graph->likelihood_ratios.end(), 0.01); // Myrthe: For testing. temporary
        return masked_graph;
    }


    size_t min_label_in_count = std::ceil(config.label_mask_in_kmer_fraction
                                    * num_in_labels);
    size_t max_label_out_count = std::floor(config.label_mask_out_kmer_fraction
                                    * num_out_labels);

    if (config.test_by_unitig) {
        sdsl::bit_vector mask(masked_graph->get_mask().size(), false);

        std::shared_ptr<DeBruijnGraph> test_graph = masked_graph;
        std::shared_ptr<CanonicalDBG> wrapped;
        if (masked_graph->get_mode() == DeBruijnGraph::PRIMARY) {
            wrapped = std::make_shared<CanonicalDBG>(masked_graph);
            test_graph = wrapped;
        }

        test_graph->call_unitigs([&](const std::string &, const auto &path) {
            uint64_t total_in_count = 0;
            uint64_t total_out_count = 0;
            auto begin = path.end();
            auto end = path.begin();
            for (auto it = path.begin(); it != path.end(); ++it) {
                node_index node = wrapped ? wrapped->get_base_node(*it) : *it;
                auto in_count = counts[node * 2];
                total_in_count += in_count;
                total_out_count += counts[node * 2 + 1];
                if (in_count) {
                    begin = std::min(it, begin);
                    end = std::max(it, end);
                }
            }

            size_t path_size = end - begin;
            if (total_in_count >= min_label_in_count*path_size && total_out_count <= max_label_out_count*path_size) {
                std::for_each(begin, end, [&](node_index node) {
                    mask[wrapped ? wrapped->get_base_node(node) : node] = true;
                });
            }
        });

        /*
        size_t kept_nodes = 0;

        call_ones(masked_graph->get_mask(), [&](node_index node) {
            uint64_t in_count = counts[node * 2];
            uint64_t out_count = counts[node * 2 + 1];

            if (in_count >= min_label_in_count && out_count <= max_label_out_count) {
                ++kept_nodes;
            } else {
                mask[node] = false;
            }
        });

        size_t num_filled = 0;
        masked_graph->call_unitigs([&](const std::string &, const auto &path) {
            auto in_mask = [&](node_index node) { return mask[node]; };
            auto begin = std::find_if(path.begin(), path.end(), in_mask);
            if (begin == path.end())
                return;

            auto end = std::find_if(path.rbegin(), path.rend(), in_mask).base();

            size_t num_in_mask = std::count_if(begin, end, in_mask);

            if (static_cast<double>(num_in_mask) / path.size() >= 0.8) {
                std::for_each(path.begin(), path.end(), [&](node_index node) {
                    if (!mask[node]) {
                        mask[node] = true;
                        ++kept_nodes;
                        ++num_filled;
                    }
                });
                return;
            }

            if (static_cast<double>(num_in_mask) / (end - begin) >= 0.8) {
                std::for_each(begin, end, [&](node_index node) {
                    if (!mask[node]) {
                        mask[node] = true;
                        ++kept_nodes;
                        ++num_filled;
                    }
                });
            }

            for (size_t i = 1; i + 1 < path.size(); ++i) {
                node_index node = path[i];
                if (!mask[node]) {
                    if (mask[path[i - 1]] && mask[path[i + 1]]) {
                        mask[node] = true;
                        ++kept_nodes;
                        ++num_filled;
                    }
                }
            }
        });
        */

        masked_graph->set_mask(new bitmap_vector(std::move(mask)));

        // logger->trace("Kept {} out of {} nodes, {} filled", kept_nodes, masked_graph->get_graph().num_nodes(), num_filled);
        masked_graph->likelihood_ratios.resize(counts.size()/2+1);
        std::fill(masked_graph->likelihood_ratios.begin(), masked_graph->likelihood_ratios.end(), 0);

        return masked_graph;
    }

    if (config.label_mask_in_unitig_fraction == 0.0
            && config.label_mask_out_unitig_fraction == 1.0
            && config.label_mask_other_unitig_fraction == 1.0) {
        logger->trace("Filtering by node");
        size_t total_nodes = masked_graph->num_nodes();
        const auto &in_mask = static_cast<const bitmap_vector &>(masked_graph->get_mask()).data();
        sdsl::bit_vector mask = in_mask;

        // TODO: make this part multithreaded
        size_t kept_nodes = 0;
        call_ones(in_mask, [&](node_index node) {
            uint64_t in_count = counts[node * 2];
            uint64_t out_count = counts[node * 2 + 1];

            if (in_count >= min_label_in_count && out_count <= max_label_out_count) {
                ++kept_nodes;
            } else {
                mask[node] = false;
            }
        });

        masked_graph->set_mask(new bitmap_vector(std::move(mask)));

        logger->trace("Kept {} out of {} nodes", kept_nodes, total_nodes);
        masked_graph->likelihood_ratios.resize(counts.size()/2+1);
        std::fill(masked_graph->likelihood_ratios.begin(), masked_graph->likelihood_ratios.end(), 0);

        return masked_graph;
    }

    logger->trace("Filtering by unitig");

    update_masked_graph_by_unitig(*masked_graph,
        [&](const std::string &, const std::vector<node_index> &path) -> Intervals {
            // return a set of intervals to keep in the graph
            size_t in_kmer_count = 0;

            size_t begin = path.size();
            size_t end = 0;
            for (size_t i = 0; i < path.size(); ++i) {
                if (counts[path[i] * 2] >= min_label_in_count) {
                    if (begin == path.size())
                        begin = i;

                    end = std::max(end, i + 1);

                    ++in_kmer_count;
                }
            }

            if (begin >= end)
                return {};

            size_t size = end - begin;
            size_t label_in_cutoff = std::ceil(config.label_mask_in_unitig_fraction * size);
            if (in_kmer_count < label_in_cutoff)
                return {};

            size_t out_kmer_count = 0;
            size_t other_kmer_count = 0;
            size_t label_out_cutoff = std::floor(config.label_mask_out_unitig_fraction * size);
            size_t other_cutoff = std::floor(config.label_mask_other_unitig_fraction * size);


            for (size_t i = begin; i < end; ++i) {
                if (counts[path[i] * 2 + 1] > max_label_out_count
                        && ++out_kmer_count > label_out_cutoff) {
                    return {};
                }

                if (check_other && other_mask[path[i]]
                        && ++other_kmer_count > other_cutoff) {
                    return {};
                }
            }

            return { std::make_pair(begin, end) };
        },
        num_threads
    );

    return masked_graph;
}


/**
 * Helpers
 */

std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                          sdsl::int_vector<> &counts, // Myrthe later: int_vector_buffer and make it a template class. Otherwise the function stays the same because initial mask should only check if a k-mer is present just once in either set.
                          sdsl::bit_vector&& mask,
                          bool add_complement,
                          size_t num_threads) {
    // counts is a double-length vector storing the in-label and out-label
    // counts interleaved
    assert(counts.size() == mask.size() * 2);

    auto masked_graph = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr,
        add_complement
            ? std::make_unique<bitmap_vector>(mask)
            : std::make_unique<bitmap_vector>(std::move(mask)),
            true,
            graph_ptr->get_mode() == DeBruijnGraph::PRIMARY ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
    );

    logger->trace("Constructed masked graph with {} nodes", masked_graph->num_nodes());

    if (add_complement) {
        logger->trace("Adding reverse complements");
        std::mutex vector_backup_mutex;
        std::atomic_thread_fence(std::memory_order_release);
        masked_graph->call_sequences([&](const std::string &seq, const std::vector<node_index> &path) {
            std::string rc_seq = seq;
            std::vector<DeBruijnGraph::node_index> rc_path = path;
            reverse_complement_seq_path(*graph_ptr, rc_seq, rc_path);

            auto it = rc_path.rbegin();
            for (size_t i = 0; i < path.size(); ++i, ++it) {
                if (*it) {
                    uint64_t in_count = atomic_fetch(counts, path[i] * 2, vector_backup_mutex, MO_RELAXED);
                    uint64_t out_count = atomic_fetch(counts, path[i] * 2 + 1, vector_backup_mutex, MO_RELAXED);
                    atomic_fetch_and_add(counts, *it * 2, in_count, vector_backup_mutex, MO_RELAXED);
                    atomic_fetch_and_add(counts, *it * 2 + 1, out_count, vector_backup_mutex, MO_RELAXED);

                    set_bit(mask.data(), *it, in_count, MO_RELAXED);
                }
            }
        }, num_threads, true);

        std::atomic_thread_fence(std::memory_order_acquire);

        masked_graph->set_mask(new bitmap_vector(std::move(mask)));

        logger->trace("Updated masked graph has {} nodes", masked_graph->num_nodes());
    }

    return masked_graph;
}


std::tuple<sdsl::int_vector<>, sdsl::bit_vector, sdsl::bit_vector, std::tuple<size_t, size_t>>
construct_diff_label_count_vector(const ColumnGenerator &generate_columns,
                                  const std::vector<std::string> &files,
                                  const tsl::hopscotch_set<Label> &labels_in,
                                  const tsl::hopscotch_set<Label> &labels_out,
                                  const DifferentialAssemblyConfig &config,
                                  size_t max_index,
                                  size_t num_labels,
                                  size_t num_threads,
                                  bool add_out_labels_to_mask) {
    bool parallel = (num_threads > 1);
    logger->trace("Allocating mask vector");

    size_t width;

    if (config.count_kmers){ // calculate the sum of column widths to derive an upper bound for the required counts vector width.
        assert(labels_in.size() > 0);
        size_t sum_widths_in = 0; size_t sum_widths_out = 0;
        for (size_t i = 0; i < files.size(); ++i) { // loopes over files
            const auto &values_fname = utils::remove_suffix(files[i],
                 annot::ColumnCompressed<>::kExtension) + annot::ColumnCompressed<>::kCountExtension;
            sdsl::int_vector_size_type col_size; uint8_t col_width;
            std::ifstream in_stream(values_fname);
            sdsl::int_vector<>::read_header(col_size, col_width, in_stream);
            auto label_encoder = annot::ColumnCompressed<>::read_label_encoder(files[i]);
            for (size_t c = 0; c < label_encoder.size(); ++c) { // loopes over kmers
                auto label = label_encoder.decode(c);
                if (labels_in.count(label)) sum_widths_in += std::pow(2, col_width);
                if (labels_out.count(label)) sum_widths_out += std::pow(2, col_width);
            }
        }
        width = std::min((size_t) 32, (size_t) std::ceil(std::log(std::max(sum_widths_in, sum_widths_out)))); // TODO, check if this is really correct. this should be 8.5 -> 9 rather then 7.
    } else{
        width = (sdsl::bits::hi(num_labels) + 1) * (1 + add_out_labels_to_mask);
    }
    assert(width);

    sdsl::bit_vector indicator(max_index + 1, false);

    logger->trace("Allocating count vector");
    // the in and out counts are stored interleaved
    sdsl::int_vector<> counts = aligned_int_vector(indicator.size() * 2, 0, width, 16); //Myrthe later: int_vector_buffer. In that case the function also to be a template too, or a template based argument, if vector_type== buffered_int_vector, then it should be a buffered (// if constexpr (std::is_same_v<vector_type, Bitmap>) { } //sdsl::int_vector_buffer<> counts = aligned_int_vector(

    logger->trace("done");

    sdsl::bit_vector other_labels(num_labels, false);
    std::mutex vector_backup_mutex;

    logger->trace("Populating count vector");
    std::atomic_thread_fence(std::memory_order_release); //Myrthe ANSWER: outer loop is over the columns, inner loop over the kmers (i indicating the rows)
            // Myrthe later: For the buffered_int_vector we would have to use operators [ ] [ ], but it might be problematic if it runs out of memory and push_back would need to be used instead.
            // Myrthe later: The second problem with buffered_int_vector might be that it is expensive, because we do random access on it. Instead a solution might be to use a piority queue of L columns, storing the smallest element as a triple (column , index, value), such that always k-mers in row x are processed before row x + 1.

    size_t in_total_kmers = 0; size_t out_total_kmers = 0;

    generate_columns([&](uint64_t offset, const Label &label, const ValueGenerator &value_generator) {
        uint8_t col_indicator = static_cast<bool>(labels_in.count(label));
        if (labels_out.count(label))
            col_indicator |= 2;

        if (!col_indicator)
            set_bit(other_labels.data(), offset, parallel, MO_RELAXED);

        ValueCallback add_in = [&](Column i, uint64_t value) {
            i = AnnotatedDBG::anno_to_graph_index(i); // indices in the annotation matrix have an offset of 1 compared to those in the graph. Since the mask should be compatible with the graph, we have to convert the index.
            assert(i != DeBruijnGraph::npos);
            set_bit(indicator.data(), i, parallel, MO_RELAXED); // set the corresponding bit to true in the mask, with the 'indicator' representing the mask.
            atomic_fetch_and_add(counts, i * 2, value, vector_backup_mutex, MO_RELAXED); // TODO Myrthe later: make sure that values do not overflow in neighbouring cells.
            in_total_kmers += value;
        };

        ValueCallback add_out = [&](Column i, uint64_t value) {
            i = AnnotatedDBG::anno_to_graph_index(i);
            assert(i != DeBruijnGraph::npos);
            if (add_out_labels_to_mask || config.test_type == "likelihoodratio")
                set_bit(indicator.data(), i, parallel, MO_RELAXED);
            atomic_fetch_and_add(counts, i * 2 + 1, value, vector_backup_mutex, MO_RELAXED);
            out_total_kmers += value;
        };

        ValueCallback add_both = [&](Column i, uint64_t value) { // theoretically a label should not be contained in both the in- and out-labels, but this is just in case.
            i = AnnotatedDBG::anno_to_graph_index(i);
            assert(i != DeBruijnGraph::npos);
            set_bit(indicator.data(), i, parallel, MO_RELAXED);
            atomic_fetch_and_add(counts, i * 2, value, vector_backup_mutex, MO_RELAXED);
            atomic_fetch_and_add(counts, i * 2 + 1, value, vector_backup_mutex, MO_RELAXED);
            out_total_kmers += value;
            in_total_kmers += value;
        };

        switch (col_indicator) {
            case 1: { value_generator(add_in); } break;
            case 2: { value_generator(add_out); } break;
            case 3: { value_generator(add_both); } break;
            default: {}
        }
    });

    size_t num_tested_kmers = 0;
    for(size_t i = 0; i < counts.size(); i+=2){
        if (counts[i] > counts[i+1])
        num_tested_kmers += 1;
    };
    logger->trace("num_tested_kmers: {}", num_tested_kmers);

    std::atomic_thread_fence(std::memory_order_acquire);
    logger->trace("done");
    assert( std::accumulate(counts.begin(), counts.end(), 0) > 0); // assert that the sum of the count vector is greater than 0

    return std::make_tuple(std::move(counts), std::move(indicator), std::move(other_labels),
                           std::make_tuple(std::move(in_total_kmers), std::move(out_total_kmers)));
}



void kmer_distribution_table(const ColumnGenerator &generate_columns,
                                  const tsl::hopscotch_set<Label> &labels_in,
                                  const tsl::hopscotch_set<Label> &labels_out,
                                  size_t max_index,
                                  MaskedDeBruijnGraph &masked_graph, // std::shared_ptr<const DeBruijnGraph> graph_ptr,
                             const DifferentialAssemblyConfig &config,
                             size_t num_threads) {

    logger->trace("Allocating count vector");
    std::vector<std::vector<int>> count_matrix(max_index*2 + 2); // outer array 'max_index + 1'*2, with arrays of initially zero lenght or of max size in labels/out_labels.

    logger->trace("Populating count matrix");
    generate_columns([&](uint64_t offset, const Label &label, const ValueGenerator &value_generator) {
        std::ignore = offset;// to overcome the error that 'offset' is set but not used.
        uint8_t col_indicator = static_cast<bool>(labels_in.count(label));
        std::cout << "col_indicator: " << std::to_string(col_indicator) << std::endl;
        if (labels_out.count(label))
            col_indicator |= 2;

        ValueCallback add_in = [&](Column i, uint64_t value) {
            std::cout << "i: " << std::to_string(i) << std::endl;
            std::cout << "value: " << std::to_string(value) << std::endl;
            count_matrix[i * 2].push_back(value);
        };

        ValueCallback add_out = [&](Column i, uint64_t value) {
            count_matrix[i * 2 + 1].push_back(value);
        };

        ValueCallback add_both = [&](Column i, uint64_t value) { // theoretically a label should not be contained in both the in- and out-labels, but this is just in case.
            count_matrix[i * 2].push_back(value);
            count_matrix[i * 2 + 1].push_back(value);

        };

        switch (col_indicator) {
            case 1: { value_generator(add_in); } break;
            case 2: { value_generator(add_out); } break;
            case 3: { value_generator(add_both); } break;
            default: {}
        }

        // Calculate the product of the poisson probabilities and compare with negative binom
    });

    if (config.test_by_unitig){         // unitig mode.
        std::vector<std::vector<int>> unitig_count_matrix;
        update_masked_graph_by_unitig(masked_graph, // intitial mask
            [&](const std::string &, const std::vector<node_index> &path)
                    -> Intervals {
                // return a set of intervals to keep in the graph
                std::vector<int> in_kmer_unitig(labels_in.size(), 0);
                std::vector<int> out_kmer_unitig(labels_out.size(), 0);
                for (size_t i = 0; i < path.size(); ++i) { // first strategy is summing counts per unitig. // The update_mask.. function takes a path. Count vector per path. QUESTION: where is the count vector created, such that unitig k-mers can be accessed sequentially?
                    assert(count_matrix.size() >= path[i] * 2);
                    if (count_matrix[path[i] * 2].size() != in_kmer_unitig.size() or count_matrix[path[i] * 2 + 1].size() != out_kmer_unitig.size()){
                        std::cout << "sizes" << std::endl;
                        std::cout << std::to_string(count_matrix[path[i] * 2].size()) << std::endl << std::flush;
                        std::cout << std::to_string(in_kmer_unitig.size()) << std::endl << std::flush;
                        std::cout << std::to_string(count_matrix[path[i] * 2 + 1].size()) << std::endl << std::flush;
                        std::cout << std::to_string(out_kmer_unitig.size()) << std::endl << std::flush;
                        std::cout << std::to_string(path[i] * 2) << std::endl << std::flush; // TODO: do I have to convert the index? It was not done in the original code
                    }else{
                        assert(count_matrix[path[i] * 2].size() == in_kmer_unitig.size());
                        assert(count_matrix[path[i] * 2 + 1].size() == out_kmer_unitig.size());
                        std::transform(in_kmer_unitig.begin(), in_kmer_unitig.end(), count_matrix[path[i] * 2].begin(), in_kmer_unitig.begin(), std::plus<int>()); // sum the values of count_matrix[x] to in_kmer_count.
                        std::transform(out_kmer_unitig.begin(), out_kmer_unitig.end(), count_matrix[path[i] * 2 + 1].begin(), out_kmer_unitig.begin(), std::plus<int>());
                    }
                }
                unitig_count_matrix.push_back(in_kmer_unitig);
                unitig_count_matrix.push_back(out_kmer_unitig); // add to matrix
                return { };
                },
            num_threads);
        count_matrix = unitig_count_matrix;
    }

    int matrix_sum = 0;
    for (size_t i = 0; i < count_matrix.size(); i+=2){
        matrix_sum += std::accumulate(count_matrix[i].begin(), count_matrix[i].end(), 0);
    }
    assert(matrix_sum > 0);

    logger->trace("Writing k-mer count matrix to table");

    std::ofstream out_stream("kmer_dist_table.py");
    out_stream << "count_matrix_in = [" ;
    for (size_t i = 0; i < std::min(count_matrix.size(), (size_t) 50000000); i+=2){
        if (count_matrix[i+1].size() > 3 or count_matrix[i].size() > 3){ // labels_out.size()/3
                std::vector<int> zeros_vector(labels_in.size()-count_matrix[i].size(), 0);
                count_matrix[i].insert(count_matrix[i].end(), zeros_vector.begin(), zeros_vector.end());
                out_stream << "[" + outstring(count_matrix[i]) << "]," <<std::endl;
            }
    }
    out_stream << "]" << std::endl;
    out_stream << "count_matrix_out = [" ;
    for (size_t i = 1; i < std::min(count_matrix.size(), (size_t) 50000000); i+=2){
        if (count_matrix[i].size() > 3 or count_matrix[i-1].size() > 3){ // labels_out.size()/3
            std::vector<int> zeros_vector(labels_out.size()-count_matrix[i].size(), 0);
            count_matrix[i].insert(count_matrix[i].end(), zeros_vector.begin(), zeros_vector.end());
            out_stream << "[" + outstring(count_matrix[i]) << "]," <<std::endl;        }
    }
    out_stream << "]"<< std::endl;
    out_stream.close();
}

void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const GetKeptIntervals &get_kept_intervals,
                                   size_t num_threads) { // traverses unitigs and takes a callback that does the filtering.
    std::atomic<uint64_t> kept_unitigs(0);
    std::atomic<uint64_t> total_unitigs(0);
    std::atomic<uint64_t> num_kept_nodes(0);
    bool parallel = num_threads > 1;

    sdsl::bit_vector mask = static_cast<const bitmap_vector&>(masked_graph.get_mask()).data();

    std::atomic_thread_fence(std::memory_order_release);

    masked_graph.call_unitigs([&](const std::string &unitig,
                                  const std::vector<node_index> &path) {
        total_unitigs.fetch_add(1, MO_RELAXED);

        size_t last = 0;
        for (const auto &pair : get_kept_intervals(unitig, path)) {
            const auto &[begin, end] = pair;
            kept_unitigs.fetch_add(1, MO_RELAXED);
            num_kept_nodes.fetch_add(end - begin, MO_RELAXED);
            for ( ; last < begin; ++last) {
                unset_bit(mask.data(), path[last], parallel, MO_RELAXED);
            }
            last = end;
        }

        for ( ; last < path.size(); ++last) {
            unset_bit(mask.data(), path[last], parallel, MO_RELAXED);
        }

    }, num_threads);
    std::atomic_thread_fence(std::memory_order_acquire);

    masked_graph.set_mask(new bitmap_vector(std::move(mask)));

    logger->trace("Kept {} out of {} unitigs with average length {} and a total of {} nodes",
                  kept_unitigs, total_unitigs,
                  static_cast<double>(num_kept_nodes + kept_unitigs * (masked_graph.get_k() - 1))
                      / kept_unitigs,
                  num_kept_nodes);
}

//%TODO A post-filtering step is done to filter out short contigs, percentage/cutoff , that is set to by default.  Short unitigs are removed afterwards
//        % TODO For each contig the confidence score is output


} // namespace graph
} // namespace mtg


