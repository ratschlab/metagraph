#include "annotated_graph_algorithm.hpp"

#include <typeinfo>
#include <mutex>
#include <variant>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/distributions/cauchy.hpp>

#include "common/logger.hpp"
#include "common/vectors/bitmap.hpp"
#include "common/vector_map.hpp"
#include "common/vectors/transpose.hpp"
#include "common/hashers/hash.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/graph_cleaning.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/int_matrix/base/int_matrix.hpp"

namespace mtg {
namespace graph {

using mtg::common::logger;

typedef AnnotatedDBG::node_index node_index;
typedef AnnotatedDBG::Annotator Annotator;
typedef AnnotatedDBG::Annotator::Label Label;
using Column = annot::matrix::BinaryMatrix::Column;
constexpr std::memory_order MO_RELAXED = std::memory_order_relaxed;

double CDF_CUTOFF = 0.95;
double HIST_CUTOFF = 1.00;
uint64_t N_BUCKETS_FOR_ESTIMATION = 3;


template<class To, class From>
std::enable_if_t<
    sizeof(To) == sizeof(From) &&
    std::is_trivially_copyable_v<From> &&
    std::is_trivially_copyable_v<To>,
    To>
// constexpr support needs compiler magic
bit_cast(const From& src) noexcept
{
    static_assert(std::is_trivially_constructible_v<To>,
        "This implementation additionally requires "
        "destination type to be trivially constructible");

    static_assert(sizeof(To) == sizeof(From));
    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

template <typename value_type, class PValStorage, typename HistGetter, typename Generator>
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, PValStorage, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                         const HistGetter &get_hist_map,
                         const Generator &generate_rows,
                         const std::vector<bool> &groups,
                         const DifferentialAssemblyConfig &config,
                         size_t num_threads = 1,
                         std::filesystem::path tmp_dir = "",
                         size_t num_parallel_files = std::numeric_limits<size_t>::max(),
                         const std::function<void()> &deallocate = []() {},
                         uint8_t max_width = 64) {
    // num_parallel_files = std::min(num_parallel_files, num_threads);
    num_parallel_files = get_num_threads();
    num_threads = get_num_threads();
    bool is_primary = graph_ptr->get_mode() == DeBruijnGraph::PRIMARY;
    bool parallel = num_parallel_files > 1;

    size_t num_labels_out = std::count(groups.begin(), groups.end(), true);
    size_t num_labels_in = groups.size() - num_labels_out;

    common::logger->trace("Graph mode: {}", is_primary ? "PRIMARY" : "other");

    common::logger->trace("Computing histogram from uncleaned counts");
    auto hists_map = get_hist_map(std::vector<size_t>(groups.size(), 1), nullptr);

    sdsl::bit_vector kept(AnnotatedDBG::graph_to_anno_index(graph_ptr->max_index() + 1), true);

    std::vector<size_t> min_counts(groups.size(), config.min_count);
    std::vector<uint64_t> check_cutoff(groups.size(), std::numeric_limits<uint64_t>::max());

    std::mutex agg_mu;
    using PairContainer = std::vector<std::pair<uint64_t, value_type>>;

    if (config.clean) {
        common::logger->trace("Cleaning count columns");

        #pragma omp parallel for num_threads(num_parallel_files)
        for (size_t j = 0; j < groups.size(); ++j) {
            // set cutoff for lower end of distribution
            auto [mean_est, nzeros_est] = estimate_ztp_mean(
                [&](const auto &callback) {
                    for (const auto &[k, c] : hists_map[j]) {
                        if (k > 0 && k <= N_BUCKETS_FOR_ESTIMATION)
                            callback(k, c);
                    }
                },
                max_width,
                0,
                N_BUCKETS_FOR_ESTIMATION
            );

            min_counts[j] = boost::math::quantile(
                boost::math::poisson(mean_est), exp(-mean_est) - CDF_CUTOFF * expm1(-mean_est)
            );
        }
    }

    auto correct_pvals = [&config](const std::vector<std::pair<double, size_t>> &m,
                                   size_t nelem,
                                   const std::vector<int64_t> &m_sums = {}) {
        size_t k = 0;
        size_t n_cutoff = 1;
        size_t acc = std::accumulate(m.begin(), m.end(), size_t(0),
                                     [](size_t sum, const auto &a) { return sum + a.second; });
        size_t total = acc;
        common::logger->trace("Calculating cutoffs for {}/{} tests", acc, nelem);
        size_t last_k = 0;
        for (size_t n = 0; n < m.size(); ++n) {
            auto [pval_min, s] = m[n];
            size_t cur_k = std::numeric_limits<size_t>::max();
            if (pval_min >= config.family_wise_error_rate) {
                cur_k = 0;
            } else if (pval_min > 0) {
                double lkd = log2(config.family_wise_error_rate) - log2(pval_min);
                if (lkd <= 64)
                    cur_k = pow(2.0, lkd);

                if (cur_k == 0) {
                    common::logger->error("k: {}\tlog k: {}\tpval_min: {}", cur_k, lkd, pval_min);
                    throw std::runtime_error("Min failed");
                }
            }

            if (cur_k != last_k && acc <= cur_k) {
                k = std::max(acc, last_k + 1);
                if (m_sums.empty()) {
                    n_cutoff = n;
                } else {
                    n_cutoff = *std::min_element(m_sums.begin() + n, m_sums.end());
                }
                break;
            }

            last_k = cur_k;
            acc -= s;
        }

        if (k == 0)
            common::logger->warn("No significant p-values achievable after correction");

        common::logger->trace("Kept {} / {}", acc, total);
        return std::make_pair(k, n_cutoff);
    };

    common::logger->trace("Updating histogram and marking discarded k-mers");

    if (groups.size())
        hists_map = get_hist_map(min_counts, &kept);

    size_t nelem = sdsl::util::cnt_one_bits(kept);

    std::vector<std::vector<std::pair<uint64_t, size_t>>> hists(groups.size());
    for (size_t j = 0; j < hists.size(); ++j) {
        hists[j] = const_cast<std::vector<std::pair<uint64_t, size_t>>&&>(hists_map[j].values_container());
        hists_map[j].clear();
        std::sort(hists[j].begin(), hists[j].end(), utils::LessFirst());
    }
    hists_map.clear();

    for (size_t j = 0; j < groups.size(); ++j) {
        size_t total_c = 0;
        for (const auto &[k, c] : hists[j]) {
            total_c += c;
        }

        if (total_c != nelem) {
            common::logger->error("FAIL {}: {} != {}", j, total_c, nelem);
            throw std::runtime_error("GGG");
        }
    }

    common::logger->trace("Computing aggregate statistics");
    int64_t in_kmers = 0;
    int64_t out_kmers = 0;
    std::vector<uint64_t> max_obs_vals(groups.size());
    std::vector<size_t> n_kmers(groups.size());
    std::vector<uint64_t> sums(groups.size());
    std::vector<uint64_t> sums_of_squares(groups.size());

    for (size_t j = 0; j < groups.size(); ++j) {
        for (const auto &[k, c] : hists[j]) {
            sums[j] += k * c;
            sums_of_squares[j] += k * k * c;
            if (k != 0)
                n_kmers[j] += c;

            max_obs_vals[j] = std::max(max_obs_vals[j], k);
        }

        if (groups[j]) {
            out_kmers += sums[j];
        } else {
            in_kmers += sums[j];
        }

        common::logger->trace("{}: n_unique: {}\tsum: {}\tmax_obs: {}\tmin_cutoff: {}\tmax_cutof: {}",
                              j, n_kmers[j], sums[j], max_obs_vals[j], min_counts[j], check_cutoff[j]);
    }

    int64_t total_kmers = in_kmers + out_kmers;

    common::logger->trace("Number of kept unique k-mers: {}\tNumber of kept k-mers: {}",
                          nelem, total_kmers);

    common::logger->trace("Allocating p-value storage");
    auto nullpval = bit_cast<uint64_t>(double(1.1));

    std::vector<PValStorage> pvals_buckets;
    std::vector<std::unique_ptr<utils::TempFile>> tmp_buckets;
    std::vector<PValStorage> pvals_min_buckets;
    std::vector<std::unique_ptr<utils::TempFile>> tmp_min_buckets;
    constexpr bool preallocated = std::is_same_v<PValStorage, std::vector<uint64_t>>;

    if constexpr(std::is_same_v<PValStorage, std::vector<uint64_t>>) {
        auto &pvals = pvals_buckets.emplace_back();
        pvals.resize(graph_ptr->max_index() + 1, nullpval);
    }

    if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
        for (size_t i = 0; i < get_num_threads() + 1; ++i) {
            auto &tmp_file = tmp_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
            pvals_buckets.emplace_back(tmp_file->name(), std::ios::out);
        }
        pvals_buckets[0].push_back(nullpval);
    }

    if (config.test_type == "gnb_exact" && config.test_by_unitig) {
        if constexpr(std::is_same_v<PValStorage, std::vector<uint64_t>>) {
            auto &pvals_min = pvals_min_buckets.emplace_back();
            pvals_min.resize(graph_ptr->max_index() + 1, nullpval);
        }

        if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
            for (size_t i = 0; i < get_num_threads() + 1; ++i) {
                auto &tmp_file = tmp_min_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
                pvals_min_buckets.emplace_back(tmp_file->name(), std::ios::out);
            }
            pvals_min_buckets[0].push_back(nullpval);
        }
    }

    std::vector<PValStorage> eff_size_buckets;
    std::vector<std::unique_ptr<utils::TempFile>> tmp_eff_buckets;
    std::vector<PValStorage> sum_buckets;
    std::vector<std::unique_ptr<utils::TempFile>> tmp_sum_buckets;
    if (config.test_by_unitig) {
        common::logger->trace("Allocating effect size and row sum storage");
        if constexpr(std::is_same_v<PValStorage, std::vector<uint64_t>>) {
            auto &eff_size = eff_size_buckets.emplace_back();
            eff_size.resize(graph_ptr->max_index() + 1);

            auto &sum = sum_buckets.emplace_back();
            sum.resize(graph_ptr->max_index() + 1);
        }

        if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
            for (size_t i = 0; i < get_num_threads() + 1; ++i) {
                auto &tmp_file_eff = tmp_eff_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
                eff_size_buckets.emplace_back(tmp_file_eff->name(), std::ios::out);

                auto &tmp_file_sum = tmp_sum_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
                sum_buckets.emplace_back(tmp_file_sum->name(), std::ios::out);
            }
            eff_size_buckets[0].push_back(0);
            sum_buckets[0].push_back(0);
        }
    }

    sdsl::bit_vector indicator_in;
    sdsl::bit_vector indicator_out;

    if (!config.test_by_unitig) {
        common::logger->trace("Allocating k-mer bitmasks");
        indicator_in = sdsl::bit_vector(graph_ptr->max_index() + 1, false);
        indicator_out = sdsl::bit_vector(graph_ptr->max_index() + 1, false);
    }

    std::function<double(int64_t, int64_t, const PairContainer&)> compute_pval;
    std::function<double(int64_t, const PairContainer&)> compute_min_pval;

    std::vector<VectorMap<uint64_t, std::pair<size_t, uint64_t>>> count_maps;

    common::logger->trace("Test: {}", config.test_type);
    if (config.test_type == "poisson_exact") {
        auto get_deviance = [&](double y, double mu) {
            y += 1e-8;
            mu += 1e-8;
            return 2 * (y * log(y/mu) - y + mu);
        };

        compute_min_pval = [&,p=static_cast<double>(in_kmers)/total_kmers,
                              mu1=static_cast<double>(in_kmers)/nelem,
                              mu2=static_cast<double>(out_kmers)/nelem,
                              get_deviance](int64_t n, const PairContainer&) {
            if (n == 0)
                return 1.0;

            auto bdist = boost::math::binomial(n, p);

            std::vector<double> devs;
            devs.reserve(n + 1);
            for (int64_t s = 0; s <= n; ++s) {
                devs.emplace_back(get_deviance(s, mu1) + get_deviance(n - s, mu2));
            }
            double max_d = *std::max_element(devs.begin(), devs.end());
            double pval = 0.0;

            int64_t s = 0;
            for ( ; s <= n; ++s) {
                if (devs[s] != max_d)
                    break;
            }
            if (s > 0)
                pval += boost::math::cdf(bdist, s - 1);

            int64_t sp = n;
            for ( ; sp >= s; --sp) {
                if (devs[sp] != max_d)
                    break;
            }

            if (sp < n)
                pval += boost::math::cdf(boost::math::complement(bdist, sp));

            return pval;
        };

        compute_pval = [&,p=static_cast<double>(in_kmers)/total_kmers,
                          mu1=static_cast<double>(in_kmers)/nelem,
                          mu2=static_cast<double>(out_kmers)/nelem,
                          get_deviance](int64_t in_sum, int64_t out_sum, const auto &row) {
            if (row.empty())
                return 1.0;

            int64_t n = in_sum + out_sum;
            auto bdist = boost::math::binomial(n, p);

            std::vector<double> devs;
            devs.reserve(n + 1);
            for (int64_t s = 0; s <= n; ++s) {
                devs.emplace_back(get_deviance(s, mu1) + get_deviance(n - s, mu2));
            }
            double pval = 0.0;

            int64_t s = 0;
            for ( ; s <= n; ++s) {
                if (devs[s] < devs[in_sum])
                    break;
            }
            if (s > 0)
                pval += boost::math::cdf(bdist, s - 1);

            int64_t sp = n;
            for ( ; sp >= s; --sp) {
                if (devs[sp] < devs[in_sum])
                    break;
            }

            if (sp < n)
                pval += boost::math::cdf(boost::math::complement(bdist, sp));

            return pval;
        };
    } else if (config.test_type == "fisher") {
        double lbase_shared = lgamma(in_kmers + 1) + lgamma(out_kmers + 1) - lgamma(total_kmers + 1);

        compute_min_pval = [&,lbase_shared](int64_t m1, const PairContainer&) {
            if (m1 == 0)
                return 1.0;

            int64_t m2 = total_kmers - m1;

            auto get_pval = [&](int64_t a) {
                double lbase = lbase_shared + lgamma(m1 + 1) + lgamma(m2 + 1);

                int64_t b = in_kmers - a;
                int64_t c = m1 - a;
                int64_t d = m2 - b;
                assert(d == out_kmers - c);

                return out_kmers < c
                    ? 1.0
                    : exp(lbase - lgamma(a + 1) - lgamma(b + 1) - lgamma(c + 1) - lgamma(d + 1));
            };

            return std::min(get_pval(0), get_pval(std::min(in_kmers, m1)));
        };

        compute_pval = [&,lbase_shared](int64_t in_sum, int64_t out_sum, const auto &row) {
            if (row.empty())
                return 1.0;

            int64_t m1 = in_sum + out_sum;
            int64_t m2 = total_kmers - m1;

            double lbase = lbase_shared + lgamma(m1 + 1) + lgamma(m2 + 1);

            int64_t a = in_sum;
            int64_t b = in_kmers - in_sum;
            int64_t c = out_sum;
            int64_t d = out_kmers - out_sum;

            return exp(lbase - lgamma(a + 1) - lgamma(b + 1) - lgamma(c + 1) - lgamma(d + 1));
        };
    } else if (config.test_type == "nbinom_exact" || config.test_type == "gnb_exact") {
        common::logger->trace("Fitting per-sample negative binomial distributions");
        auto get_rp = [&](const auto &generate,
                          const auto &add_dl,
                          const auto &add_ddl) {
            double mu = 0;
            double var = 0;
            size_t total = 0;
            generate([&](auto k, auto c) {
                mu += k * c;
                var += k * k * c;
                total += c;
            });
            mu /= total;
            double mu2 = mu * mu;
            var = (var - mu2 * total) / (total - 1);

            double r_guess = mu * mu / (var - mu);

            if (mu >= var) {
                common::logger->warn("Fit failed, falling back to Poisson: mu: {} >= var: {}", mu, var);
                return std::make_pair(0.0, 1.0);
            }

            auto get_dl = [&](double r) {
                double dl = (log(r) - log(r + mu) - boost::math::digamma(r)) * total;
                generate([&](auto k, auto c) {
                    dl += boost::math::digamma(k + r) * c;
                });
                return dl + add_dl(r);
            };

            double r_min = r_guess;
            double r_max = r_guess;

            while (true) {
                double dl_min = get_dl(r_min);
                double dl_max = get_dl(r_max);

                if (dl_min == 0) {
                    r_max = r_min;
                    break;
                }

                if (dl_max == 0) {
                    r_min = r_max;
                    break;
                }

                if (dl_min < 0)
                    r_min /= 2;

                if (dl_max > 0)
                    r_max *= 2;

                if (dl_min > 0 && dl_max < 0)
                    break;
            }

            std::tie(r_min, r_max) = boost::math::tools::bisect(
                get_dl, r_min, r_max, boost::math::tools::eps_tolerance<double>(1e-5)
            );

            double r = r_max;
            double p = r / (r + mu);
            return std::make_pair(r, p);
        };

        std::vector<std::pair<double, double>> nb_params(groups.size());
        std::vector<std::vector<std::pair<uint64_t, size_t>>::const_iterator> hist_its(groups.size());
        size_t matrix_size = 0;

        #pragma omp parallel for num_threads(num_parallel_files) reduction(+:matrix_size)
        for (size_t j = 0; j < groups.size(); ++j) {
            const auto &hist = hists[j];
            size_t total = 0;
            size_t total_cutoff = n_kmers[j] * HIST_CUTOFF;
            auto it = hist.begin();
            for ( ; it != hist.end() && total < total_cutoff; ++it) {
                if (it->first != 0)
                    total += it->second;
            }
            hist_its[j] = it;
            nb_params[j] = get_rp([&](const auto &callback) {
                std::for_each(hist.begin(), it, [&](const auto &a) {
                    const auto &[k, c] = a;
                    callback(k, c);
                });
            }, [](double) { return 0; }, [](double) { return 0; });
            const auto &[r, p] = nb_params[j];
            common::logger->trace("{}: size: {}\tmax_val: {}\tmu: {}\tvar: {}\tmle: r: {}\tp: {}",
                                  j, sums[j], (it - 1)->first, r * (1-p)/p, r*(1-p)/p/p, r, p);

            total += (hist.size() && hist[0].first == 0 ? hist[0].second : 0);

            #pragma omp atomic
            matrix_size += total;
        }

        common::logger->trace("Scaling distributions");

        double target_sum = 0.0;
        for (size_t j = 0; j < groups.size(); ++j) {
            target_sum += log(static_cast<double>(sums[j]));
        }
        target_sum = exp(target_sum / groups.size());

        common::logger->trace("Finding common p and r.");

        auto [r_map, target_p] = get_rp([&](const auto &callback) {
            for (size_t j = 0; j < groups.size(); ++j) {
                double f = target_sum / sums[j];
                std::for_each(hists[j].cbegin(), hist_its[j], [&](const auto &a) {
                    const auto &[k, c] = a;
                    callback(f * k, c);
                });
            }
        }, [](double) { return 0; }, [](double) { return 0; });

        double fit_mu = r_map * (1.0 - target_p) / target_p;
        double fit_var = fit_mu / target_p;
        common::logger->trace("Common params: r: {}\tp: {}\tmu: {}\tvar: {}",
                              r_map, target_p, fit_mu, fit_var);

        double r_in = r_map * num_labels_in;
        double r_out = r_map * num_labels_out;
        common::logger->trace("Fits: in: {}\tout: {}\tp: {}", r_in, r_out, target_p);

        if (fit_var / fit_mu - 1.0 < 1e-5)
            common::logger->warn("Fit parameters are close to a Poisson distribution");

        count_maps.resize(groups.size());

        common::logger->trace("Computing quantile maps");
        #pragma omp parallel for num_threads(num_parallel_files)
        for (size_t j = 0; j < groups.size(); ++j) {
            const auto &[r, p] = nb_params[j];

            boost::math::negative_binomial nb_out(r_map, target_p);
            double scale = target_sum / sums[j];

            auto map_value = [&](uint64_t k, double scale, const auto &old_dist, const auto &new_dist) {
                double cdf = boost::math::cdf(old_dist, k);

                if (cdf < 1.0)
                    return boost::math::quantile(new_dist, cdf);

                double ccdf = boost::math::cdf(boost::math::complement(old_dist, k));

                if (ccdf > 0.0)
                    return boost::math::quantile(boost::math::complement(new_dist, ccdf));

                return ceil(scale * k);
            };

            auto compute_map = [&](const auto &nb) {
                // ensure that 0 maps to 0
                count_maps[j][0] = std::make_pair(0, 0);

                const auto &hist = hists[j];
                for (const auto &[k, c] : hist) {
                    if (!count_maps[j].count(k)) {
                        uint64_t new_k = map_value(k, scale, nb, nb_out);
                        count_maps[j][k] = std::make_pair(new_k, c);
                    } else {
                        count_maps[j][k].second += c;
                    }
                }
            };

            if (p < 1.0) {
                compute_map(boost::math::negative_binomial(r, p));
            } else {
                compute_map(boost::math::poisson(static_cast<double>(sums[j]) * scale / nelem));
            }
        }

        common::logger->trace("Unscaled Totals: in: {}\tout: {}", in_kmers, out_kmers);
        in_kmers = 0;
        out_kmers = 0;
        for (size_t j = 0; j < groups.size(); ++j) {
            for (const auto &[k, m] : count_maps[j]) {
                const auto &[v, c] = m;

                if (groups[j]) {
                    out_kmers += v * c;
                } else {
                    in_kmers += v * c;
                }
            }
        }

        total_kmers = in_kmers + out_kmers;

        common::logger->trace("  Scaled Totals: in: {}\tout: {}", in_kmers, out_kmers);

        auto compute_min_pval_r = [](double r, double r_in, double r_out, double lscaling_base, int64_t n) {
            double half_div0 = -r_in * log2(r_in) - (n + r_out) * log2(n + r_out);
            double half_divn = -r_out * log2(r_out) - (n + r_in) * log2(n + r_in);

            double pval = 0.0;
            double lscaling = lscaling_base - lgamma(r + n) / log(2);
            if (half_div0 >= half_divn)
                pval += exp2(lscaling + (lgamma(r_out + n) - lgamma(r_out)) / log(2));

            if (half_divn >= half_div0)
                pval += exp2(lscaling + (lgamma(r_in + n) - lgamma(r_in)) / log(2));

            return std::min(1.0, pval);
        };

        auto compute_pval_r = [](double r, double r_in, double r_out, double lscaling_base, int64_t in_sum, int64_t out_sum, const PairContainer &row) {
            int64_t n = in_sum + out_sum;
            double argmin_d = r_in / r * n;
            if (in_sum == argmin_d)
                return 1.0;

            double lscaling = lscaling_base - lgamma(r + n) / log(2);

            auto get_pval = [&](const auto &get_stat) {
                double base_stat = get_stat(in_sum, out_sum,
                                            r_in != r_out && in_sum > 0 ? log2(in_sum) : 0.0,
                                            r_in != r_out && out_sum > 0 ? log2(out_sum) : 0.0);
                double pval = 0.0;
                int64_t s = 0;
                int64_t t = n;
                double ls = 0;
                double lt = log2(t);
                if (get_stat(s, t, ls, lt) >= base_stat) {
                    double base = lscaling + (lgamma(n + r_out) - lgamma(r_out)) / log(2);
                    pval += exp2(base);
                    for (++s,--t; s <= n; ++s,--t) {
                        ls = log2(s);
                        lt = t > 0 ? log2(t) : 0.0;
                        if (get_stat(s, t, ls, lt) < base_stat)
                            break;

                        base += log2(t + 1) + log2(s - 1 + r_in) - ls - log2(t + r_out);
                        pval += exp2(base);
                    }
                }

                int64_t sp = n;
                t = 0;
                ls = log2(sp);
                lt = 0;
                if (get_stat(sp, t, ls, lt) >= base_stat) {
                    double base = lscaling + (lgamma(n + r_in) - lgamma(r_in)) / log(2);
                    pval += exp2(base);
                    for (--sp,++t; sp >= s; --sp,++t) {
                        ls = sp > 0 ? log2(sp) : 0.0;
                        lt = log2(t);
                        if (get_stat(sp, t, ls, lt) < base_stat)
                            break;

                        base += log2(sp + 1) + log2(t - 1 + r_out) - lt - log2(sp + r_in);
                        pval += exp2(base);
                    }
                }

                return std::min(1.0, pval);
            };

            if (r_in == r_out) {
                return get_pval([&](int64_t s, int64_t, double, double) { return abs(argmin_d - s); });
            } else {
                return get_pval([&](int64_t s, int64_t t, double ls, double lt) {
                    return ls * s - (r_in + s)*log2(r_in + s) + lt * t - (t + r_out)*log2(t + r_out);
                });
            }
        };

        if (config.test_type == "nbinom_exact") {
            double r = r_in + r_out;
            double lscaling_base = lgamma(r) / log(2);
            compute_min_pval = [compute_min_pval_r,lscaling_base,r,r_in,r_out](int64_t n, const PairContainer&) {
                return n > 0 ? compute_min_pval_r(r, r_in, r_out, lscaling_base, n) : 1.0;
            };

            compute_pval = [compute_pval_r,lscaling_base,r,r_in,r_out](int64_t in_sum, int64_t out_sum, const PairContainer &row) {
                return row.size() ? compute_pval_r(r, r_in, r_out, lscaling_base, in_sum, out_sum, row) : 1.0;
            };
        } else {
            common::logger->trace("Enumerating row count distributions");
            std::vector<tsl::hopscotch_map<std::vector<int64_t>, double, utils::VectorHash>> vector_counts;
            {
                std::vector<std::vector<tsl::hopscotch_map<std::vector<int64_t>, double, utils::VectorHash>>> vector_counts_ms(num_threads + 1);

                generate_rows([&](uint64_t row_i, const auto &raw_row, size_t bucket_idx) {
                    int64_t in_sum = 0;
                    int64_t out_sum = 0;
                    bool in_kmer = false;
                    bool out_kmer = false;
                    std::vector<int64_t> vec;
                    if (kept[row_i]) {
                        vec.reserve(raw_row.size());
                        size_t count_in = 0;
                        size_t count_out = 0;
                        for (const auto &[j, raw_c] : raw_row) {
                            uint64_t c = raw_c;
                            if (count_maps.size())
                                c = count_maps[j].find(c)->second.first;

                            if (groups[j]) {
                                out_sum += c;
                                ++count_out;
                            } else {
                                in_sum += c;
                                ++count_in;
                            }

                            vec.emplace_back(c);
                        }

                        if (count_in + count_out >= config.min_recurrence) {
                            double out_stat = static_cast<double>(out_sum) / out_kmers * in_kmers;
                            in_kmer = count_in >= config.min_in_recurrence && count_out <= config.max_out_recurrence && in_sum > (out_kmers > 0 ? out_stat : 0.0);
                            out_kmer = count_out >= config.min_out_recurrence && count_in <= config.max_in_recurrence && in_sum < out_stat;
                        }
                    } else {
                        return;
                    }

                    if (!in_kmer && !out_kmer)
                        return;

                    size_t n = in_sum + out_sum;
                    if (n >= vector_counts_ms[bucket_idx].size())
                        vector_counts_ms[bucket_idx].resize(n + 1);

                    std::sort(vec.begin(), vec.end());
                    ++vector_counts_ms[bucket_idx][n][vec];
                });

                vector_counts = std::move(vector_counts_ms[0]);
                size_t max_size = vector_counts.size();
                for (size_t j = 1; j < vector_counts_ms.size(); ++j) {
                    max_size = std::max(max_size, vector_counts_ms[j].size());
                }
                vector_counts.resize(max_size);

                ProgressBar progress_bar(vector_counts.size() - 1, "Merging row", std::cerr, !common::get_verbose());
                #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
                for (size_t n = vector_counts.size(); n > 0; --n) {
                    for (size_t j = 1; j < vector_counts_ms.size(); ++j) {
                        if (n <= vector_counts_ms[j].size()) {
                            for (const auto &[v, c] : vector_counts_ms[j][n - 1]) {
                                vector_counts[n - 1][v] += c;
                            }
                        }
                    }
                    ++progress_bar;
                }
            }

            common::logger->trace("Fitting LN prior distribution for dispersions");
            double ln_mu = 0.0;
            double ln_var = 0.0;

            for (size_t n = 1; n < vector_counts.size(); ++n) {
                size_t total = 0;
                for (const auto &[v, c] : vector_counts[n]) {
                    total += c;
                }

                double r = target_p / (1.0 - target_p) * n / groups.size();

                #pragma omp atomic
                ln_mu += -log(r) * total;

                #pragma omp atomic
                ln_var += log(r) * log(r) * total;
            }

            size_t total = nelem * groups.size();
            ln_mu /= total;
            ln_var = ln_var / total - ln_mu * ln_mu;

            common::logger->trace("LN fit: mu: {}\tvar: {}\tE[X]: {}\tVar(X): {}\t1.0/E[X]: {}",
                                  ln_mu, ln_var,
                                  exp(ln_mu + ln_var/2), exp(ln_mu*2+ln_var)*(exp(ln_var)-1),
                                  exp(-ln_mu - ln_var/2));

            auto get_rp = [ln_mu,ln_var,gs=groups.size()](const auto &counts) {
                double mu = 0;
                double var = 0;
                for (int64_t c : counts) {
                    mu += c;
                    var += c * c;
                }
                mu /= gs;

                if (ln_var == 0) {
                    double r = exp(-ln_mu);
                    return std::make_pair(r, r / (r + mu));
                }

                double mu2 = mu * mu;
                var = var / gs - mu2;

                double a_guess = var > mu ? (var - mu) / mu2 : 1.0;

                auto get_dl = [&](double a) {
                    double r = 1.0 / a;
                    double p = r / (r + mu);
                    double lp = log(p);
                    double dl = -a - a*(log(a)-ln_mu)/ln_var - gs*lp + boost::math::digamma(r)*counts.size();
                    for (int64_t c : counts) {
                        dl -= boost::math::digamma(r + c);
                    }
                    return dl;
                };

                double a_min = a_guess;
                double a_max = a_guess;
                while (true) {
                    a_guess = (a_min + a_max) / 2;
                    double dl_min = get_dl(a_min);
                    double dl_max = get_dl(a_max);
                    if (dl_min == 0) {
                        a_max = a_min;
                        break;
                    }

                    if (dl_max == 0) {
                        a_min = a_max;
                        break;
                    }

                    if (dl_min < 0)
                        a_min /= 2;

                    if (dl_max > 0)
                        a_max *= 2;

                    if (dl_min > 0 && dl_max < 0)
                        break;
                }

                if (a_min < a_max) {
                    std::tie(a_min, a_max) = boost::math::tools::bisect(
                        get_dl,
                        a_min, a_max, boost::math::tools::eps_tolerance<double>(1e-5)
                    );
                }

                double a = a_min;
                double r = 1.0 / a;
                double p = r / (r + mu);
                return std::make_pair(r, p);
            };

            ProgressBar progress_bar(vector_counts.size() - 1, "Precomputing r's", std::cerr, !common::get_verbose());
            #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
            for (size_t n = 1; n < vector_counts.size(); ++n) {
                if (n > 0 && n % 1000 == 0)
                    progress_bar += 1000;

                for (auto it = vector_counts[n].begin(); it != vector_counts[n].end(); ++it) {
                    it.value() = get_rp(it->first).first;
                }
            }
            progress_bar += (vector_counts.size() - 1) % 1000;

            auto get_r = std::make_shared<std::function<double(int64_t, const PairContainer&)>>(
                [vc=std::move(vector_counts)](int64_t n, const auto &row) {
                    auto it = vc[n].begin();
                    if (vc[n].size() > 1) {
                        std::vector<int64_t> counts;
                        counts.reserve(row.size());
                        for (const auto &[j, c] : row) {
                            counts.emplace_back(c);
                        }
                        std::sort(counts.begin(), counts.end());
                        it = vc[n].find(counts);
                    }
                    return it->second;
                }
            );

            compute_min_pval = [compute_min_pval_r,get_r,num_labels_in,num_labels_out](int64_t n, const PairContainer &row) {
                if (row.empty())
                    return 1.0;

                double r_base = (*get_r)(n, row);
                if (r_base == 0.0)
                    return 1.0;

                double r_in = r_base * num_labels_in;
                double r_out = r_base * num_labels_out;
                double r = r_in + r_out;
                double lscaling_base = lgamma(r) / log(2);

                return compute_min_pval_r(r, r_in, r_out, lscaling_base, n);
            };

            compute_pval = [compute_pval_r,get_r,num_labels_in,num_labels_out](int64_t in_sum, int64_t out_sum, const PairContainer &row) {
                if (row.empty())
                    return 1.0;

                if (num_labels_in == num_labels_out && in_sum == out_sum)
                    return 1.0;

                int64_t n = in_sum + out_sum;
                double r_base = (*get_r)(n, row);

                if (r_base == 0.0)
                    return 1.0;

                double r_in = r_base * num_labels_in;
                double r_out = r_base * num_labels_out;
                double r = r_in + r_out;
                double lscaling_base = lgamma(r) / log(2);

                return compute_pval_r(r, r_in, r_out, lscaling_base, in_sum, out_sum, row);
            };
        }
    } else if (config.test_type == "poisson_binom") {
        // fit distribution
        std::vector<double> p;
        p.reserve(groups.size());
        for (size_t i = 0; i < groups.size(); ++i) {
            p.emplace_back(static_cast<double>(n_kmers[i]) / nelem);
        }

        common::logger->trace("p: {}", fmt::join(p, ","));

        common::logger->trace("Precomputing PMFs");
        std::vector<double> pmf_in { 1.0 };
        std::vector<double> pmf_out { 1.0 };
        std::vector<double> pmf_null { 1.0 };

        for (size_t i = 1; i <= groups.size(); ++i) {
            if (groups[i - 1]) {
                std::vector<double> pmf_out_cur(pmf_out.size() + 1);
                pmf_out_cur[0] = (1.0 - p[i - 1]) * pmf_out[0];
                pmf_out_cur[pmf_out.size()] = p[i - 1] * pmf_out.back();
                for (size_t k = 1; k < pmf_out.size(); ++k) {
                    pmf_out_cur[k] = p[i - 1] * pmf_out[k - 1] + (1.0 - p[i - 1]) * pmf_out[k];
                }
                std::swap(pmf_out_cur, pmf_out);
            } else {
                std::vector<double> pmf_in_cur(pmf_in.size() + 1);
                pmf_in_cur[0] = (1.0 - p[i - 1]) * pmf_in[0];
                pmf_in_cur[pmf_in.size()] = p[i - 1] * pmf_in.back();
                for (size_t k = 1; k < pmf_in.size(); ++k) {
                    pmf_in_cur[k] = p[i - 1] * pmf_in[k - 1] + (1.0 - p[i - 1]) * pmf_in[k];
                }
                std::swap(pmf_in_cur, pmf_in);
            }
            std::vector<double> pmf_null_cur(i + 1);
            pmf_null_cur[0] = (1.0 - p[i - 1]) * pmf_null[0];
            pmf_null_cur[pmf_null.size()] = p[i - 1] * pmf_null.back();
            for (size_t k = 1; k < pmf_null.size(); ++k) {
                pmf_null_cur[k] = p[i - 1] * pmf_null[k - 1] + (1.0 - p[i - 1]) * pmf_null[k];
            }
            std::swap(pmf_null_cur, pmf_null);
        }

        common::logger->trace("Precomputing p-values");
        std::vector<std::vector<double>> pvals(groups.size() + 1);
        pvals[0].emplace_back(1.0);

        double min_pval = 1.0;

        for (size_t n = 1; n < pvals.size(); ++n) {
            std::vector<double> probs;
            for (uint64_t s = 0; s < pmf_in.size(); ++s) {
                if (s > n)
                    break;
                uint64_t t = n - s;
                if (t < pmf_out.size()) {
                    probs.emplace_back(exp2(log2(pmf_in[s]) + log2(pmf_out[t]) - log2(pmf_null[n])));
                } else {
                    probs.emplace_back(0.0);
                }
            }

            for (uint64_t s = 0; s < probs.size(); ++s) {
                pvals[n].emplace_back(0.0);
                bool found = false;
                for (uint64_t sp = 0; sp < probs.size(); ++sp) {
                    if (pmf_in[sp] <= pmf_in[s]) {
                        found = true;
                        pvals[n][s] += probs[sp];
                    }
                }

                if (found)
                    min_pval = std::min(min_pval, pvals[n][s]);
            }
        }

        common::logger->trace("Min. p-value: {}", min_pval);
        if (min_pval >= config.family_wise_error_rate) {
            common::logger->warn("No significant p-values achievable");
            auto masked_graph_in = std::make_shared<MaskedDeBruijnGraph>(
                graph_ptr, std::make_unique<bitmap_vector>(std::move(indicator_in)), true,
                is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
            );

            auto masked_graph_out = std::make_shared<MaskedDeBruijnGraph>(
                graph_ptr, std::make_unique<bitmap_vector>(std::move(indicator_out)), true,
                is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
            );

            return std::make_tuple(masked_graph_in, masked_graph_out, PValStorage{}, nullptr);
        }

        compute_min_pval = [pvals,num_labels_in](int64_t, const PairContainer &row) {
            size_t front = row.size() - std::min(row.size(), num_labels_in);
            return std::min(pvals[row.size()][front],
                            pvals[row.size()].back());
        };

        compute_pval = [pvals,&groups](int64_t, int64_t, const auto &row) {
            uint64_t s = std::count_if(row.begin(), row.end(),
                                       [&](const auto &a) { return !groups[a.first]; });

            return pvals[row.size()][s];
        };

    } else if (config.test_type == "notest") {
        compute_min_pval = [&](int64_t n, const PairContainer&) { return n > 0 ? 0.0 : 1.0; };
        compute_pval = [&](int64_t, int64_t, const auto &row) { return row.size() ? 0.0 : 1.0; };
    } else {
        throw std::runtime_error("Test not implemented");
    }

    common::logger->trace("Computing minimum p-values");
    std::vector<std::pair<double, size_t>> m;
    std::vector<int64_t> m_sums;
    {
        std::vector<std::vector<std::pair<double, size_t>>> ms(num_threads + 1);
        std::vector<std::vector<tsl::hopscotch_map<std::vector<int64_t>, std::pair<double, size_t>, utils::VectorHash>>> ms_vecs(num_threads + 1);
        generate_rows([&](uint64_t row_i, const auto &raw_row, size_t bucket_idx) {
            int64_t in_sum = 0;
            int64_t out_sum = 0;
            int64_t out_stat_int = 0;
            bool in_kmer = false;
            bool out_kmer = false;
            PairContainer row;
            std::vector<int64_t> vals;
            if (kept[row_i]) {
                size_t count_in = 0;
                size_t count_out = 0;
                for (const auto &[j, raw_c] : raw_row) {
                    uint64_t c = raw_c;
                    if (count_maps.size())
                        c = count_maps[j].find(c)->second.first;

                    row.emplace_back(j, c);
                    vals.emplace_back(c);
                    if (groups[j]) {
                        out_sum += c;
                        ++count_out;
                    } else {
                        in_sum += c;
                        ++count_in;
                    }
                }

                if (count_in + count_out >= config.min_recurrence) {
                    double out_stat = static_cast<double>(out_sum) / out_kmers * in_kmers;
                    if (out_stat != 0)
                        out_stat_int = out_stat > 0 ? ceil(out_stat) : floor(out_stat);

                    in_kmer = count_in >= config.min_in_recurrence && count_out <= config.max_out_recurrence && in_sum > (out_kmers > 0 ? out_stat : 0.0);
                    out_kmer = count_out >= config.min_out_recurrence && count_in <= config.max_in_recurrence && in_sum < out_stat;
                }
            } else {
                return;
            }

            size_t n = config.test_type != "poisson_binom" ? in_sum + out_sum : row.size();
            if (!in_kmer && !out_kmer) {
                row.clear();
                vals.clear();
                n = 0;
            }

            double pval_min = 0;
            if (config.test_type != "gnb_exact") {
                assert(bucket_idx < ms.size());
                if (n >= ms[bucket_idx].size())
                    ms[bucket_idx].resize(n + 1, std::make_pair(1.1, 0));

                if (ms[bucket_idx][n].first == 1.1)
                    ms[bucket_idx][n].first = compute_min_pval(n, row);

                pval_min = ms[bucket_idx][n].first;
                ++ms[bucket_idx][n].second;
            } else {
                std::sort(vals.begin(), vals.end());
                if (n >= ms_vecs[bucket_idx].size())
                    ms_vecs[bucket_idx].resize(n + 1);

                auto find = ms_vecs[bucket_idx][n].find(vals);
                if (find == ms_vecs[bucket_idx][n].end()) {
                    find = ms_vecs[bucket_idx][n].try_emplace(
                        vals,
                        std::make_pair(compute_min_pval(n, row), 1)
                    ).first;
                } else {
                    ++find.value().second;
                }

                pval_min = find->second.first;
            }

            if (config.test_by_unitig) {
                bucket_idx = std::min(bucket_idx, pvals_buckets.size() - 1);
                uint64_t eff_size = bit_cast<uint64_t, int64_t>(in_sum - out_stat_int);
                if constexpr(preallocated) {
                    node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
                    eff_size_buckets[bucket_idx][node] = eff_size;
                    sum_buckets[bucket_idx][node] = n;
                    if (config.test_type == "gnb_exact")
                        pvals_min_buckets[bucket_idx][node] = bit_cast<uint64_t>(pval_min);
                } else {
                    eff_size_buckets[bucket_idx].push_back(eff_size);
                    sum_buckets[bucket_idx].push_back(n);
                    if (config.test_type == "gnb_exact")
                        pvals_min_buckets[bucket_idx].push_back(bit_cast<uint64_t>(pval_min));
                }
            }
        });

        common::logger->trace("Merging min p-value tables");
        if (config.test_type != "gnb_exact") {
            for (size_t i = 1; i < ms.size(); ++i) {
                size_t end = std::min(ms[0].size(), ms[i].size());
                for (size_t j = 0; j < end; ++j) {
                    assert(ms[0][j].first == 1.1 || ms[i][j].first == 1.1 || ms[0][j] == ms[i][j]);
                    ms[0][j].first = std::min(ms[0][j].first, ms[i][j].first);
                    ms[0][j].second += ms[i][j].second;
                }

                if (ms[i].size() > ms[0].size()) {
                    std::copy(std::make_move_iterator(ms[i].begin() + ms[0].size()),
                            std::make_move_iterator(ms[i].end()),
                            std::back_inserter(ms[0]));
                }
            }
            m = std::move(ms[0]);
            ms.resize(0);
        } else {
            std::vector<std::tuple<double, int64_t, size_t>> mn;
            for (size_t i = 0; i < ms_vecs.size(); ++i) {
                for (size_t n = 0; n < ms_vecs[i].size(); ++n) {
                    for (const auto &[vec, v] : ms_vecs[i][n]) {
                        mn.emplace_back(v.first, n, v.second);
                    }
                }
            }
            ms_vecs.resize(0);
            std::sort(mn.begin(), mn.end(), utils::GreaterFirst());

            m.reserve(mn.size());
            m_sums.reserve(mn.size());
            for (const auto &[p, s, c] : mn) {
                m.emplace_back(p, c);
                m_sums.emplace_back(s);
            }
        }
    }

    size_t k = 0;
    size_t n_cutoff = 1;

    auto combine_pvals = [dist=boost::math::cauchy()](const std::vector<double> &pvals) {
        double stat = 0.0;
        for (double pval : pvals) {
            stat += tan((0.5 - pval) * M_PI);
        }
        stat /= pvals.size();
        return boost::math::cdf(boost::math::complement(dist, stat));
    };

    std::vector<size_t> boundaries;
    if (config.test_by_unitig) {
        boundaries.reserve(eff_size_buckets.size() + 1);
        boundaries.emplace_back(0);
        for (size_t i = 0; i < eff_size_buckets.size(); ++i) {
            boundaries.emplace_back(boundaries.back() + eff_size_buckets[i].size());
        }
    }

    std::mutex pval_mu;
    std::unique_ptr<MaskedDeBruijnGraph> clean_masked_graph;
    sdsl::bit_vector unitig_start;
    if (!config.test_by_unitig) {
        std::tie(k, n_cutoff) = correct_pvals(m, nelem, m_sums);
        common::logger->trace("Picked: k: {}\tn: {}\tpval_min: {}", k, n_cutoff, config.family_wise_error_rate / k);
    } else {
        clean_masked_graph = std::make_unique<MaskedDeBruijnGraph>(
            graph_ptr,
            [&](node_index node) {
                return node != DeBruijnGraph::npos
                        && kept[AnnotatedDBG::graph_to_anno_index(node)];
            },
            true,
            is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
        );

        std::vector<std::pair<double, int64_t>> comb_msums;
        common::logger->trace("Allocating updated indicator");
        unitig_start = sdsl::bit_vector(graph_ptr->max_index() + 1, false);

        common::logger->trace("Combining min p-vals");
        clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
            std::vector<size_t> bucket_idxs;
            bucket_idxs.reserve(path.size());
            for (node_index node : path) {
                bucket_idxs.emplace_back(0);
                for (size_t i = 1; i < boundaries.size(); ++i) {
                    if (boundaries[i] > node)
                        break;

                    ++bucket_idxs.back();
                }
            }

            std::vector<double> pvals_min;
            pvals_min.reserve(path.size());

            int64_t n_sum = 0;
            {
                std::lock_guard<std::mutex> lock(pval_mu);
                for (size_t i = 0; i < path.size(); ++i) {
                    size_t bucket_idx = bucket_idxs[i];
                    size_t node_shift = path[i] - boundaries[bucket_idx];
                    size_t n = sum_buckets[bucket_idx][node_shift];
                    n_sum += n;

                    if (config.test_type == "gnb_exact") {
                        pvals_min.emplace_back(bit_cast<double, uint64_t>(pvals_min_buckets[bucket_idx][node_shift]));
                    } else {
                        pvals_min.emplace_back(m[n].first);
                    }
                }
            }

            double comb_pval = combine_pvals(pvals_min);
            if (comb_pval < config.family_wise_error_rate) {
                for (node_index node : path) {
                    set_bit(unitig_start.data(), node, parallel, MO_RELAXED);
                }
            }

            std::lock_guard<std::mutex> lock(pval_mu);
            comb_msums.emplace_back(comb_pval, n_sum);
        }, num_threads);

        std::sort(comb_msums.begin(), comb_msums.end(), utils::GreaterFirst());
        std::vector<std::pair<double, size_t>> comb_m;
        std::vector<int64_t> comb_m_sums;
        comb_m.reserve(comb_msums.size());
        comb_m_sums.reserve(comb_msums.size());
        for (const auto &[p, s] : comb_msums) {
            comb_m.emplace_back(p, 1);
            comb_m_sums.emplace_back(s);
        }

        std::tie(k, n_cutoff) = correct_pvals(comb_m, nelem, comb_m_sums);
        common::logger->trace("Picked: k: {}\tpval_min: {}", k, config.family_wise_error_rate / k);

        sum_buckets.resize(0);
        tmp_sum_buckets.resize(0);
        pvals_min_buckets.resize(0);
        tmp_min_buckets.resize(0);
    }

    if (k > 0) {
        common::logger->trace("Running differential tests");
        std::exception_ptr ex = nullptr;
        std::mutex ex_mu;
        bool keep_all_pvals = config.output_pvals || config.test_by_unitig;
        std::atomic_thread_fence(std::memory_order_release);

        generate_rows([&](uint64_t row_i, const auto &raw_row, size_t bucket_idx) {
            if (ex)
                return;

            int64_t in_sum = 0;
            int64_t out_sum = 0;
            bool in_kmer = false;
            bool out_kmer = false;
            PairContainer row;
            bucket_idx = std::min(bucket_idx, pvals_buckets.size() - 1);
            node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
            double pval = 1.1;
            if (kept[row_i]) {
                size_t count_in = 0;
                size_t count_out = 0;
                for (const auto &[j, raw_c] : raw_row) {
                    uint64_t c = raw_c;
                    if (count_maps.size())
                        c = count_maps[j].find(c)->second.first;

                    row.emplace_back(j, c);
                    if (groups[j]) {
                        out_sum += c;
                        ++count_out;
                    } else {
                        in_sum += c;
                        ++count_in;
                    }
                }

                if (count_in + count_out >= config.min_recurrence) {
                    double out_stat = static_cast<double>(out_sum) / out_kmers * in_kmers;
                    in_kmer = count_in >= config.min_in_recurrence && count_out <= config.max_out_recurrence && in_sum > (out_kmers > 0 ? out_stat : 0.0);
                    out_kmer = count_out >= config.min_out_recurrence && count_in <= config.max_in_recurrence && in_sum < out_stat;
                }
            } else {
                if (keep_all_pvals) {
                    auto &pvals = pvals_buckets[bucket_idx];
                    if constexpr(preallocated) {
                        pvals[node] = bit_cast<uint64_t>(pval);
                    } else {
                        pvals.push_back(bit_cast<uint64_t>(pval));
                    }
                }
                return;
            }

            size_t n = config.test_type != "poisson_binom" ? in_sum + out_sum : row.size();
            if (!in_kmer && !out_kmer) {
                row.clear();
                n = 0;
            }

            if (!config.output_pvals && config.test_by_unitig && !unitig_start[node]) {
                auto &pvals = pvals_buckets[bucket_idx];
                if constexpr(preallocated) {
                    pvals[node] = bit_cast<uint64_t>(pval);
                } else {
                    pvals.push_back(bit_cast<uint64_t>(pval));
                }
            }

            if (config.output_pvals
                    || (!config.test_by_unitig && n >= n_cutoff)
                    || (config.test_by_unitig && unitig_start[node])) {
                if (config.test_type != "gnb_exact" && m[n].first == 1.1) {
                    common::logger->error("in: {}\tout: {}\tn: {}", in_sum, out_sum, n);
                    throw std::runtime_error("Indexing invalid");
                }

                try {
                    pval = compute_pval(in_sum, out_sum, row);
                } catch (...) {
                    std::lock_guard<std::mutex> lock(ex_mu);
                    ex = std::current_exception();
                    return;
                }

                if (config.test_type != "gnb_exact" && pval < 1.1 && m[n].first - pval > 1e-10) {
                    common::logger->error("Min p-val estimate too high: min {} > cur {}\tn: {}\ttest: {}", m[n].first, pval, n, config.test_type);
                    throw std::runtime_error("Test failed");
                }

                if (config.test_type != "notest" && keep_all_pvals) {
                    auto &pvals = pvals_buckets[bucket_idx];
                    if constexpr(preallocated) {
                        pvals[node] = bit_cast<uint64_t>(pval);
                    } else {
                        pvals.push_back(bit_cast<uint64_t>(pval));
                    }
                }

                if (!config.test_by_unitig && in_kmer != out_kmer && pval * k < config.family_wise_error_rate) {
                    bool use_atomic = parallel && (node % 64 == 0);
                    set_bit((in_kmer ? indicator_in : indicator_out).data(), node, use_atomic, MO_RELAXED);
                }
            }
        });

        std::atomic_thread_fence(std::memory_order_acquire);

        if (ex)
            std::rethrow_exception(ex);

        unitig_start = sdsl::bit_vector();
        deallocate();

        if (!config.test_by_unitig) {
            boundaries.reserve(pvals_buckets.size() + 1);
            boundaries.emplace_back(0);
            for (size_t i = 0; i < pvals_buckets.size(); ++i) {
                boundaries.emplace_back(boundaries.back() + pvals_buckets[i].size());
            }
        }

        if (config.test_by_unitig) {
            common::logger->trace("Allocating k-mer bitmasks");
            indicator_in = sdsl::bit_vector(graph_ptr->max_index() + 1, false);
            indicator_out = sdsl::bit_vector(graph_ptr->max_index() + 1, false);

            common::logger->trace("Combining p-values within unitigs");
            std::atomic<uint64_t> num_sig_unitigs{0};

            std::atomic_thread_fence(std::memory_order_release);
            clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
                int64_t comb_eff_size = 0;
                std::vector<size_t> bucket_idxs;
                bucket_idxs.reserve(path.size());
                for (node_index node : path) {
                    bucket_idxs.emplace_back(0);
                    for (size_t i = 1; i < boundaries.size(); ++i) {
                        if (boundaries[i] > node)
                            break;

                        ++bucket_idxs.back();
                    }
                }

                // TODO: don't bother combining if the row sum is too small
                std::vector<double> pvals;
                pvals.reserve(path.size());
                {
                    std::lock_guard<std::mutex> lock(pval_mu);
                    for (size_t i = 0; i < path.size(); ++i) {
                        node_index node = path[i];
                        size_t bucket_idx = bucket_idxs[i];
                        size_t node_shift = node - boundaries[bucket_idx];

                        comb_eff_size += bit_cast<int64_t, uint64_t>(eff_size_buckets[bucket_idx][node_shift]);

                        double pval = bit_cast<double, uint64_t>(pvals_buckets[bucket_idx][node_shift]);
                        pvals.emplace_back(pval);
                    }
                }
                double comb_pval = combine_pvals(pvals);
                if (comb_pval * k < config.family_wise_error_rate && comb_eff_size != 0) {
                    ++num_sig_unitigs;
                    auto *data = comb_eff_size > 0 ? indicator_in.data() : indicator_out.data();
                    for (node_index node : path) {
                        set_bit(data, node, parallel, MO_RELAXED);
                    }
                }

                uint64_t comb_pval_enc = bit_cast<uint64_t, double>(comb_pval);

                std::lock_guard<std::mutex> lock(pval_mu);
                for (size_t i = 0; i < path.size(); ++i) {
                    size_t bucket_idx = bucket_idxs[i];
                    pvals_buckets[bucket_idx][path[i] - boundaries[bucket_idx]] = comb_pval_enc;
                }
            }, num_threads);
            std::atomic_thread_fence(std::memory_order_acquire);
            common::logger->trace("Found {} significant unitigs before correction", num_sig_unitigs);

            eff_size_buckets.resize(0);
            tmp_eff_buckets.resize(0);
        }
    } else {
        deallocate();
    }

    std::unique_ptr<utils::TempFile> tmp_file;
    PValStorage pvals;

    if constexpr(std::is_same_v<PValStorage, std::vector<uint64_t>>) {
        pvals = std::move(pvals_buckets[0]);
    }

    if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
        if (config.output_pvals) {
            common::logger->trace("Merging buckets");
            for (size_t i = 1; i < pvals_buckets.size(); ++i) {
                for (size_t j = 0; j < pvals_buckets[i].size(); ++j) {
                    pvals_buckets[0].push_back(pvals_buckets[i][j]);
                }
                pvals_buckets[i].close();
            }
        }
        pvals = std::move(pvals_buckets[0]);
        pvals_buckets.resize(0);
        tmp_file = std::move(tmp_buckets[0]);
        tmp_buckets.resize(0);
    }

    common::logger->trace("Done! Assembling contigs.");

    auto masked_graph_in = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr, std::make_unique<bitmap_vector>(std::move(indicator_in)), true,
        is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
    );

    auto masked_graph_out = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr, std::make_unique<bitmap_vector>(std::move(indicator_out)), true,
        is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
    );

    return std::make_tuple(masked_graph_in, masked_graph_out, std::move(pvals), std::move(tmp_file));
}

template <class PValStorage>
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, PValStorage, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                         const std::vector<std::string> &files,
                         const tsl::hopscotch_set<Label> &labels_in,
                         const tsl::hopscotch_set<Label> &labels_out,
                         const DifferentialAssemblyConfig &config,
                         size_t num_threads,
                         std::filesystem::path tmp_dir,
                         size_t num_parallel_files) {
    num_parallel_files = get_num_threads();
    common::logger->trace("Labels in: {}", fmt::join(labels_in, ","));
    common::logger->trace("Labels out: {}", fmt::join(labels_out, ","));

    size_t total_labels = labels_in.size() + labels_out.size();

    using ColumnValuesMem = std::vector<std::unique_ptr<const sdsl::int_vector<>>>;
    using ColumnValuesDisk = std::vector<std::unique_ptr<const sdsl::int_vector_buffer<>>>;

    std::variant<ColumnValuesMem, ColumnValuesDisk> column_values_all;

    if (config.test_by_unitig) {
        column_values_all = ColumnValuesMem(total_labels);
    } else {
        column_values_all = ColumnValuesDisk(total_labels);
    }

    return std::visit([&](auto&& column_values_all) {
        using ColumnValues = typename std::decay<decltype(column_values_all)>::type;
        using ValuesContainerPtr = typename ColumnValues::value_type;
        using ValuesContainer = typename std::decay<typename ValuesContainerPtr::element_type>::type;
        using value_type = typename ValuesContainer::value_type;
        using PairContainer = Vector<std::pair<uint64_t, value_type>>;

        std::vector<std::unique_ptr<const bit_vector>> columns_all(total_labels);
        std::vector<bool> groups(total_labels);
        annot::ColumnCompressed<>::load_columns_and_values(
            files,
            [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector>&& column, ValuesContainer&& column_values) {
                groups[offset] = !labels_in.count(label);
                columns_all[offset].reset(column.release());
                column_values_all[offset] = std::make_unique<ValuesContainer>(std::move(column_values));
            },
            num_parallel_files
        );

        uint8_t max_width = 0;
        for (const auto &col_vals : column_values_all) {
            max_width = std::max(max_width, col_vals->width());
        }

        auto generate_rows = [&](const auto &callback) {
            utils::call_rows<std::unique_ptr<const bit_vector>,
                             std::unique_ptr<const ValuesContainer>,
                             PairContainer, false>(columns_all, column_values_all, callback);
        };

        bool parallel = get_num_threads() > 1;

        return mask_nodes_by_label_dual<value_type, PValStorage>(
            graph_ptr,
            [&](const std::vector<size_t> &min_counts, sdsl::bit_vector *kept) -> std::vector<VectorMap<uint64_t, size_t>> {
                common::logger->trace("Calculating count histograms");
                std::vector<std::vector<VectorMap<uint64_t, size_t>>> hists_map_p(num_parallel_files);
                for (auto &hists_map : hists_map_p) {
                    hists_map.resize(groups.size());
                }

                std::atomic_thread_fence(std::memory_order_release);
                generate_rows([&](uint64_t row_i, const auto &row, size_t thread_id) {
                    if (row.empty()) {
                        if (kept)
                            unset_bit(kept->data(), row_i, parallel, std::memory_order_relaxed);

                        return;
                    }

                    bool found = false;
                    Vector<uint64_t> counts(groups.size());
                    for (const auto &[j, raw_c] : row) {
                        counts[j] = raw_c;
                        found |= min_counts.empty() || raw_c >= min_counts[j];
                    }

                    if (!found) {
                        if (kept)
                            unset_bit(kept->data(), row_i, parallel, std::memory_order_relaxed);

                        return;
                    }

                    for (size_t j = 0; j < counts.size(); ++j) {
                        ++hists_map_p[thread_id][j][counts[j]];
                    }
                });
                std::atomic_thread_fence(std::memory_order_acquire);

                common::logger->trace("Merging histograms");
                for (size_t thread_id = 1; thread_id < hists_map_p.size(); ++thread_id) {
                    for (size_t j = 0; j < hists_map_p[thread_id].size(); ++j) {
                        for (const auto &[k, c] : hists_map_p[thread_id][j]) {
                            hists_map_p[0][j][k] += c;
                        }
                    }
                }

                return hists_map_p[0];
            },
            generate_rows,
            groups,
            config, num_threads, tmp_dir, num_parallel_files,
            [&]() {
                columns_all.clear();
                column_values_all.clear();
            },
            max_width
        );
    }, column_values_all);
}

template
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, std::vector<uint64_t>, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<std::vector<uint64_t>>(std::shared_ptr<const DeBruijnGraph>,
                         const std::vector<std::string> &,
                         const tsl::hopscotch_set<Label> &,
                         const tsl::hopscotch_set<Label> &,
                         const DifferentialAssemblyConfig &,
                         size_t,
                         std::filesystem::path,
                         size_t);
template
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, sdsl::int_vector_buffer<64>, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<sdsl::int_vector_buffer<64>>(std::shared_ptr<const DeBruijnGraph>,
                         const std::vector<std::string> &,
                         const tsl::hopscotch_set<Label> &,
                         const tsl::hopscotch_set<Label> &,
                         const DifferentialAssemblyConfig &,
                         size_t,
                         std::filesystem::path,
                         size_t);

template <class PValStorage>
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, PValStorage, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(const AnnotatedDBG &anno_graph,
                         const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &labels_in,
                         const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &labels_out,
                         const DifferentialAssemblyConfig &config,
                         size_t num_threads,
                         std::filesystem::path tmp_dir) {
    num_threads = get_num_threads();
    common::logger->trace("Labels in: {}", fmt::join(labels_in, ","));
    common::logger->trace("Labels out: {}", fmt::join(labels_out, ","));

    size_t total_labels = labels_in.size() + labels_out.size();
    const auto &annotation = anno_graph.get_annotator();

    std::vector<bool> groups;
    groups.reserve(total_labels);
    for (const auto &label : annotation.get_label_encoder().get_labels()) {
        groups.emplace_back(!labels_in.count(label));
    }

    const auto &matrix = annotation.get_matrix();
    using Row = annot::matrix::BinaryMatrix::Row;

    if (const auto *int_matrix = dynamic_cast<const annot::matrix::IntMatrix*>(&matrix)) {
        using value_type = annot::matrix::IntMatrix::Value;
        return mask_nodes_by_label_dual<value_type, PValStorage>(
            std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr()),
            [&](const std::vector<size_t> &min_counts, sdsl::bit_vector *unmark_discarded) {
                return int_matrix->get_histograms(min_counts, unmark_discarded);
            },
            [&](const auto &callback) {
                int_matrix->call_row_values(callback, false /* ordered */);
            },
            groups,
            config, num_threads, tmp_dir, num_threads
        );
    } else {
        auto generate_bit_rows = [&](const auto &callback) {
            size_t batch_size = (matrix.num_rows() + num_threads - 1) / num_threads;
            size_t rows_per_update = 10000;
            ProgressBar progress_bar(matrix.num_rows(), "Streaming rows", std::cerr, !common::get_verbose());
            #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
            for (size_t k = 0; k < matrix.num_rows(); k += batch_size) {
                size_t begin = k;
                size_t end = std::min(begin + batch_size, matrix.num_rows());
                size_t thread_id = k / batch_size;
                std::vector<Row> row_ids;
                for (size_t row_batch_i = begin; row_batch_i < end; row_batch_i += rows_per_update) {
                    size_t row_batch_end = std::min(row_batch_i + rows_per_update, end);
                    row_ids.resize(row_batch_end - row_batch_i);
                    std::iota(row_ids.begin(), row_ids.end(), row_batch_i);
                    auto set_bits = matrix.get_rows(row_ids);
                    for (size_t i = 0; i < row_ids.size(); ++i) {
                        callback(i + row_batch_i, set_bits[i], thread_id);
                    }
                    progress_bar += row_ids.size();
                }
            }
        };

        bool parallel = get_num_threads() > 1;

        return mask_nodes_by_label_dual<uint64_t, PValStorage>(
            std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr()),
            [&](const std::vector<size_t> &, sdsl::bit_vector *kept) -> std::vector<VectorMap<uint64_t, size_t>> {
                common::logger->trace("Calculating count histograms");
                std::vector<std::vector<VectorMap<uint64_t, size_t>>> hists_map_p(num_threads);
                for (auto &hists_map : hists_map_p) {
                    hists_map.resize(groups.size());
                }
                std::atomic_thread_fence(std::memory_order_release);
                generate_bit_rows([&](uint64_t row_i, const auto &set_bits, size_t thread_id) {
                    if (set_bits.empty()) {
                        if (kept)
                            unset_bit(kept->data(), row_i, parallel, std::memory_order_relaxed);

                        return;
                    }

                    std::vector<bool> container(groups.size());
                    for (auto j : set_bits) {
                        container[j] = true;
                    }

                    for (size_t j = 0; j < container.size(); ++j) {
                        ++hists_map_p[thread_id][j][container[j]];
                    }
                });
                std::atomic_thread_fence(std::memory_order_acquire);

                common::logger->trace("Merging histograms");
                for (size_t thread_id = 1; thread_id < hists_map_p.size(); ++thread_id) {
                    for (size_t j = 0; j < hists_map_p[thread_id].size(); ++j) {
                        for (const auto &[k, c] : hists_map_p[thread_id][j]) {
                            hists_map_p[0][j][k] += c;
                        }
                    }
                }

                return hists_map_p[0];
            },
            [&](const auto &callback) {
                generate_bit_rows([&](uint64_t row_i, const auto &set_bits, size_t thread_id) {
                    Vector<std::pair<uint64_t, uint64_t>> row;
                    row.reserve(set_bits.size());
                    for (auto j : set_bits) {
                        row.emplace_back(j, 1);
                    }

                    callback(row_i, row, thread_id);
                });
            },
            groups,
            config, num_threads, tmp_dir, num_threads
        );
    }
}

template
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, std::vector<uint64_t>, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<std::vector<uint64_t>>(const AnnotatedDBG &,
                         const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &,
                         const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &,
                         const DifferentialAssemblyConfig &,
                         size_t,
                         std::filesystem::path);
template
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, sdsl::int_vector_buffer<64>, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<sdsl::int_vector_buffer<64>>(const AnnotatedDBG &,
                         const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &,
                         const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label> &,
                         const DifferentialAssemblyConfig &,
                         size_t,
                         std::filesystem::path);

} // namespace graph
} // namespace mtg


