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
    num_parallel_files = std::min(num_parallel_files, num_threads);
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
            common::logger->warn("No significant p-values achievable");

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

        common::logger->trace("{}: sum: {}\tmax_obs: {}\tmin_cutoff: {}\tmax_cutof: {}",
                              j, sums[j], max_obs_vals[j], min_counts[j], check_cutoff[j]);
    }

    int64_t total_kmers = in_kmers + out_kmers;

    common::logger->trace("Number of kept unique k-mers: {}\tNumber of kept k-mers: {}",
                          nelem, total_kmers);

    common::logger->trace("Allocating p-value storage");
    auto nullpval = bit_cast<uint64_t>(double(1.1));

    std::vector<PValStorage> pvals_buckets;
    std::vector<std::unique_ptr<utils::TempFile>> tmp_buckets;
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
    } else if (config.test_type == "nbinom_exact" || config.test_type == "zinb_exact" || config.test_type == "gnb_exact") {
        common::logger->trace("Fitting per-sample negative binomial distributions");
        auto get_rp = [&](size_t j, auto begin, auto end) {
            double mu = 0;
            double var = 0;
            size_t total = 0;
            std::for_each(begin, end, [&](const auto &a) {
                const auto &[k, c] = a;
                mu += k * c;
                var += k * k * c;
                total += c;
            });
            mu /= total;
            double mu2 = mu * mu;
            var = (var - mu2 * total) / (total - 1);

            if (mu >= var) {
                common::logger->warn("Fit {} failed, falling back to Poisson: mu: {} >= var: {}", j, mu, var);
                return std::make_pair(0.0, 1.0);
            }

            double r_guess = mu * mu / (var - mu);
            double p_guess = mu / var;
            double r = r_guess;
            try {
                auto get_dl_ddl = [&](double r) {
                    double dl = (log(r) - log(r + mu) - boost::math::digamma(r)) * total;
                    double ddl = (1.0 / r - 1.0 / (r + mu) - boost::math::trigamma(r)) * total;
                    std::for_each(begin, end, [&](const auto &a) {
                        const auto &[k, c] = a;
                        dl += boost::math::digamma(k + r) * c;
                        ddl += boost::math::trigamma(k + r) * c;
                    });
                    return std::make_pair(dl, ddl);
                };
                r = boost::math::tools::newton_raphson_iterate(get_dl_ddl, r_guess,
                                                               std::numeric_limits<double>::min(), mu * total, 30);
                auto [dl, ddl] = get_dl_ddl(r);
                if (ddl > 0) {
                    common::logger->error("Found local minimum instead: r: {}\tdl: {}\tddl: {}", r, dl, ddl);
                    throw std::domain_error("FFF");
                }
            } catch (std::exception &e) {
                common::logger->warn("Caught exception for sample {}: Falling back to initial guess", j);
                common::logger->warn("{}", e.what());
                r = r_guess;
            }

            double p = r / (r + mu);
            common::logger->trace("{}: size: {}\tmax_val: {}\tmu: {}\tvar: {}\tmoment est:\tr: {}\tp: {}\tmle: r: {}\tp: {}",
                                  j, sums[j], (end - 1)->first, mu, var, r_guess, p_guess, r, p);

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
            nb_params[j] = get_rp(j, hist.begin(), it);

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
        double real_sum = 0.0;
        double real_sum_sq = 0.0;
        for (size_t j = 0; j < groups.size(); ++j) {
            double f = target_sum / sums[j];
            common::logger->trace("{}: scale factor: {}", j, f);
            std::for_each(hists[j].cbegin(), hist_its[j], [&](const auto &a) {
                const auto &[k, c] = a;
                if (k != 0) {
                    double val = f * k;
                    real_sum += val * c;
                    real_sum_sq += val * val * c;
                }
            });
        }

        double r_guess = 0.0;
        double target_mu = 0.0;
        double target_var = 0.0;

        auto map_value = [&](uint64_t k, double scale, const auto &old_dist, const auto &new_dist) {
            double cdf = boost::math::cdf(old_dist, k);

            if (cdf < 1.0)
                return boost::math::quantile(new_dist, cdf);

            double ccdf = boost::math::cdf(boost::math::complement(old_dist, k));

            if (ccdf > 0.0)
                return boost::math::quantile(boost::math::complement(new_dist, ccdf));

            return ceil(scale * k);
        };

        // nelem * r * groups.size() / p - target_sum * groups.size() / (1 - p) = 0
        // nelem * r / p = target_sum / (1 - p)
        // r / target_mu = p / (1 - p)
        // r / target_mu = p + p * r / target_mu
        // p = r / target_mu / (r / target_mu + 1)
        // p = r / (r + target_mu)
        target_mu = real_sum / matrix_size;
        double target_mu2 = target_mu * target_mu;
        target_var = (real_sum_sq - target_mu2 * matrix_size) / (matrix_size - 1);

        if (target_var > target_mu) {
            r_guess = target_mu2 / (target_var - target_mu);
        } else {
            common::logger->error("Data is underdispersed: {} >= {}. Use poisson_exact instead", target_mu, target_var);
            throw std::domain_error("Fit failed");
        }

        auto get_dl_ddl = [&](double r) {
            double dl = (log(r) - log(r + target_mu) - boost::math::digamma(r)) * matrix_size;
            double ddl = (1.0 / r - 1.0 / (r + target_mu) - boost::math::trigamma(r)) * matrix_size;
            for (size_t j = 0; j < groups.size(); ++j) {
                double f = target_sum / sums[j];
                std::for_each(hists[j].cbegin(), hist_its[j], [&](const auto &a) {
                    const auto &[k, c] = a;
                    dl += boost::math::digamma(f * k + r) * c;
                    ddl += boost::math::trigamma(f * k + r) * c;
                });
            }
            return std::make_pair(dl, ddl);
        };
        double r_map = boost::math::tools::newton_raphson_iterate(get_dl_ddl, r_guess, 0.0, real_sum, 30);
        auto [dl, ddl] = get_dl_ddl(r_map);
        if (ddl > 0) {
            common::logger->error("Found local minimum instead: r: {}\tdl: {}\tddl: {}", r_map, dl, ddl);
            throw std::domain_error("GGG");
        }

        double target_p = r_map / (r_map + target_mu);
        double fit_mu = r_map * (1.0 - target_p) / target_p;
        double fit_var = fit_mu / target_p;
        common::logger->trace("Common params: r_guess: {}\tr: {}\tp: {}\tmu: {}\tvar: {}",
                              r_guess, r_map, target_p, fit_mu, fit_var);

        double r_in = r_map * num_labels_in;
        double r_out = r_map * num_labels_out;
        common::logger->trace("Fits: in: {}\tout: {}\tp: {}", r_in, r_out, target_p);

        double lscaling_base = lgamma(r_in + r_out) - lgamma(r_in) - lgamma(r_out);
        double mu1 = r_in * (1.0 - target_p) / target_p;
        double mu2 = r_out * (1.0 - target_p) / target_p;

        if (fit_var / fit_mu - 1.0 < 1e-5)
            common::logger->warn("Fit parameters are close to a Poisson distribution");

        count_maps.resize(groups.size());

        common::logger->trace("Computing quantile maps");
        #pragma omp parallel for num_threads(num_parallel_files)
        for (size_t j = 0; j < groups.size(); ++j) {
            const auto &[r, p] = nb_params[j];

            common::logger->trace("{}\tApproximating NB({}, {}) with NB({}, {})",
                                  j, r, p, r_map, target_p);

            boost::math::negative_binomial nb_out(r_map, target_p);
            double scale = target_sum / sums[j];

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
        int64_t total_kmers_sq = 0;
        size_t total_nonzeros = 0;
        for (size_t j = 0; j < groups.size(); ++j) {
            for (const auto &[k, m] : count_maps[j]) {
                const auto &[v, c] = m;
                if (groups[j]) {
                    out_kmers += v * c;
                } else {
                    in_kmers += v * c;
                }

                if (v > 0)
                    total_nonzeros += c;

                total_kmers_sq += v * v * c;
            }
        }
        total_kmers = in_kmers + out_kmers;
        common::logger->trace("  Scaled Totals: in: {}\tout: {}", in_kmers, out_kmers);

        int64_t in_sum = 0;
        int64_t out_sum = 0;
        for (size_t j = 0; j < groups.size(); ++j) {
            if (hists[j].size()) {
                uint64_t k = count_maps[j].find(hists[j].back().first)->second.first;
                if (groups[j]) {
                    out_sum += k;
                } else {
                    in_sum += k;
                }
            }
        }
        int64_t total_sum = in_sum + out_sum;

        if (config.test_type == "nbinom_exact") {
            auto get_deviance1 = [mu=mu1,invphi=r_in,logmuinvphi=log(mu1+r_in),logmu=log(mu1),zv=2.0*r_in*log((r_in*mu1)/r_in)](double y) {
                if (y == 0)
                    return zv;

                double logyshift = logmuinvphi - log(y + invphi);
                return 2.0 * (y * (log(y) - logmu + logyshift) + invphi * logyshift);
            };
            auto get_deviance2 = [mu=mu2,invphi=r_out,logmuinvphi=log(mu2+r_out),logmu=log(mu2),zv=2.0*r_out*log((r_out*mu2)/r_out)](double y) {
                if (y == 0)
                    return zv;

                double logyshift = logmuinvphi - log(y + invphi);
                return 2.0 * (y * (log(y) - logmu + logyshift) + invphi * logyshift);
            };

            std::vector<double> deviances1(total_sum + 1);
            std::vector<double> deviances2(total_sum + 1);
            for (int64_t s = 0; s <= total_sum; ++s) {
                deviances1[s] = get_deviance1(s);
                deviances2[s] = get_deviance2(s);
            }
            compute_min_pval = [lscaling_base,r_in,r_out,target_p,deviances1,deviances2](int64_t n, const PairContainer&) {
                if (n == 0)
                    return 1.0;

                double lscaling = lgamma(n + 1) - lgamma(r_in + r_out + n) + lscaling_base;
                std::vector<double> devs(n + 1);

                for (int64_t s = 0; s <= n; ++s) {
                    devs[s] = deviances1[s] + deviances2[n - s];
                }

                double max_d = *std::max_element(devs.begin(), devs.end());
                double pval = 0.0;
                int64_t s = 0;
                if (devs[s] == max_d) {
                    int64_t t = n;
                    double rs = r_in;
                    double rt = r_out + n;
                    double base = lscaling - lgamma(n + 1) + lgamma(rs) + lgamma(rt);
                    pval += exp(base);
                    for (++s; s <= n; ++s) {
                        if (devs[s] == max_d) {
                            --t;
                            ++rs;
                            --rt;
                            base += log(t) - log(s) + log(rs - 1) - (rt > 1 ? log(rt - 1) : lgamma(rt + 1) - lgamma(rt));
                            pval += exp(base);
                        } else {
                            break;
                        }
                    }
                }

                int64_t sp = n;
                if (devs[sp] == max_d) {
                    int64_t t = 0;
                    double rs = r_in + sp;
                    double rt = r_out;
                    double base = lscaling - lgamma(n + 1) + lgamma(rs) + lgamma(rt);
                    pval += exp(base);
                    for (--sp; sp >= s; --sp) {
                        if (devs[sp] == max_d) {
                            ++t;
                            --rs;
                            ++rt;
                            base += log(sp) - log(t) - (rs > 1 ? log(rs - 1) : lgamma(rs + 1) - lgamma(rs)) + log(rt - 1);
                            pval += exp(base);
                        } else {
                            break;
                        }
                    }
                }

                return std::min(1.0, pval);
            };

            compute_pval = [lscaling_base,r_in,r_out,target_p,deviances1,deviances2](int64_t in_sum, int64_t out_sum, const auto &row) {
                if (row.empty())
                    return 1.0;

                int64_t n = in_sum + out_sum;

                double lscaling = lgamma(n + 1) - lgamma(r_in + r_out + n) + lscaling_base;

                double in_sum_total_deviance = deviances1[in_sum] + deviances2[out_sum];

                double pval = 0.0;
                int64_t s = 0;
                int64_t t = n;
                if (deviances1[s] + deviances2[t] >= in_sum_total_deviance) {
                    double rs = r_in;
                    double rt = r_out + n;
                    double base = lscaling - lgamma(n + 1) + lgamma(rs) + lgamma(rt);
                    pval += exp(base);
                    for (++s; s <= n; ++s) {
                        if (deviances1[s] + deviances2[--t] >= in_sum_total_deviance) {
                            ++rs;
                            --rt;
                            base += log(t) - log(s) + log(rs - 1) - (rt > 1 ? log(rt - 1) : lgamma(rt + 1) - lgamma(rt));
                            pval += exp(base);
                        } else {
                            break;
                        }
                    }
                }

                int64_t sp = n;
                t = 0;
                if (deviances1[sp] + deviances2[t] >= in_sum_total_deviance) {
                    double rs = r_in + sp;
                    double rt = r_out;
                    double base = lscaling - lgamma(n + 1) + lgamma(rs) + lgamma(rt);
                    pval += exp(base);
                    for (--sp; sp >= s; --sp) {
                        if (deviances1[sp] + deviances2[++t] >= in_sum_total_deviance) {
                            --rs;
                            ++rt;
                            base += log(sp) - log(t) - (rs > 1 ? log(rs - 1) : lgamma(rs + 1) - lgamma(rs)) + log(rt - 1);
                            pval += exp(base);
                        } else {
                            break;
                        }
                    }
                }

                return std::min(1.0, pval);
            };
        } else if (config.test_type == "zinb_exact") {
            common::logger->trace("Computing merged histograms");
            std::mutex hist_mu;
            tsl::hopscotch_map<uint64_t, size_t> merged_hist_in;
            tsl::hopscotch_map<uint64_t, size_t> merged_hist_out;
            tsl::hopscotch_map<uint64_t, size_t> merged_hist_null;
            int64_t p12_sum = 0;
            int64_t p11_sum = 0;
            double m0 = 0;
            double m1 = 0;
            for (const auto &[k, m] : count_maps[0]) {
                const auto &[v, c] = m;
                m0 += v * c;
            }
            for (const auto &[k, m] : count_maps[1]) {
                const auto &[v, c] = m;
                m1 += v * c;
            }
            m0 /= nelem;
            m1 /= nelem;
            generate_rows([&](uint64_t row_i, const auto &raw_row, size_t) {
                int64_t in_sum = 0;
                int64_t out_sum = 0;
                bool in_kmer = false;
                bool out_kmer = false;
                if (kept[row_i]) {
                    size_t count_in = 0;
                    size_t count_out = 0;
                    double p0 = -m0;
                    double p1 = -m1;
                    for (const auto &[j, raw_c] : raw_row) {
                        uint64_t c = raw_c;
                        if (count_maps.size())
                            c = count_maps[j].find(c)->second.first;

                        if (j == 0) {
                            p0 = c - m0;
                        } else if (j == 1) {
                            p1 = c - m1;
                        }

                        if (groups[j]) {
                            out_sum += c;
                            ++count_out;
                        } else {
                            in_sum += c;
                            ++count_in;
                        }
                    }

                    p12_sum += p0 * p1;
                    p11_sum += p0 * p0;

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
                std::lock_guard<std::mutex> lock(hist_mu);
                ++merged_hist_in[in_sum];
                ++merged_hist_out[out_sum];
                ++merged_hist_null[n];
            }, false);

            double obs_cov = static_cast<double>(p12_sum) / (nelem - 1);
            double obs_11_cov = static_cast<double>(p11_sum) / (nelem - 1);
            common::logger->trace("Var: {}\tExp cov: {}\tObs11 cov: {}\tobs12 cov: {}",
                                  r_map * (1.0 - target_p) / target_p / target_p,
                                  (1.0 - target_p) / target_p / target_p,
                                  obs_11_cov,
                                  obs_cov);

            auto get_rp_zinb = [&](const tsl::hopscotch_map<uint64_t, size_t> &hist_in,
                                   const tsl::hopscotch_map<uint64_t, size_t> &hist_out) {
                auto get_stats = [&](const auto &hist) {
                    size_t total_nzeros = 0;
                    size_t total_zeros = 0;
                    size_t sum = 0;
                    for (const auto &[k, c] : hist) {
                        if (k != 0) {
                            sum += k * c;
                            total_nzeros += c;
                        } else {
                            total_zeros += c;
                        }
                    }
                    size_t total = total_nzeros + total_zeros;

                    return std::make_tuple(total_nzeros, total_zeros, sum, total);
                };

                auto [total_nzeros_in, total_zeros_in, sum_in, total_in] = get_stats(hist_in);
                auto [total_nzeros_out, total_zeros_out, sum_out, total_out] = get_stats(hist_out);

                auto get_r = [&](const auto &hist, size_t total, double total_nzeros, double p, double r_guess, double sum) {
                    return boost::math::tools::newton_raphson_iterate([&](double r) {
                        double dl = (log(p) - boost::math::digamma(r)) * total_nzeros + log(total_nzeros / total) * (total - total_nzeros);
                        double ddl = -boost::math::trigamma(r) * total_nzeros;
                        for (const auto &[k, c] : hist) {
                            if (k != 0) {
                                dl += boost::math::digamma(r + k) * c;
                                ddl += boost::math::trigamma(r + k) * c;
                            }
                        }

                        return std::make_pair(dl, ddl);
                    }, r_guess, std::numeric_limits<double>::min(), sum, 30);
                };

                double p = boost::math::tools::newton_raphson_iterate([&](double p) {
                    auto r_in = get_r(hist_in, total_in, total_nzeros_in, p, r_map * num_labels_in, sum_in);
                    auto r_out = get_r(hist_out, total_out, total_nzeros_out, p, r_map * num_labels_out, sum_out);

                    double dl_in = r_in * total_nzeros_in / p / (1.0 - pow(p, r_in)) - sum_in / (1.0 - p);
                    double dl_out = r_out * total_nzeros_out / p / (1.0 - pow(p, r_out)) - sum_out / (1.0 - p);
                    double ddl_in = -r_in * total_nzeros_in / p / p / (1.0 - pow(p, r_in)) / (1.0 - pow(p, r_in)) * ((1.0-pow(p, r_in)) - p * r_in * pow(p, r_in-1)) - static_cast<double>(sum_in) / (1.0-p) / (1.0-p);
                    double ddl_out = -r_out * total_nzeros_out / p / p / (1.0 - pow(p, r_out)) / (1.0 - pow(p, r_out)) * ((1.0-pow(p, r_out)) - p * r_out * pow(p, r_out-1)) - static_cast<double>(sum_out) / (1.0-p) / (1.0-p);
                    return std::make_pair(dl_in + dl_out, ddl_in + ddl_out);
                }, target_p, 0.0, 1.0, 30);

                double r_in = get_r(hist_in, total_in, total_nzeros_in, p, r_map * num_labels_in, sum_in);
                double r_out = get_r(hist_out, total_out, total_nzeros_out, p, r_map * num_labels_out, sum_out);

                double pi_in = (static_cast<double>(total_in - total_nzeros_in) / total_in - pow(p, r_in)) / (1.0 - pow(p, r_in));
                double pi_out = (static_cast<double>(total_out - total_nzeros_out) / total_out - pow(p, r_out)) / (1.0 - pow(p, r_out));

                return std::make_tuple(r_in, pi_in, r_out, pi_out, p);
            };

            auto [r_in, pi_in, r_out, pi_out, p] = get_rp_zinb(merged_hist_in, merged_hist_in);

            common::logger->trace("r: ({}, {})\tp: {}\tpi: ({}, {})",
                                  r_in, r_out, p, pi_in, pi_out);

            double mu1 = r_in * (1.0 - p) / p * (1.0 - pi_in);
            double mu2 = r_out * (1.0 - p) / p * (1.0 - pi_out);

            auto get_zv = [](double r, double mu, double pi) {
                return -2 * log(pi + (1.0 - pi) * pow((1.0-pi)*r/(mu + (1.0-pi)*r), r));
            };

            auto get_deviance1 = [r=r_in,mu=mu1,mpir=(1.0-pi_in)*r_in,zv=get_zv(r_in,mu1,pi_in)](double y) {
                if (y == 0)
                    return zv;

                return 2.0 * ((r + y) * log((mu + mpir)/(y + mpir)) + y * log(y / mu));
            };
            auto get_deviance2 = [r=r_out,mu=mu2,mpir=(1.0-pi_out)*r_out,zv=get_zv(r_out,mu2,pi_out)](double y) {
                if (y == 0)
                    return zv;

                return 2.0 * ((r + y) * log((mu + mpir)/(y + mpir)) + y * log(y / mu));
            };

            std::vector<double> deviances1(total_sum + 1);
            std::vector<double> deviances2(total_sum + 1);
            for (int64_t s = 0; s <= total_sum; ++s) {
                deviances1[s] = get_deviance1(s);
                deviances2[s] = get_deviance2(s);
            }

            auto get_lprob_zinb = [](int64_t s, double r, double p, double pi) {
                return s == 0
                    ? log(pi + (1.0 - pi) * pow(p, r))
                    : log1p(-pi) + lgamma(r + s) - lgamma(r) - lgamma(s + 1) + r * log(p) + log1p(-p) * s;
            };

            auto get_lprob_ztnb = [get_lprob_zinb](int64_t s, double r_in, double r_out, double p, double pi_in, double pi_out) {
                double pval = 0.0;
                for (int64_t ss = 0; ss <= s; ++ss) {
                    int64_t tt = s - ss;
                    pval += exp(get_lprob_zinb(ss, r_in, p, pi_in) + get_lprob_zinb(tt, r_out, p, pi_out));
                }
                return log(pval);
                // double lp1 = get_lprob_zinb(0, r_in, p, pi_in) + get_lprob_zinb(s, r_out, p, pi_out);
                // double lp2 = get_lprob_zinb(s, r_in, p, pi_in) + get_lprob_zinb(0, r_out, p, pi_out);
                // double lp3 = log1p(-pi_in) + log1p(-pi_out) + get_lprob_zinb(s, r_in + r_out, p, 0.0);
                // return log(exp(lp1) + exp(lp2) + exp(lp3));
            };

            compute_min_pval = [r_in,p_in=p,pi_in,r_out,p_out=p,pi_out,get_lprob_zinb,get_lprob_ztnb,deviances1,deviances2](int64_t n, const PairContainer&) {
                if (n == 0)
                    return 1.0;

                std::vector<double> devs(n + 1);

                for (int64_t s = 0; s <= n; ++s) {
                    devs[s] = deviances1[s] + deviances2[n - s];
                }

                double max_d = *std::max_element(devs.begin(), devs.end());
                double pval = 0.0;
                double base_lprob = get_lprob_ztnb(n, r_in, r_out, p_in, pi_in, pi_out);
                for (int64_t s = 0; s <= n; ++s) {
                    if (devs[s] == max_d)
                        pval += exp(get_lprob_zinb(s, r_in, p_in, pi_in) + get_lprob_zinb(n - s, r_out, p_out, pi_out) - base_lprob);
                }

                if (pval > 1.0) {
                    common::logger->error("fail: {}", pval);
                    throw std::runtime_error("FOO");
                }

                return pval;
            };

            compute_pval = [r_in,p_in=p,pi_in,r_out,p_out=p,pi_out,get_lprob_zinb,get_lprob_ztnb,deviances1,deviances2](int64_t in_sum, int64_t out_sum, const auto &) {
                int64_t n = in_sum + out_sum;

                if (n == 0)
                    return 1.0;

                std::vector<double> devs(n + 1);

                for (int64_t s = 0; s <= n; ++s) {
                    devs[s] = deviances1[s] + deviances2[n - s];
                }

                double target_d = devs[in_sum];
                double pval = 0.0;
                double base_lprob = get_lprob_ztnb(n, r_in, r_out, p_in, pi_in, pi_out);
                int64_t s = 0;
                for ( ; s <= n; ++s) {
                    if (devs[s] >= target_d) {
                        pval += exp(get_lprob_zinb(s, r_in, p_in, pi_in) + get_lprob_zinb(n - s, r_out, p_out, pi_out) - base_lprob);
                    } else {
                        break;
                    }
                }

                for (int64_t ns = n; ns > s; --ns) {
                    if (devs[ns] >= target_d) {
                        pval += exp(get_lprob_zinb(ns, r_in, p_in, pi_in) + get_lprob_zinb(n - ns, r_out, p_out, pi_out) - base_lprob);
                    } else {
                        break;
                    }
                }

                if (pval > 1.0) {
                    common::logger->error("fail: {}", pval);
                    throw std::runtime_error("FOO");
                }

                return pval;
            };
        } else {
            double p = target_p;

            double p1p = exp(log(p) - log1p(-p));
            size_t total = nelem * groups.size();

            double mu = static_cast<double>(total_kmers) / total;
            double mu2 = static_cast<double>(total_kmers_sq) / total;

            if (mu >= p * (mu2 - mu*mu))
                common::logger->warn("Data is underdispersed: {} >= {}", mu, p * (mu2 - mu*mu));

            double ln_var = std::max(log(p * mu2 - mu) - log(p) - 2.0 * log(mu), 0.0);
            double ln_mu = log(mu) + log(p) - log1p(-p) - ln_var / 2;

            common::logger->trace("Lognormal MoM fit: mu: {}\tvar: {}\tE[X]: {}\tVar(X): {}\t1.0/E[X]: {}",
                                  ln_mu, ln_var, exp(ln_mu + ln_var/2), exp(ln_mu*2+ln_var)*(exp(ln_var)-1), exp(-ln_mu - ln_var/2));

            auto get_r = [ln_mu,ln_var,real_sum,p1p,total,lp=log(p),sum=static_cast<double>(total_kmers),gs=groups.size()](const PairContainer &row) {
                if (ln_var == 0)
                    return exp(-ln_mu);

                auto get_l = [&](double a) {
                    double r = 1.0 / a;
                    double l = (lp * r - log(a) - pow(log(a)-ln_mu, 2.0)/2/ln_var) * gs - lgamma(r) * row.size();
                    for (const auto &[j, c] : row) {
                        l += lgamma(r + c);
                    }
                    return l;
                };

                uintmax_t max_iter = 100;
                auto [a, l] = boost::math::tools::brent_find_minima(
                    [&](double a) { return -get_l(a); },
                    1.0 / real_sum, real_sum, 30, max_iter
                );
                l *= -1;

                return 1.0 / a;
            };

            auto get_deviance = [p](double invphi, double y) {
                double mu = invphi * (1.0 - p) / p;
                if (y == 0)
                    return 2.0*invphi*log((invphi*mu)/invphi);

                double logmuinvphi = log(mu+invphi);
                double logmu = log(mu);
                double logyshift = logmuinvphi - log(y + invphi);
                return 2.0 * (y * (log(y) - logmu + logyshift) + invphi * logyshift);
            };

            compute_min_pval = [get_r,num_labels_in,num_labels_out,get_deviance,gs=groups.size()](int64_t n, const auto &row) {
                if (row.empty())
                    return 1.0;

                double r = get_r(row);
                if (r == 0.0)
                    return 1.0;

                double r_in = r * num_labels_in;
                double r_out = r * num_labels_out;
                double lscaling_base = lgamma(r_in + r_out) - lgamma(r_in) - lgamma(r_out);

                double lscaling = lgamma(n + 1) - lgamma(r_in + r_out + n) + lscaling_base;

                std::vector<double> devs(n + 1);

                for (int64_t s = 0; s <= n; ++s) {
                    devs[s] = get_deviance(r_in, s) + get_deviance(r_out, n - s);
                }

                double max_d = *std::max_element(devs.begin(), devs.end());

                double pval = 0.0;
                int64_t s = 0;
                int64_t t = n;
                if (devs[s] == max_d) {
                    double rs = r_in;
                    double rt = r_out + n;
                    double base = lscaling - lgamma(n + 1) + lgamma(rs) + lgamma(rt);
                    pval += exp(base);
                    for (++s; s <= n; ++s) {
                        --t;
                        if (devs[s] == max_d) {
                            ++rs;
                            --rt;
                            base += log(t) - log(s) + log(rs - 1) - (rt > 1 ? log(rt - 1) : lgamma(rt + 1) - lgamma(rt));
                            pval += exp(base);
                        } else {
                            break;
                        }
                    }
                }

                int64_t sp = n;
                t = 0;
                if (devs[sp] == max_d) {
                    double rs = r_in + sp;
                    double rt = r_out;
                    double base = lscaling - lgamma(n + 1) + lgamma(rs) + lgamma(rt);
                    pval += exp(base);
                    for (--sp; sp >= s; --sp) {
                        ++t;
                        if (devs[sp] == max_d) {
                            --rs;
                            ++rt;
                            base += log(sp) - log(t) - (rs > 1 ? log(rs - 1) : lgamma(rs + 1) - lgamma(rs)) + log(rt - 1);
                            pval += exp(base);
                        } else {
                            break;
                        }
                    }
                }

                return std::min(1.0, pval);
            };

            compute_pval = [get_r,num_labels_in,num_labels_out,get_deviance,gs=groups.size()](int64_t in_sum, int64_t out_sum, const auto &row) {
                if (row.empty())
                    return 1.0;

                double r = get_r(row);

                if (r == 0.0)
                    return 1.0;

                double r_in = r * num_labels_in;
                double r_out = r * num_labels_out;
                double lscaling_base = lgamma(r_in + r_out) - lgamma(r_in) - lgamma(r_out);

                int64_t n = in_sum + out_sum;

                double lscaling = lgamma(n + 1) - lgamma(r_in + r_out + n) + lscaling_base;

                std::vector<double> devs(n + 1);

                for (int64_t s = 0; s <= n; ++s) {
                    devs[s] = get_deviance(r_in, s) + get_deviance(r_out, n - s);
                }

                double in_sum_total_deviance = devs[in_sum];

                double pval = 0.0;
                int64_t s = 0;
                int64_t t = n;
                if (devs[s] >= in_sum_total_deviance) {
                    double rs = r_in;
                    double rt = r_out + n;
                    double base = lscaling - lgamma(n + 1) + lgamma(rs) + lgamma(rt);
                    pval += exp(base);
                    for (++s; s <= n; ++s) {
                        --t;
                        if (devs[s] >= in_sum_total_deviance) {
                            ++rs;
                            --rt;
                            base += log(t) - log(s) + log(rs - 1) - (rt > 1 ? log(rt - 1) : lgamma(rt + 1) - lgamma(rt));
                            pval += exp(base);
                        } else {
                            break;
                        }
                    }
                }

                int64_t sp = n;
                t = 0;
                if (devs[sp] >= in_sum_total_deviance) {
                    double rs = r_in + sp;
                    double rt = r_out;
                    double base = lscaling - lgamma(n + 1) + lgamma(rs) + lgamma(rt);
                    pval += exp(base);
                    for (--sp; sp >= s; --sp) {
                        ++t;
                        if (devs[sp] >= in_sum_total_deviance) {
                            --rs;
                            ++rt;
                            base += log(sp) - log(t) - (rs > 1 ? log(rs - 1) : lgamma(rs + 1) - lgamma(rs)) + log(rt - 1);
                            pval += exp(base);
                        } else {
                            break;
                        }
                    }
                }

                return std::min(1.0, pval);
            };
        }
    } else if (config.test_type == "notest") {
        compute_min_pval = [&](int64_t n, const PairContainer&) { return n > 0 ? 0.0 : 1.1; };
        compute_pval = [&](int64_t, int64_t, const auto &row) { return row.size() ? 0.0 : 1.1; };
    } else {
        throw std::runtime_error("Test not implemented");
    }

    common::logger->trace("Computing minimum p-values");
    std::vector<std::pair<double, size_t>> m;
    std::vector<int64_t> m_sums;
    {
        std::vector<std::vector<std::pair<double, size_t>>> ms(num_threads + 1);
        std::vector<tsl::hopscotch_map<std::vector<int64_t>, std::pair<double, size_t>, utils::VectorHash>> ms_vecs(num_threads + 1);
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

            size_t n = in_sum + out_sum;
            if (!in_kmer && !out_kmer) {
                row.clear();
                n = 0;
            }

            if (config.test_type != "gnb_exact") {
                if (n >= ms[bucket_idx].size())
                    ms[bucket_idx].resize(n + 1, std::make_pair(1.1, 0));

                if (ms[bucket_idx][n].first == 1.1)
                    ms[bucket_idx][n].first = compute_min_pval(n, row);

                ++ms[bucket_idx][n].second;
            } else {
                std::sort(vals.begin(), vals.end());
                auto find = ms_vecs[bucket_idx].find(vals);
                if (find == ms_vecs[bucket_idx].end()) {
                    ms_vecs[bucket_idx][vals] = std::make_pair(compute_min_pval(n, row), 1);
                } else {
                    ++find.value().second;
                }
            }

            if (config.test_by_unitig) {
                bucket_idx = std::min(bucket_idx, pvals_buckets.size() - 1);
                uint64_t eff_size = bit_cast<uint64_t, int64_t>(in_sum - out_stat_int);
                if constexpr(preallocated) {
                    node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
                    eff_size_buckets[bucket_idx][node] = eff_size;
                    sum_buckets[bucket_idx][node] = n;
                } else {
                    eff_size_buckets[bucket_idx].push_back(eff_size);
                    sum_buckets[bucket_idx].push_back(n);
                }
            }

        }, false);

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
            for (size_t i = 1; i < ms_vecs.size(); ++i) {
                for (const auto &[vec, v] : ms_vecs[i]) {
                    auto find = ms_vecs[0].find(vec);
                    if (find == ms_vecs[0].end()) {
                        ms_vecs[0][vec] = v;
                    } else {
                        find.value().second += v.second;
                    }
                }
            }
            ms_vecs.resize(1);

            std::vector<std::tuple<double, int64_t, size_t>> mn;
            mn.reserve(ms_vecs[0].size());
            for (const auto &[vec, v] : ms_vecs[0]) {
                mn.emplace_back(
                    v.first,
                    std::accumulate(vec.begin(), vec.end(), int64_t(0)),
                    v.second
                );
            }
            std::sort(mn.begin(), mn.end(), utils::GreaterFirst());

            m.reserve(ms_vecs[0].size());
            m_sums.reserve(ms_vecs[0].size());
            for (const auto &[p, s, c] : mn) {
                m.emplace_back(p, c);
                m_sums.emplace_back(s);
            }

            ms_vecs.resize(0);
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
                    size_t n = sum_buckets[bucket_idx][path[i] - boundaries[bucket_idx]];
                    n_sum += n;
                    pvals_min.emplace_back(m[n].first);
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

            size_t n = in_sum + out_sum;
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
        }, false);

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
                groups[offset] = labels_out.count(label);
                columns_all[offset].reset(column.release());
                column_values_all[offset] = std::make_unique<ValuesContainer>(std::move(column_values));
            },
            num_parallel_files
        );

        uint8_t max_width = 0;
        for (const auto &col_vals : column_values_all) {
            max_width = std::max(max_width, col_vals->width());
        }

        auto generate_rows = [&](const auto &callback, bool ordered) {
            if (ordered) {
                utils::call_rows<std::unique_ptr<const bit_vector>,
                                 std::unique_ptr<const ValuesContainer>,
                                 PairContainer, true>(columns_all, column_values_all,
                                                      [&](uint64_t row_i, const auto &row, size_t thread_id) {
                    callback(row_i, row, thread_id);
                });
            } else {
                utils::call_rows<std::unique_ptr<const bit_vector>,
                                 std::unique_ptr<const ValuesContainer>,
                                 PairContainer, false>(columns_all, column_values_all,
                                                       [&](uint64_t row_i, const auto &row, size_t thread_id) {
                    callback(row_i, row, thread_id);
                });
            }
        };

        bool parallel = get_num_threads() > 1;

        return mask_nodes_by_label_dual<value_type, PValStorage>(
            graph_ptr,
            [&](const std::vector<size_t> &min_counts, sdsl::bit_vector *kept) -> std::vector<VectorMap<uint64_t, size_t>> {
                common::logger->trace("Calculating count histograms");
                std::vector<VectorMap<uint64_t, size_t>> hists_map(groups.size());
                utils::call_rows<std::unique_ptr<const bit_vector>,
                                 std::unique_ptr<const ValuesContainer>,
                                 PairContainer, false>(columns_all, column_values_all,
                                                       [&](uint64_t row_i, const auto &row, size_t) {
                    if (row.empty()) {
                        if (kept)
                            unset_bit(kept->data(), row_i, parallel, std::memory_order_relaxed);

                        return;
                    }

                    bool found = false;
                    Vector<uint64_t> counts(groups.size());
                    for (const auto &[j, raw_c] : row) {
                        if (min_counts.empty() || raw_c >= min_counts[j]) {
                            counts[j] = raw_c;
                            found = true;
                        }
                    }

                    if (!found) {
                        if (kept)
                            unset_bit(kept->data(), row_i, parallel, std::memory_order_relaxed);

                        return;
                    }

                    #pragma omp critical
                    {
                        for (size_t j = 0; j < counts.size(); ++j) {
                            ++hists_map[j][counts[j]];
                        }
                    }
                });

                return hists_map;
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
    common::logger->trace("Labels in: {}", fmt::join(labels_in, ","));
    common::logger->trace("Labels out: {}", fmt::join(labels_out, ","));

    size_t total_labels = labels_in.size() + labels_out.size();
    const auto &annotation = anno_graph.get_annotator();

    std::vector<bool> groups;
    groups.reserve(total_labels);
    for (const auto &label : annotation.get_label_encoder().get_labels()) {
        groups.emplace_back(labels_out.count(label));
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
            [&](const auto &callback, bool ordered) {
                int_matrix->call_row_values(callback, ordered);
            },
            groups,
            config, num_threads, tmp_dir, num_threads
        );
    } else {
        auto generate_rows = [&](const auto &callback, bool) {
            size_t bucket_size = matrix.num_rows() / num_threads;
            ProgressBar progress_bar(matrix.num_rows(), "Streaming rows", std::cerr, !common::get_verbose());
            #pragma omp parallel for num_threads(num_threads) ordered
            for (size_t row_i = 0; row_i < matrix.num_rows(); ++row_i) {
                ++progress_bar;
                std::vector<Row> row_ids(1, row_i);
                auto set_bits = matrix.get_rows(row_ids);
                Vector<std::pair<uint64_t, uint64_t>> container;
                container.reserve(set_bits[0].size());
                for (auto j : set_bits[0]) {
                    container.emplace_back(j, 1);
                }

                #pragma omp ordered
                callback(row_i, container, row_i / bucket_size);
            }
        };

        bool parallel = get_num_threads() > 1;

        return mask_nodes_by_label_dual<uint64_t, PValStorage>(
            std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr()),
            [&](const std::vector<size_t> &, sdsl::bit_vector *kept) -> std::vector<VectorMap<uint64_t, size_t>> {
                common::logger->trace("Calculating count histograms");
                std::vector<VectorMap<uint64_t, size_t>> hists_map(groups.size());
                for (size_t j = 0; j < groups.size(); ++j) {
                    hists_map[j][0] = 0;
                    hists_map[j][1] = 0;
                }
                ProgressBar progress_bar(matrix.num_rows(), "Streaming rows", std::cerr, !common::get_verbose());
                #pragma omp parallel for num_threads(num_threads)
                for (size_t row_i = 0; row_i < matrix.num_rows(); ++row_i) {
                    ++progress_bar;
                    std::vector<Row> row_ids(1, row_i);
                    auto set_bits = matrix.get_rows(row_ids);
                    std::vector<bool> container(groups.size());
                    for (auto j : set_bits[0]) {
                        container[j] = true;
                    }

                    if (kept && set_bits.empty())
                        unset_bit(kept->data(), row_i, parallel, std::memory_order_relaxed);

                    #pragma omp critical
                    {
                    for (size_t j = 0; j < container.size(); ++j) {
                        ++hists_map[j][container[j]];
                    }
                    }
                }

                return hists_map;
            },
            generate_rows,
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


