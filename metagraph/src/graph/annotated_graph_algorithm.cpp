#include "annotated_graph_algorithm.hpp"

#include <typeinfo>
#include <mutex>
#include <variant>

#include <sdust.h>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/distributions/cauchy.hpp>

#include "common/logger.hpp"
#include "common/vectors/bitmap.hpp"
#include "common/vector_map.hpp"
#include "common/vectors/transpose.hpp"
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
    std::vector<double> medians(groups.size());
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
    if constexpr(std::is_same_v<PValStorage, std::vector<uint64_t>>) {
        auto &pvals = pvals_buckets.emplace_back();
        pvals.emplace_back(nullpval);
    }

    if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
        for (size_t i = 0; i < get_num_threads() + 1; ++i) {
            auto &tmp_file = tmp_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
            pvals_buckets.emplace_back(tmp_file->name(), std::ios::out);
        }
        pvals_buckets[0].push_back(nullpval);
    }

    std::unique_ptr<PValStorage> eff_size;
    std::unique_ptr<PValStorage> pvals_min;
    std::unique_ptr<utils::TempFile> tmp_file_eff;
    std::unique_ptr<utils::TempFile> tmp_file_min;
    if (config.test_by_unitig) {
        common::logger->trace("Allocating p-value storage");
        pvals_min = std::make_unique<PValStorage>();
        eff_size = std::make_unique<PValStorage>();
        if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
            tmp_file_min = std::make_unique<utils::TempFile>(tmp_dir);
            pvals_min = std::make_unique<sdsl::int_vector_buffer<64>>(tmp_file_min->name(), std::ios::out);

            tmp_file_eff = std::make_unique<utils::TempFile>(tmp_dir);
            eff_size = std::make_unique<sdsl::int_vector_buffer<64>>(tmp_file_eff->name(), std::ios::out);
        }

        for (size_t i = 0; i <= graph_ptr->max_index(); ++i) {
            pvals_min->push_back(nullpval);
            eff_size->push_back(0);
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
    std::function<double(int64_t)> compute_min_pval;

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
                              get_deviance](int64_t n) {
            if (n == 0)
                return 1.1;

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
                return 1.1;

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

        compute_min_pval = [&,lbase_shared](int64_t m1) {
            if (m1 == 0)
                return 1.1;

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
                return 1.1;

            int64_t m1 = in_sum + out_sum;
            int64_t m2 = total_kmers - m1;

            double lbase = lbase_shared + lgamma(m1 + 1) + lgamma(m2 + 1);

            int64_t a = in_sum;
            int64_t b = in_kmers - in_sum;
            int64_t c = out_sum;
            int64_t d = out_kmers - out_sum;

            return exp(lbase - lgamma(a + 1) - lgamma(b + 1) - lgamma(c + 1) - lgamma(d + 1));
        };
    } else if (config.test_type == "nbinom_exact") {
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
                r = boost::math::tools::newton_raphson_iterate([&](double r) {
                    double dl = (log(r) - log(r + mu) - boost::math::digamma(r)) * total;
                    double ddl = (1.0 / r - 1.0 / (r + mu) - boost::math::trigamma(r)) * total;
                    std::for_each(begin, end, [&](const auto &a) {
                        const auto &[k, c] = a;
                        dl += boost::math::digamma(k + r) * c;
                        ddl += boost::math::trigamma(k + r) * c;
                    });
                    return std::make_pair(dl, ddl);
                }, r_guess, std::numeric_limits<double>::min(), mu * total, 30);
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

        double r_map = boost::math::tools::newton_raphson_iterate([&](double r) {
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
        }, r_guess, std::numeric_limits<double>::min(), real_sum, 30);
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
                        count_maps[j][k] = std::make_pair(std::max(new_k, uint64_t(1)), c);
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
        common::logger->trace("Scaled Totals: in: {}\tout: {}", in_kmers, out_kmers);

        PairContainer best_in_row;
        PairContainer best_out_row;
        int64_t in_sum = 0;
        int64_t out_sum = 0;
        for (size_t j = 0; j < groups.size(); ++j) {
            if (hists[j].size()) {
                uint64_t k = count_maps[j].find(hists[j].back().first)->second.first;
                if (groups[j]) {
                    best_out_row.emplace_back(j, k);
                    out_sum += k;
                } else {
                    best_in_row.emplace_back(j, k);
                    in_sum += k;
                }
            }
        }
        int64_t total_sum = in_sum + out_sum;

        auto get_deviance1 = [mu=mu1+1e-8,invphi=r_in,logmuinvphi=log(mu1+1e-8+r_in),logmu=log(mu1+1e-8)](double y) {
            y += 1e-8;
            double logyshift = logmuinvphi - log(y + invphi);
            return 2.0 * (y * (log(y) - logmu + logyshift) + invphi * logyshift);
        };
        auto get_deviance2 = [mu=mu2+1e-8,invphi=r_out,logmuinvphi=log(mu2+1e-8+r_out),logmu=log(mu2+1e-8)](double y) {
            y += 1e-8;
            double logyshift = logmuinvphi - log(y + invphi);
            return 2.0 * (y * (log(y) - logmu + logyshift) + invphi * logyshift);
        };

        std::vector<double> deviances1(total_sum + 1);
        std::vector<double> deviances2(total_sum + 1);
        for (int64_t s = 0; s <= total_sum; ++s) {
            deviances1[s] = get_deviance1(s);
            deviances2[s] = get_deviance2(s);
        }

        compute_min_pval = [lscaling_base,r_in,r_out,target_p,deviances1,deviances2](int64_t n) {
            if (n == 0)
                return 1.1;

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
                return 1.1;

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

        double pval = std::min(compute_pval(in_sum, 0, best_in_row),
                               compute_pval(0, out_sum, best_out_row));
        common::logger->trace("Best achievable p-value: {}\tin_sum: {}\tout_sum: {}", pval, in_sum, out_sum);
        if (pval >= config.family_wise_error_rate) {
            common::logger->error("Best achievable p-value is too big. Use poisson_exact instead");
            throw std::domain_error("Too few samples");
        }
    } else if (config.test_type == "notest") {
        compute_min_pval = [&](int64_t n) { return n > 0 ? 0.0 : 1.1; };
        compute_pval = [&](int64_t, int64_t, const auto &row) { return row.size() ? 0.0 : 1.1; };
    } else {
        throw std::runtime_error("Test not implemented");
    }

    common::logger->trace("Computing minimum p-values");
    std::vector<std::pair<double, size_t>> m;
    {
        std::vector<std::vector<std::pair<double, size_t>>> ms(pvals_buckets.size());
        generate_rows([&](uint64_t row_i, const auto &raw_row, size_t bucket_idx) {
            int64_t in_sum = 0;
            int64_t out_sum = 0;
            int64_t out_stat_int = 0;
            bool in_kmer = false;
            bool out_kmer = false;
            PairContainer row;
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
                    if (out_stat != 0)
                        out_stat_int = out_stat > 0 ? ceil(out_stat) : floor(out_stat);

                    in_kmer = count_in >= config.min_in_recurrence && count_out <= config.max_out_recurrence && in_sum > (out_kmers > 0 ? out_stat : 0.0);
                    out_kmer = count_out >= config.min_out_recurrence && count_in <= config.max_in_recurrence && in_sum < out_stat;
                }
            }

            bucket_idx = std::min(bucket_idx, pvals_buckets.size() - 1);
            size_t n = in_sum + out_sum;
            if (!in_kmer && !out_kmer) {
                row.clear();
                n = 0;
            }

            if (n >= ms[bucket_idx].size())
                ms[bucket_idx].resize(n + 1, std::make_pair(1.1, 0));

            if (n > 0 && ms[bucket_idx][n].first == 1.1)
                ms[bucket_idx][n].first = compute_min_pval(n);

            ++ms[bucket_idx][n].second;

            // test_by_unitig
            std::ignore = out_stat_int;
            // (*eff_size)[node] = bit_cast<uint64_t, int64_t>(in_sum - out_stat_int);
            // (*pvals_min)[node] = bit_cast<uint64_t>(pval_min);
        }, false);

        common::logger->trace("Merging min p-value tables");
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
    }

    size_t k = 0;
    size_t n_cutoff = 1;
    if (!config.test_by_unitig) {
        size_t acc = std::accumulate(m.begin(), m.end(), size_t(0),
                                     [](size_t sum, const auto &a) { return sum + a.second; });
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
            if (s > 0)
                common::logger->trace("k: {}\ts: {}\tm(k): {}\tn: {}\tpval_min: {}", cur_k, s, acc, n, pval_min);

            if (cur_k != last_k && acc <= cur_k) {
                k = std::max(acc, last_k + 1);
                n_cutoff = n;
                common::logger->trace("Picked: k: {}\tn: {}\tpval_min: {}", k, n_cutoff, config.family_wise_error_rate / k);
                break;
            }

            last_k = cur_k;
            acc -= s;
        }

        if (k == 0)
            common::logger->trace("No significant k-mers found");
    }

    common::logger->trace("Running differential tests");
    std::exception_ptr ex = nullptr;
    std::mutex pval_mu;
    std::atomic_thread_fence(std::memory_order_release);

    generate_rows([&](uint64_t row_i, const auto &raw_row, size_t bucket_idx) {
        if (ex)
            return;

        int64_t in_sum = 0;
        int64_t out_sum = 0;
        bool in_kmer = false;
        bool out_kmer = false;
        PairContainer row;
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
        }

        size_t n = in_sum + out_sum;
        if (!in_kmer && !out_kmer) {
            row.clear();
            n = 0;
        }

        if (config.output_pvals || config.test_by_unitig || n >= n_cutoff) {
            if (m[n].first == 1.1) {
                common::logger->error("Min p-val too small: {}\tn: {}\ts: {}", m[n].first, n, m[n].second);
                throw std::runtime_error("Test failed");
            }

            double pval = 1.1;
            try {
                pval = compute_pval(in_sum, out_sum, row);
            } catch (...) {
                std::lock_guard<std::mutex> lock(pval_mu);
                ex = std::current_exception();
                return;
            }

            if (pval < 1.1 && m[n].first - pval > 1e-10) {
                common::logger->error("Min p-val estimate too high: min {} > cur {}\tn: {}\ttest: {}", m[n].first, pval, n, config.test_type);
                throw std::runtime_error("Test failed");
            }

            if (config.test_type != "notest" && (config.output_pvals || config.test_by_unitig)) {
                bucket_idx = std::min(bucket_idx, pvals_buckets.size() - 1);
                auto &pvals = pvals_buckets[bucket_idx];
                pvals.push_back(bit_cast<uint64_t>(pval));
            }

            if (!config.test_by_unitig && in_kmer != out_kmer && pval * k < config.family_wise_error_rate) {
                node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
                bool use_atomic = parallel && (node % 64 == 0);
                set_bit((in_kmer ? indicator_in : indicator_out).data(), node, use_atomic, MO_RELAXED);
            }
        }
    }, false);

    std::atomic_thread_fence(std::memory_order_acquire);

    if (ex)
        std::rethrow_exception(ex);

    deallocate();

    // std::vector<size_t> boundaries;
    // boundaries.reserve(pvals_buckets.size() + 1);
    // boundaries.emplace_back(0);
    // for (size_t i = 0; i < pvals_buckets.size(); ++i) {
    //     boundaries.emplace_back(boundaries.back() + pvals_buckets[i].size());
    // }

    // if (config.test_by_unitig) {
    //     common::logger->trace("Allocating k-mer bitmasks");
    //     indicator_in = sdsl::bit_vector(graph_ptr->max_index() + 1, false);
    //     indicator_out = sdsl::bit_vector(graph_ptr->max_index() + 1, false);

    //     common::logger->trace("Combining p-values within unitigs");
    //     MaskedDeBruijnGraph clean_masked_graph(
    //         graph_ptr,
    //         [&](node_index node) {
    //             return node != DeBruijnGraph::npos
    //                     && kept[AnnotatedDBG::graph_to_anno_index(node)];
    //         },
    //         true,
    //         is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
    //     );

    //     auto dist = boost::math::cauchy();

    //     std::atomic_thread_fence(std::memory_order_release);
    //     clean_masked_graph.call_unitigs([&](const std::string&, const auto &path) {
    //         double stat = 0;
    //         double stat_min = 0;
    //         int64_t comb_eff_size = 0;
    //         std::vector<size_t> bucket_idxs;
    //         bucket_idxs.reserve(path.size());
    //         for (node_index node : path) {
    //             for (size_t i = 0; i < boundaries.size(); ++i) {
    //                 if (boundaries[i] > node) {
    //                     bucket_idxs.emplace_back(i);
    //                     break;
    //                 }
    //             }
    //         }

    //         {
    //             std::lock_guard<std::mutex> lock(pval_mu);
    //             for (size_t i = 0; i < path.size(); ++i) {
    //                 node_index node = path[i];
    //                 size_t bucket_idx = bucket_idxs[i];
    //                 comb_eff_size += bit_cast<int64_t, uint64_t>((*eff_size)[node]);
    //                 double pval = bit_cast<double, uint64_t>(pvals_buckets[bucket_idx][node - boundaries[bucket_idx]]);
    //                 double pval_min = bit_cast<double, uint64_t>((*pvals_min)[node]);
    //                 stat += tan((0.5 - pval) * M_PI);
    //                 stat_min += tan((0.5 - pval_min) * M_PI);
    //             }
    //         }
    //         stat /= path.size();
    //         stat_min /= path.size();
    //         double comb_pval = boost::math::cdf(boost::math::complement(dist, stat));
    //         if (comb_pval < config.family_wise_error_rate && comb_eff_size != 0) {
    //             auto *data = comb_eff_size > 0 ? indicator_in.data() : indicator_out.data();
    //             for (node_index node : path) {
    //                 set_bit(data, node, parallel, MO_RELAXED);
    //             }
    //         }

    //         double comb_pval_min = boost::math::cdf(boost::math::complement(dist, stat_min));

    //         uint64_t k = std::numeric_limits<uint64_t>::max();
    //         if (comb_pval_min >= config.family_wise_error_rate) {
    //             k = 0;
    //         } else if (comb_pval_min > 0) {
    //             double lkd = log2(config.family_wise_error_rate) - log2(comb_pval_min);
    //             if (lkd <= 64)
    //                 k = pow(2.0, lkd);

    //             if (k == 0) {
    //                 common::logger->error("k: {}\tlog k: {}\tpval_min: {}\tpval: {}", k, lkd, comb_pval_min, comb_pval);
    //                 throw std::runtime_error("Min failed");
    //             }
    //         }

    //         uint64_t comb_pval_enc = bit_cast<uint64_t, double>(comb_pval);
    //         size_t bucket_idx = bucket_idxs[0];
    //         ++ms[bucket_idx][k];
    //         pvals_buckets[bucket_idx][path[0] - boundaries[bucket_idx]] = comb_pval_enc;
    //         if (config.output_pvals) {
    //             for (size_t i = 1; i < path.size(); ++i) {
    //                 size_t bucket_idx = bucket_idxs[i];
    //                 pvals_buckets[bucket_idx][path[i] - boundaries[bucket_idx]] = nullpval;
    //             }
    //         }
    //     }, num_threads);
    //     std::atomic_thread_fence(std::memory_order_acquire);

    //     eff_size.reset();
    //     pvals_min.reset();
    // }

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


