#include "annotated_graph_algorithm.hpp"

#include <typeinfo>
#include <mutex>

#include <sdust.h>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

#include "common/logger.hpp"
#include "common/vectors/bitmap.hpp"
#include "common/vector_map.hpp"
#include "common/vectors/transpose.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/graph_cleaning.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"

namespace mtg {
namespace graph {

using mtg::common::logger;

typedef AnnotatedDBG::node_index node_index;
typedef AnnotatedDBG::Annotator Annotator;
typedef AnnotatedDBG::Annotator::Label Label;
using Column = annot::matrix::BinaryMatrix::Column;
using PairContainer = std::vector<std::pair<uint64_t, uint64_t>>;
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
    std::vector<std::unique_ptr<const bit_vector>> columns_all(total_labels);
    std::vector<bool> groups(total_labels);

    if (config.test_by_unitig) {
        using ValuesContainer = sdsl::int_vector<>;
        std::vector<std::unique_ptr<const ValuesContainer>> column_values_all(total_labels);
        annot::ColumnCompressed<>::load_columns_and_values(
            files,
            [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> &&column, ValuesContainer&& column_values) {
                groups[offset] = labels_out.count(label);
                columns_all[offset].reset(column.release());
                column_values_all[offset] = std::make_unique<ValuesContainer>(std::move(column_values));
            },
            num_parallel_files
        );

        return mask_nodes_by_label_dual<ValuesContainer, PValStorage>(
            graph_ptr, columns_all, column_values_all, groups,
            config, num_threads, tmp_dir, num_parallel_files, true
        );
    } else {
        using ValuesContainer = sdsl::int_vector_buffer<>;
        std::vector<std::unique_ptr<const ValuesContainer>> column_values_all(total_labels);
        annot::ColumnCompressed<>::load_columns_and_values(
            files,
            [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> &&column, ValuesContainer&& column_values) {
                groups[offset] = labels_out.count(label);
                columns_all[offset].reset(column.release());
                column_values_all[offset] = std::make_unique<ValuesContainer>(std::move(column_values));
            },
            num_parallel_files
        );

        return mask_nodes_by_label_dual<ValuesContainer, PValStorage>(
            graph_ptr, columns_all, column_values_all, groups,
            config, num_threads, tmp_dir, num_parallel_files, true
        );
    }
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

template <class ValuesContainer, class PValStorage>
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, PValStorage, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                         std::vector<std::unique_ptr<const bit_vector>> &columns_all,
                         std::vector<std::unique_ptr<const ValuesContainer>> &column_values_all,
                         const std::vector<bool> &groups,
                         const DifferentialAssemblyConfig &config,
                         size_t num_threads,
                         std::filesystem::path tmp_dir,
                         size_t num_parallel_files,
                         bool deallocate) {
    num_parallel_files = std::min(num_parallel_files, num_threads);
    bool is_primary = graph_ptr->get_mode() == DeBruijnGraph::PRIMARY;
    bool parallel = num_parallel_files > 1;

    size_t num_labels_out = std::count(groups.begin(), groups.end(), true);
    size_t num_labels_in = groups.size() - num_labels_out;

    common::logger->trace("Graph mode: {}", is_primary ? "PRIMARY" : "other");

    uint8_t max_width = 0;
    for (const auto &col_vals : column_values_all) {
        max_width = std::max(max_width, col_vals->width());
    }

    sdsl::bit_vector kept(AnnotatedDBG::graph_to_anno_index(graph_ptr->max_index() + 1), true);

    std::vector<uint64_t> min_counts(groups.size(), config.min_count);
    std::vector<uint64_t> check_cutoff(groups.size(), std::numeric_limits<uint64_t>::max());

    if (config.clean && column_values_all.empty())
        common::logger->warn("Can't clean when no counts provided, skipping cleaning");

    if (config.clean && column_values_all.size()) {
        common::logger->trace("Cleaning count columns");

        #pragma omp parallel for num_threads(num_parallel_files)
        for (size_t j = 0; j < groups.size(); ++j) {
            const auto &column_values = *column_values_all[j];

            // set cutoff for lower end of distribution
            auto [mean_est, nzeros_est] = estimate_ztp_mean(
                [&](const auto &callback) {
                    for (size_t i = 0; i < column_values.size(); ++i) {
                        if (column_values[i] <= N_BUCKETS_FOR_ESTIMATION)
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
        }
    }

    common::logger->trace("Marking discarded k-mers");
    std::mutex agg_mu;
    std::vector<VectorMap<uint64_t, size_t>> hists_map(groups.size());

    if (groups.size()) {
        std::atomic_thread_fence(std::memory_order_release);
        utils::call_rows<std::unique_ptr<const bit_vector>,
                         std::unique_ptr<const ValuesContainer>,
                         PairContainer, false>(columns_all, column_values_all,
                                               [&](uint64_t row_i, const auto &row, size_t) {
            for (const auto &[j, raw_c] : row) {
                if (raw_c >= min_counts[j] && raw_c <= check_cutoff[j]) {
                    if (!config.test_by_unitig) {
                        std::lock_guard<std::mutex> lock(agg_mu);
                        sdsl::bit_vector marker(groups.size(), false);
                        for (const auto &[j, raw_c] : row) {
                            ++hists_map[j][raw_c];
                            marker[j] = true;
                        }

                        for (size_t j = 0; j < marker.size(); ++j) {
                            if (!marker[j])
                                ++hists_map[j][0];
                        }
                    }
                    return;
                }
            }

            unset_bit(kept.data(), row_i, parallel, MO_RELAXED);
        });
        std::atomic_thread_fence(std::memory_order_acquire);
    }

    size_t nelem = sdsl::util::cnt_one_bits(kept);

    sdsl::bit_vector unitig_start;
    if (config.test_by_unitig)
        unitig_start = sdsl::bit_vector(graph_ptr->max_index() + 1, true);

    std::shared_ptr<MaskedDeBruijnGraph> clean_masked_graph;
    std::atomic<size_t> num_unitigs(0);
    if (config.test_by_unitig) {
        clean_masked_graph = std::make_shared<MaskedDeBruijnGraph>(
            graph_ptr,
            [&](node_index node) {
                return node != DeBruijnGraph::npos
                        && kept[AnnotatedDBG::graph_to_anno_index(node)];
            },
            true,
            is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
        );

        std::atomic_thread_fence(std::memory_order_release);
        clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
            std::vector<uint64_t> unitig_sums(groups.size(), 0);
            std::vector<bool> all_below_cutoff(groups.size(), true);
            for (node_index node : path) {
                uint64_t row_i = AnnotatedDBG::graph_to_anno_index(node);
                for (size_t j = 0; j < groups.size(); ++j) {
                    const auto &col = *columns_all[j];
                    const auto &col_vals = *column_values_all[j];
                    if (uint64_t r = col.conditional_rank1(row_i)) {
                        uint64_t c = col_vals[r - 1];
                        unitig_sums[j] += c;
                        if (c >= min_counts[j])
                            all_below_cutoff[j] = false;
                    }
                }
            }

            unset_bit(unitig_start.data(), path[0], parallel, MO_RELAXED);
            num_unitigs.fetch_add(1, MO_RELAXED);

            std::lock_guard<std::mutex> lock(agg_mu);
            for (size_t j = 0; j < groups.size(); ++j) {
                ++hists_map[j][unitig_sums[j]];
            }
        }, num_threads);
        std::atomic_thread_fence(std::memory_order_acquire);

        nelem = num_unitigs.load();
    }

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
    std::unique_ptr<utils::TempFile> tmp_file;
    PValStorage pvals;
    if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
        tmp_file = std::make_unique<utils::TempFile>(tmp_dir);
        pvals = sdsl::int_vector_buffer<64>(tmp_file->name(), std::ios::out);
    }

    for (size_t i = 0; i <= graph_ptr->max_index(); ++i) {
        pvals.push_back(nullpval);
    }

    common::logger->trace("Allocating k-mer bitmasks");
    sdsl::bit_vector indicator_in(graph_ptr->max_index() + 1, false);
    sdsl::bit_vector indicator_out(graph_ptr->max_index() + 1, false);

    VectorMap<size_t, size_t> m;
    std::function<double(int64_t, int64_t, const PairContainer&)> compute_pval;
    std::function<double(int64_t)> compute_min_pval;

    std::vector<VectorMap<uint64_t, std::pair<size_t, uint64_t>>> count_maps;

    common::logger->trace("Test: {}", config.test_type);
    if (config.test_type == "poisson_exact") {
        compute_min_pval = [&](int64_t n) {
            if (n == 0)
                return 1.1;

            double p = static_cast<double>(in_kmers) / total_kmers;
            auto bdist = boost::math::binomial(n, p);
            auto get_deviance = [&](double y, double mu) {
                y += 1e-8;
                mu += 1e-8;
                return 2 * (y * log(y/mu) - y + mu);
            };

            double mu1 = static_cast<double>(in_kmers) / nelem;
            double mu2 = static_cast<double>(out_kmers) / nelem;

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

        compute_pval = [&](int64_t in_sum, int64_t out_sum, const auto &row) {
            if (row.empty())
                return 1.1;

            int64_t n = in_sum + out_sum;
            double p = static_cast<double>(in_kmers) / total_kmers;
            auto bdist = boost::math::binomial(n, p);
            auto get_deviance = [&](double y, double mu) {
                y += 1e-8;
                mu += 1e-8;
                return 2 * (y * log(y/mu) - y + mu);
            };

            double mu1 = static_cast<double>(in_kmers) / nelem;
            double mu2 = static_cast<double>(out_kmers) / nelem;

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

        #pragma omp parallel for num_threads(num_parallel_files)
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

            #pragma omp critical
            matrix_size += total + (hist.size() && hist[0].first == 0 ? hist[0].second : 0);
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

        auto get_deviance = [&](double y, double mu, double invphi) {
            y += 1e-8;
            mu += 1e-8;
            return 2.0 * (y * log( y/mu ) + (y + invphi) * log( (mu + invphi)/(y + invphi) ) );
        };

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

        compute_min_pval = [&,lscaling_base,r_in,r_out,target_p,get_deviance,mu1,mu2](int64_t n) {
            if (n == 0)
                return 1.1;

            double lscaling = lgamma(n + 1) + lscaling_base;
            auto get_pmf = [&](int64_t s) {
                int64_t t = n - s;
                double rs = r_in + s;
                double rt = r_out + t;
                return exp(lscaling - lgamma(s + 1) - lgamma(t + 1) + lgamma(rs) + lgamma(rt) - lgamma(rs + rt));
            };

            std::vector<double> devs;
            devs.reserve(n + 1);
            for (int64_t s = 0; s <= n; ++s) {
                devs.emplace_back(get_deviance(s, mu1, r_in) + get_deviance(n - s, mu2, r_out));
            }
            double max_d = *std::max_element(devs.begin(), devs.end());
            double pval = 0.0;
            int64_t s = 0;
            bool found = false;
            for ( ; s <= n; ++s) {
                if (devs[s] == max_d) {
                    pval += get_pmf(s);
                    found = true;
                } else {
                    break;
                }
            }

            for (int64_t sp = n; sp >= s; --sp) {
                if (devs[sp] == max_d) {
                    pval += get_pmf(sp);
                    found = true;
                } else {
                    break;
                }
            }

            if (!found) {
                common::logger->error("Devs: {}", fmt::join(devs, ","));
                throw std::runtime_error("Fail");
            }

            return std::min(1.0, pval);
        };

        compute_pval = [&,lscaling_base,r_in,r_out,target_p,get_deviance,mu1,mu2](int64_t in_sum, int64_t out_sum, const auto &row) {
            if (row.empty())
                return 1.1;

            int64_t n = in_sum + out_sum;

            double lscaling = lgamma(n + 1) + lscaling_base;
            auto get_pmf = [&](int64_t s) {
                int64_t t = n - s;
                double rs = r_in + s;
                double rt = r_out + t;
                return exp(lscaling - lgamma(s + 1) - lgamma(t + 1) + lgamma(rs) + lgamma(rt) - lgamma(rs + rt));
            };

            std::vector<double> devs;
            devs.reserve(n + 1);
            for (int64_t s = 0; s <= n; ++s) {
                devs.emplace_back(get_deviance(s, mu1, r_in) + get_deviance(n - s, mu2, r_out));
            }
            double pval = 0.0;
            int64_t s = 0;
            for ( ; s <= n; ++s) {
                if (devs[s] >= devs[in_sum]) {
                    pval += get_pmf(s);
                } else {
                    break;
                }
            }

            for (int64_t sp = n; sp >= s; --sp) {
                if (devs[sp] >= devs[in_sum]) {
                    pval += get_pmf(sp);
                } else {
                    break;
                }
            }

            return std::min(1.0, pval);
        };

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

        PairContainer best_row;
        int64_t in_sum = 0;
        int64_t out_sum = 0;
        for (size_t j = 0; j < groups.size(); ++j) {
            if (hists[j].size()) {
                uint64_t k = count_maps[j].find(hists[j].back().first)->second.first;
                best_row.emplace_back(j, k);
                if (groups[j]) {
                    out_sum += k;
                } else {
                    in_sum += k;
                }
            }
        }

        double pval = compute_pval(in_sum, out_sum, best_row);
        common::logger->trace("Best achievable p-value: {}", pval);
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

    common::logger->trace("Running differential tests");
    std::exception_ptr ex = nullptr;

    std::mutex pval_mu;
    std::atomic_thread_fence(std::memory_order_release);

    if (config.test_by_unitig) {
        clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
            if (ex)
                return;

            std::vector<int64_t> row_counts(groups.size());
            for (node_index node : path) {
                uint64_t row_i = AnnotatedDBG::graph_to_anno_index(node);
                for (size_t j = 0; j < groups.size(); ++j) {
                    const auto &col = *columns_all[j];
                    const auto &col_vals = *column_values_all[j];
                    if (uint64_t r = col.conditional_rank1(row_i))
                        row_counts[j] += col_vals[r - 1];
                }
            }

            PairContainer merged_row;
            int64_t in_sum = 0;
            int64_t out_sum = 0;
            size_t count_in = 0;
            size_t count_out = 0;
            for (size_t j = 0; j < groups.size(); ++j) {
                if (row_counts[j] > 0) {
                    if (count_maps.size())
                        row_counts[j] = count_maps[j].find(row_counts[j])->second.first;

                    merged_row.emplace_back(j, row_counts[j]);
                    if (groups[j]) {
                        out_sum += row_counts[j];
                        ++count_out;
                    } else {
                        in_sum += row_counts[j];
                        ++count_in;
                    }
                }
            }

            if (count_in + count_out < config.min_recurrence
                    || count_in < config.min_in_recurrence
                    || count_out < config.min_out_recurrence
                    || count_in > config.max_in_recurrence
                    || count_out > config.max_out_recurrence) {
                merged_row.clear();
                in_sum = 0;
                out_sum = 0;
            }

            double pval;
            double pval_min;
            try {
                pval_min = compute_min_pval(in_sum + out_sum);
                pval = compute_pval(in_sum, out_sum, merged_row);
            } catch (...) {
                std::lock_guard<std::mutex> lock(pval_mu);
                ex = std::current_exception();
                return;
            }

            if (pval >= 1.1)
                return;

            double in_stat = in_sum;
            double out_stat = static_cast<double>(out_sum) / out_kmers * in_kmers;
            bool in_kmer = (in_stat > out_stat) || (out_stat != out_stat && in_stat == in_stat);
            bool out_kmer = (in_stat < out_stat) || (in_stat != in_stat && out_stat == out_stat);
            uint64_t k = std::numeric_limits<uint64_t>::max();

            if (pval_min - pval > 1e-10) {
                common::logger->error("Min p-val estimate too high: min {} > cur {}\ttest: {}", pval_min, pval, config.test_type);
                throw std::runtime_error("Test failed");
            }

            if (pval_min >= config.family_wise_error_rate) {
                k = 0;
            } else if (pval_min > 0) {
                double lkd = log2(config.family_wise_error_rate) - log2(pval_min);
                if (lkd <= 64)
                    k = pow(2.0, lkd);

                if (k == 0) {
                    common::logger->error("k: {}\tlog k: {}\tpval_min: {}\tpval: {}", k, lkd, pval_min, pval);
                    throw std::runtime_error("Min failed");
                }
            }

            if (in_kmer != out_kmer && pval < config.family_wise_error_rate) {
                for (size_t i = 0; i < path.size(); ++i) {
                    set_bit((in_kmer ? indicator_in : indicator_out).data(), path[i], parallel, MO_RELAXED);
                }
            }

            if (config.test_type != "notest") {
                std::lock_guard<std::mutex> lock(pval_mu);
                ++m[k];
                for (size_t i = 0; i < path.size(); ++i) {
                    pvals[path[i]] = bit_cast<uint64_t>(pval);
                }
            }
        }, num_threads);
    } else {
        utils::call_rows<std::unique_ptr<const bit_vector>,
                         std::unique_ptr<const ValuesContainer>,
                         PairContainer, true>(columns_all, column_values_all,
                                              [&](uint64_t row_i, const auto &raw_row, size_t batch_id) {
            if (ex)
                return;

            int64_t in_sum = 0;
            int64_t out_sum = 0;
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

                if (count_in + count_out < config.min_recurrence
                        || count_in < config.min_in_recurrence
                        || count_out < config.min_out_recurrence
                        || count_in > config.max_in_recurrence
                        || count_out > config.max_out_recurrence) {
                    row.clear();
                    in_sum = 0;
                    out_sum = 0;
                }
            }

            node_index node = AnnotatedDBG::anno_to_graph_index(row_i);

            double pval;
            double pval_min;
            try {
                pval_min = compute_min_pval(in_sum + out_sum);

                if (pval_min >= 1.1)
                    return;

                pval = compute_pval(in_sum, out_sum, row);
            } catch (...) {
                std::lock_guard<std::mutex> lock(pval_mu);
                ex = std::current_exception();
                return;
            }

            if (pval >= 1.1)
                return;

            double in_stat = in_sum;
            double out_stat = static_cast<double>(out_sum) / out_kmers * in_kmers;
            bool in_kmer = (in_stat > out_stat) || (out_stat != out_stat && in_stat == in_stat);
            bool out_kmer = (in_stat < out_stat) || (in_stat != in_stat && out_stat == out_stat);
            uint64_t k = std::numeric_limits<uint64_t>::max();

            if (pval_min - pval > 1e-10) {
                common::logger->error("Min p-val estimate too high: min {} > cur {}\ttest: {}", pval_min, pval, config.test_type);
                throw std::runtime_error("Test failed");
            }

            if (pval_min >= config.family_wise_error_rate) {
                k = 0;
            } else if (pval_min > 0) {
                double lkd = log2(config.family_wise_error_rate) - log2(pval_min);
                if (lkd <= 64)
                    k = pow(2.0, lkd);

                if (k == 0) {
                    common::logger->error("k: {}\tlog k: {}\tpval_min: {}\tpval: {}", k, lkd, pval_min, pval);
                    throw std::runtime_error("Min failed");
                }
            }

            if (in_kmer != out_kmer && pval < config.family_wise_error_rate)
                set_bit((in_kmer ? indicator_in : indicator_out).data(), node, parallel, MO_RELAXED);

            if (config.test_type != "notest") {
                std::lock_guard<std::mutex> lock(pval_mu);
                ++m[k];
                pvals[node] = bit_cast<uint64_t>(pval);
            }
        });
    }

    std::atomic_thread_fence(std::memory_order_acquire);

    if (deallocate) {
        for (auto &col : columns_all) {
            col.reset();
        }

        for (auto &col_vals : column_values_all) {
            col_vals.reset();
        }
    }

    if (ex)
        std::rethrow_exception(ex);

    if (config.test_type != "notest") {
        uint64_t total_sig = sdsl::util::cnt_one_bits(indicator_in) + sdsl::util::cnt_one_bits(indicator_out);
        auto m_data = const_cast<std::vector<std::pair<size_t, size_t>>&&>(m.values_container());

        std::sort(m_data.begin(), m_data.end(), utils::LessFirst());
        size_t acc = std::accumulate(m_data.begin(), m_data.end(), size_t(0),
                                     [](size_t sum, const auto &a) { return sum + a.second; });
        common::logger->trace("Performed {}/{} tests", acc, nelem);
        common::logger->trace("Correcting {}/{} significant p-values", total_sig, acc);

        if (total_sig) {
            size_t k = 0;
            size_t last_k = 0;
            for (const auto &[cur_k, s] : m_data) {
                common::logger->trace("k: {}\tm(k): {}", cur_k, acc);
                if (acc <= cur_k) {
                    k = std::max(acc, last_k + 1);
                    break;
                }

                last_k = cur_k;
                acc -= s;
            }

            if (k == 0) {
                common::logger->trace("No significant k-mers found");
                sdsl::util::set_to_value(indicator_in, false);
                sdsl::util::set_to_value(indicator_out, false);
            } else {
                common::logger->trace("k: {}\talpha_corr: {}", k, config.family_wise_error_rate / k);
                for (size_t i = 0; i < indicator_in.size(); ++i) {
                    double p = bit_cast<double, uint64_t>(pvals[i]);
                    if (p < config.family_wise_error_rate && p * k >= config.family_wise_error_rate) {
                        indicator_in[i] = false;
                        indicator_out[i] = false;
                    }
                }
            }
        }
    }

    call_ones(unitig_start, [&](node_index node) { pvals[node] = nullpval; });

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

template
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, std::vector<uint64_t>, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<sdsl::int_vector_buffer<>, std::vector<uint64_t>>(std::shared_ptr<const DeBruijnGraph>,
                         std::vector<std::unique_ptr<const bit_vector>> &,
                         std::vector<std::unique_ptr<const sdsl::int_vector_buffer<>>> &,
                         const std::vector<bool> &,
                         const DifferentialAssemblyConfig &,
                         size_t,
                         std::filesystem::path,
                         size_t,
                         bool);
template
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, std::vector<uint64_t>, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<sdsl::int_vector<>, std::vector<uint64_t>>(std::shared_ptr<const DeBruijnGraph>,
                         std::vector<std::unique_ptr<const bit_vector>> &,
                         std::vector<std::unique_ptr<const sdsl::int_vector<>>> &,
                         const std::vector<bool> &,
                         const DifferentialAssemblyConfig &,
                         size_t,
                         std::filesystem::path,
                         size_t,
                         bool);

template
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, sdsl::int_vector_buffer<64>, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<sdsl::int_vector_buffer<>, sdsl::int_vector_buffer<64>>(std::shared_ptr<const DeBruijnGraph>,
                         std::vector<std::unique_ptr<const bit_vector>> &,
                         std::vector<std::unique_ptr<const sdsl::int_vector_buffer<>>> &,
                         const std::vector<bool> &,
                         const DifferentialAssemblyConfig &,
                         size_t,
                         std::filesystem::path,
                         size_t,
                         bool);
template
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, sdsl::int_vector_buffer<64>, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<sdsl::int_vector<>, sdsl::int_vector_buffer<64>>(std::shared_ptr<const DeBruijnGraph>,
                         std::vector<std::unique_ptr<const bit_vector>> &,
                         std::vector<std::unique_ptr<const sdsl::int_vector<>>> &,
                         const std::vector<bool> &,
                         const DifferentialAssemblyConfig &,
                         size_t,
                         std::filesystem::path,
                         size_t,
                         bool);

} // namespace graph
} // namespace mtg


