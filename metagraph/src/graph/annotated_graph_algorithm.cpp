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
    std::vector<VectorMap<uint64_t, size_t>> hists(groups.size());

    if (groups.size()) {
        std::atomic_thread_fence(std::memory_order_release);
        utils::call_rows<std::unique_ptr<const bit_vector>,
                         std::unique_ptr<const ValuesContainer>,
                         PairContainer, false>(columns_all, column_values_all,
                                               [&](uint64_t row_i, const auto &row, size_t) {
            for (const auto &[j, raw_c] : row) {
                if (raw_c >= min_counts[j] && raw_c <= check_cutoff[j]) {
                    std::lock_guard<std::mutex> lock(agg_mu);
                    for (const auto &[j, raw_c] : row) {
                        ++hists[j][raw_c];
                    }
                    return;
                }
            }

            unset_bit(kept.data(), row_i, parallel, MO_RELAXED);
        });
        std::atomic_thread_fence(std::memory_order_acquire);
    }

    sdsl::bit_vector unitig_start;
    if (config.test_by_unitig)
        unitig_start = sdsl::bit_vector(graph_ptr->max_index() + 1, true);

    std::shared_ptr<MaskedDeBruijnGraph> clean_masked_graph;
    std::vector<VectorMap<uint64_t, size_t>> unitig_hists(groups.size());
    size_t num_unitigs = 0;
    std::vector<uint64_t> sums_of_squares_unitigs(groups.size());
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

        clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
            unset_bit(unitig_start.data(), path[0], parallel, MO_RELAXED);
            std::vector<uint64_t> unitig_sums(groups.size());
            for (node_index node : path) {
                uint64_t row_i = AnnotatedDBG::graph_to_anno_index(node);
                for (size_t j = 0; j < groups.size(); ++j) {
                    const auto &col = *columns_all[j];
                    const auto &col_vals = *column_values_all[j];
                    if (uint64_t r = col.conditional_rank1(row_i))
                        unitig_sums[j] += col_vals[r - 1];
                }
            }

            for (size_t j = 0; j < groups.size(); ++j) {
                sums_of_squares_unitigs[j] += unitig_sums[j] * unitig_sums[j];
            }

            std::lock_guard<std::mutex> lock(agg_mu);
            ++num_unitigs;
            for (size_t j = 0; j < groups.size(); ++j) {
                ++unitig_hists[j][unitig_sums[j]];
            }
        }, num_threads);
    }

    double nelem = !config.test_by_unitig ? sdsl::util::cnt_one_bits(kept) : num_unitigs;

    common::logger->trace("Computing aggregate statistics");
    double in_kmers = 0;
    double out_kmers = 0;
    std::vector<double> medians(groups.size());
    std::vector<uint64_t> max_obs_vals(groups.size());
    std::vector<size_t> n_kmers(groups.size());
    std::vector<uint64_t> sums(groups.size());
    std::vector<uint64_t> sums_of_squares(groups.size());

    for (size_t j = 0; j < groups.size(); ++j) {
        for (const auto &[k, c] : hists[j]) {
            sums[j] += k * c;
            sums_of_squares[j] += k * k * c;
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

    double total_kmers = in_kmers + out_kmers;

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
    using RowStats = std::tuple<double, double, double, double>;
    std::function<RowStats(const PairContainer&)> compute_pval;
    std::function<RowStats(const std::vector<PairContainer>&)> compute_pval_unitig;

    common::logger->trace("Test: {}", config.test_type);
    if (config.test_type == "poisson_exact") {
        compute_pval = [&](const auto &row) {
            if (row.empty())
                return std::make_tuple(1.1, 0.0, 0.0, 1.1);

            int64_t in_sum = 0;
            int64_t out_sum = 0;

            for (const auto &[j, c] : row) {
                if (groups[j]) {
                    out_sum += c;
                } else {
                    in_sum += c;
                }
            }

            int64_t n = in_sum + out_sum;
            double pval = 1.0;
            double p_min = 1.0;

            double p = static_cast<double>(in_kmers) / total_kmers;
            auto bdist = boost::math::binomial(n, p);

            if (num_labels_in == num_labels_out) {
                p_min = boost::math::pdf(bdist, 0) + boost::math::pdf(bdist, n);
                if (in_sum != out_sum) {
                    int64_t d = std::min(in_sum, out_sum);
                    pval = boost::math::cdf(bdist, d) + boost::math::cdf(boost::math::complement(bdist, n - d - 1));
                }
            } else {
                p_min = std::min(boost::math::pdf(bdist, 0), boost::math::pdf(bdist, n));
                double d = abs(static_cast<double>(in_sum) / num_labels_in - static_cast<double>(out_sum) / num_labels_out);
                if (d > 0) {
                    // we want s s.t. | s/num_labels_in - (n-s)/num_labels_out | >= d
                    // s/num_labels_in - (n-s)/num_labels_out >= d || (n-s)/num_labels_out - s/num_labels_in >= d
                    // s/num_labels_in + s/num_labels_out >= d + n/num_labels_out || n/num_labels_out - d >= s/num_labels_in + s/num_labels_out
                    double factor = double(1.0) / num_labels_in + double(1.0) / num_labels_out;

                    double lower_d = (static_cast<double>(n) / num_labels_out - d) / factor;
                    double upper_d = (d + static_cast<double>(n) / num_labels_out) / factor;

                    // fix floating point errors
                    if (abs(lower_d - round(lower_d)) < 1e-10)
                        lower_d = round(lower_d);

                    if (abs(upper_d - round(upper_d)) < 1e-10)
                        upper_d = round(upper_d);

                    int64_t lower = floor(lower_d);
                    int64_t upper = ceil(upper_d);

                    if (lower < n && upper > 0 && lower != upper) {
                        pval = 0.0;

                        if (lower >= 0)
                            pval += boost::math::cdf(bdist, lower);

                        if (upper <= n)
                            pval += boost::math::cdf(boost::math::complement(bdist, upper - 1));
                    }
                }
            }

            if (pval > 1.0) {
                common::logger->error("{} > 1.0\t{}\t{},{}", pval, p_min, in_sum, out_sum);
                throw std::runtime_error("pval fail");
            }

            return std::make_tuple(pval,
                                   static_cast<double>(in_sum),
                                   static_cast<double>(out_sum) / out_kmers * in_kmers,
                                   p_min);
        };
    } else if (config.test_type == "fisher") {
        compute_pval = [&](const auto &row) {
            if (row.empty())
                return std::make_tuple(1.1, 0.0, 0.0, 1.1);

            int64_t in_sum = 0;
            int64_t out_sum = 0;

            for (const auto &[j, c] : row) {
                if (groups[j]) {
                    out_sum += c;
                } else {
                    in_sum += c;
                }
            }

            int64_t t = total_kmers;
            int64_t n1 = in_kmers;
            int64_t n2 = out_kmers;
            int64_t m1 = in_sum + out_sum;
            int64_t m2 = t - m1;

            double lbase = lgamma(n1 + 1) + lgamma(n2 + 1) + lgamma(m1 + 1) + lgamma(m2 + 1) - lgamma(t + 1);

            auto get_pval = [&](int64_t a) {
                int64_t b = n1 - a;
                int64_t c = m1 - a;
                int64_t d = m2 - b;
                assert(d == n2 - c);

                return b < 0 || c < 0 || d < 0
                    ? 1.0
                    : exp(lbase - lgamma(a + 1) - lgamma(b + 1) - lgamma(c + 1) - lgamma(d + 1));
            };

            double pval = get_pval(in_sum);
            double pmin = std::min(get_pval(0), get_pval(std::min(n1, m1)));

            return std::make_tuple(pval,
                                   static_cast<double>(in_sum),
                                   static_cast<double>(out_sum) / out_kmers * in_kmers,
                                   pmin);
        };
    } else if (config.test_type == "nbinom_exact") {
        common::logger->trace("Fitting per-sample negative binomial distributions");
        auto get_rp = [&](size_t j, double mu, double var, const auto &hist, double f = 1.0) {
            double var_orig = var;
            var = std::max(var, mu + 0.1);
            mu *= f;
            var *= f * f;

            double r_guess = mu * mu / (var - mu);
            double p_guess = mu / var;
            double r = r_guess;
            try {
                r = boost::math::tools::newton_raphson_iterate([&](double r) {
                    double dl = nelem * (log(r) - log(r + mu) - boost::math::digamma(r));
                    double ddl = nelem * (1.0 / r - 1.0 / (r + mu) - boost::math::trigamma(r));
                    for (const auto &[k, c] : hist) {
                        dl += boost::math::digamma(f * k + r) * c;
                        ddl += boost::math::trigamma(f * k + r) * c;
                    }
                    return std::make_pair(dl, ddl);
                }, r_guess, std::numeric_limits<double>::min(), mu * nelem, 30);
            } catch (std::exception &e) {
                common::logger->warn("Caught exception for sample {}: Falling back to initial guess", j);
                common::logger->warn("{}", e.what());
                r = r_guess;
            }

            double p = r / (r + mu);
            common::logger->trace("{}: scale: {}\tmu: {}\tvar: {}\tmoment est:\tr: {}\tp: {}\tmle: r: {}\tp: {}",
                                  j, f, mu, var_orig, r_guess, p_guess, r, p);

            return std::make_pair(r, p);
        };

        common::logger->trace("Scaling distributions");
        // geometric mean
        double target_sum = 0.0;
        for (size_t j = 0; j < groups.size(); ++j) {
            target_sum += log(static_cast<double>(sums[j]));
        }
        target_sum = exp(target_sum / groups.size());

        std::vector<std::pair<double, double>> nb_params(groups.size());
        std::vector<std::pair<double, double>> scaled_nb_params(groups.size());

        #pragma omp parallel for num_threads(num_parallel_files)
        for (size_t j = 0; j < groups.size(); ++j) {
            const auto &hist = config.test_by_unitig ? unitig_hists[j] : hists[j];
            double s = static_cast<double>(sums[j]);
            double mu = s / nelem;
            double mu2 = mu * mu;
            double var = static_cast<double>(config.test_by_unitig ? sums_of_squares_unitigs[j] : sums_of_squares[j]) / nelem - mu2;
            nb_params[j] = get_rp(j, mu, var, hist);
            scaled_nb_params[j] = get_rp(j, mu, var, hist, target_sum / s);
        }

        common::logger->trace("Finding common p");
        double p_guess = 0.0;
        for (size_t j = 0; j < groups.size(); ++j) {
            p_guess += log(scaled_nb_params[j].second);
        }
        p_guess = exp(p_guess / groups.size());

        double target_p = p_guess;
        double r_map = target_p / nelem / (1 - target_p) * target_sum;

        std::vector<VectorMap<uint64_t, std::pair<size_t, uint64_t>>> count_maps(groups.size());
        double r_in = 0;
        double r_out = 0;
        std::mutex r_mu;

        common::logger->trace("Computing quantile maps");
        #pragma omp parallel for num_threads(num_parallel_files)
        for (size_t j = 0; j < groups.size(); ++j) {
            const auto &[r, p] = nb_params[j];

            common::logger->trace("{}\tApproximating NB({}, {}) with NB({}, {})",
                                  j, r, p, r_map, target_p);

            boost::math::negative_binomial nb(r, p);
            boost::math::negative_binomial nb_out(r_map, target_p);
            double scale = exp(log(p) - log(target_p));

            // ensure that 0 maps to 0
            count_maps[j][0] = std::make_pair(0, 0);

            for (const auto &[k, c] : hists[j]) {
                if (!count_maps[j].count(k)) {
                    double cdf = boost::math::cdf(nb, k);

                    uint64_t new_k = 0;
                    if (cdf < 1.0) {
                        new_k = boost::math::quantile(nb_out, cdf);
                    } else {
                        double ccdf = boost::math::cdf(boost::math::complement(nb, k));
                        if (ccdf > 0.0) {
                            new_k = boost::math::quantile(boost::math::complement(nb_out, ccdf));
                        } else {
                            new_k = ceil(scale * k);
                        }
                    }
                    count_maps[j][k] = std::make_pair(std::max(new_k, uint64_t(1)), c);
                } else {
                    count_maps[j][k].second += c;
                }
            }

            std::lock_guard<std::mutex> lock(r_mu);
            if (groups[j]) {
                r_out += r_map;
            } else {
                r_in += r_map;
            }
        }

        common::logger->trace("Fits: in: {}\tout: {}\tp: {}", r_in, r_out, target_p);
        common::logger->trace("Unscaled Totals: in: {}\tout: {}", in_kmers, out_kmers);

        double in_kmers_adj = 0;
        double out_kmers_adj = 0;
        for (size_t j = 0; j < groups.size(); ++j) {
            for (const auto &[k, m] : count_maps[j]) {
                const auto &[v, c] = m;
                if (groups[j]) {
                    out_kmers_adj += v * c;
                } else {
                    in_kmers_adj += v * c;
                }
            }
        }

        common::logger->trace("Scaled Totals: in: {}\tout: {}", in_kmers_adj, out_kmers_adj);

        compute_pval = [&,r_in,r_out,count_maps,in_kmers_adj,out_kmers_adj](const auto &row) {
            if (row.empty())
                return std::make_tuple(1.1, 0.0, 0.0, 1.1);

            int64_t in_sum = 0;
            int64_t out_sum = 0;
            for (const auto &[j, c] : row) {
                if (groups[j]) {
                    out_sum += count_maps[j].find(c)->second.first;
                } else {
                    in_sum += count_maps[j].find(c)->second.first;
                }
            }

            int64_t n = in_sum + out_sum;
            double p_min = 1.0;
            double pval = 1.0;

            double lscaling = lgamma(n + 1) + lgamma(r_in + r_out) - lgamma(r_in) - lgamma(r_out);
            auto get_pmf = [&](int64_t s) {
                int64_t t = n - s;
                double rs = r_in + s;
                double rt = r_out + t;
                return exp(lscaling - lgamma(s + 1) - lgamma(t + 1) + lgamma(rs) + lgamma(rt) - lgamma(rs + rt));
            };

            if (num_labels_in == num_labels_out) {
                p_min = get_pmf(0) + get_pmf(n);
                if (in_sum != out_sum) {
                    lscaling = lgamma(n + 1) + lgamma(r_in + r_out) - lgamma(r_in) - lgamma(r_out);
                    int64_t d = std::min(in_sum, out_sum);
                    pval = p_min;
                    for (int64_t s = 1; s <= d; ++s) {
                        pval += get_pmf(s) + get_pmf(n - s);
                    }
                }
            } else {
                p_min = std::min(get_pmf(0), get_pmf(n));
                double d = abs(static_cast<double>(in_sum) / num_labels_in - static_cast<double>(out_sum) / num_labels_out);
                if (d > 0) {
                    // we want s s.t. | s/num_labels_in - (n-s)/num_labels_out | >= d
                    // s/num_labels_in - (n-s)/num_labels_out >= d || (n-s)/num_labels_out - s/num_labels_in >= d
                    // s/num_labels_in + s/num_labels_out >= d + n/num_labels_out || n/num_labels_out - d >= s/num_labels_in + s/num_labels_out
                    double factor = double(1.0) / num_labels_in + double(1.0) / num_labels_out;

                    double lower_d = (static_cast<double>(n) / num_labels_out - d) / factor;
                    double upper_d = (d + static_cast<double>(n) / num_labels_out) / factor;

                    // fix floating point errors
                    if (abs(lower_d - round(lower_d)) < 1e-10)
                        lower_d = round(lower_d);

                    if (abs(upper_d - round(upper_d)) < 1e-10)
                        upper_d = round(upper_d);

                    int64_t lower = floor(lower_d);
                    int64_t upper = ceil(upper_d);

                    if (lower < n && upper > 0 && lower != upper) {
                        pval = 0;
                        for (int64_t s = 0; s <= lower; ++s) {
                            pval += get_pmf(s);
                        }

                        for (int64_t s = upper; s <= n; ++s) {
                            pval += get_pmf(s);
                        }
                    }
                }
            }

            if (pval >= 1.1)
                throw std::runtime_error("pval fail");

            return std::make_tuple(std::min(pval, 1.0),
                                   static_cast<double>(in_sum),
                                   static_cast<double>(out_sum) / out_kmers_adj * in_kmers_adj,
                                   std::min(p_min, 1.0));
        };

        compute_pval_unitig = [&,r_in,r_out,count_maps,in_kmers_adj,out_kmers_adj](const auto &rows) {
            if (std::all_of(rows.begin(), rows.end(), [](const auto &row) { return row.empty(); }))
                return std::make_tuple(1.1, 0.0, 0.0, 1.1);

            int64_t in_sum = 0;
            int64_t out_sum = 0;
            for (const auto &row : rows) {
                for (const auto &[j, c] : row) {
                    if (groups[j]) {
                        out_sum += count_maps[j].find(c)->second.first;
                    } else {
                        in_sum += count_maps[j].find(c)->second.first;
                    }
                }
            }

            int64_t n = in_sum + out_sum;
            double p_min = 1.0;
            double pval = 1.0;

            double lscaling = lgamma(n + 1) + lgamma(r_in + r_out) - lgamma(r_in) - lgamma(r_out);
            auto get_pmf = [&](int64_t s) {
                int64_t t = n - s;
                double rs = r_in + s;
                double rt = r_out + t;
                return exp(lscaling - lgamma(s + 1) - lgamma(t + 1) + lgamma(rs) + lgamma(rt) - lgamma(rs + rt));
            };

            if (num_labels_in == num_labels_out) {
                p_min = get_pmf(0) + get_pmf(n);
                if (in_sum != out_sum) {
                    lscaling = lgamma(n + 1) + lgamma(r_in + r_out) - lgamma(r_in) - lgamma(r_out);
                    int64_t d = std::min(in_sum, out_sum);
                    pval = p_min;
                    for (int64_t s = 1; s <= d; ++s) {
                        pval += get_pmf(s) + get_pmf(n - s);
                    }
                }
            } else {
                p_min = std::min(get_pmf(0), get_pmf(n));
                double d = abs(static_cast<double>(in_sum) / num_labels_in - static_cast<double>(out_sum) / num_labels_out);
                if (d > 0) {
                    // we want s s.t. | s/num_labels_in - (n-s)/num_labels_out | >= d
                    // s/num_labels_in - (n-s)/num_labels_out >= d || (n-s)/num_labels_out - s/num_labels_in >= d
                    // s/num_labels_in + s/num_labels_out >= d + n/num_labels_out || n/num_labels_out - d >= s/num_labels_in + s/num_labels_out
                    double factor = double(1.0) / num_labels_in + double(1.0) / num_labels_out;

                    double lower_d = (static_cast<double>(n) / num_labels_out - d) / factor;
                    double upper_d = (d + static_cast<double>(n) / num_labels_out) / factor;

                    // fix floating point errors
                    if (abs(lower_d - round(lower_d)) < 1e-10)
                        lower_d = round(lower_d);

                    if (abs(upper_d - round(upper_d)) < 1e-10)
                        upper_d = round(upper_d);

                    int64_t lower = floor(lower_d);
                    int64_t upper = ceil(upper_d);

                    if (lower < n && upper > 0 && lower != upper) {
                        pval = 0;
                        for (int64_t s = 0; s <= lower; ++s) {
                            pval += get_pmf(s);
                        }

                        for (int64_t s = upper; s <= n; ++s) {
                            pval += get_pmf(s);
                        }
                    }
                }
            }

            if (pval >= 1.1)
                throw std::runtime_error("pval fail");

            return std::make_tuple(std::min(pval, 1.0),
                                   static_cast<double>(in_sum),
                                   static_cast<double>(out_sum) / out_kmers_adj * in_kmers_adj,
                                   std::min(p_min, 1.0));
        };
    } else if (config.test_type == "notest") {
        compute_pval = [&](const auto &row) {
            if (row.empty())
                return std::make_tuple(1.1, 0.0, 0.0, 1.1);

            int64_t in_sum = 0;
            int64_t out_sum = 0;

            for (const auto &[j, c] : row) {
                if (groups[j]) {
                    out_sum += c;
                } else {
                    in_sum += c;
                }
            }

            return std::make_tuple(0.0,
                                   static_cast<double>(in_sum),
                                   static_cast<double>(out_sum) / out_kmers * in_kmers,
                                   0.0);
        };
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

            std::vector<PairContainer> rows;
            rows.reserve(path.size());
            for (node_index node : path) {
                auto &row = rows.emplace_back();
                uint64_t row_i = AnnotatedDBG::graph_to_anno_index(node);
                for (size_t j = 0; j < groups.size(); ++j) {
                    const auto &col = *columns_all[j];
                    const auto &col_vals = *column_values_all[j];
                    if (uint64_t r = col.conditional_rank1(row_i))
                        row.emplace_back(j, col_vals[r - 1]);
                }
            }

            double pval;
            double in_stat;
            double out_stat;
            double pval_min;
            try {
                std::tie(pval, in_stat, out_stat, pval_min) = compute_pval_unitig(rows);
            } catch (...) {
                std::lock_guard<std::mutex> lock(pval_mu);
                ex = std::current_exception();
                return;
            }

            if (pval >= 1.1)
                return;

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

            PairContainer row;
            if (kept[row_i]) {
                size_t count_in = 0;
                size_t count_out = 0;
                for (const auto &[j, raw_c] : raw_row) {
                    if (uint64_t c = raw_c) {
                        row.emplace_back(j, c);
                        if (groups[j]) {
                            ++count_out;
                        } else {
                            ++count_in;
                        }
                    }
                }

                if (count_in + count_out < config.min_recurrence
                        || count_in < config.min_in_recurrence
                        || count_out < config.min_out_recurrence
                        || count_in > config.max_in_recurrence
                        || count_out > config.max_out_recurrence)
                    row.clear();
            }

            node_index node = AnnotatedDBG::anno_to_graph_index(row_i);

            double pval;
            double in_stat;
            double out_stat;
            double pval_min;
            try {
                std::tie(pval, in_stat, out_stat, pval_min) = compute_pval(row);
            } catch (...) {
                std::lock_guard<std::mutex> lock(pval_mu);
                ex = std::current_exception();
                return;
            }

            if (pval >= 1.1)
                return;

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
        size_t total_tests = acc;
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
                // TODO: fix for unitig
                common::logger->trace("Falling back to Benjamini-Yekutieli");
                VectorMap<double, double> pval_map;
                for (uint64_t p_u : pvals) {
                    auto p = bit_cast<double>(p_u);
                    if (p <= 1.0)
                        ++pval_map[p];
                }
                auto p_data = const_cast<std::vector<std::pair<double, double>>&&>(pval_map.values_container());
                std::sort(p_data.begin(), p_data.end(), utils::GreaterFirst());

                double cm = boost::math::digamma(total_tests) + 0.5772156649015329; // std::numbers::e_gamma_v<double>
                pval_map = decltype(pval_map)();
                size_t cur_rank = 0;
                for (const auto &[p, c] : p_data) {
                    cur_rank += c;
                    pval_map[p] = cm * cur_rank;
                }
                for (size_t i = 0; i < pvals.size(); ++i) {
                    double p = bit_cast<double, uint64_t>(pvals[i]);
                    if (p < config.family_wise_error_rate && p * pval_map[p] >= config.family_wise_error_rate) {
                        indicator_in[i] = false;
                        indicator_out[i] = false;
                    }
                }
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


