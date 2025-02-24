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
#include <boost/math/distributions/cauchy.hpp>

#include "common/logger.hpp"
#include "common/vectors/bitmap.hpp"
#include "common/vector_map.hpp"
#include "common/vectors/transpose.hpp"
#include "common/hashers/hash.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/graph_cleaning.hpp"
#include "graph/representation/canonical_dbg.hpp"
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

enum Group {
    IN,
    OUT,
    OTHER,
    BOTH
};


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
double get(typename PValStorage::reference&& ref) {
    return bit_cast<double>(static_cast<typename PValStorage::value_type>(ref));
}

template <class PValStorage>
void set(typename PValStorage::reference&& ref, double val) {
    ref = bit_cast<typename PValStorage::value_type>(val);
}

template <class PValStorage>
void push_back(PValStorage &v, double val) {
    v.push_back(bit_cast<typename PValStorage::value_type>(val));
}

template <typename value_type, class PValStorage, typename HistGetter, typename Generator>
std::tuple<std::shared_ptr<DeBruijnGraph>, std::shared_ptr<DeBruijnGraph>, PValStorage, std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                         const HistGetter &get_hist_map,
                         const Generator &generate_rows,
                         const std::vector<Group> &groups,
                         const DifferentialAssemblyConfig &config,
                         size_t num_threads = 1,
                         std::filesystem::path tmp_dir = "",
                         size_t num_parallel_files = std::numeric_limits<size_t>::max(),
                         const std::function<void()> &deallocate = []() {},
                         uint8_t max_width = 64) {
    if (auto canonical = std::dynamic_pointer_cast<const CanonicalDBG>(graph_ptr))
        graph_ptr = canonical->get_graph_ptr();

    num_parallel_files = get_num_threads();
    num_threads = get_num_threads();
    bool is_primary = graph_ptr->get_mode() == DeBruijnGraph::PRIMARY;

    size_t num_labels_both = std::count(groups.begin(), groups.end(), Group::BOTH);
    size_t num_labels_out = std::count(groups.begin(), groups.end(), Group::OUT);
    size_t num_labels_in = std::count(groups.begin(), groups.end(), Group::IN);
    num_labels_in += num_labels_both;
    num_labels_out += num_labels_both;

    common::logger->trace("Graph mode: {}", is_primary ? "PRIMARY" : "other");

    common::logger->trace("Computing histogram from uncleaned counts");
    auto hists_map = get_hist_map(std::vector<size_t>(groups.size(), 1), nullptr);

    sdsl::bit_vector kept_bv(AnnotatedDBG::graph_to_anno_index(graph_ptr->max_index() + 1), true);

    std::vector<size_t> min_counts(groups.size(), config.min_count);
    std::vector<uint64_t> check_cutoff(groups.size(), std::numeric_limits<uint64_t>::max());

    // std::mutex agg_mu;
    // using PairContainer = std::vector<std::pair<uint64_t, value_type>>;

    if (config.clean) {
        common::logger->trace("Cleaning count columns");

        #pragma omp parallel for num_threads(num_parallel_files)
        for (size_t j = 0; j < groups.size(); ++j) {
            if (groups[j] == Group::OTHER) {
                min_counts[j] = std::numeric_limits<size_t>::max();
                continue;
            }

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

    auto correct_pvals = [&config](const std::vector<std::pair<double, size_t>> &m) {
        size_t total = std::accumulate(m.begin(), m.end(), size_t(0),
                                       [](size_t sum, const auto &a) { return sum + a.second; });
        size_t acc = 0;

        for (const auto &[pval_min, s] : m) {
            if (s == 0)
                continue;

            size_t c = floor(config.family_wise_error_rate / pval_min + 1);
            size_t nsig = total - acc;
            // there are nsig k-mers s.t. pval_min < alpha / c
            if (nsig > 0 && c >= nsig)
                return std::make_pair(c, nsig);

            acc += s;
        }

        return std::make_pair<size_t, size_t>(std::numeric_limits<size_t>::max(),
                                              std::numeric_limits<size_t>::max());
    };

    common::logger->trace("Updating histogram and marking discarded k-mers");

    if (groups.size())
        hists_map = get_hist_map(min_counts, &kept_bv);

    bit_vector_stat kept(std::move(kept_bv));
    size_t nelem = kept.num_set_bits();

    std::unique_ptr<MaskedDeBruijnGraph> clean_masked_graph;
    if (config.test_type != "notest" && config.test_by_unitig) {
        clean_masked_graph = std::make_unique<MaskedDeBruijnGraph>(
            graph_ptr,
            [&](node_index node) {
                return node != DeBruijnGraph::npos
                        && kept[AnnotatedDBG::graph_to_anno_index(node)];
            },
            true,
            is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
        );
    }

    std::vector<std::vector<std::pair<uint64_t, size_t>>> hists(groups.size());
    for (size_t j = 0; j < hists.size(); ++j) {
        hists[j] = const_cast<std::vector<std::pair<uint64_t, size_t>>&&>(hists_map[j].values_container());
        std::sort(hists[j].begin(), hists[j].end(), utils::LessFirst());
        hists_map[j].clear();
    }
    hists_map.clear();

    for (size_t j = 0; j < groups.size(); ++j) {
        size_t total_c = 0;
        for (const auto &[k, c] : hists[j]) {
            total_c += c;
        }

        assert(total_c == nelem);

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

    int64_t total_kmers = 0;
    for (size_t j = 0; j < groups.size(); ++j) {
        for (const auto &[k, c] : hists[j]) {
            sums[j] += k * c;
            if (k != 0)
                n_kmers[j] += c;

            max_obs_vals[j] = std::max(max_obs_vals[j], k);
        }

        if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
            out_kmers += sums[j];

        if (groups[j] == Group::IN || groups[j] == Group::BOTH)
            in_kmers += sums[j];

        total_kmers += sums[j];

        common::logger->trace("{}: n_unique: {}\tsum: {}\tmax_obs: {}\tmin_cutoff: {}\tmax_cutof: {}",
                              j, n_kmers[j], sums[j], max_obs_vals[j], min_counts[j], check_cutoff[j]);
    }

    common::logger->trace("Number of kept unique k-mers: {}\tNumber of kept k-mers: {}",
                          nelem, total_kmers);

    // common::logger->trace("Allocating p-value storage");
    // auto nullpval = bit_cast<uint64_t>(double(1.1));

    // std::vector<PValStorage> mean_buckets;
    // std::vector<std::unique_ptr<utils::TempFile>> tmp_buckets;
    // if constexpr(std::is_same_v<PValStorage, std::vector<uint64_t>>) {
    //     auto &means = mean_buckets.emplace_back();
    //     means.resize(graph_ptr->max_index() + 1, 0);
    // }

    // if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
    //     for (size_t i = 0; i < get_num_threads() + 1; ++i) {
    //         auto &tmp_file = tmp_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
    //         mean_buckets.emplace_back(tmp_file->name(), std::ios::out);
    //     }
    //     mean_buckets[0].push_back(0);
    // }

    std::vector<PValStorage> pvals_buckets;
    std::vector<std::unique_ptr<utils::TempFile>> tmp_buckets;
    std::vector<PValStorage> pvals_min_buckets;
    std::vector<std::unique_ptr<utils::TempFile>> tmp_min_buckets;
    std::vector<PValStorage> eff_size_buckets;
    std::vector<std::unique_ptr<utils::TempFile>> tmp_eff_buckets;
    // constexpr bool preallocated = std::is_same_v<PValStorage, std::vector<uint64_t>>;

    if (config.test_type != "notest" && config.test_by_unitig) {
        if constexpr(std::is_same_v<PValStorage, std::vector<uint64_t>>) {
            auto &pvals = pvals_buckets.emplace_back();
            pvals.reserve(graph_ptr->max_index());

            auto &pvals_min = pvals_min_buckets.emplace_back();
            pvals_min.reserve(graph_ptr->max_index());

            auto &eff_size = eff_size_buckets.emplace_back();
            eff_size.reserve(graph_ptr->max_index());
        }

        if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
            for (size_t i = 0; i < get_num_threads() + 1; ++i) {
                auto &tmp_file = tmp_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
                pvals_buckets.emplace_back(tmp_file->name(), std::ios::out);
            }

            for (size_t i = 0; i < get_num_threads() + 1; ++i) {
                auto &tmp_file = tmp_min_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
                pvals_min_buckets.emplace_back(tmp_file->name(), std::ios::out);
            }

            for (size_t i = 0; i < get_num_threads() + 1; ++i) {
                auto &tmp_file_eff = tmp_eff_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
                eff_size_buckets.emplace_back(tmp_file_eff->name(), std::ios::out);
            }
        }
    }

    // if ((config.test_type == "gnb_exact" || config.test_type == "lnb_exact") && config.test_by_unitig) {

    // }


    // std::vector<PValStorage> sum_buckets;
    // std::vector<std::unique_ptr<utils::TempFile>> tmp_sum_buckets;
    // if (config.test_by_unitig) {
    //     common::logger->trace("Allocating effect size and row sum storage");

    // }

    sdsl::bit_vector indicator_in;
    sdsl::bit_vector indicator_out;

    // if (!config.test_by_unitig) {
        common::logger->trace("Allocating k-mer bitmasks");
        indicator_in = sdsl::bit_vector(graph_ptr->max_index() + 1, false);
        indicator_out = sdsl::bit_vector(graph_ptr->max_index() + 1, false);
    // }

    // std::function<double(int64_t, int64_t, const PairContainer&)> compute_pval;
    // std::function<double(int64_t, const PairContainer&)> compute_min_pval;

    // std::vector<VectorMap<uint64_t, std::pair<size_t, uint64_t>>> count_maps;

    common::logger->trace("Test: {}\tby unitig: {}", config.test_type, config.test_by_unitig);

    // prefilter rows
    auto generate_clean_rows = [&](const auto &callback) {
        generate_rows([&](uint64_t row_i, const auto &row, size_t bucket_idx) {
            if (!kept[row_i])
                return;

            size_t count_in = 0;
            size_t count_out = 0;
            size_t total_count = 0;
            for (const auto &[j, c] : row) {
                if (min_counts.size() && c < min_counts[j])
                    continue;

                if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                    ++count_out;

                if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                    ++count_in;

                total_count += groups[j] != Group::OTHER;
            }

            if (total_count < config.min_recurrence)
                return;

            bool in_kmer = count_in >= config.min_in_recurrence && count_in <= config.max_in_recurrence;
            bool out_kmer = count_out >= config.min_out_recurrence && count_out <= config.max_out_recurrence;

            if (in_kmer || out_kmer)
                callback(row_i, row, bucket_idx);
        });
    };

    // precompute negative binomial fits for poisson_bayes and nbinom_exact tests
    std::vector<std::pair<double, double>> nb_params(groups.size());
    std::pair<double, double> nb_params_a;
    std::pair<double, double> nb_params_b;
    std::pair<double, double> nb_params_null;
    double nb_base = 0.0;
    if (config.test_type == "poisson_bayes" || config.test_type == "nbinom_exact") {
        common::logger->trace("Fitting per-sample negative binomial distributions");
        auto get_rp = [&](const auto &generate) {
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
                return dl;
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
                get_dl, r_min, r_max, boost::math::tools::eps_tolerance<double>(5)
            );

            double r = r_max;
            double p = r / (r + mu);
            return std::make_pair(r, p);
        };

        #pragma omp parallel for num_threads(num_parallel_files)
        for (size_t j = 0; j < groups.size(); ++j) {
            const auto &hist = hists[j];
            if (hist.size()) {
                nb_params[j] = get_rp([&](const auto &callback) {
                    for (const auto &[k, c] : hist) {
                        callback(k, c);
                    }
                });
                const auto &[r, p] = nb_params[j];
                common::logger->trace("{}: size: {}\tmax_val: {}\tmu: {}\tvar: {}\tmle: r: {}\tp: {}",
                                    j, sums[j], (hist.end() - 1)->first, r * (1-p)/p, r*(1-p)/p/p, r, p);
            }
        }

        if (config.test_type == "nbinom_exact") {
            double mu_a = 0.0;
            double mu_b = 0.0;
            double var_a = 0.0;
            double var_b = 0.0;
            for (size_t j = 0; j < groups.size(); ++j) {
                const auto &[r, p] = nb_params[j];
                double mu = r * (1.0 - p) / p;
                double var = mu / p;
                if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
                    mu_a += mu;
                    var_a += var;
                }

                if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
                    mu_b += mu;
                    var_b += var;
                }
            }

            nb_params_a.first = mu_a * mu_a / (var_a - mu_a);
            nb_params_a.second = mu_a / var_a;
            common::logger->trace("In: r: {}\tp: {}", nb_params_a.first, nb_params_a.second);

            nb_params_b.first = mu_b * mu_b / (var_b - mu_b);
            nb_params_b.second = mu_b / var_b;
            common::logger->trace("Out: r: {}\tp: {}", nb_params_b.first, nb_params_b.second);

            double mu_null = mu_a + mu_b;
            double var_null = var_a + var_b; // the two are independent, so their covariance is 0
            nb_params_null.first = mu_null * mu_null / (var_null - mu_null);
            nb_params_null.second = mu_null / var_null;
            common::logger->trace("Null: r: {}\tp: {}", nb_params_null.first, nb_params_null.second);

            nb_base = (lgamma(nb_params_null.first) - lgamma(nb_params_a.first) - lgamma(nb_params_b.first)) / log(2.0);
            nb_base += nb_params_a.first*log2(nb_params_a.second) + nb_params_b.first*log2(nb_params_b.second)-nb_params_null.first*log2(nb_params_null.second);
        }
    }

    // precompute p-values for poisson_binom test
    std::vector<std::vector<double>> pb_pvals(num_labels_in + num_labels_out + 1);
    std::vector<double> mid_points(num_labels_in + num_labels_out + 1);
    if (config.test_type == "poisson_binom") {
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
            if (groups[i - 1] == Group::OUT || groups[i - 1] == Group::BOTH) {
                std::vector<double> pmf_out_cur(pmf_out.size() + 1);
                pmf_out_cur[0] = (1.0 - p[i - 1]) * pmf_out[0];
                pmf_out_cur[pmf_out.size()] = p[i - 1] * pmf_out.back();
                for (size_t k = 1; k < pmf_out.size(); ++k) {
                    pmf_out_cur[k] = p[i - 1] * pmf_out[k - 1] + (1.0 - p[i - 1]) * pmf_out[k];
                }
                std::swap(pmf_out_cur, pmf_out);
            }

            if (groups[i - 1] == Group::IN || groups[i - 1] == Group::BOTH) {
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

        if (pmf_in.size() != num_labels_in + 1) {
            common::logger->error("PMF in wrong: {} != {}", pmf_in.size(), num_labels_in + 1);
            throw std::domain_error("");
        }

        if (pmf_out.size() != num_labels_out + 1) {
            common::logger->error("PMF out wrong: {} != {}", pmf_out.size(), num_labels_out + 1);
            throw std::domain_error("");
        }

        if (pmf_null.size() != num_labels_in + num_labels_out + 1) {
            common::logger->error("PMF null wrong: {} != {}", pmf_null.size(), num_labels_in + num_labels_out + 1);
            throw std::domain_error("");
        }

        common::logger->trace("Precomputing p-values");
        pb_pvals[0].emplace_back(1.0);

        double min_pval = 1.0;

        for (size_t n = 1; n < pb_pvals.size(); ++n) {
            std::vector<double> probs;
            size_t front = n - std::min(n, num_labels_out);
            for (uint64_t s = 0; s < pmf_in.size(); ++s) {
                if (s > n)
                    break;
                uint64_t t = n - s;
                if (s < pmf_in.size() && t < pmf_out.size()) {
                    if (s < front) {
                        common::logger->error("Attempting non-zero p-value in impossible configuration: {},{}", n, s);
                        throw std::domain_error("");
                    }
                    probs.emplace_back(exp2(log2(pmf_in[s]) + log2(pmf_out[t]) - log2(pmf_null[n])));
                } else {
                    probs.emplace_back(0.0);
                }
            }

            for (uint64_t s = 0; s < probs.size(); ++s) {
                pb_pvals[n].emplace_back(0.0);
                if (s > n)
                    break;
                uint64_t t = n - s;
                if (s >= pmf_in.size() || t >= pmf_out.size())
                    continue;
                for (uint64_t sp = 0; sp < probs.size(); ++sp) {
                    if (sp > n)
                        break;
                    uint64_t tp = n - sp;
                    if (sp >= pmf_in.size() || tp >= pmf_out.size())
                        continue;
                    if (probs[sp] <= probs[s])
                        pb_pvals[n][s] += probs[sp];
                }
            }

            for (uint64_t s = 0; s < front; ++s) {
                if (pb_pvals[n][s] > 0) {
                    common::logger->error("Non-zero p-value in impossible configuration: {},{}", n, s);
                    throw std::domain_error("");
                }
            }
            double local_min_pval = std::min(pb_pvals[n][front],
                                            pb_pvals[n].back());
            min_pval = std::min(min_pval, local_min_pval);

            size_t max_prob_counts = 0;
            double max_prob_pos = 0;
            // double max_prob = -1;

            for (uint64_t s = front; s < pb_pvals[n].size(); ++s) {
                if (pb_pvals[n][s] < local_min_pval) {
                    common::logger->error("Min p-value not at boundary: {},{}: {} < {}", n, s, pb_pvals[n][s], local_min_pval);
                    throw std::domain_error("");
                }
            }
            if (front + 1 == pb_pvals[n].size()) {
                max_prob_counts = 1;
                max_prob_pos = front;
            } else if (front + 2 == pb_pvals[n].size()) {
                max_prob_counts = 2;
                max_prob_pos = front + front + 1;
            } else {
                for (uint64_t s = front + 1; s < pb_pvals[n].size(); ++s) {
                    if (probs[s] == probs[s-1]) {
                        max_prob_counts = 2;
                        max_prob_pos = s + s;
                        break;
                    }

                    if (s + 1 < pb_pvals[n].size()) {
                        if (probs[s-1] < probs[s] && probs[s] < probs[s+1])
                            continue;

                        if (probs[s-1] > probs[s] && probs[s] > probs[s+1])
                            continue;

                        if (probs[s] > probs[s-1] && probs[s] > probs[s+1]) {
                            if (probs[s] - probs[s-1] == probs[s]-probs[s+1]) {
                                max_prob_counts = 1;
                                max_prob_pos = s;
                            } else if (probs[s] - probs[s-1] > probs[s]-probs[s+1]) {
                                max_prob_counts = 2;
                                max_prob_pos = s + s + 1;
                            } else {
                                max_prob_counts = 2;
                                max_prob_pos = s + s - 1;
                            }
                        }
                    }
                }
            }

            if (max_prob_counts == 0) {
                max_prob_counts = 2;
                max_prob_pos = front - 1 + front;
            }

            mid_points[n] = static_cast<double>(max_prob_pos) / max_prob_counts;
            common::logger->trace("Midpoint: n: {}\tmp: {}\t{}\t{}",
                                    n, mid_points[n],
                                    fmt::join(probs,","),
                                    fmt::join(pb_pvals[n],","));
        }

        common::logger->trace("Min. p-value: {}", min_pval);
        if (min_pval >= config.family_wise_error_rate) {
            common::logger->warn("No significant p-values achievable");
            auto masked_graph_in = std::make_shared<MaskedDeBruijnGraph>(
                graph_ptr, [](node_index) { return false; }, true,
                is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
            );

            auto masked_graph_out = std::make_shared<MaskedDeBruijnGraph>(
                graph_ptr, [](node_index) { return false; }, true,
                is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
            );

            return std::make_tuple(masked_graph_in, masked_graph_out, PValStorage{}, nullptr);
        }
    }

    auto combine_pvals = [dist=boost::math::cauchy()](const std::vector<double> &pvals) {
        double stat = 0.0;
        for (double pval : pvals) {
            stat += tan((0.5 - pval) * M_PI);
        }
        stat /= pvals.size();
        return boost::math::cdf(boost::math::complement(dist, stat));
    };

    // precompute minimal attainable p-values
    std::vector<std::pair<double, size_t>> m;
    std::vector<size_t> offsets;
    if (config.test_type != "notest") {
        common::logger->trace("Computing minimum p-values");
        VectorMap<double, size_t> m_vec;
        std::vector<VectorMap<double, size_t>> ms(num_threads + 1);
        generate_clean_rows([&](uint64_t row_i, const auto &row, uint64_t bucket_idx) {
            double pval = 1.0;
            if (config.test_type == "poisson_binom") {
                size_t n = 0;
                for (const auto &[j, c] : row) {
                    n += static_cast<int64_t>(groups[j] != Group::OTHER) + (groups[j] == Group::BOTH);
                }
                size_t front = n - std::min(n, num_labels_out);
                pval = std::min(pb_pvals[n][front], pb_pvals[n].back());
            } else if (config.test_type == "poisson_exact" || config.test_type == "poisson_bayes") {
                double mu1 = 0.0;
                double mu2 = 0.0;
                if (config.test_type == "poisson_exact") {
                    mu1 = static_cast<double>(in_kmers) / nelem;
                    mu2 = static_cast<double>(out_kmers) / nelem;
                }
                int64_t n = 0;
                for (const auto &[j, c] : row) {
                    n += c;
                    if (config.test_type == "poisson_bayes") {
                        const auto &[r, p] = nb_params[j];
                        double lambda_j = sqrt((1.0 - p) * (c + r - 1));
                        if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                            mu2 += lambda_j;

                        if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                            mu1 += lambda_j;
                    }
                }

                if (n > 0) {
                    double p = mu1 / (mu1 + mu2);
                    boost::math::binomial bdist(n, p);
                    double pval0 = boost::math::pdf(bdist, 0);
                    double pvaln = boost::math::pdf(bdist, n);
                    pval = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                }
            } else if (config.test_type == "nbinom_exact") {
                int64_t n = 0;
                for (const auto &[j, c] : row) {
                    n += c;
                }

                if (n > 0) {
                    auto [r_n, p_n] = nb_params_null;
                    auto [r_a, p_a] = nb_params_a;
                    auto [r_b, p_b] = nb_params_b;

                    double base = nb_base + (lgamma(n + 1) - lgamma(r_n + n)) / log(2.0) - n*log2(1.0-p_n);
                    auto get_pval = [&](int64_t s) {
                        int64_t t = n - s;
                        double sbase = (lgamma(r_a+s)+lgamma(r_b+t)-lgamma(s+1)-lgamma(t+1)) / log(2.0);
                        return exp2(base + sbase +s*log2(1.0-p_a) + t*log2(1.0-p_b));
                    };

                    double pval0 = get_pval(0);
                    double pvaln = get_pval(n);
                    pval = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                }
            }

            if (!config.test_by_unitig) {
                ++ms[bucket_idx][pval];
            } else {
                push_back(pvals_min_buckets[bucket_idx], pval);
            }
        });

        if (config.test_by_unitig) {
            common::logger->trace("Combining minimal attainable p-values for unitigs");
            offsets.reserve(num_threads + 1);
            offsets.emplace_back(0);
            for (const auto &bucket : pvals_min_buckets) {
                offsets.emplace_back(offsets.back() + bucket.size());
            }
            std::mutex m_vec_mu;
            clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
                std::vector<double> pvals;
                pvals.reserve(path.size());
                {
                    std::lock_guard<std::mutex> lock(m_vec_mu);
                    for (node_index node : path) {
                        size_t row = AnnotatedDBG::graph_to_anno_index(node);
                        size_t row_rank = kept.rank1(row);
                        size_t bucket_idx = std::lower_bound(offsets.begin(), offsets.end(), row_rank) - offsets.begin() - 1;
                        size_t offset = row_rank - offsets[bucket_idx];
                        pvals.emplace_back(get<PValStorage>(pvals_min_buckets[bucket_idx][offset]));
                    }
                }
                double comp_min_pval = combine_pvals(pvals);
                size_t row = AnnotatedDBG::graph_to_anno_index(path[0]);
                size_t row_rank = kept.rank1(row);
                size_t bucket_idx = std::lower_bound(offsets.begin(), offsets.end(), row_rank) - offsets.begin() - 1;
                size_t offset = row_rank - offsets[bucket_idx];
                ++ms[bucket_idx][comp_min_pval];
                std::lock_guard<std::mutex> lock(m_vec_mu);
                set<PValStorage>(pvals_min_buckets[bucket_idx][offset], comp_min_pval);
            });
        }

        std::swap(m_vec, ms[0]);
        for (size_t i = 1; i < ms.size(); ++i) {
            for (const auto &[pval, c] : ms[i]) {
                m_vec[pval] += c;
            }
        }

        m = const_cast<std::vector<std::pair<double, size_t>>&&>(m_vec.values_container());
        std::sort(m.rbegin(), m.rend());
    }

    // determine cutoffs for multiple testing correction
    auto [k_min, k] = correct_pvals(m);

    common::logger->trace("Picked: k: {}\tn: {} / {}\tpval_min: {}", k_min, k, nelem, config.family_wise_error_rate / k);

    common::logger->trace("Running differential tests");
    std::atomic_thread_fence(std::memory_order_release);
    generate_clean_rows([&](uint64_t row_i, const auto &row, uint64_t bucket_idx) {
        if (config.test_type == "notest") {
            size_t count_in = 0;
            size_t count_out = 0;
            for (const auto &[j, c] : row) {
                if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                    ++count_out;

                if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                    ++count_in;
            }

            bool in_kmer = count_in >= config.min_in_recurrence && count_in <= config.max_in_recurrence;
            bool out_kmer = count_out >= config.min_out_recurrence && count_out <= config.max_out_recurrence;

            node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
            if (in_kmer)
                set_bit(indicator_in.data(), node, true, std::memory_order_relaxed);

            if (out_kmer)
                set_bit(indicator_out.data(), node, true, std::memory_order_relaxed);

            return;
        }

        auto set_pval = [&](double pval, double eff_size) {
            if (!config.test_by_unitig) {
                if (pval * k < config.family_wise_error_rate) {
                    node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
                    if (eff_size > 0) {
                        set_bit(indicator_in.data(), node, true, std::memory_order_relaxed);
                    } else if (eff_size < 0) {
                        set_bit(indicator_out.data(), node, true, std::memory_order_relaxed);
                    }
                }
            } else {
                push_back(pvals_buckets[bucket_idx], pval);
                push_back(eff_size_buckets[bucket_idx], eff_size);
            }
        };

        double eff_size = 0.0;
        double pval = 1.0;
        if (config.test_type == "poisson_binom") {
            size_t n = 0;
            uint64_t s = 0;
            for (const auto &[j, c] : row) {
                n += static_cast<int64_t>(groups[j] != Group::OTHER) + (groups[j] == Group::BOTH);
                s += (groups[j] == Group::IN);
            }
            size_t front = n - std::min(n, num_labels_out);
            double min_pval = std::min(pb_pvals[n][front], pb_pvals[n].back());

            if (!config.test_by_unitig && min_pval * k_min >= config.family_wise_error_rate)
                return;

            pval = pb_pvals[n][s];
            eff_size = static_cast<double>(s) - mid_points[n];
        } else if (config.test_type == "poisson_exact" || config.test_type == "poisson_bayes") {
            double mu1 = 0.0;
            double mu2 = 0.0;
            if (config.test_type == "poisson_exact") {
                mu1 = static_cast<double>(in_kmers) / nelem;
                mu2 = static_cast<double>(out_kmers) / nelem;
            }
            int64_t n = 0;
            int64_t in_sum = 0;
            int64_t out_sum = 0;
            // sdsl::bit_vector found(num_labels_in + num_labels_out);
            for (const auto &[j, c] : row) {
                n += c;
                // double lambda_j = 0.0;
                // if (config.test_type == "poisson_bayes") {
                //     const auto &[r, p] = nb_params[j];
                //     lambda_j = (1.0 - p) * (r - 1.0 + c);
                //     found[j] = true;
                // }
                if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
                    // mu2 += lambda_j;
                    out_sum += c;
                }

                if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
                    // mu1 += lambda_j;
                    in_sum += c;
                }
            }
            if (config.test_type == "poisson_bayes") {
                size_t unitig_id = kmer_to_unitig[AnnotatedDBG::anno_to_graph_index(row_i)];
                auto &[unitig_size, unitig] = unitig_sums[unitig_id];
                for (size_t j = 0; j < num_labels_in + num_labels_out; ++j) {
                    int64_t unitig_sum = unitig[j];
                    const auto &[r, p] = nb_params[j];
                    double lambda_j = unitig_sum > 0 || r > 1
                        ? (r - 1 + unitig_sum) / (unitig_size + p / (1.0 - p))
                        : 0.0;
                    if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                        mu2 += lambda_j;

                    if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                        mu1 += lambda_j;
                }
                // for (size_t j = 0; j < found.size(); ++j) {
                //     if (!found[j]) {
                //         const auto &[r, p] = nb_params[j];
                //         if (r > 1.0) {
                //             double lambda_j = (1.0 - p) * (r - 1);
                //             if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                //                 mu2 += lambda_j;

                //             if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                //                 mu1 += lambda_j;
                //         }
                //     }
                // }
            }

            double p = mu1 / (mu1 + mu2);
            boost::math::binomial bdist(n, p);
            if (n != 0) {
                double pval0 = boost::math::pdf(bdist, 0);
                double pvaln = boost::math::pdf(bdist, n);
                double min_pval = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                if (!config.test_by_unitig && min_pval * k_min >= config.family_wise_error_rate)
                    return;

                auto get_deviance = [&](double y, double mu) {
                    y += 1e-8;
                    mu += 1e-8;
                    return 2 * (y * log(y/mu) - y + mu);
                };

                std::vector<double> devs;
                devs.reserve(n + 1);
                for (int64_t s = 0; s <= n; ++s) {
                    devs.emplace_back(get_deviance(s, mu1) + get_deviance(n - s, mu2));
                }

                pval = 0.0;

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

                if (!config.test_by_unitig && pval * k >= config.family_wise_error_rate)
                    return;
            }

            eff_size = static_cast<double>(in_sum) - boost::math::mode(bdist);
        } else if (config.test_type == "nbinom_exact") {
            int64_t n = 0;
            int64_t in_sum = 0;
            int64_t out_sum = 0;
            for (const auto &[j, c] : row) {
                n += c;
                if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
                    out_sum += c;
                }

                if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
                    in_sum += c;
                }
            }

            if (n > 0) {
                auto [r_n, p_n] = nb_params_null;
                auto [r_a, p_a] = nb_params_a;
                auto [r_b, p_b] = nb_params_b;

                double midpoint = r_a * n / (r_a + r_b);

                eff_size = static_cast<double>(in_sum) - midpoint;
                double in_dist = abs(eff_size);

                double base = nb_base + (lgamma(n + 1) - lgamma(r_n + n)) / log(2.0) - n*log2(1.0-p_n);
                auto get_pval = [&](int64_t s) {
                    int64_t t = n - s;
                    double sbase = (lgamma(r_a+s)+lgamma(r_b+t)-lgamma(s+1)-lgamma(t+1)) / log(2.0);
                    return exp2(base + sbase +s*log2(1.0-p_a) + t*log2(1.0-p_b));
                };

                {
                    double pval0 = get_pval(0);
                    double pvaln = get_pval(n);
                    double min_pval = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                    if (!config.test_by_unitig && min_pval * k >= config.family_wise_error_rate)
                        return;
                }

                pval = 0.0;
                for (int64_t s = 0; s <= n; ++s) {
                    double dist = abs(midpoint - s);
                    if (dist >= in_dist)
                        pval += get_pval(s);
                }
            }
        }

        set_pval(pval, eff_size);
    });

    if (config.test_by_unitig) {
        common::logger->trace("Combining p-values for unitigs");
        std::mutex m_vec_mu;
        clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
            double pval_min = 1.0;
            {
                size_t row = AnnotatedDBG::graph_to_anno_index(path[0]);
                size_t row_rank = kept.rank1(row);
                size_t bucket_idx = std::lower_bound(offsets.begin(), offsets.end(), row_rank) - offsets.begin() - 1;
                size_t offset = row_rank - offsets[bucket_idx];
                std::lock_guard<std::mutex> lock(m_vec_mu);
                pval_min = get<PValStorage>(pvals_min_buckets[bucket_idx][offset]);
            }

            if (pval_min * k_min >= config.family_wise_error_rate)
                return;

            std::vector<double> pvals;
            pvals.reserve(path.size());
            double comb_eff_size = 0.0;
            {
                std::lock_guard<std::mutex> lock(m_vec_mu);
                for (node_index node : path) {
                    size_t row = AnnotatedDBG::graph_to_anno_index(node);
                    size_t row_rank = kept.rank1(row);
                    size_t bucket_idx = std::lower_bound(offsets.begin(), offsets.end(), row_rank) - offsets.begin() - 1;
                    size_t offset = row_rank - offsets[bucket_idx];
                    pvals.emplace_back(get<PValStorage>(pvals_buckets[bucket_idx][offset]));
                    comb_eff_size += get<PValStorage>(eff_size_buckets[bucket_idx][offset]);
                }
            }
            double comp_pval = combine_pvals(pvals);
            if (comp_pval * k >= config.family_wise_error_rate || comb_eff_size == 0.0)
                return;

            auto &indicator = comb_eff_size > 0 ? indicator_in : indicator_out;
            for (node_index node : path) {
                set_bit(indicator.data(), node, true, std::memory_order_relaxed);
            }
        });

        pvals_min_buckets.resize(0);
        tmp_min_buckets.resize(0);
        pvals_buckets.resize(0);
        tmp_buckets.resize(0);
        eff_size_buckets.resize(0);
        tmp_eff_buckets.resize(0);
    }

    std::atomic_thread_fence(std::memory_order_acquire);

    // common::logger->trace("Computing minimum p-values");
    // std::vector<std::pair<double, size_t>> m;
    // std::vector<int64_t> m_sums;
    // {
    //     std::vector<std::vector<std::pair<double, size_t>>> ms(num_threads + 1);
    //     std::vector<std::vector<tsl::hopscotch_map<std::vector<int64_t>, std::pair<double, size_t>, utils::VectorHash>>> ms_vecs(num_threads + 1);
    //     generate_rows([&](uint64_t row_i, const auto &raw_row, size_t bucket_idx) {
    //         int64_t in_sum = 0;
    //         int64_t out_sum = 0;
    //         int64_t out_stat_int = 0;
    //         bool in_kmer = false;
    //         bool out_kmer = false;
    //         PairContainer row;
    //         std::vector<int64_t> vals;
    //         if (kept[row_i]) {
    //             size_t count_in = 0;
    //             size_t count_out = 0;
    //             size_t total_count = 0;
    //             for (const auto &[j, raw_c] : raw_row) {
    //                 if (min_counts.size() && raw_c < min_counts[j])
    //                     continue;

    //                 uint64_t c = raw_c;
    //                 if (count_maps.size())
    //                     c = count_maps[j].find(c)->second.first;

    //                 row.emplace_back(j, c);
    //                 vals.emplace_back(c);
    //                 if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
    //                     out_sum += c;
    //                     ++count_out;
    //                 }

    //                 if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
    //                     in_sum += c;
    //                     ++count_in;
    //                 }

    //                 total_count += groups[j] != Group::OTHER;
    //             }

    //             if (total_count >= config.min_recurrence) {
    //                 double out_stat = static_cast<double>(out_sum) / out_kmers * in_kmers;
    //                 if (out_stat != 0)
    //                     out_stat_int = out_stat > 0 ? ceil(out_stat) : floor(out_stat);

    //                 in_kmer = count_in >= config.min_in_recurrence && count_out <= config.max_out_recurrence && in_sum > (out_kmers > 0 ? out_stat : 0.0);
    //                 out_kmer = count_out >= config.min_out_recurrence && count_in <= config.max_in_recurrence && in_sum < out_stat;
    //             }
    //         } else {
    //             return;
    //         }

    //         size_t n = config.test_type != "poisson_binom" ? in_sum + out_sum : row.size();
    //         if (!in_kmer && !out_kmer) {
    //             row.clear();
    //             vals.clear();
    //             n = 0;
    //         }

    //         double pval_min = 0;
    //         if (config.test_type != "gnb_exact" && config.test_type != "lnb_exact") {
    //             assert(bucket_idx < ms.size());
    //             if (n >= ms[bucket_idx].size())
    //                 ms[bucket_idx].resize(n + 1, std::make_pair(1.1, 0));

    //             if (ms[bucket_idx][n].first == 1.1)
    //                 ms[bucket_idx][n].first = compute_min_pval(n, row);

    //             pval_min = ms[bucket_idx][n].first;
    //             ++ms[bucket_idx][n].second;
    //         } else {
    //             std::sort(vals.begin(), vals.end());
    //             if (n >= ms_vecs[bucket_idx].size())
    //                 ms_vecs[bucket_idx].resize(n + 1);

    //             auto find = ms_vecs[bucket_idx][n].find(vals);
    //             if (find == ms_vecs[bucket_idx][n].end()) {
    //                 find = ms_vecs[bucket_idx][n].try_emplace(
    //                     vals,
    //                     std::make_pair(compute_min_pval(n, row), 1)
    //                 ).first;
    //             } else {
    //                 ++find.value().second;
    //             }

    //             pval_min = find->second.first;
    //         }

    //         if (config.test_by_unitig) {
    //             bucket_idx = std::min(bucket_idx, pvals_buckets.size() - 1);
    //             uint64_t eff_size = bit_cast<uint64_t, int64_t>(in_sum - out_stat_int);
    //             if constexpr(preallocated) {
    //                 node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
    //                 eff_size_buckets[bucket_idx][node] = eff_size;
    //                 sum_buckets[bucket_idx][node] = n;
    //                 if (config.test_type == "gnb_exact" || config.test_type == "lnb_exact")
    //                     pvals_min_buckets[bucket_idx][node] = bit_cast<uint64_t>(pval_min);
    //             } else {
    //                 eff_size_buckets[bucket_idx].push_back(eff_size);
    //                 sum_buckets[bucket_idx].push_back(n);
    //                 if (config.test_type == "gnb_exact" || config.test_type == "lnb_exact")
    //                     pvals_min_buckets[bucket_idx].push_back(bit_cast<uint64_t>(pval_min));
    //             }
    //         }
    //     });

    //     common::logger->trace("Merging min p-value tables");
    //     if (config.test_type != "gnb_exact" && config.test_type != "lnb_exact") {
    //         for (size_t i = 1; i < ms.size(); ++i) {
    //             size_t end = std::min(ms[0].size(), ms[i].size());
    //             for (size_t j = 0; j < end; ++j) {
    //                 assert(ms[0][j].first == 1.1 || ms[i][j].first == 1.1 || ms[0][j] == ms[i][j]);
    //                 ms[0][j].first = std::min(ms[0][j].first, ms[i][j].first);
    //                 ms[0][j].second += ms[i][j].second;
    //             }

    //             if (ms[i].size() > ms[0].size()) {
    //                 std::copy(std::make_move_iterator(ms[i].begin() + ms[0].size()),
    //                         std::make_move_iterator(ms[i].end()),
    //                         std::back_inserter(ms[0]));
    //             }
    //         }
    //         m = std::move(ms[0]);
    //         ms.resize(0);
    //     } else {
    //         std::vector<std::tuple<double, int64_t, size_t>> mn;
    //         for (size_t i = 0; i < ms_vecs.size(); ++i) {
    //             for (size_t n = 0; n < ms_vecs[i].size(); ++n) {
    //                 for (const auto &[vec, v] : ms_vecs[i][n]) {
    //                     mn.emplace_back(v.first, n, v.second);
    //                 }
    //             }
    //         }
    //         ms_vecs.resize(0);
    //         std::sort(mn.begin(), mn.end(), utils::GreaterFirst());

    //         m.reserve(mn.size());
    //         m_sums.reserve(mn.size());
    //         for (const auto &[p, s, c] : mn) {
    //             m.emplace_back(p, c);
    //             m_sums.emplace_back(s);
    //         }
    //     }
    // }

    // size_t k = 0;
    // size_t n_cutoff = 1;

    // auto combine_pvals = [dist=boost::math::cauchy()](const std::vector<double> &pvals) {
    //     double stat = 0.0;
    //     for (double pval : pvals) {
    //         stat += tan((0.5 - pval) * M_PI);
    //     }
    //     stat /= pvals.size();
    //     return boost::math::cdf(boost::math::complement(dist, stat));
    // };

    // std::vector<size_t> boundaries;
    // if (config.test_by_unitig) {
    //     boundaries.reserve(eff_size_buckets.size() + 1);
    //     boundaries.emplace_back(0);
    //     for (size_t i = 0; i < eff_size_buckets.size(); ++i) {
    //         boundaries.emplace_back(boundaries.back() + eff_size_buckets[i].size());
    //     }
    // }

    // std::mutex pval_mu;
    // std::unique_ptr<MaskedDeBruijnGraph> clean_masked_graph;
    // sdsl::bit_vector unitig_start;
    // if (!config.test_by_unitig) {
    //     std::tie(k, n_cutoff) = correct_pvals(m, nelem, m_sums);
    //     common::logger->trace("Picked: k: {}\tn: {}\tpval_min: {}", k, n_cutoff, config.family_wise_error_rate / k);
    // } else {
    //     clean_masked_graph = std::make_unique<MaskedDeBruijnGraph>(
    //         graph_ptr,
    //         [&](node_index node) {
    //             return node != DeBruijnGraph::npos
    //                     && kept[AnnotatedDBG::graph_to_anno_index(node)];
    //         },
    //         true,
    //         is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
    //     );

    //     std::vector<std::pair<double, int64_t>> comb_msums;
    //     common::logger->trace("Allocating updated indicator");
    //     unitig_start = sdsl::bit_vector(graph_ptr->max_index() + 1, false);

    //     common::logger->trace("Combining min p-vals");
    //     clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
    //         std::vector<size_t> bucket_idxs;
    //         bucket_idxs.reserve(path.size());
    //         for (node_index node : path) {
    //             bucket_idxs.emplace_back(0);
    //             for (size_t i = 1; i < boundaries.size(); ++i) {
    //                 if (boundaries[i] > node)
    //                     break;

    //                 ++bucket_idxs.back();
    //             }
    //         }

    //         std::vector<double> pvals_min;
    //         pvals_min.reserve(path.size());

    //         int64_t n_sum = 0;
    //         {
    //             std::lock_guard<std::mutex> lock(pval_mu);
    //             for (size_t i = 0; i < path.size(); ++i) {
    //                 size_t bucket_idx = bucket_idxs[i];
    //                 size_t node_shift = path[i] - boundaries[bucket_idx];
    //                 size_t n = sum_buckets[bucket_idx][node_shift];
    //                 n_sum += n;

    //                 if (config.test_type == "gnb_exact" || config.test_type == "lnb_exact") {
    //                     pvals_min.emplace_back(bit_cast<double, uint64_t>(pvals_min_buckets[bucket_idx][node_shift]));
    //                 } else {
    //                     pvals_min.emplace_back(m[n].first);
    //                 }
    //             }
    //         }

    //         double comb_pval = *std::min_element(pvals_min.begin(), pvals_min.end());
    //         if (comb_pval < config.family_wise_error_rate) {
    //             for (node_index node : path) {
    //                 set_bit(unitig_start.data(), node, parallel, MO_RELAXED);
    //             }
    //         }

    //         std::lock_guard<std::mutex> lock(pval_mu);
    //         comb_msums.emplace_back(comb_pval, n_sum);
    //     }, num_threads);

    //     std::sort(comb_msums.begin(), comb_msums.end(), utils::GreaterFirst());
    //     std::vector<std::pair<double, size_t>> comb_m;
    //     std::vector<int64_t> comb_m_sums;
    //     comb_m.reserve(comb_msums.size());
    //     comb_m_sums.reserve(comb_msums.size());
    //     for (const auto &[p, s] : comb_msums) {
    //         comb_m.emplace_back(p, 1);
    //         comb_m_sums.emplace_back(s);
    //     }

    //     std::tie(k, n_cutoff) = correct_pvals(comb_m, nelem, comb_m_sums);
    //     common::logger->trace("Picked: k: {}\tpval_min: {}", k, config.family_wise_error_rate / k);

    //     sum_buckets.resize(0);
    //     tmp_sum_buckets.resize(0);
    //     pvals_min_buckets.resize(0);
    //     tmp_min_buckets.resize(0);
    // }

    // if (k > 0) {
    //     common::logger->trace("Running differential tests");
    //     std::exception_ptr ex = nullptr;
    //     std::mutex ex_mu;
    //     bool keep_all_pvals = config.output_pvals || config.test_by_unitig;
    //     std::atomic_thread_fence(std::memory_order_release);

    //     generate_rows([&](uint64_t row_i, const auto &raw_row, size_t bucket_idx) {
    //         if (ex)
    //             return;

    //         int64_t in_sum = 0;
    //         int64_t out_sum = 0;
    //         bool in_kmer = false;
    //         bool out_kmer = false;
    //         PairContainer row;
    //         bucket_idx = std::min(bucket_idx, pvals_buckets.size() - 1);
    //         node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
    //         double pval = 1.1;
    //         if (kept[row_i]) {
    //             size_t count_in = 0;
    //             size_t count_out = 0;
    //             size_t total_count = 0;
    //             for (const auto &[j, raw_c] : raw_row) {
    //                 if (min_counts.size() && raw_c < min_counts[j])
    //                     continue;

    //                 uint64_t c = raw_c;
    //                 if (count_maps.size())
    //                     c = count_maps[j].find(c)->second.first;

    //                 row.emplace_back(j, c);
    //                 if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
    //                     out_sum += c;
    //                     ++count_out;
    //                 }

    //                 if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
    //                     in_sum += c;
    //                     ++count_in;
    //                 }

    //                 total_count += groups[j] != Group::OTHER;
    //             }

    //             if (total_count >= config.min_recurrence) {
    //                 double out_stat = static_cast<double>(out_sum) / out_kmers * in_kmers;
    //                 in_kmer = count_in >= config.min_in_recurrence && count_out <= config.max_out_recurrence && in_sum > (out_kmers > 0 ? out_stat : 0.0);
    //                 out_kmer = count_out >= config.min_out_recurrence && count_in <= config.max_in_recurrence && in_sum < out_stat;
    //             }
    //         } else {
    //             if (keep_all_pvals) {
    //                 auto &pvals = pvals_buckets[bucket_idx];
    //                 if constexpr(preallocated) {
    //                     pvals[node] = bit_cast<uint64_t>(pval);
    //                 } else {
    //                     pvals.push_back(bit_cast<uint64_t>(pval));
    //                 }
    //             }
    //             return;
    //         }

    //         size_t n = config.test_type != "poisson_binom" ? in_sum + out_sum : row.size();

    //         if (!in_kmer && !out_kmer) {
    //             row.clear();
    //             n = 0;
    //         }

    //         if (!config.output_pvals && config.test_by_unitig && !unitig_start[node]) {
    //             auto &pvals = pvals_buckets[bucket_idx];
    //             if constexpr(preallocated) {
    //                 pvals[node] = bit_cast<uint64_t>(pval);
    //             } else {
    //                 pvals.push_back(bit_cast<uint64_t>(pval));
    //             }
    //         }

    //         if (config.output_pvals
    //                 || (!config.test_by_unitig && n >= n_cutoff)
    //                 || (config.test_by_unitig && unitig_start[node])) {
    //             if (config.test_type != "gnb_exact" && config.test_type != "lnb_exact" && m[n].first == 1.1) {
    //                 common::logger->error("in: {}\tout: {}\tn: {}", in_sum, out_sum, n);
    //                 throw std::runtime_error("Indexing invalid");
    //             }

    //             try {
    //                 pval = compute_pval(in_sum, out_sum, row);
    //             } catch (...) {
    //                 std::lock_guard<std::mutex> lock(ex_mu);
    //                 ex = std::current_exception();
    //                 return;
    //             }

    //             if (config.test_type != "gnb_exact" && config.test_type != "lnb_exact" && pval < 1.1 && m[n].first - pval > 1e-10) {
    //                 common::logger->error("Min p-val estimate too high: min {} > cur {}\tn: {}\ttest: {}", m[n].first, pval, n, config.test_type);
    //                 throw std::runtime_error("Test failed");
    //             }

    //             if (config.test_type != "notest" && keep_all_pvals) {
    //                 auto &pvals = pvals_buckets[bucket_idx];
    //                 if constexpr(preallocated) {
    //                     pvals[node] = bit_cast<uint64_t>(pval);
    //                 } else {
    //                     pvals.push_back(bit_cast<uint64_t>(pval));
    //                 }
    //             }

    //             if (!config.test_by_unitig && in_kmer != out_kmer && pval * k < config.family_wise_error_rate) {
    //                 bool use_atomic = parallel && (node % 64 == 0);
    //                 set_bit((in_kmer ? indicator_in : indicator_out).data(), node, use_atomic, MO_RELAXED);
    //             }
    //         }
    //     });

    //     std::atomic_thread_fence(std::memory_order_acquire);

    //     if (ex)
    //         std::rethrow_exception(ex);

    //     unitig_start = sdsl::bit_vector();
    //     deallocate();

    //     if (!config.test_by_unitig) {
    //         boundaries.reserve(pvals_buckets.size() + 1);
    //         boundaries.emplace_back(0);
    //         for (size_t i = 0; i < pvals_buckets.size(); ++i) {
    //             boundaries.emplace_back(boundaries.back() + pvals_buckets[i].size());
    //         }
    //     }

    //     if (config.test_by_unitig) {
    //         common::logger->trace("Allocating k-mer bitmasks");
    //         indicator_in = sdsl::bit_vector(graph_ptr->max_index() + 1, false);
    //         indicator_out = sdsl::bit_vector(graph_ptr->max_index() + 1, false);

    //         common::logger->trace("Combining p-values within unitigs");
    //         std::atomic<uint64_t> num_sig_unitigs{0};

    //         std::atomic_thread_fence(std::memory_order_release);
    //         clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
    //             int64_t comb_eff_size = 0;
    //             std::vector<size_t> bucket_idxs;
    //             bucket_idxs.reserve(path.size());
    //             for (node_index node : path) {
    //                 bucket_idxs.emplace_back(0);
    //                 for (size_t i = 1; i < boundaries.size(); ++i) {
    //                     if (boundaries[i] > node)
    //                         break;

    //                     ++bucket_idxs.back();
    //                 }
    //             }

    //             // TODO: don't bother combining if the row sum is too small
    //             std::vector<double> pvals;
    //             pvals.reserve(path.size());
    //             {
    //                 std::lock_guard<std::mutex> lock(pval_mu);
    //                 for (size_t i = 0; i < path.size(); ++i) {
    //                     node_index node = path[i];
    //                     size_t bucket_idx = bucket_idxs[i];
    //                     size_t node_shift = node - boundaries[bucket_idx];

    //                     comb_eff_size += bit_cast<int64_t, uint64_t>(eff_size_buckets[bucket_idx][node_shift]);

    //                     double pval = bit_cast<double, uint64_t>(pvals_buckets[bucket_idx][node_shift]);
    //                     pvals.emplace_back(pval);
    //                 }
    //             }
    //             double comb_pval = combine_pvals(pvals);
    //             if (comb_pval * k < config.family_wise_error_rate && comb_eff_size != 0) {
    //                 ++num_sig_unitigs;
    //                 auto *data = comb_eff_size > 0 ? indicator_in.data() : indicator_out.data();
    //                 for (node_index node : path) {
    //                     set_bit(data, node, parallel, MO_RELAXED);
    //                 }
    //             }

    //             uint64_t comb_pval_enc = bit_cast<uint64_t, double>(comb_pval);

    //             std::lock_guard<std::mutex> lock(pval_mu);
    //             for (size_t i = 0; i < path.size(); ++i) {
    //                 size_t bucket_idx = bucket_idxs[i];
    //                 pvals_buckets[bucket_idx][path[i] - boundaries[bucket_idx]] = comb_pval_enc;
    //             }
    //         }, num_threads);
    //         std::atomic_thread_fence(std::memory_order_acquire);
    //         common::logger->trace("Found {} significant unitigs before correction", num_sig_unitigs);

    //         eff_size_buckets.resize(0);
    //         tmp_eff_buckets.resize(0);
    //     }
    // } else {
    //     deallocate();
    // }

    std::unique_ptr<utils::TempFile> tmp_file;
    PValStorage pvals;

    // if constexpr(std::is_same_v<PValStorage, std::vector<uint64_t>>) {
    //     pvals = std::move(pvals_buckets[0]);
    // }

    // if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
    //     if (config.output_pvals) {
    //         common::logger->trace("Merging buckets");
    //         for (size_t i = 1; i < pvals_buckets.size(); ++i) {
    //             for (size_t j = 0; j < pvals_buckets[i].size(); ++j) {
    //                 pvals_buckets[0].push_back(pvals_buckets[i][j]);
    //             }
    //             pvals_buckets[i].close();
    //         }
    //     }
    //     pvals = std::move(pvals_buckets[0]);
    //     pvals_buckets.resize(0);
    //     tmp_file = std::move(tmp_buckets[0]);
    //     tmp_buckets.resize(0);
    // }

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

    size_t total_labels = 0;
    for (const auto &file : files) {
        total_labels += annot::ColumnCompressed<>::read_num_labels(file);
    }

    using ColumnValuesMem = std::vector<std::unique_ptr<const sdsl::int_vector<>>>;
    using ColumnValuesDisk = std::vector<std::unique_ptr<const sdsl::int_vector_buffer<>>>;

    std::variant<ColumnValuesMem, ColumnValuesDisk> column_values_all;

    if (config.test_by_unitig) {
        column_values_all = ColumnValuesMem();
    } else {
        column_values_all = ColumnValuesDisk();
    }

    return std::visit([&](auto&& column_values_all) {
        using ColumnValues = typename std::decay<decltype(column_values_all)>::type;
        using ValuesContainerPtr = typename ColumnValues::value_type;
        using ValuesContainer = typename std::decay<typename ValuesContainerPtr::element_type>::type;
        using value_type = typename ValuesContainer::value_type;
        using PairContainer = Vector<std::pair<uint64_t, value_type>>;

        std::vector<std::unique_ptr<const bit_vector>> columns_all(total_labels);
        std::vector<Group> groups(total_labels);
        if (files.empty()) {
            throw std::runtime_error("No files provided");
        }

        if (std::filesystem::exists(files[0] + ".counts")) {
            column_values_all.resize(total_labels);
            annot::ColumnCompressed<>::load_columns_and_values(
                files,
                [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector>&& column, ValuesContainer&& column_values) {
                    bool is_in = labels_in.count(label);
                    bool is_out = labels_out.count(label);
                    if (is_in) {
                        groups[offset] = !is_out ? Group::IN : Group::BOTH;
                    } else {
                        groups[offset] = is_out ? Group::OUT : Group::OTHER;
                    }
                    columns_all[offset].reset(column.release());
                    column_values_all[offset] = std::make_unique<ValuesContainer>(std::move(column_values));
                },
                num_parallel_files
            );
        } else {
            annot::ColumnCompressed<>::merge_load(
                files,
                [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector>&& column) {
                    bool is_in = labels_in.count(label);
                    bool is_out = labels_out.count(label);
                    if (is_in) {
                        groups[offset] = !is_out ? Group::IN : Group::BOTH;
                    } else {
                        groups[offset] = is_out ? Group::OUT : Group::OTHER;
                    }
                    columns_all[offset].reset(column.release());
                },
                num_parallel_files
            );
        }

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

    std::vector<Group> groups;
    groups.reserve(total_labels);
    for (const auto &label : annotation.get_label_encoder().get_labels()) {
        bool is_in = labels_in.count(label);
        bool is_out = labels_out.count(label);
        if (is_in) {
            groups.emplace_back(!is_out ? Group::IN : Group::BOTH);
        } else {
            groups.emplace_back(is_out ? Group::OUT : Group::OTHER);
        }
    }

    auto graph_ptr = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());

    const auto &matrix = annotation.get_matrix();
    using Row = annot::matrix::BinaryMatrix::Row;

    if (const auto *int_matrix = dynamic_cast<const annot::matrix::IntMatrix*>(&matrix)) {
        using value_type = annot::matrix::IntMatrix::Value;
        return mask_nodes_by_label_dual<value_type, PValStorage>(
            graph_ptr,
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
            size_t n = matrix.num_rows();
            size_t batch_size = (n + num_threads - 1) / num_threads;
            size_t rows_per_update = 10000;
            ProgressBar progress_bar(n, "Streaming rows", std::cerr, !common::get_verbose());
            #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
            for (size_t k = 0; k < n; k += batch_size) {
                size_t begin = k;
                size_t end = std::min(begin + batch_size, n);
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
            graph_ptr,
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


