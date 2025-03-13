#include "annotated_graph_algorithm.hpp"

#include <typeinfo>
#include <mutex>
#include <variant>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/hypergeometric.hpp>

#include <sdust.h>

#include "common/logger.hpp"
#include "common/vectors/bitmap.hpp"
#include "common/vector_map.hpp"
#include "common/vector_set.hpp"
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

enum Group { IN, OUT, OTHER, BOTH };


template <class To, class From>
std::enable_if_t<sizeof(To) == sizeof(From) && std::is_trivially_copyable_v<From>
                         && std::is_trivially_copyable_v<To>,
                 To>
// constexpr support needs compiler magic
bit_cast(const From& src) noexcept {
    static_assert(std::is_trivially_constructible_v<To>,
                  "This implementation additionally requires "
                  "destination type to be trivially constructible");

    static_assert(sizeof(To) == sizeof(From));
    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

#if !_PROTEIN_GRAPH
inline bool is_low_complexity(std::string_view s, int T = 20, int W = 64) {
    int n;
    std::unique_ptr<uint64_t, decltype(std::free)*> r {
        sdust(0, (const uint8_t*)s.data(), s.size(), T, W, &n), std::free
    };
    return n > 0;
}
#else
inline bool is_low_complexity(std::string_view, int = 20, int = 64) {
    // TODO: implement a checker here
    return false;
}
#endif

template <class PValStorage>
double get(typename PValStorage::reference&& ref) {
    return bit_cast<double>(static_cast<typename PValStorage::value_type>(ref));
}

template <class PValStorage>
void set(typename PValStorage::reference&& ref, double val) {
    ref = bit_cast<typename PValStorage::value_type>(val);
}

template <class PValStorage>
void push_back(PValStorage& v, double val) {
    v.push_back(bit_cast<typename PValStorage::value_type>(val));
}

template <typename Generator>
std::tuple<long double, long double, long double, long double>
get_rp(const Generator& generate) {
    long double mu = 0;
    long double var = 0;
    size_t total = 0;
    generate([&](auto k, auto c) {
        mu += k * c;
        var += k * k * c;
        total += c;
    });
    mu /= total;
    long double mu2 = mu * mu;
    var = (var - mu2 * total) / (total - 1);

    if (mu >= var) {
        common::logger->warn("Fit failed, falling back to Poisson: mu: {} >= var: {}", mu,
                             var);
        return { 0.0, 1.0, mu, var };
    }

    auto get_l = [&](long double r) {
        long double p = r / (r + mu);
        long double l = -lgammal(r) * total + mu * total * log1pl(-p) + total * r * logl(p);
        generate([&](auto k, auto c) { l += lgammal(k + r) * c; });
        return l;
    };

    auto get_dl_r = [&](long double r) {
        long double dl = (log(r) - log(r + mu) - boost::math::digamma(r)) * total;
        generate([&](auto k, auto c) { dl += boost::math::digamma(k + r) * c; });
        return dl;
    };

    auto get_dl_a = [&](long double a) {
        long double dl = mu * mu * total / (1.0 + a * mu) + total * logl(1 + a * mu) / a / a
                + total / a / a * boost::math::digamma(1.0 / a);
        generate([&](auto k, auto c) {
            dl -= boost::math::digamma(k + 1.0 / a) / a / a * c;
        });
        return dl;
    };

    long double r1 = 0.0;
    long double r2 = 0.0;
    try {
        auto [r_min, r_max] = boost::math::tools::bisect(
                get_dl_r, std::numeric_limits<double>::min(), 1.0,
                boost::math::tools::eps_tolerance<long double>(5));
        r1 = (r_min + r_max) / 2.0;
    } catch (boost::wrapexcept<boost::math::evaluation_error>& e) {
    }

    try {
        auto [a_min, a_max] = boost::math::tools::bisect(
                get_dl_a, std::numeric_limits<double>::min(), 1.0,
                boost::math::tools::eps_tolerance<long double>(5));
        r2 = (1.0 / a_min + 1.0 / a_max) / 2.0;
    } catch (boost::wrapexcept<boost::math::evaluation_error>& e) {
    }

    long double r = 0.0;
    if (r1 == 0 && r2 == 0) {
        throw std::runtime_error("Failed to fit");
    } else if (r1 == 0) {
        r = r2;
    } else if (r2 == 0) {
        r = r1;
    } else {
        r = get_l(r1) > get_l(r2) ? r1 : r2;
    }

    long double p = r / (r + mu);
    // long double p = 0.1;
    // r = p * mu / (1.0 - p);
    return { r, p, mu, var };
};

template <class G1, class G2>
std::pair<long double, long double> mann_whitneyu(const G1& generate_a, const G2& generate_b) {
    std::vector<long double> x;
    std::vector<long double> y;
    std::vector<long double> xy;
    generate_a([&](auto c) {
        x.emplace_back(c);
        xy.emplace_back(c);
    });
    generate_b([&](auto c) {
        y.emplace_back(c);
        xy.emplace_back(c);
    });

    long double n1 = x.size();
    long double n2 = y.size();
    VectorMap<long double, size_t> counts;
    for (long double c : xy) {
        ++counts[c];
    }
    if (counts.size() <= 1)
        return std::make_pair(1.0, 0.0);

    std::vector<std::pair<long double, size_t>> ranked_counts = counts.values_container();
    std::sort(ranked_counts.begin(), ranked_counts.end());
    tsl::hopscotch_map<long double, long double> count_to_rank;
    size_t cur_rank = 0;
    for (const auto& [v, c] : ranked_counts) {
        size_t last_rank = cur_rank;
        cur_rank += c;
        long double r = 0.0;
        for (size_t i = last_rank + 1; i <= cur_rank; ++i) {
            r += i;
        }
        r /= c;
        count_to_rank[v] = r;
    }

    long double R1 = 0.0;
    for (long double v : x) {
        R1 += count_to_rank[v];
    }
    long double R2 = 0.0;
    for (long double v : y) {
        R2 += count_to_rank[v];
    }

    long double U1 = R1 - n1 * (n1 + 1) / 2.0;
    long double U2 = n1 * n2 - U1;
    long double U = std::max(U1, U2);

    long double mu = n1 * n2 / 2.0;
    long double n = n1 + n2;
    long double tie = 0.0;
    for (const auto& [v, c] : ranked_counts) {
        tie += c * c * c - c;
    }
    long double s = sqrtl(n1 * n2 / 12 * ((n + 1) - tie / (n * (n - 1))));
    long double numerator = U - mu - 0.5;
    long double z = numerator / s;
    boost::math::normal dist;
    long double p = boost::math::cdf(boost::math::complement(dist, z)) * 2;
    p = std::min(p, 1.0L);

    long double eff_size = R1 / n1 - R2 / n2;
    return std::make_pair(p, eff_size);
}

template <class G1, class G2>
std::pair<long double, long double> brunner_munzel(const G1& generate_a,
                                                   const G2& generate_b) {
    std::vector<long double> x;
    std::vector<long double> y;
    std::vector<long double> xy;
    generate_a([&](auto c) {
        x.emplace_back(c);
        xy.emplace_back(c);
    });
    generate_b([&](auto c) {
        y.emplace_back(c);
        xy.emplace_back(c);
    });

    long double nx = x.size();
    long double ny = y.size();
    VectorMap<long double, size_t> counts;
    for (long double c : xy) {
        ++counts[c];
    }
    if (counts.size() <= 1)
        return std::make_pair(1.0, 0.0);

    std::vector<std::pair<long double, size_t>> ranked_counts = counts.values_container();
    std::sort(ranked_counts.begin(), ranked_counts.end());
    tsl::hopscotch_map<long double, long double> count_to_rank;
    size_t cur_rank = 0;
    for (const auto& [v, c] : ranked_counts) {
        size_t last_rank = cur_rank;
        cur_rank += c;
        long double r = 0.0;
        for (size_t i = last_rank + 1; i <= cur_rank; ++i) {
            r += i;
        }
        r /= c;
        count_to_rank[v] = r;
    }

    std::vector<long double> rankc;
    rankc.reserve(xy.size());
    for (long double v : xy) {
        rankc.emplace_back(count_to_rank[v]);
    }

    std::vector<long double> rankcx;
    rankcx.reserve(x.size());
    for (long double v : x) {
        rankcx.emplace_back(count_to_rank[v]);
    }
    std::vector<long double> rankcy;
    rankcy.reserve(y.size());
    for (long double v : y) {
        rankcy.emplace_back(count_to_rank[v]);
    }

    long double rankcx_mean = std::accumulate(rankcx.begin(), rankcx.end(), 0.0L) / nx;
    long double rankcy_mean = std::accumulate(rankcy.begin(), rankcy.end(), 0.0L) / ny;

    counts.clear();
    for (long double c : x) {
        ++counts[c];
    }
    ranked_counts = counts.values_container();
    std::sort(ranked_counts.begin(), ranked_counts.end());
    count_to_rank.clear();
    cur_rank = 0;
    for (const auto& [v, c] : ranked_counts) {
        size_t last_rank = cur_rank;
        cur_rank += c;
        long double r = 0.0;
        for (size_t i = last_rank + 1; i <= cur_rank; ++i) {
            r += i;
        }
        r /= c;
        count_to_rank[v] = r;
    }
    std::vector<long double> rankx;
    rankx.reserve(x.size());
    for (long double v : x) {
        rankx.emplace_back(count_to_rank[v]);
    }


    counts.clear();
    for (long double c : y) {
        ++counts[c];
    }
    ranked_counts = counts.values_container();
    std::sort(ranked_counts.begin(), ranked_counts.end());
    count_to_rank.clear();
    cur_rank = 0;
    for (const auto& [v, c] : ranked_counts) {
        size_t last_rank = cur_rank;
        cur_rank += c;
        long double r = 0.0;
        for (size_t i = last_rank + 1; i <= cur_rank; ++i) {
            r += i;
        }
        r /= c;
        count_to_rank[v] = r;
    }
    std::vector<long double> ranky;
    ranky.reserve(y.size());
    for (long double v : x) {
        ranky.emplace_back(count_to_rank[v]);
    }

    long double rankx_mean = std::accumulate(rankx.begin(), rankx.end(), 0.0L) / nx;
    long double ranky_mean = std::accumulate(ranky.begin(), ranky.end(), 0.0L) / ny;

    long double Sx = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        Sx += pow(rankcx[i] - rankx[i] - rankcx_mean + rankx_mean, 2.0);
    }
    Sx /= nx - 1;
    long double Sy = 0.0;
    for (size_t i = 0; i < y.size(); ++i) {
        Sy += pow(rankcy[i] - ranky[i] - rankcy_mean + ranky_mean, 2.0);
    }
    Sy /= ny - 1;

    long double wbfn = nx * ny * (rankcy_mean - rankcx_mean);

    if (wbfn == 0.0)
        return std::make_pair(1.0, 0.0);

    long double wbfn_denom = (nx + ny) * sqrtl(nx * Sx + ny * Sy);
    if (wbfn_denom == 0.0)
        return std::make_pair(0.0, -wbfn);

    wbfn /= wbfn_denom;
    long double p = 1.0;
    long double eff_size = -wbfn;

    if (nx + ny >= 50) {
        p = boost::math::cdf(boost::math::complement(boost::math::normal(), abs(wbfn))) * 2.0;
    } else if (nx + ny >= 20) {
        long double df_denom = pow(nx * Sx, 2.0) / (nx - 1);
        df_denom += pow(ny * Sy, 2.0) / (ny - 1);
        if (df_denom == 0.0) {
            p = boost::math::cdf(boost::math::complement(boost::math::normal(), abs(wbfn)))
                    * 2.0;
        } else {
            long double df_numer = pow(nx * Sx + ny * Sy, 2.0);
            long double df = df_numer / df_denom;
            p = boost::math::cdf(
                        boost::math::complement(boost::math::students_t(df), abs(wbfn)))
                    * 2.0;
        }
    } else {
        throw std::runtime_error("Too few samples");
    }

    return std::make_pair(p, eff_size);
}

template <typename value_type, class PValStorage, typename HistGetter, typename Generator>
std::tuple<std::shared_ptr<DeBruijnGraph>,
           std::shared_ptr<DeBruijnGraph>,
           PValStorage,
           std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(
        std::shared_ptr<const DeBruijnGraph> graph_ptr,
        const HistGetter& get_hist_map,
        const Generator& generate_rows,
        const std::vector<Group>& groups,
        const DifferentialAssemblyConfig& config,
        size_t num_threads = 1,
        std::filesystem::path tmp_dir = "",
        size_t num_parallel_files = std::numeric_limits<size_t>::max(),
        const std::function<void()>& deallocate = []() {},
        uint8_t max_width = 64) {
    if (auto canonical = std::dynamic_pointer_cast<const CanonicalDBG>(graph_ptr))
        graph_ptr = canonical->get_graph_ptr();

    num_parallel_files = get_num_threads();
    num_threads = get_num_threads();
    bool is_primary = graph_ptr->get_mode() == DeBruijnGraph::PRIMARY;

    auto combine_pvals
            = [dist = boost::math::cauchy()](const std::vector<long double>& pvals) {
                  long double stat = 0.0;
                  for (long double pval : pvals) {
                      stat += tan((0.5 - pval) * M_PI);
                  }
                  stat /= pvals.size();
                  return boost::math::cdf(boost::math::complement(dist, stat));
              };

    size_t num_labels_both = std::count(groups.begin(), groups.end(), Group::BOTH);
    size_t num_labels_out = std::count(groups.begin(), groups.end(), Group::OUT);
    size_t num_labels_in = std::count(groups.begin(), groups.end(), Group::IN);
    num_labels_in += num_labels_both;
    num_labels_out += num_labels_both;

    common::logger->trace("Graph mode: {}", is_primary ? "PRIMARY" : "other");

    sdsl::bit_vector kept_bv(AnnotatedDBG::graph_to_anno_index(graph_ptr->max_index() + 1),
                             true);

    std::vector<size_t> min_counts(groups.size(), config.min_count);
    std::vector<uint64_t> check_cutoff(groups.size(), std::numeric_limits<uint64_t>::max());

    // std::mutex agg_mu;
    // using PairContainer = std::vector<std::pair<uint64_t, value_type>>;

    if (config.clean) {
        common::logger->trace("Computing histogram from uncleaned counts");
        auto hists_map = get_hist_map(std::vector<size_t>(groups.size(), 1), nullptr);
        common::logger->trace("Cleaning count columns");

#pragma omp parallel for num_threads(num_parallel_files)
        for (size_t j = 0; j < groups.size(); ++j) {
            if (groups[j] == Group::OTHER) {
                min_counts[j] = std::numeric_limits<size_t>::max();
                continue;
            }

            // set cutoff for lower end of distribution
            auto [mean_est, nzeros_est] = estimate_ztp_mean(
                    [&](const auto& callback) {
                        for (const auto& [k, c] : hists_map[j]) {
                            if (k > 0 && k <= N_BUCKETS_FOR_ESTIMATION)
                                callback(k, c);
                        }
                    },
                    max_width, 0, N_BUCKETS_FOR_ESTIMATION);

            min_counts[j]
                    = boost::math::quantile(boost::math::poisson(mean_est),
                                            exp(-mean_est) - CDF_CUTOFF * expm1(-mean_est));
        }
    }

    auto correct_pvals = [&config](const std::vector<std::pair<long double, size_t>>& m) {
        size_t total
                = std::accumulate(m.begin(), m.end(), size_t(0),
                                  [](size_t sum, const auto& a) { return sum + a.second; });
        size_t acc = 0;

        for (size_t n = 0; n < m.size(); ++n) {
            const auto& [pval_min, s] = m[n];
            if (s == 0)
                continue;

            size_t c = floor(config.family_wise_error_rate / pval_min + 1);
            size_t nsig = total - acc;
            // there are nsig k-mers s.t. pval_min < alpha / c
            if (nsig > 0 && c >= nsig)
                return std::make_tuple(c, nsig, n);

            acc += s;
        }

        return std::make_tuple(std::numeric_limits<size_t>::max(),
                               std::numeric_limits<size_t>::max(),
                               std::numeric_limits<size_t>::max());
    };

    common::logger->trace("Marking discarded k-mers");
    auto hists_map = get_hist_map(min_counts, &kept_bv);
    std::unique_ptr<MaskedDeBruijnGraph> clean_masked_graph;
    // bit_vector_stat kept(std::move(kept_bv));
    bit_vector_stat kept(kept_bv);

    clean_masked_graph = std::make_unique<MaskedDeBruijnGraph>(
            graph_ptr,
            [&](node_index node) {
                return node != DeBruijnGraph::npos
                        && kept[AnnotatedDBG::graph_to_anno_index(node)];
            },
            true, is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC);

    // common::logger->trace("Filtering low-complexity k-mers");
    // std::atomic_thread_fence(std::memory_order_release);
    // clean_masked_graph->call_sequences([&](const std::string &contig, const auto &path) {
    //     for (size_t i = 0; i < path.size(); ++i) {
    //         if (is_low_complexity(std::string_view(contig.c_str() + i, graph_ptr->get_k())))
    //             unset_bit(kept_bv.data(), AnnotatedDBG::graph_to_anno_index(path[i]), true, std::memory_order_relaxed);
    //     }
    // }, num_threads);
    // std::atomic_thread_fence(std::memory_order_acquire);
    // kept = bit_vector_stat(std::move(kept_bv));


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

    // std::vector<PValStorage> pvals_buckets;
    // std::vector<std::unique_ptr<utils::TempFile>> tmp_buckets;
    // std::vector<PValStorage> pvals_min_buckets;
    // std::vector<std::unique_ptr<utils::TempFile>> tmp_min_buckets;
    // std::vector<PValStorage> eff_size_buckets;
    // std::vector<std::unique_ptr<utils::TempFile>> tmp_eff_buckets;
    // // constexpr bool preallocated = std::is_same_v<PValStorage, std::vector<uint64_t>>;

    // if (config.test_type != "notest" && config.test_by_unitig) {
    //     if constexpr(std::is_same_v<PValStorage, std::vector<uint64_t>>) {
    //         auto &pvals = pvals_buckets.emplace_back();
    //         pvals.reserve(graph_ptr->max_index());

    //         auto &pvals_min = pvals_min_buckets.emplace_back();
    //         pvals_min.reserve(graph_ptr->max_index());

    //         auto &eff_size = eff_size_buckets.emplace_back();
    //         eff_size.reserve(graph_ptr->max_index());
    //     }

    //     if constexpr(std::is_same_v<PValStorage, sdsl::int_vector_buffer<64>>) {
    //         for (size_t i = 0; i < get_num_threads() + 1; ++i) {
    //             auto &tmp_file = tmp_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
    //             pvals_buckets.emplace_back(tmp_file->name(), std::ios::out);
    //         }

    //         for (size_t i = 0; i < get_num_threads() + 1; ++i) {
    //             auto &tmp_file = tmp_min_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
    //             pvals_min_buckets.emplace_back(tmp_file->name(), std::ios::out);
    //         }

    //         for (size_t i = 0; i < get_num_threads() + 1; ++i) {
    //             auto &tmp_file_eff = tmp_eff_buckets.emplace_back(std::make_unique<utils::TempFile>(tmp_dir));
    //             eff_size_buckets.emplace_back(tmp_file_eff->name(), std::ios::out);
    //         }
    //     }
    // }

    // if ((config.test_type == "gnb_exact" || config.test_type == "lnb_exact") && config.test_by_unitig) {

    // }


    // std::vector<PValStorage> sum_buckets;
    // std::vector<std::unique_ptr<utils::TempFile>> tmp_sum_buckets;
    // if (config.test_by_unitig) {
    //     common::logger->trace("Allocating effect size and row sum storage");

    // }

    sdsl::bit_vector indicator_in;
    sdsl::bit_vector indicator_out;

    common::logger->trace("Allocating k-mer bitmasks");
    indicator_in = sdsl::bit_vector(graph_ptr->max_index() + 1, false);
    indicator_out = sdsl::bit_vector(graph_ptr->max_index() + 1, false);

    common::logger->trace("Test: {}\tby unitig: {}", config.test_type, config.test_by_unitig);

    // prefilter rows
    auto generate_clean_rows = [&](const auto& callback) {
        generate_rows([&](uint64_t row_i, const auto& row, size_t bucket_idx) {
            if (!kept[row_i])
                return;

            size_t count_in = 0;
            size_t count_out = 0;
            size_t total_count = 0;
            for (const auto& [j, c] : row) {
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

            bool in_kmer = count_in >= config.min_in_recurrence
                    && count_in <= config.max_in_recurrence;
            bool out_kmer = count_out >= config.min_out_recurrence
                    && count_out <= config.max_out_recurrence;

            if (in_kmer || out_kmer) {
                callback(row_i, row, bucket_idx);
                // std::vector<std::pair<size_t, size_t>> new_row;
                // for (const auto &[j, c] : row) {
                //     if (c < 1024)
                //         new_row.emplace_back(j, c);
                // }
                // callback(row_i, new_row, bucket_idx);
            }
        });
    };

    // common::logger->trace("Marking low-complexity k-mers");
    // std::atomic_thread_fence(std::memory_order_release);
    // clean_masked_graph->call_sequences([&](const std::string &contig, const auto &path) {
    //     for (size_t i = 0; i < path.size(); ++i) {
    //         std::string_view kmer(contig.c_str() + i, clean_masked_graph->get_k());
    //         if (is_low_complexity(kmer)) {
    //             auto row = AnnotatedDBG::graph_to_anno_index(path[i]);
    //             unset_bit(kept_bv.data(), row, true, std::memory_order_relaxed);
    //         }
    //     }
    // }, num_threads);
    // std::atomic_thread_fence(std::memory_order_acquire);

    // common::logger->trace("Discarding counts from low-complexity k-mers");
    // {
    //     std::vector<std::vector<VectorMap<uint64_t, size_t>>> hist_maps_m(num_threads +
    //     1); for (auto &hist_maps : hist_maps_m) {
    //         hist_maps.resize(groups.size());
    //     }
    //     std::atomic_thread_fence(std::memory_order_release);
    //     generate_clean_rows([&](auto row_i, const auto &row, size_t bucket_idx) {
    //         if (fetch_bit(kept_bv.data(), row_i, true, std::memory_order_acquire)) {
    //             if (row.empty()) {
    //                 unset_bit(kept_bv.data(), row_i, true, std::memory_order_release);
    //             } else {
    //                 sdsl::bit_vector found(groups.size());
    //                 for (const auto &[j, c] : row) {
    //                     found[j] = true;
    //                     ++hist_maps_m[bucket_idx][j][c];
    //                 }
    //                 for (size_t j = 0; j < found.size(); ++j) {
    //                     if (!found[j])
    //                         ++hist_maps_m[bucket_idx][j][0];
    //                 }
    //             }
    //         }
    //     });
    //     std::atomic_thread_fence(std::memory_order_acquire);
    //     std::swap(hists_map, hist_maps_m[0]);
    //     for (size_t i = 1; i < hist_maps_m.size(); ++i) {
    //         for (size_t j = 0; j < hist_maps_m[i].size(); ++j) {
    //             for (const auto &[k, c] : hist_maps_m[i][j]) {
    //                 hists_map[j][k] += c;
    //             }
    //         }
    //     }
    // }

    // std::atomic_thread_fence(std::memory_order_release);
    // generate_rows([&](uint64_t row_i, const auto &row, size_t) {
    //     if (kept[row_i] && !kept_bv[row_i]) {
    //         sdsl::bit_vector found(groups.size());
    //         for (const auto &[j, c] : row) {
    //             found[j] = true;
    //             auto find = hists_map[j].find(c);
    //             assert(find != hists_map[j].end());
    //             uint64_t prev = __atomic_fetch_sub(&find.value(), 1,
    //             std::memory_order_relaxed); std::ignore = prev; assert(prev);
    //         }
    //         for (size_t j = 0; j < found.size(); ++j) {
    //             if (!found[j]) {
    //                 auto find = hists_map[j].find(0);
    //                 assert(find != hists_map[j].end());
    //                 uint64_t prev = __atomic_fetch_sub(&find.value(), 1,
    //                 std::memory_order_relaxed); std::ignore = prev; assert(prev);
    //             }
    //         }
    //     }
    // });
    // std::atomic_thread_fence(std::memory_order_acquire);
    // size_t old_nelem = kept.num_set_bits();
    // kept = bit_vector_stat(std::move(kept_bv));

    size_t nelem = kept.num_set_bits();
    // common::logger->trace("Kept {} / {} high-complexity k-mers", nelem, old_nelem);

    for (size_t j = 0; j < groups.size(); ++j) {
        size_t total_c = 0;
        for (const auto& [k, c] : hists_map[j]) {
            total_c += c;
        }

        assert(total_c == nelem);

        if (total_c != nelem) {
            common::logger->error("{}: {} != {}", j, total_c, nelem);
            throw std::runtime_error("FAIL");
        }
    }

    common::logger->trace("Computing aggregate statistics");
    int64_t in_kmers = 0;
    int64_t out_kmers = 0;
    int64_t in_sq_kmers = 0;
    int64_t out_sq_kmers = 0;
    int64_t total_sq_kmers = 0;
    std::vector<uint64_t> max_obs_vals(groups.size());
    std::vector<size_t> n_kmers(groups.size());
    std::vector<uint64_t> sums(groups.size());
    std::vector<uint64_t> sq_sums(groups.size());
    uint64_t max_in_obs_val = 0.0;
    uint64_t max_out_obs_val = 0.0;

    int64_t total_kmers = 0;
    for (size_t j = 0; j < groups.size(); ++j) {
        for (const auto& [k, c] : hists_map[j]) {
            sums[j] += k * c;
            sq_sums[j] += k * k * c;
            if (k != 0)
                n_kmers[j] += c;

            max_obs_vals[j] = std::max(max_obs_vals[j], k);
        }

        if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
            out_kmers += sums[j];
            out_sq_kmers += sq_sums[j];
            max_out_obs_val += max_obs_vals[j];
        }

        if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
            in_kmers += sums[j];
            in_sq_kmers += sq_sums[j];
            max_in_obs_val += max_obs_vals[j];
        }

        total_kmers += sums[j];
        total_sq_kmers += sq_sums[j];

        common::logger->trace(
                "{}: n_unique: {}\tsum: {}\tmax_obs: {}\tmin_cutoff: {}\tmax_cutof: {}",
                j, n_kmers[j], sums[j], max_obs_vals[j], min_counts[j], check_cutoff[j]);
    }

    common::logger->trace("Number of kept unique k-mers: {}\tNumber of kept k-mers: {}",
                          nelem, total_kmers);

    std::vector<std::vector<std::pair<uint64_t, size_t>>> hists(groups.size());
    for (size_t j = 0; j < hists.size(); ++j) {
        hists[j] = const_cast<std::vector<std::pair<uint64_t, size_t>>&&>(
                hists_map[j].values_container());
        std::sort(hists[j].begin(), hists[j].end(), utils::LessFirst());
        hists_map[j].clear();
    }
    hists_map.resize(0);

    std::vector<size_t> kmer_to_unitig;
    std::vector<std::tuple<size_t, size_t, size_t>> counts;
    // std::vector<std::pair<size_t, std::vector<uint64_t>>> agg_counts;
    if (config.test_by_unitig && config.test_type != "notest") {
        clean_masked_graph = std::make_unique<MaskedDeBruijnGraph>(
                graph_ptr,
                [&](node_index node) {
                    return node != DeBruijnGraph::npos
                            && kept[AnnotatedDBG::graph_to_anno_index(node)];
                },
                true, is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC);

        common::logger->trace("Associating k-mers to unitigs");
        kmer_to_unitig.resize(graph_ptr->max_index() + 1);
        // std::atomic_thread_fence(std::memory_order_release);
        // std::atomic<size_t> num_unitigs{0};
        std::mutex mu;
        clean_masked_graph->call_unitigs(
                [&](const std::string&, const auto& path) {
                    size_t unitig_id = 0;
                    {
                        unitig_id = counts.size();
                        std::lock_guard<std::mutex> lock(mu);
                        counts.emplace_back(path.size(), 0, 0);
                    }
                    // size_t unitig_id = num_unitigs.fetch_add(1, std::memory_order_relaxed);
                    for (node_index node : path) {
                        kmer_to_unitig[node] = unitig_id;
                    }
                },
                num_threads);
        // std::atomic_thread_fence(std::memory_order_acquire);
        common::logger->trace("Assembled {} unitigs", counts.size());

        common::logger->trace("Aggregating counts for unitigs");
        generate_clean_rows([&](uint64_t row_i, const auto& row, size_t) {
            node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
            size_t unitig_id = kmer_to_unitig[node];
            auto& [len, in_count, out_count] = counts[unitig_id];
            for (const auto& [j, c] : row) {
                if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                    out_count += c;

                if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                    in_count += c;
            }
        });

        common::logger->trace("Computing unitig count histogram");
        std::vector<VectorMap<size_t, size_t>> unitig_hists_map(2);
        for (const auto& [len, in_count, out_count] : counts) {
            size_t avg_in_count = std::ceil(static_cast<double>(in_count) / len);
            size_t avg_out_count = std::ceil(static_cast<double>(out_count) / len);
            ++unitig_hists_map[0][avg_in_count];
            ++unitig_hists_map[1][avg_out_count];
        }

        common::logger->trace("Fitting negative binomials");
        auto [r_a, p_a, mu_a, var_a] = get_rp([&](const auto& callback) {
            for (const auto& [k, c] : unitig_hists_map[0]) {
                callback(k, c);
            }
        });
        auto [r_b, p_b, mu_b, var_b] = get_rp([&](const auto& callback) {
            for (const auto& [k, c] : unitig_hists_map[1]) {
                callback(k, c);
            }
        });
        auto [r_n, p_n, mu_n, var_n] = get_rp([&](const auto& callback) {
            for (const auto& [k, c] : unitig_hists_map[0]) {
                auto find = unitig_hists_map[1].find(k);
                if (find != unitig_hists_map[1].end()) {
                    callback(k, c + find->second);
                } else {
                    callback(k, c);
                }
            }

            for (const auto& [k, c] : unitig_hists_map[1]) {
                if (!unitig_hists_map[0].count(k))
                    callback(k, c);
            }
        });

        common::logger->trace("In: r: {}\tp: {}\t mu: {}\tsample var: {}\tvar: {}", r_a,
                              p_a, mu_a, var_a, mu_a / p_a);
        common::logger->trace("Out: r: {}\tp: {}", r_b, p_b, mu_b, var_b, mu_b / p_b);
        common::logger->trace("Null: r: {}\tp: {}", r_n, p_n, mu_n, var_n, mu_n / p_n);

        long double nb_base = (lgammal(r_n) - lgammal(r_a) - lgammal(r_b)) / logl(2.0)
                + r_a * log2l(p_a) + r_b * log2l(p_b) - r_n * log2l(p_n);

        common::logger->trace("Allocating p-value storage");
        std::vector<std::tuple<long double, Group, size_t>> pvals(
                counts.size(), std::make_tuple(1.0, Group::OTHER, counts.size()));

        common::logger->trace("Calculating p-values");
        size_t rows_per_update = 10000;
        ProgressBar progress_bar(pvals.size(), "Streaming unitigs", std::cerr,
                                 !common::get_verbose());
#pragma omp parallel for num_threads(num_threads) schedule(static)
        for (size_t i = 0; i < counts.size(); ++i) {
            if (i > 0 && (i % rows_per_update) == 0)
                progress_bar += rows_per_update;

            const auto& [len, in_sum, out_sum] = counts[i];
            size_t n = in_sum + out_sum;

            long double eff_size = 0.0;
            auto& [pval, group, idx] = pvals[i];
            idx = i;

            long double midpoint = r_a * n / (r_a + r_b);

            auto get_deviance_a = [&](long double s) {
                if (s == 0)
                    return r_a * logl(p_a) * -2.0;

                return ((logl(s) - log1pl(-p_a)) * s - (r_a + s) * logl(r_a + s)
                        + r_a * (logl(r_a) - logl(p_a)))
                        * 2.0;
            };

            auto get_deviance_b = [&](long double t) {
                if (t == 0)
                    return r_b * logl(p_b) * -2.0;

                return ((logl(t) - log1pl(-p_b)) * t - (r_b + t) * logl(r_b + t)
                        + r_b * (logl(r_b) - logl(p_b)))
                        * 2.0;
            };

            auto get_deviance = [&](long double s, long double t) {
                return get_deviance_a(s) + get_deviance_b(t);
            };

            long double min_dev
                    = get_deviance(midpoint, static_cast<long double>(n) - midpoint);
            long double in_dev = get_deviance(in_sum, out_sum);

            if (in_dev > min_dev) {
                eff_size = static_cast<long double>(in_sum) - midpoint;
                if (eff_size > 0) {
                    group = Group::IN;
                } else if (eff_size < 0) {
                    group = Group::OUT;
                }
                pval = 0.0;
                idx = i;
                long double base = nb_base
                        + (lgammal(n + 1) - lgammal(r_n + n) - n * log1pl(-p_n)) / logl(2.0);
                double l1pa = log1pl(-p_a) / logl(2.0);
                double l1pb = log1pl(-p_b) / logl(2.0);

                {
                    size_t s = 0;
                    size_t t = n;
                    long double sbase = (lgammal(r_a + s) + lgammal(r_b + t)
                                         - lgammal(s + 1) - lgammal(t + 1))
                                    / logl(2.0)
                            + s * l1pa + t * l1pb;
                    pval += exp2(base + sbase);
                    for (++s; s <= n; ++s) {
                        long double dev = get_deviance(s, t - 1);
                        if (dev >= in_dev) {
                            sbase += log2l(r_a + s) - log2l(r_b + t) - log2l(s + 1)
                                    + log2l(t + 1) + l1pa - l1pb;
                            --t;
                            pval += exp2(base + sbase);
                            if (pval >= config.family_wise_error_rate)
                                break;
                        } else {
                            break;
                        }
                    }
                }
                if (pval < config.family_wise_error_rate && get_deviance(n, 0) >= in_dev) {
                    size_t s = n;
                    size_t t = 0;
                    long double sbase = (lgammal(r_a + s) + lgammal(r_b + t)
                                         - lgammal(s + 1) - lgammal(t + 1))
                                    / logl(2.0)
                            + s * l1pa + t * l1pb;
                    pval += exp2(base + sbase);
                    for (++t; t <= n; ++t) {
                        long double dev = get_deviance(s - 1, t);
                        if (dev >= in_dev) {
                            sbase -= log2l(r_a + s) - log2l(r_b + t) - log2l(s + 1)
                                    + log2l(t + 1) + l1pa - l1pb;
                            --s;
                            pval += exp2(base + sbase);
                            if (pval >= config.family_wise_error_rate)
                                break;
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        progress_bar += (pvals.size() % rows_per_update);

        common::logger->trace("Correcting p-values");
        common::logger->trace("Sorting {} / {} p-vals", pvals.size(), counts.size());
        std::sort(pvals.begin(), pvals.end());
        long double harm = 1.0;
        // long double harm = 0.0;
        // for (size_t i = 1; i <= nelem; ++i) {
        //     harm += 1.0 / i;
        // }

        common::logger->trace("Selecing cut-off");
        size_t k = 0;
        for (size_t i = 0; i < pvals.size(); ++i) {
            const auto& [pval, group, idx] = pvals[i];
            if (pval >= config.family_wise_error_rate)
                break;

            if (pval <= config.family_wise_error_rate * (i + 1) / pvals.size() / harm)
                k = std::max(i + 1, k);
        }

        common::logger->trace("Found {} / {} significant p-values. Minimum p-value is {}",
                              k, pvals.size(), std::get<0>(pvals[0]));

        for (size_t i = k; i < pvals.size(); ++i) {
            std::get<1>(pvals[i]) = Group::BOTH;
        }

        common::logger->trace("Reordering p-values");
        std::sort(pvals.begin(), pvals.end(), [&](const auto& a, const auto& b) {
            return std::get<2>(a) < std::get<2>(b);
        });

        common::logger->trace("Marking significant p-values");
        std::atomic_thread_fence(std::memory_order_release);
        generate_clean_rows([&](uint64_t row_i, const auto& row, size_t) {
            node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
            size_t unitig_id = kmer_to_unitig[node];
            const auto& [pval, group, idx] = pvals[unitig_id];
            if (group == Group::IN)
                set_bit(indicator_in.data(), node, std::memory_order_relaxed);

            if (group == Group::OUT)
                set_bit(indicator_out.data(), node, std::memory_order_relaxed);
        });
        std::atomic_thread_fence(std::memory_order_acquire);
    } else {
        // precompute negative binomial fits for poisson_bayes and nbinom_exact tests
        std::vector<std::pair<long double, long double>> nb_params(groups.size());
        std::pair<long double, long double> nb_params_a;
        std::pair<long double, long double> nb_params_b;
        std::pair<long double, long double> nb_params_null;
        long double nb_base = 0.0;
        // long double alpha = in_kmers;
        // long double beta = out_kmers;
        // long double hg_base = (lgammal(alpha + beta) - lgammal(alpha) - lgammal(beta)) / logl(2.0);
        if (config.test_type == "nbinom_exact") {
            common::logger->trace("Fitting per-sample negative binomial distributions");


#pragma omp parallel for num_threads(num_parallel_files)
            for (size_t j = 0; j < groups.size(); ++j) {
                const auto& hist = hists[j];
                if (hist.size()) {
                    auto [r, p, mu, var] = get_rp([&](const auto& callback) {
                        for (const auto& [k, c] : hist) {
                            callback(k, c);
                        }
                    });
                    nb_params[j] = std::make_pair(r, p);
                    common::logger->trace(
                            "{}: size: {}\tmax_val: {}\tsample mean: {}\tsample var: "
                            "{}\tmu: {}\tvar: {}\tmle: r: {}\tp: {}",
                            j, sums[j], (hist.end() - 1)->first, mu, var, r * (1 - p) / p,
                            r * (1 - p) / p / p, r, p);
                }
            }

            if (config.test_type == "nbinom_exact") {
                // based on this method for combining gamma distributions
                // https://stats.stackexchange.com/questions/72479/generic-sum-of-gamma-random-variables
                long double num_a = 0.0;
                long double num_b = 0.0;

                long double denom_a = 0.0;
                long double denom_b = 0.0;

                for (size_t j = 0; j < groups.size(); ++j) {
                    const auto& [r, p] = nb_params[j];
                    long double theta = (1.0 - p) / p;
                    if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
                        num_b += theta * r;
                        denom_b += theta * theta * r;
                    }

                    if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
                        num_a += theta * r;
                        denom_a += theta * theta * r;
                    }
                }

                long double r_a = num_a * num_a / denom_a;
                long double theta_a = num_a / r_a;

                long double r_b = num_b * num_b / denom_b;
                long double theta_b = num_b / r_b;

                long double r_null = pow(num_a + num_a, 2.0) / (denom_a + denom_b);
                long double theta_null = (num_a + num_b) / r_null;

                nb_params_a.first = r_a;
                nb_params_b.first = r_b;
                nb_params_null.first = r_null;

                // theta = (1 - p) / p
                // ptheta = 1 - p
                // ptheta + p = 1
                // p = 1/(1+theta)

                nb_params_a.second = 1.0 / (theta_a + 1.0);
                nb_params_b.second = 1.0 / (theta_b + 1.0);
                nb_params_null.second = 1.0 / (theta_null + 1.0);
                // long double mu_a = 0.0;
                // long double mu_b = 0.0;
                // long double var_a = 0.0;
                // long double var_b = 0.0;
                // for (size_t j = 0; j < groups.size(); ++j) {
                //     const auto &[r, p] = nb_params[j];
                //     long double mu = r * (1.0 - p) / p;
                //     long double var = mu / p;
                //     if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
                //         mu_a += mu;
                //         var_a += var;
                //     }

                //     if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
                //         mu_b += mu;
                //         var_b += var;
                //     }
                // }

                // nb_params_a.first = mu_a * mu_a / (var_a - mu_a);
                // nb_params_a.second = mu_a / var_a;
                common::logger->trace("In: r: {}\tp: {}", nb_params_a.first,
                                      nb_params_a.second);

                // nb_params_b.first = mu_b * mu_b / (var_b - mu_b);
                // nb_params_b.second = mu_b / var_b;
                common::logger->trace("Out: r: {}\tp: {}", nb_params_b.first,
                                      nb_params_b.second);

                // long double mu_null = mu_a + mu_b;
                // long double var_null = var_a + var_b; // the two are independent, so
                // their covariance is 0 nb_params_null.first = mu_null * mu_null /
                // (var_null - mu_null); nb_params_null.second = mu_null / var_null;
                common::logger->trace("Null: r: {}\tp: {}", nb_params_null.first,
                                      nb_params_null.second);

                nb_base = (lgammal(nb_params_null.first) - lgammal(nb_params_a.first)
                           - lgammal(nb_params_b.first))
                        / logl(2.0);
                nb_base += nb_params_a.first * log2l(nb_params_a.second)
                        + nb_params_b.first * log2l(nb_params_b.second)
                        - nb_params_null.first * log2l(nb_params_null.second);
            }
        }

        using RowGenerator
                = std::function<void(const std::function<void(size_t, uint64_t)>&)>;
        using GetLambda
                = std::function<std::pair<long double, long double>(size_t, const RowGenerator&)>;
        GetLambda get_lambda;

        std::pair<long double, long double> pln_params_a;
        std::pair<long double, long double> pln_params_b;
        std::pair<long double, long double> pln_params_null;
        if (config.test_type == "poisson_bayes") {
            auto get_params = [&](long double kmers, long double sq_kmers) {
                long double mu = static_cast<long double>(kmers) / nelem;
                long double mu2 = mu * mu;
                long double var = sq_kmers / nelem - mu2;
                long double log_sig2 = log((var - mu) / mu2 + 1);
                long double log_mu = log(mu) - log_sig2 / 2.0;
                return std::make_pair(log_mu, log_sig2);
            };

            pln_params_a = get_params(in_kmers, in_sq_kmers);
            pln_params_b = get_params(out_kmers, out_sq_kmers);
            pln_params_null = get_params(total_kmers, total_sq_kmers);
            common::logger->trace("Lognormal MoM\tin:\tmu: {}\tvar: {}",
                                  pln_params_a.first, pln_params_a.second);
            common::logger->trace("Lognormal MoM\tout:\tmu: {}\tvar: {}",
                                  pln_params_b.first, pln_params_b.second);
            common::logger->trace("Lognormal MoM\tnull:\tmu: {}\tvar: {}",
                                  pln_params_null.first, pln_params_null.second);

            // common::logger->trace("Precomputing lambdas");

            // #pragma omp parallel for num_threads(num_threads)
            // for (size_t j = 0; j < groups.size(); ++j) {
            //     auto [mu, sig2] = groups[j] == Group::IN ? pln_params_a : pln_params_b;
            //     for (const auto &[k, c] : hists[j]) {
            //         auto [lambda_j_min, lambda_j_max] = boost::math::tools::bisect(
            //             [&](long double lambda) {
            //                 return (mu - logl(lambda)) / sig2 - lambda + k - 1;
            //             },
            //             std::numeric_limits<long double>::min(),
            //             static_cast<long double>(max_obs_vals[j]),
            //             boost::math::tools::eps_tolerance<long double>(5)
            //         );
            //         lambdas[j][k] = (lambda_j_min + lambda_j_max) / 2;
            //     }
            // }
        }

        if (config.test_type == "poisson_exact") {
            get_lambda = [&](size_t, const RowGenerator&) {
                return std::make_pair(static_cast<long double>(in_kmers) / nelem,
                                      static_cast<long double>(out_kmers) / nelem);
            };
        }

        std::vector<std::pair<tsl::hopscotch_map<size_t, long double>,
                              tsl::hopscotch_map<size_t, long double>>>
                lambdas(num_threads + 1);
        if (config.test_type == "poisson_bayes") {
            get_lambda = [&](size_t bucket_idx, const RowGenerator& generate)
                    -> std::pair<long double, long double> {
                size_t in_sum = 0;
                size_t out_sum = 0;
                generate([&](size_t j, uint64_t c) {
                    Group gp = groups[j];

                    if (gp == Group::OUT || gp == Group::BOTH)
                        out_sum += c;

                    if (gp == Group::IN || gp == Group::BOTH)
                        in_sum += c;
                });

                auto& [in_lambdas, out_lambdas] = lambdas[bucket_idx];
                auto in_found = in_lambdas.find(in_sum);
                auto out_found = out_lambdas.find(out_sum);

                if (in_found == in_lambdas.end()) {
                    auto [mu, sig2] = pln_params_a;
                    auto [lambda_j_min, lambda_j_max] = boost::math::tools::bisect(
                            [&](long double lambda) {
                                return (mu - logl(lambda)) / sig2 - lambda + in_sum - 1;
                            },
                            std::numeric_limits<long double>::min(),
                            static_cast<long double>(max_in_obs_val),
                            boost::math::tools::eps_tolerance<long double>(5));
                    long double lambda = (lambda_j_min + lambda_j_max) / 2;
                    in_found = in_lambdas.try_emplace(in_sum, lambda).first;
                }

                if (out_found == out_lambdas.end()) {
                    auto [mu, sig2] = pln_params_b;
                    auto [lambda_j_min, lambda_j_max] = boost::math::tools::bisect(
                            [&](long double lambda) {
                                return (mu - logl(lambda)) / sig2 - lambda + out_sum - 1;
                            },
                            std::numeric_limits<long double>::min(),
                            static_cast<long double>(max_out_obs_val),
                            boost::math::tools::eps_tolerance<long double>(5));
                    long double lambda = (lambda_j_min + lambda_j_max) / 2;
                    out_found = out_lambdas.try_emplace(out_sum, lambda).first;
                }

                return std::make_pair(in_found->second, out_found->second);
            };
        }

        // precompute p-values for poisson_binom test
        std::vector<std::vector<long double>> pb_pvals(num_labels_in + num_labels_out + 1);
        std::vector<long double> mid_points(num_labels_in + num_labels_out + 1);
        if (config.test_type == "poisson_binom") {
            // fit distribution
            std::vector<long double> p;
            p.reserve(groups.size());
            for (size_t i = 0; i < groups.size(); ++i) {
                p.emplace_back(static_cast<long double>(n_kmers[i]) / nelem);
            }

            common::logger->trace("p: {}", fmt::join(p, ","));

            common::logger->trace("Precomputing PMFs");
            std::vector<long double> pmf_in { 1.0L };
            std::vector<long double> pmf_out { 1.0L };
            std::vector<long double> pmf_null { 1.0L };

            for (size_t i = 1; i <= groups.size(); ++i) {
                if (groups[i - 1] == Group::OUT || groups[i - 1] == Group::BOTH) {
                    std::vector<long double> pmf_out_cur(pmf_out.size() + 1);
                    // pmf_out_cur[0] = (1.0L - p[i - 1]) * pmf_out[0];
                    pmf_out_cur[0] = expl(log1pl(-p[i - 1]) + logl(pmf_out[0]));
                    pmf_out_cur[pmf_out.size()] = p[i - 1] * pmf_out.back();
                    for (size_t k = 1; k < pmf_out.size(); ++k) {
                        pmf_out_cur[k] = p[i - 1] * pmf_out[k - 1]
                                + (1.0L - p[i - 1]) * pmf_out[k];
                    }
                    std::swap(pmf_out_cur, pmf_out);
                }

                if (groups[i - 1] == Group::IN || groups[i - 1] == Group::BOTH) {
                    std::vector<long double> pmf_in_cur(pmf_in.size() + 1);
                    // pmf_in_cur[0] = (1.0L - p[i - 1]) * pmf_in[0];
                    pmf_in_cur[0] = expl(log1pl(-p[i - 1]) + logl(pmf_in[0]));
                    pmf_in_cur[pmf_in.size()] = p[i - 1] * pmf_in.back();
                    for (size_t k = 1; k < pmf_in.size(); ++k) {
                        pmf_in_cur[k]
                                = p[i - 1] * pmf_in[k - 1] + (1.0L - p[i - 1]) * pmf_in[k];
                    }
                    std::swap(pmf_in_cur, pmf_in);
                }
                std::vector<long double> pmf_null_cur(i + 1);
                // pmf_null_cur[0] = (1.0L - p[i - 1]) * pmf_null[0];
                pmf_null_cur[0] = expl(log1pl(-p[i - 1]) + logl(pmf_null[0]));
                pmf_null_cur[pmf_null.size()] = p[i - 1] * pmf_null.back();
                for (size_t k = 1; k < pmf_null.size(); ++k) {
                    pmf_null_cur[k]
                            = p[i - 1] * pmf_null[k - 1] + (1.0L - p[i - 1]) * pmf_null[k];
                }
                std::swap(pmf_null_cur, pmf_null);
            }

            if (pmf_in.size() != num_labels_in + 1) {
                common::logger->error("PMF in wrong: {} != {}", pmf_in.size(),
                                      num_labels_in + 1);
                throw std::domain_error("");
            }

            if (pmf_out.size() != num_labels_out + 1) {
                common::logger->error("PMF out wrong: {} != {}", pmf_out.size(),
                                      num_labels_out + 1);
                throw std::domain_error("");
            }

            if (pmf_null.size() != num_labels_in + num_labels_out + 1) {
                common::logger->error("PMF null wrong: {} != {}", pmf_null.size(),
                                      num_labels_in + num_labels_out + 1);
                throw std::domain_error("");
            }

            common::logger->trace("Precomputing p-values");
            pb_pvals[0].emplace_back(1.0);

            long double min_pval = 1.0;
            // long double last_local_min_pval = 1.0;

            for (size_t n = 1; n < pb_pvals.size(); ++n) {
                std::vector<long double> probs;
                size_t front = n - std::min(n, num_labels_out);
                for (uint64_t s = 0; s < pmf_in.size(); ++s) {
                    if (s > n)
                        break;
                    uint64_t t = n - s;
                    if (s < pmf_in.size() && t < pmf_out.size()) {
                        if (s < front) {
                            common::logger->error(
                                    "Attempting non-zero p-value in impossible "
                                    "configuration: {},{}",
                                    n, s);
                            throw std::domain_error("");
                        }
                        probs.emplace_back(exp2l(log2l(pmf_in[s]) + log2l(pmf_out[t])
                                                 - log2l(pmf_null[n])));
                    } else {
                        probs.emplace_back(0.0);
                    }
                }

                long double sum_probs = std::accumulate(probs.begin(), probs.end(), 0.0L);
                if (abs(sum_probs - 1.0L) > 1e-5) {
                    common::logger->error("Sum of probs for n={} = {}", n, sum_probs);
                    throw std::runtime_error("Fail");
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
                        common::logger->error(
                                "Non-zero p-value in impossible configuration: {},{}", n, s);
                        throw std::domain_error("");
                    }
                }
                long double local_min_pval
                        = std::min(pb_pvals[n][front], pb_pvals[n].back());
                min_pval = std::min(min_pval, local_min_pval);

                size_t max_prob_counts = 0;
                long double max_prob_pos = 0;
                // double max_prob = -1;

                for (uint64_t s = front; s < pb_pvals[n].size(); ++s) {
                    if (pb_pvals[n][s] < local_min_pval) {
                        common::logger->error(
                                "Min p-value not at boundary: {},{}: {} < {}", n, s,
                                pb_pvals[n][s], local_min_pval);
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
                        if (probs[s] == probs[s - 1]) {
                            max_prob_counts = 2;
                            max_prob_pos = s + s;
                            break;
                        }

                        if (s + 1 < pb_pvals[n].size()) {
                            if (probs[s - 1] < probs[s] && probs[s] < probs[s + 1])
                                continue;

                            if (probs[s - 1] > probs[s] && probs[s] > probs[s + 1])
                                continue;

                            if (probs[s] > probs[s - 1] && probs[s] > probs[s + 1]) {
                                if (probs[s] - probs[s - 1] == probs[s] - probs[s + 1]) {
                                    max_prob_counts = 1;
                                    max_prob_pos = s;
                                } else if (probs[s] - probs[s - 1] > probs[s] - probs[s + 1]) {
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

                mid_points[n] = static_cast<long double>(max_prob_pos) / max_prob_counts;
                common::logger->trace("Midpoint: n: {}\tmp: {}\t{}\t{}", n, mid_points[n],
                                      fmt::join(probs, ","), fmt::join(pb_pvals[n], ","));

                // if (local_min_pval > last_local_min_pval) {
                //     common::logger->error("Min p-vals failed: n-1: {} -> {}\tn: {} -> {}",
                //                           n - 1, last_local_min_pval, n, local_min_pval);
                //     throw std::runtime_error("Min p-val failed");
                // }
                // last_local_min_pval = local_min_pval;
            }

            common::logger->trace("Min. p-value: {}", min_pval);
            if (min_pval >= config.family_wise_error_rate) {
                common::logger->warn("No significant p-values achievable");
                auto masked_graph_in = std::make_shared<MaskedDeBruijnGraph>(
                        graph_ptr, [](node_index) { return false; }, true,
                        is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC);

                auto masked_graph_out = std::make_shared<MaskedDeBruijnGraph>(
                        graph_ptr, [](node_index) { return false; }, true,
                        is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC);

                return std::make_tuple(masked_graph_in, masked_graph_out, PValStorage {},
                                       nullptr);
            }
        }

        // std::pair<long double, long double> bb_params_a;
        // std::pair<long double, long double> bb_params_b;
        // std::pair<long double, long double> bb_params_null;
        std::vector<std::tuple<long double, long double, long double>> bb_bases;
        // long double bb_alpha = 0.0;
        // long double bb_beta = 0.0;
        if (config.test_type == "bb_hypergeometric") {
            std::vector<std::vector<std::tuple<size_t, long double, long double>>> probs_m(
                    num_threads + 1);
            generate_clean_rows([&](auto row_i, const auto& row, size_t bucket_idx) {
                size_t n = 0;
                size_t in_sum = 0;
                for (const auto& [j, c] : row) {
                    if (groups[j] == Group::IN)
                        in_sum += c;

                    n += c;
                }
                long double p = static_cast<long double>(in_sum) / n;
                if (probs_m[bucket_idx].size() <= n)
                    probs_m[bucket_idx].resize(n + 1);

                auto& [c, prob_sum, sq_prob_sum] = probs_m[bucket_idx][n];
                prob_sum += p;
                sq_prob_sum += p * p;
                ++c;
            });

            size_t max_size = 0;
            for (const auto& probs : probs_m) {
                max_size = std::max(max_size, probs.size());
            }
            std::vector<std::tuple<size_t, long double, long double>> probs(max_size);
            for (const auto& prob : probs_m) {
                for (size_t n = 0; n < prob.size(); ++n) {
                    auto& [total_c, total_prob_sum, total_sq_prob_sum] = probs[n];
                    const auto& [c, prob_sum, sq_prob_sum] = prob[n];
                    total_c += c;
                    total_prob_sum += prob_sum;
                    total_sq_prob_sum += sq_prob_sum;
                }
            }
            bb_bases.resize(max_size);

            // convert prob sums and sq prob sums to mean and variance
            common::logger->trace("Fitting beta distributions");
            ProgressBar progress_bar(probs.size(), "Fitting", std::cerr,
                                     !common::get_verbose());
#pragma omp parallel for num_threads(num_threads) schedule(static)
            for (size_t n = 0; n < probs.size(); ++n) {
                ++progress_bar;
                const auto& [c, prob_sum, sq_prob_sum] = probs[n];

                if (c == 0)
                    continue;

                auto& [bb_alpha, bb_beta, bb_base] = bb_bases[n];
                long double mn = prob_sum / c;
                long double nu = n; // using the mean-sample size parametrization
                // long double var = sq_prob_sum / c - mn * mn;
                // long double nu = mn * (1.0 - mn) / var - 1.0; // MoM estimate
                bb_alpha = nu * mn;
                bb_beta = nu * (1.0 - mn);
                long double bb_null = bb_alpha + bb_beta;

                bb_base = (lgammal(bb_null) - lgammal(bb_alpha) - lgammal(bb_beta))
                        / logl(2.0);

                common::logger->trace("n: {}\talpha: {}\tbeta: {}", n, bb_alpha, bb_beta);
            }


            // common::logger->trace("Beta:\tmean: {}\tvar: {}\talpha: {}\tbeta: {}", mn,
            //                       var, bb_alpha, bb_beta);
            // nb_base += bb_alpha*log2l(0.1) + bb_beta*log2l(0.1)-bb_null*log2l(0.1);

            // auto get_ab = [&](long double n, long double n2) {
            //     long double mu = n / nelem;
            //     long double var = n2 / nelem - mu * mu;
            //     long double a = mu * (n * mu - mu * mu - var) / n / (var - mu + mu * mu
            //     / n); long double b = a * (n / mu - 1); return std::make_pair(a, b);
            // };

            // bb_params_a = get_ab(in_kmers, in_sq_kmers);
            // bb_params_b = get_ab(out_kmers, out_sq_kmers);
            // bb_params_null = get_ab(total_kmers, total_sq_kmers);

            // const auto &[a_a, b_a] = bb_params_a;
            // const auto &[a_b, b_b] = bb_params_b;
            // const auto &[a_n, b_n] = bb_params_null;

            // long double nx = static_cast<long double>(in_kmers);
            // long double ny = static_cast<long double>(out_kmers);

            // bb_base = (lgammal(nx + 1) + lgammal(ny + 1) - lgammal(nx + ny + 1)
            //             - lgammal(a_a + b_a + nx) - lgammal(a_b + b_b + ny)
            //             + lgammal(a_n + b_n + nx + ny) + lgammal(a_n) + lgammal(b_n) - lgammal(a_n + b_n)
            //             - lgammal(a_a) - lgammal(b_a) + lgammal(a_a + b_a)
            //             - lgammal(a_b) - lgammal(b_b) + lgammal(a_b + b_b)) / logl(2.0);
        }

        std::pair<long double, long double> gp_params_a;
        std::pair<long double, long double> gp_params_b;
        std::vector<std::vector<long double>> gp_pvals;
        // std::pair<long double, long double> gp_params_n;
        // long double gp_base = 0.0;
        if (config.test_type == "gpoisson") {
            auto& [lambda_a, omega_a] = gp_params_a;
            auto& [lambda_b, omega_b] = gp_params_b;
            // auto &[lambda_n, omega_n] = gp_params_n;

            lambda_a = static_cast<long double>(in_kmers) / nelem;
            lambda_b = static_cast<long double>(out_kmers) / nelem;
            // lambda_n = static_cast<long double>(total_kmers) / nelem;

            common::logger->trace("Calculating row sums");
            tsl::hopscotch_map<size_t, size_t> row_sums_a;
            tsl::hopscotch_map<size_t, size_t> row_sums_b;
            VectorSet<size_t> row_sums_n;
            // tsl::hopscotch_map<size_t, size_t> row_sums_n;
            {
                std::vector<tsl::hopscotch_map<size_t, size_t>> row_sums_a_m(num_threads + 1);
                std::vector<tsl::hopscotch_map<size_t, size_t>> row_sums_b_m(num_threads + 1);
                std::vector<VectorSet<size_t>> row_sums_n_m(num_threads + 1);
                generate_clean_rows([&](auto, const auto& row, size_t bucket_idx) {
                    size_t in_sum = 0;
                    size_t out_sum = 0;
                    for (const auto& [j, c] : row) {
                        if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                            in_sum += c;

                        if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                            out_sum += c;
                    }
                    ++row_sums_a_m[bucket_idx][in_sum];
                    ++row_sums_b_m[bucket_idx][out_sum];
                    row_sums_n_m[bucket_idx].emplace(in_sum + out_sum);
                    // ++row_sums_n_m[bucket_idx][in_sum + out_sum];
                });

                row_sums_a = std::move(row_sums_a_m[0]);
                row_sums_b = std::move(row_sums_b_m[0]);
                row_sums_n = std::move(row_sums_n_m[0]);

                for (size_t i = 1; i <= num_threads; ++i) {
                    for (const auto& [k, c] : row_sums_a_m[i]) {
                        row_sums_a[k] += c;
                    }
                    for (const auto& [k, c] : row_sums_b_m[i]) {
                        row_sums_b[k] += c;
                    }
                    for (size_t k : row_sums_n_m[i]) {
                        row_sums_n.emplace(k);
                    }
                    // for (const auto &[k, c] : row_sums_n_m[i]) {
                    //     row_sums_n[k] += c;
                    // }
                }
            }

            long double var_a
                    = static_cast<long double>(in_sq_kmers) / nelem - lambda_a * lambda_a;
            long double var_b
                    = static_cast<long double>(out_sq_kmers) / nelem - lambda_b * lambda_b;
            // long double var_n = static_cast<long double>(total_sq_kmers) / nelem - lambda_n * lambda_n;

            auto get_omega = [&](long double lambda, long double var, const auto& generate) {
                long double val = sqrtl(lambda / var);
                if (1.0 <= val) {
                    common::logger->trace("Invalid omega: {}\tlambda: {}\tvar: {}",
                                          1.0 - val, lambda, var);
                    throw std::runtime_error("Fit fail");
                }

                return 1.0 - val;
                // auto [omega_min, omega_max] = boost::math::tools::bisect([&](long double omega) {
                //     long double dl = nelem / (1.0 - omega);
                //     generate([&](int64_t k, auto c) {
                //         dl += (static_cast<long double>(k) - lambda) * (k - 1) / ((1.0 - omega) * lambda + k * omega) * c;
                //     });
                //     return dl;
                // }, 0.0, 1.0 - std::numeric_limits<float>::min(), boost::math::tools::eps_tolerance<double>(5));
                // return (omega_min + omega_max) / 2.0;
            };

            omega_a = get_omega(lambda_a, var_a, [&](const auto& callback) {
                for (const auto& [k, c] : row_sums_a) {
                    callback(k, c);
                }
            });
            omega_b = get_omega(lambda_b, var_b, [&](const auto& callback) {
                for (const auto& [k, c] : row_sums_b) {
                    callback(k, c);
                }
            });
            // omega_n = get_omega(lambda_n, var_n, [&](const auto &callback) {
            //     for (const auto &[k, c] : row_sums_n) {
            //         callback(k, c);
            //     }
            // });

            common::logger->trace(
                    "Generalized Poisson In:\tlambda: {}\tomega: {}\tsample variance: "
                    "{}\tvar: {}",
                    lambda_a, omega_a, var_a, lambda_a / pow(1.0 - omega_a, 2.0));
            common::logger->trace(
                    "Generalized Poisson Out:\tlambda: {}\tomega: {}\tsample variance: "
                    "{}\tvar: {}",
                    lambda_b, omega_b, var_b, lambda_b / pow(1.0 - omega_b, 2.0));
            // common::logger->trace("Generalized Poisson Null:\tlambda: {}\tomega: {}\tsample variance: {}\tvar: {}",
            //                       lambda_n, omega_n, var_n, lambda_n / pow(1.0 - omega_n, 2.0));

            long double gp_base = (log1pl(-omega_a) + log1pl(-omega_b)
                                   - (1.0 - omega_a) * lambda_a - (1.0 - omega_b) * lambda_b)
                            / logl(2.0)
                    + log2l(lambda_a) + log2l(lambda_b);

            gp_pvals.resize(*std::max_element(row_sums_n.begin(), row_sums_n.end()) + 1);

            size_t rows_per_update = 1000;
            ProgressBar progress_bar(row_sums_n.size(), "Precomputing p-values",
                                     std::cerr, !common::get_verbose());
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
            for (size_t i = 0; i < row_sums_n.size(); ++i) {
                if (i > 0 && (i % rows_per_update) == 0)
                    progress_bar += rows_per_update;

                size_t n = row_sums_n.values_container()[i];
                size_t s = 0;
                size_t t = n;
                long double lprob = gp_base
                        + (-lgammal(s + 1) - lgamma(t + 1) - omega_a * s - omega_b * t)
                                / logl(2.0);
                // auto get_prob = [&](int64_t s) {
                //     int64_t t = n - s;
                //     long double lprob = gp_base + (-lgammal(s + 1) - lgamma(t + 1) - omega_a * s - omega_b * t) / logl(2.0);
                //     if (s == 0) {
                //         lprob -= log1pl(-omega_a) * lambda_a / logl(2.0);
                //     } else if (s > 1) {
                //         lprob += log2l((1.0 - omega_a) * lambda_a + omega_a * s) * (s - 1);
                //     }
                //     if (t == 0) {
                //         lprob -= log1pl(-omega_b) * lambda_b / logl(2.0);
                //     } else if (t > 1) {
                //         lprob += log2l((1.0 - omega_b) * lambda_b + omega_b * t) * (t - 1);
                //     }
                //     return exp2l(lprob);
                // };
                std::vector<long double> probs;
                probs.reserve(n + 1);
                long double cur_lprob = lprob - log1pl(-omega_a) * lambda_a / logl(2.0)
                        + log2l((1.0 - omega_b) * lambda_b + omega_b * t) * (t - 1);
                probs.emplace_back(exp2l(cur_lprob));
                long double prob_sum = probs.back();
                for (size_t s = 1; s <= n; ++s) {
                    int64_t t = n - s;
                    lprob += log2l(t + 1) - log2l(s) - omega_a + omega_b;
                    long double cur_lprob = lprob;
                    if (s == 0) {
                        cur_lprob -= log1pl(-omega_a) * lambda_a / logl(2.0);
                    } else if (s > 1) {
                        cur_lprob += log2l((1.0 - omega_a) * lambda_a + omega_a * s)
                                * (s - 1);
                    }
                    if (t == 0) {
                        cur_lprob -= log1pl(-omega_b) * lambda_b / logl(2.0);
                    } else if (t > 1) {
                        cur_lprob += log2l((1.0 - omega_b) * lambda_b + omega_b * t)
                                * (t - 1);
                    }
                    probs.emplace_back(exp2l(cur_lprob));
                    prob_sum += probs.back();
                }

                auto& pvals_n = gp_pvals[n];
                pvals_n.reserve(n + 1);
                for (size_t s = 0; s <= n; ++s) {
                    for (size_t sp = 0; sp <= n; ++sp) {
                        if (probs[sp] <= probs[s])
                            pvals_n[s] += probs[sp];
                    }
                    pvals_n[s] /= prob_sum;
                }
            }
            progress_bar += (row_sums_n.size() % rows_per_update);
        }

        long double pig_base_a = 0.0;
        long double pig_base_b = 0.0;
        long double pig_base_n = 0.0;
        long double lambda_a = 0.0;
        long double lambda_b = 0.0;
        long double lambda_n = 0.0;
        if (config.test_type == "pig") {
            // long double mu_a = static_cast<long double>(in_kmers) / nelem;
            // long double mu_b = static_cast<long double>(out_kmers) / nelem;

            // long double var_a = static_cast<long double>(in_sq_kmers) / nelem - mu_a * mu_a;
            // long double var_b = static_cast<long double>(out_sq_kmers) / nelem - mu_b * mu_b;

            // common::logger->trace("PIG\tIn:\tmu: {}\tvar: {}\tvar/mu: {}\tOut:\tmu: {}\tvar: {}\tvar/mu: {}",
            //                       mu_a, var_a, var_a/mu_a,
            //                       mu_b, var_b, var_b/mu_b);

            auto get_lp = [&](long double x, long double mu, long double lambda) {
                long double phi = 1.0L / lambda;
                int64_t i;
                long double p, pi1m, pi2m;
                long double twophi = phi + phi;

                long double A, B;
                long double mu2 = mu * mu;
                long double twophimu2 = twophi * mu2;

                p = (1.0L - sqrtl(1.0L + twophimu2)) / phi / mu; /* log p[0] */
                if (x == 0.0) {
                    if (p > 0) {
                        common::logger->error("mu: {}\tlambda: {}\tx: {}\tlogp: {}", mu,
                                              lambda, x, p);
                        throw std::runtime_error("Fail");
                    }
                    return p;
                }

                pi2m = expl(p); /* p[i - 2] = p[0]*/
                p = logl(mu) + p - log1pl(twophimu2) / 2.0L; /* log p[1] */
                if (x == 1.0) {
                    if (p > 0) {
                        common::logger->error("mu: {}\tlambda: {}\tx: {}\tlogp: {}", mu,
                                              lambda, x, p);
                        throw std::runtime_error("Fail");
                    }
                    return p;
                }

                pi1m = expl(p); /* p[i - 1] = p[1] */
                A = 1.0L / (1.0L + 1.0L / twophimu2); /* constant in first term */
                B = mu2 / (1.0L + twophimu2); /* constant in second term */
                for (i = 2; i <= x; i++) {
                    p = A * (1.0L - 1.5L / i) * pi1m + (B * pi2m) / (i * (i - 1));
                    pi2m = pi1m;
                    pi1m = p;
                }

                if (p > 1.0) {
                    common::logger->error("mu: {}\tlambda: {}\tx: {}\tlogp: {}", mu,
                                          lambda, x, logl(p));
                    throw std::runtime_error("Fail");
                }

                return logl(p);
            };

            std::vector<long double> lambdas(groups.size());
#pragma omp parallel for num_threads(num_threads)
            for (size_t i = 0; i < groups.size(); ++i) {
                long double mu = static_cast<long double>(sums[i]) / nelem;
                long double lambda = 0.0;
                long double val = 0.0;
                try {
                    std::tie(lambda, val) = boost::math::tools::brent_find_minima(
                            [&](long double lambda) {
                                long double val = 0.0;
                                for (const auto& [k, c] : hists[i]) {
                                    val += get_lp(k, mu, lambda) * c;
                                }
                                return -val;
                            },
                            std::numeric_limits<long double>::min(), 1.0L, 5);
                } catch (boost::wrapexcept<boost::math::evaluation_error>&) {
                }
                try {
                    auto [phi, val2] = boost::math::tools::brent_find_minima(
                            [&](long double phi) {
                                long double val = 0.0;
                                for (const auto& [k, c] : hists[i]) {
                                    val += get_lp(k, mu, 1.0L / phi) * c;
                                }
                                return -val;
                            },
                            std::numeric_limits<long double>::min(), 1.0L, 5);
                    if (val == 0.0 || val2 < val)
                        lambda = 1.0L / phi;
                } catch (boost::wrapexcept<boost::math::evaluation_error>&) {
                }
                if (val == 0.0 || lambda == 0.0) {
                    throw std::runtime_error("FAILS");
                }

                long double var = static_cast<long double>(sq_sums[i]) / nelem - mu * mu;
                common::logger->trace(
                        "i: {}\tmu: {}\tsample var: {}\tlambda: {}\tvar: {}", i, mu, var,
                        lambda, mu + mu * mu * mu / lambda);
                lambdas[i] = lambda;
            }

            size_t min_sum_a = std::numeric_limits<size_t>::max();
            size_t min_sum_b = std::numeric_limits<size_t>::max();
            size_t min_elem_a = 0;
            size_t min_elem_b = 0;
            for (size_t i = 0; i < lambdas.size(); ++i) {
                if (groups[i] == Group::IN && sums[i] < min_sum_a) {
                    min_sum_a = sums[i];
                    min_elem_a = i;
                }

                if (groups[i] == Group::OUT && sums[i] < min_sum_b) {
                    min_sum_b = sums[i];
                    min_elem_b = i;
                }
            }
            size_t min_sum_n = std::min(min_sum_a, min_sum_b);
            size_t min_elem_n = min_sum_a < min_sum_b ? min_elem_a : min_elem_b;

            long double min_lambda_a = lambdas[min_elem_a];
            long double min_lambda_b = lambdas[min_elem_b];
            long double min_lambda_n = lambdas[min_elem_n];

            lambda_a = 0.0;
            lambda_b = 0.0;
            lambda_n = 0.0;
            for (size_t i = 0; i < lambdas.size(); ++i) {
                if (groups[i] == Group::IN)
                    lambda_a += pow(static_cast<long double>(sums[i]) / min_sum_a, 2.0);

                if (groups[i] == Group::OUT)
                    lambda_b += pow(static_cast<long double>(sums[i]) / min_sum_b, 2.0);

                lambda_n += pow(static_cast<long double>(sums[i]) / min_sum_n, 2.0);
            }
            lambda_a *= min_lambda_a;
            lambda_b *= min_lambda_b;
            lambda_n *= min_lambda_n;

            common::logger->trace("PIG:\tIN:\tmu: {}\tlambda: {}",
                                  static_cast<long double>(in_kmers) / nelem, lambda_a);
            common::logger->trace("PIG:\tOUT:\tmu: {}\tlambda: {}",
                                  static_cast<long double>(out_kmers) / nelem, lambda_b);
            common::logger->trace("PIG:\tNULL:\tmu: {}\tlambda: {}",
                                  static_cast<long double>(total_kmers) / nelem, lambda_n);

            // long double phi_a = (var_a - mu_a) / mu_a / mu_a / mu_a;
            // long double phi_b = (var_b - mu_b) / mu_b / mu_b / mu_b;


            // long double lambda_a = 1.0 / phi_a;
            // long double lambda_b = 1.0 / phi_b;

            // common::logger->trace("PIG\tIn:\tlambda: {}\tOut:\tlambda: {}", lambda_a, lambda_b);

            // long double min_mu = std::min(mu_a, mu_b);
            // long double w_a = mu_a / min_mu;
            // long double w_b = mu_b / min_mu;
            // long double min_lambda = min_mu == mu_a ? lambda_a : lambda_b;

            // long double mu_n = mu_a + mu_b;
            // long double lambda_n = min_lambda * pow(w_a + w_b, 2.0);
            // // throw std::runtime_error("Check");

            // common::logger->trace("PIG\tNull:\tmu: {}\tlambda: {}", mu_n, lambda_n);

            // pig_base_a = 0.5 * (logl(2.0) - logl(M_PI) + logl(lambda_a)) + lambda_a / mu_a;
            // pig_base_b = 0.5 * (logl(2.0) - logl(M_PI) + logl(lambda_b)) + lambda_b / mu_b;
            // pig_base_n = 0.5 * (logl(2.0) - logl(M_PI) + logl(lambda_n)) + lambda_n / mu_n;


            // long double var_n = static_cast<long double>(total_sq_kmers) / nelem - mu_n
            // * mu_n; long double lambda = mu_n * mu_n * mu_n / (var_n - mu_n);

            // common::logger->trace("PIG: lambda: {}", lambda);
        }

        bool fdr = config.test_type == "nbinom_exact" || config.test_type == "pig"
                || config.test_type == "mwu" || config.test_type == "bm"
                || config.test_type == "bb_hypergeometric" || config.test_type == "cmh";

        // precompute minimal attainable p-values
        std::vector<std::pair<long double, size_t>> m;
        std::vector<size_t> offsets;
        if (config.test_type != "notest" && !fdr) {
            common::logger->trace("Computing minimum p-values");
            std::vector<std::vector<std::pair<long double, size_t>>> ms(num_threads + 1);
            std::vector<VectorMap<long double, size_t>> ms_map(num_threads + 1);
            std::vector<std::vector<std::vector<long double>>> gp_pvals_m(num_threads + 1);
            std::mutex mu;
            generate_clean_rows([&](uint64_t row_i, const auto& row, uint64_t bucket_idx) {
                long double pval = 1.0;
                size_t n = 0;
                if (config.test_type == "poisson_binom") {
                    for (const auto& [j, c] : row) {
                        n += static_cast<int64_t>(groups[j] != Group::OTHER)
                                + (groups[j] == Group::BOTH);
                    }
                    if (ms[bucket_idx].size() <= n)
                        ms[bucket_idx].resize(n + 1);

                    if (ms[bucket_idx][n].second == 0) {
                        size_t front = n - std::min(n, num_labels_out);
                        pval = std::min(pb_pvals[n][front], pb_pvals[n].back());
                        ms[bucket_idx][n].first = pval;
                        ms[bucket_idx][n].second = 1;
                    } else {
                        ++ms[bucket_idx][n].second;
                    }

                } else if (config.test_type == "poisson_exact") {
                    sdsl::bit_vector found(num_labels_in + num_labels_out);
                    size_t n = 0;
                    for (const auto& [j, c] : row) {
                        n += c;
                        found[j] = true;
                    }

                    auto [mu1, mu2] = get_lambda(bucket_idx, [&](const auto& callback) {
                        for (const auto& [j, c] : row) {
                            callback(j, c);
                        }
                        for (size_t j = 0; j < found.size(); ++j) {
                            if (!found[j])
                                callback(j, 0);
                        }
                    });

                    if (n > 0) {
                        if (ms[bucket_idx].size() <= n)
                            ms[bucket_idx].resize(n + 1);

                        if (ms[bucket_idx][n].second == 0) {
                            long double p = mu1 / (mu1 + mu2);
                            boost::math::binomial bdist(n, p);
                            long double pval0 = boost::math::pdf(bdist, 0);
                            long double pvaln = boost::math::pdf(bdist, n);
                            pval = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                            ms[bucket_idx][n].first = pval;
                            ms[bucket_idx][n].second = 1;
                        } else {
                            ++ms[bucket_idx][n].second;
                        }
                    }
                } else if (config.test_type == "nbinom_exact") {
                    for (const auto& [j, c] : row) {
                        n += c;
                    }

                    if (n > 0) {
                        if (ms[bucket_idx].size() <= n)
                            ms[bucket_idx].resize(n + 1);

                        if (ms[bucket_idx][n].second == 0) {
                            auto [r_n, p_n] = nb_params_null;
                            auto [r_a, p_a] = nb_params_a;
                            auto [r_b, p_b] = nb_params_b;

                            long double l21p = log1p(-p_n) / log2l(2.0);
                            long double base = nb_base
                                    + (lgammal(n + 1) - lgammal(r_n + n)) / log2l(2.0)
                                    - n * l21p;
                            auto get_pval = [&](int64_t s) {
                                int64_t t = n - s;
                                long double sbase = (lgammal(r_a + s) + lgammal(r_b + t)
                                                     - lgammal(s + 1) - lgammal(t + 1)
                                                     + s * log1pl(-p_a) + t * log1pl(-p_b))
                                        / logl(2.0);
                                return exp2l(base + sbase);
                            };

                            long double pval0 = get_pval(0);
                            long double pvaln = get_pval(n);
                            pval = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                            ms[bucket_idx][n].first = pval;
                            ms[bucket_idx][n].second = 1;
                        } else {
                            ++ms[bucket_idx][n].second;
                        }
                    }
                } else if (config.test_type == "hypergeometric"
                           || config.test_type == "poisson_bayes") {
                    size_t n = 0;
                    for (const auto& [j, c] : row) {
                        n += c;
                    }
                    if (n > 0) {
                        // if (ms[bucket_idx].size() <= n)
                        //     ms[bucket_idx].resize(n + 1);

                        // if (ms[bucket_idx][n].second == 0) {
                        sdsl::bit_vector found(groups.size());
                        auto [mu1, mu2] = get_lambda(bucket_idx, [&](const auto& callback) {
                            for (const auto& [j, c] : row) {
                                callback(j, c);
                            }
                            for (size_t j = 0; j < found.size(); ++j) {
                                if (!found[j])
                                    callback(j, 0);
                            }
                        });
                        boost::math::binomial_distribution dist(static_cast<long double>(n),
                                                                mu1 / (mu1 + mu2));
                        // boost::math::hypergeometric_distribution dist(mu1 * nelem, n, (mu1 + mu2) * nelem);
                        long double pval0 = boost::math::pdf(dist, 0);
                        long double pvaln = boost::math::pdf(dist, n);
                        pval = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                        ++ms_map[bucket_idx][pval];
                        // ms[bucket_idx][n].first = pval;
                        // ms[bucket_idx][n].second = 1;
                        // } else {
                        //    ++ms[bucket_idx][n].second;
                        // }
                    }
                } else if (config.test_type == "pig") {
                    size_t n = 0;
                    for (const auto& [j, c] : row) {
                        n += c;
                    }

                    if (n > 0) {
                        if (ms[bucket_idx].size() <= n)
                            ms[bucket_idx].resize(n + 1);

                        if (ms[bucket_idx][n].second == 0) {
                            // https://search.r-project.org/CRAN/refmans/actuar/html/PoissonInverseGaussian.html
                            long double mu_a = static_cast<long double>(in_kmers) / nelem;
                            long double mu_b = static_cast<long double>(out_kmers) / nelem;

                            // long double var_a = static_cast<long double>(in_sq_kmers) /
                            // nelem - mu_a * mu_a; long double var_b = static_cast<long
                            // double>(out_sq_kmers) / nelem - mu_b * mu_b;

                            // long double phi_a = (var_a - mu_a) / mu_a / mu_a / mu_a;
                            // long double phi_b = (var_b - mu_b) / mu_b / mu_b / mu_b;
                            // long double lambda_a = 1.0 / phi_a;
                            // long double lambda_b = 1.0 / phi_b;

                            auto get_lp = [&](long double x, long double mu,
                                              long double lambda, long double base) {
                                long double phi = 1.0L / lambda;
                                int64_t i;
                                long double p, pi1m, pi2m;
                                long double twophi = phi + phi;

                                long double A, B;
                                long double mu2 = mu * mu;
                                long double twophimu2 = twophi * mu2;

                                p = (1.0L - sqrtl(1.0L + twophimu2)) / phi
                                        / mu; /* log p[0] */
                                if (x == 0.0) {
                                    if (p > 0) {
                                        common::logger->error(
                                                "mu: {}\tlambda: {}\tn: {}\tx: {}\tlogp: "
                                                "{}",
                                                mu, lambda, n, x, p);
                                        throw std::runtime_error("Fail");
                                    }
                                    return p;
                                }

                                pi2m = expl(p); /* p[i - 2] = p[0]*/
                                p = logl(mu) + p - log1pl(twophimu2) / 2.0L; /* log p[1] */
                                if (x == 1.0) {
                                    if (p > 0) {
                                        common::logger->error(
                                                "mu: {}\tlambda: {}\tn: {}\tx: {}\tlogp: "
                                                "{}",
                                                mu, lambda, n, x, p);
                                        throw std::runtime_error("Fail");
                                    }
                                    return p;
                                }

                                pi1m = expl(p); /* p[i - 1] = p[1] */
                                A = 1.0L / (1.0L + 1.0L / twophimu2); /* constant in first term */
                                B = mu2 / (1.0L + twophimu2); /* constant in second term */
                                for (i = 2; i <= x; i++) {
                                    p = A * (1.0L - 1.5L / i) * pi1m
                                            + (B * pi2m) / (i * (i - 1));
                                    pi2m = pi1m;
                                    pi1m = p;
                                }

                                if (p > 1.0) {
                                    common::logger->error(
                                            "mu: {}\tlambda: {}\tn: {}\tx: {}\tlogp: {}",
                                            mu, lambda, n, x, logl(p));
                                    throw std::runtime_error("Fail");
                                }

                                return logl(p);
                                // long double lprob = 0.0;
                                // if (x == 0) {
                                //     // https://github.com/cran/actuar/blob/acf6c71ebef1374230b1c24789ce4a925134142f/src/poisinvgauss.c#L66
                                //     long double mu2 = mu * mu;
                                //     long double phi = 1.0 / lambda;
                                //     long double twophi = phi + phi;
                                //     long double twophimu2 = twophi * mu2;

                                //     lprob = (1.0 - sqrtl(1.0 + twophimu2))/phi/mu;
                                // } else {
                                //     lprob = base - lgamma(x + 1)
                                //             - 0.5*(x - 0.5) * (logl(2.0) - logl(lambda)
                                //             + logl(1.0 + lambda / 2 / mu / mu))
                                //             + logl(cyl_bessel_kl(x - 0.5, sqrtl(lambda
                                //             * 2 * (1.0 + lambda / 2 / mu / mu))));
                                // }
                                // if (lprob > 0) {
                                //     common::logger->error("mu: {}\tlambda: {}\tx:
                                //     {}\tn: {}\tlog(pval): {}\tpval: {}",
                                //         mu, lambda, x, n, lprob, expl(lprob));
                                //     throw std::runtime_error("pval > 1");
                                // }
                                // return lprob;
                            };

                            // long double min_mu = std::min(mu_a, mu_b);
                            // long double w_a = mu_a / min_mu;
                            // long double w_b = mu_b / min_mu;
                            // long double min_lambda = min_mu == mu_a ? lambda_a : lambda_b;

                            long double mu_n = mu_a + mu_b;
                            // long double lambda_n = min_lambda * pow(w_a + w_b, 2.0);

                            long double base_lprob = get_lp(n, mu_n, lambda_n, pig_base_n);
                            auto get_lprob = [&](int64_t s) {
                                return get_lp(s, mu_a, lambda_a, pig_base_a)
                                        + get_lp(n - s, mu_b, lambda_b, pig_base_b)
                                        - base_lprob;
                            };

                            long double lpval0 = get_lprob(0);
                            long double lpvaln = get_lprob(n);

                            pval = expl(std::min(lpval0, lpvaln)) * (1 + (lpval0 == lpvaln));
                            ms[bucket_idx][n].first = pval;
                            ms[bucket_idx][n].second = 1;
                        } else {
                            ++ms[bucket_idx][n].second;
                        }
                    }
                } else if (config.test_type == "gpoisson") {
                    size_t n = 0;
                    for (const auto& [j, c] : row) {
                        n += c;
                    }

                    if (n > 0) {
                        if (ms[bucket_idx].size() <= n) {
                            gp_pvals_m[bucket_idx].resize(n + 1);
                            ms[bucket_idx].resize(n + 1);
                        }

                        if (ms[bucket_idx][n].second == 0) {
                            // const auto &[lambda_a, omega_a] = gp_params_a;
                            // const auto &[lambda_b, omega_b] = gp_params_b;
                            // size_t s = 0;
                            // size_t t = n;
                            // long double lprob = gp_base + (-lgammal(s + 1) - lgamma(t + 1) - omega_a * s - omega_b * t) / logl(2.0);
                            // // auto get_prob = [&](int64_t s) {
                            // //     int64_t t = n - s;
                            // //     long double lprob = gp_base + (-lgammal(s + 1) - lgamma(t + 1) - omega_a * s - omega_b * t) / logl(2.0);
                            // //     if (s == 0) {
                            // //         lprob -= log1pl(-omega_a) * lambda_a / logl(2.0);
                            // //     } else if (s > 1) {
                            // //         lprob += log2l((1.0 - omega_a) * lambda_a + omega_a * s) * (s - 1);
                            // //     }
                            // //     if (t == 0) {
                            // //         lprob -= log1pl(-omega_b) * lambda_b / logl(2.0);
                            // //     } else if (t > 1) {
                            // //         lprob += log2l((1.0 - omega_b) * lambda_b + omega_b * t) * (t - 1);
                            // //     }
                            // //     return exp2l(lprob);
                            // // };
                            // std::vector<long double> probs;
                            // probs.reserve(n + 1);
                            // long double cur_lprob = lprob - log1pl(-omega_a) * lambda_a / logl(2.0)
                            //                             + log2l((1.0 - omega_b) * lambda_b + omega_b * t) * (t - 1);
                            // probs.emplace_back(exp2l(cur_lprob));
                            // long double prob_sum = probs.back();
                            // for (size_t s = 1; s <= n; ++s) {
                            //     int64_t t = n - s;
                            //     lprob += log2l(t + 1) - log2l(s) - omega_a + omega_b;
                            //     long double cur_lprob = lprob;
                            //     if (s == 0) {
                            //         cur_lprob -= log1pl(-omega_a) * lambda_a / logl(2.0);
                            //     } else if (s > 1) {
                            //         cur_lprob += log2l((1.0 - omega_a) * lambda_a + omega_a * s) * (s - 1);
                            //     }
                            //     if (t == 0) {
                            //         cur_lprob -= log1pl(-omega_b) * lambda_b / logl(2.0);
                            //     } else if (t > 1) {
                            //         cur_lprob += log2l((1.0 - omega_b) * lambda_b + omega_b * t) * (t - 1);
                            //     }
                            //     probs.emplace_back(exp2l(cur_lprob));
                            //     prob_sum += probs.back();
                            // }

                            const auto& pvals_n = gp_pvals[n];
                            // pvals_n.reserve(n + 1);
                            // for (size_t s = 0; s <= n; ++s) {
                            //     for (size_t sp = 0; sp <= n; ++sp) {
                            //         if (probs[sp] <= probs[s])
                            //             pvals_n[s] += probs[sp];
                            //     }
                            //     pvals_n[s] /= prob_sum;
                            // }

                            pval = std::min(pvals_n[0], pvals_n[n]);
                            ms[bucket_idx][n].first = pval;
                            ms[bucket_idx][n].second = 1;
                        } else {
                            ++ms[bucket_idx][n].second;
                        }
                    }
                    // } else if (config.test_type == "bb_hypergeometric") {
                    //     size_t n = 0;
                    //     for (const auto& [j, c] : row) {
                    //         n += c;
                    //     }
                    //     if (n > 0) {
                    //         if (ms[bucket_idx].size() <= n)
                    //             ms[bucket_idx].resize(n + 1);

                    //         if (ms[bucket_idx][n].second == 0) {
                    //             const auto& [a_a, b_a] = bb_params_a;
                    //             const auto& [a_b, b_b] = bb_params_b;
                    //             const auto& [a_n, b_n] = bb_params_null;
                    //             long double nx = static_cast<long double>(in_kmers);
                    //             long double ny = static_cast<long double>(out_kmers);
                    //             long double base
                    //                     = (lgammal(n + 1) + lgammal(nx + ny - n + 1)
                    //                        - lgammal(a_n + n) - lgammal(b_n + nx + ny - n))
                    //                             / logl(2.0)
                    //                     + bb_base;
                    //             int64_t s = 0;
                    //             long double pval0
                    //                     = (-lgammal(s + 1) - lgammal(nx - s + 1)
                    //                        - lgammal(n - s + 1) - lgammal(ny - n + s + 1)
                    //                        + lgammal(a_a + s) + lgammal(b_a + nx - s)
                    //                        + lgammal(a_b + n - s) + lgammal(b_b + ny - n + s))
                    //                             / logl(2.0)
                    //                     + base;
                    //             s = n;
                    //             long double pvaln
                    //                     = (-lgammal(s + 1) - lgammal(nx - s + 1)
                    //                        - lgammal(n - s + 1) - lgammal(ny - n + s + 1)
                    //                        + lgammal(a_a + s) + lgammal(b_a + nx - s)
                    //                        + lgammal(a_b + n - s) + lgammal(b_b + ny - n + s))
                    //                             / logl(2.0)
                    //                     + base;
                    //             pval = exp2l(std::min(pval0, pvaln)) * (1 + (pval0 == pvaln));
                    //             ms[bucket_idx][n].first = pval;
                    //             ms[bucket_idx][n].second = 1;
                    //         } else {
                    //             ++ms[bucket_idx][n].second;
                    //         }
                    //     }
                }

                // if (config.test_by_unitig) {
                //     node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
                //     size_t unitig_id = kmer_to_unitig[node];
                //     std::lock_guard<std::mutex> lock(mu);
                //     ++counts_map[unitig_id][n];
                // }
            });

            // if (config.test_by_unitig) {
            //     common::logger->trace("Combining minimal attainable p-values for unitigs");
            //     offsets.reserve(num_threads + 1);
            //     offsets.emplace_back(0);
            //     for (const auto &bucket : pvals_min_buckets) {
            //         offsets.emplace_back(offsets.back() + bucket.size());
            //     }
            //     std::mutex m_vec_mu;
            //     clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
            //         std::vector<long double> pvals;
            //         pvals.reserve(path.size());
            //         {
            //             std::lock_guard<std::mutex> lock(m_vec_mu);
            //             for (node_index node : path) {
            //                 size_t row = AnnotatedDBG::graph_to_anno_index(node);
            //                 size_t row_rank = kept.rank1(row);
            //                 size_t bucket_idx = std::lower_bound(offsets.begin(), offsets.end(), row_rank) - offsets.begin() - 1;
            //                 size_t offset = row_rank - offsets[bucket_idx];
            //                 pvals.emplace_back(get<PValStorage>(pvals_min_buckets[bucket_idx][offset]));
            //             }
            //         }
            //         double comp_min_pval = combine_pvals(pvals);
            //         size_t row = AnnotatedDBG::graph_to_anno_index(path[0]);
            //         size_t row_rank = kept.rank1(row);
            //         size_t bucket_idx = std::lower_bound(offsets.begin(), offsets.end(), row_rank) - offsets.begin() - 1;
            //         size_t offset = row_rank - offsets[bucket_idx];
            //         ++ms[bucket_idx][comp_min_pval];
            //         std::lock_guard<std::mutex> lock(m_vec_mu);
            //         set<PValStorage>(pvals_min_buckets[bucket_idx][offset], comp_min_pval);
            //     });
            // }

            if (config.test_type == "poisson_bayes") {
                VectorMap<long double, size_t> m_map = std::move(ms_map[0]);
                for (size_t i = 1; i < ms_map.size(); ++i) {
                    for (const auto& [pval, c] : ms_map[i]) {
                        m_map[pval] += c;
                    }
                }
                m = const_cast<std::vector<std::pair<long double, size_t>>&&>(
                        m_map.values_container());
                std::sort(m.rbegin(), m.rend());
            } else {
                std::swap(m, ms[0]);
                // std::swap(gp_pvals, gp_pvals_m[0]);
                size_t max_size = m.size();
                for (size_t i = 1; i < ms.size(); ++i) {
                    max_size = std::max(max_size, ms[i].size());
                }
                m.resize(max_size);
                // gp_pvals.resize(max_size);
                for (size_t i = 1; i < ms.size(); ++i) {
                    for (size_t n = 0; n < ms[i].size(); ++i) {
                        const auto& [pval, c] = ms[i][n];
                        m[n].first = pval;
                        m[n].second += c;
                        // if (gp_pvals[n].empty() && gp_pvals_m[i][n].size())
                        //     std::swap(gp_pvals[n], gp_pvals_m[i][n]);
                    }
                }

                if (config.test_type == "poisson_binom")
                    std::sort(m.rbegin(), m.rend());

                // if (config.test_type == "gpoisson") {
                //     m.erase(std::remove_if(m.begin(), m.end(), [&](const auto &a) {
                //     return !a.second; }), m.end()); std::sort(m.rbegin(), m.rend());
                // }
            }

            if (m.size() > 1) {
                auto it = m.begin();
                auto jt = it + 1;
                while (jt != m.end()) {
                    while (it->second == 0 && jt != m.end()) {
                        ++it;
                        ++jt;
                    }
                    while (jt != m.end() && jt->second == 0) {
                        ++jt;
                    }
                    if (jt == m.end())
                        break;

                    if (it->first - 1.0L > 1e-5L) {
                        common::logger->error("pval it > 1.0: {} ({})", it->first,
                                              it->first - 1.0L);
                        throw std::runtime_error("p-val fail");
                    }

                    if (jt->first - 1.0L > 1e-5L) {
                        common::logger->error("pval jt > 1.0: {} ({})", jt->first,
                                              jt->first - 1.0L);
                        throw std::runtime_error("p-val fail");
                    }

                    if (jt->first - it->first > 1e-5) {
                        common::logger->error(
                                "min p-vals not sorted: {} vs. {}\t{} !>= {}",
                                it - m.begin(), jt - m.begin(), it->first, jt->first);
                        throw std::runtime_error("p-val fail");
                    }
                    ++it;
                    ++jt;
                }
            }

            // if (config.test_by_unitig) {
            //     common::logger->trace("Correcting for monotigs");
            //     for (const auto &count_map : counts_map) {
            //         for (const auto &[n, c] : count_map) {
            //             if (c > 1)
            //                 m[n].second -= c - 1;
            //         }
            //     }
            // }
        }

        // determine cutoffs for multiple testing correction
        auto [k_min, k, n_cutoff] = correct_pvals(m);
        if (fdr) {
            k_min = 1;
            k = 1;
            n_cutoff = 1;
        }
        // n_cutoff = 1;
        // k = 1;

        common::logger->trace("Picked: k: {}\tn: {} / {}\tpval_min: {}\tn_cutoff: {}",
                              k_min, k, nelem, config.family_wise_error_rate / k, n_cutoff);

        common::logger->trace("Running differential tests");
        std::vector<std::pair<long double, size_t>> nb_pvals;
        if (fdr)
            nb_pvals.resize(kept.num_set_bits(), std::make_pair(1.0, kept.size()));

        std::atomic_thread_fence(std::memory_order_release);
        generate_clean_rows([&](uint64_t row_i, const auto& row, uint64_t bucket_idx) {
            if (config.test_type == "notest") {
                size_t count_in = 0;
                size_t count_out = 0;
                for (const auto& [j, c] : row) {
                    if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                        ++count_out;

                    if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                        ++count_in;
                }

                bool in_kmer = count_in >= config.min_in_recurrence
                        && count_in <= config.max_in_recurrence;
                bool out_kmer = count_out >= config.min_out_recurrence
                        && count_out <= config.max_out_recurrence;

                node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
                if (in_kmer)
                    set_bit(indicator_in.data(), node, true, std::memory_order_relaxed);

                if (out_kmer)
                    set_bit(indicator_out.data(), node, true, std::memory_order_relaxed);

                return;
            }

            auto set_pval = [&](long double pval, long double eff_size) {
                // if (!config.test_by_unitig) {
                if (pval * k < config.family_wise_error_rate) {
                    if (eff_size == 0) {
                        common::logger->error("Effect size 0 when p-value is {}", pval);
                        throw std::runtime_error("Pval failure");
                    }

                    node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
                    if (eff_size > 0) {
                        set_bit(indicator_in.data(), node, true, std::memory_order_relaxed);
                    } else if (eff_size < 0) {
                        set_bit(indicator_out.data(), node, true, std::memory_order_relaxed);
                    }
                }
                // } else {
                //     push_back(pvals_buckets[bucket_idx], pval);
                //     push_back(eff_size_buckets[bucket_idx], eff_size);
                // }
            };

            long double eff_size = 0.0;
            long double pval = 1.0;
            if (config.test_type == "poisson_binom") {
                size_t n = 0;
                uint64_t s = 0;
                for (const auto& [j, c] : row) {
                    n += static_cast<int64_t>(groups[j] != Group::OTHER)
                            + (groups[j] == Group::BOTH);
                    s += (groups[j] == Group::IN);
                }

                // if (n < n_cutoff)
                //     return;

                size_t front = n - std::min(n, num_labels_out);
                double min_pval = std::min(pb_pvals[n][front], pb_pvals[n].back());

                if (min_pval * k_min >= config.family_wise_error_rate)
                    return;

                pval = pb_pvals[n][s];
                eff_size = static_cast<double>(s) - mid_points[n];
                if (fdr) {
                    size_t nb_idx = kept.rank1(row_i) - 1;
                    nb_pvals[nb_idx].second = row_i;
                    if (pval < config.family_wise_error_rate) {
                        nb_pvals[nb_idx].first = pval;
                    }
                }
            } else if (config.test_type == "poisson_exact") {
                size_t n = 0;
                sdsl::bit_vector found(num_labels_in + num_labels_out);
                int64_t in_sum = 0;
                int64_t out_sum = 0;
                for (const auto& [j, c] : row) {
                    n += c;
                    found[j] = true;
                    if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                        out_sum += c;

                    if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                        in_sum += c;
                }

                auto [mu1, mu2] = get_lambda(bucket_idx, [&](const auto& callback) {
                    for (const auto& [j, c] : row) {
                        callback(j, c);
                    }
                    for (size_t j = 0; j < found.size(); ++j) {
                        if (!found[j])
                            callback(j, 0);
                    }
                });

                if (n < n_cutoff)
                    return;

                long double p = mu1 / (mu1 + mu2);
                boost::math::binomial bdist(n, p);
                if (n != 0) {
                    // double pval0 = boost::math::pdf(bdist, 0);
                    // double pvaln = boost::math::pdf(bdist, n);
                    // double min_pval = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                    // if (!config.test_by_unitig && min_pval * k_min >= config.family_wise_error_rate)
                    //     return;

                    auto get_deviance = [&](long double y, long double mu) {
                        y += 1e-8;
                        mu += 1e-8;
                        return 2 * (y * log(y / mu) - y + mu);
                    };

                    std::vector<long double> devs;
                    devs.reserve(n + 1);
                    for (size_t s = 0; s <= n; ++s) {
                        devs.emplace_back(get_deviance(s, mu1) + get_deviance(n - s, mu2));
                    }

                    pval = 0.0;

                    size_t s = 0;
                    for (; s <= n; ++s) {
                        if (devs[s] < devs[in_sum])
                            break;
                    }
                    if (s > 0)
                        pval += boost::math::cdf(bdist, s - 1);

                    size_t sp = n;
                    for (; sp >= s; --sp) {
                        if (devs[sp] < devs[in_sum])
                            break;
                    }

                    if (sp < n)
                        pval += boost::math::cdf(boost::math::complement(bdist, sp));

                    if (pval * k >= config.family_wise_error_rate)
                        return;

                    eff_size = static_cast<long double>(in_sum) - boost::math::mode(bdist);
                }

            } else if (config.test_type == "pig") {
                size_t n = 0;
                int64_t in_sum = 0;
                int64_t out_sum = 0;
                for (const auto& [j, c] : row) {
                    n += c;
                    if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
                        out_sum += c;
                    }

                    if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
                        in_sum += c;
                    }
                }

                if (n < n_cutoff)
                    return;

                if (n > 0) {
                    long double mu_a = static_cast<long double>(in_kmers) / nelem;
                    long double mu_b = static_cast<long double>(out_kmers) / nelem;

                    // long double var_a = static_cast<long double>(in_sq_kmers) / nelem - mu_a * mu_a;
                    // long double var_b = static_cast<long double>(out_sq_kmers) / nelem - mu_b * mu_b;

                    // long double phi_a = (var_a - mu_a) / mu_a / mu_a / mu_a;
                    // long double phi_b = (var_b - mu_b) / mu_b / mu_b / mu_b;
                    // long double lambda_a = 1.0 / phi_a;
                    // long double lambda_b = 1.0 / phi_b;

                    auto get_lprobs = [&](long double mu, long double lambda) {
                        // https://github.com/cran/actuar/blob/acf6c71ebef1374230b1c24789ce4a925134142f/src/poisinvgauss.c#L66
                        long double phi = 1.0L / lambda;
                        uint64_t i;
                        long double p, pi1m, pi2m;
                        long double twophi = phi + phi;

                        long double A, B;
                        long double mu2 = mu * mu;
                        long double twophimu2 = twophi * mu2;

                        p = (1.0L - sqrtl(1.0L + twophimu2)) / phi / mu; /* log p[0] */
                        std::vector<long double> lprobs;
                        lprobs.reserve(n);
                        lprobs.emplace_back(p);
                        // if (x == 0.0) {
                        //     if (p > 0) {
                        //         common::logger->error("mu: {}\tlambda: {}\tn: {}\tx:
                        //         {}\tlogp: {}", mu, lambda, n, x, p); throw
                        //         std::runtime_error("Fail");
                        //     }
                        //     return p;
                        // }

                        pi2m = expl(p); /* p[i - 2] = p[0]*/
                        p = logl(mu) + p - log1pl(twophimu2) / 2.0L; /* log p[1] */
                        lprobs.emplace_back(p);
                        // if (x == 1.0) {
                        //     if (p > 0) {
                        //         common::logger->error("mu: {}\tlambda: {}\tn: {}\tx:
                        //         {}\tlogp: {}", mu, lambda, n, x, p); throw
                        //         std::runtime_error("Fail");
                        //     }
                        //     return p;
                        // }

                        pi1m = expl(p); /* p[i - 1] = p[1] */
                        A = 1.0L / (1.0L + 1.0L / twophimu2); /* constant in first term */
                        B = mu2 / (1.0L + twophimu2); /* constant in second term */
                        for (i = 2; i <= n; i++) {
                            p = A * (1.0L - 1.5L / i) * pi1m + (B * pi2m) / (i * (i - 1));
                            pi2m = pi1m;
                            pi1m = p;
                            lprobs.emplace_back(logl(p));
                        }

                        return lprobs;

                        // if (p > 1.0) {
                        //     common::logger->error("mu: {}\tlambda: {}\tn: {}\tx:
                        //     {}\tlogp: {}", mu, lambda, n, x, logl(p)); throw
                        //     std::runtime_error("Fail");
                        // }

                        // return logl(p);
                    };

                    // long double min_mu = std::min(mu_a, mu_b);
                    // long double w_a = mu_a / min_mu;
                    // long double w_b = mu_b / min_mu;
                    // long double min_lambda = min_mu == mu_a ? lambda_a : lambda_b;

                    long double mu_n = mu_a + mu_b;
                    // long double lambda_n = min_lambda * pow(w_a + w_b, 2.0);

                    auto get_deviance
                            = [&](long double x, long double mu, long double lambda) {
                                  return lambda * lambda * pow(x - mu, 2.0) / x / mu / mu;
                              };

                    auto get_total_deviance = [&](size_t s) {
                        if (s == 0 || s == n)
                            return std::numeric_limits<long double>::max();

                        return get_deviance(s, mu_a, lambda_a)
                                + get_deviance(n - s, mu_b, lambda_b);
                    };

                    long double base_lprob = get_lprobs(mu_n, lambda_n).back();
                    auto lprobs_a = get_lprobs(mu_a, lambda_a);
                    auto lprobs_b = get_lprobs(mu_b, lambda_b);
                    auto get_lprob = [&](int64_t s) {
                        return lprobs_a[s] + lprobs_b[n - s] - base_lprob;
                        // return get_lp(s, mu_a, lambda_a) + get_lp(n - s, mu_b, lambda_b) - base_lprob;
                    };

                    long double in_lprob = get_lprob(in_sum);
                    long double in_deviance = get_total_deviance(in_sum);
                    pval = 0.0;
                    eff_size = static_cast<long double>(in_sum) - mu_a / mu_n * n;
                    for (size_t s = 0; s <= n; ++s) {
                        long double dev = get_total_deviance(s);
                        if (dev >= in_deviance) {
                            long double lprob = get_lprob(s);
                            pval += expl(lprob);
                            if (pval >= config.family_wise_error_rate)
                                break;
                        }
                    }
                    // size_t s = 0;
                    // for ( ; s <= n; ++s) {
                    //     long double lprob = get_lprob(s);
                    //     if (lprob <= in_lprob) {
                    //         pval += expl(lprob);
                    //     } else {
                    //         break;
                    //     }
                    // }
                    // for (size_t t = 0; t <= n; ++t) {
                    //     size_t sp = n - t;
                    //     if (sp == s)
                    //         break;
                    //     long double lprob = get_lprob(sp);
                    //     if (lprob <= in_lprob) {
                    //         pval += expl(lprob);
                    //     } else {
                    //         break;
                    //     }
                    // }

                    if (pval == 0.0) {
                        common::logger->trace(
                                "In: ({},{})\tOut: ({},{})\tNull: ({},{})\tn: "
                                "{}\tin_sum: {}\tin_lprob: {} = {} + {} - {}",
                                mu_a, lambda_a, mu_b, lambda_b, mu_n, lambda_n, n, in_sum,
                                in_lprob, lprobs_a[in_sum], lprobs_b[out_sum],
                                // get_lp(in_sum, mu_a, lambda_a),
                                // get_lp(out_sum, mu_b, lambda_b),
                                base_lprob);
                        throw std::runtime_error("Fail calc p-val");
                    }

                    size_t nb_idx = kept.rank1(row_i) - 1;
                    nb_pvals[nb_idx].first = pval;
                    nb_pvals[nb_idx].second = row_i;
                }
            } else if (config.test_type == "gpoisson") {
                size_t n = 0;
                size_t in_sum = 0;
                size_t out_sum = 0;
                for (const auto& [j, c] : row) {
                    n += c;
                    if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                        out_sum += c;

                    if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                        in_sum += c;
                }

                if (n < n_cutoff)
                    return;

                if (n != 0) {
                    size_t mode = 0;
                    size_t nmode = 0;
                    const auto& pvals = gp_pvals[n];
                    for (size_t s = 0; s <= n; ++s) {
                        if (pvals[s] >= 1.0) {
                            mode += s;
                            ++nmode;
                        }
                    }
                    long double midpoint = static_cast<long double>(mode) / nmode;
                    eff_size = static_cast<long double>(in_sum) - midpoint;
                    pval = pvals[in_sum];
                    // const auto &[lambda_a, omega_a] = gp_params_a;
                    // const auto &[lambda_b, omega_b] = gp_params_b;
                    // const auto &cprobs = prob_totals[n];
                    // size_t mode = 0;
                    // long double last_prob = probs[0];
                    // long double in_prob = probs[0];
                    // for (size_t i = 1; i < probs.size(); ++i) {
                    //     long double cur_prob = probs[i] - probs[i - 1];
                    //     if (cur_prob > last_prob) {
                    //         last_prob = cur_prob;
                    //         mode = i;
                    //     }
                    //     if (i == in_sum)
                    //         in_prob = cur_prob;
                    // }
                    // pval = 0.0;
                    // eff_size = static_cast<ssize_t>(in_sum) - mode;
                    // long double in_prob = probs[in_sum];
                    // for (size_t s = 0; s <= mode; ++s) {
                    //     if (probs[s] <= in_prob)
                    //         pval +=
                    // }
                    // const auto &[lambda_n, omega_n] = gp_params_n;
                    // long double base = gp_base + lgammal(n + 1) + omega_n * n - logl((1.0 - omega_n)*lambda_n + omega_n*n) * (n - 1);
                    // auto get_lprob = [&](int64_t s) {
                    //     int64_t t = n - s;
                    //     long double lprob = base - lgammal(s + 1) - lgammal(t + 1) - omega_a * s - omega_b * t;
                    //     if (s == 0) {
                    //         lprob -= logl(lambda_a) + log1pl(-omega_a);
                    //     } else if (s > 1) {
                    //         lprob += logl((1.0 - omega_a) * lambda_a + omega_a*s) * (s - 1);
                    //     }
                    //     if (t == 0) {
                    //         lprob -= logl(lambda_b) + log1pl(-omega_b);
                    //     } else if (t > 1) {
                    //         lprob += logl((1.0 - omega_b) * lambda_b + omega_b*t) * (t - 1);
                    //     }
                    //     return lprob;
                    // };
                    // // {
                    // //     long double lpval0 = get_lprob(0);
                    // //     long double lpvaln = get_lprob(n);
                    // //     long double min_pval = expl(std::min(lpval0, lpvaln)) * (1 + (lpval0 == lpvaln));
                    // //     if (min_pval * k_min >= config.family_wise_error_rate)
                    // //         return;
                    // // }
                    // pval = 0.0;
                    // eff_size = static_cast<long double>(in_sum) - lambda_a / lambda_n;
                    // long double inlp = get_lprob(in_sum);
                    // for (size_t s = 0; s <= n; ++s) {
                    //     long double lp = get_lprob(s);
                    //     if (lp <= inlp)
                    //         pval += expl(lp);

                    //     if (pval >= k * config.family_wise_error_rate)
                    //         break;
                    // }
                }
            } else if (config.test_type == "nbinom_exact") {
                size_t nb_idx = kept.rank1(row_i) - 1;
                size_t n = 0;
                int64_t in_sum = 0;
                int64_t out_sum = 0;
                for (const auto& [j, c] : row) {
                    n += c;
                    if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
                        out_sum += c;
                    }

                    if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
                        in_sum += c;
                    }
                }

                // if (n < n_cutoff)
                //     return;

                if (n > 0) {
                    // long double min_pval = 1.0;
                    auto [r_a, p_a] = nb_params_a;
                    auto [r_b, p_b] = nb_params_b;
                    auto [r_n, p_n] = nb_params_null;
                    long double min_pval = 1.0;
                    {
                        long double l21p = log1p(-p_n) / log2l(2.0);
                        long double base = nb_base
                                + (lgammal(n + 1) - lgammal(r_n + n)) / log2l(2.0)
                                - n * l21p;
                        auto get_pval = [&](int64_t s) {
                            int64_t t = n - s;
                            long double sbase = (lgammal(r_a + s) + lgammal(r_b + t)
                                                 - lgammal(s + 1) - lgammal(t + 1)
                                                 + s * log1pl(-p_a) + t * log1pl(-p_b))
                                    / logl(2.0);
                            return exp2l(base + sbase);
                        };

                        long double pval0 = get_pval(0);
                        long double pvaln = get_pval(n);
                        min_pval = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                    }

                    if (min_pval < config.family_wise_error_rate) {
                        nb_pvals[nb_idx].second = row_i;
                        long double midpoint = r_a * n / (r_a + r_b);

                        auto get_deviance_a = [&](long double s) {
                            if (s == 0)
                                return r_a * logl(p_a) * -2.0;

                            return ((logl(s) - log1pl(-p_a)) * s - (r_a + s) * logl(r_a + s)
                                    + r_a * (logl(r_a) - logl(p_a)))
                                    * 2.0;
                        };

                        auto get_deviance_b = [&](long double t) {
                            if (t == 0)
                                return r_b * logl(p_b) * -2.0;

                            return ((logl(t) - log1pl(-p_b)) * t - (r_b + t) * logl(r_b + t)
                                    + r_b * (logl(r_b) - logl(p_b)))
                                    * 2.0;
                        };

                        auto get_deviance = [&](long double s, long double t) {
                            return get_deviance_a(s) + get_deviance_b(t);
                        };

                        long double min_dev
                                = get_deviance(midpoint,
                                               static_cast<long double>(n) - midpoint);
                        long double in_dev = get_deviance(in_sum, out_sum);

                        if (in_dev > min_dev) {
                            eff_size = static_cast<long double>(in_sum) - midpoint;
                            pval = 0.0;
                            long double base = nb_base
                                    + (lgammal(n + 1) - lgammal(r_n + n) - n * log1pl(-p_n))
                                            / logl(2.0);
                            double l1pa = log1pl(-p_a) / logl(2.0);
                            double l1pb = log1pl(-p_b) / logl(2.0);

                            {
                                size_t s = 0;
                                size_t t = n;
                                long double sbase = (lgammal(r_a + s) + lgammal(r_b + t)
                                                     - lgammal(s + 1) - lgammal(t + 1))
                                                / logl(2.0)
                                        + s * l1pa + t * l1pb;
                                pval += exp2(base + sbase);
                                for (++s; s <= n; ++s) {
                                    long double dev = get_deviance(s, t - 1);
                                    if (dev >= in_dev) {
                                        sbase += log2l(r_a + s) - log2l(r_b + t)
                                                - log2l(s + 1) + log2l(t + 1) + l1pa - l1pb;
                                        --t;
                                        pval += exp2(base + sbase);
                                        if (pval >= config.family_wise_error_rate)
                                            break;
                                    } else {
                                        break;
                                    }
                                }
                            }
                            if (pval < config.family_wise_error_rate
                                && get_deviance(n, 0) >= in_dev) {
                                size_t s = n;
                                size_t t = 0;
                                long double sbase = (lgammal(r_a + s) + lgammal(r_b + t)
                                                     - lgammal(s + 1) - lgammal(t + 1))
                                                / logl(2.0)
                                        + s * l1pa + t * l1pb;
                                pval += exp2(base + sbase);
                                for (++t; t <= n; ++t) {
                                    long double dev = get_deviance(s - 1, t);
                                    if (dev >= in_dev) {
                                        sbase -= log2l(r_a + s) - log2l(r_b + t)
                                                - log2l(s + 1) + log2l(t + 1) + l1pa - l1pb;
                                        --s;
                                        pval += exp2(base + sbase);
                                        if (pval >= config.family_wise_error_rate)
                                            break;
                                    } else {
                                        break;
                                    }
                                }
                            }
                            nb_pvals[nb_idx].first = pval;
                        }
                    } else {
                        nb_pvals[nb_idx].second = row_i;
                        nb_pvals[nb_idx].first = 1.1;
                    }
                }
            } else if (config.test_type == "hypergeometric"
                       || config.test_type == "poisson_bayes") {
                size_t n = 0;
                int64_t in_sum = 0;
                int64_t out_sum = 0;
                for (const auto& [j, c] : row) {
                    n += c;
                    if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                        out_sum += c;

                    if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                        in_sum += c;
                }
                if (n != 0) {
                    sdsl::bit_vector found(groups.size());
                    auto [mu1, mu2] = get_lambda(bucket_idx, [&](const auto& callback) {
                        for (const auto& [j, c] : row) {
                            callback(j, c);
                        }
                        for (size_t j = 0; j < found.size(); ++j) {
                            if (!found[j])
                                callback(j, 0);
                        }
                    });
                    boost::math::binomial_distribution dist(static_cast<long double>(n),
                                                            mu1 / (mu1 + mu2));
                    // boost::math::hypergeometric_distribution dist(mu1 * nelem, n, (mu1 + mu2) * nelem);
                    eff_size = static_cast<double>(in_sum) - boost::math::mode(dist);

                    long double pval0 = boost::math::pdf(dist, 0);
                    long double pvaln = boost::math::pdf(dist, n);
                    long double pval_min = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                    if (pval_min * k_min >= config.family_wise_error_rate)
                        return;

                    long double in_prob = boost::math::pdf(dist, in_sum);
                    if (eff_size > 0) {
                        pval = boost::math::cdf(boost::math::complement(dist, in_sum - 1))
                                + boost::math::cdf(
                                        dist,
                                        floor(boost::math::quantile(
                                                boost::math::complement(dist, in_prob))));
                    } else if (eff_size < 0) {
                        pval = boost::math::cdf(dist, in_sum)
                                + boost::math::cdf(
                                        dist, ceil(boost::math::quantile(dist, in_prob)));
                    }

                    // double in_sz = abs(eff_size);

                    // pval = 0.0;

                    // size_t s = 0;
                    // for ( ; s <= n; ++s) {
                    //     double sz = abs(static_cast<double>(s) - mean);
                    //     if (sz < in_sz)
                    //         break;
                    // }
                    // if (s > 0)
                    //     pval += boost::math::cdf(dist, s - 1);

                    // size_t sp = n;
                    // for ( ; sp >= s; --sp) {
                    //     double sz = abs(static_cast<double>(sp) - mean);
                    //     if (sz < in_sz)
                    //         break;
                    // }

                    // if (sp < n)
                    //     pval += boost::math::cdf(boost::math::complement(dist, sp));

                    if (pval * k >= config.family_wise_error_rate)
                        return;
                }
            } else if (config.test_type == "mwu") {
                auto generate_a = [&](const auto& callback) {
                    sdsl::bit_vector found(groups.size());
                    for (const auto& [j, c] : row) {
                        found[j] = true;
                        if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                            callback(static_cast<long double>(c) / sums[j]);
                    }
                    for (size_t j = 0; j < found.size(); ++j) {
                        if (!found[j] && (groups[j] == Group::IN || groups[j] == Group::BOTH))
                            callback(0.0L);
                    }
                };
                auto generate_b = [&](const auto& callback) {
                    sdsl::bit_vector found(groups.size());
                    for (const auto& [j, c] : row) {
                        found[j] = true;
                        if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                            callback(static_cast<long double>(c) / sums[j]);
                    }
                    for (size_t j = 0; j < found.size(); ++j) {
                        if (!found[j] && (groups[j] == Group::OUT || groups[j] == Group::BOTH))
                            callback(0.0L);
                    }
                };

                std::tie(pval, eff_size) = mann_whitneyu(generate_a, generate_b);

                if (pval < config.family_wise_error_rate) {
                    size_t nb_idx = kept.rank1(row_i) - 1;
                    nb_pvals[nb_idx].first = pval;
                    nb_pvals[nb_idx].second = row_i;
                }
            } else if (config.test_type == "bm") {
                auto generate_a = [&](const auto& callback) {
                    sdsl::bit_vector found(groups.size());
                    for (const auto& [j, c] : row) {
                        found[j] = true;
                        if (groups[j] == Group::IN || groups[j] == Group::BOTH)
                            callback(static_cast<long double>(c) / sums[j]);
                    }
                    for (size_t j = 0; j < found.size(); ++j) {
                        if (!found[j] && (groups[j] == Group::IN || groups[j] == Group::BOTH))
                            callback(0.0L);
                    }
                };
                auto generate_b = [&](const auto& callback) {
                    sdsl::bit_vector found(groups.size());
                    for (const auto& [j, c] : row) {
                        found[j] = true;
                        if (groups[j] == Group::OUT || groups[j] == Group::BOTH)
                            callback(static_cast<long double>(c) / sums[j]);
                    }
                    for (size_t j = 0; j < found.size(); ++j) {
                        if (!found[j] && (groups[j] == Group::OUT || groups[j] == Group::BOTH))
                            callback(0.0L);
                    }
                };

                std::tie(pval, eff_size) = brunner_munzel(generate_a, generate_b);

                if (pval < config.family_wise_error_rate) {
                    size_t nb_idx = kept.rank1(row_i) - 1;
                    nb_pvals[nb_idx].first = pval;
                    nb_pvals[nb_idx].second = row_i;
                }

            } else if (config.test_type == "cmh") {
                // size_t n = 0;
                int64_t in_sum = 0;
                int64_t out_sum = 0;
                for (const auto& [j, c] : row) {
                    // n += c;
                    if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
                        out_sum += c;
                    }

                    if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
                        in_sum += c;
                    }
                }
                eff_size = static_cast<long double>(in_sum) / in_kmers
                        - static_cast<long double>(out_sum) / out_kmers;
                // // int64_t a = in_sum;
                // int64_t b = in_kmers - in_sum;
                // int64_t c = out_sum;
                // // int64_t d = out_kmers - out_sum;
                // // int64_t N = total_kmers;
                // if (b + c >= 25) {
                //     long double stat = static_cast<long double>(pow(b - c, 2.0)) / (b + c);
                //     pval = boost::math::cdf(
                //             boost::math::complement(boost::math::chi_squared(1), stat));
                // } else if (b >= c) {
                //     pval = boost::math::cdf(
                //             boost::math::complement(boost::math::binomial(b + c, 0.5), b - 1));
                // } else {
                //     throw std::runtime_error("CMH not implemented for this case");
                // }
                int64_t t = total_kmers;
                int64_t n1 = in_kmers;
                int64_t n2 = out_kmers;
                int64_t m1 = in_sum + out_sum;
                int64_t m2 = t - m1;
                assert(m2 >= 0);

                if (m2 == 0) {
                    common::logger->error(
                            "This row has all counts. Too few data points. Use fisher "
                            "instead");
                    throw std::domain_error("Test failed");
                }

                double lbase
                        = log(n1) + log(n2) + log(m1) + log(m2) - log(t) * 2 - log(t - 1);
                double factor = static_cast<double>(n1 * m1) / t;

                auto get_chi_stat = [&](int64_t a) {
                    int64_t b = n1 - a;
                    int64_t c = m1 - a;
                    int64_t d = m2 - b;
                    assert(d == n2 - c);

                    if (b < 0 || c < 0 || d < 0)
                        return 0.0;

                    double a_shift = static_cast<double>(a) - factor;

                    if (a_shift > 0)
                        return exp(log(a_shift) * 2.0 - lbase);

                    return a_shift * a_shift / exp(lbase);
                };

                double chi_stat = get_chi_stat(in_sum);
                if (chi_stat <= 0) {
                    common::logger->error(
                            "Test statistic {} <= 0, too few data points. Use fisher "
                            "instead.\t{}\t{},{},{},{}\t{}\t{}",
                            chi_stat, in_sum, n1, n2, m1, m2, t, factor);
                    throw std::domain_error("Test failed");
                }

                double max_chi_stat
                        = std::max(get_chi_stat(0), get_chi_stat(std::min(n1, m1)));
                if (max_chi_stat <= 0) {
                    common::logger->error(
                            "Best test statistic {} <= 0, too few data points. Use "
                            "fisher instead.\t{},{},{},{}\t{}\t{}",
                            max_chi_stat, n1, n2, m1, m2, t, factor);
                    throw std::domain_error("Test failed");
                }

                pval = boost::math::cdf(
                        boost::math::complement(boost::math::chi_squared(1), chi_stat));
                if (pval < config.family_wise_error_rate) {
                    size_t nb_idx = kept.rank1(row_i) - 1;
                    nb_pvals[nb_idx].first = pval;
                    nb_pvals[nb_idx].second = row_i;
                }
            } else if (config.test_type == "bb_hypergeometric") {
                size_t nb_idx = kept.rank1(row_i) - 1;
                size_t n = 0;
                int64_t in_sum = 0;
                int64_t out_sum = 0;
                for (const auto& [j, c] : row) {
                    n += c;
                    if (groups[j] == Group::OUT || groups[j] == Group::BOTH) {
                        out_sum += c;
                    }

                    if (groups[j] == Group::IN || groups[j] == Group::BOTH) {
                        in_sum += c;
                    }
                }

                // if (n < n_cutoff)
                //     return;

                if (n > 0) {
                    // long double min_pval = 1.0;
                    const auto& [r_a, r_b, bb_base] = bb_bases[n];
                    long double r_n = r_a + r_b;
                    long double p_a = 0.1;
                    long double p_b = 0.1;
                    long double p_n = 0.1;
                    long double min_pval = 1.0;
                    {
                        long double l21p = log1p(-p_n) / log2l(2.0);
                        long double base = bb_base
                                + (lgammal(n + 1) - lgammal(r_n + n)) / log2l(2.0)
                                - n * l21p;
                        auto get_pval = [&](int64_t s) {
                            int64_t t = n - s;
                            long double sbase = (lgammal(r_a + s) + lgammal(r_b + t)
                                                 - lgammal(s + 1) - lgammal(t + 1)
                                                 + s * log1pl(-p_a) + t * log1pl(-p_b))
                                    / logl(2.0);
                            return exp2l(base + sbase);
                        };

                        long double pval0 = get_pval(0);
                        long double pvaln = get_pval(n);
                        min_pval = std::min(pval0, pvaln) * (1 + (pval0 == pvaln));
                    }

                    if (min_pval < config.family_wise_error_rate) {
                        nb_pvals[nb_idx].second = row_i;
                        long double midpoint = r_a * n / (r_a + r_b);

                        auto get_deviance_a = [&](long double s) {
                            if (s == 0)
                                return r_a * logl(p_a) * -2.0;

                            return ((logl(s) - log1pl(-p_a)) * s - (r_a + s) * logl(r_a + s)
                                    + r_a * (logl(r_a) - logl(p_a)))
                                    * 2.0;
                        };

                        auto get_deviance_b = [&](long double t) {
                            if (t == 0)
                                return r_b * logl(p_b) * -2.0;

                            return ((logl(t) - log1pl(-p_b)) * t - (r_b + t) * logl(r_b + t)
                                    + r_b * (logl(r_b) - logl(p_b)))
                                    * 2.0;
                        };

                        auto get_deviance = [&](long double s, long double t) {
                            return get_deviance_a(s) + get_deviance_b(t);
                        };

                        long double min_dev
                                = get_deviance(midpoint,
                                               static_cast<long double>(n) - midpoint);
                        long double in_dev = get_deviance(in_sum, out_sum);

                        if (in_dev > min_dev) {
                            eff_size = static_cast<long double>(in_sum) - midpoint;
                            pval = 0.0;
                            long double base = bb_base
                                    + (lgammal(n + 1) - lgammal(r_n + n) - n * log1pl(-p_n))
                                            / logl(2.0);
                            double l1pa = log1pl(-p_a) / logl(2.0);
                            double l1pb = log1pl(-p_b) / logl(2.0);

                            {
                                size_t s = 0;
                                size_t t = n;
                                long double sbase = (lgammal(r_a + s) + lgammal(r_b + t)
                                                     - lgammal(s + 1) - lgammal(t + 1))
                                                / logl(2.0)
                                        + s * l1pa + t * l1pb;
                                pval += exp2(base + sbase);
                                for (++s; s <= n; ++s) {
                                    long double dev = get_deviance(s, t - 1);
                                    if (dev >= in_dev) {
                                        sbase += log2l(r_a + s) - log2l(r_b + t)
                                                - log2l(s + 1) + log2l(t + 1) + l1pa - l1pb;
                                        --t;
                                        pval += exp2(base + sbase);
                                        if (pval >= config.family_wise_error_rate)
                                            break;
                                    } else {
                                        break;
                                    }
                                }
                            }
                            if (pval < config.family_wise_error_rate
                                && get_deviance(n, 0) >= in_dev) {
                                size_t s = n;
                                size_t t = 0;
                                long double sbase = (lgammal(r_a + s) + lgammal(r_b + t)
                                                     - lgammal(s + 1) - lgammal(t + 1))
                                                / logl(2.0)
                                        + s * l1pa + t * l1pb;
                                pval += exp2(base + sbase);
                                for (++t; t <= n; ++t) {
                                    long double dev = get_deviance(s - 1, t);
                                    if (dev >= in_dev) {
                                        sbase -= log2l(r_a + s) - log2l(r_b + t)
                                                - log2l(s + 1) + log2l(t + 1) + l1pa - l1pb;
                                        --s;
                                        pval += exp2(base + sbase);
                                        if (pval >= config.family_wise_error_rate)
                                            break;
                                    } else {
                                        break;
                                    }
                                }
                            }
                            nb_pvals[nb_idx].first = pval;
                        }
                    } else {
                        nb_pvals[nb_idx].second = row_i;
                        nb_pvals[nb_idx].first = 1.1;
                    }
                }
            }

            set_pval(pval, eff_size);
        });

        if (fdr) {
            // if (false) {
            common::logger->trace("Correcting p-vals");
            nb_pvals.erase(std::remove_if(nb_pvals.begin(), nb_pvals.end(),
                                          [&](const auto& a) { return a.first > 1.0; }),
                           nb_pvals.end());
            size_t num_tests = nb_pvals.size();
            common::logger->trace("Discarded {} / {} untestable hypotheses",
                                  nelem - num_tests, nelem);
            nb_pvals.erase(std::remove_if(nb_pvals.begin(), nb_pvals.end(),
                                          [&](const auto& a) {
                                              return a.second == kept.size();
                                          }),
                           nb_pvals.end());
            common::logger->trace("Sorting {} / {} p-vals", nb_pvals.size(), num_tests);
            std::sort(nb_pvals.begin(), nb_pvals.end());

            // size_t k = 0;
            // for (size_t i = 0; i < nb_pvals.size(); ++i) {
            //     const auto& [pval, row_i] = nb_pvals[i];
            //     if (pval >= config.family_wise_error_rate)
            //         break;

            //     if (pval * (num_tests - i) >= config.family_wise_error_rate) {
            //         node_index node =
            //         AnnotatedDBG::anno_to_graph_index(nb_pvals[i].second);
            //         indicator_in[node] = false; indicator_out[node] = false;
            //     } else {
            //         ++k;
            //     }
            // }
            // common::logger->trace(
            //         "Found {} / {} significant p-values. Minimum p-value is {}", k,
            //         nb_pvals.size(), nb_pvals[0].first);

            long double harm = 0.0;
            for (size_t i = 1; i <= num_tests; ++i) {
                harm += 1.0 / i;
            }

            common::logger->trace("Selecting cut-off");
            size_t k = 0;
            for (size_t i = 0; i < nb_pvals.size(); ++i) {
                if (nb_pvals[i].first >= config.family_wise_error_rate)
                    break;

                if (nb_pvals[i].first
                    <= config.family_wise_error_rate * (i + 1) / num_tests / harm) {
                    k = std::max(i + 1, k);
                    node_index node = AnnotatedDBG::anno_to_graph_index(nb_pvals[i].second);
                    if (!indicator_in[node] && !indicator_out[node]) {
                        common::logger->error(
                                "Node {} not marked, but p-value {} significant", node,
                                nb_pvals[i].first);
                        throw std::runtime_error("Index fail");
                    }
                }
            }

            common::logger->trace(
                    "Found {} / {} significant p-values. Minimum p-value is {}", k,
                    nb_pvals.size(), nb_pvals[0].first);
            for (size_t i = k; i < nb_pvals.size(); ++i) {
                const auto& [pval, row_i] = nb_pvals[i];
                if (nb_pvals[i].first >= config.family_wise_error_rate)
                    break;

                node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
                if (!indicator_in[node] && !indicator_out[node]) {
                    common::logger->error(
                            "Node {} not marked, but original p-value {} significant",
                            node, nb_pvals[i].first);
                    throw std::runtime_error("Index fail");
                }
                indicator_in[node] = false;
                indicator_out[node] = false;
            }
        }
        // if (config.test_by_unitig) {
        //     common::logger->trace("Combining p-values for unitigs");
        //     std::mutex m_vec_mu;
        //     clean_masked_graph->call_unitigs([&](const std::string&, const auto &path) {
        //         long double pval_min = 1.0;
        //         {
        //             size_t row = AnnotatedDBG::graph_to_anno_index(path[0]);
        //             size_t row_rank = kept.rank1(row);
        //             size_t bucket_idx = std::lower_bound(offsets.begin(), offsets.end(), row_rank) - offsets.begin() - 1;
        //             size_t offset = row_rank - offsets[bucket_idx];
        //             std::lock_guard<std::mutex> lock(m_vec_mu);
        //             pval_min = get<PValStorage>(pvals_min_buckets[bucket_idx][offset]);
        //         }

        //         if (pval_min * k_min >= config.family_wise_error_rate)
        //             return;

        //         std::vector<long double> pvals;
        //         pvals.reserve(path.size());
        //         long double comb_eff_size = 0.0;
        //         {
        //             std::lock_guard<std::mutex> lock(m_vec_mu);
        //             for (node_index node : path) {
        //                 size_t row = AnnotatedDBG::graph_to_anno_index(node);
        //                 size_t row_rank = kept.rank1(row);
        //                 size_t bucket_idx = std::lower_bound(offsets.begin(), offsets.end(), row_rank) - offsets.begin() - 1;
        //                 size_t offset = row_rank - offsets[bucket_idx];
        //                 pvals.emplace_back(get<PValStorage>(pvals_buckets[bucket_idx][offset]));
        //                 comb_eff_size += get<PValStorage>(eff_size_buckets[bucket_idx][offset]);
        //             }
        //         }
        //         long double comp_pval = combine_pvals(pvals);
        //         if (comp_pval * k >= config.family_wise_error_rate || comb_eff_size == 0.0)
        //             return;

        //         auto &indicator = comb_eff_size > 0 ? indicator_in : indicator_out;
        //         for (node_index node : path) {
        //             set_bit(indicator.data(), node, true, std::memory_order_relaxed);
        //         }
        //     });

        //     pvals_min_buckets.resize(0);
        //     tmp_min_buckets.resize(0);
        //     pvals_buckets.resize(0);
        //     tmp_buckets.resize(0);
        //     eff_size_buckets.resize(0);
        //     tmp_eff_buckets.resize(0);
        // }
        std::atomic_thread_fence(std::memory_order_acquire);
    }

    // common::logger->trace("Computing minimum p-values");
    // std::vector<std::pair<double, size_t>> m;
    // std::vector<int64_t> m_sums;
    // {
    //     std::vector<std::vector<std::pair<double, size_t>>> ms(num_threads + 1);
    //     std::vector<std::vector<tsl::hopscotch_map<std::vector<int64_t>,
    //     std::pair<double, size_t>, utils::VectorHash>>> ms_vecs(num_threads + 1);
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
    //                 double out_stat = static_cast<double>(out_sum) / out_kmers *
    //                 in_kmers; if (out_stat != 0)
    //                     out_stat_int = out_stat > 0 ? ceil(out_stat) : floor(out_stat);

    //                 in_kmer = count_in >= config.min_in_recurrence && count_out <=
    //                 config.max_out_recurrence && in_sum > (out_kmers > 0 ? out_stat :
    //                 0.0); out_kmer = count_out >= config.min_out_recurrence && count_in
    //                 <= config.max_in_recurrence && in_sum < out_stat;
    //             }
    //         } else {
    //             return;
    //         }

    //         size_t n = config.test_type != "poisson_binom" ? in_sum + out_sum :
    //         row.size(); if (!in_kmer && !out_kmer) {
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
    //                 assert(ms[0][j].first == 1.1 || ms[i][j].first == 1.1 || ms[0][j]
    //                 == ms[i][j]); ms[0][j].first = std::min(ms[0][j].first,
    //                 ms[i][j].first); ms[0][j].second += ms[i][j].second;
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
    //     common::logger->trace("Picked: k: {}\tn: {}\tpval_min: {}", k, n_cutoff,
    //     config.family_wise_error_rate / k);
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
    //                 double out_stat = static_cast<double>(out_sum) / out_kmers *
    //                 in_kmers; in_kmer = count_in >= config.min_in_recurrence &&
    //                 count_out <= config.max_out_recurrence && in_sum > (out_kmers > 0 ?
    //                 out_stat : 0.0); out_kmer = count_out >= config.min_out_recurrence
    //                 && count_in <= config.max_in_recurrence && in_sum < out_stat;
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

    //             if (config.test_type != "gnb_exact" && config.test_type != "lnb_exact"
    //             && pval < 1.1 && m[n].first - pval > 1e-10) {
    //                 common::logger->error("Min p-val estimate too high: min {} > cur
    //                 {}\tn: {}\ttest: {}", m[n].first, pval, n, config.test_type); throw
    //                 std::runtime_error("Test failed");
    //             }

    //             if (config.test_type != "notest" && keep_all_pvals) {
    //                 auto &pvals = pvals_buckets[bucket_idx];
    //                 if constexpr(preallocated) {
    //                     pvals[node] = bit_cast<uint64_t>(pval);
    //                 } else {
    //                     pvals.push_back(bit_cast<uint64_t>(pval));
    //                 }
    //             }

    //             if (!config.test_by_unitig && in_kmer != out_kmer && pval * k <
    //             config.family_wise_error_rate) {
    //                 bool use_atomic = parallel && (node % 64 == 0);
    //                 set_bit((in_kmer ? indicator_in : indicator_out).data(), node,
    //                 use_atomic, MO_RELAXED);
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
            is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC);

    auto masked_graph_out = std::make_shared<MaskedDeBruijnGraph>(
            graph_ptr, std::make_unique<bitmap_vector>(std::move(indicator_out)), true,
            is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC);

    return std::make_tuple(masked_graph_in, masked_graph_out, std::move(pvals),
                           std::move(tmp_file));
}

template <class PValStorage>
std::tuple<std::shared_ptr<DeBruijnGraph>,
           std::shared_ptr<DeBruijnGraph>,
           PValStorage,
           std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                         const std::vector<std::string>& files,
                         const tsl::hopscotch_set<Label>& labels_in,
                         const tsl::hopscotch_set<Label>& labels_out,
                         const DifferentialAssemblyConfig& config,
                         size_t num_threads,
                         std::filesystem::path tmp_dir,
                         size_t num_parallel_files) {
    num_parallel_files = get_num_threads();
    common::logger->trace("Labels in: {}", fmt::join(labels_in, ","));
    common::logger->trace("Labels out: {}", fmt::join(labels_out, ","));

    size_t total_labels = 0;
    for (const auto& file : files) {
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

    return std::visit(
            [&](auto&& column_values_all) {
                using ColumnValues = typename std::decay<decltype(column_values_all)>::type;
                using ValuesContainerPtr = typename ColumnValues::value_type;
                using ValuesContainer =
                        typename std::decay<typename ValuesContainerPtr::element_type>::type;
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
                            [&](uint64_t offset, const Label& label,
                                std::unique_ptr<bit_vector>&& column,
                                ValuesContainer&& column_values) {
                                bool is_in = labels_in.count(label);
                                bool is_out = labels_out.count(label);
                                if (is_in) {
                                    groups[offset] = !is_out ? Group::IN : Group::BOTH;
                                } else {
                                    groups[offset] = is_out ? Group::OUT : Group::OTHER;
                                }
                                columns_all[offset].reset(column.release());
                                column_values_all[offset]
                                        = std::make_unique<ValuesContainer>(
                                                std::move(column_values));
                            },
                            num_parallel_files);
                } else {
                    annot::ColumnCompressed<>::merge_load(
                            files,
                            [&](uint64_t offset, const Label& label,
                                std::unique_ptr<bit_vector>&& column) {
                                bool is_in = labels_in.count(label);
                                bool is_out = labels_out.count(label);
                                if (is_in) {
                                    groups[offset] = !is_out ? Group::IN : Group::BOTH;
                                } else {
                                    groups[offset] = is_out ? Group::OUT : Group::OTHER;
                                }
                                columns_all[offset].reset(column.release());
                            },
                            num_parallel_files);
                }

                uint8_t max_width = 0;
                for (const auto& col_vals : column_values_all) {
                    max_width = std::max(max_width, col_vals->width());
                }

                auto generate_rows = [&](const auto& callback) {
                    utils::call_rows<std::unique_ptr<const bit_vector>,
                                     std::unique_ptr<const ValuesContainer>, PairContainer, false>(
                            columns_all, column_values_all, callback);
                };

                bool parallel = get_num_threads() > 1;

                return mask_nodes_by_label_dual<value_type, PValStorage>(
                        graph_ptr,
                        [&](const std::vector<size_t>& min_counts, sdsl::bit_vector* kept)
                                -> std::vector<VectorMap<uint64_t, size_t>> {
                            common::logger->trace("Calculating count histograms");
                            std::vector<std::vector<VectorMap<uint64_t, size_t>>> hists_map_p(
                                    num_parallel_files);
                            for (auto& hists_map : hists_map_p) {
                                hists_map.resize(groups.size());
                            }

                            std::atomic_thread_fence(std::memory_order_release);
                            generate_rows([&](uint64_t row_i, const auto& row,
                                              size_t thread_id) {
                                if (row.empty()) {
                                    if (kept)
                                        unset_bit(kept->data(), row_i, parallel,
                                                  std::memory_order_relaxed);

                                    return;
                                }

                                bool found = false;
                                Vector<uint64_t> counts(groups.size());
                                for (const auto& [j, raw_c] : row) {
                                    if (min_counts.empty() || raw_c >= min_counts[j]) {
                                        counts[j] = raw_c;
                                        found = true;
                                    }
                                }

                                if (!found) {
                                    if (kept)
                                        unset_bit(kept->data(), row_i, parallel,
                                                  std::memory_order_relaxed);

                                    return;
                                }

                                for (size_t j = 0; j < counts.size(); ++j) {
                                    ++hists_map_p[thread_id][j][counts[j]];
                                }
                            });
                            std::atomic_thread_fence(std::memory_order_acquire);

                            common::logger->trace("Merging histograms");
                            for (size_t thread_id = 1; thread_id < hists_map_p.size();
                                 ++thread_id) {
                                for (size_t j = 0; j < hists_map_p[thread_id].size(); ++j) {
                                    for (const auto& [k, c] : hists_map_p[thread_id][j]) {
                                        hists_map_p[0][j][k] += c;
                                    }
                                }
                            }

                            return hists_map_p[0];
                        },
                        generate_rows, groups, config, num_threads, tmp_dir,
                        num_parallel_files,
                        [&]() {
                            columns_all.clear();
                            column_values_all.clear();
                        },
                        max_width);
            },
            column_values_all);
}

template std::tuple<std::shared_ptr<DeBruijnGraph>,
                    std::shared_ptr<DeBruijnGraph>,
                    std::vector<uint64_t>,
                    std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<std::vector<uint64_t>>(std::shared_ptr<const DeBruijnGraph>,
                                                const std::vector<std::string>&,
                                                const tsl::hopscotch_set<Label>&,
                                                const tsl::hopscotch_set<Label>&,
                                                const DifferentialAssemblyConfig&,
                                                size_t,
                                                std::filesystem::path,
                                                size_t);
template std::tuple<std::shared_ptr<DeBruijnGraph>,
                    std::shared_ptr<DeBruijnGraph>,
                    sdsl::int_vector_buffer<64>,
                    std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<sdsl::int_vector_buffer<64>>(std::shared_ptr<const DeBruijnGraph>,
                                                      const std::vector<std::string>&,
                                                      const tsl::hopscotch_set<Label>&,
                                                      const tsl::hopscotch_set<Label>&,
                                                      const DifferentialAssemblyConfig&,
                                                      size_t,
                                                      std::filesystem::path,
                                                      size_t);

template <class PValStorage>
std::tuple<std::shared_ptr<DeBruijnGraph>,
           std::shared_ptr<DeBruijnGraph>,
           PValStorage,
           std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(
        const AnnotatedDBG& anno_graph,
        const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label>& labels_in,
        const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label>& labels_out,
        const DifferentialAssemblyConfig& config,
        size_t num_threads,
        std::filesystem::path tmp_dir) {
    num_threads = get_num_threads();
    common::logger->trace("Labels in: {}", fmt::join(labels_in, ","));
    common::logger->trace("Labels out: {}", fmt::join(labels_out, ","));

    size_t total_labels = labels_in.size() + labels_out.size();
    const auto& annotation = anno_graph.get_annotator();

    std::vector<Group> groups;
    groups.reserve(total_labels);
    for (const auto& label : annotation.get_label_encoder().get_labels()) {
        bool is_in = labels_in.count(label);
        bool is_out = labels_out.count(label);
        if (is_in) {
            groups.emplace_back(!is_out ? Group::IN : Group::BOTH);
        } else {
            groups.emplace_back(is_out ? Group::OUT : Group::OTHER);
        }
    }

    auto graph_ptr
            = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());

    const auto& matrix = annotation.get_matrix();
    using Row = annot::matrix::BinaryMatrix::Row;

    if (const auto* int_matrix = dynamic_cast<const annot::matrix::IntMatrix*>(&matrix)) {
        using value_type = annot::matrix::IntMatrix::Value;
        return mask_nodes_by_label_dual<value_type, PValStorage>(
                graph_ptr,
                [&](const std::vector<size_t>& min_counts, sdsl::bit_vector* unmark_discarded) {
                    return int_matrix->get_histograms(min_counts, unmark_discarded);
                },
                [&](const auto& callback) {
                    int_matrix->call_row_values(callback, false /* ordered */);
                },
                groups, config, num_threads, tmp_dir, num_threads);
    } else {
        auto generate_bit_rows = [&](const auto& callback) {
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
                for (size_t row_batch_i = begin; row_batch_i < end;
                     row_batch_i += rows_per_update) {
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
                [&](const std::vector<size_t>&,
                    sdsl::bit_vector* kept) -> std::vector<VectorMap<uint64_t, size_t>> {
                    common::logger->trace("Calculating count histograms");
                    std::vector<std::vector<VectorMap<uint64_t, size_t>>> hists_map_p(
                            num_threads);
                    for (auto& hists_map : hists_map_p) {
                        hists_map.resize(groups.size());
                    }
                    std::atomic_thread_fence(std::memory_order_release);
                    generate_bit_rows(
                            [&](uint64_t row_i, const auto& set_bits, size_t thread_id) {
                                if (set_bits.empty()) {
                                    if (kept)
                                        unset_bit(kept->data(), row_i, parallel,
                                                  std::memory_order_relaxed);

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
                            for (const auto& [k, c] : hists_map_p[thread_id][j]) {
                                hists_map_p[0][j][k] += c;
                            }
                        }
                    }

                    return hists_map_p[0];
                },
                [&](const auto& callback) {
                    generate_bit_rows(
                            [&](uint64_t row_i, const auto& set_bits, size_t thread_id) {
                                Vector<std::pair<uint64_t, uint64_t>> row;
                                row.reserve(set_bits.size());
                                for (auto j : set_bits) {
                                    row.emplace_back(j, 1);
                                }

                                callback(row_i, row, thread_id);
                            });
                },
                groups, config, num_threads, tmp_dir, num_threads);
    }
}

template std::tuple<std::shared_ptr<DeBruijnGraph>,
                    std::shared_ptr<DeBruijnGraph>,
                    std::vector<uint64_t>,
                    std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<std::vector<uint64_t>>(
        const AnnotatedDBG&,
        const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label>&,
        const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label>&,
        const DifferentialAssemblyConfig&,
        size_t,
        std::filesystem::path);
template std::tuple<std::shared_ptr<DeBruijnGraph>,
                    std::shared_ptr<DeBruijnGraph>,
                    sdsl::int_vector_buffer<64>,
                    std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual<sdsl::int_vector_buffer<64>>(
        const AnnotatedDBG&,
        const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label>&,
        const tsl::hopscotch_set<typename AnnotatedDBG::Annotator::Label>&,
        const DifferentialAssemblyConfig&,
        size_t,
        std::filesystem::path);

} // namespace graph
} // namespace mtg
