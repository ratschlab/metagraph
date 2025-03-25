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
#include "graph/representation/hash/dbg_sshash.hpp"
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

long double combine_pvals(const std::vector<long double>& pvals) {
    static const boost::math::cauchy dist;
    long double stat = 0.0;
    for (long double pval : pvals) {
        stat += tan((0.5 - pval) * M_PI);
    }
    stat /= pvals.size();
    return boost::math::cdf(boost::math::complement(dist, stat));
};


struct PairVectorHash {
    template <class Vector>
    std::size_t operator()(const Vector& vector) const {
        uint64_t hash = 0;
        for (const auto& [a, b] : vector) {
            hash ^= a + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            hash ^= b + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return static_cast<std::size_t>(hash);
    }
};

template <typename value_type, class PValStorage, typename PairContainer, typename HistGetter, typename Generator, typename UnitigGenerator>
std::tuple<std::shared_ptr<DeBruijnGraph>,
           std::shared_ptr<DeBruijnGraph>,
           PValStorage,
           std::unique_ptr<utils::TempFile>>
mask_nodes_by_label_dual(
        std::shared_ptr<const DeBruijnGraph> graph_ptr,
        const HistGetter& get_hist_map,
        const Generator& generate_rows,
        const UnitigGenerator& generate_unitigs,
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
    bit_vector_stat kept(kept_bv);

    clean_masked_graph = std::make_unique<MaskedDeBruijnGraph>(
            graph_ptr,
            [&](node_index node) {
                return node != DeBruijnGraph::npos
                        && kept[AnnotatedDBG::graph_to_anno_index(node)];
            },
            true, is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC);

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
            }
        });
    };

    size_t nelem = kept.num_set_bits();

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
    clean_masked_graph = std::make_unique<MaskedDeBruijnGraph>(
            graph_ptr,
            [&](node_index node) {
                return node != DeBruijnGraph::npos
                        && kept[AnnotatedDBG::graph_to_anno_index(node)];
            },
            true, is_primary ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC);
    if (config.test_by_unitig && config.test_type != "notest") {
        std::vector<std::tuple<long double, size_t, Group, node_index>> nb_pvals;
        nb_pvals.resize(nelem,
                        std::make_tuple(1.0L, graph_ptr->max_index() + 1, Group::OTHER,
                                        graph_ptr->max_index() + 1));

        generate_clean_rows([&](auto row_i, const auto& row, size_t) {
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

            size_t idx = kept.rank1(row_i) - 1;
            auto& [pval, monotig_id, group, node] = nb_pvals[idx];

            long double eff_size;

            std::tie(pval, eff_size) = mann_whitneyu(generate_a, generate_b);
            if (eff_size != 0.0) {
                node = AnnotatedDBG::anno_to_graph_index(row_i);
                if (eff_size > 0) {
                    group = Group::IN;
                } else if (eff_size < 0) {
                    group = Group::OUT;
                }
            }
        });

        uint64_t num_tests = 0;
        std::atomic_thread_fence(std::memory_order_release);
        clean_masked_graph->call_unitigs(
                [&](const std::string&, const auto& path) {
                    std::vector<long double> pvals;
                    pvals.reserve(path.size());

                    Group last_group = Group::BOTH;
                    size_t last_i = 0;

                    std::vector<size_t> ranks;
                    ranks.reserve(path.size());
                    for (node_index node : path) {
                        size_t row_i = AnnotatedDBG::graph_to_anno_index(node);
                        ranks.emplace_back(kept.rank1(row_i) - 1);
                    }

                    for (size_t i = 0; i < path.size(); ++i) {
                        size_t idx = ranks[i];
                        const auto& [pval, monotig_id, group, stored_node] = nb_pvals[idx];
                        if (group != last_group) {
                            if (pvals.size()) {
                                uint64_t monotig_id
                                        = __atomic_fetch_add(&num_tests, 1,
                                                             std::memory_order_relaxed);
                                long double comb_pval = combine_pvals(pvals);
                                for (size_t j = last_i; j < i; ++j) {
                                    size_t jdx = ranks[j];
                                    std::get<0>(nb_pvals[jdx]) = comb_pval;
                                    std::get<1>(nb_pvals[jdx]) = monotig_id;
                                    if (comb_pval < config.family_wise_error_rate) {
                                        if (last_group == Group::IN) {
                                            set_bit(indicator_in.data(), path[j], true,
                                                    std::memory_order_relaxed);
                                        } else if (last_group == Group::OUT) {
                                            set_bit(indicator_out.data(), path[j], true,
                                                    std::memory_order_relaxed);
                                        }
                                    }
                                }
                            }
                            pvals.clear();
                            last_group = group;
                            last_i = i;
                        }
                        pvals.emplace_back(pval);
                    }
                    if (pvals.size()) {
                        uint64_t monotig_id = __atomic_fetch_add(&num_tests, 1,
                                                                 std::memory_order_relaxed);
                        long double comb_pval = combine_pvals(pvals);
                        for (size_t j = last_i; j < path.size(); ++j) {
                            size_t jdx = ranks[j];
                            std::get<0>(nb_pvals[jdx]) = comb_pval;
                            std::get<1>(nb_pvals[jdx]) = monotig_id;
                            if (comb_pval < config.family_wise_error_rate) {
                                if (last_group == Group::IN) {
                                    set_bit(indicator_in.data(), path[j], true,
                                            std::memory_order_relaxed);
                                } else if (last_group == Group::OUT) {
                                    set_bit(indicator_out.data(), path[j], true,
                                            std::memory_order_relaxed);
                                }
                            }
                        }
                    }
                },
                num_threads);
        std::atomic_thread_fence(std::memory_order_acquire);
        common::logger->trace("Found {} monotigs", num_tests);
        std::atomic_thread_fence(std::memory_order_acquire);
        std::sort(nb_pvals.begin(), nb_pvals.end());

        long double harm = 0.0;
        for (size_t i = 1; i <= num_tests; ++i) {
            harm += 1.0 / i;
        }

        common::logger->trace("Selecting cut-off");
        size_t k = 0;
        uint64_t last_unitig_id = graph_ptr->max_index() + 1;
        size_t cur_count = 0;
        for (size_t i = 0; i < nb_pvals.size(); ++i) {
            const auto& [pval, unitig_id, group, node] = nb_pvals[i];
            if (pval >= config.family_wise_error_rate)
                break;

            if (unitig_id != last_unitig_id) {
                last_unitig_id = unitig_id;
                ++cur_count;
            }

            if (pval <= config.family_wise_error_rate * cur_count / num_tests / harm) {
                k = std::max(i + 1, k);
            }
        }

        common::logger->trace("Found {} / {} significant p-values. Minimum p-value is {}",
                              k, nb_pvals.size(), std::get<0>(nb_pvals[0]));

        for (size_t i = k; i < nb_pvals.size(); ++i) {
            const auto& [pval, unitig_id, group, node] = nb_pvals[i];
            if (pval >= config.family_wise_error_rate)
                break;

            if (!indicator_in[node] && !indicator_out[node]) {
                common::logger->error(
                        "Node {} not marked, but original p-value {} significant", node,
                        pval);
                throw std::runtime_error("Index fail");
            }
            indicator_in[node] = false;
            indicator_out[node] = false;
        }
    } else {
        // precompute negative binomial fits for poisson_bayes and nbinom_exact tests
        std::vector<std::pair<long double, long double>> nb_params(groups.size());
        std::pair<long double, long double> nb_params_a;
        std::pair<long double, long double> nb_params_b;
        std::pair<long double, long double> nb_params_null;
        long double nb_base = 0.0;
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

        if (config.test_type == "poisson_exact") {
            get_lambda = [&](size_t, const RowGenerator&) {
                return std::make_pair(static_cast<long double>(in_kmers) / nelem,
                                      static_cast<long double>(out_kmers) / nelem);
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

        bool fdr = config.test_by_unitig || config.test_type == "nbinom_exact"
                || config.test_type == "mwu" || config.test_type == "cmh";

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
                }
            });

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
                    }
                }

                if (config.test_type == "poisson_binom")
                    std::sort(m.rbegin(), m.rend());
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
        }

        // determine cutoffs for multiple testing correction
        auto [k_min, k, n_cutoff] = correct_pvals(m);
        if (fdr) {
            k_min = 1;
            k = 1;
            n_cutoff = 1;
        }

        common::logger->trace("Picked: k: {}\tn: {} / {}\tpval_min: {}\tn_cutoff: {}",
                              k_min, k, nelem, config.family_wise_error_rate / k, n_cutoff);

        common::logger->trace("Running differential tests");
        std::vector<std::pair<long double, size_t>> nb_pvals;
        std::vector<long double> eff_sizes;
        if (fdr) {
            nb_pvals.resize(kept.num_set_bits(), std::make_pair(1.0, kept.size()));
            eff_sizes.resize(kept.num_set_bits());
        }

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
                if (fdr) {
                    size_t nb_idx = kept.rank1(row_i) - 1;
                    nb_pvals[nb_idx].first = pval;
                    nb_pvals[nb_idx].second = row_i;
                    eff_sizes[nb_idx] = eff_size;
                }
                // if (!config.test_by_unitig) {
                if (pval * k < config.family_wise_error_rate) {
                    if (eff_size == 0) {
                        common::logger->error("Effect size 0 when p-value is {}", pval);
                        throw std::runtime_error("Pval failure");
                    }

                    if (!config.test_by_unitig) {
                        node_index node = AnnotatedDBG::anno_to_graph_index(row_i);
                        if (eff_size > 0) {
                            set_bit(indicator_in.data(), node, true,
                                    std::memory_order_relaxed);
                        } else if (eff_size < 0) {
                            set_bit(indicator_out.data(), node, true,
                                    std::memory_order_relaxed);
                        }
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

                if (!fdr) {
                    size_t front = n - std::min(n, num_labels_out);
                    double min_pval = std::min(pb_pvals[n][front], pb_pvals[n].back());

                    if (min_pval * k_min >= config.family_wise_error_rate)
                        return;
                }

                pval = pb_pvals[n][s];
                eff_size = static_cast<double>(s) - mid_points[n];
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

            } else if (config.test_type == "nbinom_exact") {
                // size_t nb_idx = kept.rank1(row_i) - 1;
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

                    // if (min_pval < config.family_wise_error_rate) {
                    // nb_pvals[nb_idx].second = row_i;
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
                        // nb_pvals[nb_idx].first = pval;
                    }
                    // } else {
                    //     nb_pvals[nb_idx].second = row_i;
                    //     nb_pvals[nb_idx].first = 1.1;
                    // }
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

                // if (pval < config.family_wise_error_rate) {
                //     size_t nb_idx = kept.rank1(row_i) - 1;
                //     nb_pvals[nb_idx].first = pval;
                //     nb_pvals[nb_idx].second = row_i;
                // }
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
                // if (pval < config.family_wise_error_rate) {
                //     size_t nb_idx = kept.rank1(row_i) - 1;
                //     nb_pvals[nb_idx].first = pval;
                //     nb_pvals[nb_idx].second = row_i;
                // }
            }

            set_pval(pval, eff_size);
        });

        std::atomic_thread_fence(std::memory_order_acquire);

        if (fdr) {
            size_t num_tests = nb_pvals.size();
            std::vector<std::pair<long double, size_t>> sorted_unitig_counts;
            if (config.test_by_unitig) {
                VectorMap<long double, size_t> pval_to_count;
                std::mutex mu;
                common::logger->trace("Combining p-values within unitigs");
                num_tests = 0;
                std::atomic_thread_fence(std::memory_order_release);
                clean_masked_graph->call_unitigs(
                        [&](const std::string&, const auto& path) {
                            __atomic_fetch_add(&num_tests, 1, std::memory_order_relaxed);
                            std::vector<long double> pvals;
                            pvals.reserve(path.size());
                            long double agg_eff_size = 0.0;
                            for (node_index node : path) {
                                auto row_i = AnnotatedDBG::graph_to_anno_index(node);
                                size_t nb_idx = kept.rank1(row_i) - 1;
                                long double pval = nb_pvals[nb_idx].first;
                                long double eff_size = eff_sizes[nb_idx];
                                pvals.emplace_back(pval);
                                agg_eff_size += eff_size;
                            }
                            long double agg_pval = combine_pvals(pvals);
                            for (node_index node : path) {
                                auto row_i = AnnotatedDBG::graph_to_anno_index(node);
                                size_t nb_idx = kept.rank1(row_i) - 1;
                                nb_pvals[nb_idx].first = agg_pval;
                                if (agg_pval < config.family_wise_error_rate) {
                                    if (agg_eff_size > 0) {
                                        set_bit(indicator_in.data(), node, true,
                                                std::memory_order_relaxed);
                                    } else {
                                        set_bit(indicator_out.data(), node, true,
                                                std::memory_order_relaxed);
                                    }
                                }
                            }
                            std::lock_guard<std::mutex> lock(mu);
                            pval_to_count[agg_pval] += path.size();
                        },
                        num_threads);
                std::atomic_thread_fence(std::memory_order_acquire);
                sorted_unitig_counts
                        = const_cast<std::vector<std::pair<long double, size_t>>&&>(
                                pval_to_count.values_container());
                std::sort(sorted_unitig_counts.begin(), sorted_unitig_counts.end());
                for (size_t i = 1; i < sorted_unitig_counts.size(); ++i) {
                    sorted_unitig_counts[i].second += sorted_unitig_counts[i - 1].second;
                }
            }

            // if (false) {
            common::logger->trace("Correcting p-vals");

            common::logger->trace("Discarded {} / {} untestable hypotheses",
                                  nelem - num_tests, nelem);
            nb_pvals.erase(std::remove_if(nb_pvals.begin(), nb_pvals.end(),
                                          [&](const auto& a) {
                                              return a.second == kept.size();
                                          }),
                           nb_pvals.end());
            common::logger->trace("Sorting {} / {} p-vals", nb_pvals.size(), num_tests);
            std::sort(nb_pvals.begin(), nb_pvals.end());

            long double harm = 0.0;
            for (size_t i = 1; i <= num_tests; ++i) {
                harm += 1.0 / i;
            }

            common::logger->trace("Selecting cut-off");
            size_t k = 0;
            auto it = sorted_unitig_counts.begin();
            for (size_t i = 0; i < nb_pvals.size(); ++i) {
                if (nb_pvals[i].first >= config.family_wise_error_rate)
                    break;

                if (it != sorted_unitig_counts.end() && i == it->second)
                    ++it;

                size_t cur_count = it == sorted_unitig_counts.end()
                        ? i
                        : it - sorted_unitig_counts.begin();
                if (nb_pvals[i].first <= config.family_wise_error_rate * (cur_count + 1)
                            / num_tests / harm) {
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
    }

    std::unique_ptr<utils::TempFile> tmp_file;
    PValStorage pvals;

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

                auto generate_unitigs = [&](const DeBruijnGraph&, const auto&) {
                    throw std::runtime_error(
                            "Unitigs not implemented for column annotations");
                };

                bool parallel = get_num_threads() > 1;

                return mask_nodes_by_label_dual<value_type, PValStorage, PairContainer>(
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
                        generate_rows, generate_unitigs, groups, config, num_threads,
                        tmp_dir, num_parallel_files,
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
        using PairContainer = Vector<std::pair<uint64_t, value_type>>;
        // std::mutex mu;

        auto generate_unitigs = [&](const DeBruijnGraph& graph, const auto& callback) {
            graph.call_unitigs(
                    [&](const std::string&, const auto& path) {
                        std::vector<annot::matrix::BinaryMatrix::Row> rows;
                        rows.reserve(path.size());
                        for (auto node : path) {
                            rows.emplace_back(AnnotatedDBG::graph_to_anno_index(node));
                        }
                        auto row_values = int_matrix->get_row_values(rows);
                        callback(path, row_values);
                    },
                    num_threads);
        };

        return mask_nodes_by_label_dual<value_type, PValStorage, PairContainer>(
                graph_ptr,
                [&](const std::vector<size_t>& min_counts, sdsl::bit_vector* unmark_discarded) {
                    return int_matrix->get_histograms(min_counts, unmark_discarded);
                },
                [&](const auto& callback) {
                    int_matrix->call_row_values(callback, false /* ordered */);
                },
                generate_unitigs, groups, config, num_threads, tmp_dir, num_threads);
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

        using value_type = annot::matrix::IntMatrix::Value;
        using PairContainer = Vector<std::pair<uint64_t, value_type>>;
        return mask_nodes_by_label_dual<uint64_t, PValStorage, PairContainer>(
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
                [&](const DeBruijnGraph&, const auto&) {
                    throw std::runtime_error(
                            "Unitigs not implemented for binary annotations");
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
