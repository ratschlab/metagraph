#include "annotated_graph_algorithm.hpp"

#include "common/logger.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vectors/bitmap.hpp"
#include "graph/representation/masked_graph.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"

// https://github.com/tlemane/kmdiff/blob/6a56ce6f20abbf63928a19ebbfecb1968efd4cd3/include/kmdiff/model.hpp
/*****************************************************************************
 *   kmdiff
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>

#include <alglibinternal.h>
#include <spdlog/spdlog.h>
#include <specialfunctions.h>

#define KMTRICKS_PUBLIC
#include <kmtricks/utils.hpp>

#include <kmdiff/range.hpp>
#include <kmdiff/kmer.hpp>
#include <kmdiff/log_factorial_table.hpp>
#include <kmdiff/correction.hpp>
#include <kmdiff/imodel.hpp>


// Myrthe: the following code calculates the p-value for a single k-mer across the samples. All comes together in `process`.
namespace kmdiff {

    template <size_t MAX_C>
    class Model
    {
        using count_t = typename km::selectC<MAX_C>::type;

    public:
        Model() {}

        template <typename Iterable>
        double compute_mean(Iterable&& v)
        {
            return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        }

        template <typename Iterable>
        std::tuple<double, size_t> compute_mean_e(Iterable& v)
        {
            size_t sum = 0;
            size_t positive = 0;
            for (auto& e : v)
            {
                if (e > 0) positive++;
                sum += e;
            }
            double mean = sum / static_cast<double>(v.size());
            return std::make_tuple(mean, positive);
        }

        template <typename Iterable>
        std::tuple<double, size_t> compute_sum_e(Iterable& v)
        {
            double sum = 0;
            size_t positive = 0;
            for (auto& e : v)
            {
                if (e > 0) positive++;
                sum += e;
            }
            return std::make_tuple(sum, positive);
        }

        template <typename Iterable>
        std::tuple<double, double> compute_mean_sd(Iterable& v)
        {
            double mean = compute_mean(v);
            double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
            double sd = std::sqrt(sq_sum / v.size() - mean * mean);
            return std::make_tuple(mean, sd);
        }

    };

    template <size_t MAX_C>
    class PoissonLikelihood : public IModel<MAX_C>, public Model<MAX_C>
    {
        /*
          Based on https://github.com/atifrahman/HAWK
          https://doi.org/10.7554/eLife.32920.001
          https://doi.org/10.1371/journal.pone.0245058
        */

        using count_t = typename km::selectC<MAX_C>::type;

    public:
        PoissonLikelihood(
                size_t nb_controls,
                size_t nb_cases,
                const std::vector<size_t>& total_controls,
                const std::vector<size_t>& total_cases,
                size_t preload)
                : m_nb_controls(nb_controls),
                  m_nb_cases(nb_cases),
                  m_total_kmer_controls(total_controls),
                  m_total_kmer_cases(total_cases),
                  m_preload(preload)
        {
        }

        void configure(const std::string& config) override {}
    private:
        double log_factorial(int k)
        {
            double res = 0;
            while (k > 1)
            {
                res += log(k);
                k--;
            }
            return res;
        }

        double poisson_prob(int k, double lambda)
        {
            if (lambda <= 0) return 0;
            if (k < 0) k = 0;
            return (-lambda + (k * log(lambda) - m_lf_table[k]));
        }

    public:

        std::tuple<double, Significance, double, double>
        process(const Range<count_t>& controls, const Range<count_t>& cases) override
        {
            auto [mean_control, positive_controls] = this->compute_sum_e(controls);
            auto [mean_case, positive_cases] = this->compute_sum_e(cases);

            double mean = (mean_control + mean_case) / static_cast<double>(m_sum_controls + m_sum_cases);

            double null_hypothesis = 0;
            double alt_hypothesis = 0;

            alt_hypothesis += poisson_prob(mean_control, mean_control);
            alt_hypothesis += poisson_prob(mean_case, mean_case);

            null_hypothesis += poisson_prob(mean_control, mean * m_sum_controls);
            null_hypothesis += poisson_prob(mean_case, mean * m_sum_cases);

            double likelihood_ratio = alt_hypothesis - null_hypothesis;

            if (likelihood_ratio < 0) likelihood_ratio = 0;
            double p_value = alglib::chisquarecdistribution(1, 2 * likelihood_ratio);

            Significance sign;

            mean_control = mean_control * m_sum_cases / m_sum_controls;

            if (mean_control < mean_case)
                sign = Significance::CASE;
            else if (mean_control > mean_case)
                sign = Significance::CONTROL;
            else
                sign = Significance::NO;

            return std::make_tuple(p_value, sign, mean_control, mean_case);
        }

        std::tuple<double, Significance, double, double>
        process_myrthe(double in_sum, double out_sum, double in_total_kmers, double out_total_kmers) override // previous is equivalent to mean_control,  mean_case, latter is equivalent to m_sum_controls, m_sum_cases
        {
            auto mean_control = out_sum; auto mean_case = in_sum;
            double mean = (mean_control + mean_case) / static_cast<double>(m_sum_controls + m_sum_cases);

            double null_hypothesis = 0;
            double alt_hypothesis = 0;

            alt_hypothesis += poisson_prob(mean_control, mean_control);
            alt_hypothesis += poisson_prob(mean_case, mean_case);

            null_hypothesis += poisson_prob(mean_control, mean * m_sum_controls);
            null_hypothesis += poisson_prob(mean_case, mean * m_sum_cases);

            double likelihood_ratio = alt_hypothesis - null_hypothesis;

            if (likelihood_ratio < 0) likelihood_ratio = 0;
            double p_value = alglib::chisquarecdistribution(1, 2 * likelihood_ratio);

            Significance sign;

            mean_control = mean_control * m_sum_cases / m_sum_controls;

            if (mean_control < mean_case)
                sign = Significance::CASE;
            else if (mean_control > mean_case)
                sign = Significance::CONTROL;
            else
                sign = Significance::NO;

            return std::make_tuple(p_value, sign, mean_control, mean_case);
        }

    private:
        size_t m_nb_controls{0};
        size_t m_nb_cases{0};

        const std::vector<size_t>& m_total_kmer_controls;
        const std::vector<size_t>& m_total_kmer_cases;

        size_t m_sum_controls{
                std::accumulate(m_total_kmer_controls.begin(), m_total_kmer_controls.end(), 0ULL)};

        size_t m_sum_cases{std::accumulate(m_total_kmer_cases.begin(), m_total_kmer_cases.end(), 0ULL)};

        size_t m_preload;
        LogFactorialTable m_lf_table{m_preload};
    };

} // end of namespace kmdiff
