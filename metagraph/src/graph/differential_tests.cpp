#include <cassert>
#include <math.h>
#include "differential_tests.hpp"

// https://github.com/tlemane/kmdiff/blob/6a56ce6f20abbf63928a19ebbfecb1968efd4cd3/include/kmdiff/log_factorial_table.hpp
// cpp https://github.com/tlemane/kmdiff/blob/6a56ce6f20abbf63928a19ebbfecb1968efd4cd3/src/log_factorial_table.cpp
LogFactorialTable::LogFactorialTable(size_t size)
        : m_size(size)
{
    m_table.reserve(m_size);
    for (size_t i=0; i<size; i++)
        m_table.push_back(log_factorial(i));
}

double LogFactorialTable::log_factorial(size_t k)
{
    double res = 0;
    while (k > 1)
    {
        res += log(k);
        k--;
    }
    return res;
}


DifferentialTest::DifferentialTest(double family_wise_error_rate, size_t total_hypotheses, size_t preload_table_size,
                                   size_t in_total_kmers, size_t out_total_kmers) :
        preload_table_size(preload_table_size), // similar to the width of the counts vector, this should be the maximum joint coverage over the in_labels resp. out_labels
        family_wise_error_rate(family_wise_error_rate),
        total_hypotheses(total_hypotheses),
        out_total_kmers(out_total_kmers),
        in_total_kmers(in_total_kmers)
{}


//double DifferentialTest::chi_sqrd_1_cdf(double x)// DOF = 1
//{
//    assert(x > 0); // assert that it is positive.
//    double igf = lower_incomplete_gamma_function(x/0.5);
//    if(std::isnan(igf) or std::isinf(igf) or igf <= 1e-8) return 1e-14;
//    return 1.0 - igf / gamma;
//}
//
//double DifferentialTest::lower_incomplete_gamma_function(double x) // lower incomplete_gamma_function lower incomplete. with s = 1/2
//{ // https://en.wikipedia.org/wiki/Incomplete_gamma_function#Lower_incomplete_gamma_function
//    double Sc = 2 * std::sqrt(x) * exp(-x);
//    double s = 0.5; double sum = 1.0; double nominator = 1.0; double denominator = 1.0;
//    for(int I = 0; I < 200; I++)
//    {
//        s++;
//        nominator *= x;
//        denominator *= s;
//        sum += (nominator / denominator);
//    }
//    return sum * Sc;
//}

double DifferentialTest::poisson_prob(int k, double lambda) // https://en.wikipedia.org/wiki/Poisson_distribution
{
    if (lambda <= 0) return 0;
    if (k < 0) k = 0;
    return (-lambda + (k * std::log(lambda) - lf_table[k]));
}

//m be the total number of null hypotheses
// family-wise error rate (FWER)
bool DifferentialTest::bonferroni_correction(double& pvalue) //  &family_wise_error_rate, std::size_t & total_hypotheses, --> can also be class parameters.
{
    return pvalue < (family_wise_error_rate / total_hypotheses);
}

// adapted code from kmdiff.
// previous is equivalent to mean_control,  mean_case, latter is equivalent to m_sum_controls, m_sum_cases
std::tuple<double, bool> DifferentialTest::likelihood_ratio_test(double in_sum, double out_sum)
{
    double mean = (out_sum + in_sum) / static_cast<double>(out_total_kmers + in_total_kmers);

    double null_hypothesis = 0;
    double alt_hypothesis = 0;

    alt_hypothesis += poisson_prob(out_sum, out_sum);
    alt_hypothesis += poisson_prob(in_sum, in_sum);

    null_hypothesis += poisson_prob(out_sum, mean * out_total_kmers);
    null_hypothesis += poisson_prob(in_sum, mean * in_total_kmers);

    double likelihood_ratio = alt_hypothesis - null_hypothesis;

    if (likelihood_ratio < 0) likelihood_ratio = 0;
    double p_value = Chi2PLookupTable.getPValue(2 * likelihood_ratio, 1);

    auto out_sum_normalized = out_sum * in_total_kmers / out_total_kmers;

    bool sign = false;
    if (out_sum_normalized < in_sum) sign = true;

    return std::make_tuple(p_value, sign);
}

