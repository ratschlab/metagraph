#include <math.h>
#include <string>
#include <numbers>
#include <iostream>
#include <cassert>

#include "differential_tests.hpp"
#include "lookup_table_chisqrd_cdf.cpp"


// https://github.com/tlemane/kmdiff/blob/6a56ce6f20abbf63928a19ebbfecb1968efd4cd3/include/kmdiff/log_factorial_table.hpp
// cpp https://github.com/tlemane/kmdiff/blob/6a56ce6f20abbf63928a19ebbfecb1968efd4cd3/src/log_factorial_table.cpp
namespace mtg{
//LogFactorialTable::LogFactorialTable(size_t size)
//        : m_size(size)
//{
//    m_table.reserve(m_size);
//    for (size_t i=0; i<size; i++)
//        m_table.push_back(log_factorial(i));
//}
//
//double LogFactorialTable::log_factorial(size_t k) // TODO this can be calculated more efficiently.
//{
//    double res = 0;
//    while (k > 1)
//    {
//        res += log(k);
//        k--;
//    }
//    return res;
//}
//
//double LogFactorialTable::approximate_log_factorial(size_t k) {
//    return k * log(k) - k * log2e_v;
//}


//double DifferentialTest::chi_sqrd_1_cdf(double x)// DOF = 1 // https://www.codeproject.com/Articles/432194/How-to-Calculate-the-Chi-Squared-P-Value
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

//double DifferentialTest::poisson_pmf(int k, double lambda) // https://en.wikipedia.org/wiki/Poisson_distribution and https://elifesciences.org/articles/32920#appendix-1
//{
//    if (lambda <= 0) return 0;
//    if (k < 0) k = 0;
//    return (-lambda + (k * std::log(lambda)) - lf_table[k]); // One does not need to subtract the factorial of k because it divides away in the likelihood ratio.
//}
//
//double DifferentialTest::binomial_pmf(int k, double lambda) // https://en.wikipedia.org/wiki/Poisson_distribution and https://elifesciences.org/articles/32920#appendix-1
//{
//    if (lambda <= 0) return 0;
//    if (k < 0) k = 0;
//    return (-lambda + (k * std::log(lambda)) - lf_table[k]); // One does not need to subtract the factorial of k because it divides away in the likelihood ratio.
//}
//


DifferentialTest::DifferentialTest(double family_wise_error_rate, size_t total_hypotheses, size_t ,
                                   size_t in_total_kmers, size_t out_total_kmers) : //preload_table_size(preload_table_size), // similar to the width of the counts vector, this should be the maximum joint coverage over the in_labels resp. out_labels
      family_wise_error_rate(family_wise_error_rate),
      total_hypotheses(total_hypotheses),
      out_total_kmers(out_total_kmers),
      in_total_kmers(in_total_kmers)
{}

double DifferentialTest::poisson_prob(int k, double lambda) // https://en.wikipedia.org/wiki/Poisson_distribution and https://elifesciences.org/articles/32920#appendix-1
{
    if (lambda <= 0) return 0;
    if (k < 0) k = 0;
    return (-lambda + (k * std::log(lambda))); // One does not need to subtract the factorial of k because it divides away in the likelihood ratio.
}



bool DifferentialTest::bonferroni_correction(double& pvalue)
{
    return pvalue < (family_wise_error_rate / total_hypotheses); //  store "family_wise_error_rate / total_hypotheses" as a class variable, corrected_pvalue_threshold (2.4e-8)
}


double DifferentialTest::lrt_threshold(){ // Binary search to find the threshold on the liklihood ratio test
    double corrected_pvalue = family_wise_error_rate/total_hypotheses; // Take into account the Bonferonni multiple testing correction
    if (family_wise_error_rate ==1)// Myrthe: ugly temporair method to put multiple testing off, and simply use a pvalue of 1 instead.
        return (double) 3840 / 1000 / 2;//corrected_pvalue = 0.05;
    Chi2PLookup lookup;
    int low = 0;
    int high = lookup.divisor * lookup.cutoff[0] - 1;
    while (low <= high) {
        int mid = (low + high) >> 1; // the mean of low and high
        double current_pvalue = lookup.getPValue(mid / lookup.divisor, 1);
        if (current_pvalue > corrected_pvalue)
        {low = mid + 1;}
        else
        {high = mid - 1;}
    } // this binary search is not perfect yet: for 0.05 mid = 4000, wheras it should be more like 3840
    return (double)  low / lookup.divisor / 2.0; // Note that 2*LRT is Chi-squared distributed.
}


// adapted code from kmdiff.
// previous is equivalent to mean_control,  mean_case, latter is equivalent to m_sum_controls, m_sum_cases
std::tuple<bool, double> DifferentialTest::likelihood_ratio_test(double in_sum, double out_sum)
{
    double mean = (out_sum + in_sum) / static_cast<double>(out_total_kmers + in_total_kmers);

    double alt_hypothesis = poisson_prob(out_sum, out_sum) +
            poisson_prob(in_sum, in_sum); // K2 equals N2Theta2 because Theta2 = K2/N2 (?). Shouldn't Theta2 be = K2/|background labels| ?

    double null_hypothesis = poisson_prob(out_sum, mean * out_total_kmers) +
            poisson_prob(in_sum, mean * in_total_kmers);

    double likelihood_ratio = alt_hypothesis - null_hypothesis;

    double out_sum_normalized = (double) out_sum * in_total_kmers / out_total_kmers;
    bool sign = false;
    if (out_sum_normalized < in_sum) sign = true;
    if (out_sum_normalized >= in_sum // This checks if the rate for the foreground group is higher compared to the background.
        and likelihood_ratio > likelihood_ratio_threshold)
        std::cout << sign;
    if (out_sum_normalized < in_sum // This checks if the rate for the foreground group is higher compared to the background.
        and likelihood_ratio > likelihood_ratio_threshold)
        return std::make_tuple(true, likelihood_ratio);

//
//
//    if (likelihood_ratio < 0) likelihood_ratio = 0;
//    Chi2PLookup Chi2PLookupTable; // generated with https://github.com/MoseleyBioinformaticsLab/chi2plookup
//    double pvalue = Chi2PLookupTable.getPValue(2 * likelihood_ratio, 1);
//
//

//    bool efficient_test = (likelihood_ratio > likelihood_ratio_threshold);
//    bool normal_test = (sign and bonferroni_correction(pvalue));
//    std::ignore = efficient_test;
//    std::ignore = normal_test;
//
//    if ((likelihood_ratio > likelihood_ratio_threshold) != (sign and bonferroni_correction(pvalue))){
//        std::cout  << "efficient test: " + std::to_string(likelihood_ratio > likelihood_ratio_threshold) << std::endl;  // TODO test this instead of the code below.
//        std::cout << "normal test: " + std::to_string(sign and bonferroni_correction(pvalue)) << std::endl; // The efficient test always returns false it seems.
//    }
//    assert((likelihood_ratio > likelihood_ratio_threshold) == (sign and bonferroni_correction(pvalue)));
//
//    if (sign)
//        if (bonferroni_correction(pvalue))
//            return std::make_tuple(true, likelihood_ratio);

    return std::make_tuple(false, likelihood_ratio);
}

} // namespace mtg
