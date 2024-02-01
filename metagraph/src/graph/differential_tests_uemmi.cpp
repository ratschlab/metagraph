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

DifferentialTest::DifferentialTest(double family_wise_error_rate, size_t total_hypotheses, size_t ,
                                   size_t in_total_kmers, size_t out_total_kmers, size_t in_unique_kmers, size_t out_unique_kmers): 
    family_wise_error_rate(family_wise_error_rate),
    total_hypotheses(total_hypotheses),
    in_total_kmers(in_total_kmers),
    out_total_kmers(out_total_kmers), 
    in_unique_kmers(in_unique_kmers),
    out_unique_kmers(out_unique_kmers)
    {}

// loglikelihood of poisson dist.
double DifferentialTest::poisson_prob(int k, double lambda) // https://en.wikipedia.org/wiki/Poisson_distribution and https://elifesciences.org/articles/32920#appendix-1
{
    if (lambda <= 0) return 0;
    if (k < 0) k = 0;
    return (-lambda + (k * std::log(lambda))); // One does not need to subtract the factorial of k because it divides away in the likelihood ratio.
}

// TODO: implement negative binomial 
// loglikelihood of negative binomial dist.
// p = succes probability = avg occurance of kmers in sampels (either in and out or both together), 
// r = number of failures = number of kmers in out group, 
// k = x = number of successes = number kmers in in group
// needs counts for all samples not just sum or avg


// not in use at the moment
bool DifferentialTest::bonferroni_correction(double& pvalue)
{
    return pvalue < (family_wise_error_rate / total_hypotheses); 
}

// return likelihood ratio thereshold corresponding to corrected alpha value
double DifferentialTest::lrt_threshold(){ // Binary search to find the threshold on the liklihood ratio test 
    double corrected_pvalue = family_wise_error_rate/total_hypotheses; // Take into account the Bonferonni multiple testing correction
    std::cout << "total_hypotheses " << total_hypotheses << "\n";
    if (family_wise_error_rate ==1) // TODO: change how to handle turning off mutliple testing correction. currtenly used like a flag, has no statistical meaning
        return (double) 3840 / 1000 / 2; //corrected_pvalue = 0.05;
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
    } // get index where pvalue is closest to corrected pvalue, from that index get the corresponding likelihood ratio threshold
    return (double)  low / lookup.divisor / 2.0; // Note that 2*LRT is Chi-squared distributed.
}


// adapted code from kmdiff.
// previous is equivalent to mean_control,  mean_case, latter is equivalent to m_sum_controls, m_sum_cases
std::tuple<bool, double> DifferentialTest::likelihood_ratio_test(double in_sum, double out_sum)
{
    // add pseudocounts
    in_sum += 1;
    out_sum += 1;
    double mean = (out_sum + in_sum)*(out_total_kmers+in_total_kmers)/ static_cast<double>(out_unique_kmers + in_unique_kmers);
    // mean = (out_sum + in_sum)/ static_cast<double>(out_total_kmers + in_total_kmers);
    // double lambda = (in_sum/in_total_kmers + out_sum/out_total_kmers)/2 * ((out_total_kmers+in_total_kmers)/2);
    std::cout << "out sum " << out_sum << "\n";
    std::cout << "in sum " << in_sum << "\n";
    std::cout << "mean " << mean << "\n";
    std::cout << "out_total_kmers " << out_total_kmers << "\n";
    std::cout << "in_total_kmers " << in_total_kmers << "\n";
    double alt_hypothesis = poisson_prob(out_sum, out_sum) +
            poisson_prob(in_sum, in_sum); // K2 equals N2Theta2 because Theta2 = K2/N2 (?). Shouldn't Theta2 be = K2/|background labels| ?

    double null_hypothesis = poisson_prob(out_sum, mean * out_total_kmers) +
            poisson_prob(in_sum, mean * in_total_kmers);

    double likelihood_ratio = alt_hypothesis - null_hypothesis;
    double out_sum_normalized = (double) out_sum * in_total_kmers / out_total_kmers;  
    if ((out_sum_normalized < in_sum) && (likelihood_ratio > likelihood_ratio_threshold)){
        return std::make_tuple(true, likelihood_ratio);
    }
    return std::make_tuple(false, likelihood_ratio);
}

} // namespace mtg
