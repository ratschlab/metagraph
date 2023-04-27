#include <math.h>
#include <string>
#include <numbers>
#include <iostream>
#include <cassert>

#include "differential_tests.hpp"
#include "lookup_table_chisqrd_cdf.cpp"



namespace mtg{
DifferentialTest::DifferentialTest(double family_wise_error_rate, size_t total_hypotheses,
                                   double pvalue, size_t in_total_kmers, size_t out_total_kmers) :
      family_wise_error_rate(family_wise_error_rate),
      total_hypotheses(total_hypotheses),
      pvalue(pvalue),
      out_total_kmers(out_total_kmers),
      in_total_kmers(in_total_kmers)
{}

double DifferentialTest::poisson_prob(int k, double lambda)
{
    if (lambda <= 0) return 0;
    if (k < 0) k = 0;
    return (-lambda + (k * std::log(lambda))); // One does not need to subtract the factorial of k, because it divides away in the likelihood ratio.
}


double DifferentialTest::lrt_threshold(){ // Binary search to find the threshold on the liklihood ratio test
    if (family_wise_error_rate) // Take into account the Bonferonni multiple testing correction
        pvalue = family_wise_error_rate/total_hypotheses;
    Chi2PLookup lookup; // Lookup table was created with the library
    int low = 0;
    int high = lookup.divisor * lookup.cutoff[0] - 1;
    while (low <= high) {
        int mid = (low + high) >> 1; // the mean of low and high
        double current_pvalue = lookup.getPValue(mid / lookup.divisor, 1);
        if (current_pvalue > pvalue)
        {low = mid + 1;}
        else
        {high = mid - 1;}
    } // TODO this binary search is not perfect yet: for 0.05 mid = 4000, whereas it should be more like 3840
    return (double)  low / lookup.divisor / 2.0; // Note that 2*LRT is Chi-squared distributed.
}



std::tuple<bool, double> DifferentialTest::likelihood_ratio_test(double in_sum, double out_sum)
{ // This code has been adapted from kmdiff.
    double mean = (out_sum + in_sum) / static_cast<double>(out_total_kmers + in_total_kmers);

    double alt_hypothesis = poisson_prob(out_sum, out_sum) +
            poisson_prob(in_sum, in_sum);

    double null_hypothesis = poisson_prob(out_sum, mean * out_total_kmers) +
            poisson_prob(in_sum, mean * in_total_kmers);

    double likelihood_ratio = alt_hypothesis - null_hypothesis;

    double out_sum_normalized = (double) out_sum * in_total_kmers / out_total_kmers;
    if (out_sum_normalized < in_sum // This checks if the rate for the foreground group is higher compared to the background.
        and likelihood_ratio > likelihood_ratio_threshold)
        return std::make_tuple(true, likelihood_ratio);

    return std::make_tuple(false, likelihood_ratio);
}

} // namespace mtg
