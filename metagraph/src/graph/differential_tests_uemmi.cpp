#include <math.h>
#include <string>
// #include <numbers>
#include <iostream>
#include <cassert>
#include <algorithm>
#include "common/logger.hpp"

#include "differential_tests.hpp"
#include "lookup_table_chisqrd_cdf.hpp"


// https://github.com/tlemane/kmdiff/blob/6a56ce6f20abbf63928a19ebbfecb1968efd4cd3/include/kmdiff/log_factorial_table.hpp
// cpp https://github.com/tlemane/kmdiff/blob/6a56ce6f20abbf63928a19ebbfecb1968efd4cd3/src/log_factorial_table.cpp
namespace mtg{

DifferentialTest::DifferentialTest(double family_wise_error_rate, size_t total_hypotheses, size_t ,
                                   size_t in_total_kmers, size_t out_total_kmers): 
    family_wise_error_rate(family_wise_error_rate),
    total_hypotheses(total_hypotheses),
    in_total_kmers(in_total_kmers),
    out_total_kmers(out_total_kmers)
    {}

// loglikelihood of poisson dist.
double DifferentialTest::poisson_prob(int k, double lambda) // https://en.wikipedia.org/wiki/Poisson_distribution and https://elifesciences.org/articles/32920#appendix-1
{
    if (lambda <= 0) return 0;
    if (k < 0) k = 0;
    return (-lambda + (k * std::log(lambda))); // One does not need to subtract the factorial of k because it divides away in the likelihood ratio.
}

double DifferentialTest::get_t_test_alpha(int df, double alpha=0.05){
    alpha = alpha/2; // convert to one sided test
    alpha = alpha/total_hypotheses; // bonferroni correction
    common::logger->trace("alpha: {}", alpha);
    double critical_value = t_table.getCriticalValue(alpha, df);
    common::logger->trace("Critical value: {}", critical_value);
    return critical_value;
}

std::vector<double> DifferentialTest::get_midranks(std::vector<double> in_counts, int size_in_counts){
    std::vector<double> mid_ranks;
    std::sort(in_counts.begin(), in_counts.end());
    double curr_value = in_counts[0];
    int tied = 1;
    for (int i = 1; i<size_in_counts; i++){
        if (in_counts[i] == curr_value){
            tied++;
            if (i == size_in_counts-1){
                i += 1;
                for (int j = i-tied; j<i; j++){
                    mid_ranks.push_back((i-tied) + (tied+1)/2.0);
                }
            }
        }
        else{
            if (i == size_in_counts-1){
                for (int j = i-tied; j<i; j++){
                    mid_ranks.push_back((i-tied) + (tied+1)/2.0);
                }
                mid_ranks.push_back(i+1);
            }
            else{
                for (int j = i-tied; j<i; j++){
                    mid_ranks.push_back((i-tied) + (tied+1)/2.0);
                }
                curr_value = in_counts[i];
                tied = 1;
            }
            
        }
    }
    return mid_ranks;
}

int DifferentialTest::get_df_approx(std::vector<double> in_counts, std::vector<double> out_counts){
    int m = in_counts.size();
    int n = out_counts.size();
    double var_in = get_var(in_counts, m);
    double var_out = get_var(out_counts, n);
    if (var_in == 0 && var_out == 0){
        return std::min(m, n) - 1;
    }
    double df_approx = (var_in*var_in/m + var_out*var_out/n)*(var_in*var_in/m + var_out*var_out/n) / (((var_in*var_in/m)*(var_in*var_in/m)/(m - 1)) + ((var_out*var_out/n)*(var_out*var_out/n)/(n - 1)));
    return int(std::floor(df_approx));
}

int DifferentialTest::get_df_conservative(std::vector<double> in_counts, std::vector<double> out_counts){
    int m = in_counts.size();
    int n = out_counts.size();
    return std::min(m, n) - 1;
}

double DifferentialTest::get_var(std::vector<double> counts, int n){
    double mean = std::accumulate(counts.begin(), counts.end(), 0.0) / n;
    double mean_of_sqaures = std::inner_product(counts.begin(), counts.end(), counts.begin(), 0)/n;
    double var = mean_of_sqaures - mean*mean;
    return var;
}

std::tuple<bool, double> DifferentialTest::brunner_munzel_test(std::vector<double> in_counts, std::vector<double> out_counts){
    std::vector<double> mid_ranks_in;
    std::vector<double> mid_ranks_out;
    int m = in_counts.size();
    int n = out_counts.size();
    int N = m + n;
    // populate mid_ranks_in and mid_ranks_out
    mid_ranks_in = get_midranks(in_counts, m);
    mid_ranks_out = get_midranks(out_counts, n);
    // get means of midranks
    double mean_mid_rank_in = std::accumulate(mid_ranks_in.begin(), mid_ranks_in.end(), 0.0) / mid_ranks_in.size();
    double mean_mid_rank_out = std::accumulate(mid_ranks_out.begin(), mid_ranks_out.end(), 0.0) / mid_ranks_out.size();
    // get variances
    double var_in = get_var(in_counts, m);
    double var_out = get_var(out_counts, n);
    if (var_in == 0 && var_out == 0){
        var_0++;
        var_in = 1;
        var_out = 1;
    }         
    // get test statistic
    double b = (mean_mid_rank_out - mean_mid_rank_in) / (N * std::sqrt(var_in/(m*n*n) + var_out/(m*m*n)));
    // get degrees of freedom 
    if (df_precalc == -1){
        common::logger->trace("df_precalc not set, calculating df_precalc");
        df_precalc = get_df_conservative(in_counts, out_counts);
    }
    // get t statistic alpha value
    if (alpha_precalc == -1){
        common::logger->trace("alpha_precalc not set, calculating alpha_precalc");
        alpha_precalc = get_t_test_alpha(df_precalc, 0.05);
    }
    if( std::abs(b) < alpha_precalc){
        return std::tuple(true, b);
    }
    else{
        return std::tuple(false, b);
    }
}


std::tuple<std::vector<int>, std::vector<double>> sorted_pvalues(){
    std::vector<int> nodes_sorted;
    std::vector<double> pvalues_sorted;

    return std::make_tuple(nodes_sorted, pvalues_sorted);
}

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

int DifferentialTest::benjamini_yekutieli(std::vector<std::pair<double, int>> likelihood_ratios, double lrt_threshold){
    int m = likelihood_ratios.size();
    double gamma = 0.5772156649; // euler-mascheroni constant
    double c = std::log(m)  + gamma + 1/(2*m); // taylor expansion for sum(1/i) from 1 to m
    int max_i = 0;
    for (int i = 0; i <= m-1; i++){
        double pvalue = likelihood_ratios[i].first;
        double qvalue = lrt_threshold * (i+1) / (m * c);
        if (pvalue <= qvalue){
            max_i = i;
        } 
    }
    return max_i;
}


// adapted code from kmdiff.
// previous is equivalent to mean_control,  mean_case, latter is equivalent to m_sum_controls, m_sum_cases
std::tuple<bool, double> DifferentialTest::likelihood_ratio_test(double in_sum, double out_sum)
{
    // add pseudocounts
    in_sum += 1;
    out_sum += 1;
    double mean = (out_sum + in_sum)/ static_cast<double>(out_total_kmers + in_total_kmers);
    double alt_hypothesis = poisson_prob(out_sum, out_sum) +
            poisson_prob(in_sum, in_sum); // K2 equals N2Theta2 because Theta2 = K2/N2 (?). Shouldn't Theta2 be = K2/|background labels| ?

    double null_hypothesis = poisson_prob(out_sum, mean * out_total_kmers) +
            poisson_prob(in_sum, mean * in_total_kmers);

    double likelihood_ratio = alt_hypothesis - null_hypothesis;
    // double out_sum_normalized = (double) out_sum * in_total_kmers / out_total_kmers;  
    if (likelihood_ratio > likelihood_ratio_threshold){
        return std::make_tuple(true, likelihood_ratio);
    }
    return std::make_tuple(false, likelihood_ratio);
}

} // namespace mtg
