#include "graph_cleaning.hpp"

#include <cmath>

#include "common/logger.hpp"

using mtg::common::logger;

bool is_unreliable_unitig(const std::vector<SequenceGraph::node_index> &path,
                          const NodeWeights &node_weights,
                          uint64_t min_median_abundance) {
    assert(path.size());

    if (min_median_abundance <= 1)
        return false;

    uint64_t num_weak_kmers = 0;

    for (auto node : path) {
        num_weak_kmers += node_weights[node] < min_median_abundance;
    }

    assert(num_weak_kmers <= path.size());

    // check if median is smaller than the threshold
    return num_weak_kmers * 2 > path.size();
}


int cleaning_pick_kmer_threshold(const uint64_t *kmer_covg, size_t arrlen,
                                 double *alpha_est_ptr, double *beta_est_ptr,
                                 double *false_pos_ptr, double *false_neg_ptr);

uint64_t estimate_min_kmer_abundance(const DeBruijnGraph &graph,
                                     const NodeWeights &node_weights,
                                     int fallback_cutoff,
                                     uint64_t num_singleton_kmers) {
    std::vector<uint64_t> hist;
    graph.call_nodes([&](auto i) {
        uint64_t kmer_count = node_weights[i];
        assert(kmer_count && "All k-mers in graph must have non-zero counts");
        while (kmer_count >= hist.size()) {
            hist.push_back(0);
        }
        hist[kmer_count]++;
    });

    hist.resize(std::max((uint64_t)hist.size(), (uint64_t)10), 0);

    if (num_singleton_kmers) {
        logger->info("The count for singleton k-mers in histogram is reset from {} to {}",
                  hist[1], num_singleton_kmers);
        hist[1] = num_singleton_kmers;
    }

    double alpha_est_ptr, beta_est_ptr, false_pos_ptr, false_neg_ptr;
    int cutoff = cleaning_pick_kmer_threshold(hist.data(), hist.size(),
                                              &alpha_est_ptr, &beta_est_ptr,
                                              &false_pos_ptr, &false_neg_ptr);

    if (cutoff != -1)
        return cutoff;
    if (fallback_cutoff == -1) {
        logger->error("Cannot estimate expected minimum k-mer abundance "
                "and fallback is disabled (--fallback -1). Terminating.");
        std::exit(129);
    }
    logger->warn("Cannot estimate expected minimum k-mer abundance. "
            "Using fallback value: {}", fallback_cutoff);

    return fallback_cutoff;
}


/** The next routines were copied from:
 * https://github.com/mcveanlab/mccortex/blob/97aba198d632ee98ac1aa496db33d1a7a8cb7e51/src/tools/clean_graph.c

The MIT License (MIT)

Copyright (c) 2014  Isaac Turner <turner.isaac@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

 */
// Find cutoff by finding first coverage level where errors make up less than
// `fdr` of total coverage
// returns -1 if not found
static inline int pick_cutoff_with_fdr_thresh(const double *e_covg,
                                              const uint64_t *kmer_covg,
                                              size_t arrlen, double fdr)
{
  size_t i;
  for(i = 1; i < arrlen; i++) {
    // printf(" %zu: %f %zu test: %f < %f\n", i, e_covg[i], kmer_covg[i],
    //                                       e_covg[i] / kmer_covg[i], fdr);
    if(e_covg[i] / kmer_covg[i] <= fdr) {
      return i;
    }
  }
  return -1;
}

// Get highest cutoff where false-positives < false-negatives
// i.e. proportion of real kmers we are removing is less than the
//      proportion of bad kmers we are keeping
// returns -1 if not found
static inline int pick_cutoff_FP_lt_FN(const double *e_covg, double e_total,
                                      const uint64_t *kmer_covg, uint64_t d_total,
                                      size_t arrlen)
{
  size_t i;
  // for(i = 1; i < arrlen; i++) { printf("  %zu", kmer_covg[i]); } printf("\n");
  // for(i = 1; i < arrlen; i++) { printf("  %f", e_covg[i]); } printf("\n");
  // printf(" e_total: %f d_total: %f\n", e_total, (double)d_total);
  double e_rem = e_total, d_rem = d_total;
  double e_sum = 0, d_sum = 0;
  for(i = 1; i < arrlen; i++) {
    e_sum += e_covg[i];
    d_sum += kmer_covg[i];
    e_rem -= e_covg[i];
    d_rem -= kmer_covg[i];
    // printf(" %zu: e_total: %f d_total: %f\n", i, e_rem, d_rem);
    if(1-e_sum/d_sum > e_rem/d_rem) {
      return i;
    }
  }
  return -1;
}

static inline int pick_cutoff_loss_vs_error(const double *e_covg,
                                            double e_total,
                                            const uint64_t *kmer_covg,
                                            size_t arrlen)
{
  size_t i;
  // for(i = 1; i < arrlen; i++) { printf("  %zu", kmer_covg[i]); } printf("\n");
  // for(i = 1; i < arrlen; i++) { printf("  %f", e_covg[i]); } printf("\n");
  // printf(" e_total: %f d_total: %f\n", e_total, (double)d_total);
  double e_rem = e_total;
  double e_sum = 0, d_sum = 0;
  for(i = 1; i < arrlen; i++) {
    e_sum += e_covg[i];
    d_sum += kmer_covg[i];
    e_rem -= e_covg[i];
    double lost_seq = (d_sum-e_sum);
    double rem_err = e_rem;
    // printf(" %zu: e_total: %f d_total: %f\n", i, e_rem, d_rem);
    if(lost_seq > rem_err) return i;
  }
  return -1;
}

static inline void cutoff_get_FP_FN(const double *e_covg, double e_total,
                                    const uint64_t *kmer_covg, uint64_t d_total,
                                    size_t cutoff,
                                    double *false_pos, double *false_neg)
{
  size_t i;
  double e_rem = e_total, d_rem = d_total;
  double e_sum = 0, d_sum = 0;
  for(i = 1; i < cutoff; i++) {
    e_sum += e_covg[i];
    d_sum += kmer_covg[i];
    e_rem -= e_covg[i];
    d_rem -= kmer_covg[i];
  }
  *false_pos = 1-e_sum/d_sum;
  *false_neg = e_rem/d_rem;
}

// Check if at least `frac_covg_kept` coverage is kept when using threshold
static inline bool is_cutoff_good(const uint64_t *kmer_covg, size_t arrlen,
                                  size_t cutoff, double frac_covg_kept)
{
  uint64_t kmers_below = 0, kmers_above = 0;
  size_t i;
  for(i = 0;      i < cutoff; i++) kmers_below += kmer_covg[i]*i;
  for(i = cutoff; i < arrlen; i++) kmers_above += kmer_covg[i]*i;

  // At least 20% of kmers should be kept
  return !arrlen || // any cutoff is good if no kmers
         ((double)kmers_above/(kmers_below+kmers_above) >= frac_covg_kept);
}

/**
 * Pick a cleaning threshold from kmer coverage histogram. Assumes low coverage
 * kmers are all due to error. Fits a poisson with a gamma distributed mean.
 * Then chooses a cleaning threshold such than FDR (uncleaned kmers) occur at a
 * rate of < the FDR paramater.
 *
 * Translated from Gil McVean's initial proposed method in R code
 *
 * @param kmer_covg Histogram of kmer counts at coverages 1,2,.. arrlen-1
 * @param arrlen    Length of array kmer_covg
 * @param alpha_est_ptr If not NULL, used to return estimate for alpha
 * @param beta_est_ptr  If not NULL, used to return estimate for beta
 * @return -1 if no cut-off satisfies FDR, otherwise returns coverage cutoff
 */
int cleaning_pick_kmer_threshold(const uint64_t *kmer_covg, size_t arrlen,
                                 double *alpha_est_ptr, double *beta_est_ptr,
                                 double *false_pos_ptr, double *false_neg_ptr)
{
  assert(arrlen >= 10);
  assert(kmer_covg[0] == 0 && "Shouldn't see any kmers with coverage zero");

  size_t i, min_a_est_idx = 0;
  double r1, r2, rr, min_a_est = std::numeric_limits<double>::max(), tmp;
  double aa, faa, a_est, b_est, c0;

  if (mtg::common::get_verbose())
  {
    std::cout << "k-mer count histogram:\n";
    for(i = 1; i < arrlen; i++)
    {
        if (kmer_covg[i] != 0) {
            std::cout << i << ": " << kmer_covg[i] << ", ";
        }
    }
    std::cout << std::endl;
  }

  r1 = (double)kmer_covg[2] / kmer_covg[1];
  r2 = (double)kmer_covg[3] / kmer_covg[2];
  rr = r2 / r1;

  // printf("r1: %.2f r2: %.2f rr: %.2f\n", r1, r2, rr);

  // iterate aa = { 0.01, 0.02, ..., 1.99, 2.00 }
  // find aa value that minimises abs(faa-rr)
  for(i = 1; i <= 200; i++)
  {
    aa = i*0.01;
    faa = std::tgamma(aa) * std::tgamma(aa+2) / (2 * std::pow(std::tgamma(aa+1), 2));
    tmp = std::fabs(faa-rr);
    if(tmp < min_a_est) { min_a_est = tmp; min_a_est_idx = i; }
  }

  // a_est, b_est are estimates for alpha, beta of gamma distribution
  a_est = min_a_est_idx*0.01;
  b_est = std::tgamma(a_est + 1.0) / (r1 * std::tgamma(a_est)) - 1.0;
  b_est = std::max(b_est, 1.); // Avoid beta values <1
  c0 = kmer_covg[1] * std::pow(b_est/(1+b_est),-a_est);

  if(alpha_est_ptr) *alpha_est_ptr = a_est;
  if(beta_est_ptr)  *beta_est_ptr  = b_est;

  // printf("min_a_est_idx: %zu\n", min_a_est_idx);
  // printf("a_est: %f b_est %f c0: %f\n", a_est, b_est, c0);

  // keep coverage estimates on the stack - this should be ok
  double e_covg_tmp, e_covg[arrlen];
  double e_total = 0;
  uint64_t d_total = 0;

  // Calculate some values here for speed
  double log_b_est          = std::log(b_est);
  double log_one_plus_b_est = std::log(1 + b_est);
  double lgamma_a_est       = std::lgamma(a_est);

  // note: lfactorial(x) = lgamma(x+1)

  for(i = 1; i < arrlen; i++)
  {
    e_covg_tmp = a_est * log_b_est - lgamma_a_est - std::lgamma(i)
                   + std::lgamma(a_est + i - 1)
                   - (a_est + i - 1) * log_one_plus_b_est;
    e_covg[i] = std::exp(e_covg_tmp) * c0;
    e_total += e_covg[i];
    d_total += kmer_covg[i];
  }

  // for(i = 1; i < MIN2(arrlen,100); i++)
  //   printf("  %zu: %f %zu\n", i, e_covg[i], (size_t)kmer_covg[i]);

  int cutoff = -1;

  // Find cutoff by finding first coverage level where errors make up less than
  // 0.1% of total coverage
  cutoff = pick_cutoff_with_fdr_thresh(e_covg, kmer_covg, arrlen, 0.001);
  // printf("A cutoff: %i\n", cutoff);

  // Pick highest cutoff that keeps FP < FN
  if(cutoff < 0)
    cutoff = pick_cutoff_FP_lt_FN(e_covg, e_total, kmer_covg, d_total, arrlen);

  if(cutoff < 0)
    cutoff = pick_cutoff_loss_vs_error(e_covg, e_total, kmer_covg, arrlen);

  // printf("B cutoff: %i\n", cutoff);

  if(cutoff < 0) return -1;

  // printf("C cutoff: %i\n", cutoff);

  // Check cutoff keeps at least 20% of coverage
  // (WGS should be much higher, Exome sequencing needs low cutoff)
  if(!is_cutoff_good(kmer_covg, arrlen, cutoff, 0.2)) return -1;

  // printf("D cutoff: %i\n", cutoff);

  // Calculate FP,FN rates
  if(false_pos_ptr || false_neg_ptr) {
    double false_pos = 0, false_neg = 0;
    cutoff_get_FP_FN(e_covg, e_total, kmer_covg, d_total, cutoff,
                     &false_pos, &false_neg);
    // printf("  FP: %f, FN: %f\n", false_pos, false_neg);
    if(false_pos_ptr) *false_pos_ptr = false_pos;
    if(false_neg_ptr) *false_neg_ptr = false_neg;
  }

  // printf(" kmers_above : %zu / (%zu + %zu) = %f\n",
  //        kmers_above, kmers_below, kmers_above,
  //        (double)kmers_above/(kmers_below+kmers_above));

  // printf("cutoff: %i\n", cutoff);

  // printf(" cutoff: %zu fdr: %f fdr_limit: %f good: %i\n",
  //        cutoff, fdr, fdr_limit, (int)good_cutoff);

  return cutoff;
}
