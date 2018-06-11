#ifndef __READS_FILTERING_HPP__
#define __READS_FILTERING_HPP__

#include "utils.hpp"


typedef std::function<void(const std::string&)> CallbackRead;


std::vector<bool> filter_reads(std::function<void(CallbackRead)> generate_reads,
                               size_t k,
                               size_t noise_kmer_frequency,
                               bool verbose,
                               utils::ThreadPool *thread_pool,
                               const std::string &kmc_base = "");


#endif // __READS_FILTERING_HPP__
