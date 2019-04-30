#ifndef __READS_FILTERING_HPP__
#define __READS_FILTERING_HPP__

#include "threading.hpp"


typedef std::function<void(const std::string&)> CallbackRead;


std::vector<bool> filter_reads(std::function<void(CallbackRead)> generate_reads,
                               size_t k,
                               size_t max_unreliable_abundance,
                               size_t unreliable_kmers_threshold,
                               bool verbose,
                               ThreadPool *thread_pool,
                               const std::string &kmc_base = "");


#endif // __READS_FILTERING_HPP__
