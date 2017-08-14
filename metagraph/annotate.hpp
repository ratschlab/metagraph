#ifndef __ANNOTATE_HPP__
#define __ANNOTATE_HPP__

#include <cstdint>
#include <vector>
#include <set>
#include <string>

#include "kseq.h"
#include "dbg_succinct_libmaus.hpp"

namespace annotate {

    void annotate_seq(DBG_succ* G, kstring_t &seq, kstring_t &label, uint64_t start=0, uint64_t end=0, pthread_mutex_t* anno_mutex=NULL);

    void annotate_kmer(DBG_succ* G, std::string &kmer, uint32_t &label, uint64_t &previous, pthread_mutex_t* anno_mutex, bool ignore=false);

    std::vector<uint32_t> classify_path(DBG_succ* G, std::vector<uint64_t> path);

    std::set<uint32_t> classify_read(DBG_succ* G, kstring_t &read, uint64_t max_distance);

    uint32_t insert_new_combination(std::vector<uint32_t> &combination_vector, std::vector<uint32_t> &new_combination);

    std::vector<uint32_t> add_to_combination(std::vector<uint32_t> &old_combination, uint32_t value);

    std::vector<uint32_t> get_curr_combination(std::vector<uint32_t> &combinations, uint32_t idx);

    void annotationToScreen(DBG_succ* G);
}
#endif
