#ifndef __DATATYPES_HPP__
#define __DATATYPES_HPP__

#include <functional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <pthread.h>
#include <assert.h>
#include "kseq.h"

class DBG_succ;

struct HitInfo {
    uint64_t rl;
    uint64_t ru;
    uint64_t str_pos;
    uint64_t graph_pos;
    uint64_t distance;
    std::string cigar;
    std::vector<uint64_t> path;

    HitInfo(uint64_t rl_, uint64_t ru_, uint64_t str_pos_, uint64_t graph_pos_, uint64_t distance_, std::string cigar_, std::vector<uint64_t> path_):
        rl(rl_),
        ru(ru_),
        str_pos(str_pos_),
        graph_pos(graph_pos_),
        distance(distance_),
        cigar(cigar_),
        path(path_) {}

    HitInfo(const HitInfo& other):
        rl(other.rl),
        ru(other.ru),
        str_pos(other.str_pos),
        graph_pos(other.graph_pos),
        distance(other.distance),
        cigar(other.cigar),
        path(other.path) {}
};

class HitInfoCompare {
    bool is_reverse;
public:
    HitInfoCompare(const bool& is_reverse_ = false) {
        is_reverse = is_reverse_;
    }

    bool operator() (const HitInfo& lhs, const HitInfo& rhs) const {
        if (is_reverse) return (lhs.distance < rhs.distance);
        else return (lhs.distance > rhs.distance);
    }
};

struct AnnotationSet {
    std::set<uint32_t> annotation;
};

struct AnnotationHash {
    std::unordered_map<uint32_t, uint32_t> bitmap;
    uint32_t last = 1;
    
    uint32_t knuth_hash(const uint32_t i) const {
        return i * UINT32_C(2654435761);
    };

    std::uint32_t operator()(const std::vector<uint32_t> &a) const {
        std::vector<uint32_t>::const_iterator it = a.begin();
        uint32_t h1 = knuth_hash(*it);
        it++;
        if (a.size() > 1) {
            for (; it != a.end(); ++it) {
                uint32_t h2 = knuth_hash(*it);
                h1 = h1 ^ (h2 << 1);
            }
        }
        return h1;
    };

    std::uint32_t bithash(const std::vector<uint32_t> &a) {
        uint32_t h1 = 0;
        for (std::vector<uint32_t>::const_iterator it = a.begin(); it!=a.end();++it) {
            std::unordered_map<uint32_t, uint32_t>::iterator c = bitmap.find(*it);
            if (c != bitmap.end()) {
                h1 += c->second;
            } else {
                bitmap[*it] = last;
                h1 += last;
                last <<= 1;
            }
        }
        assert(h1 || !a.size());
        return h1;
    };


    /*std::uint32_t operator()(const std::set<uint32_t> &a) const {
        std::set<uint32_t>::iterator it = a.begin();
        //std::uint32_t h1 = (uint32_t) std::hash<uint32_t>{}(*it);
        uint32_t h1 = knuth_hash(*it);
        it++;
        if (a.size() > 1) {
            for (; it != a.end(); ++it) {
                //std::uint32_t h2 = (uint32_t) std::hash<uint32_t>{}(*it);
                uint32_t h2 = knuth_hash(*it);
                h1 = h1 ^ (h2 << 1);
                //h1 = (h1 * h2);
            }
        }
        return h1;
    }*/
};

struct ParallelMergeContainer {
    std::vector<std::pair<uint64_t, uint64_t> > ref_bins;
    std::vector<std::vector<std::pair<uint64_t, uint64_t> > > bins;
    std::vector<DBG_succ*> result;
    std::vector<DBG_succ*> graphs;
    unsigned int idx;
    unsigned int k;
    unsigned int bins_done;
    unsigned int first = 0;
    unsigned int last = 0;

    /* Helper function to subset the bins to the chunk
     * computed in the current distributed compute.
     */
    void subset_bins(unsigned int idx, unsigned int total, unsigned int bins_per_part){

        //std::cerr << "ref bins " << ref_bins.size() << " total " << total << " per part " << bins_per_part << std::endl;
        assert(ref_bins.size() == (total * bins_per_part));

        std::vector< std::pair<uint64_t, uint64_t> > new_ref_bins;
        //std::cerr << "min: " << binsize_min << " max: " << binsize_max << " thresh: " << threshold << " total: " << total << std::endl;

        /*size_t start, end;
        if (idx < threshold) {
            start = binsize_max * idx;
            end = (idx == (total - 1)) ? ref_bins.size() : binsize_max * (idx + 1);
        } else {
            start = (threshold * binsize_max) + ((idx - threshold) * binsize_min);
            end = (idx == (total - 1)) ? ref_bins.size() : (threshold * binsize_max) + ((idx - threshold + 1) * binsize_min);
        }*/
        size_t start = idx * bins_per_part;
        size_t end = (idx + 1) * bins_per_part;

        if (start > 0)
            first = ref_bins.at(start - 1).second;
        if (end < ref_bins.size())
            last = ref_bins.at(end).second;

        for (size_t i = start; i < end; i++) {
            new_ref_bins.push_back(ref_bins.at(i));
        }

        ref_bins = new_ref_bins;
    }

    void print_bins() {
        for (size_t ii = 0; ii < bins.size(); ii++) {
            std::cerr << "graph " << ii + 1 << std::endl;
            for (size_t i = 0; i < bins.at(ii).size(); i++) 
                std::cerr << bins.at(ii).at(i).first << " - " << bins.at(ii).at(i).second << std::endl;
        }
    }

    /* Show an overview of the distribution of merging bin 
     * sizes.
     */
    void get_bin_stats() {
        size_t min_bin = 0, max_bin = 0, total_bin = 0;

        size_t cum_size;
        for (size_t i = 0; i < ref_bins.size(); ++i) {
            cum_size = 0;
            for (size_t ii = 0; ii < bins.size(); ii++) {
                cum_size += (bins.at(ii).at(i).first == 0) ? 0 : bins.at(ii).at(i).second - bins.at(ii).at(i).first + 1;
            }
            if (cum_size > 0) {
                min_bin = (min_bin == 0) ? cum_size : std::min(min_bin, cum_size);
                max_bin = (max_bin == 0) ? cum_size : std::max(max_bin, cum_size);
            }
            total_bin += cum_size;
        }

        std::cout << std::endl;
        std::cout << "Total number of bins: " << ref_bins.size() << std::endl;
        std::cout << "Total size: " << total_bin << std::endl;
        std::cout << "Smallest bin: " << min_bin << std::endl;
        std::cout << "Largest bin: " << max_bin << std::endl;
        std::cout << "Average bin size: " << total_bin / ref_bins.size() << std::endl << std::endl;
    }
};


struct ParallelMergeContainer_deprecated {
    std::vector<std::pair<uint64_t, uint64_t> > bins_g1;
    std::vector<std::pair<uint64_t, uint64_t> > bins_g2;
    std::vector<DBG_succ*> result;
    DBG_succ* graph1;
    DBG_succ* graph2;
    unsigned int idx;
    unsigned int k;
    unsigned int bins_done;

    /* Helper function to rebalance the bins for
     * a somewhat equal distribution of sizes.
     */
    void rebalance_bins(uint64_t target_bins) {

        std::vector<uint64_t> combined_bins;
        uint64_t total_sum = 0;
        size_t size1, size2;
        for (size_t i = 0; i < bins_g1.size(); ++i) {
            size1 = (bins_g1.at(i).first == 0) ? 0 : bins_g1.at(i).second - bins_g1.at(i).first + 1;
            size2 = (bins_g2.at(i).first == 0) ? 0 : bins_g2.at(i).second - bins_g2.at(i).first + 1;
            combined_bins.push_back(size1 + size2);
            total_sum += combined_bins.back();
            //std::cerr << "bin 1: " << bins_g1.at(i).first << " - " << bins_g1.at(i).second << " size: " << size1 << " --- " << "bin 2: " << bins_g2.at(i).first << " - " << bins_g2.at(i).second << " size: " << size2 << " total: " << size1 + size2 << std::endl; 
        }
        uint64_t target_bin_size = (total_sum / target_bins) + 1;

        std::vector<std::pair<uint64_t, uint64_t> > new_bins_g1;
        std::vector<std::pair<uint64_t, uint64_t> > new_bins_g2;

        uint64_t start_g1 = 0, start_g2 = 0, end_g1 = 0, end_g2 = 0;
        uint64_t cum_sum = 0;
        for (size_t i = 0; i < combined_bins.size(); ++i) {
            cum_sum += combined_bins.at(i);
            if (start_g1 == 0 && bins_g1.at(i).first > 0) {
                start_g1 = bins_g1.at(i).first;
                end_g1 = bins_g1.at(i).second;
            }
            if (start_g2 == 0 && bins_g2.at(i).first > 0) {
                start_g2 = bins_g2.at(i).first;
                end_g2 = bins_g2.at(i).second;
            }
            end_g1 = std::max(end_g1, bins_g1.at(i).second);
            end_g2 = std::max(end_g2, bins_g2.at(i).second);

            if (cum_sum >= target_bin_size) {
                new_bins_g1.push_back(std::make_pair(start_g1, end_g1));
                new_bins_g2.push_back(std::make_pair(start_g2, end_g2));
                cum_sum = 0;
                start_g1 = start_g2 = end_g1 = end_g2 = 0;
            }
        }
        if ((start_g1 > 0) || (start_g2 > 0)) {
            new_bins_g1.push_back(std::make_pair(start_g1, end_g1));
            new_bins_g2.push_back(std::make_pair(start_g2, end_g2));
        }

        bins_g1 = new_bins_g1;
        bins_g2 = new_bins_g2;
    }


    /* Show an overview of the distribution of merging bin 
     * sizes.
     */
    void get_bin_stats() {
        size_t min_bin = 0, max_bin = 0, total_bin = 0;
        size_t curr_size, size1, size2;
        for (size_t i = 0; i < bins_g1.size(); ++i) {
            size1 = (bins_g1.at(i).first == 0) ? 0 : bins_g1.at(i).second - bins_g1.at(i).first + 1;
            size2 = (bins_g2.at(i).first == 0) ? 0 : bins_g2.at(i).second - bins_g2.at(i).first + 1;
            curr_size = (size1 + size2);
            if (curr_size > 0) {
                min_bin = (min_bin == 0) ? curr_size : std::min(min_bin, curr_size);
                max_bin = (max_bin == 0) ? curr_size : std::max(max_bin, curr_size);
            }
            total_bin += curr_size;
        }

        std::cout << std::endl;
        std::cout << "Total number of bins: " << bins_g1.size() << std::endl;
        std::cout << "Total size: " << total_bin << std::endl;
        std::cout << "Smallest bin: " << min_bin << std::endl;
        std::cout << "Largest bin: " << max_bin << std::endl;
        std::cout << "Average bin size: " << total_bin / bins_g1.size() << std::endl << std::endl;
    }

    /* Helper function to subset the bins to the chunk
     * computed in the current distributed compute.
     */
    void subset_bins(unsigned int idx, unsigned int total) {
        if (bins_g1.size() < total)
            std::cerr << "total number of possible bins is " << bins_g1.size() << std::endl
                      << "current choice: " << total << std::endl
                      << "Please use at max " << bins_g1.size() << " total bins" << std::endl;
        assert(bins_g1.size() >= total);
        // augment size of last bin until end of bins
        size_t binsize_min = bins_g1.size() / total;
        size_t binsize_max = (bins_g1.size() + total - 1) / total;
        size_t threshold = bins_g1.size() - (total * binsize_min);

        std::vector<std::pair<uint64_t, uint64_t> > new_bins_g1;
        std::vector<std::pair<uint64_t, uint64_t> > new_bins_g2;

        size_t start;
        size_t end;
        if (idx < threshold) {
            start = binsize_max * idx;
            end = (idx == (total - 1)) ? bins_g1.size() : binsize_max * (idx + 1);
        } else {
            start = (threshold * binsize_max) + ((idx - threshold) * binsize_min);
            end = (idx == (total - 1)) ? bins_g1.size() : (threshold * binsize_max) + ((idx - threshold + 1) * binsize_min);
        }

        for (size_t i = start; i < end; i++) {
            new_bins_g1.push_back(bins_g1.at(i));
            new_bins_g2.push_back(bins_g2.at(i));
        }

        bins_g1 = new_bins_g1;
        bins_g2 = new_bins_g2;
    }
};

struct ParallelAnnotateContainer {
    kstring_t* seq;
    kstring_t* label;
    DBG_succ* graph;
    uint64_t idx;
    uint64_t binsize;
    uint64_t total_bins; 
    //pthread_mutex_t* anno_mutex;
};

#endif
