#ifndef __DATATYPES_HPP__
#define __DATATYPES_HPP__

#include <functional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <pthread.h>
#include <assert.h>
#include "kseq.h"

#include <sdsl/wavelet_trees.hpp>

#include <libmaus2/bitbtree/bitbtree.hpp>
#include <libmaus2/wavelet/DynamicWaveletTree.hpp>

class DBG_succ;

struct ColorNode {
    
    ColorNode() {
        left = NULL;
        right = NULL;
        idx = 0;
    }

    ColorNode* left;
    ColorNode* right;
    uint64_t idx;
    
};

/*class ColorTree {
  
    private:
    ColorNode* root;

    public:
    ColorTree() {
        root = new ColorNode();
    }

    void add_color(sdsl::vector_rrr<63>* b) {
        for (size_t i = 0; i < b->size(); ++i) {
            if (b[i]) {
                curr = left ? left : (new ColorNode());
            } else {
                curr = right ? right : (new ColorNode());
            }
        }
    }

    uint64_t color_index(sdsl::vector_rrr<63>* b) {

    }

};
*/

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

class bit_vector {

public:

    virtual ~bit_vector() {};

    virtual uint64_t size() = 0;
    virtual void set(size_t id, bool val) = 0;
    virtual void setBitQuick(size_t id, bool val) = 0;
    //virtual bool const& operator[](size_t id) const = 0;
    virtual bool operator[](size_t id) = 0;
    virtual void insertBit(size_t id, bool val) = 0;
    virtual void deleteBit(size_t id) = 0;
    virtual void serialise(std::ostream &out) = 0;
    virtual void deserialise(std::istream &in) = 0;
    virtual uint64_t select1(size_t id) = 0;
    virtual uint64_t rank1(size_t id) = 0;
};


class bit_vector_dyn: public bit_vector, public libmaus2::bitbtree::BitBTree<6, 64> {

public:
    bit_vector_dyn() : libmaus2::bitbtree::BitBTree<6, 64>() {
    }

    bit_vector_dyn(size_t size, bool def) : libmaus2::bitbtree::BitBTree<6, 64>(size, def) {
    }

    ~bit_vector_dyn() {
    }

    bit_vector_dyn(std::vector<uint8_t> &v) : libmaus2::bitbtree::BitBTree<6, 64>(v.size(), false) {
        for (size_t i = 0; i < v.size(); ++i) {
            if (v.at(i))
                this->set(i, true);
        }
    }

    bit_vector_dyn(std::vector<bool> &v) : libmaus2::bitbtree::BitBTree<6, 64>(v.size(), false) {
        for (size_t i = 0; i < v.size(); ++i) {
            if (v.at(i))
                this->set(i, true);
        }
    }


    bit_vector_dyn(bit_vector* v) : libmaus2::bitbtree::BitBTree<6, 64>(v->size(), false) {
        for (size_t i = 0; i < v->size(); ++i) {
            if (v->operator[](i))
                this->set(i, true);
        }
    }

    bit_vector_dyn(std::istream &in) {
        this->deserialise(in);
    }

    uint64_t size() {
        return this->libmaus2::bitbtree::BitBTree<6, 64>::size();
    }

    void set(size_t id, bool val) {
        this->libmaus2::bitbtree::BitBTree<6, 64>::set(id, val);
    }
    
    void setBitQuick(size_t id, bool val) {
        this->libmaus2::bitbtree::BitBTree<6, 64>::setBitQuick(id, val);
    }
    
    //bool const& operator[](size_t id) const {
    //    return this->libmaus2::bitbtree::BitBTree<6, 64>::operator[](id);
    //}

    bool operator[](size_t id) {
        return this->libmaus2::bitbtree::BitBTree<6, 64>::operator[](id);
    }

    void insertBit(size_t id, bool val) {
        this->libmaus2::bitbtree::BitBTree<6, 64>::insertBit(id, val);
    }

    void deleteBit(size_t id) {
        this->libmaus2::bitbtree::BitBTree<6, 64>::deleteBit(id);
    }

    void deserialise(std::istream &in) {
        this->libmaus2::bitbtree::BitBTree<6, 64>::deserialise(in);
    }

    void serialise(std::ostream &out) {
        this->libmaus2::bitbtree::BitBTree<6, 64>::serialise(out);
    }

    uint64_t select1(size_t id) {
        return this->libmaus2::bitbtree::BitBTree<6, 64>::select1(id);
    }

    uint64_t rank1(size_t id) {
        return this->libmaus2::bitbtree::BitBTree<6, 64>::rank1(id);
    }
};

class bit_vector_stat: public bit_vector, public sdsl::bit_vector {

private:
    sdsl::rank_support_v5<> rk;
    sdsl::select_support_mcl<> slct;
    bool update_rs;
    void init_rs() {
        rk = sdsl::rank_support_v5<>(this);
        slct = sdsl::select_support_mcl<>(this);
        update_rs = false;
    }
public:
    bit_vector_stat(size_t size, bool def) : sdsl::bit_vector(size, def) {
        init_rs();
    }

    bit_vector_stat(bit_vector* V) : sdsl::bit_vector(V->size(), 0) {
        for (size_t i = 0; i < V->size(); ++i) {
            if (V->operator[](i))
                this->sdsl::bit_vector::operator[](i) = true;
        }
        init_rs();
    }

    bit_vector_stat(std::istream &in) {
        this->deserialise(in);
        init_rs();
    }

    bit_vector_stat() : sdsl::bit_vector() {
        init_rs();
    }

    ~bit_vector_stat() {
    }

    void set(size_t id, bool val) {
        this->sdsl::bit_vector::operator[](id) = val;
        update_rs = true;
    }

    void setBitQuick(size_t id, bool val) {
        this->sdsl::bit_vector::operator[](id) = val;
        update_rs = true;
    }

    //bool const& operator[](size_t id) const {
    //    return this->operator[](id);
    //}

    bool operator[](size_t id) {
        update_rs = true;
        return this->sdsl::bit_vector::operator[](id);
    }

    void insertBit(size_t id, bool val) {
        this->resize(this->size()+1);
        if (this->size() > 1)
            std::copy_backward(this->begin()+id,(this->end())-1,this->end());
        set(id, val);
        update_rs = true;
    }
    void deleteBit(size_t id) {
        if (this->size() > 1)
            std::copy(this->begin()+id+1,this->end(),this->begin()+id);
        this->resize(this->size()-1);
        update_rs = true;
    }
    void deserialise(std::istream &in) {
        this->load(in);
    }
    void serialise(std::ostream &out) {
        this->serialize(out);
    }
    uint64_t select1(size_t id) {
        if (update_rs)
            init_rs();
        //compensating for libmaus weirdness
        //id++;
        //size_t maxrank = rk(this->size());
        //if (id > maxrank) {
            //TODO: should this line ever be reached?
        //    return this->size();
        //}
        return slct(id + 1);
    }
    uint64_t rank1(size_t id) {
        if (update_rs)
            init_rs();
        //the rank method in SDSL does not include id in the count
        return rk(id >= this->size() ? this->size() : id + 1);
    }
    uint64_t size() {
        return this->sdsl::bit_vector::size();
    }
};


class wavelet_tree {
    
public:

    virtual ~wavelet_tree() {};

    virtual uint64_t size() = 0;
    //virtual void deserialise(std::istream &in) = 0;
    virtual void serialise(std::ostream &out) = 0;
    virtual void insert(uint64_t val, size_t id) = 0;
    virtual void remove(size_t id) = 0;
    //virtual void set(size_t id, uint64_t val) = 0;
    virtual uint64_t rank(uint64_t c, uint64_t i) = 0;
    virtual uint64_t select(uint64_t c, uint64_t i) = 0;
    //virtual uint64_t const& operator[](size_t id) const = 0;
    virtual uint64_t operator[](size_t id) = 0;
    //virtual uint64_t get(size_t id) = 0;
    // this only makes sense when implemented in the dynamic part
    virtual bool get_bit_raw(uint64_t id) = 0;
    
};


class wavelet_tree_stat: public wavelet_tree, public sdsl::int_vector<> {

private:
    sdsl::wt_int<> wwt;
    size_t logsigma;
    bool update_rs;
    void init_wt() {
        std::cout << "Initializing WT index" << std::endl;
        this->resize(n);
        sdsl::construct_im(wwt, *this);
        update_rs=false;
    }
    size_t n;

public:
    wavelet_tree_stat(size_t logsigma, size_t size, uint64_t def) 
        : sdsl::int_vector<>(2 * size + 1, def, 1<< logsigma) {
        n = size;
        update_rs = true;
    }

    wavelet_tree_stat(wavelet_tree* T, size_t logsigma)
        : sdsl::int_vector<>(T->size(), 0, logsigma) {
        n = T->size();
        for (size_t i = 0; i < n; ++i)
            this->sdsl::int_vector<>::operator[](i) = T->operator[](i);
        update_rs = true;
    }

    wavelet_tree_stat(size_t logsigma)
        : wavelet_tree_stat(logsigma, 0, 0) {
    }

    wavelet_tree_stat(std::istream &in) {
        this->deserialise(in);
        // we assume the wwt has been build before serialization
        update_rs = false;
    }

    ~wavelet_tree_stat() {
    }

    uint64_t size() {
        return n;
    }

    void deserialise(std::istream &in) {
        wwt.load(in);
        this->load(in);
        n = this->sdsl::int_vector<>::size();
    }

    void serialise(std::ostream &out) {
        this->resize(n);
        if (update_rs)
            init_wt();
        wwt.serialize(out);
        this->serialize(out);
    }

    void insert(uint64_t val, size_t id) {
        if (n == this->size()) {
            this->resize(2*n+1);
        }
        n++;
        if (this->size() > 1)
            std::copy_backward(this->begin()+id,this->begin()+n-1,this->begin()+n);
        this->sdsl::int_vector<>::operator[](id) = val;
        update_rs = true;
    }

    void remove(size_t id) {
        if (this->size() > 1)
            std::copy(this->begin()+id+1,this->begin()+n,this->begin()+id);
        n--;
        update_rs = true;
    }

    uint64_t rank(uint64_t c, uint64_t i) {
        if (update_rs)
            init_wt();
        return wwt.rank(i >= wwt.size() ? wwt.size() : i+1, c);
    }

    uint64_t select(uint64_t c, uint64_t i) {

        if (update_rs)
            init_wt();
        i++;
        uint64_t maxrank = wwt.rank(wwt.size(), c);
        if (i > maxrank) {
            //TODO: should this line ever be reached?
            return wwt.size();
        }
        return wwt.select(i, c);
    }

    /*sdsl::int_vector<>::reference operator[](const size_t id) {
        update_rs = true;
        return sdsl::int_vector<>::operator[](id);
    }*/

    //uint64_t const& operator[](size_t id) const {
    //    return this->operator[](id);
    //}

    uint64_t operator[](size_t id) {
//        update_rs = true;
        return this->sdsl::int_vector<>::operator[](id);
    }

/*    void set(size_t id, uint64_t val) {
        update_rs = true;
        this->sdsl::int_vector<>::operator[](id) = val;
    }

    uint64_t get(size_t id) {
        return this->sdsl::int_vector<>::operator[](id);
    }
    */

    bool get_bit_raw(uint64_t id) {
       throw std::logic_error("Not Implemented"); 
       return id;
    }

};


//Child of DynamicWaveletTree allowing for construction from an int vector
class wavelet_tree_dyn : public wavelet_tree, public libmaus2::wavelet::DynamicWaveletTree<6, 64> {

private:

    libmaus2::bitbtree::BitBTree<6, 64>* makeTree(std::vector<uint8_t> &W_stat, size_t const b, unsigned int parallel=1) {

        size_t const n = W_stat.size();

        // compute total offsets for the individual bins
        std::vector<uint64_t> offsets((1ull << (b-1)) - 1, 0);
        #pragma omp parallel num_threads(parallel)
        {
            uint64_t v, m, o;
            #pragma omp for
            for (size_t i = 0; i < n; ++i) {
                m = (1ull << (b - 1));
                v = (uint64_t) W_stat.at(i);
                o = 0;
                for (uint64_t ib = 1; ib < b; ++ib) {
                    bool const bit = m & v;
                    if (!bit) {
                        #pragma omp critical
                        offsets.at(o) += 1;
                    }
                    o = 2*o + 1 + bit;
                    m >>= 1;
                }
            }
        }
        //bit_vector *tmp = new bit_vector_dyn(n * b, false);
        libmaus2::bitbtree::BitBTree<6, 64> *tmp = new libmaus2::bitbtree::BitBTree<6, 64>(n * b, false);
        std::vector<uint64_t> upto_offsets ((1ull << (b - 1)) - 1, 0);

        #pragma omp parallel num_threads(parallel)
        {
            uint64_t m, v, o, p, co;
            bool bit;
            #pragma omp for
            for (uint64_t i = 0; i < n; ++i) {
                m = (1ull << (b - 1));
                v = (uint64_t) W_stat.at(i);
                o = 0;
                p = i;
                co = 0;
                for (uint64_t ib = 0; ib < b - 1; ++ib) {
                    bit = m & v;
                    if (bit) {
                        #pragma omp critical
                        tmp->setBitQuick(ib * n + p + co, true);
                        co += offsets.at(o);
                        p -= upto_offsets.at(o);
                    } else {
                        p -= (p - upto_offsets.at(o)); 
                        upto_offsets.at(o) += 1;
                    }
                    //std::cerr << "o: " << o << " offset[o]: " << offsets.at(o) << std::endl;
                    o = 2*o + 1 + bit;
                    m >>= 1;
                }
                bit = m & v;
                if (bit) {
                    //std::cerr << "b - 1: " << b - 1 << " n: " << n << " p: " << p << " co: " << co << std::endl;
                    #pragma omp critical
                    tmp->setBitQuick((b - 1) * n + p + co, true); 
                }
            }
        }
        return tmp;
    }

    libmaus2::bitbtree::BitBTree<6, 64>* makeTree(wavelet_tree* W_stat, size_t const b, unsigned int parallel=1) {

        size_t const n = W_stat->size();

        // compute total offsets for the individual bins
        std::vector<uint64_t> offsets((1ull << (b-1)) - 1, 0);
        #pragma omp parallel num_threads(parallel)
        {
            uint64_t v, m, o;
            #pragma omp for
            for (size_t i = 0; i < n; ++i) {
                m = (1ull << (b - 1));
                v = (uint64_t) W_stat->operator[](i);
                o = 0;
                for (uint64_t ib = 1; ib < b; ++ib) {
                    bool const bit = m & v;
                    if (!bit) {
                        #pragma omp critical
                        offsets.at(o) += 1;
                    }
                    o = 2*o + 1 + bit;
                    m >>= 1;
                }
            }
        }
        //bit_vector *tmp = new bit_vector_dyn(n * b, false);
        libmaus2::bitbtree::BitBTree<6, 64> *tmp = new libmaus2::bitbtree::BitBTree<6, 64>(n * b, false);
        std::vector<uint64_t> upto_offsets ((1ull << (b - 1)) - 1, 0);

        #pragma omp parallel num_threads(parallel)
        {
            uint64_t m, v, o, p, co;
            bool bit;
            #pragma omp for
            for (uint64_t i = 0; i < n; ++i) {
                m = (1ull << (b - 1));
                v = (uint64_t) W_stat->operator[](i);
                o = 0;
                p = i;
                co = 0;
                for (uint64_t ib = 0; ib < b - 1; ++ib) {
                    bit = m & v;
                    if (bit) {
                        #pragma omp critical
                        tmp->setBitQuick(ib * n + p + co, true);
                        co += offsets.at(o);
                        p -= upto_offsets.at(o);
                    } else {
                        p -= (p - upto_offsets.at(o)); 
                        upto_offsets.at(o) += 1;
                    }
                    //std::cerr << "o: " << o << " offset[o]: " << offsets.at(o) << std::endl;
                    o = 2*o + 1 + bit;
                    m >>= 1;
                }
                bit = m & v;
                if (bit) {
                    //std::cerr << "b - 1: " << b - 1 << " n: " << n << " p: " << p << " co: " << co << std::endl;
                    #pragma omp critical
                    tmp->setBitQuick((b - 1) * n + p + co, true); 
                }
            }
        }
        return tmp;
    }

public:
    wavelet_tree_dyn(size_t b)
        : libmaus2::wavelet::DynamicWaveletTree<6, 64>(b) {
    }

    wavelet_tree_dyn(std::istream &in)
        : libmaus2::wavelet::DynamicWaveletTree<6, 64>(in) {
    }

    wavelet_tree_dyn(std::vector<uint8_t> &W_stat, size_t b, unsigned int parallel=1)
        : libmaus2::wavelet::DynamicWaveletTree<6, 64>(makeTree(W_stat, b, parallel), b, W_stat.size()) {
    }

    wavelet_tree_dyn(wavelet_tree* W_stat, size_t b, unsigned int parallel=1)
        : libmaus2::wavelet::DynamicWaveletTree<6, 64>(makeTree(W_stat, b, parallel), b, W_stat->size()) {
    }

    //wavelet_tree_dyn(std::vector<uint8_t> &W_stat, size_t b, unsigned int parallel=1)
    //    : libmaus2::wavelet::DynamicWaveletTree<6, 64>(dynamic_cast<libmaus2::bitbtree::BitBTree<6, 64>*>(makeTree(W_stat, b, parallel)), b, W_stat.size()) {
    //}


    wavelet_tree_dyn(libmaus2::bitbtree::BitBTree<6, 64>* bt, size_t b, size_t n)
        : libmaus2::wavelet::DynamicWaveletTree<6, 64>(bt, b, n) {
    }
    //wavelet_tree_dyn(bit_vector* bt, size_t b, size_t n)
    //    : libmaus2::wavelet::DynamicWaveletTree<6, 64>(dynamic_cast<libmaus2::bitbtree::BitBTree<6, 64>*>(bt), b, n) {
    //}

    ~wavelet_tree_dyn() {}

    uint64_t size() {
        return this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::size();
    }

    //void deserialise(std::istream &in) {
    //    this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::deserialise(in);
    //}

    void serialise(std::ostream &out) {
        this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::serialise(out);
    }

    void insert(uint64_t val, size_t id) {
        this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::insert(val, id);
    }

    void remove(size_t id) {
        this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::remove(id);
    }

    uint64_t rank(uint64_t c, uint64_t i) {
        return this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::rank(c, i);
    }

    uint64_t select(uint64_t c, uint64_t i) {
        return this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::select(c, i);
    }

   // uint64_t const& operator[](size_t id) const {
   //     return this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::operator[](id);
   // }

    uint64_t operator[](size_t id) {
        return this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::operator[](id);
    }

//    uint64_t get(size_t id) {
//        return this->operator[](id);
//    }

//    void set(size_t id, uint64_t val) {
//        this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::operator[](id) = val;
//    }

    bool get_bit_raw(uint64_t id) {
        return this->libmaus2::wavelet::DynamicWaveletTree<6, 64>::R->operator[](id);
    }

};


#endif
