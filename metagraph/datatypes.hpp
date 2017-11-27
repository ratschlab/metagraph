#ifndef __DATATYPES_HPP__
#define __DATATYPES_HPP__

#include <pthread.h>
#include <assert.h>
#include <functional>
#include <set>
#include <unordered_map>
#include <unordered_set>

#define __STDC_FORMAT_MACROS 1
#include <inttypes.h>

#include <sdsl/wavelet_trees.hpp>
#include <libmaus2/bitbtree/bitbtree.hpp>
#include <libmaus2/wavelet/DynamicWaveletTree.hpp>

#include "kseq.h"


struct HitInfo {
    uint64_t rl;
    uint64_t ru;
    uint64_t str_pos;
    uint64_t graph_pos;
    uint64_t distance;
    std::string cigar;
    std::vector<uint64_t> path;
};


class HitInfoCompare {
  public:
    explicit HitInfoCompare(bool is_reverse) : is_reverse_(is_reverse) {}
    HitInfoCompare() : HitInfoCompare::HitInfoCompare(false) {}

    bool operator()(const HitInfo &lhs, const HitInfo &rhs) const {
        return is_reverse_ ? lhs.distance < rhs.distance
                           : lhs.distance > rhs.distance;
    }
  private:
    bool is_reverse_;
};


class bit_vector {
  public:
    virtual ~bit_vector() {};

    virtual uint64_t size() const = 0;
    virtual void set(size_t id, bool val) = 0;
    // Can be redefined, e.g. without rebalancing
    virtual void setBitQuick(size_t id, bool val) { set(id, val); };
    virtual bool operator[](size_t id) const = 0;
    virtual void insertBit(size_t id, bool val) = 0;
    virtual void deleteBit(size_t id) = 0;
    virtual void serialise(std::ostream &out) const = 0;
    virtual void deserialise(std::istream &in) = 0;
    virtual uint64_t select1(size_t id) const = 0;
    virtual uint64_t rank1(size_t id) const = 0;
    friend inline std::ostream& operator<<(std::ostream &os, const bit_vector& bv);
  protected:
    virtual void print(std::ostream& os) const = 0;
};

inline std::ostream& operator<<(std::ostream& os, const bit_vector& bv) {
    bv.print(os);
    return os;
}

class bit_vector_dyn : public bit_vector {
  public:
    bit_vector_dyn() : vector_() {}

    bit_vector_dyn(size_t size, bool def) : vector_(size, def) {}

    explicit bit_vector_dyn(std::vector<bool> &v) : vector_(v.size(), false) {
        for (size_t i = 0; i < v.size(); ++i) {
            if (v.at(i))
                this->set(i, true);
        }
    }

    explicit bit_vector_dyn(const bit_vector_dyn &v) : vector_(v.vector_) {}

    explicit bit_vector_dyn(const bit_vector &v) : vector_(v.size(), false) {
        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i])
                this->set(i, true);
        }
    }

    explicit bit_vector_dyn(std::istream &in) {
        this->deserialise(in);
    }

    uint64_t size() const {
        return vector_.size();
    }

    void set(size_t id, bool val) {
        vector_.set(id, val);
    }

    void setBitQuick(size_t id, bool val) {
        vector_.setBitQuick(id, val);
    }

    bool operator[](size_t id) const {
        return vector_[id];
    }

    void insertBit(size_t id, bool val) {
        vector_.insertBit(id, val);
    }

    void deleteBit(size_t id) {
        vector_.deleteBit(id);
    }

    void deserialise(std::istream &in) {
        vector_.deserialise(in);
    }

    void serialise(std::ostream &out) const {
        vector_.serialise(out);
    }

    uint64_t select1(size_t id) const {
        return vector_.select1(id);
    }

    uint64_t rank1(size_t id) const {
        return vector_.rank1(id);
    }

  protected:
    void print(std::ostream& os) const {
        os << static_cast<libmaus2::bitbtree::BitBTree<6, 64> >(vector_);
    }

  private:
    libmaus2::bitbtree::BitBTree<6, 64> vector_;
};


class bit_vector_stat : public bit_vector {
  public:
    bit_vector_stat() : vector_() {}

    bit_vector_stat(size_t size, bool def) : vector_(size, def) {}

    bit_vector_stat(const bit_vector_stat &other) : vector_(other.vector_) {}

    bit_vector_stat(const bit_vector &V) : vector_(V.size(), 0) {
        for (size_t i = 0; i < V.size(); ++i) {
            if (V[i])
                vector_[i] = 1;
        }
    }

    bit_vector_stat(std::istream &in) {
        deserialise(in);
    }

    void set(size_t id, bool val) {
        vector_[id] = val;
        requires_update_ = true;
    }

    bool operator[](size_t id) const {
        return vector_[id];
    }

    void insertBit(size_t id, bool val) {
        vector_.resize(size() + 1);
        if (vector_.size() > 1)
            std::copy_backward(vector_.begin() + id, vector_.end() - 1, vector_.end());
        set(id, val);
        requires_update_ = true;
    }
    void deleteBit(size_t id) {
        if (vector_.size() > 1)
            std::copy(vector_.begin() + id + 1, vector_.end(), vector_.begin() + id);
        vector_.resize(vector_.size() - 1);
        requires_update_ = true;
    }
    void deserialise(std::istream &in) {
        vector_.load(in);
        requires_update_ = true;
    }
    void serialise(std::ostream &out) const {
        vector_.serialize(out);
    }
    uint64_t select1(size_t id) const {
        if (requires_update_)
            const_cast<bit_vector_stat *>(this)->init_rs();
        return slct(id + 1);
    }
    uint64_t rank1(size_t id) const {
        if (requires_update_)
            const_cast<bit_vector_stat *>(this)->init_rs();
        //the rank method in SDSL does not include id in the count
        return rk(id >= this->size() ? this->size() : id + 1);
    }
    uint64_t size() const {
        return vector_.size();
    }
    
  protected:
    void print(std::ostream& os) const {
        os << static_cast<sdsl::bit_vector>(vector_);
    }


  private:
    sdsl::bit_vector vector_;

    // maintain rank/select operations
    sdsl::rank_support_v5<> rk;
    sdsl::select_support_mcl<> slct;
    bool requires_update_ = true;

    void init_rs() {
        rk = sdsl::rank_support_v5<>(&vector_);
        slct = sdsl::select_support_mcl<>(&vector_);
        requires_update_ = false;
    }
};


class wavelet_tree {
  public:
    virtual ~wavelet_tree() {};

    virtual uint64_t size() const = 0;
    //virtual void deserialise(std::istream &in) = 0;
    virtual void serialise(std::ostream &out) const = 0;
    virtual void insert(uint64_t val, size_t id) = 0;
    virtual void remove(size_t id) = 0;
    //virtual void set(size_t id, uint64_t val) = 0;
    virtual uint64_t rank(uint64_t c, uint64_t i) const = 0;
    virtual uint64_t select(uint64_t c, uint64_t i) const = 0;
    //virtual uint64_t const& operator[](size_t id) const = 0;
    virtual uint64_t operator[](size_t id) const = 0;
    //virtual uint64_t get(size_t id) = 0;
    // this only makes sense when implemented in the dynamic part
    virtual bool get_bit_raw(uint64_t id) const = 0;
    friend inline std::ostream& operator<<(std::ostream &os, const wavelet_tree& wt);
  protected:
    virtual void print(std::ostream& os) const = 0;
};

inline std::ostream& operator<<(std::ostream &os, const wavelet_tree& wt) {
    wt.print(os);
    return os;
}


class wavelet_tree_stat : public wavelet_tree {
  public:
    wavelet_tree_stat(size_t logsigma_, size_t size, uint64_t def)
            : int_vector_(2 * size + 1, def, 1 << logsigma_), n_(size) {}

    wavelet_tree_stat(wavelet_tree *T, size_t logsigma_)
            : int_vector_(T->size(), 0, logsigma_), n_(T->size()) {
        for (size_t i = 0; i < n_; ++i) {
            int_vector_[i] = T->operator[](i);
        }
    }

    wavelet_tree_stat(size_t logsigma_)
            : wavelet_tree_stat(logsigma_, 0, 0) {}

    wavelet_tree_stat(std::istream &in) { deserialise(in); }

    uint64_t size() const { return n_; }

    void deserialise(std::istream &in) {
        wwt_.load(in);
        int_vector_.load(in);
        n_ = int_vector_.size();
        // we assume the wwt_ has been build before serialization
        requires_update_ = false;
    }

    void serialise(std::ostream &out) const {
        if (requires_update_)
            const_cast<wavelet_tree_stat*>(this)->init_wt();
        int_vector_.serialize(out);
        wwt_.serialize(out);
    }

    void insert(uint64_t val, size_t id) {
        if (n_ == size()) {
            int_vector_.resize(2 * n_ + 1);
        }
        n_++;
        if (size() > 1)
            std::copy_backward(int_vector_.begin() + id,
                               int_vector_.begin() + n_ - 1,
                               int_vector_.begin() + n_);
        int_vector_[id] = val;
        requires_update_ = true;
    }

    void remove(size_t id) {
        if (this->size() > 1)
            std::copy(int_vector_.begin() + id + 1,
                      int_vector_.begin() + n_,
                      int_vector_.begin() + id);
        n_--;
        requires_update_ = true;
    }

    uint64_t rank(uint64_t c, uint64_t i) const {
        if (requires_update_)
            const_cast<wavelet_tree_stat*>(this)->init_wt();
        return wwt_.rank(i >= wwt_.size() ? wwt_.size() : i+1, c);
    }

    uint64_t select(uint64_t c, uint64_t i) const {

        if (requires_update_)
            const_cast<wavelet_tree_stat*>(this)->init_wt();
        i++;
        uint64_t maxrank = wwt_.rank(wwt_.size(), c);
        if (i > maxrank) {
            //TODO: should this line ever be reached?
            return wwt_.size();
        }
        return wwt_.select(i, c);
    }

    /*sdsl::int_vector<>::reference operator[](const size_t id) {
        requires_update_ = true;
        return sdsl::int_vector<>::operator[](id);
    }*/

    //uint64_t const& operator[](size_t id) const {
    //    return this->operator[](id);
    //}

    uint64_t operator[](size_t id) const {
//        requires_update_ = true;
        return int_vector_[id];
    }

/*    void set(size_t id, uint64_t val) {
        requires_update_ = true;
        int_vector_.operator[](id) = val;
    }

    uint64_t get(size_t id) {
        return int_vector_.operator[](id);
    }
    */

    bool get_bit_raw(uint64_t id) const {
       throw std::logic_error("Not Implemented");
       return id;
    }

  protected:
    void print(std::ostream& os) const {
        os << wwt_;
    }

  private:
    sdsl::int_vector<> int_vector_;
    sdsl::wt_int<> wwt_;
    size_t logsigma_;
    bool requires_update_ = true;
    size_t n_;

    void init_wt() {
        std::cout << "Initializing WT index" << std::endl;
        int_vector_.resize(n_);
        sdsl::construct_im(wwt_, &int_vector_);
        requires_update_ = false;
    }
};


//Child of DynamicWaveletTree allowing for construction from an int vector
class wavelet_tree_dyn : public wavelet_tree {
  public:
    wavelet_tree_dyn(size_t b)
        : wavelet_tree_(b) {
    }

    wavelet_tree_dyn(std::istream &in)
        : wavelet_tree_(in) {
    }

    wavelet_tree_dyn(const std::vector<uint8_t> &W_stat, size_t b, unsigned int parallel=1)
        : wavelet_tree_(initialize_tree(W_stat, b, parallel), b, W_stat.size()) {
    }

    wavelet_tree_dyn(const wavelet_tree &W_stat, size_t b, unsigned int parallel=1)
        : wavelet_tree_(initialize_tree(W_stat, b, parallel), b, W_stat.size()) {
    }

    //wavelet_tree_dyn(std::vector<uint8_t> &W_stat, size_t b, unsigned int parallel=1)
    //    : wavelet_tree_(dynamic_cast<libmaus2::bitbtree::BitBTree<6, 64>*>(makeTree(W_stat, b, parallel)), b, W_stat.size()) {
    //}


    wavelet_tree_dyn(libmaus2::bitbtree::BitBTree<6, 64> *bt, size_t b, size_t n)
        : wavelet_tree_(bt, b, n) {
    }

    //wavelet_tree_dyn(bit_vector* bt, size_t b, size_t n)
    //    : wavelet_tree_(dynamic_cast<libmaus2::bitbtree::BitBTree<6, 64>*>(bt), b, n) {
    //}

    uint64_t size() const {
        return wavelet_tree_.size();
    }

    //void deserialise(std::istream &in) {
    //    wavelet_tree_.deserialise(in);
    //}

    void serialise(std::ostream &out) const {
        wavelet_tree_.serialise(out);
    }

    void insert(uint64_t val, size_t id) {
        wavelet_tree_.insert(val, id);
    }

    void remove(size_t id) {
        wavelet_tree_.remove(id);
    }

    uint64_t rank(uint64_t c, uint64_t i) const {
        return wavelet_tree_.rank(c, i);
    }

    uint64_t select(uint64_t c, uint64_t i) const {
        return wavelet_tree_.select(c, i);
    }

   // uint64_t const& operator[](size_t id) const {
   //     return wavelet_tree_[id];
   // }

    uint64_t operator[](size_t id) const {
        return wavelet_tree_[id];
    }

//    uint64_t get(size_t id) {
//        return this->operator[](id);
//    }

//    void set(size_t id, uint64_t val) {
//        wavelet_tree_[id] = val;
//    }

    bool get_bit_raw(uint64_t id) const {
        return wavelet_tree_.R->operator[](id);
    }
  
  protected:
    void print(std::ostream& os) const {
        os << wavelet_tree_;
    }


  private:
    libmaus2::wavelet::DynamicWaveletTree<6, 64> wavelet_tree_;

    template <class Vector>
    libmaus2::bitbtree::BitBTree<6, 64>* initialize_tree(const Vector &W_stat, size_t const b, unsigned int parallel=1) {
        size_t const n = W_stat.size();

        // compute total offsets for the individual bins
        std::vector<uint64_t> offsets((1ull << (b-1)) - 1, 0);
        #pragma omp parallel num_threads(parallel)
        {
            uint64_t v, m, o;
            #pragma omp for
            for (size_t i = 0; i < n; ++i) {
                m = (1ull << (b - 1));
                v = static_cast<uint64_t>(W_stat[i]);
                o = 0;
                for (uint64_t ib = 1; ib < b; ++ib) {
                    bool const bit = m & v;
                    if (!bit) {
                        #pragma omp critical
                        offsets.at(o) += 1;
                    }
                    o = 2 * o + 1 + bit;
                    m >>= 1;
                }
            }
        }

        auto *tmp = new libmaus2::bitbtree::BitBTree<6, 64>(n * b, false);
        std::vector<uint64_t> upto_offsets((1ull << (b - 1)) - 1, 0);

        #pragma omp parallel num_threads(parallel)
        {
            uint64_t m, v, o, p, co;
            bool bit;
            #pragma omp for
            for (uint64_t i = 0; i < n; ++i) {
                m = (1ull << (b - 1));
                v = (uint64_t) W_stat[i];
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
                    o = 2 * o + 1 + bit;
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
};

#endif // __DATATYPES_HPP__
