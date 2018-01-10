#ifndef __BIT_VECTOR_HPP__
#define __BIT_VECTOR_HPP__

#include <cassert>

#define __STDC_FORMAT_MACROS 1
#include <cinttypes>

#include <sdsl/wavelet_trees.hpp>
#include <libmaus2/bitbtree/bitbtree.hpp>


class bit_vector {
  public:
    virtual ~bit_vector() {};

    // Computes the number of set bits in the subarray indexed by [0,1,...,id]
    virtual uint64_t rank1(uint64_t id) const = 0;
    // Returns the i-th set bit, starting from 1
    virtual uint64_t select1(uint64_t i) const = 0;
    virtual uint64_t size() const = 0;
    virtual void set(uint64_t id, bool val) = 0;
    // Can be redefined, e.g. without rebalancing
    virtual void setBitQuick(uint64_t id, bool val) { set(id, val); };
    virtual bool operator[](uint64_t id) const = 0;
    virtual void insertBit(uint64_t id, bool val) = 0;
    virtual void deleteBit(uint64_t id) = 0;
    virtual void serialise(std::ostream &out) const = 0;
    virtual bool deserialise(std::istream &in) = 0;
    virtual std::vector<bool> to_vector() const {
        std::vector<bool> result(size());
        for (uint64_t i = 0; i < size(); ++i) {
            result[i] = operator[](i);
        }
        return result;
    }

    friend inline std::ostream& operator<<(std::ostream &os,
                                           const bit_vector &bv);

  protected:
    virtual void print(std::ostream &os) const = 0;
};

inline std::ostream& operator<<(std::ostream &os, const bit_vector &bv) {
    bv.print(os);
    return os;
}


class bit_vector_dyn : public bit_vector {
  public:
    bit_vector_dyn() : vector_() {}

    bit_vector_dyn(uint64_t size, bool value) : vector_(size, value) {}

    explicit bit_vector_dyn(const std::vector<bool> &v) : vector_(v.size(), false) {
        for (uint64_t i = 0; i < v.size(); ++i) {
            if (v.at(i))
                this->set(i, true);
        }
    }

    explicit bit_vector_dyn(const bit_vector_dyn &v) : vector_(v.vector_) {}

    explicit bit_vector_dyn(const bit_vector &v) : vector_(v.size(), false) {
        for (uint64_t i = 0; i < v.size(); ++i) {
            if (v[i])
                this->set(i, true);
        }
    }

    explicit bit_vector_dyn(std::istream &in) {
        this->deserialise(in);
    }

    uint64_t rank1(uint64_t id) const {
        return vector_.rank1(id);
    }

    uint64_t select1(uint64_t id) const {
        assert(id > 0 && size() > 0 && id <= rank1(size() - 1));
        return vector_.select1(id - 1);
    }

    uint64_t size() const {
        return vector_.size();
    }

    void set(uint64_t id, bool val) {
        vector_.set(id, val);
    }

    void setBitQuick(uint64_t id, bool val) {
        vector_.setBitQuick(id, val);
    }

    bool operator[](uint64_t id) const {
        assert(id < size());
        return vector_[id];
    }

    void insertBit(uint64_t id, bool val) {
        vector_.insertBit(id, val);
    }

    void deleteBit(uint64_t id) {
        assert(size() > 0);
        vector_.deleteBit(id);
    }

    bool deserialise(std::istream &in) {
        vector_.deserialise(in);
        return true;
    }

    void serialise(std::ostream &out) const {
        vector_.serialise(out);
    }

  protected:
    void print(std::ostream& os) const { os << vector_; }

  private:
    libmaus2::bitbtree::BitBTree<6, 64> vector_;
};


class bit_vector_stat : public bit_vector {
  public:
    bit_vector_stat() : vector_() {}

    bit_vector_stat(uint64_t size, bool value) : vector_(size, value) {}

    bit_vector_stat(const bit_vector_stat &other) : vector_(other.vector_) {}

    bit_vector_stat(const bit_vector &other)
            : vector_(other.size(), 0) {
        for (uint64_t i = 0; i < other.size(); ++i) {
            if (other[i])
                vector_[i] = 1;
        }
    }

    bit_vector_stat(const std::vector<bool> &other)
            : vector_(other.size(), 0) {
        for (uint64_t i = 0; i < other.size(); ++i) {
            if (other.at(i))
                vector_[i] = 1;
        }
    }


    bit_vector_stat(std::initializer_list<bool> init) : vector_(init) {}

    bit_vector_stat(std::istream &in) {
        deserialise(in);
    }

    uint64_t rank1(uint64_t id) const {
        if (requires_update_)
            const_cast<bit_vector_stat*>(this)->init_rs();
        //the rank method in SDSL does not include id in the count
        return rk_(id >= this->size() ? this->size() : id + 1);
    }

    uint64_t select1(uint64_t id) const {
        assert(id > 0 && size() > 0 && id <= rank1(size() - 1));
        if (requires_update_)
            const_cast<bit_vector_stat*>(this)->init_rs();
        return slct_(id);
    }

    void set(uint64_t id, bool val) {
        vector_[id] = val;
        requires_update_ = true;
    }

    bool operator[](uint64_t id) const {
        assert(id < size());
        return vector_[id];
    }

    void insertBit(uint64_t id, bool val) {
        vector_.resize(size() + 1);
        if (vector_.size() > 1)
            std::copy_backward(vector_.begin() + id, vector_.end() - 1, vector_.end());
        set(id, val);
        requires_update_ = true;
    }

    void deleteBit(uint64_t id) {
        assert(size() > 0);
        if (vector_.size() > 1)
            std::copy(vector_.begin() + id + 1, vector_.end(), vector_.begin() + id);
        vector_.resize(vector_.size() - 1);
        requires_update_ = true;
    }

    bool deserialise(std::istream &in) {
        vector_.load(in);
        requires_update_ = true;
        return true;
    }

    void serialise(std::ostream &out) const {
        vector_.serialize(out);
    }

    uint64_t size() const {
        return vector_.size();
    }

  protected:
    void print(std::ostream &os) const { os << vector_; }

  private:
    sdsl::bit_vector vector_;

    // maintain rank/select operations
    sdsl::rank_support_v5<> rk_;
    sdsl::select_support_mcl<> slct_;
    bool requires_update_ = true;

    void init_rs() {
        rk_ = sdsl::rank_support_v5<>(&vector_);
        slct_ = sdsl::select_support_mcl<>(&vector_);
        requires_update_ = false;
    }
};


#endif // __BIT_VECTOR_HPP__
