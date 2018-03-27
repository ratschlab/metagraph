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
    virtual void set(uint64_t id, bool val) = 0;
    // Can be redefined, e.g. without rebalancing
    virtual void setBitQuick(uint64_t id, bool val) { set(id, val); };
    virtual bool operator[](uint64_t id) const = 0;
    virtual void insertBit(uint64_t id, bool val) = 0;
    virtual void deleteBit(uint64_t id) = 0;
    virtual void serialise(std::ostream &out) const = 0;
    virtual bool deserialise(std::istream &in) = 0;
    virtual uint64_t size() const = 0;
    virtual uint64_t get_num_set_bits() const { return rank1(size() - 1); }

    virtual std::vector<bool> to_vector() const;
    virtual void print(std::ostream &os) const = 0;
};

std::ostream& operator<<(std::ostream &os, const bit_vector &bv);


class bit_vector_dyn : public bit_vector {
  public:
    bit_vector_dyn() : vector_() {}
    bit_vector_dyn(uint64_t size, bool value) : vector_(size, value) {}
    explicit bit_vector_dyn(const std::vector<uint64_t> &v, size_t num_bits);
    explicit bit_vector_dyn(const std::vector<bool> &v);
    explicit bit_vector_dyn(const bit_vector_dyn &v);
    explicit bit_vector_dyn(const bit_vector &v);
    explicit bit_vector_dyn(std::istream &in);

    uint64_t rank1(uint64_t id) const;
    uint64_t select1(uint64_t id) const;

    void set(uint64_t id, bool val);
    void setBitQuick(uint64_t id, bool val) override;
    bool operator[](uint64_t id) const;

    void insertBit(uint64_t id, bool val);
    void deleteBit(uint64_t id);

    bool deserialise(std::istream &in);
    void serialise(std::ostream &out) const;

    uint64_t size() const { return vector_.size(); }

    void print(std::ostream& os) const { os << vector_; }

  private:
    libmaus2::bitbtree::BitBTree<6, 64> vector_;
};


class bit_vector_stat : public bit_vector {
  public:
    bit_vector_stat() {}
    bit_vector_stat(uint64_t size, bool value);
    explicit bit_vector_stat(const bit_vector_stat &other);
    explicit bit_vector_stat(const bit_vector &other);
    explicit bit_vector_stat(const std::vector<bool> &other);
    explicit bit_vector_stat(std::initializer_list<bool> init);
    explicit bit_vector_stat(std::istream &in);

    uint64_t rank1(uint64_t id) const;
    uint64_t select1(uint64_t id) const;

    void set(uint64_t id, bool val);
    bool operator[](uint64_t id) const;

    void insertBit(uint64_t id, bool val);
    void deleteBit(uint64_t id);

    bool deserialise(std::istream &in);
    void serialise(std::ostream &out) const;

    uint64_t size() const { return vector_.size(); }
    uint64_t get_num_set_bits() const override { return num_set_bits_; }

  protected:
    void print(std::ostream &os) const { os << vector_; }

  private:
    sdsl::bit_vector vector_;
    uint64_t num_set_bits_ = 0;

    // maintain rank/select operations
    sdsl::rank_support_v5<> rk_;
    sdsl::select_support_mcl<> slct_;
    bool requires_update_ = true;

    void init_rs();
};


#endif // __BIT_VECTOR_HPP__
