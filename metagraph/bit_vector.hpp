#ifndef __BIT_VECTOR_HPP__
#define __BIT_VECTOR_HPP__

#include <cstdint>

#include <sdsl/wavelet_trees.hpp>
#include <libmaus2/bitbtree/bitbtree.hpp>

class bit_vector {
  public:
    typedef std::vector<bool> BoolVector;

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

    // this invalidates the current object
    template <class Vector>
    Vector convert_to();

    virtual BoolVector to_vector() const;
};

std::ostream& operator<<(std::ostream &os, const bit_vector &bv);

class bit_vector_dyn : public bit_vector {
  public:
    bit_vector_dyn() : vector_() {}
    bit_vector_dyn(uint64_t size, bool value) : vector_(size, value) {}

    template <class BinaryVector>
    explicit bit_vector_dyn(const BinaryVector &v);
    explicit bit_vector_dyn(std::initializer_list<bool> init);

    explicit bit_vector_dyn(const bit_vector_dyn &v);
    bit_vector_dyn(bit_vector_dyn&& v);

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

  private:
    libmaus2::bitbtree::BitBTree<6, 64> vector_;

    // constructor from packed bits
    explicit bit_vector_dyn(const std::vector<uint64_t> &v, size_t num_bits);
};


class bit_vector_stat : public bit_vector {
  public:
    bit_vector_stat() {}
    bit_vector_stat(uint64_t size, bool value);
    explicit bit_vector_stat(const bit_vector_stat &other);
    explicit bit_vector_stat(const std::vector<bool> &other);
    explicit bit_vector_stat(std::initializer_list<bool> init);

    bit_vector_stat(bit_vector_stat&& other);
    bit_vector_stat(sdsl::bit_vector&& vector);
    bit_vector_stat& operator=(sdsl::bit_vector&& vector);

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

  private:
    sdsl::bit_vector vector_;
    uint64_t num_set_bits_ = 0;

    // maintain rank/select operations
    sdsl::rank_support_v5<> rk_;
    sdsl::select_support_mcl<> slct_;
    bool requires_update_ = true;

    void init_rs();
};

class bit_vector_small : public bit_vector {
    typedef sdsl::sd_vector<> bv_type;

  public:
    bit_vector_small(uint64_t size = 0, bool value = false);
    explicit bit_vector_small(const bit_vector_small &other);
    explicit bit_vector_small(const std::vector<bool> &other);
    explicit bit_vector_small(std::initializer_list<bool> init);

    bit_vector_small(bit_vector_small&& vector);
    bit_vector_small(sdsl::bit_vector&& vector);
    bit_vector_small& operator=(sdsl::bit_vector&& vector);

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

  private:
    bit_vector_small(const sdsl::bit_vector &bv, uint64_t num_set_bits);
    sdsl::sd_vector<> vector_;
    uint64_t num_set_bits_ = 0;

    static sdsl::bit_vector invert(const sdsl::bit_vector &other, uint64_t num_set_bits);

    // maintain rank/select operations
    typename decltype(vector_)::rank_0_type rk_;
    typename decltype(vector_)::select_0_type slct_;
};


#endif // __BIT_VECTOR_HPP__
