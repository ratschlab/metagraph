#ifndef __BIT_VECTOR_HPP__
#define __BIT_VECTOR_HPP__

#include <cstdint>
#include <mutex>
#include <atomic>

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
    virtual bool load(std::istream &in) = 0;
    virtual void serialize(std::ostream &out) const = 0;
    virtual uint64_t size() const = 0;
    virtual uint64_t num_set_bits() const { return rank1(size() - 1); }

    // FYI: This function invalidates the current object
    template <class Vector>
    Vector convert_to();

    virtual std::vector<bool> to_vector() const;
};

std::ostream& operator<<(std::ostream &os, const bit_vector &bv);


class bit_vector_dyn : public bit_vector {
  public:
    explicit bit_vector_dyn(uint64_t size = 0, bool value = 0);
    template <class BitVector>
    explicit bit_vector_dyn(const BitVector &vector);
    bit_vector_dyn(std::initializer_list<bool> init);

    uint64_t rank1(uint64_t id) const;
    uint64_t select1(uint64_t id) const;

    void set(uint64_t id, bool val);
    void setBitQuick(uint64_t id, bool val) override;
    bool operator[](uint64_t id) const;

    void insertBit(uint64_t id, bool val);
    void deleteBit(uint64_t id);

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    uint64_t size() const { return vector_.size(); }

  private:
    bit_vector_dyn(const std::vector<uint64_t> &bits_packed, size_t num_bits);

    libmaus2::bitbtree::BitBTree<6, 64> vector_;
};


class bit_vector_stat : public bit_vector {
    friend bit_vector;

  public:
    explicit bit_vector_stat(uint64_t size = 0, bool value = 0);
    explicit bit_vector_stat(const std::vector<bool> &vector);
    explicit bit_vector_stat(const bit_vector_stat &other);
    bit_vector_stat(sdsl::bit_vector&& vector);
    bit_vector_stat(bit_vector_stat&& other);
    bit_vector_stat(std::initializer_list<bool> init);

    bit_vector_stat& operator=(const bit_vector_stat &other);
    bit_vector_stat& operator=(bit_vector_stat&& other);

    uint64_t rank1(uint64_t id) const;
    uint64_t select1(uint64_t id) const;

    void set(uint64_t id, bool val);
    bool operator[](uint64_t id) const;

    void insertBit(uint64_t id, bool val);
    void deleteBit(uint64_t id);

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    uint64_t size() const { return vector_.size(); }
    uint64_t num_set_bits() const override { return num_set_bits_; }

  private:
    void init_rs();

    sdsl::bit_vector vector_;
    uint64_t num_set_bits_ = 0;

    // maintain rank/select operations
    sdsl::rank_support_v5<> rk_;
    sdsl::select_support_mcl<> slct_;
    std::atomic_bool requires_update_ { true };
    std::mutex mu_;
};


class bit_vector_small : public bit_vector {
  public:
    explicit bit_vector_small(uint64_t size = 0, bool value = false);
    template <class BitVector>
    explicit bit_vector_small(const BitVector &vector);
    explicit bit_vector_small(const bit_vector_small &other);
    bit_vector_small(bit_vector_small&& other);
    bit_vector_small(std::initializer_list<bool> init);

    bit_vector_small& operator=(const bit_vector_small &other);
    bit_vector_small& operator=(bit_vector_small&& other);

    uint64_t rank1(uint64_t id) const;
    uint64_t select1(uint64_t id) const;

    void set(uint64_t id, bool val);
    bool operator[](uint64_t id) const;

    void insertBit(uint64_t id, bool val);
    void deleteBit(uint64_t id);

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    uint64_t size() const { return vector_.size(); }

    std::vector<bool> to_vector() const override;
    void add_to(std::vector<bool> *other) const;

    bool is_inverted() const { return inverted_; }

  private:
    bool inverted_;
    sdsl::sd_vector<> vector_;
    sdsl::sd_vector<>::rank_1_type rk1_;
    sdsl::sd_vector<>::select_1_type slct1_;
    sdsl::sd_vector<>::select_0_type slct0_;
};


#endif // __BIT_VECTOR_HPP__
