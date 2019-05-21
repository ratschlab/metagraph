#ifndef __BIT_VECTOR_HPP__
#define __BIT_VECTOR_HPP__

#include <cstdint>
#include <mutex>
#include <atomic>

#include <sdsl/wavelet_trees.hpp>
#include <libmaus2/bitbtree/bitbtree.hpp>

#include "bitmap.hpp"


// Bitmap with rank and select operations ( rank/select dictionary )
class bit_vector : public bitmap {
  public:
    virtual ~bit_vector() {};

    // Computes the number of set bits in the subarray indexed by [0,1,...,id]
    virtual uint64_t rank1(uint64_t id) const = 0;
    virtual uint64_t rank0(uint64_t id) const;
    // Returns the i-th set bit, starting from 1
    virtual uint64_t select1(uint64_t i) const = 0;

    virtual uint64_t next1(uint64_t id) const = 0;
    virtual uint64_t prev1(uint64_t id) const = 0;

    virtual void insert_bit(uint64_t id, bool val) = 0;
    virtual void delete_bit(uint64_t id) = 0;
    virtual void set(uint64_t id, bool val) override = 0;
    // Can be redefined, e.g. without rebalancing
    virtual void setBitQuick(uint64_t id, bool val) { set(id, val); };

    virtual bool operator[](uint64_t id) const override = 0;
    virtual uint64_t get_int(uint64_t id, uint32_t width) const override = 0;

    virtual bool load(std::istream &in) = 0;
    virtual void serialize(std::ostream &out) const = 0;

    virtual uint64_t size() const override = 0;
    virtual uint64_t num_set_bits() const override { return rank1(size() - 1); }

    /*
        Convert vector to other types:
            - bit_vector_small
            - bit_vector_smart
            - bit_vector_dyn
            - bit_vector_stat
            - bit_vector_sd
            - bit_vector_rrr<3>
            - bit_vector_rrr<8>
            - bit_vector_rrr<15>
            - bit_vector_rrr<31>
            - bit_vector_rrr<63>
            - bit_vector_rrr<127>
            - bit_vector_rrr<255>
            - sdsl::bit_vector

        FYI: This function invalidates the current object
    */
    template <class Vector>
    Vector convert_to();

    template <class Vector>
    Vector copy_to() const;

    virtual std::unique_ptr<bit_vector> copy() const = 0;

    virtual sdsl::bit_vector to_vector() const;
};

std::ostream& operator<<(std::ostream &os, const bit_vector &bv);


class bit_vector_dyn : public bit_vector {
  public:
    explicit bit_vector_dyn(uint64_t size = 0, bool value = 0);
    explicit bit_vector_dyn(const sdsl::bit_vector &vector);
    bit_vector_dyn(std::initializer_list<bool> init);

    virtual std::unique_ptr<bit_vector> copy() const override;

    uint64_t rank1(uint64_t id) const override;
    uint64_t select1(uint64_t id) const override;
    uint64_t select0(uint64_t id) const;

    uint64_t next1(uint64_t id) const override;
    uint64_t prev1(uint64_t id) const override;

    void set(uint64_t id, bool val) override;
    void setBitQuick(uint64_t id, bool val) override;
    bool operator[](uint64_t id) const override;
    uint64_t get_int(uint64_t id, uint32_t width) const override;

    void insert_bit(uint64_t id, bool val) override;
    void delete_bit(uint64_t id) override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    uint64_t size() const override { return vector_.size(); }
    uint64_t num_set_bits() const override { return vector_.count1(); }

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const override;

  private:
    libmaus2::bitbtree::BitBTree<6, 64> vector_;
};


class bit_vector_stat : public bit_vector {
    friend bit_vector;

  public:
    explicit bit_vector_stat(uint64_t size = 0, bool value = 0);
    explicit bit_vector_stat(const sdsl::bit_vector &vector) noexcept
          : bit_vector_stat(sdsl::bit_vector(vector)) {}
    explicit bit_vector_stat(const bit_vector_stat &other);
    bit_vector_stat(const std::function<void(const VoidCall<uint64_t>&)> &call_ones,
                    uint64_t size);
    bit_vector_stat(sdsl::bit_vector&& vector) noexcept;
    bit_vector_stat(bit_vector_stat&& other) noexcept;
    bit_vector_stat(std::initializer_list<bool> init);

    bit_vector_stat& operator=(const bit_vector_stat &other);
    bit_vector_stat& operator=(bit_vector_stat&& other) noexcept;

    virtual std::unique_ptr<bit_vector> copy() const override;

    uint64_t rank1(uint64_t id) const override;
    uint64_t select1(uint64_t id) const override;

    uint64_t next1(uint64_t id) const override;
    uint64_t prev1(uint64_t id) const override;

    void set(uint64_t id, bool val) override;
    bool operator[](uint64_t id) const override;
    uint64_t get_int(uint64_t id, uint32_t width) const override;

    void insert_bit(uint64_t id, bool val) override;
    void delete_bit(uint64_t id) override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    uint64_t size() const override { return vector_.size(); }
    uint64_t num_set_bits() const override { return num_set_bits_; }

    sdsl::bit_vector to_vector() const override { return vector_; }

    virtual void add_to(sdsl::bit_vector *other) const override;

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const override;

    const sdsl::bit_vector& data() const { return vector_; }

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


class bit_vector_sd : public bit_vector {
  public:
    explicit bit_vector_sd(uint64_t size = 0, bool value = false);
    explicit bit_vector_sd(const sdsl::bit_vector &vector);
    explicit bit_vector_sd(const bit_vector_sd &other);

    bit_vector_sd(bit_vector_sd&& other) noexcept;
    bit_vector_sd(std::initializer_list<bool> init);
    bit_vector_sd(const std::function<void(const VoidCall<uint64_t>&)> &call_ones,
                  uint64_t size,
                  uint64_t num_set_bits);

    bit_vector_sd& operator=(const bit_vector_sd &other);
    bit_vector_sd& operator=(bit_vector_sd&& other) noexcept;

    virtual std::unique_ptr<bit_vector> copy() const override;

    uint64_t rank1(uint64_t id) const override;
    uint64_t select0(uint64_t id) const;
    uint64_t select1(uint64_t id) const override;

    uint64_t next1(uint64_t id) const override;
    uint64_t prev1(uint64_t id) const override;

    void set(uint64_t id, bool val) override;
    bool operator[](uint64_t id) const override;
    uint64_t get_int(uint64_t id, uint32_t width) const override;

    void insert_bit(uint64_t id, bool val) override;
    void delete_bit(uint64_t id) override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    uint64_t size() const override { return vector_.size(); }

    sdsl::bit_vector to_vector() const override;

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const override;

    bool is_inverted() const { return inverted_; }

    const sdsl::sd_vector<>& data() const { return vector_; }

  private:
    bool inverted_;
    sdsl::sd_vector<> vector_;
    sdsl::sd_vector<>::rank_1_type rk1_;
    sdsl::sd_vector<>::select_1_type slct1_;
    sdsl::sd_vector<>::select_0_type slct0_;
};


// default block size: 63
template <size_t log_block_size = 63>
class bit_vector_rrr : public bit_vector {
    static constexpr auto kBlockSize = log_block_size;

  public:
    explicit bit_vector_rrr(uint64_t size = 0, bool value = false);
    explicit bit_vector_rrr(const sdsl::bit_vector &vector);
    explicit bit_vector_rrr(const bit_vector_rrr &other);

    bit_vector_rrr(bit_vector_rrr&& other) noexcept;
    bit_vector_rrr(std::initializer_list<bool> init);

    bit_vector_rrr& operator=(const bit_vector_rrr &other);
    bit_vector_rrr& operator=(bit_vector_rrr&& other) noexcept;

    virtual std::unique_ptr<bit_vector> copy() const override;

    uint64_t rank1(uint64_t id) const override;
    uint64_t select0(uint64_t id) const;
    uint64_t select1(uint64_t id) const override;

    uint64_t next1(uint64_t id) const override;
    uint64_t prev1(uint64_t id) const override;

    void set(uint64_t id, bool val) override;
    bool operator[](uint64_t id) const override;
    uint64_t get_int(uint64_t id, uint32_t width) const override;

    void insert_bit(uint64_t id, bool val) override;
    void delete_bit(uint64_t id) override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    uint64_t size() const override { return vector_.size(); }

    sdsl::bit_vector to_vector() const override;

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const override;

    const sdsl::rrr_vector<kBlockSize>& data() const { return vector_; }

  private:
    typename sdsl::rrr_vector<kBlockSize> vector_;
    typename sdsl::rrr_vector<kBlockSize>::rank_1_type rk1_;
    typename sdsl::rrr_vector<kBlockSize>::select_1_type slct1_;
    typename sdsl::rrr_vector<kBlockSize>::select_0_type slct0_;
};


class bit_vector_adaptive : public bit_vector {
    friend bit_vector;

  public:
    virtual ~bit_vector_adaptive() {};

    virtual uint64_t rank1(uint64_t id) const override final { return vector_->rank1(id); }
    virtual uint64_t select1(uint64_t id) const override final { return vector_->select1(id); }

    virtual uint64_t next1(uint64_t id) const override final { return vector_->next1(id); }
    virtual uint64_t prev1(uint64_t id) const override final { return vector_->prev1(id); }

    virtual void set(uint64_t id, bool val) override final { vector_->set(id, val); }
    virtual bool operator[](uint64_t id) const override final { return (*vector_)[id]; }
    virtual uint64_t get_int(uint64_t id, uint32_t width) const override final {
        return vector_->get_int(id, width);
    }

    virtual void insert_bit(uint64_t id, bool val) override final { vector_->insert_bit(id, val); }
    virtual void delete_bit(uint64_t id) override final { vector_->delete_bit(id); }

    virtual bool load(std::istream &in) override final;
    virtual void serialize(std::ostream &out) const override final;

    virtual uint64_t size() const override final { return vector_->size(); }

    virtual sdsl::bit_vector to_vector() const override final { return vector_->to_vector(); }

    virtual void add_to(sdsl::bit_vector *other) const override final { vector_->add_to(other); }

    virtual void call_ones_in_range(uint64_t begin, uint64_t end,
                                    const VoidCall<uint64_t> &callback) const override final {
        vector_->call_ones_in_range(begin, end, callback);
    }

    enum VectorCode {
    // FUI: don't change order of the variables!
    // Add new ones to the end, if any.
    // Otherwise, serialized vectors will be not loadable
        RRR_VECTOR = 0,
        SD_VECTOR,
        STAT_VECTOR
    };

    typedef VectorCode (*DefineRepresentation)(uint64_t /* size */,
                                               uint64_t /* num_set_bits */);

    VectorCode representation_tag() const;

    const bit_vector& data() const { return *vector_; }

  protected:
    bit_vector_adaptive() {}

    explicit bit_vector_adaptive(const bit_vector_adaptive &other);
    bit_vector_adaptive(bit_vector_adaptive&& other) = default;
    bit_vector_adaptive& operator=(const bit_vector_adaptive &other);
    bit_vector_adaptive& operator=(bit_vector_adaptive&& other) = default;

    std::unique_ptr<bit_vector> vector_;
};


/**
 * static hybrid vector:
 *      the internal representation is defined in constructor.
 */
template <bit_vector_adaptive::DefineRepresentation optimal_representation>
class bit_vector_adaptive_stat : public bit_vector_adaptive {
    friend bit_vector;

  public:
    explicit bit_vector_adaptive_stat(uint64_t size = 0, bool value = false);

    template <class BitVector>
    explicit bit_vector_adaptive_stat(const BitVector &vector)
      : bit_vector_adaptive_stat(static_cast<bit_vector&&>(bit_vector_stat(vector))) {}

    bit_vector_adaptive_stat(const VoidCall<const VoidCall<uint64_t>&> &call_ones,
                             uint64_t size,
                             uint64_t num_set_bits);

    bit_vector_adaptive_stat(std::initializer_list<bool> init)
      : bit_vector_adaptive_stat(sdsl::bit_vector(init)) {}

    virtual std::unique_ptr<bit_vector> copy() const override final;

  private:
    explicit bit_vector_adaptive_stat(const bit_vector &vector);
    bit_vector_adaptive_stat(bit_vector&& vector);
};

/**
 * hybrid vector: the smallest representation
 * combines:
 *    - bit_vector_sd
 *    - bit_vector_rrr
 */
bit_vector_adaptive::VectorCode
smallest_representation(uint64_t size, uint64_t num_set_bits);

typedef bit_vector_adaptive_stat<smallest_representation> bit_vector_small;

/**
 * hybrid vector: a good tradeoff between the speed and size
 * combines:
 *    - bit_vector_sd
 *    - bit_vector_stat
 */
bit_vector_adaptive::VectorCode
smart_representation(uint64_t size, uint64_t num_set_bits);

typedef bit_vector_adaptive_stat<smart_representation> bit_vector_smart;

/**
 * Predict space taken by the vector with given parameters
 * In bits.
 */
template <class BitVector>
uint64_t predict_size(uint64_t size, uint64_t num_set_bits);

#endif // __BIT_VECTOR_HPP__
