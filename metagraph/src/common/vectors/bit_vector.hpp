#ifndef __BIT_VECTOR_HPP__
#define __BIT_VECTOR_HPP__

#include <cstdint>
#include <mutex>
#include <atomic>

#include <sdsl/int_vector.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/hyb_vector.hpp>
#include <sdsl/bit_vector_il.hpp>
#include <dynamic.hpp>

#include "bitmap.hpp"


// Bitmap with rank and select operations ( rank/select dictionary )
class bit_vector : public bitmap {
  public:
    virtual ~bit_vector() {}

    // Computes the number of set bits in the subarray indexed by [0,1,...,id]
    virtual uint64_t rank1(uint64_t id) const = 0;
    virtual uint64_t rank0(uint64_t id) const;
    // Returns the i-th set bit, starting from 1
    virtual uint64_t select1(uint64_t i) const = 0;
    // Query bit and rank
    virtual std::pair<bool, uint64_t> inverse_select(uint64_t id) const {
        return std::make_pair(operator[](id), rank1(id));
    }
    // Query bit and return rank if the bit is set, otherwise return zero
    virtual uint64_t conditional_rank1(uint64_t id) const {
        if (operator[](id)) {
            return rank1(id);
        } else {
            return 0;
        }
    }

    virtual uint64_t next1(uint64_t id) const = 0;
    virtual uint64_t prev1(uint64_t id) const = 0;

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
            - bit_vector_hyb<>
            - bit_vector_il<>
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

    virtual sdsl::bit_vector to_vector() const = 0;

    virtual void add_to(sdsl::bit_vector *other) const override;
};

std::ostream& operator<<(std::ostream &os, const bit_vector &bv);


class bit_vector_dyn : public bit_vector {
  public:
    explicit bit_vector_dyn(uint64_t size = 0, bool value = 0);
    explicit bit_vector_dyn(const sdsl::bit_vector &vector);
    bit_vector_dyn(std::initializer_list<bool> init);

    std::unique_ptr<bit_vector> copy() const override;

    uint64_t rank1(uint64_t id) const override;
    uint64_t select1(uint64_t id) const override;
    uint64_t select0(uint64_t id) const;

    uint64_t next1(uint64_t id) const override;
    uint64_t prev1(uint64_t id) const override;

    bool operator[](uint64_t id) const override;
    uint64_t get_int(uint64_t id, uint32_t width) const override;

    void insert_bit(uint64_t id, bool val);
    void delete_bit(uint64_t id);
    void set(uint64_t id, bool val);

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    uint64_t size() const override { return vector_.size(); }

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const override;

    sdsl::bit_vector to_vector() const override;

  private:
    dyn::suc_bv vector_;
};


class bit_vector_stat : public bit_vector {
    friend bit_vector;

  public:
    explicit bit_vector_stat(uint64_t size = 0, bool value = 0);
    explicit bit_vector_stat(const sdsl::bit_vector &vector) noexcept;
    bit_vector_stat(const sdsl::bit_vector &vector, uint64_t num_set_bits);
    explicit bit_vector_stat(const bit_vector_stat &other);
    bit_vector_stat(const std::function<void(const VoidCall<uint64_t>&)> &call_ones,
                    uint64_t size);
    bit_vector_stat(sdsl::bit_vector&& vector) noexcept;
    bit_vector_stat(sdsl::bit_vector&& vector, uint64_t num_set_bits);
    bit_vector_stat(bit_vector_stat&& other) noexcept;
    bit_vector_stat(std::initializer_list<bool> init);

    bit_vector_stat& operator=(const bit_vector_stat &other);
    bit_vector_stat& operator=(bit_vector_stat&& other) noexcept;

    std::unique_ptr<bit_vector> copy() const override;

    uint64_t rank1(uint64_t id) const override;
    uint64_t select1(uint64_t id) const override;

    uint64_t next1(uint64_t id) const override;
    uint64_t prev1(uint64_t id) const override;

    bool operator[](uint64_t id) const override;
    uint64_t get_int(uint64_t id, uint32_t width) const override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    uint64_t size() const override { return vector_.size(); }
    uint64_t num_set_bits() const override { return num_set_bits_; }

    sdsl::bit_vector to_vector() const override { return vector_; }

    void add_to(sdsl::bit_vector *other) const override;

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const override;

    const sdsl::bit_vector& data() const { return vector_; }

  private:
    void init_rs() const;

    sdsl::bit_vector vector_;
    uint64_t num_set_bits_ = 0;

    // maintain rank/select operations
    mutable sdsl::rank_support_v5<> rk_;
    mutable sdsl::select_support_mcl<> slct_;
};


class bit_vector_sd : public bit_vector {
  public:
    explicit bit_vector_sd(uint64_t size = 0, bool value = false);
    explicit bit_vector_sd(const sdsl::bit_vector &vector);
    explicit bit_vector_sd(const sdsl::bit_vector &vector, uint64_t num_set_bits);
    explicit bit_vector_sd(const bit_vector_sd &other);

    bit_vector_sd(bit_vector_sd&& other) noexcept;
    bit_vector_sd(std::initializer_list<bool> init);
    bit_vector_sd(const std::function<void(const VoidCall<uint64_t>&)> &call_ones,
                  uint64_t size,
                  uint64_t num_set_bits);

    bit_vector_sd& operator=(const bit_vector_sd &other);
    bit_vector_sd& operator=(bit_vector_sd&& other) noexcept;

    std::unique_ptr<bit_vector> copy() const override;

    uint64_t rank1(uint64_t id) const override;
    uint64_t select0(uint64_t id) const;
    uint64_t select1(uint64_t id) const override;

    uint64_t next1(uint64_t id) const override;
    uint64_t prev1(uint64_t id) const override;

    bool operator[](uint64_t id) const override;
    uint64_t get_int(uint64_t id, uint32_t width) const override;

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


template <uint32_t block_rate = 16>
class bit_vector_hyb : public bit_vector {
  public:
    explicit bit_vector_hyb(uint64_t size = 0, bool value = false);
    explicit bit_vector_hyb(const sdsl::bit_vector &vector);
    explicit bit_vector_hyb(const sdsl::bit_vector &vector, uint64_t num_set_bits);
    explicit bit_vector_hyb(const bit_vector_hyb &other);

    bit_vector_hyb(bit_vector_hyb&& other) noexcept;
    bit_vector_hyb(std::initializer_list<bool> init);

    bit_vector_hyb& operator=(const bit_vector_hyb &other);
    bit_vector_hyb& operator=(bit_vector_hyb&& other) noexcept;

    std::unique_ptr<bit_vector> copy() const override;

    uint64_t rank1(uint64_t id) const override;
    uint64_t select0(uint64_t id) const;
    uint64_t select1(uint64_t id) const override;

    uint64_t next1(uint64_t id) const override;
    uint64_t prev1(uint64_t id) const override;

    bool operator[](uint64_t id) const override;
    uint64_t get_int(uint64_t id, uint32_t width) const override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    uint64_t size() const override { return vector_.size(); }

    sdsl::bit_vector to_vector() const override;

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const override;

    const sdsl::hyb_vector<block_rate>& data() const { return vector_; }

  private:
    sdsl::hyb_vector<block_rate> vector_;
    typename sdsl::hyb_vector<block_rate>::rank_1_type rk1_;
    typename sdsl::hyb_vector<block_rate>::select_1_type slct1_;
    typename sdsl::hyb_vector<block_rate>::select_0_type slct0_;
};


template <uint32_t block_size = 512>
class bit_vector_il : public bit_vector {
    static constexpr auto kBlockSize = block_size;

  public:
    explicit bit_vector_il(uint64_t size = 0, bool value = false);
    explicit bit_vector_il(const sdsl::bit_vector &vector);
    explicit bit_vector_il(const sdsl::bit_vector &vector, uint64_t num_set_bits);
    explicit bit_vector_il(const bit_vector_il &other);

    bit_vector_il(bit_vector_il&& other) noexcept;
    bit_vector_il(std::initializer_list<bool> init);

    bit_vector_il& operator=(const bit_vector_il &other);
    bit_vector_il& operator=(bit_vector_il&& other) noexcept;

    std::unique_ptr<bit_vector> copy() const override;

    uint64_t rank1(uint64_t id) const override;
    uint64_t select0(uint64_t id) const;
    uint64_t select1(uint64_t id) const override;

    uint64_t next1(uint64_t id) const override;
    uint64_t prev1(uint64_t id) const override;

    bool operator[](uint64_t id) const override;
    uint64_t get_int(uint64_t id, uint32_t width) const override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    uint64_t size() const override { return vector_.size(); }

    sdsl::bit_vector to_vector() const override;

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const override;

    const sdsl::bit_vector_il<kBlockSize>& data() const { return vector_; }

  private:
    sdsl::bit_vector_il<kBlockSize> vector_;
    typename sdsl::bit_vector_il<kBlockSize>::rank_1_type rk1_;
    typename sdsl::bit_vector_il<kBlockSize>::select_1_type slct1_;
    typename sdsl::bit_vector_il<kBlockSize>::select_0_type slct0_;
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

    std::unique_ptr<bit_vector> copy() const override;

    uint64_t rank1(uint64_t id) const override;
    uint64_t select0(uint64_t id) const;
    uint64_t select1(uint64_t id) const override;
    std::pair<bool, uint64_t> inverse_select(uint64_t id) const override;
    uint64_t conditional_rank1(uint64_t id) const override;

    uint64_t next1(uint64_t id) const override;
    uint64_t prev1(uint64_t id) const override;

    bool operator[](uint64_t id) const override;
    uint64_t get_int(uint64_t id, uint32_t width) const override;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    uint64_t size() const override { return vector_.size(); }

    sdsl::bit_vector to_vector() const override;

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const override;

    const sdsl::rrr_vector<kBlockSize>& data() const { return vector_; }

  private:
    sdsl::rrr_vector<kBlockSize> vector_;
    typename sdsl::rrr_vector<kBlockSize>::rank_1_type rk1_;
    typename sdsl::rrr_vector<kBlockSize>::select_1_type slct1_;
    typename sdsl::rrr_vector<kBlockSize>::select_0_type slct0_;
};


class bit_vector_adaptive : public bit_vector {
    friend bit_vector;

  public:
    virtual ~bit_vector_adaptive() {}

    virtual uint64_t rank1(uint64_t id) const override final { return vector_->rank1(id); }
    virtual uint64_t select1(uint64_t id) const override final { return vector_->select1(id); }
    virtual std::pair<bool, uint64_t> inverse_select(uint64_t id) const override final {
        return vector_->inverse_select(id);
    }
    virtual uint64_t conditional_rank1(uint64_t id) const override final {
        return vector_->conditional_rank1(id);
    }

    virtual uint64_t next1(uint64_t id) const override final { return vector_->next1(id); }
    virtual uint64_t prev1(uint64_t id) const override final { return vector_->prev1(id); }

    virtual bool operator[](uint64_t id) const override final { return (*vector_)[id]; }
    virtual uint64_t get_int(uint64_t id, uint32_t width) const override final {
        return vector_->get_int(id, width);
    }

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

    explicit bit_vector_adaptive_stat(const bit_vector &vector);
    explicit bit_vector_adaptive_stat(const sdsl::bit_vector &vector);

    bit_vector_adaptive_stat(bit_vector&& vector);
    bit_vector_adaptive_stat(sdsl::bit_vector&& vector)
      : bit_vector_adaptive_stat(bit_vector_stat(std::move(vector))) {}

    bit_vector_adaptive_stat(const VoidCall<const VoidCall<uint64_t>&> &call_ones,
                             uint64_t size,
                             uint64_t num_set_bits);

    bit_vector_adaptive_stat(std::initializer_list<bool> init)
      : bit_vector_adaptive_stat(sdsl::bit_vector(init)) {}

    std::unique_ptr<bit_vector> copy() const override final;
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

// indexes are distinct and sorted
sdsl::bit_vector subvector(const bit_vector &col,
                           const std::vector<uint64_t> &indexes);

#endif // __BIT_VECTOR_HPP__
