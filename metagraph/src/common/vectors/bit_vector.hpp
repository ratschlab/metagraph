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

#include "vector_algorithm.hpp"
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


template <class bv_type,
          class rank_1_type,
          class select_1_type,
          class select_0_type>
class bit_vector_sdsl : public bit_vector {
    template<typename>
    struct is_rrr : std::false_type {};
    template<uint16_t t_bs, class t_rac, uint16_t t_k>
    struct is_rrr<sdsl::rrr_vector<t_bs, t_rac, t_k>> : std::true_type {};

    template<typename>
    struct bv_traits {};

    template<uint32_t k_sblock_rate>
    struct bv_traits<sdsl::hyb_vector<k_sblock_rate>> {
        // hyb_vector doesn't support select
        static constexpr size_t MAX_ITER_BIT_VECTOR = 1000;
        static constexpr size_t SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR = 1000;
    };
    template<uint32_t t_bs>
    struct bv_traits<sdsl::bit_vector_il<t_bs>> {
        static constexpr size_t MAX_ITER_BIT_VECTOR = 1000;
        static constexpr size_t SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR = 1000;
    };
    template<uint16_t t_bs, class t_rac, uint16_t t_k>
    struct bv_traits<sdsl::rrr_vector<t_bs, t_rac, t_k>> {
        static constexpr size_t MAX_ITER_BIT_VECTOR = 1;
        // sequential word access for bit_vector_rrr per bit is 21 times faster than select
        static constexpr size_t SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR = 21;
    };

  public:
    explicit bit_vector_sdsl(uint64_t size = 0, bool value = false)
      : bit_vector_sdsl(sdsl::bit_vector(size, value)) {}
    explicit bit_vector_sdsl(const sdsl::bit_vector &vector)
      : vector_(vector), rk1_(&vector_), slct1_(&vector_), slct0_(&vector_) {}
    explicit bit_vector_sdsl(const bit_vector_sdsl &other) { *this = other; }

    bit_vector_sdsl(sdsl::bit_vector&& other) noexcept
      : vector_(std::move(other)), rk1_(&vector_), slct1_(&vector_), slct0_(&vector_) {}
    bit_vector_sdsl(bit_vector_sdsl&& other) noexcept { *this = std::move(other); }
    bit_vector_sdsl(std::initializer_list<bool> init)
      : bit_vector_sdsl(sdsl::bit_vector(init)) {}

    bit_vector_sdsl& operator=(const bit_vector_sdsl &other) {
        vector_ = other.vector_;
        rk1_ = other.rk1_;
        rk1_.set_vector(&vector_);
        slct1_ = other.slct1_;
        slct1_.set_vector(&vector_);
        slct0_ = other.slct0_;
        slct0_.set_vector(&vector_);
        return *this;
    }
    bit_vector_sdsl& operator=(bit_vector_sdsl&& other) noexcept {
        vector_ = std::move(other.vector_);
        rk1_ = std::move(other.rk1_);
        rk1_.set_vector(&vector_);
        slct1_ = std::move(other.slct1_);
        slct1_.set_vector(&vector_);
        slct0_ = std::move(other.slct0_);
        slct0_.set_vector(&vector_);
        return *this;
    }

    std::unique_ptr<bit_vector> copy() const override {
        return std::make_unique<bit_vector_sdsl>(*this);
    }

    uint64_t rank1(uint64_t id) const override {
        //the rank method in SDSL does not include id in the count
        return rk1_(id >= this->size() ? this->size() : id + 1);
    }
    uint64_t select0(uint64_t id) const {
        assert(id > 0 && size() > 0 && id <= size() - num_set_bits());
        assert(num_set_bits() == rank1(size() - 1));
        return slct0_(id);
    }
    uint64_t select1(uint64_t id) const override {
        assert(id > 0 && size() > 0 && id <= num_set_bits());
        assert(num_set_bits() == rank1(size() - 1));
        return slct1_(id);
    }
    std::pair<bool, uint64_t> inverse_select(uint64_t id) const override {
        if constexpr(is_rrr<bv_type>{}) {
            if constexpr(bv_type::block_size != 15) {
                std::pair<bool, uint64_t> pair = vector_.inverse_select(id);
                pair.second += pair.first;
                return pair;
            }
        }
        // TODO: implement inverse_select for other bit vectors as well
        return bit_vector::inverse_select(id);
    }
    uint64_t conditional_rank1(uint64_t id) const override {
        if constexpr(is_rrr<bv_type>{}) {
            if constexpr(bv_type::block_size != 15) {
                return vector_.conditional_rank(id);
            }
        }
        return bit_vector::conditional_rank1(id);
    }

    uint64_t next1(uint64_t id) const override {
        assert(id < size());
        return ::next1(*this, id, bv_traits<bv_type>::MAX_ITER_BIT_VECTOR);
    }
    uint64_t prev1(uint64_t id) const override {
        assert(id < size());
        return ::prev1(*this, id, bv_traits<bv_type>::MAX_ITER_BIT_VECTOR);
    }

    bool operator[](uint64_t id) const override {
        assert(id < size());
        return vector_[id];
    }
    uint64_t get_int(uint64_t id, uint32_t width) const override {
        return vector_.get_int(id, width);
    }

    bool load(std::istream &in) override {
        if (!in.good())
            return false;

        try {
            vector_.load(in);
            if (!in.good())
                return false;
            rk1_ = decltype(rk1_)(&vector_);
            slct1_ = decltype(slct1_)(&vector_);
            slct0_ = decltype(slct0_)(&vector_);
            return true;
        } catch (const std::bad_alloc &exception) {
            std::cerr << "ERROR: Not enough memory to load "
                      << typeid(bv_type).name() << std::endl;
            return false;
        } catch (...) {
            return false;
        }
    }
    void serialize(std::ostream &out) const override {
        vector_.serialize(out);

        if (!out.good())
            throw std::ofstream::failure("Error when dumping bit_vector_rrr");
    }

    uint64_t size() const override { return vector_.size(); }

    sdsl::bit_vector to_vector() const override {
        if (num_set_bits()
                < size() / bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR) {
            // sparse
            sdsl::bit_vector vector(size(), 0);
            uint64_t max_rank = size() ? rank1(size() - 1) : 0;
            for (uint64_t i = 1; i <= max_rank; ++i) {
                vector[slct1_(i)] = 1;
            }
            return vector;

        } else if ((size() - num_set_bits())
                < size() / bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR) {
            // dense
            sdsl::bit_vector vector(size(), 1);
            uint64_t max_rank = size() ? rank0(size() - 1) : 0;
            for (uint64_t i = 1; i <= max_rank; ++i) {
                vector[slct0_(i)] = 0;
            }
            return vector;

        } else {
            // moderate density
            sdsl::bit_vector vector(size());

            uint64_t i;
            for (i = 0; i + 64 <= vector.size(); i += 64) {
                vector.set_int(i, vector_.get_int(i));
            }
            if (i < vector.size())
                vector.set_int(i, vector_.get_int(i, vector.size() - i), vector.size() - i);

            return vector;
        }
    }

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const override {
        assert(begin <= end);
        assert(end <= size());

        if (num_set_bits()
                < size() / bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR) {
            // sparse
            uint64_t num_ones = end ? rank1(end - 1) : 0;
            for (uint64_t r = begin ? rank1(begin - 1) + 1 : 1; r <= num_ones; ++r) {
                callback(select1(r));
            }
        } else if ((size() - num_set_bits())
                    < size() / bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR) {
            // dense
            uint64_t one_pos = 0;
            uint64_t zero_pos = 0;
            uint64_t num_zeros = end ? rank0(end - 1) : 0;
            for (uint64_t r = begin ? rank0(begin - 1) + 1 : 1; r <= num_zeros; ++r) {
                zero_pos = select0(r);
                while (one_pos < zero_pos) {
                    callback(one_pos++);
                }
                one_pos++;
            }
            while (one_pos < end) {
                callback(one_pos++);
            }
        } else {
            // moderate density
            assert(begin <= end);
            assert(end <= size());
            ::call_ones(vector_, begin, end, callback);
        }
    }

    const bv_type& data() const { return vector_; }

  private:
    bv_type vector_;
    rank_1_type rk1_;
    select_1_type slct1_;
    select_0_type slct0_;
};

template <uint32_t k_sblock_rate = 16>
using bit_vector_hyb
    = bit_vector_sdsl<sdsl::hyb_vector<k_sblock_rate>,
                      typename sdsl::hyb_vector<k_sblock_rate>::rank_1_type,
                      typename sdsl::hyb_vector<k_sblock_rate>::select_1_type,
                      typename sdsl::hyb_vector<k_sblock_rate>::select_0_type>;

template <uint32_t block_size = 512>
using bit_vector_il
    = bit_vector_sdsl<sdsl::bit_vector_il<block_size>,
                      typename sdsl::bit_vector_il<block_size>::rank_1_type,
                      typename sdsl::bit_vector_il<block_size>::select_1_type,
                      typename sdsl::bit_vector_il<block_size>::select_0_type>;

template <uint32_t block_size = 63>
using bit_vector_rrr
    = bit_vector_sdsl<sdsl::rrr_vector<block_size>,
                      typename sdsl::rrr_vector<block_size>::rank_1_type,
                      typename sdsl::rrr_vector<block_size>::select_1_type,
                      typename sdsl::rrr_vector<block_size>::select_0_type>;


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
