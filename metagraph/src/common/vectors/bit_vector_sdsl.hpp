#ifndef __BIT_VECTOR_SDSL_HPP__
#define __BIT_VECTOR_SDSL_HPP__

#include <cmath>
#include <cstdint>

#include <sdsl/int_vector.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/hyb_vector.hpp>
#include <sdsl/bit_vector_il.hpp>

#include "common/serialization.hpp"
#include "vector_algorithm.hpp"
#include "bit_vector.hpp"


template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
class bit_vector_sdsl : public bit_vector {
    friend bit_vector;

    template <typename>
    struct is_rrr : std::false_type {};
    template <uint16_t t_bs, class t_rac, uint16_t t_k>
    struct is_rrr<sdsl::rrr_vector<t_bs, t_rac, t_k>> : std::true_type {};

    template <typename>
    struct bv_traits;

  public:
    explicit bit_vector_sdsl(uint64_t size = 0, bool value = false)
      : bit_vector_sdsl(sdsl::bit_vector(size, value)) {}
    explicit bit_vector_sdsl(const sdsl::bit_vector &vector)
      : vector_(vector), rk1_(&vector_), slct1_(&vector_), slct0_(&vector_),
        num_set_bits_(rk1_(vector_.size())) {}
    explicit bit_vector_sdsl(const bit_vector_sdsl &other) { *this = other; }

    bit_vector_sdsl(sdsl::bit_vector&& vector)
      : vector_(std::move(vector)), rk1_(&vector_), slct1_(&vector_), slct0_(&vector_),
        num_set_bits_(rk1_(vector_.size())) {}
    bit_vector_sdsl(bit_vector_sdsl&& other) { *this = std::move(other); }
    bit_vector_sdsl(std::initializer_list<bool> init)
      : bit_vector_sdsl(sdsl::bit_vector(init)) {}

    inline bit_vector_sdsl& operator=(const bit_vector_sdsl &other);
    inline bit_vector_sdsl& operator=(bit_vector_sdsl&& other) noexcept;

    std::unique_ptr<bit_vector> copy() const override {
        return std::make_unique<bit_vector_sdsl>(*this);
    }

    inline uint64_t rank1(uint64_t id) const override;
    inline uint64_t select0(uint64_t id) const override;
    inline uint64_t select1(uint64_t id) const override;
    inline std::pair<bool, uint64_t> inverse_select(uint64_t id) const override;
    inline uint64_t conditional_rank1(uint64_t id) const override;

    inline uint64_t next1(uint64_t id) const override;
    inline uint64_t prev1(uint64_t id) const override;

    inline bool operator[](uint64_t id) const override;
    inline uint64_t get_int(uint64_t id, uint32_t width) const override;

    inline bool load(std::istream &in) override;
    inline void serialize(std::ostream &out) const override;

    inline uint64_t size() const override { return vector_.size(); }
    inline uint64_t num_set_bits() const override { return num_set_bits_; }

    inline void call_ones_in_range(uint64_t begin, uint64_t end,
                                   const VoidCall<uint64_t> &callback) const override;

    inline void add_to(sdsl::bit_vector *other) const override;

    inline sdsl::bit_vector to_vector() const override;

    inline const bv_type& data() const { return vector_; }

    /**
     * Predict space taken by the vector with given its parameters in bits.
     */
    static inline uint64_t predict_size(uint64_t size, uint64_t num_set_bits) {
       return bv_traits<bv_type>::predict_size(size, num_set_bits)
                + bv_traits<rank_1_type>::predict_size(size, num_set_bits)
                + bv_traits<select_1_type>::predict_size(size, num_set_bits)
                + bv_traits<select_0_type>::predict_size(size, num_set_bits);
    }

  private:
    bv_type vector_;
    rank_1_type rk1_;
    select_1_type slct1_;
    select_0_type slct0_;
    uint64_t num_set_bits_;
};


template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>&
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::operator=(const bit_vector_sdsl &other) {
    vector_ = other.vector_;
    rk1_ = other.rk1_;
    rk1_.set_vector(&vector_);
    slct1_ = other.slct1_;
    slct1_.set_vector(&vector_);
    slct0_ = other.slct0_;
    slct0_.set_vector(&vector_);
    num_set_bits_ = other.num_set_bits_;
    return *this;
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>&
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::operator=(bit_vector_sdsl&& other) noexcept {
    vector_ = std::move(other.vector_);
    rk1_ = std::move(other.rk1_);
    rk1_.set_vector(&vector_);
    slct1_ = std::move(other.slct1_);
    slct1_.set_vector(&vector_);
    slct0_ = std::move(other.slct0_);
    slct0_.set_vector(&vector_);
    num_set_bits_ = other.num_set_bits_;
    return *this;
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::rank1(uint64_t id) const {
    //the rank method in SDSL does not include id in the count
    return rk1_(id >= this->size() ? this->size() : id + 1);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::select0(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= size() - num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));
    return slct0_(id);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));
    return slct1_(id);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
std::pair<bool, uint64_t>
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::inverse_select(uint64_t id) const {
    if constexpr(is_rrr<bv_type>{}) {
        // TODO: implement inverse_select for sdsl::rrr_vector<15>
        if constexpr(bv_type::block_size != 15) {
            std::pair<bool, uint64_t> pair = vector_.inverse_select(id);
            pair.second += pair.first;
            return pair;
        }
    }
    // TODO: implement inverse_select for other bit vectors as well
    return bit_vector::inverse_select(id);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::conditional_rank1(uint64_t id) const {
    if constexpr(is_rrr<bv_type>{}) {
        // TODO: implement conditional_rank for sdsl::rrr_vector<15>
        if constexpr(bv_type::block_size != 15) {
            return vector_.conditional_rank(id);
        }
    }
    return bit_vector::conditional_rank1(id);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::next1(uint64_t id) const {
    assert(id < size());

    if constexpr(std::is_same_v<bv_type, sdsl::bit_vector>) {
        auto next = next_bit(vector_, id, bv_traits<bv_type>::MAX_ITER_BIT_VECTOR);
        if (next < vector_.size())
            return next;

        if (vector_.size() - id <= bv_traits<bv_type>::MAX_ITER_BIT_VECTOR)
            return size();

        uint64_t rk = rank1(id) + 1;
        return rk <= num_set_bits() ? select1(rk) : size();
    } else {
        return ::next1(*this, id, bv_traits<bv_type>::MAX_ITER_BIT_VECTOR);
    }
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::prev1(uint64_t id) const {
    assert(id < size());

    if constexpr(std::is_same_v<bv_type, sdsl::bit_vector>) {
        auto prev = prev_bit(vector_, id, bv_traits<bv_type>::MAX_ITER_BIT_VECTOR);
        if (prev <= id)
            return prev;

        if (id < bv_traits<bv_type>::MAX_ITER_BIT_VECTOR)
            return size();

        uint64_t rk = rank1(id);
        return rk ? select1(rk) : size();
    } else {
        return ::prev1(*this, id, bv_traits<bv_type>::MAX_ITER_BIT_VECTOR);
    }
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
bool
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
uint64_t
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::get_int(uint64_t id, uint32_t width) const {
    assert(id + width <= size());
    return vector_.get_int(id, width);
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
bool
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);
        if constexpr(std::is_same_v<bv_type, sdsl::bit_vector>) {
            load_number(in);
        }
        rk1_.load(in, &vector_);
        slct1_.load(in, &vector_);
        slct0_.load(in, &vector_);
        num_set_bits_ = rk1_(vector_.size());
        return in.good();

    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load "
                  << typeid(bv_type).name() << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
void
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::serialize(std::ostream &out) const {
    vector_.serialize(out);
    if constexpr(std::is_same_v<bv_type, sdsl::bit_vector>) {
        serialize_number(out, num_set_bits_);
    }
    rk1_.serialize(out);
    slct1_.serialize(out);
    slct0_.serialize(out);

    if (!out.good())
        throw std::ofstream::failure("Error when dumping bit_vector");
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
sdsl::bit_vector
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::to_vector() const {
    if constexpr(std::is_same_v<bv_type, sdsl::bit_vector>) {
        return vector_;
    } else {
        return bit_vector::to_vector_adaptive(
            bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR
        );
    }
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
void
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::call_ones_in_range(uint64_t begin, uint64_t end,
                     const VoidCall<uint64_t> &callback) const {
    assert(begin <= end);
    assert(end <= size());

    if constexpr(std::is_same_v<bv_type, sdsl::bit_vector>) {
        ::call_ones(vector_, begin, end, callback);
    } else {
        bit_vector::call_ones_adaptive(begin, end, callback,
            bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR
        );
    }
}

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
void
bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());

    if constexpr(std::is_same_v<bv_type, sdsl::bit_vector>) {
        *other |= vector_;
    } else {
        bit_vector::add_to_adaptive(other,
            bv_traits<bv_type>::SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR
        );
    }
}


template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
template <typename type>
struct bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::bv_traits {
    static inline uint64_t predict_size(uint64_t /*size*/, uint64_t /*num_set_bits*/) {
        return sizeof(type) * 8;
    }
};

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
template <uint32_t k_sblock_rate>
struct bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::bv_traits<sdsl::hyb_vector<k_sblock_rate>> {
    // hyb_vector doesn't support select
    static constexpr size_t MAX_ITER_BIT_VECTOR = 1000;
    static constexpr size_t SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR = 1000;

    static inline uint64_t predict_size(uint64_t /*size*/, uint64_t /*num_set_bits*/) {
        throw std::runtime_error(std::string("Error: unknown space taken for this bit_vector")
                                    + typeid(bv_type).name());
    }
};

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
template <uint32_t t_bs>
struct bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::bv_traits<sdsl::bit_vector_il<t_bs>> {
    static constexpr size_t MAX_ITER_BIT_VECTOR = 1000;
    static constexpr size_t SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR = 1000;

    static inline uint64_t predict_size(uint64_t size, uint64_t /*num_set_bits*/) {
        uint64_t num_blocks = (size + 64) / 64 + (size + t_bs) / t_bs + 1;
        return 4 * 64
            + num_blocks * 64 + 64
            + (num_blocks > 1024 * 64
                ? std::min(1024ULL, 1ULL << sdsl::bits::hi((size + t_bs) / t_bs)) * 64 + 64
                : 64ULL);
    }
};

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
template <uint16_t t_bs, class t_rac, uint16_t t_k>
struct bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::bv_traits<sdsl::rrr_vector<t_bs, t_rac, t_k>> {
    static constexpr size_t MAX_ITER_BIT_VECTOR = 1;
    static constexpr size_t SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR = 2;

    static inline double logbinomial(double n, double m) {
        return (lgamma(n + 1)
                    - lgamma(m + 1)
                    - lgamma(n - m + 1)) / log(2);
    }

    static inline uint64_t predict_size(uint64_t size, uint64_t num_set_bits) {
        // TODO: correct the formula for block_size = 15
        uint64_t bt = (size + t_bs) / t_bs;
        uint64_t blocks = bt * logbinomial(t_bs, t_bs * (double)num_set_bits / size);
        return (bt * (sdsl::bits::hi(t_bs) + 1) + 63) / 64 * 64 + 64 + 8
            + (blocks + 63) / 64 * 64 + 64
            + ((bt + t_k - 1) / t_k
                    * (sdsl::bits::hi(blocks) + 1) + 63) / 64 * 64 + 64 + 8
            + (((bt + t_k - 1) / t_k + ((size % (t_k * t_bs)) > 0))
                    * (sdsl::bits::hi(num_set_bits) + 1) + 63) / 64 * 64 + 64 + 8;
    }
};

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
template <uint8_t t_width>
struct bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::bv_traits<sdsl::int_vector<t_width>> {
    static constexpr size_t MAX_ITER_BIT_VECTOR = 1000;
    static constexpr size_t SEQ_BITWISE_WORD_ACCESS_VS_SELECT_FACTOR = 1000;

    static inline uint64_t predict_size(uint64_t size, uint64_t /*num_set_bits*/) {
        return sizeof(sdsl::int_vector<t_width>) * 8 + (size + 64) / 64 * 64;
    }
};

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
template <uint8_t t_b, uint8_t t_pat_len>
struct bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::bv_traits<sdsl::rank_support_v5<t_b, t_pat_len>> {
    static inline uint64_t predict_size(uint64_t size, uint64_t /*num_set_bits*/) {
        return sizeof(sdsl::rank_support_v5<t_b, t_pat_len>) * 8
            + footprint_rank_support_v5(size);
    }
};

template <class bv_type, class rank_1_type, class select_1_type, class select_0_type>
template <uint8_t t_b, uint8_t t_pat_len>
struct bit_vector_sdsl<bv_type, rank_1_type, select_1_type, select_0_type>
::bv_traits<sdsl::select_support_mcl<t_b, t_pat_len>> {
    static inline uint64_t predict_size(uint64_t size, uint64_t num_set_bits) {
        return sizeof(sdsl::select_support_mcl<t_b, t_pat_len>) * 8
            + footprint_select_support_mcl(size, num_set_bits);
    }
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

using bit_vector_stat
    = bit_vector_sdsl<sdsl::bit_vector,
                      sdsl::rank_support_v5<1>,
                      sdsl::select_support_mcl<1>,
                      sdsl::select_support_scan<0>>;

#endif // __BIT_VECTOR_SDSL_HPP__
