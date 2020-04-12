#ifndef __VECTOR_ALGORITHM_HPP__
#define __VECTOR_ALGORITHM_HPP__

#include <functional>
#include <cassert>

#include <sdsl/int_vector.hpp>
#include <sdsl/select_support_scan.hpp>

class ThreadPool;
class bit_vector;


sdsl::bit_vector to_sdsl(const std::vector<bool> &vector);
sdsl::bit_vector to_sdsl(const std::vector<uint8_t> &vector);

template <class Vector>
sdsl::int_vector<> pack_vector(const Vector &vector, uint8_t bits_per_number) {
    if constexpr(std::is_same_v<Vector, sdsl::int_vector<>>) {
        if (bits_per_number == vector.width())
            return vector;
    }

    sdsl::int_vector<> packed(vector.size(), 0, bits_per_number);
    for (uint64_t i = 0; i < vector.size(); ++i) {
        packed[i] = vector[i];
    }
    return packed;
}

sdsl::int_vector<> pack_vector(sdsl::int_vector<>&& vector,
                               uint8_t bits_per_number);


template <class Bitmap, class Callback>
void call_ones(const Bitmap &vector,
               uint64_t begin, uint64_t end,
               Callback callback) {
    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for (; i < end && i & 0x3F; ++i) {
        if (vector[i])
            callback(i);
    }
    uint64_t word;
    for (uint64_t j = i + 64; j <= end; j += 64) {
        word = vector.get_int(i, 64);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for (; i < j; ++i) {
            if (vector[i])
                callback(i);
        }
    }
    for (; i < end; ++i) {
        if (vector[i])
            callback(i);
    }
}

template <class Bitmap, class Callback>
void call_ones(const Bitmap &vector, Callback callback) {
    call_ones(vector, 0, vector.size(), callback);
}

template <class Bitmap, class Callback>
void call_zeros(const Bitmap &vector,
                uint64_t begin, uint64_t end,
                Callback callback) {
    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for (; i < end && i & 0x3F; ++i) {
        if (!vector[i])
            callback(i);
    }
    uint64_t word;
    for (uint64_t j = i + 64; j <= end; j += 64) {
        word = ~vector.get_int(i, 64);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for (; i < j; ++i) {
            if (!vector[i])
                callback(i);
        }
    }
    for (; i < end; ++i) {
        if (!vector[i])
            callback(i);
    }
}

template <class Bitmap, class Callback>
void call_zeros(const Bitmap &vector, Callback callback) {
    call_zeros(vector, 0, vector.size(), callback);
}

uint64_t count_ones(const sdsl::bit_vector &vector, uint64_t begin, uint64_t end);

uint64_t inner_prod(const sdsl::bit_vector &first, const sdsl::bit_vector &second);

void compute_or(const std::vector<const bit_vector *> &columns,
                sdsl::bit_vector *result,
                ThreadPool &thread_pool);

// assumes that all bits that are set in |column| are set in |reference| too.
sdsl::bit_vector generate_subindex(const bit_vector &column,
                                   const sdsl::bit_vector &reference,
                                   uint64_t reference_num_set_bits,
                                   ThreadPool &thread_pool);

// Apply the bitwise AND of vector with right-shifts of itself. Only works for
// values of offset < 64
sdsl::bit_vector autocorrelate(const sdsl::bit_vector &vector, uint8_t offset);


// Call (uint64_t index, uint64_t value) for each non-zero value in [begin, end)
template <class Callback>
void call_nonzeros(const sdsl::int_vector<> &vector,
                   uint64_t begin, uint64_t end,
                   Callback callback) {
    if (begin >= end)
        return;

    uint64_t i = begin;
    switch (vector.width()) {
        case 64:
            std::for_each(reinterpret_cast<const uint64_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint64_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        case 32:
            std::for_each(reinterpret_cast<const uint32_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint32_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        case 16:
            std::for_each(reinterpret_cast<const uint16_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint16_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        case 8:
            std::for_each(reinterpret_cast<const uint8_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint8_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        default:
            // vector.width() is not a power of two
            assert(vector.width() < 64);

            size_t begin_word = begin * vector.width() / 64;
            size_t end_word = (end * vector.width() + 63) / 64;

            for (uint64_t w = begin_word, it = begin; w < end_word; ++w) {
                if (!vector.data()[w])
                    continue;

                it = std::max(it, w * 64 / vector.width());

                auto it_end = std::min(end,
                    ((w + 1) * 64 + vector.width() - 1) / vector.width());

                for (; it < it_end; ++it) {
                    if (vector[it])
                        callback(it, vector[it]);
                }
            }
    }
}

// Call (uint64_t index, uint64_t value) for each non-zero value in |vector|.
template <class Callback>
void call_nonzeros(const sdsl::int_vector<> &vector, Callback callback) {
    return call_nonzeros(vector, 0, vector.size(), callback);
}

template <typename BitVector>
inline uint64_t next1(const BitVector &v,
                      uint64_t pos,
                      size_t num_steps) {
    assert(pos < v.size());

    for (size_t t = 1; t < num_steps; ++t, ++pos) {
        if (pos == v.size() || v[pos])
            return pos;
    }
    if (pos == v.size())
        return pos;

    uint64_t rk;

    if (num_steps >= 1) {
        auto pair = v.inverse_select(pos);
        if (pair.first)
            return pos;

        rk = pair.second + 1;

    } else {
        rk = pos ? v.rank1(pos - 1) + 1 : 1;
    }

    return rk <= v.num_set_bits() ? v.select1(rk) : v.size();
}

template <typename BitVector>
inline uint64_t prev1(const BitVector &v,
                      uint64_t pos,
                      size_t num_steps) {
    assert(pos < v.size());

    for (size_t t = 1; t < num_steps; ++t, --pos) {
        if (v[pos])
            return pos;

        if (pos == 0)
            return v.size();
    }

    uint64_t rk;

    if (num_steps >= 1) {
        auto pair = v.inverse_select(pos);
        if (pair.first)
            return pos;

        rk = pair.second;

    } else {
        rk = v.rank1(pos);
    }

    return rk ? v.select1(rk) : v.size();
}


// taken from https://github.com/xxsds/sdsl-lite/blob/master/include/sdsl/util.hpp
// this function has been modified.
template <class t_int_vec>
typename t_int_vec::size_type
next_bit(const t_int_vec &v,
         uint64_t idx,
         uint64_t max_steps = std::numeric_limits<uint64_t>::max()) {
    assert(idx < v.bit_size());

    uint64_t pos  = idx >> 6;
    uint64_t node = v.data()[pos];
    node >>= (idx & 0x3F);
    if (node)
        return std::min(idx + sdsl::bits::lo(node), v.bit_size());

    uint64_t end = idx + std::min(max_steps, v.bit_size() - idx);
    for (++pos; (pos << 6) < end; ++pos) {
        if (v.data()[pos])
            return std::min((pos << 6) | sdsl::bits::lo(v.data()[pos]), v.bit_size());
    }
    return v.bit_size();
}

// taken from https://github.com/xxsds/sdsl-lite/blob/master/include/sdsl/util.hpp
// this function has been modified.
template <class t_int_vec>
typename t_int_vec::size_type
prev_bit(const t_int_vec &v,
         uint64_t idx,
         uint64_t max_steps = std::numeric_limits<uint64_t>::max()) {
    assert(idx < v.bit_size());

    int64_t pos  = idx >> 6;
    uint64_t node = v.data()[pos];
    node <<= 63 - (idx & 0x3F);
    if (node)
        return idx - (63 - sdsl::bits::hi(node));

    // last position to visit: 0 or (idx + 1 - max_steps)
    int64_t r_end_word = ((idx + 1 - std::min(idx + 1, max_steps)) >> 6) - 1;
    assert(r_end_word >= -1);
    for (--pos; pos > r_end_word; --pos) {
        if (v.data()[pos])
            return (pos << 6) | sdsl::bits::hi(v.data()[pos]);
    }
    return v.bit_size();
}

/**
 * Atomic bit fetching, setting, and unsetting on packed vectors.
 * fetch_and_* return the old values. The default memorder __ATOMIC_SEQ_CST
 * enforces the ordering of writes and reads across threads. See
 * https://gcc.gnu.org/onlinedocs/gcc/_005f_005fatomic-Builtins.html
 * for other options.
 */

template <class t_int_vec>
inline bool atomic_fetch_and_set_bit(t_int_vec &v,
                                     size_t i,
                                     bool atomic = false,
                                     int memorder = __ATOMIC_SEQ_CST) {
    // these assume that the underlying vector contains packed 64-bit integers
    static_assert(sizeof(*v.data()) == 8);

    if (atomic) {
        return (__atomic_fetch_or(&v.data()[i >> 6],
                                  1llu << (i & 0x3F),
                                  memorder) >> (i & 0x3F)) & 1;
    } else {
        uint64_t *word = &v.data()[i >> 6];
        if (!((*word >> (i & 0x3F)) & 1)) {
            *word |= (1llu << (i & 0x3F));
            return false;
        } else {
            return true;
        }
    }
}

template <class t_int_vec>
inline bool atomic_fetch_and_unset_bit(t_int_vec &v,
                                       size_t i,
                                       bool atomic = false,
                                       int memorder = __ATOMIC_SEQ_CST) {
    // these assume that the underlying vector contains packed 64-bit integers
    static_assert(sizeof(*v.data()) == 8);

    if (atomic) {
        return (__atomic_fetch_and(&v.data()[i >> 6],
                                   ~(1llu << (i & 0x3F)),
                                   memorder) >> (i & 0x3F)) & 1;
    } else {
        uint64_t *word = &v.data()[i >> 6];
        if ((*word >> (i & 0x3F)) & 1) {
            *word &= ~(1llu << (i & 0x3F));
            return true;
        } else {
            return false;
        }
    }
}

template <class t_int_vec>
inline bool atomic_fetch_bit(t_int_vec &v,
                             size_t i,
                             bool atomic = false,
                             int memorder = __ATOMIC_SEQ_CST) {
    // these assume that the underlying vector contains packed 64-bit integers
    static_assert(sizeof(*v.data()) == 8);

    return atomic
        ? ((__atomic_load_n(&v.data()[i >> 6], memorder) >> (i & 0x3F)) & 1)
        : ((v.data()[i >> 6] >> (i & 0x3F)) & 1);
}

template <class t_int_vec>
inline void atomic_set_bit(t_int_vec &v,
                           size_t i,
                           bool atomic = false,
                           int memorder = __ATOMIC_SEQ_CST) {
    // these assume that the underlying vector contains packed 64-bit integers
    static_assert(sizeof(*v.data()) == 8);

    if (atomic) {
        __atomic_or_fetch(&v.data()[i >> 6], 1llu << (i & 0x3F), memorder);
    } else {
        v.data()[i >> 6] |= (1llu << (i & 0x3F));
    }
}

template <class t_int_vec>
inline void atomic_unset_bit(t_int_vec &v,
                             size_t i,
                             bool atomic = false,
                             int memorder = __ATOMIC_SEQ_CST) {
    // these assume that the underlying vector contains packed 64-bit integers
    static_assert(sizeof(*v.data()) == 8);

    if (atomic) {
        __atomic_and_fetch(&v.data()[i >> 6], ~(1llu << (i & 0x3F)), memorder);
    } else {
        v.data()[i >> 6] &= ~(1llu << (i & 0x3F));
    }
}

namespace sdsl {

// based on sdsl::select_support_scan
template <uint8_t t_b = 1, uint8_t t_pat_len = 1>
class select_support_scan_offset : public select_support_scan<t_b, t_pat_len> {
  public:
    using typename select_support_scan<t_b, t_pat_len>::size_type;

    explicit select_support_scan_offset(const bit_vector *v = nullptr)
          : select_support_scan<t_b, t_pat_len>::select_support_scan(v) {}

    select_support_scan_offset(const select_support_scan<t_b,t_pat_len> &ss)
          : select_support_scan<t_b, t_pat_len>::select_support_scan(ss.m_v) {}

    inline size_type select_offset(size_type i, size_type offset = 0) const {
        using trait = select_support_trait<t_b, t_pat_len>;
        const uint64_t *data = this->m_v->data() + (offset >> 6);
        size_type word_pos = offset >> 6;
        size_type word_off = offset & 0x3F;
        uint64_t carry = trait::init_carry(data, word_pos);
        size_type args = trait::args_in_the_first_word(*data, word_off, carry);
        if (args >= i) {
            return (word_pos << 6)
                + trait::ith_arg_pos_in_the_first_word(*data, i, word_off, carry);
        }
        word_pos++;
        size_type sum_args = args;
        carry = trait::get_carry(*data);
        uint64_t old_carry = carry;
        args = trait::args_in_the_word(*(++data), carry);
        while (sum_args + args < i) {
            sum_args += args;
            assert(data + 1 < this->m_v->data() + (this->m_v->capacity() >> 6));
            old_carry = carry;
            args = trait::args_in_the_word(*(++data), carry);
            word_pos++;
        }
        return (word_pos << 6)
            + trait::ith_arg_pos_in_the_word(*data, i - sum_args, old_carry);
    }
};

} // namespace sdsl

// Predict the memory footprint in bits for sdsl::sd_vector<>
uint64_t footprint_sd_vector(uint64_t size, uint64_t num_set_bits);

// Predict the memory footprint in bits for sdsl::select_support_mcl<>
uint64_t footprint_select_support_mcl(uint64_t size, uint64_t num_set_bits);

// Predict the memory footprint in bits for sdsl::rank_support_v5<>
uint64_t footprint_rank_support_v5(uint64_t size);

#endif // __VECTOR_ALGORITHM_HPP__
