#ifndef __VECTOR_ALGORITHM_HPP__
#define __VECTOR_ALGORITHM_HPP__

#include <functional>
#include <cassert>
#include <cstdlib>

#include <sdsl/int_vector.hpp>
#include <sdsl/select_support_scan.hpp>

class ThreadPool;
class bit_vector;


sdsl::bit_vector to_sdsl(const std::vector<bool> &vector);
sdsl::bit_vector to_sdsl(const std::vector<uint8_t> &vector);

template <class Vector>
inline sdsl::int_vector<> pack_vector(const Vector &vector, uint8_t bits_per_number) {
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

sdsl::int_vector<> pack_vector(sdsl::int_vector<>&& vector, uint8_t bits_per_number);


/**
 * Atomic bit fetching, setting, and unsetting on packed vectors.
 * fetch_and_* return the old values. The default memory order __ATOMIC_SEQ_CST
 * enforces the ordering of writes and reads across threads. See
 * https://gcc.gnu.org/onlinedocs/gcc/_005f_005fatomic-Builtins.html
 * for more details.
 */

inline bool fetch_and_set_bit(uint64_t *v,
                              uint64_t i,
                              bool atomic = false,
                              int mo = __ATOMIC_SEQ_CST);

inline bool fetch_and_unset_bit(uint64_t *v,
                                uint64_t i,
                                bool atomic = false,
                                int mo = __ATOMIC_SEQ_CST);

inline bool fetch_bit(const uint64_t *v,
                      uint64_t i,
                      bool atomic = false,
                      int mo = __ATOMIC_SEQ_CST);

inline void set_bit(uint64_t *v,
                    uint64_t i,
                    bool atomic = false,
                    int mo = __ATOMIC_SEQ_CST);

inline void unset_bit(uint64_t *v,
                      uint64_t i,
                      bool atomic = false,
                      int mo = __ATOMIC_SEQ_CST);

inline uint64_t atomic_fetch(const sdsl::int_vector<> &vector, uint64_t i,
                             std::mutex &backup_mutex,
                             int mo = __ATOMIC_SEQ_CST);

inline uint64_t atomic_fetch_and_add(sdsl::int_vector<> &vector, uint64_t i,
                                     uint64_t count,
                                     std::mutex &backup_mutex,
                                     int mo = __ATOMIC_SEQ_CST);

inline uint64_t atomic_exchange(sdsl::int_vector<> &vector, uint64_t i, uint64_t val,
                                std::mutex &backup_mutex,
                                int mo = __ATOMIC_SEQ_CST);


template <class Bitmap, class Callback>
inline void call_ones(const Bitmap &vector,
                      uint64_t begin, uint64_t end,
                      Callback callback);

template <class Callback>
inline void call_ones(const sdsl::bit_vector &vector,
                      uint64_t begin, uint64_t end,
                      Callback callback,
                      bool atomic,
                      int mo = __ATOMIC_SEQ_CST);

template <class Bitmap, class Callback>
inline void call_ones(const Bitmap &vector, Callback callback) {
    call_ones(vector, 0, vector.size(), callback);
}

template <class Callback>
inline void call_ones(const sdsl::bit_vector &vector,
                      Callback callback,
                      bool atomic,
                      int mo = __ATOMIC_SEQ_CST) {
    call_ones(vector, 0, vector.size(), callback, atomic, mo);
}

template <class Bitmap, class Callback>
inline void call_zeros(const Bitmap &vector,
                       uint64_t begin, uint64_t end,
                       Callback callback);

template <class Callback>
inline void call_zeros(const sdsl::bit_vector &vector,
                       uint64_t begin, uint64_t end,
                       Callback callback,
                       bool atomic,
                       int mo = __ATOMIC_SEQ_CST);

template <class Bitmap, class Callback>
inline void call_zeros(const Bitmap &vector, Callback callback) {
    call_zeros(vector, 0, vector.size(), callback);
}

template <class Callback>
inline void call_zeros(const sdsl::bit_vector &vector,
                       Callback callback,
                       bool atomic,
                       int mo = __ATOMIC_SEQ_CST) {
    call_zeros(vector, 0, vector.size(), callback, atomic, mo);
}

uint64_t count_ones(const sdsl::bit_vector &vector, uint64_t begin, uint64_t end);

uint64_t inner_prod(const sdsl::bit_vector &first, const sdsl::bit_vector &second);

void compute_or(const std::vector<const bit_vector *> &columns,
                sdsl::bit_vector *result,
                ThreadPool &thread_pool);

// Call this version only for sparse vectors (with the density about 1% or less).
// The buffer must have capacity to store 3 x (number of set bits in all columns)
// 64-bit integers.
std::unique_ptr<bit_vector> compute_or(const std::vector<const bit_vector *> &columns,
                                       uint64_t *buffer,
                                       ThreadPool &thread_pool);

// Assumes that all bits that are set in |column| are set in |reference| too
sdsl::bit_vector generate_subindex(const bit_vector &column,
                                   const sdsl::bit_vector &reference,
                                   uint64_t reference_num_set_bits,
                                   ThreadPool &thread_pool);
// Assumes that all bits that are set in |column| are set in |reference| too
sdsl::bit_vector generate_subindex(const bit_vector &column,
                                   const bit_vector &reference,
                                   ThreadPool &thread_pool);

// Apply the bitwise AND of vector with right-shifts of itself. Only works for
// values of offset < 64
sdsl::bit_vector autocorrelate(const sdsl::bit_vector &vector, uint8_t offset);


// Call (uint64_t index, uint64_t value) for each non-zero value in [begin, end)
template <class Callback>
inline void call_nonzeros(const sdsl::int_vector<> &vector,
                          uint64_t begin, uint64_t end,
                          Callback callback);

// Call (uint64_t index, uint64_t value) for each non-zero value in |vector|.
template <class Callback>
void call_nonzeros(const sdsl::int_vector<> &vector, Callback callback) {
    return call_nonzeros(vector, 0, vector.size(), callback);
}

template <typename BitVector>
inline uint64_t next1(const BitVector &v, uint64_t pos, size_t num_steps);

template <typename BitVector>
inline uint64_t prev1(const BitVector &v, uint64_t pos, size_t num_steps);


// taken from https://github.com/xxsds/sdsl-lite/blob/master/include/sdsl/util.hpp
// this function has been modified.
template <class t_int_vec>
typename t_int_vec::size_type
next_bit(const t_int_vec &v,
         uint64_t idx,
         uint64_t max_steps = std::numeric_limits<uint64_t>::max());

// taken from https://github.com/xxsds/sdsl-lite/blob/master/include/sdsl/util.hpp
// this function has been modified.
template <class t_int_vec>
typename t_int_vec::size_type
prev_bit(const t_int_vec &v,
         uint64_t idx,
         uint64_t max_steps = std::numeric_limits<uint64_t>::max());


// Return an int_vector whose limbs are aligned to alignment bytes.
// This is useful for algorithms requiring access to larger words in the underlying
// vector (atomic increment, SIMD, etc.)
template <int t_width = 0>
inline sdsl::int_vector<t_width> aligned_int_vector(size_t size = 0, uint64_t val = 0,
                                                    uint8_t width = 0,
                                                    size_t alignment = 8);

inline sdsl::bit_vector aligned_bit_vector(size_t size = 0,
                                           bool val = 0,
                                           size_t alignment = 8) {
    return aligned_int_vector<1>(size, val, 1, alignment);
}


// Predict the memory footprint in bits for sdsl::sd_vector<>
uint64_t footprint_sd_vector(uint64_t size, uint64_t num_set_bits);

// Predict the memory footprint in bits for sdsl::select_support_mcl<>
uint64_t footprint_select_support_mcl(uint64_t size, uint64_t num_set_bits);

// Predict the memory footprint in bits for sdsl::rank_support_v5<>
uint64_t footprint_rank_support_v5(uint64_t size);


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


////////////////////////////////////////////////////////////////
//                      IMPLEMENTATIONS                       //
////////////////////////////////////////////////////////////////

inline bool fetch_and_set_bit(uint64_t *v, uint64_t i, bool atomic, int mo) {
    const uint64_t mask = (1llu << (i & 0x3F));

    if (atomic) {
        return __atomic_fetch_or(&v[i >> 6], mask, mo) & mask;
    } else {
        uint64_t &word = v[i >> 6];
        if (word & mask) {
            return true;
        } else {
            word |= mask;
            return false;
        }
    }
}

inline bool fetch_and_unset_bit(uint64_t *v, uint64_t i, bool atomic, int mo) {
    const uint64_t mask = (1llu << (i & 0x3F));

    if (atomic) {
        return __atomic_fetch_and(&v[i >> 6], ~mask, mo) & mask;
    } else {
        uint64_t &word = v[i >> 6];
        if (word & mask) {
            word &= ~mask;
            return true;
        } else {
            return false;
        }
    }
}

inline bool fetch_bit(const uint64_t *v, uint64_t i, bool atomic, int mo) {
    return atomic
        ? ((__atomic_load_n(&v[i >> 6], mo) >> (i & 0x3F)) & 1)
        : ((v[i >> 6] >> (i & 0x3F)) & 1);
}

inline void set_bit(uint64_t *v, uint64_t i, bool atomic, int mo) {
    if (atomic) {
        __atomic_or_fetch(&v[i >> 6], 1llu << (i & 0x3F), mo);
    } else {
        v[i >> 6] |= (1llu << (i & 0x3F));
    }
}

inline void unset_bit(uint64_t *v, uint64_t i, bool atomic, int mo) {
    if (atomic) {
        __atomic_and_fetch(&v[i >> 6], ~(1llu << (i & 0x3F)), mo);
    } else {
        v[i >> 6] &= ~(1llu << (i & 0x3F));
    }
}

inline uint64_t atomic_fetch(const sdsl::int_vector<> &vector,
                             uint64_t i,
                             std::mutex &backup_mutex,
                             int mo) {
    size_t width = vector.width();
    size_t bit_pos = i * width;
    uint64_t mask = (1llu << width) - 1;
    if (width + 7 > 64) {
        // there is no way to reliably modify without a mutex
        std::lock_guard<std::mutex> lock(backup_mutex);
        return vector[i];
    } else if ((bit_pos & 0x3F) + width <= 64) {
        // element fits in an aligned word
        const uint64_t *word = &vector.data()[bit_pos >> 6];
        return (__atomic_load_n(word, mo) >> (bit_pos & 0x3F)) & mask;
#if defined(MODE_TI) && defined(__CX16__)
    } else if ((bit_pos & 0x7F) + width <= 128) {
        // element fits in an aligned double word
        const __uint128_t *word = &reinterpret_cast<const __uint128_t*>(vector.data())[bit_pos >> 7];
        return (__atomic_load_n(word, mo) >> (bit_pos & 0x7F)) & mask;
#endif
    } else {
        uint8_t shift = bit_pos & 0x7;
        const uint8_t *word = &reinterpret_cast<const uint8_t*>(vector.data())[bit_pos >> 3];
        if (shift + width <= 8) {
            // read from a byte
            return (__atomic_load_n(word, mo) >> shift) & mask;
        } else if (shift + width <= 16) {
            // unaligned read from two bytes
            return (__atomic_load_n((const uint16_t*)word, mo) >> shift) & mask;
        } else if (shift + width <= 32) {
            // unaligned read from four bytes
            return (__atomic_load_n((const uint32_t*)word, mo) >> shift) & mask;
        } else if (shift + width <= 64) {
            // unaligned read from eight bytes
            return (__atomic_load_n((const uint64_t*)word, mo) >> shift) & mask;
        } else {
            assert(false);
        }
    }

    return 0;
}

inline uint64_t atomic_fetch_and_add(sdsl::int_vector<> &vector,
                                     uint64_t i,
                                     uint64_t val,
                                     std::mutex &backup_mutex,
                                     int mo) {
    size_t width = vector.width();
    size_t bit_pos = i * width;
    uint64_t mask = (1llu << width) - 1;
    if (width + 7 > 64) {
        // there is no way to reliably modify without a mutex
        std::lock_guard<std::mutex> lock(backup_mutex);
        uint64_t old_val = vector[i];
        vector[i] += val;
        return old_val;
    } else if ((bit_pos & 0x3F) + width <= 64) {
        // element fits in an aligned word
        uint64_t *word = &vector.data()[bit_pos >> 6];
        uint8_t shift = bit_pos & 0x3F;
        return (__atomic_fetch_add(word, val << shift, mo) >> shift) & mask;
#if defined(MODE_TI) && defined(__CX16__)
    } else if ((bit_pos & 0x7F) + width <= 128) {
        // element fits in an aligned double word
        __uint128_t *word = &reinterpret_cast<__uint128_t*>(vector.data())[bit_pos >> 7];
        uint8_t shift = bit_pos & 0x7F;
        // TODO: GCC only generates cmpxchg16b instruction only with older __sync functions
        return (__sync_fetch_and_add(word, __uint128_t(val) << shift) >> shift) & mask;
#endif
    } else {
        uint8_t shift = bit_pos & 0x7;
        uint8_t *word = &reinterpret_cast<uint8_t*>(vector.data())[bit_pos >> 3];
        if (shift + width <= 8) {
            // read from a byte
            return (__atomic_fetch_add(word, val << shift, mo) >> shift) & mask;
        } else if (shift + width <= 16) {
            // unaligned read from two bytes
            return (__atomic_fetch_add((uint16_t*)word, val << shift, mo) >> shift) & mask;
        } else if (shift + width <= 32) {
            // unaligned read from four bytes
            return (__atomic_fetch_add((uint32_t*)word, val << shift, mo) >> shift) & mask;
        } else if (shift + width <= 64) {
            // unaligned read from eight bytes
            return (__atomic_fetch_add((uint64_t*)word, val << shift, mo) >> shift) & mask;
        } else {
            assert(false);
        }
    }

    return 0;
}

inline uint64_t atomic_exchange(sdsl::int_vector<> &vector,
                                uint64_t i,
                                uint64_t val,
                                std::mutex &backup_mutex,
                                int mo) {
    size_t width = vector.width();
    size_t bit_pos = i * width;
    uint64_t mask = (1llu << width) - 1;
    if (width + 7 > 64) {
        // there is no way to reliably modify without a mutex
        std::lock_guard<std::mutex> lock(backup_mutex);
        uint64_t old_val = vector[i];
        vector[i] = val;
        return old_val;
    } else if ((bit_pos & 0x3F) + width <= 64) {
        // element fits in an aligned word
        uint64_t *word = &vector.data()[bit_pos >> 6];
        uint8_t shift = bit_pos & 0x3F;
        uint64_t desired;
        uint64_t exp = *word;
        uint64_t inv_mask = ~(mask << shift);
        val <<= shift;
        do {
            desired = val | (inv_mask & exp);
        } while (!__atomic_compare_exchange(word, &exp, &desired, true, mo, __ATOMIC_RELAXED));
        return (exp >> shift) & mask;
#if defined(MODE_TI) && defined(__CX16__)
    } else if ((bit_pos & 0x7F) + width <= 128) {
        // element fits in an aligned double word
        __uint128_t *word = &reinterpret_cast<__uint128_t*>(vector.data())[bit_pos >> 7];
        uint8_t shift = bit_pos & 0x7F;
        __uint128_t desired;
        __uint128_t exp;
        __uint128_t inv_mask = ~(__uint128_t(mask) << shift);
        __uint128_t big_val = __uint128_t(val) << shift;
        // TODO: GCC only generates cmpxchg16b instruction only with older __sync functions
        do {
            exp = *word;
            desired = big_val | (inv_mask & exp);
        } while (!__sync_bool_compare_and_swap(word, exp, desired));
        return (exp >> shift) & mask;
#endif
    } else {
        uint8_t shift = bit_pos & 0x7;
        uint8_t *word = &reinterpret_cast<uint8_t*>(vector.data())[bit_pos >> 3];
        if (shift + width <= 8) {
            // read from a byte
            uint8_t desired;
            uint8_t exp = *word;
            uint8_t inv_mask = ~(mask << shift);
            val <<= shift;
            do {
                desired = val | (inv_mask & exp);
            } while (!__atomic_compare_exchange(word, &exp, &desired, true, mo, __ATOMIC_RELAXED));
            return (exp >> shift) & mask;
        } else if (shift + width <= 16) {
            // unaligned read from two bytes
            uint16_t *this_word = (uint16_t*)word;
            uint16_t desired;
            uint16_t exp = *this_word;
            uint16_t inv_mask = ~(mask << shift);
            val <<= shift;
            do {
                desired = val | (inv_mask & exp);
            } while (!__atomic_compare_exchange(this_word, &exp, &desired, true, mo, __ATOMIC_RELAXED));
            return (exp >> shift) & mask;
        } else if (shift + width <= 32) {
            // unaligned read from four bytes
            uint32_t *this_word = (uint32_t*)word;
            uint32_t desired;
            uint32_t exp = *this_word;
            uint32_t inv_mask = ~(mask << shift);
            val <<= shift;
            do {
                desired = val | (inv_mask & exp);
            } while (!__atomic_compare_exchange(this_word, &exp, &desired, true, mo, __ATOMIC_RELAXED));
            return (exp >> shift) & mask;
        } else if (shift + width <= 64) {
            // unaligned read from eight bytes
            uint64_t *this_word = (uint64_t*)word;
            uint64_t desired;
            uint64_t exp = *this_word;
            uint64_t inv_mask = ~(mask << shift);
            val <<= shift;
            do {
                desired = val | (inv_mask & exp);
            } while (!__atomic_compare_exchange(this_word, &exp, &desired, true, mo, __ATOMIC_RELAXED));
            return (exp >> shift) & mask;
        } else {
            assert(false);
        }
    }

    return 0;
}


template <class Bitmap, class Callback>
void call_ones(const Bitmap &vector, uint64_t begin, uint64_t end, Callback callback) {
    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for ( ; i < end && (i & 0x3F); ++i) {
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

        for ( ; i < j; ++i) {
            if (vector[i])
                callback(i);
        }
    }
    for ( ; i < end; ++i) {
        if (vector[i])
            callback(i);
    }
}

template <class Callback>
void call_ones(const sdsl::bit_vector &vector,
               uint64_t begin, uint64_t end,
               Callback callback,
               bool atomic,
               int mo) {
    if (!atomic) {
        call_ones(vector, begin, end, callback);
        return;
    }

    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for ( ; i < end && (i & 0x3F); ++i) {
        if (fetch_bit(vector.data(), i, true))
            callback(i);
    }
    uint64_t word;
    for (uint64_t j = i + 64; j <= end; j += 64) {
        word = __atomic_load_n(&vector.data()[i >> 6], mo);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for ( ; i < j; ++i) {
            if (fetch_bit(vector.data(), i, true))
                callback(i);
        }
    }
    for ( ; i < end; ++i) {
        if (fetch_bit(vector.data(), i, true))
            callback(i);
    }
}

template <class Bitmap, class Callback>
void call_zeros(const Bitmap &vector, uint64_t begin, uint64_t end, Callback callback) {
    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for ( ; i < end && (i & 0x3F); ++i) {
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

        for ( ; i < j; ++i) {
            if (!vector[i])
                callback(i);
        }
    }
    for ( ; i < end; ++i) {
        if (!vector[i])
            callback(i);
    }
}

template <class Callback>
void call_zeros(const sdsl::bit_vector &vector,
                uint64_t begin, uint64_t end,
                Callback callback,
                bool atomic,
                int mo) {
    if (!atomic) {
        call_zeros(vector, begin, end, callback);
        return;
    }

    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for ( ; i < end && (i & 0x3F); ++i) {
        if (!fetch_bit(vector.data(), i, true))
            callback(i);
    }
    uint64_t word;
    for (uint64_t j = i + 64; j <= end; j += 64) {
        word = ~__atomic_load_n(&vector.data()[i >> 6], mo);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for ( ; i < j; ++i) {
            if (!fetch_bit(vector.data(), i, true))
                callback(i);
        }
    }
    for ( ; i < end; ++i) {
        if (!fetch_bit(vector.data(), i, true))
            callback(i);
    }
}

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

template <typename BitVector>
inline uint64_t next1(const BitVector &v, uint64_t pos, size_t num_steps) {
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
inline uint64_t prev1(const BitVector &v, uint64_t pos, size_t num_steps) {
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
inline typename t_int_vec::size_type
next_bit(const t_int_vec &v, uint64_t idx, uint64_t max_steps) {
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
inline typename t_int_vec::size_type
prev_bit(const t_int_vec &v, uint64_t idx, uint64_t max_steps) {
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

template <int t_width>
inline sdsl::int_vector<t_width>
aligned_int_vector(size_t size, uint64_t val, uint8_t width, size_t alignment) {
    assert(__builtin_popcountll(alignment) == 1);
    assert(alignment >= 8);

    // This is a dirty hack to allow for reallocating an int_vector<t_width>'s
    // underlying storage
    struct int_vector_access {
        typename sdsl::int_vector<t_width>::size_type m_size;
        uint64_t *m_data;
        typename sdsl::int_vector<t_width>::int_width_type m_width;
    };
    static_assert(sizeof(sdsl::int_vector<t_width>) == sizeof(int_vector_access));
    assert(!t_width || t_width == width);

    sdsl::int_vector<t_width> v;
    auto &v_cast = reinterpret_cast<int_vector_access&>(v);
    v_cast.m_size = size * width;
    v_cast.m_width = width;
    free(v_cast.m_data);

    // Round up to the nearest multiple of |alignment| bytes
    size_t capacity_bytes = ((((v_cast.m_size + 7) >> 3) + alignment - 1) / alignment) * alignment;
    if (posix_memalign((void**)&v_cast.m_data, alignment, capacity_bytes) || !v_cast.m_data)
        throw std::bad_alloc();

    memset(v_cast.m_data, 0, capacity_bytes);

    if (val)
        sdsl::util::set_to_value(v, val);

    return v;
}

#endif // __VECTOR_ALGORITHM_HPP__
