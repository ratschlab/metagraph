#ifndef __VECTOR_ALGORITHM_HPP__
#define __VECTOR_ALGORITHM_HPP__

#include <functional>
#include <cassert>

#include <sdsl/int_vector.hpp>
#include <sdsl/select_support_scan.hpp>

class ThreadPool;
class bit_vector;
class bit_vector_stat;


sdsl::bit_vector to_sdsl(const std::vector<bool> &vector);
sdsl::bit_vector to_sdsl(const std::vector<uint8_t> &vector);

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
        word = vector.get_int(i);
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
        word = ~vector.get_int(i);
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
                                   const bit_vector_stat &reference);

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

#endif // __VECTOR_ALGORITHM_HPP__
