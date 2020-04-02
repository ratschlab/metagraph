#ifndef __BIT_VECTOR_HPP__
#define __BIT_VECTOR_HPP__

#include <cstdint>

#include <sdsl/int_vector.hpp>

#include "bitmap.hpp"


// Bitmap with rank and select operations ( rank/select dictionary )
class bit_vector : public bitmap {
  public:
    virtual ~bit_vector() {}

    // Computes the number of set bits in the subarray indexed by [0,1,...,id]
    virtual uint64_t rank1(uint64_t id) const = 0;
    virtual uint64_t rank0(uint64_t id) const {
        return std::min(id + 1, size()) - rank1(id);
    }
    // Returns the i-th set bit, starting from 1
    virtual uint64_t select1(uint64_t i) const = 0;
    // Returns the i-th unset bit, starting from 1
    virtual uint64_t select0(uint64_t i) const = 0;
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

  protected:
    // routines using density-based adaptive iteration with select and access
    void call_ones_adaptive(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback,
                            double WORD_ACCESS_VS_SELECT_FACTOR) const;

    void add_to_adaptive(sdsl::bit_vector *other,
                         double WORD_ACCESS_VS_SELECT_FACTOR) const;

    sdsl::bit_vector to_vector_adaptive(double WORD_ACCESS_VS_SELECT_FACTOR) const;
};

std::ostream& operator<<(std::ostream &os, const bit_vector &bv);

#endif // __BIT_VECTOR_HPP__
