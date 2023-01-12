#ifndef __SD_VECTOR_DISK__
#define __SD_VECTOR_DISK__

#include <type_traits>

#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/util.hpp>
#include <sdsl/iterators.hpp>

// Based on https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/sd_vector.hpp

namespace sdsl {

// forward declaration needed for friend declaration
template<uint8_t t_b          = 1,
         class t_hi_bit_vector= bit_vector,
         class t_select_1     = typename t_hi_bit_vector::select_1_type,
         class t_select_0     = typename t_hi_bit_vector::select_0_type>
class rank_support_sd_disk;  // in sd_vector

// forward declaration needed for friend declaration
template<uint8_t t_b          = 1,
         class t_hi_bit_vector= bit_vector,
         class t_select_1     = typename t_hi_bit_vector::select_1_type,
         class t_select_0     = typename t_hi_bit_vector::select_0_type>
class select_support_sd_disk;  // in sd_vector

//! Class for in-place construction of sd_vector from a strictly increasing sequence
/*! \par Building an sd_vector will clear the builder.
 */
class sd_vector_disk_builder
{
        template<typename, typename, typename>
        friend class sd_vector_disk;

    public:
        typedef bit_vector::size_type size_type;

    private:
        size_type m_size, m_capacity;
        uint8_t   m_wl;
        size_type m_tail, m_items;
        size_type m_last_high, m_highpos;

        std::string filename;
        size_t file_offset;
        int_vector_buffer<> m_low;
        bit_vector   m_high;

    public:
        sd_vector_disk_builder() :
            m_size(0), m_capacity(0),
            m_wl(0),
            m_tail(0), m_items(0),
            m_last_high(0), m_highpos(0) {}

        //! Constructor
        /*! \param n Vector size.
         *  \param m The number of 1-bits.
         *  \param filename The buffer filename to where the bitmap will eventually be stored.
         */
        sd_vector_disk_builder(size_type n, size_type m, const std::string &fname, size_t offset = 0) :
            m_size(n), m_capacity(m),
            m_wl(0),
            m_tail(0), m_items(0),
            m_last_high(0), m_highpos(0)
        {
            if(m_capacity > m_size)
                throw std::runtime_error("sd_vector_builder: requested capacity is larger than vector size.");

            size_type logm = bits::hi(m_capacity) + 1;
            const size_type logn = bits::hi(m_size) + 1;
            if(logm == logn)
            {
                --logm; // to ensure logn-logm > 0
                assert(logn - logm > 0);
            }
            m_wl = logn - logm;

            std::ofstream out(filename, (file_offset ? std::ios::in : std::ios::openmode(0))|std::ios::binary);
            out.seekp(file_offset);
            write_member(m_size, out);
            write_member(m_wl, out);
            out.close();

            filename = fname;
            file_offset = offset;
            m_low = int_vector_buffer<>(filename, std::ios::out, 1024 * 1024, m_wl, false,
                                        file_offset + sizeof(m_size) + sizeof(m_wl));
            m_high = bit_vector(m_capacity + (1ULL << logm), 0);
        }

        inline size_type size() const { return m_size; }
        inline size_type capacity() const { return m_capacity; }
        inline size_type tail() const { return m_tail; }
        inline size_type items() const { return m_items; }

        //! Set a bit to 1.
        /*! \param i The position of the bit.
         *  \par The position must be strictly greater than for the previous call.
         */
        inline void set(size_type i)
        {
            assert(i >= m_tail && i < m_size);
            assert(m_items < m_capacity);

            size_type cur_high = i >> m_wl;
            m_highpos += (cur_high - m_last_high);
            m_last_high = cur_high;
            m_low[m_items++] = i; // int_vector truncates the most significant logm bits
            m_high[m_highpos++] = 1;  // write 1 for the entry
            m_tail = i + 1;
        }
};


//! A bit vector which compresses very sparse populated bit vectors by
// representing the positions of 1 by the Elias-Fano representation for non-decreasing sequences
/*!
 * \par Other implementations of this data structure:
 *  - the sdarray of Okanohara and Sadakane
 *  - Sebastiano Vigna implemented a elias_fano class in this sux library.
 *
 * \par References
 *  - P. Elias: ,,Efficient storage and retrieval by content and address of static files'',
 *              Journal of the ACM, 1974
 *  - R. Fano: ,,On the number of bits required to implement an associative memory''.
 *             Memorandum 61. Computer Structures Group, Project MAC, MIT, 1971
 *  - D. Okanohara, K. Sadakane: ,,Practical Entropy-Compressed Rank/Select Dictionary'',
 *             Proceedings of ALENEX 2007.
 *
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<class t_hi_bit_vector = bit_vector,
         class t_select_1     = typename t_hi_bit_vector::select_1_type,
         class t_select_0     = typename t_hi_bit_vector::select_0_type>
class sd_vector_disk
{
    public:
        typedef bit_vector::size_type                   size_type;
        typedef size_type                               value_type;
        typedef bit_vector::difference_type             difference_type;
        typedef random_access_const_iterator<sd_vector_disk> iterator;
        typedef iterator                                const_iterator;
        typedef bv_tag                                  index_category;
        typedef t_select_0                              select_0_support_type;
        typedef t_select_1                              select_1_support_type;

        typedef rank_support_sd_disk<0, t_hi_bit_vector, select_1_support_type, select_0_support_type> rank_0_type;
        typedef rank_support_sd_disk<1, t_hi_bit_vector, select_1_support_type, select_0_support_type> rank_1_type;
        typedef select_support_sd_disk<0, t_hi_bit_vector, select_1_support_type, select_0_support_type> select_0_type;
        typedef select_support_sd_disk<1, t_hi_bit_vector, select_1_support_type, select_0_support_type> select_1_type;

        typedef t_hi_bit_vector hi_bit_vector_type;
    private:
        // we need this variables to represent the m ones of the original bit vector of size n
        size_type m_size = 0;  // length of the original bit vector
        uint8_t   m_wl   = 0;  // log n - log m, where n is the length of the original bit vector
        // and m is the number of ones in the bit vector, wl is the abbreviation
        // for ,,width (of) low (part)''

        int_vector_buffer<>   m_low;           // vector for the least significant bits of the positions of the m ones
        hi_bit_vector_type    m_high;          // bit vector that represents the most significant bit in permuted order
        select_1_support_type m_high_1_select; // select support for the ones in m_high
        select_0_support_type m_high_0_select; // select support for the zeros in m_high

    public:
        const uint8_t&               wl            = m_wl;
        const hi_bit_vector_type&    high          = m_high;
        const int_vector_buffer<>&   low           = m_low;
        const select_1_support_type& high_1_select = m_high_1_select;
        const select_0_support_type& high_0_select = m_high_0_select;

        sd_vector_disk() { }

        sd_vector_disk(sd_vector_disk_builder& builder)
        {
            if (builder.items() < builder.capacity()) {
                throw std::runtime_error("sd_vector_disk: the builder is not full.");
            } else if (builder.items() > builder.capacity()) {
                throw std::runtime_error("sd_vector_disk: builder overflow.");
            }

            m_size = builder.m_size;
            m_wl = builder.m_wl;
            // flush m_low to file
            builder.m_low.close();

            if constexpr(std::is_same<hi_bit_vector_type, bit_vector>::value) {
                m_high.swap(builder.m_high);
            } else {
                util::assign(m_high, builder.m_high);
            }
            util::init_support(m_high_1_select, &m_high);
            util::init_support(m_high_0_select, &m_high);

            std::ofstream out(builder.filename, std::ios::in | std::ios::ate | std::ios::binary);
            int64_t wb = (builder.capacity()*m_wl+7)/8;
            size_t m_low_endpos = builder.file_offset + sizeof(m_size) + sizeof(m_wl) + 9 + wb + (wb%8 ? 8-wb%8 : 0);
            if ((size_t)out.tellp() != m_low_endpos)
                throw std::runtime_error("sd_vector_disk: bad offset calculation");

            m_high.serialize(out);
            m_high_1_select.serialize(out);
            m_high_0_select.serialize(out);

            // open m_low in read-only
            m_low = int_vector_buffer<>(builder.filename, std::ios::in, 1024 * 1024, m_wl, false,
                                        builder.file_offset + sizeof(m_size) + sizeof(m_wl));

            builder = sd_vector_disk_builder();
        }

        //! Accessing the i-th element of the original bit_vector
        /*! \param i An index i with \f$ 0 \leq i < size()  \f$.
        *   \return The i-th bit of the original bit_vector
        *   \par Time complexity
        *           \f$ \Order{t_{select0} + n/m} \f$, where m equals the number of zeros
        *    \par Remark
         *         The time complexity can be easily improved to
        *            \f$\Order{t_{select0}+\log(n/m)}\f$
        *        by using binary search in the second step.
        */
        value_type operator[](size_type i)const
        {
            size_type high_val = (i >> (m_wl));
            size_type sel_high = m_high_0_select(high_val + 1);
            size_type rank_low = sel_high - high_val;
            if (0 == rank_low)
                return 0;
            size_type val_low = i & bits::lo_set[ m_wl ]; // extract the low m_wl = log n -log m bits
            --sel_high; --rank_low;
            while (m_high[sel_high] and m_low[rank_low] > val_low) {
                if (sel_high > 0) {
                    --sel_high; --rank_low;
                } else
                    return 0;
            }
            return m_high[sel_high] and m_low[rank_low] == val_low;
        }

        //! Get the integer value of the binary string of length len starting at position idx.
        /*! \param idx Starting index of the binary representation of the integer.
         *  \param len Length of the binary representation of the integer. Default value is 64.
         *  \returns The integer value of the binary string of length len starting at position idx.
         *
         *  \pre idx+len-1 in [0..size()-1]
         *  \pre len in [1..64]
         */
        uint64_t get_int(size_type idx, const uint8_t len=64) const
        {
            uint64_t i = idx+len-1;
            uint64_t high_val = (i >> (m_wl));
            uint64_t sel_high = m_high_0_select(high_val + 1);
            uint64_t rank_low = sel_high - high_val;
            if (0 == rank_low)
                return 0;
            size_type val_low = i & bits::lo_set[ m_wl ]; // extract the low m_wl = log n -log m bits
            --sel_high; --rank_low;
            while (m_high[sel_high] and m_low[rank_low] > val_low) {
                if (sel_high > 0) {
                    --sel_high; --rank_low;
                } else
                    return 0;
            }
            uint64_t res = 0;
            while (true) {
                while (!m_high[sel_high]) {
                    if (sel_high > 0 and(high_val << m_wl) >=idx) {
                        --sel_high; --high_val;
                    } else {
                        return res;
                    }
                }
                while (m_high[sel_high]) {
                    uint64_t val = (high_val << m_wl) + m_low[rank_low];
                    if (val >= idx) {
                        res |= 1ULL<<(val-idx);
                    } else {
                        return res;
                    }
                    if (sel_high > 0) {
                        --sel_high; --rank_low;
                    } else {
                        return res;
                    }
                }
            }
        }

        //! Swap method
        void swap(sd_vector_disk& v)
        {
            if (this != &v) {
                std::swap(m_size, v.m_size);
                std::swap(m_wl, v.m_wl);
                m_low.swap(v.m_low);
                m_high.swap(v.m_high);
                util::swap_support(m_high_1_select, v.m_high_1_select, &m_high, &v.m_high);
                util::swap_support(m_high_0_select, v.m_high_0_select, &m_high, &v.m_high);
            }
        }

        //! Returns the size of the original bit vector.
        size_type size()const
        {
            return m_size;
        }

        sd_vector_disk& operator=(const sd_vector_disk& v)
        {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        sd_vector_disk& operator=(sd_vector_disk&& v)
        {
            if (this != &v) {
                m_size = v.m_size;
                m_wl   = v.m_wl;
                m_low  = std::move(v.m_low);
                m_high = std::move(v.m_high);
                m_high_1_select = std::move(v.m_high_1_select);
                m_high_1_select.set_vector(&m_high);
                m_high_0_select = std::move(v.m_high_0_select);
                m_high_0_select.set_vector(&m_high);
            }
            return *this;
        }

        //! Serializes the data structure into the given ostream
        // size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        // {
        //     throw std::runtime_error("Serialization for sd_vector_disk is not supported, use the constructor")
        // }

        //! Loads the data structure from the given istream.
        void load(const std::string &filename, size_t file_offset = 0)
        {
            std::ifstream in(filename, std::ios::binary);
            in.seekg(file_offset);
            read_member(m_size, in);
            read_member(m_wl, in);
            size_t m_low_offset = in.tellg();
            // skip m_low in file
            uint64_t m_low_size;
            uint8_t m_low_width;
            int_vector<>::read_header(m_low_size, m_low_width, in);
            int64_t wb = (m_low_size*m_wl+7)/8;
            assert(m_low_width == m_wl);
            in.seekg(m_low_offset + 9 + wb + (wb%8 ? 8-wb%8 : 0));
            m_high.load(in);
            m_high_1_select.load(in, &m_high);
            m_high_0_select.load(in, &m_high);

            m_low = int_vector_buffer<>(filename, std::ios::in, 1024 * 1024, m_wl, false, m_low_offset);
        }

        iterator begin() const
        {
            return iterator(this, 0);
        }

        iterator end() const
        {
            return iterator(this, size());
        }
};

template<uint8_t t_b>
struct rank_support_sd_disk_trait {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r,size_type)
    {
        return r;
    }
};

template<>
struct rank_support_sd_disk_trait<0> {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, size_type n)
    {
        return n - r;
    }
};

//! Rank data structure for sd_vector_disk
/*! \tparam t_b             Bit pattern.
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<uint8_t t_b, class t_hi_bit_vector, class t_select_1, class t_select_0>
class rank_support_sd_disk
{
        static_assert(t_b == 1u or t_b == 0u , "rank_support_sd_disk: bit pattern must be `0` or `1`");
    public:
        typedef bit_vector::size_type size_type;
        typedef sd_vector_disk<t_hi_bit_vector, t_select_1, t_select_0> bit_vector_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const bit_vector_type* m_v;

    public:

        explicit rank_support_sd_disk(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type rank(size_type i)const
        {
            assert(m_v != nullptr);
            assert(i <= m_v->size());
            // split problem in two parts:
            // (1) find  >=
            size_type high_val = (i >> (m_v->wl));
            size_type sel_high = m_v->high_0_select(high_val + 1);
            size_type rank_low = sel_high - high_val; //
            if (0 == rank_low)
                return rank_support_sd_disk_trait<t_b>::adjust_rank(0, i);
            size_type val_low = i & bits::lo_set[ m_v->wl ];
            // now since rank_low > 0 => sel_high > 0
            do {
                if (!sel_high)
                    return rank_support_sd_disk_trait<t_b>::adjust_rank(0, i);
                --sel_high; --rank_low;
            } while (m_v->high[sel_high] and m_v->low[rank_low] >= val_low);
            return rank_support_sd_disk_trait<t_b>::adjust_rank(rank_low+1, i);
        }

        size_type operator()(size_type i)const
        {
            return rank(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        rank_support_sd_disk& operator=(const rank_support_sd_disk& rs)
        {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(rank_support_sd_disk&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            return serialize_empty_object(out, v, name, this);
        }
};

template<uint8_t t_b, class t_sd_vec>
struct select_support_sd_disk_trait {
    typedef bit_vector::size_type size_type;
    static size_type select(size_type i, const t_sd_vec* v)
    {
        return v->low[i-1] +  // lower part of the number
               ((v->high_1_select(i) + 1 - i)  << (v->wl));  // upper part
        //^-number of 0 before the i-th 1-^    ^-shift by wl
    }
};

template<class t_sd_vec>
struct select_support_sd_disk_trait<0, t_sd_vec> {
    typedef bit_vector::size_type size_type;
    static size_type select(size_type i, const t_sd_vec* v)
    {
        auto ones  = v->low.size();
        assert(0 < i and i <= v->size() - ones);
        size_type lb = 1, rb = ones+1;
        size_type r0 = 0;
        size_type pos = (size_type)-1;
        // rb exclusive
        // invariant: rank0(select_1(rb)) >= i
        while (lb < rb) {
            auto mid = lb + (rb-lb)/2;
            auto x = select_support_sd_disk_trait<1, t_sd_vec>::select(mid, v);
            auto rank0 = x + 1 - mid;
            if (rank0 >= i) {
                rb = mid;
            } else {
                r0 = rank0;
                pos = x;
                lb = mid + 1;
            }
        }
        return pos + i - r0;
    }
};

//! Select data structure for sd_vector_disk
/*! \tparam t_b             Bit pattern.
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<uint8_t t_b, class t_hi_bit_vector, class t_select_1, class t_select_0>
class select_support_sd_disk
{
    public:
        typedef bit_vector::size_type size_type;
        typedef sd_vector_disk<t_hi_bit_vector, t_select_1, t_select_0> bit_vector_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const bit_vector_type* m_v;
    public:

        explicit select_support_sd_disk(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select(size_type i)const
        {
            return select_support_sd_disk_trait<t_b, bit_vector_type>::select(i, m_v);
        }

        size_type operator()(size_type i)const
        {
            return select(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        select_support_sd_disk& operator=(const select_support_sd_disk& ss)
        {
            if (this != &ss) {
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(select_support_sd_disk&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            return serialize_empty_object(out, v, name, this);
        }
};

} // namespace sdsl

#endif // __SD_VECTOR_DISK__
