#ifndef __SD_VECTOR_BUILDER_DISK__
#define __SD_VECTOR_BUILDER_DISK__

#include <filesystem>
#include <type_traits>

#include <sdsl/sd_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/int_vector_mapper.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/util.hpp>


namespace sdsl {


// returns the number of bytes taken on disk by int_vector<0> of this size
template <uint8_t t_width = 0>
inline uint64_t int_vector_space(uint64_t size, uint8_t width) {
    uint64_t wb = (size * width + 7) / 8;
    return (t_width ? 8 : 9) + wb + (wb%8 ? 8-wb%8 : 0);
}

//! Class for in-place construction of sd_vector from a strictly increasing sequence
/*! \par Building an sd_vector will clear the builder.
 */
template<class t_select_1     = typename bit_vector::select_1_type,
         class t_select_0     = typename bit_vector::select_0_type>
class sd_vector_builder_disk
{
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
        int_vector_buffer<1> m_high;
        bool m_flushed = false;

    public:
        //! Constructor
        /*! \param n Vector size.
         *  \param m The number of 1-bits.
         *  \param filename The buffer filename to where the bitmap will eventually be stored.
         *  \param offset The offset position in the file where the bitmap will eventually be written.
         */
        sd_vector_builder_disk(size_type n, size_type m, const std::string &fname, size_t offset = 0) :
            m_size(n), m_capacity(m),
            m_wl(0),
            m_tail(0), m_items(0),
            m_last_high(0), m_highpos(0), filename(fname), file_offset(offset)
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
            uint64_t m_low_offset = file_offset + sizeof(m_size) + sizeof(m_wl);
            uint64_t m_high_offset = m_low_offset + int_vector_space<>(m_capacity, m_wl);
            // fill the file with 0-bytes for `m_low` and `m_high`
            std::filesystem::resize_file(filename, m_high_offset + int_vector_space<1>(m_capacity + (1ULL << logm), 1));
            // write headers for `m_low` and `m_high`
            out.open(filename, std::ios::in | std::ios::binary);
            out.seekp(m_low_offset);
            int_vector<0>::write_header(m_capacity, m_wl, out);
            out.seekp(m_high_offset);
            bit_vector::write_header(m_capacity + (1ULL << logm), 1, out);
            out.close();

            m_low = int_vector_buffer<>(filename, std::ios::in|std::ios::out, 1024 * 1024, m_wl, false, m_low_offset);
            m_high = int_vector_buffer<1>(filename, std::ios::in|std::ios::out, 1024 * 1024, 1, false, m_high_offset);
        }

        ~sd_vector_builder_disk()
        {
            if (!m_flushed)
                finish();
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

        void finish()
        {
            if (m_flushed) {
                std::runtime_error("sd_vector_disk: already flushed.");
            } else if (items() < capacity()) {
                throw std::runtime_error("sd_vector_disk: the builder is not full.");
            } else if (items() > capacity()) {
                throw std::runtime_error("sd_vector_disk: builder overflow.");
            }

            // flush to file
            m_low.close();
            m_high.close();

            // open m_low and m_high in int_vector_mapper for read-only
            uint64_t low_offset = file_offset + sizeof(m_size) + sizeof(m_wl);

            uint64_t high_offset = low_offset + int_vector_space<>(capacity(), m_wl);
            int_vector_mapper<1,std::ios_base::in> high_mapper(filename, false, false, high_offset);
            const bit_vector &high = high_mapper.wrapper();

            std::ofstream out(filename, std::ios::in | std::ios::ate | std::ios::binary);
            if ((size_t)out.tellp() != high_offset + int_vector_space<1>(high.size(), 1))
                throw std::runtime_error("sd_vector_disk: bad offset calculation");

            {
                t_select_1 high_1_select;
                util::init_support(high_1_select, &high);
                high_1_select.serialize(out);
            }
            {
                t_select_0 high_0_select;
                util::init_support(high_0_select, &high);
                high_0_select.serialize(out);
            }

            m_flushed = true;
        }
};

} // namespace sdsl

#endif // __SD_VECTOR_BUILDER_DISK__
