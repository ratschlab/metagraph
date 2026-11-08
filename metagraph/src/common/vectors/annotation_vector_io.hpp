#ifndef MTG_ANNOTATION_VECTOR_IO_HPP
#define MTG_ANNOTATION_VECTOR_IO_HPP
#include <cstdint>
#include <cstring>
#include <istream>
#include <ostream>
#include <span>
#include <stdexcept>
#include <sdsl/packed_vector_view.hpp>
#include <sdsl/shared_mapper.hpp>
namespace mtg::detail {
// Read and write the existing MetaGraph annotation-vector format:
// a 64-bit bit count, an optional width byte, then the packed payload.
// SDSL provides the shared packed view and ownership conversion.
template <uint8_t Width> class annotation_vector_io : public sdsl::packed_vector_view<Width> {
    using base = sdsl::packed_vector_view<Width>;
    using vector_type = sdsl::int_vector<Width>;

  public:
    annotation_vector_io() = default;
    using base::base;
    const base &packed() const { return *this; }
    void load(std::istream &in) {
        base next;
        if (auto *stream = dynamic_cast<sdsl::mmap_ifstream *>(&in)) {
            auto owner = stream->get_mmap_context();
            auto pos = in.tellg();
            if (pos < 0 || !owner)
                throw std::ios_base::failure("Invalid mapped vector position");
            uint64_t offset = pos, header = 8 + (Width == 0), bits;
            if (offset > owner->file_size_bytes() || header > owner->file_size_bytes() - offset)
                throw std::ios_base::failure("Packed vector header outside file");
            std::memcpy(&bits, owner->data() + offset, 8);
            uint8_t width = Width ? Width : owner->data()[offset + 8];
            std::span<const std::byte> file(reinterpret_cast<const std::byte *>(owner->data()),
                                            owner->file_size_bytes());
            next = base::from_bytes(owner, file.subspan(offset + header), bits, width);
            in.seekg(header + next.bytes().size(), std::ios_base::cur);
        } else {
            vector_type value;
            value.load(in);
            next = base(std::move(value));
        }
        if (!in)
            throw std::ios_base::failure("Cannot load packed vector view");
        static_cast<base &>(*this) = std::move(next);
    }
    void serialize(std::ostream &out) const {
        auto bits = this->bit_size();
        out.write(reinterpret_cast<const char *>(&bits), 8);
        if constexpr (Width == 0) {
            auto w = this->width();
            out.write(reinterpret_cast<const char *>(&w), 1);
        }
        auto bytes = this->bytes();
        if (!bytes.empty())
            out.write(reinterpret_cast<const char *>(bytes.data()), bytes.size());
    }
};
} // namespace mtg::detail
#endif
