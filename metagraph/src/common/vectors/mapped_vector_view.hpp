#ifndef MTG_MAPPED_VECTOR_VIEW_HPP
#define MTG_MAPPED_VECTOR_VIEW_HPP

#include <cassert>
#include <memory>
#include <ranges>
#include <sdsl/int_vector.hpp>

namespace mtg::detail {

// Share immutable packed storage; materialize at an owning boundary.
// Iterators do not retain the owner; this intentionally is not a borrowed_range.
template <uint8_t Width>
class mapped_vector_view : public std::ranges::view_interface<mapped_vector_view<Width>> {
    using vector_type = sdsl::int_vector<Width>;
    vector_type mapped_;
    std::shared_ptr<sdsl::mmap_context> context_;
    uint64_t offset_ = 0;
    std::shared_ptr<vector_type> storage_;
    const uint64_t *words_ = nullptr;

  public:
    mapped_vector_view() = default;
    mapped_vector_view(const mapped_vector_view &other)
        : context_(other.context_), offset_(other.offset_), storage_(other.storage_), words_(other.words_) {
        if (context_) mapped_.load(context_, offset_);
    }
    mapped_vector_view(mapped_vector_view &&) = default;
    mapped_vector_view &operator=(mapped_vector_view &&) = default;
    mapped_vector_view &operator=(const mapped_vector_view &other) {
        if (this != &other) *this = mapped_vector_view(other);
        return *this;
    }
    explicit mapped_vector_view(vector_type value) {
        if constexpr (Width != 0) {
            if (value.empty()) return;
        }
        // An ordinary vector argument can itself carry a read-only mapping.
        // load() validates embedded bounds for zero-copy; this path owns its value.
        storage_ = value.is_mapped() ? std::make_shared<vector_type>(value)
                                    : std::make_shared<vector_type>(std::move(value));
        words_ = storage_->data();
    }

    const vector_type &packed() const {
        if (context_) return mapped_;
        if (storage_) return *storage_;
        static const vector_type empty;
        return empty;
    }
    auto begin() const { return packed().begin(); }
    auto end() const { return packed().end(); }
    uint64_t size() const { return packed().size(); }
    uint8_t width() const { return packed().width(); }
    uint64_t operator[](uint64_t i) const {
        assert(i < size());
        const uint8_t bits = Width ? Width : width();
        return sdsl::bits::read_int(words_ + ((i * bits) >> 6), (i * bits) & 63, bits);
    }
    uint64_t get_int(uint64_t bit, uint8_t length) const {
        assert(bit + length <= packed().bit_size());
        return length ? sdsl::bits::read_int(words_ + (bit >> 6), bit & 63, length) : 0;
    }
    bool is_mapped() const { return static_cast<bool>(context_); }

    vector_type to_owned() const { return packed(); }
    vector_type into_owned() && {
        if (!context_ && storage_ && storage_.use_count() == 1) {
            auto result = std::move(*storage_);
            storage_.reset();
            words_ = nullptr;
            return result;
        }
        return to_owned();
    }

    void load(std::istream &in) {
        mapped_vector_view next;
        if (auto *stream = dynamic_cast<sdsl::mmap_ifstream *>(&in)) {
            const auto position = in.tellg();
            if (position < 0) throw std::ios_base::failure("Invalid mapped vector position");
            next.context_ = stream->get_mmap_context();
            next.offset_ = static_cast<uint64_t>(position);
            const auto bytes = next.mapped_.load(next.context_, next.offset_);
            in.seekg(bytes, std::ios_base::cur);
        } else {
            next.storage_ = std::make_shared<vector_type>();
            next.storage_->load(in);
        }
        if (!in) throw std::ios_base::failure("Cannot load packed vector view");
        next.words_ = next.packed().data();
        *this = std::move(next);
    }
    void serialize(std::ostream &out) const { packed().serialize(out); }
};

static_assert(std::ranges::view<mapped_vector_view<1>>);
static_assert(std::ranges::sized_range<const mapped_vector_view<1>>);
static_assert(!std::ranges::borrowed_range<mapped_vector_view<1>>);

} // namespace mtg::detail
#endif
