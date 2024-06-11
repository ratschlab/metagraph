#include "bitmap.hpp"

#include "common/algorithms.hpp"
#include "bit_vector_sdsl.hpp"
#include "bit_vector_adaptive.hpp"
#include "vector_algorithm.hpp"


// since std::set needs 32 bytes per element
// switch when n < 64m + 32 * 8 * m \approx (1 << 8)m
const size_t bitmap_adaptive::kMaxNumIndicesLogRatio = 8;
const size_t bitmap_adaptive::kNumIndicesMargin = 2;
const size_t bitmap_adaptive::kRowCutoff = 1'000'000;

const size_t SPARSE_DENSE_FACTOR = 128;


////////////////////////////////////////////////////////////////
//                           bitmap                           //
////////////////////////////////////////////////////////////////

bool bitmap::operator==(const bitmap &other) const {
    if (size() != other.size() || num_set_bits() != other.num_set_bits())
        return false;

    const bitmap *bitmap_p[] = { this, &other };
    const sdsl::bit_vector *bv_p[] = { nullptr, nullptr };

    for (int i : { 0, 1 }) {
        while (dynamic_cast<const bitmap_adaptive*>(bitmap_p[i])) {
            bitmap_p[i] = &dynamic_cast<const bitmap_adaptive&>(*bitmap_p[i]).data();
        }
        if (dynamic_cast<const bitmap_vector*>(bitmap_p[i])) {
            bv_p[i] = &dynamic_cast<const bitmap_vector&>(*bitmap_p[i]).data();
        } else if (dynamic_cast<const bit_vector_stat*>(bitmap_p[i])) {
            bv_p[i] = &dynamic_cast<const bit_vector_stat&>(*bitmap_p[i]).data();
        }
    }

    if (bv_p[0] && bv_p[1]) {
        return *bv_p[0] == *bv_p[1];
    } else if (dynamic_cast<const bitmap_set*>(bitmap_p[0])
                && dynamic_cast<const bitmap_set*>(bitmap_p[1])) {
        return dynamic_cast<const bitmap_set&>(*bitmap_p[0]).size()
                == dynamic_cast<const bitmap_set&>(*bitmap_p[1]).size()
            && dynamic_cast<const bitmap_set&>(*bitmap_p[0]).data()
                == dynamic_cast<const bitmap_set&>(*bitmap_p[1]).data();
    }

    // for very sparse vectors
    if (num_set_bits() * SPARSE_DENSE_FACTOR < size()) {
        bool equal = true;
        call_ones([&other,&equal](auto i) { if (equal && !other[i]) equal = false; });
        return equal;
    }

    // TODO: implement call_zeros and add the same thing for very dense vectors

    uint64_t i;
    const uint64_t end = size();
    for (i = 0; i + 64 <= end; i += 64) {
        if (get_int(i, 64) != other.get_int(i, 64))
            return false;
    }
    if (i < size()) {
        if (get_int(i, size() - i) != other.get_int(i, size() - i))
            return false;
    }
    return true;
}

void bitmap::call_ones(const VoidCall<uint64_t> &callback) const {
    call_ones_in_range(0, size(), callback);
}


////////////////////////////////////////////////////////////////
//                       bitmap dynamic                       //
////////////////////////////////////////////////////////////////

bitmap_dyn& bitmap_dyn::operator|=(const bitmap &other) {
    assert(size() == other.size());
    other.call_ones([this](auto i) { this->set(i, true); });
    return *this;
}


////////////////////////////////////////////////////////////////
// bitmap_set: std::set<uint64_t> storing indices of set bits //
////////////////////////////////////////////////////////////////

bitmap_set::bitmap_set(uint64_t size, bool value)
      : size_(size) {
    if (value) {
        for (uint64_t i = 0; i < size; ++i) {
            bits_.insert(bits_.end(), i);
        }
    }
}

bitmap_set::bitmap_set(uint64_t size, const std::set<uint64_t> &bits)
      : size_(size), bits_(bits) {
    assert(bits_.lower_bound(size_) == bits_.end());
}

bitmap_set::bitmap_set(uint64_t size, std::initializer_list<uint64_t> init)
      : size_(size),
        bits_(init) {
    assert(bits_.lower_bound(size_) == bits_.end());
}

bitmap_set::bitmap_set(uint64_t size, std::set<uint64_t>&& bits) noexcept
      : size_(size), bits_(std::move(bits)) {
    assert(bits_.lower_bound(size_) == bits_.end());
}

void bitmap_set::set(uint64_t id, bool val) {
    if (val) {
        bits_.insert(id);
    } else {
        bits_.erase(id);
    }
}

uint64_t bitmap_set::get_int(uint64_t id, uint32_t width) const {
    sdsl::bit_vector vect(64, false);
    const auto ub = bits_.lower_bound(id + width);
    for (auto lb = bits_.lower_bound(id); lb != ub; ++lb) {
        assert(*lb >= id);
        vect[*lb - id] = true;
    }

    return *vect.data();
}

void bitmap_set::insert_zeros(const std::vector<uint64_t> &pos) {
    assert(std::is_sorted(pos.begin(), pos.end()));
    assert(!pos.size() || pos.back() < size_ + pos.size());

    size_ += pos.size();

    std::set<uint64_t> bits;
    uint64_t offset = 0;

    for (auto i : bits_) {
        while (offset < pos.size() && i + offset >= pos[offset]) {
            ++offset;
        }
        bits.emplace_hint(bits.end(), i + offset);
    }
    assert(bits.size() == bits_.size());

    bits_ = std::move(bits);
}

void bitmap_set::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());
    call_ones([other](auto i) { (*other)[i] = true; });
}

void bitmap_set::call_ones_in_range(uint64_t begin, uint64_t end,
                                    const VoidCall<uint64_t> &callback) const {
    std::for_each(bits_.lower_bound(begin), bits_.lower_bound(end), callback);
}


/////////////////////////////////////////////////////////////
// bitmap_vector: sdsl::bit_vector                         //
/////////////////////////////////////////////////////////////

bitmap_vector::bitmap_vector(uint64_t size, bool value)
      : num_set_bits_(value ? size : 0), bit_vector_(size, value) {}

bitmap_vector::bitmap_vector(const sdsl::bit_vector &vector)
      : num_set_bits_(sdsl::util::cnt_one_bits(vector)), bit_vector_(vector) {}

bitmap_vector::bitmap_vector(std::initializer_list<bool> init)
      : bitmap_vector(sdsl::bit_vector(init)) {}

bitmap_vector::bitmap_vector(sdsl::bit_vector&& vector) noexcept
      : num_set_bits_(sdsl::util::cnt_one_bits(vector)),
        bit_vector_(std::move(vector)) {}

void bitmap_vector::set(uint64_t id, bool val) {
    if (bit_vector_[id] == val)
        return;

    if (val) {
        num_set_bits_++;
    } else {
        num_set_bits_--;
    }

    bit_vector_[id] = val;
}

void bitmap_vector::insert_zeros(const std::vector<uint64_t> &pos) {
    utils::insert(&bit_vector_, pos, 0);
}

bitmap_vector& bitmap_vector::operator|=(const bitmap &other) {
    assert(size() == other.size());

    if (dynamic_cast<const bitmap_vector*>(&other)) {
        bit_vector_ |= dynamic_cast<const bitmap_vector &>(other).data();
        num_set_bits_ = sdsl::util::cnt_one_bits(bit_vector_);

    } else if (dynamic_cast<const bitmap_adaptive*>(&other)) {
        *this |= dynamic_cast<const bitmap_adaptive &>(other).data();

    } else if (dynamic_cast<const bit_vector_stat*>(&other)) {
        bit_vector_ |= dynamic_cast<const bit_vector_stat &>(other).data();
        num_set_bits_ = sdsl::util::cnt_one_bits(bit_vector_);

    } else if (dynamic_cast<const bit_vector_adaptive*>(&other)) {
        *this |= dynamic_cast<const bit_vector_adaptive &>(other).data();

    } else {
        other.call_ones([this](auto i) { this->set(i, true); });
    }

    return *this;
}

void bitmap_vector::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());
    *other |= bit_vector_;
}

void bitmap_vector::call_ones_in_range(uint64_t begin, uint64_t end,
                                       const VoidCall<uint64_t> &callback) const {
    ::call_ones(bit_vector_, begin, end, callback);
}


/////////////////////////////////////////////////////////////
// bitmap_adaptive: set and vector                         //
/////////////////////////////////////////////////////////////

bitmap_adaptive::bitmap_adaptive(uint64_t size, bool value) {
    if (value || size < kRowCutoff) {
        bitmap_.reset(new bitmap_vector(size, value));
    } else {
        bitmap_.reset(new bitmap_set(size, value));
    }
}

bitmap_adaptive::bitmap_adaptive(const sdsl::bit_vector &vector)
      : bitmap_(new bitmap_vector(vector)) {}

bitmap_adaptive::bitmap_adaptive(uint64_t size, const std::set<uint64_t> &bits)
      : bitmap_(new bitmap_set(size, bits)) {}

bitmap_adaptive::bitmap_adaptive(std::initializer_list<bool> bitmap)
      : bitmap_adaptive(sdsl::bit_vector(bitmap)) {}

bitmap_adaptive::bitmap_adaptive(uint64_t size, std::initializer_list<uint64_t> bits)
      : bitmap_adaptive(size, std::set<uint64_t>(bits)) {}

bitmap_adaptive::bitmap_adaptive(sdsl::bit_vector&& vector) noexcept
      : bitmap_(new bitmap_vector(std::move(vector))) {}

bitmap_adaptive::bitmap_adaptive(uint64_t size, std::set<uint64_t>&& bits) noexcept
      : bitmap_(new bitmap_set(size, std::move(bits))) {}

void bitmap_adaptive::insert_zeros(const std::vector<uint64_t> &pos) {
    bitmap_->insert_zeros(pos);
    check_switch();
}

bitmap_adaptive& bitmap_adaptive::operator|=(const bitmap &other) {
    assert(size() == other.size());

    if (other.num_set_bits() >= (size() >> kMaxNumIndicesLogRatio)) {
        to_bit_vector();
    }

    *bitmap_ |= other;

    check_switch();

    return *this;
}

void bitmap_adaptive::set(uint64_t id, bool val) {
    bitmap_->set(id, val);
    check_switch();
}

void bitmap_adaptive::check_switch() {
    if (size() < kRowCutoff
            || bitmap_->num_set_bits()
                    >= (bitmap_->size() >> kMaxNumIndicesLogRatio)) {
        to_bit_vector();
    } else if (size() >= kRowCutoff
            && bitmap_->num_set_bits()
                    < (bitmap_->size() >> (kMaxNumIndicesLogRatio + kNumIndicesMargin))) {
        to_set();
    }
}

void bitmap_adaptive::to_bit_vector() {
    if (dynamic_cast<bitmap_set*>(bitmap_.get())) {
        sdsl::bit_vector vector(bitmap_->size(), false);
        bitmap_->add_to(&vector);
        bitmap_.reset(new bitmap_vector(std::move(vector)));
    }
}

void bitmap_adaptive::to_set() {
    if (dynamic_cast<bitmap_vector*>(bitmap_.get())) {
        std::set<uint64_t> bits;
        bitmap_->call_ones([&](auto i) { bits.insert(bits.end(), i); });
        bitmap_.reset(new bitmap_set(bitmap_->size(), std::move(bits)));
    }
}


////////////////////////////////////////////////////////////////
//                         bitmap_lazy                        //
////////////////////////////////////////////////////////////////

bitmap_lazy::bitmap_lazy(size_t size, bool value)
      : bitmap_lazy([value](auto) { return value; }, size, size * value) {}

bitmap_lazy::bitmap_lazy(std::function<bool(uint64_t)>&& callback,
                         size_t size,
                         size_t num_set_bits) noexcept
      : is_bit_set_(std::move(callback)),
        size_(size),
        num_set_bits_(num_set_bits) {}

uint64_t bitmap_lazy::get_int(uint64_t id, uint32_t width) const {
    assert(width <= sizeof(uint64_t) * 8);

    uint64_t word = 0;
    if (!width)
        return word;

    uint64_t i = id + width - 1;
    while (i != id) {
        word = (word << 1) | operator[](i--);
    }

    return (word << 1) | operator[](id);
}

// TODO: remove num_set_bits() from bitmap and from bitmap_lazy
uint64_t bitmap_lazy::num_set_bits() const {
    if (num_set_bits_ == static_cast<size_t>(-1)) {
        throw std::runtime_error("Number of set bits not pre-computed");
    }

    return num_set_bits_;
}

void bitmap_lazy::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());
    call_ones([other](auto i) { (*other)[i] = true; });
}

void bitmap_lazy::call_ones_in_range(uint64_t begin, uint64_t end,
                                     const VoidCall<uint64_t> &callback) const {
    for (uint64_t i = begin; i < end; ++i) {
        if (operator[](i))
            callback(i);
    }
}


////////////////////////////////////////////////////////////////
//                         bitmap_generator                   //
////////////////////////////////////////////////////////////////

bitmap_generator::bitmap_generator(size_t size, bool value)
      : bitmap_generator([value,size](const auto &index_callback) {
                             if (value) {
                                 for (size_t i = 0; i < size; ++i) {
                                     index_callback(i);
                                 }
                             }
                         },
                         size, size * value) {}

bitmap_generator
::bitmap_generator(std::function<void(const VoidCall<uint64_t>&)>&& generator,
                   size_t size,
                   size_t num_set_bits) noexcept
      : size_(size), num_set_bits_(num_set_bits), generator_(std::move(generator)) {}

void bitmap_generator::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());
    call_ones([other](auto i) { (*other)[i] = true; });
}

void bitmap_generator::call_ones_in_range(uint64_t begin, uint64_t end,
                                          const VoidCall<uint64_t> &callback) const {
    if (begin || end != size_)
        throw std::runtime_error("This is inefficient and should not be called");

    generator_(callback);
}
