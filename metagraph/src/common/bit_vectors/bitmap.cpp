#include "bitmap.hpp"
#include "utils.hpp"


// since std::set needs 32 bytes per element
// switch when n < 64m + 32 * 8 * m \approx (1 << 8)m
const size_t bitmap_adaptive::kMaxNumIndicesLogRatio = 8;
const size_t bitmap_adaptive::kNumIndicesMargin = 2;
const size_t bitmap_adaptive::kRowCutoff = 1'000'000;


sdsl::bit_vector to_sdsl(const std::vector<bool> &vector) {
    sdsl::bit_vector result(vector.size(), 0);
    for (size_t i = 0; i < vector.size(); ++i) {
        if (vector[i])
            result[i] = 1;
    }
    return result;
}

void call_ones(const sdsl::bit_vector &vector,
               uint64_t begin, uint64_t end,
               const VoidCall<uint64_t> &callback) {
    assert(begin <= end);
    assert(end <= vector.size());

    if (sdsl::util::cnt_one_bits(vector) > vector.size() / 2) {
        //TODO: benchmark to check if this actually makes it faster
        for (uint64_t i = begin; i < end; ++i) {
            if (vector[i])
                callback(i);
        }
        return;
    }

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

void call_ones(const sdsl::bit_vector &vector,
               const VoidCall<uint64_t> &callback) {
    call_ones(vector, 0, vector.size(), callback);
}

void call_zeros(const sdsl::bit_vector &vector,
                uint64_t begin, uint64_t end,
                const VoidCall<uint64_t> &callback) {
    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for (; i < end && i & 0x3F; ++i) {
        if (vector[i])
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

void call_zeros(const sdsl::bit_vector &vector,
               const VoidCall<uint64_t> &callback) {
    call_zeros(vector, 0, vector.size(), callback);
}

uint64_t inner_prod(const sdsl::bit_vector &first,
                    const sdsl::bit_vector &second) {
    assert(first.size() == second.size());

    if (first.empty())
        return 0;

    const uint64_t *first_data = first.data();
    const uint64_t *second_data = second.data();

    uint64_t result = sdsl::bits::cnt((*first_data) & (*second_data));

    for (typename sdsl::bit_vector::size_type i = 1; i < (first.capacity() >> 6); ++i) {
        result += sdsl::bits::cnt(*(++first_data) & *(++second_data));
    }
    if (first.bit_size() & 0x3F) {
        result -= sdsl::bits::cnt((*first_data) & (*second_data)
                                    & (~sdsl::bits::lo_set[first.bit_size() & 0x3F]));
    }
    return result;
}


////////////////////////////////////////////////////////////////
//                           bitmap                           //
////////////////////////////////////////////////////////////////

void bitmap::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());
    call_ones([other](auto i) { (*other)[i] = true; });
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
    size_ += pos.size();
    utils::insert_default_values(pos, &bits_);
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
    utils::insert_default_values(pos, &bit_vector_);
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
        bitmap_->call_ones([&](auto i) { vector[i] = true; });
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
