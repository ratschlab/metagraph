#include "bitmap.hpp"
#include "utils.hpp"


// since std::set needs 32 bytes per element
// switch when n < 64m + 32 * 8 * m \approx (1 << 8)m
const size_t bitmap_adaptive::kMaxNumIndicesLogRatio = 8;
const size_t bitmap_adaptive::kRowCutoff = 1'000'000;


void call_ones(const sdsl::bit_vector &vector,
               const std::function<void(uint64_t)> &callback) {
    uint64_t j = 64;
    uint64_t i = 0;
    uint64_t word;
    for (; j <= vector.size(); j += 64) {
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
    for (; i < vector.size(); ++i) {
        if (vector[i])
            callback(i);
    }
}

void call_zeros(const sdsl::bit_vector &vector,
                const std::function<void(uint64_t)> &callback) {
    uint64_t j = 64;
    uint64_t i = 0;
    uint64_t word;
    for (; j <= vector.size(); j += 64) {
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
    for (; i < vector.size(); ++i) {
        if (!vector[i])
            callback(i);
    }
}

uint64_t count_num_set_bits(const sdsl::bit_vector &vector) {
    uint64_t count = 0;
    uint64_t i = 0;
    for (; i + 64 <= vector.size(); i += 64) {
        count += sdsl::bits::cnt(vector.get_int(i));
    }
    for (; i < vector.size(); ++i) {
        if (vector[i])
            count++;
    }
    return count;
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

bool bitmap_set::operator[](uint64_t id) const {
    return bits_.count(id);
}

uint64_t bitmap_set::get_int(uint64_t id, uint32_t width) const {
    sdsl::bit_vector vect(64, false);
    const auto ub = bits_.upper_bound(id + width);
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

void bitmap_set::call_ones(const std::function<void(uint64_t)> &callback) const {
    std::for_each(bits_.begin(), bits_.end(), callback);
}


/////////////////////////////////////////////////////////////
// bitmap_vector: sdsl::bit_vector                         //
/////////////////////////////////////////////////////////////

void bitmap_vector::call_ones(const std::function<void(uint64_t)> &callback) const {
    ::call_ones(bit_vector_, callback);
}

bitmap_vector
::bitmap_vector(uint64_t size, bool value)
      : num_set_bits_(value ? size : 0), bit_vector_(size, value) {}

bitmap_vector
::bitmap_vector(const sdsl::bit_vector &vector)
      : num_set_bits_(count_num_set_bits(vector)), bit_vector_(vector) {}

bitmap_vector
::bitmap_vector(std::initializer_list<bool> init)
      : bitmap_vector(sdsl::bit_vector(init)) {}

bitmap_vector
::bitmap_vector(sdsl::bit_vector&& vector) noexcept
      : num_set_bits_(count_num_set_bits(vector)),
        bit_vector_(std::move(vector)) {}

void bitmap_vector::set(uint64_t id, bool val) {
    if (bit_vector_[id] != val) {
        if (val) {
            num_set_bits_++;
        } else {
            num_set_bits_--;
        }

        bit_vector_[id] = val;
    }
}

void bitmap_vector::insert_zeros(const std::vector<uint64_t> &pos) {
    utils::insert_default_values(pos, &bit_vector_);
}


/////////////////////////////////////////////////////////////
// bitmap_adaptive: set and vector                         //
/////////////////////////////////////////////////////////////

bitmap_adaptive
::bitmap_adaptive(uint64_t size, bool value) {
    if (value || size < kRowCutoff) {
        bitmap_.reset(new bitmap_vector(size, value));
    } else {
        bitmap_.reset(new bitmap_set(size, value));
    }
}

bitmap_adaptive
::bitmap_adaptive(const sdsl::bit_vector &vector)
      : bitmap_(new bitmap_vector(vector)) {}

bitmap_adaptive
::bitmap_adaptive(uint64_t size, const std::set<uint64_t> &bits)
      : bitmap_(new bitmap_set(size, bits)) {}

bitmap_adaptive
::bitmap_adaptive(std::initializer_list<bool> bitmap)
      : bitmap_adaptive(sdsl::bit_vector(bitmap)) {}

bitmap_adaptive
::bitmap_adaptive(uint64_t size, std::initializer_list<uint64_t> bits)
      : bitmap_adaptive(size, std::set<uint64_t>(bits)) {}

bitmap_adaptive
::bitmap_adaptive(sdsl::bit_vector&& vector) noexcept
      : bitmap_(new bitmap_vector(std::move(vector))) {}

bitmap_adaptive
::bitmap_adaptive(uint64_t size, std::set<uint64_t>&& bits) noexcept
      : bitmap_(new bitmap_set(size, std::move(bits))) {}

void bitmap_adaptive::insert_zeros(const std::vector<uint64_t> &pos) {
    bitmap_->insert_zeros(pos);
    check_switch();
}

void bitmap_adaptive::set(uint64_t id, bool val) {
    bitmap_->set(id, val);
    check_switch();
}

void bitmap_adaptive::check_switch() {
    if (size() < kRowCutoff
        || bitmap_->num_set_bits() >= (bitmap_->size() >> kMaxNumIndicesLogRatio)) {
        to_bit_vector();
    } else {
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
