#include "bitmap.hpp"
#include "utils.hpp"


// since std::set needs 32 bytes per element
// switch when n < 64m + 32 * 8 * m \approx (1 << 8)m
const size_t kMaxNumIndicesLogRatio = 8;
const size_t kRowCutoff = 1'000'000;

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
            bits_.emplace_hint(bits_.end(), i);
        }
    }
}

bitmap_set::bitmap_set(uint64_t size, const std::set<uint64_t> &bits)
      : size_(size), bits_(bits) {}

bitmap_set::bitmap_set(const bitmap_set &bits) : size_(bits.size_), bits_(bits.bits_) {}

bitmap_set::bitmap_set(uint64_t size, std::initializer_list<uint64_t> init)
      : size_(size),
        bits_(init) { }

bitmap_set::bitmap_set(uint64_t size, std::set<uint64_t>&& bits) noexcept
      : size_(size), bits_(std::move(bits)) {}

bitmap_set::bitmap_set(bitmap_set&& bits) noexcept
      : size_(bits.size_), bits_(std::move(bits.bits_)) {}

void bitmap_set::set(uint64_t id, bool val) {
    auto find = bits_.find(id);
    if (find == bits_.end()) {
        if (val)
            bits_.insert(id);
    } else {
        if (!val)
            bits_.erase(find);
    }
}

bool bitmap_set::operator[](uint64_t id) const {
    return bits_.find(id) != bits_.end();
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
// bitmap_vector: sdsl::bit_vector with mutexes //
/////////////////////////////////////////////////////////////

void bitmap_vector::call_ones(const std::function<void(uint64_t)> &callback) const {
    ::call_ones(bit_vector_, callback);
}

void bitmap_vector::reset_mutex_pool(uint64_t pool_size) {
    bit_vector_mutexes_.clear();
    bit_vector_mutexes_.reserve(pool_size);
    for (size_t i = 0; i < pool_size; ++i) {
        bit_vector_mutexes_.emplace_back(new std::mutex());
    }
}

void bitmap_vector::call_critical(const std::function<void(bitmap_vector *vector)> &callback) {
    // lock all mutexes
    std::vector<std::unique_ptr<std::lock_guard<std::mutex>>> vector_locks;
    vector_locks.reserve(bit_vector_mutexes_.size());
    for (const auto &mutex : bit_vector_mutexes_) {
        assert(mutex.get());
        vector_locks.emplace_back(new std::lock_guard<std::mutex>(*mutex));
    }

    callback(this);
}

void bitmap_vector::call_critical(const std::function<void(const bitmap_vector *vector)> &callback) const {
    // lock all mutexes
    std::vector<std::unique_ptr<std::lock_guard<std::mutex>>> vector_locks;
    vector_locks.reserve(bit_vector_mutexes_.size());
    for (const auto &mutex : bit_vector_mutexes_) {
        assert(mutex.get());
        vector_locks.emplace_back(new std::lock_guard<std::mutex>(*mutex));
    }

    callback(this);
}

bitmap_vector
::bitmap_vector(uint64_t size, bool value, uint64_t pool_size)
      : num_set_bits_(value ? size : 0),
        bit_vector_(size, value) {
    reset_mutex_pool(pool_size);
}

bitmap_vector
::bitmap_vector(const sdsl::bit_vector &vector, uint64_t pool_size)
      : num_set_bits_(count_num_set_bits(vector)),
        bit_vector_(vector) {
    reset_mutex_pool(pool_size);
}

bitmap_vector
::bitmap_vector(const bitmap_vector &vector) {
    vector.call_critical([&](const auto *vector) {
        num_set_bits_ = vector->num_set_bits_;
        bit_vector_ = vector->bit_vector_;
    });
    reset_mutex_pool(vector.bit_vector_mutexes_.size());
}

bitmap_vector
::bitmap_vector(std::initializer_list<bool> init, uint64_t pool_size)
      : bitmap_vector(sdsl::bit_vector(init), pool_size) {}

bitmap_vector
::bitmap_vector(sdsl::bit_vector&& vector, uint64_t pool_size) noexcept
      : num_set_bits_(count_num_set_bits(vector)),
        bit_vector_(std::move(vector)) {
    reset_mutex_pool(pool_size);
}

bitmap_vector
::bitmap_vector(bitmap_vector&& vector) noexcept {
    vector.call_critical([&](auto *vector) {
        num_set_bits_ = vector->num_set_bits_;
        bit_vector_ = std::move(vector->bit_vector_);
    });
    reset_mutex_pool(vector.bit_vector_mutexes_.size());
}

void bitmap_vector::set(uint64_t id, bool val) {
    if (bit_vector_[id] != val) {
        assert(bit_vector_mutexes_[id % bit_vector_mutexes_.size()].get());
        std::lock_guard<std::mutex> lock(
            *bit_vector_mutexes_[id % bit_vector_mutexes_.size()]
        );

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
::bitmap_adaptive(uint64_t size, bool value, uint64_t pool_size)
      : pool_size_(pool_size) {
    if (value) {
        bitmap_.reset(new bitmap_vector(size, value, pool_size_));
    } else {
        bitmap_.reset(new bitmap_set(size, value));
    }
}

bitmap_adaptive
::bitmap_adaptive(const sdsl::bit_vector &vector, uint64_t pool_size)
      : pool_size_(pool_size),
        bitmap_(new bitmap_vector(vector, pool_size_)) {}

bitmap_adaptive
::bitmap_adaptive(uint64_t size,
                  const std::set<uint64_t> &bits,
                  uint64_t pool_size)
      : pool_size_(pool_size),
        bitmap_(new bitmap_set(size, bits)) {}

bitmap_adaptive
::bitmap_adaptive(std::initializer_list<bool> init, uint64_t pool_size)
      : bitmap_adaptive(sdsl::bit_vector(init), pool_size) {}

bitmap_adaptive
::bitmap_adaptive(uint64_t size,
                  std::initializer_list<uint64_t> init,
                  uint64_t pool_size)
      : bitmap_adaptive(size, std::set<uint64_t>(init), pool_size) {}

bitmap_adaptive
::bitmap_adaptive(sdsl::bit_vector&& vector, uint64_t pool_size) noexcept
      : pool_size_(pool_size),
        bitmap_(new bitmap_vector(std::move(vector), pool_size)) {}

bitmap_adaptive
::bitmap_adaptive(uint64_t size,
                  std::set<uint64_t>&& bits,
                  uint64_t pool_size) noexcept
      : pool_size_(pool_size),
        bitmap_(new bitmap_set(size, std::move(bits))) {}

void bitmap_adaptive::set(uint64_t id, bool val) {
    if (val
        && bitmap_->size() > kRowCutoff
        && bitmap_->num_set_bits() >= (bitmap_->size() >> kMaxNumIndicesLogRatio))
        to_bit_vector();

    bitmap_->set(id, val);
}

void bitmap_adaptive::to_bit_vector() {
    if (dynamic_cast<bitmap_set*>(bitmap_.get())) {
        sdsl::bit_vector vector(bitmap_->size(), false);
        bitmap_->call_ones([&](auto i) { vector[i] = true; });
        bitmap_.reset(new bitmap_vector(std::move(vector), pool_size_));
    }
}

void bitmap_adaptive::to_set() {
    if (dynamic_cast<bitmap_vector*>(bitmap_.get())) {
        std::set<uint64_t> bits;
        bitmap_->call_ones([&](auto i) { bits.emplace_hint(bits.end(), i); });
        bitmap_.reset(new bitmap_set(bitmap_->size(), std::move(bits)));
    }
}
