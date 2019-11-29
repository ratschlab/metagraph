#include "bit_vector.hpp"

#include <cmath>
#include <type_traits>
#include <cassert>

#include "utils/serialization.hpp"

const double STAT_BITS_PER_BIT_IF_SPARSE = 1.07;
// TODO: why is this so large? Check the sdsl implementation
const double STAT_BITS_PER_BIT_IF_DENSE = 1.51;

// TODO: run benchmarks and optimize these parameters
const size_t MAX_ITER_BIT_VECTOR_STAT = 1000;
const size_t MAX_ITER_BIT_VECTOR_DYN = 0;
const size_t MAX_ITER_BIT_VECTOR_SD = 10;
const size_t MAX_ITER_BIT_VECTOR_RRR = 10;


uint64_t bit_vector::rank0(uint64_t id) const {
    return std::min(id + 1, size()) - rank1(id);
}

sdsl::bit_vector bit_vector::to_vector() const {
    sdsl::bit_vector result(size(), false);
    call_ones([&result](auto i) { result[i] = true; });
    return result;
}

std::ostream& operator<<(std::ostream &os, const bit_vector &bv) {
    return os << bv.to_vector();
}

template <class Vector>
Vector bit_vector::convert_to() {
    static_assert(!std::is_same<Vector, bit_vector_smart>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_small>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_adaptive>::value, "");

    if (dynamic_cast<Vector*>(this)) {
        // to the same type, no conversion
        return dynamic_cast<Vector&&>(*this);

    } else if (dynamic_cast<bit_vector_stat*>(this)) {
        // stat -> anything else
        return Vector(std::move(dynamic_cast<bit_vector_stat*>(this)->vector_));

    } else if (dynamic_cast<bit_vector_adaptive*>(this)) {
        // adaptive(x) -> anything
        return dynamic_cast<bit_vector_adaptive*>(this)->vector_->convert_to<Vector>();

    } else {
        // anything -> anything (slower: with full reconstruction)
        return Vector(to_vector());
    }
}
template bit_vector_dyn bit_vector::convert_to<bit_vector_dyn>();
template bit_vector_stat bit_vector::convert_to<bit_vector_stat>();
template bit_vector_sd bit_vector::convert_to<bit_vector_sd>();
template bit_vector_rrr<3> bit_vector::convert_to<bit_vector_rrr<3>>();
template bit_vector_rrr<8> bit_vector::convert_to<bit_vector_rrr<8>>();
template bit_vector_rrr<15> bit_vector::convert_to<bit_vector_rrr<15>>();
template bit_vector_rrr<31> bit_vector::convert_to<bit_vector_rrr<31>>();
template bit_vector_rrr<63> bit_vector::convert_to<bit_vector_rrr<63>>();
template bit_vector_rrr<127> bit_vector::convert_to<bit_vector_rrr<127>>();
template bit_vector_rrr<255> bit_vector::convert_to<bit_vector_rrr<255>>();
template sdsl::bit_vector bit_vector::convert_to<sdsl::bit_vector>();
template<> bit_vector_small bit_vector::convert_to() {
    return bit_vector_small(std::move(*this));
}
template<> bit_vector_smart bit_vector::convert_to() {
    return bit_vector_smart(std::move(*this));
}

template <class Vector>
Vector bit_vector::copy_to() const {
    static_assert(!std::is_same<Vector, bit_vector_smart>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_small>::value, "");
    static_assert(!std::is_same<Vector, bit_vector_adaptive>::value, "");

    if (dynamic_cast<const Vector*>(this)) {
        // copy to the same type, no conversion
        return Vector(dynamic_cast<const Vector&>(*this));

    } else if (dynamic_cast<const bit_vector_stat*>(this)) {
        // copy stat -> anything else
        auto bv = dynamic_cast<const bit_vector_stat*>(this)->vector_;
        return Vector(std::move(bv));

    } else if (dynamic_cast<const bit_vector_adaptive*>(this)) {
        // copy adaptive(x) -> anything
        return dynamic_cast<const bit_vector_adaptive*>(this)->vector_->copy_to<Vector>();

    } else {
        // anything -> anything (slower: with full reconstruction)
        return Vector(to_vector());
    }
}
template bit_vector_dyn bit_vector::copy_to<bit_vector_dyn>() const;
template bit_vector_stat bit_vector::copy_to<bit_vector_stat>() const;
template bit_vector_sd bit_vector::copy_to<bit_vector_sd>() const;
template bit_vector_rrr<3> bit_vector::copy_to<bit_vector_rrr<3>>() const;
template bit_vector_rrr<8> bit_vector::copy_to<bit_vector_rrr<8>>() const;
template bit_vector_rrr<15> bit_vector::copy_to<bit_vector_rrr<15>>() const;
template bit_vector_rrr<31> bit_vector::copy_to<bit_vector_rrr<31>>() const;
template bit_vector_rrr<63> bit_vector::copy_to<bit_vector_rrr<63>>() const;
template bit_vector_rrr<127> bit_vector::copy_to<bit_vector_rrr<127>>() const;
template bit_vector_rrr<255> bit_vector::copy_to<bit_vector_rrr<255>>() const;
template sdsl::bit_vector bit_vector::copy_to<sdsl::bit_vector>() const;
template<> bit_vector_small bit_vector::copy_to() const {
    return bit_vector_small(*this);
}
template<> bit_vector_smart bit_vector::copy_to() const {
    return bit_vector_smart(*this);
}



// taken from https://github.com/xxsds/sdsl-lite/blob/master/include/sdsl/util.hpp
// this function has been modified.
template <class t_int_vec>
typename t_int_vec::size_type
next_bit(const t_int_vec &v,
         uint64_t idx,
         uint64_t max_steps = std::numeric_limits<uint64_t>::max()) {
    assert(idx < v.bit_size());

    uint64_t pos  = idx >> 6;
    uint64_t node = v.data()[pos];
    node >>= (idx & 0x3F);
    if (node)
        return std::min(idx + sdsl::bits::lo(node), v.bit_size());

    uint64_t end = idx + std::min(max_steps, v.bit_size() - idx);
    for (++pos; (pos << 6) < end; ++pos) {
        if (v.data()[pos])
            return std::min((pos << 6) | sdsl::bits::lo(v.data()[pos]), v.bit_size());
    }
    return v.bit_size();
}

// taken from https://github.com/xxsds/sdsl-lite/blob/master/include/sdsl/util.hpp
// this function has been modified.
template <class t_int_vec>
typename t_int_vec::size_type
prev_bit(const t_int_vec &v,
         uint64_t idx,
         uint64_t max_steps = std::numeric_limits<uint64_t>::max()) {
    assert(idx < v.bit_size());

    int64_t pos  = idx >> 6;
    uint64_t node = v.data()[pos];
    node <<= 63 - (idx & 0x3F);
    if (node)
        return idx - (63 - sdsl::bits::hi(node));

    // last position to visit: 0 or (idx + 1 - max_steps)
    int64_t r_end_word = ((idx + 1 - std::min(idx + 1, max_steps)) >> 6) - 1;
    assert(r_end_word >= -1);
    for (--pos; pos > r_end_word; --pos) {
        if (v.data()[pos])
            return (pos << 6) | sdsl::bits::hi(v.data()[pos]);
    }
    return v.bit_size();
}

template <typename BitVector>
inline uint64_t next1(const BitVector &v,
                      uint64_t pos,
                      size_t num_steps) {
    assert(pos < v.size());

    if (v[pos])
        return pos;

    for (size_t t = 1; t < num_steps; ++t) {
        if (pos + t == v.size() || v[pos + t])
            return pos + t;
    }

    uint64_t rk = v.rank1(pos) + 1;
    return rk <= v.num_set_bits() ? v.select1(rk) : v.size();
}

template <typename BitVector>
inline uint64_t prev1(const BitVector &v,
                      uint64_t pos,
                      size_t num_steps) {
    assert(pos < v.size());

    for (size_t t = 0; t < num_steps; ++t, --pos) {
        if (v[pos])
            return pos;

        if (pos == 0)
            return v.size();
    }

    uint64_t rk = v.rank1(pos);
    return rk ? v.select1(rk) : v.size();
}


/////////////////////////////
// bit_vector_dyn, DYNAMIC //
/////////////////////////////

bit_vector_dyn::bit_vector_dyn(uint64_t size, bool value) {
    uint64_t i;
    for (i = 0; i + 64 <= size; i += 64) {
        vector_.push_word(value ? ~uint64_t(0) : 0, 64);
    }
    if (i < size)
        vector_.push_word(value ? (uint64_t(1) << (size - i)) - 1 : 0, size - i);
}

bit_vector_dyn::bit_vector_dyn(const sdsl::bit_vector &v) {
    uint64_t i;
    for (i = 0; i + 64 <= v.size(); i += 64) {
        vector_.push_word(v.get_int(i, 64), 64);
    }
    if (i < v.size())
        vector_.push_word(v.get_int(i, v.size() - i), v.size() - i);
}

bit_vector_dyn::bit_vector_dyn(std::initializer_list<bool> init)
      : bit_vector_dyn(sdsl::bit_vector(init)) {}

std::unique_ptr<bit_vector> bit_vector_dyn::copy() const {
    return std::make_unique<bit_vector_dyn>(*this);
}

uint64_t bit_vector_dyn::rank1(uint64_t id) const {
    return vector_.rank1(std::min(id + 1, size()));
}

uint64_t bit_vector_dyn::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= rank1(size() - 1));
    return vector_.select1(id - 1);
}

uint64_t bit_vector_dyn::select0(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= rank0(size() - 1));
    return vector_.select0(id - 1);
}

uint64_t bit_vector_dyn::next1(uint64_t pos) const {
    assert(pos < size());

    return ::next1(*this, pos, MAX_ITER_BIT_VECTOR_DYN);
}

uint64_t bit_vector_dyn::prev1(uint64_t pos) const {
    assert(pos < size());

    return ::prev1(*this, pos, MAX_ITER_BIT_VECTOR_DYN);
}

void bit_vector_dyn::set(uint64_t id, bool val) {
    vector_.set(id, val);
}

bool bit_vector_dyn::operator[](uint64_t id) const {
    assert(id < size());
    return vector_.at(id);
}

uint64_t bit_vector_dyn::get_int(uint64_t id, uint32_t width) const {
    assert(id + width <= size());
    assert(width);
    uint64_t word = 0;
    for (int64_t pos = id + width - 1; pos >= static_cast<int64_t>(id); --pos) {
        word = (word << 1) + vector_.at(pos);
    }
    return word;
}

void bit_vector_dyn::insert_bit(uint64_t id, bool val) {
    assert(id <= size());
    vector_.insert(id, val);
}

void bit_vector_dyn::delete_bit(uint64_t id) {
    assert(id < size());
    vector_.remove(id);
}

bool bit_vector_dyn::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        //TODO: catch reading errors
        vector_.load(in);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_sd." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void bit_vector_dyn::serialize(std::ostream &out) const {
    vector_.serialize(out);
}

void bit_vector_dyn::call_ones_in_range(uint64_t begin, uint64_t end,
                                        const VoidCall<uint64_t> &callback) const {
    assert(begin <= end);
    assert(end <= size());

    if (2 * num_set_bits() < size()) {
        // sparse
        uint64_t num_ones = end ? rank1(end - 1) : 0;
        for (uint64_t r = begin ? rank1(begin - 1) + 1 : 1; r <= num_ones; ++r) {
            callback(select1(r));
        }
    } else {
        // dense
        uint64_t one_pos = 0;
        uint64_t zero_pos = 0;
        uint64_t num_zeros = end ? rank0(end - 1) : 0;
        for (uint64_t r = begin ? rank0(begin - 1) + 1 : 1; r <= num_zeros; ++r) {
            zero_pos = select0(r);
            while (one_pos < zero_pos) {
                callback(one_pos++);
            }
            one_pos++;
        }
        while (one_pos < end) {
            callback(one_pos++);
        }
    }
}


///////////////////////////////////////////////
// bit_vector_stat, sdsl rank-select support //
///////////////////////////////////////////////

bit_vector_stat::bit_vector_stat(uint64_t size, bool value)
      : vector_(size, value) {
    if (value)
        num_set_bits_ = size;
}

bit_vector_stat::bit_vector_stat(const sdsl::bit_vector &vector) noexcept
      : bit_vector_stat(sdsl::bit_vector(vector)) {}

bit_vector_stat::bit_vector_stat(const sdsl::bit_vector &vector, uint64_t num_set_bits)
      : bit_vector_stat(sdsl::bit_vector(vector), num_set_bits) {}

bit_vector_stat::bit_vector_stat(const bit_vector_stat &other) {
    *this = other;
}

bit_vector_stat
::bit_vector_stat(const std::function<void(const VoidCall<uint64_t>&)> &call_ones,
                  uint64_t size)
      : vector_(size, false),
        num_set_bits_(0) {
    call_ones([&](uint64_t pos) {
        assert(pos < size);
        vector_[pos] = true;
        num_set_bits_++;
    });
}

bit_vector_stat::bit_vector_stat(sdsl::bit_vector&& vector) noexcept
      : vector_(std::move(vector)),
        num_set_bits_(sdsl::util::cnt_one_bits(vector_)) {}

bit_vector_stat::bit_vector_stat(sdsl::bit_vector&& vector, uint64_t num_set_bits)
      : vector_(std::move(vector)),
        num_set_bits_(num_set_bits) {
    assert(num_set_bits_ == sdsl::util::cnt_one_bits(vector_));
}

bit_vector_stat::bit_vector_stat(bit_vector_stat&& other) noexcept {
    *this = std::move(other);
}

bit_vector_stat::bit_vector_stat(std::initializer_list<bool> init)
      : bit_vector_stat(sdsl::bit_vector(init)) {}

bit_vector_stat& bit_vector_stat::operator=(const bit_vector_stat &other) {
    vector_ = other.vector_;
    num_set_bits_ = other.num_set_bits_;
    if (!other.requires_update_) {
        rk_ = other.rk_;
        rk_.set_vector(&vector_);
        slct_ = other.slct_;
        slct_.set_vector(&vector_);
        requires_update_ = false;
    }
    return *this;
}

bit_vector_stat& bit_vector_stat::operator=(bit_vector_stat&& other) noexcept {
    vector_ = std::move(other.vector_);
    num_set_bits_ = other.num_set_bits_;
    if (!other.requires_update_) {
        rk_ = std::move(other.rk_);
        rk_.set_vector(&vector_);
        slct_ = std::move(other.slct_);
        slct_.set_vector(&vector_);
        requires_update_ = false;
    }
    return *this;
}

std::unique_ptr<bit_vector> bit_vector_stat::copy() const {
    return std::make_unique<bit_vector_stat>(*this);
}

uint64_t bit_vector_stat::rank1(uint64_t id) const {
    if (requires_update_)
        init_rs();
    //the rank method in SDSL does not include id in the count
    return rk_(id >= this->size() ? this->size() : id + 1);
}

uint64_t bit_vector_stat::select1(uint64_t id) const {
    assert(id > 0 && size() > 0);

    if (requires_update_)
        init_rs();

    assert(id <= num_set_bits_);
    return slct_(id);
}

uint64_t bit_vector_stat::next1(uint64_t pos) const {
    assert(pos < size());

    auto next = next_bit(vector_, pos, MAX_ITER_BIT_VECTOR_STAT);
    if (next < vector_.size())
        return next;

    if (vector_.size() - pos <= MAX_ITER_BIT_VECTOR_STAT)
        return size();

    uint64_t rk = rank1(pos) + 1;
    return rk <= num_set_bits() ? select1(rk) : size();
}

uint64_t bit_vector_stat::prev1(uint64_t pos) const {
    assert(pos < size());

    auto prev = prev_bit(vector_, pos, MAX_ITER_BIT_VECTOR_STAT);
    if (prev <= pos)
        return prev;

    if (pos < MAX_ITER_BIT_VECTOR_STAT)
        return size();

    uint64_t rk = rank1(pos);
    return rk ? select1(rk) : size();
}

bool bit_vector_stat::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
}

uint64_t bit_vector_stat::get_int(uint64_t id, uint32_t width) const {
    return vector_.get_int(id, width);
}

void bit_vector_stat::set(uint64_t id, bool val) {
    if (vector_[id] == val)
        return;

    if (val) {
        num_set_bits_++;
    } else {
        num_set_bits_--;
    }

    vector_[id] = val;
    requires_update_ = true;
}

void bit_vector_stat::insert_bit(uint64_t id, bool val) {
    assert(id <= size());

    if (val)
        num_set_bits_++;

    vector_.resize(size() + 1);
    std::copy_backward(vector_.begin() + id, vector_.end() - 1, vector_.end());

    vector_[id] = val;
    requires_update_ = true;
}

void bit_vector_stat::delete_bit(uint64_t id) {
    assert(size() > 0);
    assert(id < size());

    if (vector_[id])
        num_set_bits_--;

    std::copy(vector_.begin() + id + 1, vector_.end(), vector_.begin() + id);

    vector_.resize(vector_.size() - 1);
    requires_update_ = true;
}

void bit_vector_stat::serialize(std::ostream &out) const {
    vector_.serialize(out);

    if (requires_update_)
        init_rs();

    serialize_number(out, num_set_bits_);
    rk_.serialize(out);
    slct_.serialize(out);
}

bool bit_vector_stat::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);

        try {
            num_set_bits_ = load_number(in);
            rk_.load(in, &vector_);
            slct_.load(in, &vector_);
            requires_update_ = false;
        } catch (...) {
            std::cerr << "Warning: Loading from file without bit_vector rank"
                      << " and select support dumped. Reserialize to"
                      << " make the loading faster." << std::endl;

            num_set_bits_ = sdsl::util::cnt_one_bits(vector_);
            requires_update_ = true;
            init_rs();
        }
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_stat." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void bit_vector_stat::add_to(sdsl::bit_vector *other) const {
    assert(other);
    assert(other->size() == size());
    *other |= vector_;
}

void bit_vector_stat::call_ones_in_range(uint64_t begin, uint64_t end,
                                         const VoidCall<uint64_t> &callback) const {
    assert(begin <= end);
    assert(end <= size());
    ::call_ones(vector_, begin, end, callback);
}

void bit_vector_stat::init_rs() const {
    std::unique_lock<std::mutex> lock(mu_);

    if (!requires_update_)
        return;

    rk_ = sdsl::rank_support_v5<>(&vector_);
    slct_ = sdsl::select_support_mcl<>(&vector_);

    requires_update_ = false;

    assert(num_set_bits_ == (size() ? rank1(size() - 1) : 0));
    assert(num_set_bits_ == num_set_bits());
}


////////////////////////////////////////////////////////////////
// bit_vector_sd, sdsl compressed with rank-select support    //
////////////////////////////////////////////////////////////////

bit_vector_sd::bit_vector_sd(uint64_t size, bool value)
      : inverted_(value && size) {
    sdsl::sd_vector_builder builder(size, 0);
    vector_ = decltype(vector_)(builder);
    rk1_ = decltype(rk1_)(&vector_);
    slct1_ = decltype(slct1_)(&vector_);
    slct0_ = decltype(slct0_)(&vector_);
}

bit_vector_sd::bit_vector_sd(const sdsl::bit_vector &vector)
      : bit_vector_sd(vector, sdsl::util::cnt_one_bits(vector)) {}

bit_vector_sd::bit_vector_sd(const sdsl::bit_vector &vector, uint64_t num_set_bits) {
    assert(num_set_bits == sdsl::util::cnt_one_bits(vector));

    // check if it needs to be inverted
    if (num_set_bits <= vector.size() / 2) {
        // vector is sparse, no need to invert
        vector_ = sdsl::sd_vector<>(vector);
        inverted_ = false;
    } else {
        // invert
        sdsl::sd_vector_builder builder(vector.size(), vector.size() - num_set_bits);
        ::call_zeros(vector, [&builder](uint64_t i) { builder.set(i); });
        vector_ = sdsl::sd_vector<>(builder);
        inverted_ = true;
    }
    slct0_ = decltype(slct0_)(&vector_);
    slct1_ = decltype(slct1_)(&vector_);
    rk1_ = decltype(rk1_)(&vector_);
}

bit_vector_sd::bit_vector_sd(const bit_vector_sd &other) {
    *this = other;
}

bit_vector_sd::bit_vector_sd(bit_vector_sd&& other) noexcept {
    *this = std::move(other);
}

bit_vector_sd
::bit_vector_sd(const std::function<void(const VoidCall<uint64_t>&)> &call_ones,
                uint64_t size,
                uint64_t num_set_bits)
      : inverted_(num_set_bits > size / 2) {
    sdsl::sd_vector_builder builder(size,
                                    !inverted_
                                    ? num_set_bits
                                    : size - num_set_bits);
    if (inverted_) {
        uint64_t last_pos = 0;
        call_ones([&](uint64_t pos) {
            while (last_pos < pos) {
                builder.set(last_pos++);
            }
            ++last_pos;
        });
        while (last_pos < size) {
            builder.set(last_pos++);
        }
    } else {
        call_ones([&](uint64_t pos) { builder.set(pos); });
    }

    assert(builder.items() == builder.capacity());

    vector_ = sdsl::sd_vector<>(builder);
    slct0_ = decltype(slct0_)(&vector_);
    slct1_ = decltype(slct1_)(&vector_);
    rk1_ = decltype(rk1_)(&vector_);
}

bit_vector_sd::bit_vector_sd(std::initializer_list<bool> init)
      : bit_vector_sd(sdsl::bit_vector(init)) {}

bit_vector_sd& bit_vector_sd::operator=(const bit_vector_sd &other) {
    inverted_ = other.inverted_;
    vector_ = other.vector_;
    rk1_ = other.rk1_;
    rk1_.set_vector(&vector_);
    slct1_ = other.slct1_;
    slct1_.set_vector(&vector_);
    slct0_ = other.slct0_;
    slct0_.set_vector(&vector_);
    return *this;
}

bit_vector_sd& bit_vector_sd::operator=(bit_vector_sd&& other) noexcept {
    inverted_ = other.inverted_;
    vector_ = std::move(other.vector_);
    rk1_ = std::move(other.rk1_);
    rk1_.set_vector(&vector_);
    slct1_ = std::move(other.slct1_);
    slct1_.set_vector(&vector_);
    slct0_ = std::move(other.slct0_);
    slct0_.set_vector(&vector_);
    return *this;
}

std::unique_ptr<bit_vector> bit_vector_sd::copy() const {
    return std::make_unique<bit_vector_sd>(*this);
}

uint64_t bit_vector_sd::rank1(uint64_t id) const {
    //the rank method in SDSL does not include id in the count
    size_t idx = id >= this->size() ? this->size() : id + 1;
    return !inverted_ ? rk1_(idx)
                      : idx - rk1_(idx);
}

uint64_t bit_vector_sd::select0(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= rank0(size() - 1));
    assert(num_set_bits() == rank1(size() - 1));

    return !inverted_ ? slct0_(id) : slct1_(id);
}

uint64_t bit_vector_sd::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));

    return !inverted_ ? slct1_(id) : slct0_(id);
}

uint64_t bit_vector_sd::next1(uint64_t pos) const {
    assert(pos < size());

    return ::next1(*this, pos, MAX_ITER_BIT_VECTOR_SD);
}

uint64_t bit_vector_sd::prev1(uint64_t pos) const {
    assert(pos < size());

    return ::prev1(*this, pos, MAX_ITER_BIT_VECTOR_SD);
}

void bit_vector_sd::set(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

bool bit_vector_sd::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id] != inverted_;
}

uint64_t bit_vector_sd::get_int(uint64_t id, uint32_t width) const {
    if (inverted_)
        return ~vector_.get_int(id, width) & sdsl::bits::lo_set[width];

    return vector_.get_int(id, width);
}

void bit_vector_sd::insert_bit(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

void bit_vector_sd::delete_bit(uint64_t) {
    throw std::runtime_error("Not supported");
}

bool bit_vector_sd::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);
        inverted_ = in.get();
        if (!in.good())
            return false;
        rk1_ = decltype(rk1_)(&vector_);
        slct1_ = decltype(slct1_)(&vector_);
        slct0_ = decltype(slct0_)(&vector_);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_sd." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void bit_vector_sd::serialize(std::ostream &out) const {
    vector_.serialize(out);
    if (!out.put(inverted_).good())
        throw std::ofstream::failure("Error when dumping bit_vector_sd");
}

sdsl::bit_vector bit_vector_sd::to_vector() const {
    sdsl::bit_vector vector(size(), inverted_);
    uint64_t max_rank = rk1_(size());
    for (uint64_t i = 1; i <= max_rank; ++i) {
        vector[slct1_(i)] = !inverted_;
    }
    return vector;
}

void bit_vector_sd::call_ones_in_range(uint64_t begin, uint64_t end,
                                       const VoidCall<uint64_t> &callback) const {
    assert(begin <= end);
    assert(end <= size());

    if (!inverted_) {
        // sparse
        uint64_t num_ones = rk1_(end);
        for (uint64_t i = rk1_(begin) + 1; i <= num_ones; ++i) {
            callback(slct1_(i));
        }
    } else {
        // dense, vector_ keeps positions of zeros
        uint64_t one_pos = 0;
        uint64_t zero_pos = 0;
        uint64_t num_zeros = rk1_(end);
        for (uint64_t r = rk1_(begin) + 1; r <= num_zeros; ++r) {
            zero_pos = slct1_(r);
            while (one_pos < zero_pos) {
                callback(one_pos++);
            }
            one_pos++;
        }
        while (one_pos < end) {
            callback(one_pos++);
        }
    }
}


////////////////////////////////////////////////////////////////
//  bit_vector_rrr, sdsl compressed with rank-select support  //
////////////////////////////////////////////////////////////////

template <size_t log_block_size>
bit_vector_rrr<log_block_size>
::bit_vector_rrr(uint64_t size, bool value)
      : vector_(sdsl::bit_vector(size, value)),
        rk1_(&vector_),
        slct1_(&vector_),
        slct0_(&vector_) {}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>
::bit_vector_rrr(const sdsl::bit_vector &vector)
      : vector_(vector),
        rk1_(&vector_),
        slct1_(&vector_),
        slct0_(&vector_) {}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>
::bit_vector_rrr(const bit_vector_rrr<log_block_size> &other) {
    *this = other;
}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>
::bit_vector_rrr(bit_vector_rrr<log_block_size>&& other) noexcept {
    *this = std::move(other);
}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>
::bit_vector_rrr(std::initializer_list<bool> init)
      : bit_vector_rrr(sdsl::bit_vector(init)) {}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>&
bit_vector_rrr<log_block_size>::operator=(const bit_vector_rrr<log_block_size> &other) {
    vector_ = other.vector_;
    rk1_ = other.rk1_;
    rk1_.set_vector(&vector_);
    slct1_ = other.slct1_;
    slct1_.set_vector(&vector_);
    slct0_ = other.slct0_;
    slct0_.set_vector(&vector_);
    return *this;
}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>&
bit_vector_rrr<log_block_size>::operator=(bit_vector_rrr<log_block_size>&& other) noexcept {
    vector_ = std::move(other.vector_);
    rk1_ = std::move(other.rk1_);
    rk1_.set_vector(&vector_);
    slct1_ = std::move(other.slct1_);
    slct1_.set_vector(&vector_);
    slct0_ = std::move(other.slct0_);
    slct0_.set_vector(&vector_);
    return *this;
}

template <size_t log_block_size>
std::unique_ptr<bit_vector> bit_vector_rrr<log_block_size>::copy() const {
    return std::make_unique<bit_vector_rrr<log_block_size>>(*this);
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>::rank1(uint64_t id) const {
    //the rank method in SDSL does not include id in the count
    return rk1_(id >= this->size() ? this->size() : id + 1);
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>::select0(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= size() - num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));

    return slct0_(id);
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));

    return slct1_(id);
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>::next1(uint64_t pos) const {
    assert(pos < size());

    return ::next1(*this, pos, MAX_ITER_BIT_VECTOR_RRR);
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>::prev1(uint64_t pos) const {
    assert(pos < size());

    return ::prev1(*this, pos, MAX_ITER_BIT_VECTOR_RRR);
}

template <size_t log_block_size>
void bit_vector_rrr<log_block_size>::set(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

template <size_t log_block_size>
bool bit_vector_rrr<log_block_size>::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>
::get_int(uint64_t id, uint32_t width) const {
    return vector_.get_int(id, width);
}

template <size_t log_block_size>
void bit_vector_rrr<log_block_size>::insert_bit(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

template <size_t log_block_size>
void bit_vector_rrr<log_block_size>::delete_bit(uint64_t) {
    throw std::runtime_error("Not supported");
}

template <size_t log_block_size>
bool bit_vector_rrr<log_block_size>::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);
        if (!in.good())
            return false;
        rk1_ = decltype(rk1_)(&vector_);
        slct1_ = decltype(slct1_)(&vector_);
        slct0_ = decltype(slct0_)(&vector_);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_rrr." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

template <size_t log_block_size>
void bit_vector_rrr<log_block_size>::serialize(std::ostream &out) const {
    vector_.serialize(out);

    if (!out.good())
        throw std::ofstream::failure("Error when dumping bit_vector_rrr");
}

template <size_t log_block_size>
sdsl::bit_vector bit_vector_rrr<log_block_size>::to_vector() const {
    if (2 * num_set_bits() < size()) {
        // sparse
        sdsl::bit_vector vector(size(), 0);
        uint64_t max_rank = size() ? rank1(size() - 1) : 0;
        for (uint64_t i = 1; i <= max_rank; ++i) {
            vector[slct1_(i)] = 1;
        }
        return vector;
    } else {
        // dense
        sdsl::bit_vector vector(size(), 1);
        uint64_t max_rank = size() ? rank0(size() - 1) : 0;
        for (uint64_t i = 1; i <= max_rank; ++i) {
            vector[slct0_(i)] = 0;
        }
        return vector;
    }
}

template <size_t log_block_size>
void bit_vector_rrr<log_block_size>
::call_ones_in_range(uint64_t begin, uint64_t end,
                     const VoidCall<uint64_t> &callback) const {
    assert(begin <= end);
    assert(end <= size());

    if (2 * num_set_bits() < size()) {
        // sparse
        uint64_t num_ones = end ? rank1(end - 1) : 0;
        for (uint64_t r = begin ? rank1(begin - 1) + 1 : 1; r <= num_ones; ++r) {
            callback(select1(r));
        }
    } else {
        // dense
        uint64_t one_pos = 0;
        uint64_t zero_pos = 0;
        uint64_t num_zeros = end ? rank0(end - 1) : 0;
        for (uint64_t r = begin ? rank0(begin - 1) + 1 : 1; r <= num_zeros; ++r) {
            zero_pos = select0(r);
            while (one_pos < zero_pos) {
                callback(one_pos++);
            }
            one_pos++;
        }
        while (one_pos < end) {
            callback(one_pos++);
        }
    }
}

template class bit_vector_rrr<3>;
template class bit_vector_rrr<8>;
template class bit_vector_rrr<15>;
template class bit_vector_rrr<31>;
template class bit_vector_rrr<63>;
template class bit_vector_rrr<127>;
template class bit_vector_rrr<255>;


////////////////////////////////////////////
// bit_vector_adaptive abstract interface //
////////////////////////////////////////////

bit_vector_adaptive::bit_vector_adaptive(const bit_vector_adaptive &other)
      : vector_(other.vector_->copy()) {}

bit_vector_adaptive&
bit_vector_adaptive::operator=(const bit_vector_adaptive &other) {
    vector_ = other.vector_->copy();
    return *this;
}

bit_vector_adaptive::VectorCode
bit_vector_adaptive::representation_tag() const {
    if (dynamic_cast<const bit_vector_sd*>(vector_.get())) {
        return VectorCode::SD_VECTOR;

    } else if (dynamic_cast<const bit_vector_rrr<>*>(vector_.get())) {
        return VectorCode::RRR_VECTOR;

    } else if (dynamic_cast<const bit_vector_stat*>(vector_.get())) {
        return VectorCode::STAT_VECTOR;

    } else {
        throw std::runtime_error("Unsupported type");
    }
}

bool bit_vector_adaptive::load(std::istream &in) {
    switch (load_number(in)) {
        case VectorCode::SD_VECTOR:
            vector_.reset(new bit_vector_sd());
            break;
        case VectorCode::RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>());
            break;
        case VectorCode::STAT_VECTOR:
            vector_.reset(new bit_vector_stat());
            break;
        default:
            return false;
    }
    return vector_->load(in);
}

void bit_vector_adaptive::serialize(std::ostream &out) const {
    serialize_number(out, representation_tag());
    vector_->serialize(out);
}


////////////////////////////////////////////
// bit_vector_adaptive_stat               //
////////////////////////////////////////////

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(uint64_t size, bool value) {
    switch (optimal_representation(size, size * static_cast<uint64_t>(value))) {
        case SD_VECTOR:
            vector_.reset(new bit_vector_sd(size, value));
            break;
        case RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>(size, value));
            break;
        case STAT_VECTOR:
            vector_.reset(new bit_vector_stat(size, value));
            break;
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(const bit_vector &vector) {
    switch (optimal_representation(vector.size(), vector.num_set_bits())) {
        case SD_VECTOR:
            vector_.reset(new bit_vector_sd(vector.copy_to<bit_vector_sd>()));
            break;
        case RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>(vector.copy_to<bit_vector_rrr<>>()));
            break;
        case STAT_VECTOR:
            vector_.reset(new bit_vector_stat(vector.copy_to<bit_vector_stat>()));
            break;
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(const sdsl::bit_vector &vector) {
    uint64_t num_set_bits = sdsl::util::cnt_one_bits(vector);

    switch (optimal_representation(vector.size(), num_set_bits)) {
        case SD_VECTOR:
            vector_.reset(new bit_vector_sd(vector, num_set_bits));
            break;
        case RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>(vector));
            break;
        case STAT_VECTOR:
            vector_.reset(new bit_vector_stat(vector, num_set_bits));
            break;
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(bit_vector&& vector) {
    switch (optimal_representation(vector.size(), vector.num_set_bits())) {
        case SD_VECTOR:
            vector_.reset(new bit_vector_sd(vector.convert_to<bit_vector_sd>()));
            break;
        case RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>(vector.convert_to<bit_vector_rrr<>>()));
            break;
        case STAT_VECTOR:
            vector_.reset(new bit_vector_stat(vector.convert_to<bit_vector_stat>()));
            break;
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
bit_vector_adaptive_stat<optimal_representation>
::bit_vector_adaptive_stat(const VoidCall<const VoidCall<uint64_t>&> &call_ones,
                           uint64_t size,
                           uint64_t num_set_bits) {
    switch (optimal_representation(size, num_set_bits)) {
        case SD_VECTOR:
            vector_.reset(new bit_vector_sd(call_ones, size, num_set_bits));
            break;
        case RRR_VECTOR:
            vector_.reset(new bit_vector_rrr<>(
                bit_vector_stat(call_ones, size).convert_to<bit_vector_rrr<>>()
            ));
            break;
        case STAT_VECTOR:
            vector_.reset(new bit_vector_stat(call_ones, size));
            break;
    }
}

template <bit_vector_adaptive::DefineRepresentation optimal_representation>
std::unique_ptr<bit_vector>
bit_vector_adaptive_stat<optimal_representation>
::copy() const {
    auto copy = new bit_vector_adaptive_stat();
    copy->vector_ = vector_->copy();
    return std::unique_ptr<bit_vector> { copy };
}

template class bit_vector_adaptive_stat<smallest_representation>;
template class bit_vector_adaptive_stat<smart_representation>;


///////////////////////////////////////////////////
// bit_vector_small, the smallest representation //
///////////////////////////////////////////////////

bit_vector_adaptive::VectorCode
smallest_representation(uint64_t size, uint64_t num_set_bits) {
    assert(num_set_bits <= size);

    return predict_size<bit_vector_sd>(size, num_set_bits)
            < predict_size<bit_vector_rrr<>>(size, num_set_bits)
                  ? bit_vector_adaptive::VectorCode::SD_VECTOR
                  : bit_vector_adaptive::VectorCode::RRR_VECTOR;
}

////////////////////////////////////////////////
// bit_vector_smart, fast/small hybrid method //
////////////////////////////////////////////////

bit_vector_adaptive::VectorCode
smart_representation(uint64_t size, uint64_t num_set_bits) {
    assert(num_set_bits <= size);

    return predict_size<bit_vector_sd>(size, num_set_bits)
            < predict_size<bit_vector_stat>(size, num_set_bits)
                  ? bit_vector_adaptive::VectorCode::SD_VECTOR
                  : bit_vector_adaptive::VectorCode::STAT_VECTOR;
}


double logbinomial(uint64_t n, uint64_t m) {
    return (lgamma(n + 1)
                - lgamma(m + 1)
                - lgamma(n - m + 1)) / log(2);
}

double entropy(double q) {
    assert(q >= 0);
    assert(q <= 1);

    if (q == 0 || q == 1)
        return 0;

    return q * log2(q) + (1 - q) * log2(1 - q);
}

// TODO: write unit tests for these. Check if approximately equals to the serialized dumps
uint64_t bv_space_taken_rrr(uint64_t size, uint64_t num_set_bits, uint8_t block_size) {
    return std::ceil(logbinomial(size, num_set_bits))
            + (size + block_size - 1) / block_size * std::ceil(log2(block_size + 1));
}

uint64_t bv_space_taken_sd(uint64_t size, uint64_t num_set_bits) {
    num_set_bits = std::min(num_set_bits, size - num_set_bits);
    return std::ceil((log2(size + 1) - log2(num_set_bits + 1) + 3) * num_set_bits);
}

uint64_t bv_space_taken_stat(uint64_t size, uint64_t num_set_bits) {
    return STAT_BITS_PER_BIT_IF_SPARSE * size
            + (STAT_BITS_PER_BIT_IF_DENSE - STAT_BITS_PER_BIT_IF_SPARSE) * num_set_bits;
}

template <class VectorType>
uint64_t predict_size(uint64_t size, uint64_t num_set_bits) {
    assert(size >= num_set_bits);

    if constexpr(std::is_base_of<bit_vector_stat, VectorType>::value)
        return bv_space_taken_stat(size, num_set_bits);

    if constexpr(std::is_base_of<bit_vector_sd, VectorType>::value)
        return bv_space_taken_sd(size, num_set_bits);

    if constexpr(std::is_base_of<bit_vector_rrr<15>, VectorType>::value)
        return bv_space_taken_rrr(size, num_set_bits, 15);

    if constexpr(std::is_base_of<bit_vector_rrr<31>, VectorType>::value)
        return bv_space_taken_rrr(size, num_set_bits, 31);

    if constexpr(std::is_base_of<bit_vector_rrr<63>, VectorType>::value)
        return bv_space_taken_rrr(size, num_set_bits, 63);

    if constexpr(std::is_base_of<bit_vector_rrr<127>, VectorType>::value)
        return bv_space_taken_rrr(size, num_set_bits, 127);

    if constexpr(std::is_base_of<bit_vector_rrr<255>, VectorType>::value)
        return bv_space_taken_rrr(size, num_set_bits, 255);

    if constexpr(std::is_base_of<bit_vector_small, VectorType>::value)
        return std::min(bv_space_taken_sd(size, num_set_bits),
                        predict_size<bit_vector_rrr<>>(size, num_set_bits));

    if constexpr(std::is_base_of<bit_vector_smart, VectorType>::value)
        return std::min(bv_space_taken_sd(size, num_set_bits),
                        predict_size<bit_vector_stat>(size, num_set_bits));

    throw std::runtime_error(std::string("Error: unknown space taken for this bit_vector")
                                 + typeid(VectorType).name());
}

// template uint64_t predict_size<bit_vector_dyn>(uint64_t, uint64_t);
template uint64_t predict_size<bit_vector_stat>(uint64_t, uint64_t);
template uint64_t predict_size<bit_vector_sd>(uint64_t, uint64_t);
// template uint64_t predict_size<bit_vector_rrr<3>>(uint64_t, uint64_t);
// template uint64_t predict_size<bit_vector_rrr<8>>(uint64_t, uint64_t);
template uint64_t predict_size<bit_vector_rrr<15>>(uint64_t, uint64_t);
template uint64_t predict_size<bit_vector_rrr<31>>(uint64_t, uint64_t);
template uint64_t predict_size<bit_vector_rrr<63>>(uint64_t, uint64_t);
template uint64_t predict_size<bit_vector_rrr<127>>(uint64_t, uint64_t);
template uint64_t predict_size<bit_vector_rrr<255>>(uint64_t, uint64_t);
template uint64_t predict_size<bit_vector_small>(uint64_t, uint64_t);
template uint64_t predict_size<bit_vector_smart>(uint64_t, uint64_t);


// indexes are distinct and sorted
sdsl::bit_vector subvector(const bit_vector &col,
                           const std::vector<uint64_t> &indexes) {
    assert(indexes.size() <= col.size());

    sdsl::bit_vector shrinked(indexes.size(), 0);

    uint64_t max_rank = col.num_set_bits();
    if (!max_rank)
        return shrinked;

    // the case of uncompressed vector
    if (dynamic_cast<const bit_vector_stat *>(&col)) {
        for (size_t j = 0; j < indexes.size(); ++j) {
            if (col[indexes[j]])
                shrinked[j] = true;
        }
        return shrinked;
    }

    uint64_t cur_rank = 1;
    uint64_t next_pos = col.select1(1);

    for (size_t j = 0; j < indexes.size(); ++j) {
        if (indexes[j] < next_pos)
            continue;

        if (indexes[j] == next_pos) {
            shrinked[j] = true;
            continue;
        }

        // indexes[j] > next_pos
        if (col[indexes[j]]) {
            shrinked[j] = true;
            continue;
        }

        // we found a zero, update next_pos
        cur_rank = col.rank1(indexes[j]) + 1;
        if (cur_rank > max_rank)
            return shrinked;

        next_pos = col.select1(cur_rank);
    }

    return shrinked;
}
