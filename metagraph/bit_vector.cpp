#include "bit_vector.hpp"

#include <cassert>
#include <libmaus2/bitio/putBit.hpp>


bit_vector::BoolVector bit_vector::to_vector() const {
    BoolVector result(size());
    for (uint64_t i = 0; i < size(); ++i) {
        result[i] = operator[](i);
    }
    return result;
}

std::ostream& operator<<(std::ostream &os, const bit_vector &bv) {
    for (bool a : bv.to_vector()) {
        os << a;
    }
    return os;
}

template <class Vector>
Vector bit_vector::convert_to() {
    if (dynamic_cast<Vector*>(this)) {
        return dynamic_cast<Vector&&>(*this);
    } else if (dynamic_cast<bit_vector_small*>(this)) {
        sdsl::bit_vector bv(size());
        for (uint64_t i = 1; i <= get_num_set_bits(); ++i) {
            bv[select1(i)] = 1;
        }
        return Vector(std::move(bv));
    } else {
        sdsl::bit_vector bv(size());
        for (uint64_t i = 0; i < size(); ++i) {
            if (operator[](i))
                bv[i] = 1;
        }
        return Vector(std::move(bv));
    }
}
template bit_vector_dyn   bit_vector::convert_to<bit_vector_dyn>();
template bit_vector_stat  bit_vector::convert_to<bit_vector_stat>();
template bit_vector_small bit_vector::convert_to<bit_vector_small>();

/////////////////////////////
// bit_vector_dyn, libmaus //
/////////////////////////////

template <class BitVector>
std::vector<uint64_t> pack_bits(const BitVector &v) {
    std::vector<uint64_t> bits((v.size() + 63) / 64);
    for (size_t i = 0; i < v.size(); ++i) {
        libmaus2::bitio::putBit(bits.data(), i, v[i]);
    }
    return bits;
}

bit_vector_dyn::bit_vector_dyn(const std::vector<uint64_t> &v, size_t num_bits)
      : vector_(num_bits, v.data()) {}

template <class BinaryVector>
bit_vector_dyn::bit_vector_dyn(const BinaryVector &v)
      : bit_vector_dyn(pack_bits(v), v.size()) {}
template bit_vector_dyn::bit_vector_dyn(const BoolVector&);
template bit_vector_dyn::bit_vector_dyn(const sdsl::bit_vector&);

bit_vector_dyn::bit_vector_dyn(const bit_vector_dyn &v)
      : vector_(v.vector_) {}

bit_vector_dyn::bit_vector_dyn(std::initializer_list<bool> init)
      : bit_vector_dyn(pack_bits(std::vector<bool>(init)), init.size()) {}

bit_vector_dyn::bit_vector_dyn(bit_vector_dyn&& v)
      : vector_(std::move(v.vector_)) {}

uint64_t bit_vector_dyn::rank1(uint64_t id) const {
    return vector_.rank1(id);
}

uint64_t bit_vector_dyn::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= rank1(size() - 1));
    return vector_.select1(id - 1);
}

void bit_vector_dyn::set(uint64_t id, bool val) {
    vector_.set(id, val);
}

void bit_vector_dyn::setBitQuick(uint64_t id, bool val) {
    vector_.setBitQuick(id, val);
}

bool bit_vector_dyn::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
}

void bit_vector_dyn::insertBit(uint64_t id, bool val) {
    vector_.insertBit(id, val);
}

void bit_vector_dyn::deleteBit(uint64_t id) {
    assert(size() > 0);
    vector_.deleteBit(id);
}

bool bit_vector_dyn::deserialise(std::istream &in) {
    if (!in.good())
        return false;

    //TODO: catch reading errors
    vector_.deserialise(in);
    return true;
}

void bit_vector_dyn::serialise(std::ostream &out) const {
    vector_.serialise(out);
}

///////////////////////////////////////////////
// bit_vector_stat, sdsl rank-select support //
///////////////////////////////////////////////

bit_vector_stat::bit_vector_stat(uint64_t size, bool value)
      : vector_(size, value) {
    if (value)
        num_set_bits_ = size;
}

bit_vector_stat::bit_vector_stat(const bit_vector_stat &other)
      : vector_(other.vector_), num_set_bits_(other.num_set_bits_) {}

bit_vector_stat::bit_vector_stat(const std::vector<bool> &other)
      : vector_(other.size(), 0) {
    for (uint64_t i = 0; i < other.size(); ++i) {
        if (other.at(i)) {
            vector_[i] = 1;
            num_set_bits_++;
        }
    }
}

bit_vector_stat::bit_vector_stat(std::initializer_list<bool> init)
      : vector_(init) {
    num_set_bits_ = std::count(vector_.begin(), vector_.end(), 1);
}

bit_vector_stat::bit_vector_stat(bit_vector_stat&& other)
      : vector_(std::move(other.vector_)), num_set_bits_(other.num_set_bits_) {}

bit_vector_stat::bit_vector_stat(sdsl::bit_vector&& vector)
      : vector_(std::move(vector)) {
    num_set_bits_ = std::count(vector_.begin(), vector_.end(), 1);
}

bit_vector_stat& bit_vector_stat::operator=(sdsl::bit_vector&& vector) {
    vector_ = std::move(vector);
    num_set_bits_ = std::count(vector_.begin(), vector_.end(), 1);
    return *this;
}

uint64_t bit_vector_stat::rank1(uint64_t id) const {
    if (requires_update_)
        const_cast<bit_vector_stat*>(this)->init_rs();
    //the rank method in SDSL does not include id in the count
    return rk_(id >= this->size() ? this->size() : id + 1);
}

uint64_t bit_vector_stat::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits_);
    assert(num_set_bits_ == rank1(size() - 1));

    if (requires_update_)
        const_cast<bit_vector_stat*>(this)->init_rs();
    return slct_(id);
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

bool bit_vector_stat::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
}

void bit_vector_stat::insertBit(uint64_t id, bool val) {
    assert(id <= size());

    if (val)
        num_set_bits_++;

    vector_.resize(size() + 1);
    std::copy_backward(vector_.begin() + id, vector_.end() - 1, vector_.end());

    vector_[id] = val;
    requires_update_ = true;
}

void bit_vector_stat::deleteBit(uint64_t id) {
    assert(size() > 0);
    assert(id < size());

    if (vector_[id])
        num_set_bits_--;

    std::copy(vector_.begin() + id + 1, vector_.end(), vector_.begin() + id);

    vector_.resize(vector_.size() - 1);
    requires_update_ = true;
}

bool bit_vector_stat::deserialise(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);
        /*
        rk_.load(in);
        slct_.load(in);

        rk_.set_vector(&vector_);
        slct_.set_vector(&vector_);
        */
        num_set_bits_ = std::count(vector_.begin(), vector_.end(), 1);
        requires_update_ = true;
        init_rs();
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_stat." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void bit_vector_stat::serialise(std::ostream &out) const {
    vector_.serialize(out);
    /*
    if (requires_update_)
        init_rs();

    rk_.serialize(out);
    slct_.serialize(out);
    */
}

void bit_vector_stat::init_rs() {
    rk_ = sdsl::rank_support_v5<>(&vector_);
    slct_ = sdsl::select_support_mcl<>(&vector_);
    requires_update_ = false;
}

////////////////////////////////////////////////////////////////
// bit_vector_small, sdsl compressed with rank-select support //
////////////////////////////////////////////////////////////////


sdsl::bit_vector bit_vector_small::invert(const sdsl::bit_vector &other, uint64_t num_set_bits) {
    sdsl::bit_vector bv(other.size());
    sdsl::bit_vector::select_0_type slct_0(&other);
    uint64_t num_unset_bits = other.size() - num_set_bits;
    for (uint64_t i = 1; i <= num_unset_bits; ++i) {
        bv[slct_0(i)] = 1;
    }
    return bv;
}

bit_vector_small::bit_vector_small(uint64_t size, bool value)
      : vector_(sdsl::bit_vector(size, !value)),
        num_set_bits_(size * value),
        rk_(&vector_),
        slct_(&vector_) {}

bit_vector_small::bit_vector_small(const bit_vector_small &other)
      : vector_(other.vector_), num_set_bits_(other.num_set_bits_), rk_(&vector_), slct_(&vector_) {}

bit_vector_small::bit_vector_small(const std::vector<bool> &other)
      : num_set_bits_(0) {
    sdsl::bit_vector vector(other.size(), 0);
    for (uint64_t i = 0; i < other.size(); ++i) {
        if (!other[i]) {
            vector[i] = 1;
        } else {
            num_set_bits_++;
        }
    }
    vector_ = decltype(vector_)(vector);
    rk_ = decltype(rk_)(&vector_);
    slct_ = decltype(slct_)(&vector_);
}

bit_vector_small::bit_vector_small(std::initializer_list<bool> init)
        : num_set_bits_(0) {
    sdsl::bit_vector vector(init.size(), 0);
    uint64_t i = 0;
    for (bool a : init) {
        if (!a) {
            vector[i] = 1;
        } else {
            num_set_bits_++;
        }
        ++i;
    }
    vector_ = decltype(vector_)(vector);
    rk_ = decltype(rk_)(&vector_);
    slct_ = decltype(slct_)(&vector_);
}

bit_vector_small::bit_vector_small(bit_vector_small&& vector)
      : vector_(std::move(vector_)),
        num_set_bits_(vector.num_set_bits_),
        rk_(&vector_),
        slct_(&vector_) {}

bit_vector_small::bit_vector_small(sdsl::bit_vector&& vector) {
    *this = std::move(vector);
}

bit_vector_small& bit_vector_small::operator=(sdsl::bit_vector&& vector) {
    num_set_bits_ = std::count(vector.begin(), vector.end(), 1);
    vector_ = invert(std::move(vector), num_set_bits_);
    rk_ = decltype(rk_)(&vector_);
    slct_ = decltype(slct_)(&vector_);
    return *this;
}

uint64_t bit_vector_small::rank1(uint64_t id) const {
    //the rank method in SDSL does not include id in the count
    return rk_(id >= this->size() ? this->size() : id + 1);
}

uint64_t bit_vector_small::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits_);
    assert(num_set_bits_ == rank1(size() - 1));

    return slct_(id);
}

void bit_vector_small::set(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

bool bit_vector_small::operator[](uint64_t id) const {
    assert(id < size());
    return !vector_[id];
}

void bit_vector_small::insertBit(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

void bit_vector_small::deleteBit(uint64_t) {
    throw std::runtime_error("Not supported");
}

bool bit_vector_small::deserialise(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);
        num_set_bits_ = std::count(vector_.begin(), vector_.end(), 0);
        rk_ = decltype(rk_)(&vector_);
        slct_ = decltype(slct_)(&vector_);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_small." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void bit_vector_small::serialise(std::ostream &out) const {
    vector_.serialize(out);
}
