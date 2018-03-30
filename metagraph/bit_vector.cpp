#include "bit_vector.hpp"

#include <cassert>
#include <libmaus2/bitio/putBit.hpp>


std::vector<bool> bit_vector::to_vector() const {
    std::vector<bool> result(size());
    for (uint64_t i = 0; i < size(); ++i) {
        result[i] = operator[](i);
    }
    return result;
}

std::ostream& operator<<(std::ostream &os, const bit_vector &bv) {
    bv.print(os);
    return os;
}


/////////////////////////////
// bit_vector_dyn, libmaus //
/////////////////////////////

std::vector<uint64_t> pack_bits(const std::vector<bool> &v) {
    std::vector<uint64_t> bits((v.size() + 63) / 64);
    for (size_t i = 0; i < v.size(); ++i) {
        libmaus2::bitio::putBit(bits.data(), i, v[i]);
    }
    return bits;
}

bit_vector_dyn::bit_vector_dyn(const std::vector<uint64_t> &v, size_t num_bits)
      : vector_(num_bits, v.data()) {}

bit_vector_dyn::bit_vector_dyn(const std::vector<bool> &v)
      : bit_vector_dyn(pack_bits(v), v.size()) {}

bit_vector_dyn::bit_vector_dyn(const bit_vector_dyn &v)
      : vector_(v.vector_) {}

bit_vector_dyn::bit_vector_dyn(const bit_vector &v)
      : bit_vector_dyn(v.to_vector()) {}

bit_vector_dyn::bit_vector_dyn(std::istream &in) {
    if (!deserialise(in)) {
        throw "ERROR: deserialisation error";
    }
}

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

bit_vector_stat::bit_vector_stat(const bit_vector &other)
      : vector_(other.size(), 0) {
    for (uint64_t i = 0; i < other.size(); ++i) {
        if (other[i])
            vector_[i] = 1;
    }
    num_set_bits_ = other.get_num_set_bits();
}

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

bit_vector_stat::bit_vector_stat(std::istream &in) {
    if (!deserialise(in)) {
        throw "ERROR: deserialisation error";
    }
}

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
        num_set_bits_ = std::count(vector_.begin(), vector_.end(), 1);
        requires_update_ = true;
        return true;
    } catch (...) {
        return false;
    }
}

void bit_vector_stat::serialise(std::ostream &out) const {
    vector_.serialize(out);
}

void bit_vector_stat::init_rs() {
    rk_ = sdsl::rank_support_v5<>(&vector_);
    slct_ = sdsl::select_support_mcl<>(&vector_);
    requires_update_ = false;
}
