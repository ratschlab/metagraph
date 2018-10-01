#include "bit_vector.hpp"

#include <cassert>
#include <libmaus2/bitio/putBit.hpp>

#include "serialization.hpp"


std::vector<bool> bit_vector::to_vector() const {
    std::vector<bool> result(size());
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
        return Vector(dynamic_cast<Vector&&>(*this));
    } else if (dynamic_cast<bit_vector_small*>(this)) {
        sdsl::bit_vector bv(size(), 0);
        uint64_t max_rank = num_set_bits();
        for (uint64_t i = 1; i <= max_rank; ++i) {
            bv[select1(i)] = 1;
        }
        return Vector(std::move(bv));
    } else if (dynamic_cast<bit_vector_stat*>(this)) {
        return Vector(std::move(dynamic_cast<bit_vector_stat*>(this)->vector_));
    } else {
        sdsl::bit_vector bv(size(), 0);
        for (uint64_t i = 0; i < size(); ++i) {
            if (operator[](i))
                bv[i] = 1;
        }
        return Vector(std::move(bv));
    }
}
template bit_vector_dyn bit_vector::convert_to<bit_vector_dyn>();
template bit_vector_stat bit_vector::convert_to<bit_vector_stat>();
template bit_vector_small bit_vector::convert_to<bit_vector_small>();
template sdsl::bit_vector bit_vector::convert_to<sdsl::bit_vector>();

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

bit_vector_dyn::bit_vector_dyn(const std::vector<uint64_t> &bits_packed,
                               size_t num_bits)
      : vector_(num_bits, bits_packed.data()) {}

bit_vector_dyn::bit_vector_dyn(uint64_t size, bool value)
      : vector_(size, value) {}

template <class BitVector>
bit_vector_dyn::bit_vector_dyn(const BitVector &v)
      : bit_vector_dyn(pack_bits(v), v.size()) {}

template bit_vector_dyn::bit_vector_dyn(const std::vector<bool> &);
template bit_vector_dyn::bit_vector_dyn(const sdsl::bit_vector &);

bit_vector_dyn::bit_vector_dyn(std::initializer_list<bool> init)
      : bit_vector_dyn(pack_bits(std::vector<bool>(init)), init.size()) {}

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

bool bit_vector_dyn::load(std::istream &in) {
    if (!in.good())
        return false;

    //TODO: catch reading errors
    vector_.deserialise(in);
    return true;
}

void bit_vector_dyn::serialize(std::ostream &out) const {
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

bit_vector_stat::bit_vector_stat(const std::vector<bool> &v)
      : vector_(v.size(), 0) {
    for (uint64_t i = 0; i < v.size(); ++i) {
        if (v[i]) {
            vector_[i] = 1;
            num_set_bits_++;
        }
    }
}

bit_vector_stat::bit_vector_stat(const bit_vector_stat &other) {
    *this = other;
}

bit_vector_stat::bit_vector_stat(sdsl::bit_vector&& vector)
      : vector_(std::move(vector)) {
    num_set_bits_ = std::count(vector_.begin(), vector_.end(), 1);
}

bit_vector_stat::bit_vector_stat(bit_vector_stat&& other) {
    *this = std::move(other);
}

bit_vector_stat::bit_vector_stat(std::initializer_list<bool> init)
      : vector_(init) {
    num_set_bits_ = std::count(vector_.begin(), vector_.end(), 1);
}

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

bit_vector_stat& bit_vector_stat::operator=(bit_vector_stat&& other) {
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

bool bit_vector_stat::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
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

void bit_vector_stat::serialize(std::ostream &out) const {
    vector_.serialize(out);

    if (requires_update_)
        const_cast<bit_vector_stat*>(this)->init_rs();

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

            num_set_bits_ = std::count(vector_.begin(), vector_.end(), 1);
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

void bit_vector_stat::init_rs() {
    std::unique_lock<std::mutex> lock(mu_);

    if (!requires_update_)
        return;

    rk_ = sdsl::rank_support_v5<>(&vector_);
    slct_ = sdsl::select_support_mcl<>(&vector_);

    requires_update_ = false;
}

////////////////////////////////////////////////////////////////
// bit_vector_small, sdsl compressed with rank-select support //
////////////////////////////////////////////////////////////////

template <class BitVector>
bool invert_vector_if_dense(const BitVector &vector, sdsl::bit_vector *result) {
    assert(result);

    result->resize(vector.size());

    uint64_t num_set_bits = std::count(vector.begin(), vector.end(), true);
    if (num_set_bits <= vector.size() / 2) {
        // vector is sparse, no need to invert
        std::copy(vector.begin(), vector.end(), result->begin());
        return false;
    } else {
        // invert
        for (uint64_t i = 0; i < vector.size(); ++i) {
            (*result)[i] = !vector[i];
        }
        return true;
    }
}

bit_vector_small::bit_vector_small(uint64_t size, bool value)
      : inverted_(value && size),
        vector_(sdsl::bit_vector(size, 0)),
        rk1_(&vector_) {
    if (!inverted_) {
        slct1_ = decltype(slct1_)(&vector_);
    } else {
        slct0_ = decltype(slct0_)(&vector_);
    }
}

template <class BitVector>
bit_vector_small::bit_vector_small(const BitVector &vector) {
    sdsl::bit_vector bv;
    inverted_ = invert_vector_if_dense(vector, &bv);
    vector_ = decltype(vector_)(bv);
    rk1_ = decltype(rk1_)(&vector_);
    if (!inverted_) {
        slct1_ = decltype(slct1_)(&vector_);
    } else {
        slct0_ = decltype(slct0_)(&vector_);
    }
}

template bit_vector_small::bit_vector_small(const std::vector<bool> &);
template bit_vector_small::bit_vector_small(const sdsl::bit_vector &);

bit_vector_small::bit_vector_small(const bit_vector_small &other) {
    *this = other;
}

bit_vector_small::bit_vector_small(bit_vector_small&& other) {
    *this = std::move(other);
}

bit_vector_small::bit_vector_small(std::initializer_list<bool> init)
      : bit_vector_small(sdsl::bit_vector(init)) {}

bit_vector_small& bit_vector_small::operator=(const bit_vector_small &other) {
    inverted_ = other.inverted_;
    vector_ = other.vector_;
    rk1_ = other.rk1_;
    rk1_.set_vector(&vector_);
    if (!inverted_) {
        slct1_ = other.slct1_;
        slct1_.set_vector(&vector_);
    } else {
        slct0_ = other.slct0_;
        slct0_.set_vector(&vector_);
    }
    return *this;
}

bit_vector_small& bit_vector_small::operator=(bit_vector_small&& other) {
    inverted_ = other.inverted_;
    vector_ = std::move(other.vector_);
    rk1_ = std::move(other.rk1_);
    rk1_.set_vector(&vector_);
    if (!inverted_) {
        slct1_ = std::move(other.slct1_);
        slct1_.set_vector(&vector_);
    } else {
        slct0_ = std::move(other.slct0_);
        slct0_.set_vector(&vector_);
    }
    return *this;
}

uint64_t bit_vector_small::rank1(uint64_t id) const {
    //the rank method in SDSL does not include id in the count
    size_t idx = id >= this->size() ? this->size() : id + 1;
    return !inverted_ ? rk1_(idx)
                      : idx - rk1_(idx);
}

uint64_t bit_vector_small::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));

    return !inverted_ ? slct1_(id) : slct0_(id);
}

void bit_vector_small::set(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

bool bit_vector_small::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id] != inverted_;
}

void bit_vector_small::insertBit(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

void bit_vector_small::deleteBit(uint64_t) {
    throw std::runtime_error("Not supported");
}

bool bit_vector_small::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);
        inverted_ = in.get();
        if (!in.good())
            return false;
        rk1_ = decltype(rk1_)(&vector_);
        if (!inverted_) {
            slct1_ = decltype(slct1_)(&vector_);
        } else {
            slct0_ = decltype(slct0_)(&vector_);
        }
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_small." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void bit_vector_small::serialize(std::ostream &out) const {
    vector_.serialize(out);
    if (!out.put(inverted_).good())
        throw std::ofstream::failure("Error when dumping bit_vector_small");
}

std::vector<bool> bit_vector_small::to_vector() const {
    std::vector<bool> vector(size(), 0);
    add_to(&vector);
    return vector;
}

void bit_vector_small::add_to(std::vector<bool> *other) const {
    assert(other->size() == size());
    uint64_t max_rank = num_set_bits();
    for (uint64_t i = 1; i <= max_rank; ++i) {
        (*other)[select1(i)] = 1;
    }
}
