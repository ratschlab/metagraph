#include "wavelet_tree.hpp"

#include <cassert>
#include "serialization.hpp"


// TODO: run benchmarks and optimize these parameters
const size_t MAX_ITER_WAVELET_TREE_STAT = 1000;
const size_t MAX_ITER_WAVELET_TREE_DYN = 10;
const size_t MAX_ITER_WAVELET_TREE_SMALL = 10;

typedef libmaus2::bitbtree::BitBTree<6, 64> BitBTree;

/////////////////////////////////
// wavelet_tree shared methods //
/////////////////////////////////

template <class Vector>
sdsl::int_vector<> pack_vector(const Vector &vector, uint8_t bits_per_number) {
    sdsl::int_vector<> packed(vector.size(), 0, bits_per_number);
    for (uint64_t i = 0; i < vector.size(); ++i) {
        packed[i] = vector[i];
    }
    return packed;
}

template <>
sdsl::int_vector<> pack_vector(const sdsl::int_vector<> &vector,
                               uint8_t bits_per_number) {
    if (bits_per_number == vector.width())
        return vector;

    sdsl::int_vector<> packed(vector.size(), 0, bits_per_number);
    for (uint64_t i = 0; i < vector.size(); ++i) {
        packed[i] = vector[i];
    }
    return packed;
}

template <class WaveletTree>
WaveletTree wavelet_tree::convert_to() {
    if (dynamic_cast<WaveletTree*>(this)) {
        return WaveletTree(dynamic_cast<WaveletTree&&>(*this));
    } else if (dynamic_cast<wavelet_tree_small*>(this)
                && typeid(WaveletTree) == typeid(wavelet_tree_stat)) {
        return WaveletTree(
            logsigma(),
            std::move(dynamic_cast<wavelet_tree_small*>(this)->wwt_)
        );
    } else if (dynamic_cast<wavelet_tree_stat*>(this)
                && !dynamic_cast<wavelet_tree_stat*>(this)->requires_update_
                && typeid(WaveletTree) == typeid(wavelet_tree_small)) {
        return WaveletTree(
            logsigma(),
            std::move(dynamic_cast<wavelet_tree_stat*>(this)->wwt_)
        );
    } else {
        uint8_t bits_per_number = logsigma();

        sdsl::int_vector<> vector;
        if (dynamic_cast<wavelet_tree_stat*>(this)) {
            vector = std::move(dynamic_cast<wavelet_tree_stat*>(this)->int_vector_);
            vector.resize(dynamic_cast<wavelet_tree_stat*>(this)->n_);
        } else if (dynamic_cast<wavelet_tree_fast*>(this)) {
            vector = std::move(dynamic_cast<wavelet_tree_fast*>(this)->int_vector_);
        } else {
            vector = to_vector();
        }
        clear();

        return WaveletTree(bits_per_number, std::move(vector));
    }
}

template wavelet_tree_dyn wavelet_tree::convert_to<wavelet_tree_dyn>();
template wavelet_tree_stat wavelet_tree::convert_to<wavelet_tree_stat>();
template wavelet_tree_small wavelet_tree::convert_to<wavelet_tree_small>();
template wavelet_tree_fast wavelet_tree::convert_to<wavelet_tree_fast>();


template <typename Vector>
inline uint64_t next(const Vector &v,
                     uint64_t pos,
                     uint64_t value,
                     size_t num_steps) {
    assert(pos < v.size());

    if (v[pos] == value)
        return pos;

    for (size_t t = 1; t < num_steps; ++t) {
        if (pos + t == v.size() || v[pos + t] == value)
            return pos + t;
    }

    uint64_t rk = v.rank(value, pos) + 1;
    return rk <= v.rank(value, v.size() - 1)
            ? v.select(value, rk)
            : v.size();
}

template <typename Vector>
inline uint64_t prev(const Vector &v,
                     uint64_t pos,
                     uint64_t value,
                     size_t num_steps) {
    assert(pos < v.size());

    for (size_t t = 0; t < num_steps; ++t, --pos) {
        if (v[pos] == value)
            return pos;

        if (pos == 0)
            return v.size();
    }

    uint64_t rk = v.rank(value, pos);
    return rk ? v.select(value, rk)
              : v.size();
}

///////////////////////////////////////////////////////////
// wavelet_tree_stat sdsl rank/select, int_vector access //
///////////////////////////////////////////////////////////

wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma,
                                     uint64_t size, uint64_t value)
      : int_vector_(size, value, logsigma), n_(size) {}

template <class Vector>
wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma, const Vector &vector)
      : int_vector_(pack_vector(vector, logsigma)), n_(int_vector_.size()) {}

template
wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma,
                                     const sdsl::int_vector<> &vector);
template
wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma,
                                     const std::vector<uint8_t> &vector);
template
wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma,
                                     const std::vector<uint64_t> &vector);
template
wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma,
                                     const std::vector<int> &vector);

wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma,
                                     sdsl::int_vector<>&& vector) {
    if (vector.width() == logsigma) {
        int_vector_ = std::move(vector);
    } else {
        int_vector_ = pack_vector(vector, logsigma);
    }
    n_ = int_vector_.size();
}

wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma, sdsl::wt_huff<>&& wwt)
      : int_vector_(pack_vector(wwt, logsigma)),
        wwt_(std::move(wwt)),
        requires_update_(false),
        n_(wwt_.size()) {}

wavelet_tree_stat::wavelet_tree_stat(const wavelet_tree_stat &other) {
    *this = other;
}

wavelet_tree_stat::wavelet_tree_stat(wavelet_tree_stat&& other) noexcept {
    *this = std::move(other);
}

wavelet_tree_stat& wavelet_tree_stat::operator=(const wavelet_tree_stat &other) {
    int_vector_ = other.int_vector_;
    n_ = other.n_;
    if (!other.requires_update_) {
        requires_update_ = false;
        wwt_ = other.wwt_;
    }
    return *this;
}

wavelet_tree_stat& wavelet_tree_stat::operator=(wavelet_tree_stat&& other) noexcept {
    int_vector_ = std::move(other.int_vector_);
    n_ = other.n_;
    if (!other.requires_update_) {
        requires_update_ = false;
        wwt_ = std::move(other.wwt_);
    }
    return *this;
}

bool wavelet_tree_stat::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        int_vector_.load(in);
        wwt_.load(in);
        n_ = int_vector_.size();
        // we assume the wwt_ has been build before serialization
        requires_update_ = false;
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load wavelet_tree_stat." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void wavelet_tree_stat::serialize(std::ostream &out) const {
    if (requires_update_)
        const_cast<wavelet_tree_stat*>(this)->init_wt();

    int_vector_.serialize(out);
    wwt_.serialize(out);
}

void wavelet_tree_stat::set(uint64_t id, uint64_t val) {
    if (int_vector_[id] == val)
        return;

    if (!requires_update_) {
        std::unique_lock<std::mutex> lock(mu_);
        wwt_ = decltype(wwt_)();
        requires_update_ = true;
    }
    int_vector_[id] = val;
}

void wavelet_tree_stat::insert(uint64_t id, uint64_t val) {
    assert(id <= size());

    if (!requires_update_) {
        requires_update_ = true;
        wwt_ = decltype(wwt_)();
    }
    if (n_ == size()) {
        int_vector_.resize(2 * n_ + 1);
    }
    n_++;
    if (size() > 1)
        std::copy_backward(int_vector_.begin() + id,
                           int_vector_.begin() + n_ - 1,
                           int_vector_.begin() + n_);
    int_vector_[id] = val;
}

void wavelet_tree_stat::remove(uint64_t id) {
    assert(id < size());

    if (!requires_update_) {
        requires_update_ = true;
        wwt_ = decltype(wwt_)();
    }
    if (this->size() > 1)
        std::copy(int_vector_.begin() + id + 1,
                  int_vector_.begin() + n_,
                  int_vector_.begin() + id);
    n_--;
}

uint64_t wavelet_tree_stat::rank(uint64_t c, uint64_t i) const {
    if (requires_update_)
        const_cast<wavelet_tree_stat*>(this)->init_wt();

    return wwt_.rank(std::min(i + 1, size()), c);
}

uint64_t wavelet_tree_stat::select(uint64_t c, uint64_t i) const {
    assert(i > 0 && size() > 0);

    if (requires_update_)
        const_cast<wavelet_tree_stat*>(this)->init_wt();

    assert(i <= rank(c, size() - 1));
    return wwt_.select(i, c);
}

uint64_t wavelet_tree_stat::operator[](uint64_t id) const {
    assert(id < size());
    return int_vector_[id];
}

uint64_t wavelet_tree_stat::next(uint64_t pos, uint64_t value) const {
    assert(pos < size());

    return ::next(*this, pos, value, MAX_ITER_WAVELET_TREE_STAT);
}

uint64_t wavelet_tree_stat::prev(uint64_t pos, uint64_t value) const {
    assert(pos < size());

    return ::prev(*this, pos, value, MAX_ITER_WAVELET_TREE_STAT);
}

void wavelet_tree_stat::clear() {
    int_vector_ = decltype(int_vector_)(0, 0, int_vector_.width());
    wwt_ = decltype(wwt_)();
    requires_update_ = false;
    n_ = 0;
}

sdsl::int_vector<> wavelet_tree_stat::to_vector() const {
    sdsl::int_vector<> vector = int_vector_;
    vector.resize(n_);
    return vector;
}

void wavelet_tree_stat::init_wt() {
    std::unique_lock<std::mutex> lock(mu_);

    if (!requires_update_)
        return;

    int_vector_.resize(n_);
    wwt_ = decltype(wwt_)(int_vector_);

    requires_update_ = false;
}


///////////////////////////////////////////
// wavelet_tree_dyn libmaus2 rank/select //
///////////////////////////////////////////

template <class Vector>
BitBTree* initialize_tree(const Vector &vector, uint64_t b) {
    uint64_t n = vector.size();

    // compute total offsets for the individual bins
    std::vector<uint64_t> offsets((1ull << (b - 1)) - 1, 0);
    uint64_t v, m, o;
    for (uint64_t i = 0; i < n; ++i) {
        m = (1ull << (b - 1));
        v = static_cast<uint64_t>(vector[i]);
        o = 0;
        for (uint64_t ib = 1; ib < b; ++ib) {
            bool const bit = m & v;
            if (!bit) {
                offsets.at(o) += 1;
            }
            o = 2 * o + 1 + bit;
            m >>= 1;
        }
    }

    auto *tmp = new BitBTree(n * b, false);
    std::vector<uint64_t> upto_offsets((1ull << (b - 1)) - 1, 0);

    uint64_t p, co;
    bool bit;
    for (uint64_t i = 0; i < n; ++i) {
        m = (1ull << (b - 1));
        v = (uint64_t) vector[i];
        o = 0;
        p = i;
        co = 0;
        for (uint64_t ib = 0; ib < b - 1; ++ib) {
            bit = m & v;
            if (bit) {
                tmp->setBitQuick(ib * n + p + co, true);
                co += offsets.at(o);
                p -= upto_offsets.at(o);
            } else {
                p -= (p - upto_offsets.at(o));
                upto_offsets.at(o) += 1;
            }
            //std::cerr << "o: " << o << " offset[o]: " << offsets.at(o) << std::endl;
            o = 2 * o + 1 + bit;
            m >>= 1;
        }
        bit = m & v;
        if (bit) {
            //std::cerr << "b - 1: " << b - 1 << " n: " << n << " p: " << p << " co: " << co << std::endl;
            tmp->setBitQuick((b - 1) * n + p + co, true);
        }
    }

    return tmp;
}

wavelet_tree_dyn::wavelet_tree_dyn(uint8_t logsigma)
      : wwt_(new typename decltype(wwt_)::element_type(logsigma)) {}

template <class Vector>
wavelet_tree_dyn::wavelet_tree_dyn(uint8_t logsigma, const Vector &vector) {
    wwt_.reset(new typename decltype(wwt_)::element_type(
        initialize_tree(vector, logsigma), logsigma, vector.size()
    ));
}

template
wavelet_tree_dyn::wavelet_tree_dyn(uint8_t, const sdsl::int_vector<> &);

template
wavelet_tree_dyn::wavelet_tree_dyn(uint8_t, const std::vector<uint8_t> &);

template
wavelet_tree_dyn::wavelet_tree_dyn(uint8_t, const std::vector<int> &);

template
wavelet_tree_dyn::wavelet_tree_dyn(uint8_t, const std::vector<uint64_t> &);

wavelet_tree_dyn::wavelet_tree_dyn(const wavelet_tree_dyn &other) {
    *this = other;
}

wavelet_tree_dyn::wavelet_tree_dyn(wavelet_tree_dyn&& other) noexcept {
    *this = std::move(other);
}

// TODO: copy constructor not defined in libmaus2
wavelet_tree_dyn& wavelet_tree_dyn::operator=(const wavelet_tree_dyn &other) {
    wwt_.reset(new typename decltype(wwt_)::element_type(
        initialize_tree(*other.wwt_, other.wwt_->b), other.wwt_->b, other.size()
    ));
    return *this;
}

wavelet_tree_dyn& wavelet_tree_dyn::operator=(wavelet_tree_dyn&& other) noexcept {
    wwt_ = std::move(other.wwt_);
    return *this;
}

uint64_t wavelet_tree_dyn::rank(uint64_t c, uint64_t i) const {
    return size() > 0
            ? wwt_->rank(c, std::min(i, size() - 1))
            : 0;
}

uint64_t wavelet_tree_dyn::select(uint64_t c, uint64_t i) const {
    assert(i > 0 && size() > 0 && i <= rank(c, size() - 1));
    return wwt_->select(c, i - 1);
}

uint64_t wavelet_tree_dyn::operator[](uint64_t id) const {
    assert(id < size());
    return (*wwt_)[id];
}

uint64_t wavelet_tree_dyn::next(uint64_t pos, uint64_t value) const {
    assert(pos < size());

    return ::next(*this, pos, value, MAX_ITER_WAVELET_TREE_DYN);
}

uint64_t wavelet_tree_dyn::prev(uint64_t pos, uint64_t value) const {
    assert(pos < size());

    return ::prev(*this, pos, value, MAX_ITER_WAVELET_TREE_DYN);
}

void wavelet_tree_dyn::set(uint64_t id, uint64_t val) {
    remove(id);
    insert(id, val);
}

void wavelet_tree_dyn::insert(uint64_t id, uint64_t val) {
    assert(id <= size());
    wwt_->insert(val, id);
}

void wavelet_tree_dyn::remove(uint64_t id) {
    assert(id < size());
    wwt_->remove(id);
}

void wavelet_tree_dyn::serialize(std::ostream &out) const {
    wwt_->serialise(out);
}

bool wavelet_tree_dyn::load(std::istream &in) {
    if (!in.good())
        return false;

    //TODO: catch reading errors
    wwt_.reset(new libmaus2::wavelet::DynamicWaveletTree<6, 64>(in));
    return true;
}

void wavelet_tree_dyn::clear() {
    wwt_.reset(new libmaus2::wavelet::DynamicWaveletTree<6, 64>(wwt_->b));
}

sdsl::int_vector<> wavelet_tree_dyn::to_vector() const {
    size_t b = wwt_->b;
    size_t n = size();

    std::vector<uint64_t> offsets(1ull << b, 0);
    std::queue<uint64_t> blocks;
    std::queue<uint64_t> new_blocks;

    blocks.push(n);

    size_t pos = 0;
    uint64_t o = 0;

    for (size_t ib = 0; ib < b; ++ib) {
        while (!blocks.empty()) {
            uint64_t cnt = blocks.front();
            blocks.pop();
            uint64_t epos = pos + cnt;
            for ( ; pos < epos; ++pos) {
                offsets.at(o) += !wwt_->R->operator[](pos);
            }
            if (ib < b - 1) {
                new_blocks.push(offsets.at(o));
                new_blocks.push(cnt - offsets.at(o));
            }
            o++;
        }
        if (ib < b - 1)
            blocks.swap(new_blocks);
    }
    //std::cerr << "R size: " << R->size() << std::endl;

    sdsl::int_vector<> vector(n, 0, wwt_->b);

    bool bit;
    std::vector<uint64_t> upto_offsets((1ull << (b - 1)) - 1, 0);
    uint64_t p, co, v, m;
    for (size_t i = 0; i < n; ++i) {
        m = 1ull << (b - 1);
        //v = (uint64_t) vector.at(i);
        v = 0;
        o = 0;
        p = i;
        co = 0;
        for (size_t ib = 0; ib < b - 1; ++ib) {
            bit = wwt_->R->operator[](ib * n + p + co);
            if (bit) {
                v |= m;
                co += offsets.at(o);
                p -= upto_offsets.at(o);
            } else {
                p -= (p - upto_offsets.at(o));
                upto_offsets.at(o) += 1;
            }
            o = 2 * o + 1 + bit;
            m >>= 1;
        }
        bit = wwt_->R->operator[]((b - 1) * n + p + co);
        if (bit) {
            v |= m;
        }
        vector[i] = v;
    }
    return vector;
}


////////////////////////////////////////////////////
// wavelet_tree_small sdsl rank/select, wt access //
////////////////////////////////////////////////////

template <class Vector>
wavelet_tree_small::wavelet_tree_small(uint8_t logsigma, const Vector &vector)
      : wwt_(pack_vector(vector, logsigma)), logsigma_(logsigma) {}

template
wavelet_tree_small::wavelet_tree_small(uint8_t logsigma,
                                       const sdsl::int_vector<> &vector);
template
wavelet_tree_small::wavelet_tree_small(uint8_t logsigma,
                                       const std::vector<uint8_t> &vector);
template
wavelet_tree_small::wavelet_tree_small(uint8_t logsigma,
                                       const std::vector<uint64_t> &vector);
template
wavelet_tree_small::wavelet_tree_small(uint8_t logsigma,
                                       const std::vector<int> &vector);

wavelet_tree_small::wavelet_tree_small(uint8_t logsigma, const sdsl::wt_huff<> &wwt)
      : wwt_(wwt), logsigma_(logsigma) {}

wavelet_tree_small::wavelet_tree_small(uint8_t logsigma, sdsl::wt_huff<>&& wwt)
      : wwt_(std::move(wwt)), logsigma_(logsigma) {}

uint64_t wavelet_tree_small::rank(uint64_t c, uint64_t i) const {
    return wwt_.rank(std::min(i + 1, size()), c);
}

uint64_t wavelet_tree_small::select(uint64_t c, uint64_t i) const {
    assert(i > 0 && size() > 0);
    assert(i <= rank(c, size() - 1));
    return wwt_.select(i, c);
}

uint64_t wavelet_tree_small::operator[](uint64_t id) const {
    assert(id < size());
    return wwt_[id];
}

uint64_t wavelet_tree_small::next(uint64_t pos, uint64_t value) const {
    assert(pos < size());

    return ::next(*this, pos, value, MAX_ITER_WAVELET_TREE_SMALL);
}

uint64_t wavelet_tree_small::prev(uint64_t pos, uint64_t value) const {
    assert(pos < size());

    return ::prev(*this, pos, value, MAX_ITER_WAVELET_TREE_SMALL);
}

void wavelet_tree_small::set(uint64_t, uint64_t) {
    throw std::runtime_error("Not supported");
}

void wavelet_tree_small::insert(uint64_t, uint64_t) {
    throw std::runtime_error("Not supported");
}

void wavelet_tree_small::remove(uint64_t) {
    throw std::runtime_error("Not supported");
}

void wavelet_tree_small::serialize(std::ostream &out) const {
    wwt_.serialize(out);
}

bool wavelet_tree_small::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        wwt_.load(in);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load wavelet_tree_small" << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

sdsl::int_vector<> wavelet_tree_small::to_vector() const {
    sdsl::int_vector<> vector(wwt_.size(), 0, logsigma_);
    std::copy(wwt_.begin(), wwt_.end(), vector.begin());
    return vector;
}


////////////////////////////////////////////////////
// wavelet_tree_fast based on bit_vectors, static //
////////////////////////////////////////////////////

wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, uint64_t size, uint64_t value)
      : int_vector_(size, value, logsigma), bitmaps_(1 << logsigma) {
    bitmaps_.at(value) = bit_vector_stat(size, 1);
}

template <class Vector>
wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, const Vector &vector)
      : int_vector_(pack_vector(vector, logsigma)), bitmaps_(1 << logsigma) {

    for (size_t i = 0; i < vector.size(); ++i) {

        uint64_t value = vector[i];

        assert(value < static_cast<uint64_t>(bitmaps_.size()));

        if (!bitmaps_[value].size())
            bitmaps_[value] = bit_vector_stat(vector.size(), 0);

        bitmaps_[value].set(i, 1);
    }
}

template wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, const sdsl::int_vector<> &vector);
template wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, const std::vector<uint8_t> &vector);
template wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, const std::vector<uint64_t> &vector);
template wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, const std::vector<int> &vector);

uint64_t wavelet_tree_fast::rank(uint64_t c, uint64_t i) const {
    assert(c < bitmaps_.size());
    return bitmaps_[c].size() ? bitmaps_[c].rank1(i) : 0;
}

uint64_t wavelet_tree_fast::select(uint64_t c, uint64_t i) const {
    assert(i > 0 && size() > 0);
    assert(c < bitmaps_.size());
    assert(i <= rank(c, size() - 1));
    return bitmaps_[c].select1(i);
}

uint64_t wavelet_tree_fast::operator[](uint64_t id) const {
    assert(id < size());
    return int_vector_[id];
}

uint64_t wavelet_tree_fast::next(uint64_t pos, uint64_t value) const {
    assert(pos < size());

    return ::next(*this, pos, value, MAX_ITER_WAVELET_TREE_STAT);
}

uint64_t wavelet_tree_fast::prev(uint64_t pos, uint64_t value) const {
    assert(pos < size());

    return ::prev(*this, pos, value, MAX_ITER_WAVELET_TREE_STAT);
}

void wavelet_tree_fast::set(uint64_t id, uint64_t val) {
    assert(id < int_vector_.size());

    if (int_vector_[id] == val)
        return;

    int_vector_[id] = val;

    for (size_t c = 0; c < bitmaps_.size(); ++c) {
        if (c != val && bitmaps_[c].size())
            bitmaps_[c].set(id, 0);
    }

    if (!bitmaps_[val].size())
        bitmaps_[val] = bit_vector_stat(int_vector_.size(), 0);

    bitmaps_[val].set(id, 1);
}

void wavelet_tree_fast::insert(uint64_t id, uint64_t val) {
    assert(id <= int_vector_.size());

    int_vector_.resize(int_vector_.size() + 1);

    std::copy_backward(int_vector_.begin() + id,
                       int_vector_.end() - 1,
                       int_vector_.end());
    int_vector_[id] = val;

    assert(val < bitmaps_.size());

    for (size_t c = 0; c < bitmaps_.size(); ++c) {
        if (c != val && bitmaps_[c].size())
            bitmaps_[c].insert_bit(id, 0);
    }

    if (!bitmaps_[val].size()) {
        bitmaps_[val] = bit_vector_stat(int_vector_.size(), 0);
        bitmaps_[val].set(id, 1);
    } else {
        bitmaps_[val].insert_bit(id, 1);
    }
}

void wavelet_tree_fast::remove(uint64_t id) {
    assert(id < int_vector_.size());

    std::copy(int_vector_.begin() + id + 1,
              int_vector_.end(),
              int_vector_.begin() + id);

    int_vector_.resize(int_vector_.size() - 1);

    for (size_t c = 0; c < bitmaps_.size(); ++c) {
        if (bitmaps_[c].size())
            bitmaps_[c].delete_bit(id);
    }
}

void wavelet_tree_fast::serialize(std::ostream &out) const {
    int_vector_.serialize(out);
    for (const auto &v : bitmaps_) {
        v.serialize(out);
    }
}

bool wavelet_tree_fast::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        int_vector_.load(in);

        bitmaps_.resize(1 << logsigma());

        for (auto &v : bitmaps_) {
            if (!v.load(in))
                return false;
        }

        return true;

    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load wavelet_tree_fast" << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void wavelet_tree_fast::clear() {
    int_vector_.resize(0);
    bitmaps_.assign(bitmaps_.size(), bit_vector_stat());
}
