#include "wavelet_tree.hpp"

#include <cassert>


/////////////////////////////////
// wavelet_tree shared methods //
/////////////////////////////////

template <class Vector>
sdsl::int_vector<> dump_int_vector(const Vector &v, uint64_t logsigma) {
    sdsl::int_vector<> iv(v.size(), 0, logsigma);

    for (uint64_t i = 0; i < v.size(); ++i) {
        iv[i] = v[i];
    }

    return iv;
}

template <class WaveletTree>
WaveletTree wavelet_tree::convert_to(uint64_t logsigma) {
    if (dynamic_cast<WaveletTree*>(this)) {
        return WaveletTree(dynamic_cast<WaveletTree&&>(*this));
    }
    return WaveletTree(logsigma, dump_int_vector(*this, logsigma));
}
template wavelet_tree_dyn wavelet_tree::convert_to(uint64_t);

template <>
wavelet_tree_stat wavelet_tree::convert_to(uint64_t logsigma) {
    if (dynamic_cast<wavelet_tree_stat*>(this)) {
        return wavelet_tree_stat(dynamic_cast<wavelet_tree_stat&&>(*this));
    } else if (dynamic_cast<wavelet_tree_small*>(this)) {
        return wavelet_tree_stat(logsigma, std::move(dynamic_cast<wavelet_tree_small&&>(*this).wwt_));
    }
    return wavelet_tree_stat(logsigma, dump_int_vector(*this, logsigma));
}

template <>
wavelet_tree_small wavelet_tree::convert_to(uint64_t logsigma){
    if (dynamic_cast<wavelet_tree_small*>(this)) {
        return wavelet_tree_small(dynamic_cast<wavelet_tree_small&&>(*this));
    } else if (dynamic_cast<wavelet_tree_stat*>(this)) {
        auto&& input = dynamic_cast<wavelet_tree_stat&&>(*this);
        if (input.requires_update_) {
            input.wwt_ = decltype(input.wwt_)();
            input.init_wt();
        }

        input.int_vector_ = decltype(input.int_vector_)();
        return wavelet_tree_small(logsigma, std::move(input.wwt_));
    }
    return wavelet_tree_small(logsigma, dump_int_vector(*this, logsigma));
}

///////////////////////////////////////////////////////////
// wavelet_tree_stat sdsl rank/select, int_vector access //
///////////////////////////////////////////////////////////

wavelet_tree_stat::wavelet_tree_stat(uint64_t logsigma,
                                     uint64_t size, uint64_t value)
      : int_vector_(size, value, logsigma), n_(size) {}

template <class Vector>
wavelet_tree_stat::wavelet_tree_stat(uint64_t logsigma, const Vector &vector)
      : int_vector_(dump_int_vector(vector, logsigma)), n_(int_vector_.size()) {}

template
wavelet_tree_stat
::wavelet_tree_stat<std::vector<uint8_t>>(uint64_t logsigma,
                                          const std::vector<uint8_t> &vector);

template
wavelet_tree_stat
::wavelet_tree_stat<std::vector<uint64_t>>(uint64_t logsigma,
                                           const std::vector<uint64_t> &vector);

template
wavelet_tree_stat
::wavelet_tree_stat<std::vector<int>>(uint64_t logsigma,
                                      const std::vector<int> &vector);

wavelet_tree_stat::wavelet_tree_stat(const wavelet_tree_stat &other)
      : int_vector_(other.int_vector_),
        wwt_(other.wwt_),
        requires_update_(false),
        n_(other.n_) {}

wavelet_tree_stat::wavelet_tree_stat(uint64_t logsigma, wt_type&& wwt)
      : int_vector_(dump_int_vector(wwt, logsigma)),
        wwt_(std::move(wwt)),
        requires_update_(false),
        n_(wwt_.size()) {}

wavelet_tree_stat::wavelet_tree_stat(wavelet_tree_stat&& other)
      : int_vector_(std::move(other.int_vector_)),
        wwt_(std::move(other.wwt_)),
        requires_update_(false),
        n_(other.n_) {}

wavelet_tree_stat::wavelet_tree_stat(sdsl::int_vector<>&& vector) {
    *this = std::move(vector);
}

wavelet_tree_stat& wavelet_tree_stat::operator=(sdsl::int_vector<>&& vector) {
    wwt_ = decltype(wwt_)();
    int_vector_ = std::move(vector);
    n_ = int_vector_.size();
    requires_update_ = true;
    return *this;
}

bool wavelet_tree_stat::deserialise(std::istream &in) {
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

void wavelet_tree_stat::serialise(std::ostream &out) const {
    if (requires_update_)
        const_cast<wavelet_tree_stat*>(this)->init_wt();

    int_vector_.serialize(out);
    wwt_.serialize(out);
}

void wavelet_tree_stat::insert(uint64_t id, uint64_t val) {
    wwt_ = decltype(wwt_)();
    if (n_ == size()) {
        int_vector_.resize(2 * n_ + 1);
    }
    n_++;
    if (size() > 1)
        std::copy_backward(int_vector_.begin() + id,
                           int_vector_.begin() + n_ - 1,
                           int_vector_.begin() + n_);
    int_vector_[id] = val;
    requires_update_ = true;
}

void wavelet_tree_stat::remove(uint64_t id) {
    wwt_ = decltype(wwt_)();
    if (this->size() > 1)
        std::copy(int_vector_.begin() + id + 1,
                  int_vector_.begin() + n_,
                  int_vector_.begin() + id);
    n_--;
    requires_update_ = true;
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

std::vector<uint64_t> wavelet_tree_stat::to_vector() const {
    std::vector<uint64_t> result(size());
    for (uint64_t i = 0; i < size(); ++i) {
        result[i] = int_vector_[i];
    }
    return result;
}

void wavelet_tree_stat::set(uint64_t id, uint64_t val) {
    requires_update_ = true;
    int_vector_[id] = val;
}

void wavelet_tree_stat::init_wt() {
    int_vector_.resize(n_);

    // release obsolete wavelet tree
    wwt_ = decltype(wwt_)();
    // initialize new wavelet tree
    wwt_ = decltype(wwt_)(int_vector_);
    requires_update_ = false;
}

////////////////////////////////////////////////////
// wavelet_tree_small sdsl rank/select, wt access //
////////////////////////////////////////////////////

wavelet_tree_small::wavelet_tree_small(uint64_t logsigma,
                                       uint64_t size, uint64_t value) {
    sdsl::int_vector<> int_vector_(size, value, logsigma);
    //TODO fix constructor in sdsl
    wwt_ = decltype(wwt_)(int_vector_);
}

template <class Vector>
wavelet_tree_small::wavelet_tree_small(uint64_t logsigma, const Vector &vector) {
    auto iv = dump_int_vector(vector, logsigma);
    //TODO fix constructor in sdsl
    wwt_ = decltype(wwt_)(iv);
}

template
wavelet_tree_small
::wavelet_tree_small<std::vector<uint8_t>>(uint64_t logsigma,
                                           const std::vector<uint8_t> &vector);

template
wavelet_tree_small
::wavelet_tree_small<std::vector<uint64_t>>(uint64_t logsigma,
                                            const std::vector<uint64_t> &vector);

template
wavelet_tree_small
::wavelet_tree_small<std::vector<int>>(uint64_t logsigma,
                                       const std::vector<int> &vector);

wavelet_tree_small::wavelet_tree_small(const wavelet_tree_small &other)
      : wwt_(other.wwt_) {}

wavelet_tree_small::wavelet_tree_small(wt_type&& wwt)
      : wwt_(std::move(wwt)) {}

wavelet_tree_small::wavelet_tree_small(wavelet_tree_small&& other)
      : wwt_(std::move(other.wwt_)) {}

wavelet_tree_small::wavelet_tree_small(sdsl::int_vector<>&& vector) {
    *this = std::move(vector);
}

wavelet_tree_small& wavelet_tree_small::operator=(sdsl::int_vector<>&& vector) {
    //TODO: this should be std::move, fix line 257 in wt_pc.hpp
    wwt_ = decltype(wwt_)();
    wwt_ = decltype(wwt_)(vector);
    return *this;
}

bool wavelet_tree_small::deserialise(std::istream &in) {
    if (!in.good())
        return false;

    try {
        wwt_.load(in);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load wavelet_tree_small." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void wavelet_tree_small::serialise(std::ostream &out) const {
    wwt_.serialize(out);
}

void wavelet_tree_small::insert(uint64_t, uint64_t) {
    throw std::runtime_error("Not supported\n");
}

void wavelet_tree_small::remove(uint64_t) {
    throw std::runtime_error("Not supported\n");
}

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

std::vector<uint64_t> wavelet_tree_small::to_vector() const {
    return std::vector<uint64_t>(wwt_.begin(), wwt_.end());
}

void wavelet_tree_small::set(uint64_t, uint64_t) {
    throw std::runtime_error("Not supported\n");
}

///////////////////////////////////////////
// wavelet_tree_dyn libmaus2 rank/select //
///////////////////////////////////////////

template <class Vector>
libmaus2::bitbtree::BitBTree<6, 64>* initialize_tree(const Vector &W_stat,
                                                     const uint64_t b) {
    uint64_t const n = W_stat.size();

    // compute total offsets for the individual bins
    std::vector<uint64_t> offsets((1ull << (b - 1)) - 1, 0);
    uint64_t v, m, o;
    for (uint64_t i = 0; i < n; ++i) {
        m = (1ull << (b - 1));
        v = static_cast<uint64_t>(W_stat[i]);
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

    auto *tmp = new libmaus2::bitbtree::BitBTree<6, 64>(n * b, false);
    std::vector<uint64_t> upto_offsets((1ull << (b - 1)) - 1, 0);

    uint64_t p, co;
    bool bit;
    for (uint64_t i = 0; i < n; ++i) {
        m = (1ull << (b - 1));
        v = (uint64_t) W_stat[i];
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


template <class Vector>
wavelet_tree_dyn::wavelet_tree_dyn(uint64_t logsigma,
                                   const Vector &W_stat)
    : wwt_(initialize_tree(W_stat, logsigma), logsigma, W_stat.size()) {}

template
wavelet_tree_dyn
::wavelet_tree_dyn<std::vector<uint8_t>>(uint64_t logsigma,
                                         const std::vector<uint8_t> &W_stat);

template
wavelet_tree_dyn
::wavelet_tree_dyn<std::vector<int>>(uint64_t logsigma,
                                     const std::vector<int> &W_stat);

template
wavelet_tree_dyn
::wavelet_tree_dyn<std::vector<uint64_t>>(uint64_t logsigma,
                                          const std::vector<uint64_t> &W_stat);

// TODO: copy constructor not defined in libmaus2
wavelet_tree_dyn::wavelet_tree_dyn(const wavelet_tree_dyn &other)
      : wwt_(initialize_tree(other.wwt_, other.wwt_.b), other.wwt_.b, other.size()) {}

wavelet_tree_dyn::wavelet_tree_dyn(wavelet_tree_dyn&& other)
      : wwt_(std::move(other.wwt_)) {}

bool wavelet_tree_dyn::deserialise(std::istream &in) {
    if (!in.good())
        return false;

    //TODO: catch reading errors
    wwt_.R =
        libmaus2::wavelet::DynamicWaveletTree<6, 64>::loadBitBTree(in);
    const_cast<uint64_t&>(wwt_.b) =
        libmaus2::util::NumberSerialisation::deserialiseNumber(in);
    wwt_.n =
        libmaus2::util::NumberSerialisation::deserialiseNumber(in);
    return true;
}

void wavelet_tree_dyn::serialise(std::ostream &out) const {
    wwt_.serialise(out);
}

void wavelet_tree_dyn::insert(uint64_t id, uint64_t val) {
    wwt_.insert(val, id);
}

void wavelet_tree_dyn::remove(uint64_t id) {
    wwt_.remove(id);
}

uint64_t wavelet_tree_dyn::rank(uint64_t c, uint64_t i) const {
    return size() > 0
            ? wwt_.rank(c, std::min(i, size() - 1))
            : 0;
}

uint64_t wavelet_tree_dyn::select(uint64_t c, uint64_t i) const {
    assert(i > 0 && size() > 0 && i <= rank(c, size() - 1));
    return wwt_.select(c, i - 1);
}

uint64_t wavelet_tree_dyn::operator[](uint64_t id) const {
    assert(id < size());
    return wwt_[id];
}

void wavelet_tree_dyn::set(uint64_t id, uint64_t val) {
    remove(id);
    insert(id, val);
}

std::vector<uint64_t> wavelet_tree_dyn::to_vector() const {
    std::vector<uint64_t> W_stat;

    size_t n = size();

    const size_t b = wwt_.b;

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
                offsets.at(o) += !get_bit_raw(pos);
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

    bool bit;
    std::vector<uint64_t> upto_offsets((1ull << (b - 1)) - 1, 0);
    uint64_t p, co, v, m;
    for (size_t i = 0; i < n; ++i) {
        m = 1ull << (b - 1);
        //v = (uint64_t) W_stat.at(i);
        v = 0;
        o = 0;
        p = i;
        co = 0;
        for (size_t ib = 0; ib < b - 1; ++ib) {
            bit = get_bit_raw(ib * n + p + co);
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
        bit = get_bit_raw((b - 1) * n + p + co);
        if (bit) {
            v |= m;
        }
        W_stat.push_back(v);
    }
    return W_stat;
}

bool wavelet_tree_dyn::get_bit_raw(uint64_t id) const {
    return wwt_.R->operator[](id);
}
