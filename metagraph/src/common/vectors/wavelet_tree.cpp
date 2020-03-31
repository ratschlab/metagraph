#include "wavelet_tree.hpp"

#include <cassert>

#include "common/algorithms.hpp"
#include "common/serialization.hpp"

typedef wavelet_tree::TAlphabet TAlphabet;

const size_t MAX_ITER_WAVELET_TREE_STAT = 1000;
const size_t MAX_ITER_WAVELET_TREE_DYN = 0;
const size_t MAX_ITER_WAVELET_TREE_SMALL = 20;


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

bool wavelet_tree::operator==(const wavelet_tree &other) const {
    if (size() != other.size())
        return false;

    const wavelet_tree *bitmap_p[] = { this, &other };
    const sdsl::int_vector<> *bv_p[] = { nullptr, nullptr };

    for (int i : { 0, 1 }) {
        if (dynamic_cast<const wavelet_tree_stat*>(bitmap_p[i])) {
            bv_p[i] = &dynamic_cast<const wavelet_tree_stat&>(*bitmap_p[i]).data();
        } else if (dynamic_cast<const wavelet_tree_fast*>(bitmap_p[i])) {
            bv_p[i] = &dynamic_cast<const wavelet_tree_fast&>(*bitmap_p[i]).data();
        }
    }

    if (bv_p[0] && bv_p[1])
        return *bv_p[0] == *bv_p[1];

    const uint64_t end = size();
    for (uint64_t i = 0; i < end; ++i) {
        if ((*this)[i] != other[i])
            return false;
    }
    return true;
}

template <class WaveletTree>
WaveletTree wavelet_tree::convert_to() {
    if (dynamic_cast<WaveletTree*>(this)) {
        return WaveletTree(dynamic_cast<WaveletTree&&>(*this));

    } else if (dynamic_cast<wavelet_tree_small*>(this)
                && std::is_same_v<WaveletTree, wavelet_tree_stat>) {
        return WaveletTree(
            logsigma(),
            std::move(dynamic_cast<wavelet_tree_small*>(this)->wwt_)
        );

    } else if (dynamic_cast<wavelet_tree_stat*>(this)
                && std::is_same_v<WaveletTree, wavelet_tree_small>) {
        return WaveletTree(
            logsigma(),
            std::move(dynamic_cast<wavelet_tree_stat*>(this)->wwt_)
        );

    } else {
        sdsl::int_vector<> vector;
        if (auto *wt_stat = dynamic_cast<wavelet_tree_stat*>(this)) {
            vector = std::move(wt_stat->int_vector_);
        } else if (auto *wt_fast = dynamic_cast<wavelet_tree_fast*>(this)) {
            vector = std::move(wt_fast->int_vector_);
        } else {
            vector = to_vector();
        }
        clear();

        return WaveletTree(vector.width(), std::move(vector));
    }
}

template wavelet_tree_dyn wavelet_tree::convert_to<wavelet_tree_dyn>();
template wavelet_tree_stat wavelet_tree::convert_to<wavelet_tree_stat>();
template wavelet_tree_small wavelet_tree::convert_to<wavelet_tree_small>();
template wavelet_tree_fast wavelet_tree::convert_to<wavelet_tree_fast>();


template <typename Vector>
inline uint64_t next(const Vector &v,
                     uint64_t pos,
                     TAlphabet c,
                     size_t num_steps) {
    assert(pos < v.size());

    if (v[pos] == c)
        return pos;

    for (size_t t = 1; t < num_steps; ++t) {
        if (pos + t == v.size() || v[pos + t] == c)
            return pos + t;
    }

    uint64_t rk = v.rank(c, pos) + 1;
    return rk <= v.count(c)
            ? v.select(c, rk)
            : v.size();
}

template <typename Vector>
inline uint64_t prev(const Vector &v,
                     uint64_t pos,
                     TAlphabet c,
                     size_t num_steps) {
    assert(pos < v.size());

    for (size_t t = 0; t < num_steps; ++t, --pos) {
        if (v[pos] == c)
            return pos;

        if (pos == 0)
            return v.size();
    }

    uint64_t rk = v.rank(c, pos);
    return rk ? v.select(c, rk)
              : v.size();
}

///////////////////////////////////////////////////////////
// wavelet_tree_stat sdsl rank/select, int_vector access //
///////////////////////////////////////////////////////////

wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma,
                                     uint64_t size, TAlphabet c)
      : int_vector_(size, c, logsigma),
        wwt_(int_vector_),
        count_(1 << logsigma, 0) {
    count_[c] = size;
}

template <class Vector>
wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma, const Vector &vector)
      : int_vector_(pack_vector(vector, logsigma)),
        wwt_(int_vector_),
        count_(1 << logsigma) {
    for (TAlphabet c = 0; c < count_.size(); ++c) {
        count_[c] = rank(c, size());
    }
}

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
    wwt_ = decltype(wwt_)(int_vector_),

    count_.resize(1 << this->logsigma());
    for (TAlphabet c = 0; c < count_.size(); ++c) {
        count_[c] = rank(c, size());
    }
}

wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma, sdsl::wt_huff<>&& wwt)
      : int_vector_(pack_vector(wwt, logsigma)),
        wwt_(std::move(wwt)),
        count_(1 << logsigma) {
    for (TAlphabet c = 0; c < count_.size(); ++c) {
        count_[c] = rank(c, size());
    }
}

wavelet_tree_stat::wavelet_tree_stat(const wavelet_tree_stat &other) {
    *this = other;
}

wavelet_tree_stat::wavelet_tree_stat(wavelet_tree_stat&& other) noexcept {
    *this = std::move(other);
}

wavelet_tree_stat& wavelet_tree_stat::operator=(const wavelet_tree_stat &other) {
    int_vector_ = other.int_vector_;
    wwt_ = other.wwt_;
    count_ = other.count_;
    return *this;
}

wavelet_tree_stat& wavelet_tree_stat::operator=(wavelet_tree_stat&& other) noexcept {
    int_vector_ = std::move(other.int_vector_);
    wwt_ = std::move(other.wwt_);
    count_ = std::move(other.count_);
    return *this;
}

bool wavelet_tree_stat::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        int_vector_.load(in);
        wwt_.load(in);

        count_.resize(1 << logsigma());
        for (TAlphabet c = 0; c < count_.size(); ++c) {
            count_[c] = rank(c, size());
        }

        return true;

    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load wavelet_tree_stat." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void wavelet_tree_stat::serialize(std::ostream &out) const {
    int_vector_.serialize(out);
    wwt_.serialize(out);
}

uint64_t wavelet_tree_stat::rank(TAlphabet c, uint64_t i) const {
    assert(c < (1llu << logsigma()));
    return wwt_.rank(std::min(i + 1, size()), c);
}

uint64_t wavelet_tree_stat::select(TAlphabet c, uint64_t i) const {
    assert(i > 0 && size() > 0);
    assert(i <= rank(c, size() - 1));
    assert(c < (1llu << logsigma()));
    return wwt_.select(i, c);
}

TAlphabet wavelet_tree_stat::operator[](uint64_t i) const {
    assert(i < size());
    return int_vector_[i];
}

uint64_t wavelet_tree_stat::next(uint64_t pos, TAlphabet c) const {
    assert(pos < size());
    assert(c < (1llu << logsigma()));

    return ::next(*this, pos, c, MAX_ITER_WAVELET_TREE_STAT);
}

uint64_t wavelet_tree_stat::prev(uint64_t pos, TAlphabet c) const {
    assert(pos < size());
    assert(c < (1llu << logsigma()));

    return ::prev(*this, pos, c, MAX_ITER_WAVELET_TREE_STAT);
}

void wavelet_tree_stat::clear() {
    int_vector_ = decltype(int_vector_)(0, 0, int_vector_.width());
    wwt_ = decltype(wwt_)();
}


////////////////////////////////////////////////
// wavelet_tree_dyn xxsds/DYNAMIC rank/select //
////////////////////////////////////////////////

wavelet_tree_dyn::wavelet_tree_dyn(uint8_t logsigma)
      : dwt_(1ull << logsigma) {}

template <class Vector>
wavelet_tree_dyn::wavelet_tree_dyn(uint8_t logsigma, const Vector &vector)
      : dwt_(1ull << logsigma) {
    dwt_.push_many(1ull << logsigma, vector);
}

template wavelet_tree_dyn::wavelet_tree_dyn(uint8_t, const sdsl::int_vector<> &);
template wavelet_tree_dyn::wavelet_tree_dyn(uint8_t, const std::vector<uint8_t> &);
template wavelet_tree_dyn::wavelet_tree_dyn(uint8_t, const std::vector<int> &);
template wavelet_tree_dyn::wavelet_tree_dyn(uint8_t, const std::vector<uint64_t> &);

uint64_t wavelet_tree_dyn::rank(TAlphabet c, uint64_t i) const {
    assert(c < (1llu << logsigma()));
    return size() > 0
            ? dwt_.rank(std::min(i + 1, size()), c)
            : 0;
}

uint64_t wavelet_tree_dyn::select(TAlphabet c, uint64_t i) const {
    assert(i > 0 && size() > 0 && i <= rank(c, size() - 1));
    assert(c < (1llu << logsigma()));
    return dwt_.select(i - 1, c);
}

TAlphabet wavelet_tree_dyn::operator[](uint64_t i) const {
    assert(i < size());
    return dwt_.at(i);
}

uint64_t wavelet_tree_dyn::next(uint64_t pos, TAlphabet c) const {
    assert(pos < size());
    assert(c < (1llu << logsigma()));

    return ::next(*this, pos, c, MAX_ITER_WAVELET_TREE_DYN);
}

uint64_t wavelet_tree_dyn::prev(uint64_t pos, TAlphabet c) const {
    assert(pos < size());
    assert(c < (1llu << logsigma()));

    return ::prev(*this, pos, c, MAX_ITER_WAVELET_TREE_DYN);
}

void wavelet_tree_dyn::set(uint64_t i, TAlphabet c) {
    remove(i);
    insert(i, c);
}

void wavelet_tree_dyn::insert(uint64_t i, TAlphabet c) {
    assert(i <= size());
    dwt_.insert(i, c);
}

void wavelet_tree_dyn::remove(uint64_t i) {
    assert(i < size());
    dwt_.remove(i);
}

uint8_t wavelet_tree_dyn::logsigma() const {
    return dwt_.alphabet_size() ? sdsl::bits::hi(dwt_.alphabet_size() - 1) + 1 : 1;
}

void wavelet_tree_dyn::serialize(std::ostream &out) const {
    dwt_.serialize(out);
}

bool wavelet_tree_dyn::load(std::istream &in) {
    if (!in.good())
        return false;

    //TODO: catch reading errors
    dwt_.load(in);
    return in.good();
}

void wavelet_tree_dyn::clear() {
    dwt_ = decltype(dwt_)(dwt_.alphabet_size());
}

sdsl::int_vector<> wavelet_tree_dyn::to_vector() const {
    sdsl::int_vector<> vector(dwt_.size(), 0, logsigma());
    for(size_t i = 0; i < vector.size(); ++i) {
        vector[i] = dwt_.at(i);
    }
    return vector;
}


////////////////////////////////////////////////////
// wavelet_tree_small sdsl rank/select, wt access //
////////////////////////////////////////////////////

template <class Vector>
wavelet_tree_small::wavelet_tree_small(uint8_t logsigma, const Vector &vector)
      : wwt_(pack_vector(vector, logsigma)), logsigma_(logsigma), count_(1 << logsigma) {
    for (TAlphabet c = 0; c < count_.size(); ++c) {
        count_[c] = rank(c, size());
    }
}

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
      : wwt_(wwt), logsigma_(logsigma), count_(1 << logsigma) {
    for (TAlphabet c = 0; c < count_.size(); ++c) {
        count_[c] = rank(c, size());
    }
}

wavelet_tree_small::wavelet_tree_small(uint8_t logsigma, sdsl::wt_huff<>&& wwt)
      : wwt_(std::move(wwt)), logsigma_(logsigma), count_(1 << logsigma) {
    for (TAlphabet c = 0; c < count_.size(); ++c) {
        count_[c] = rank(c, size());
    }
}

uint64_t wavelet_tree_small::rank(TAlphabet c, uint64_t i) const {
    assert(c < (1llu << logsigma()));
    return wwt_.rank(std::min(i + 1, size()), c);
}

uint64_t wavelet_tree_small::select(TAlphabet c, uint64_t i) const {
    assert(i > 0 && size() > 0);
    assert(i <= rank(c, size() - 1));
    assert(c < (1llu << logsigma()));
    return wwt_.select(i, c);
}

TAlphabet wavelet_tree_small::operator[](uint64_t i) const {
    assert(i < size());
    return wwt_[i];
}

uint64_t wavelet_tree_small::next(uint64_t pos, TAlphabet c) const {
    assert(pos < size());
    assert(c < (1llu << logsigma()));

    return ::next(*this, pos, c, MAX_ITER_WAVELET_TREE_SMALL);
}

uint64_t wavelet_tree_small::prev(uint64_t pos, TAlphabet c) const {
    assert(pos < size());
    assert(c < (1llu << logsigma()));

    return ::prev(*this, pos, c, MAX_ITER_WAVELET_TREE_SMALL);
}

void wavelet_tree_small::serialize(std::ostream &out) const {
    wwt_.serialize(out);
}

bool wavelet_tree_small::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        wwt_.load(in);

        count_.resize(1 << logsigma());
        for (TAlphabet c = 0; c < count_.size(); ++c) {
            count_[c] = rank(c, size());
        }

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

wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, uint64_t size, TAlphabet c)
      : int_vector_(size, c, logsigma), bitmaps_(1 << logsigma) {
    bitmaps_.at(c) = bit_vector_stat(size, 1);
}

template <class Vector>
wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, const Vector &vector)
      : int_vector_(pack_vector(vector, logsigma)) {

    std::vector<sdsl::bit_vector> bitmaps(1 << logsigma);

    for (size_t i = 0; i < vector.size(); ++i) {

        TAlphabet c = vector[i];

        assert(c < static_cast<uint64_t>(bitmaps.size()));

        if (!bitmaps[c].size())
            bitmaps[c] = sdsl::bit_vector(vector.size(), 0);

        bitmaps[c][i] = 1;
    }

    for (auto &bv : bitmaps) {
        bitmaps_.emplace_back(std::move(bv));
    }
}

template wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, const sdsl::int_vector<> &vector);
template wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, const std::vector<uint8_t> &vector);
template wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, const std::vector<uint64_t> &vector);
template wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, const std::vector<int> &vector);

uint64_t wavelet_tree_fast::rank(TAlphabet c, uint64_t i) const {
    assert(c < bitmaps_.size());
    return bitmaps_[c].size() ? bitmaps_[c].rank1(i) : 0;
}

uint64_t wavelet_tree_fast::select(TAlphabet c, uint64_t i) const {
    assert(i > 0 && size() > 0);
    assert(c < bitmaps_.size());
    assert(i <= rank(c, size() - 1));
    return bitmaps_[c].select1(i);
}

TAlphabet wavelet_tree_fast::operator[](uint64_t i) const {
    assert(i < size());
    return int_vector_[i];
}

uint64_t wavelet_tree_fast::next(uint64_t pos, TAlphabet c) const {
    assert(pos < size());
    assert(c < (1llu << logsigma()));

    return ::next(*this, pos, c, MAX_ITER_WAVELET_TREE_STAT);
}

uint64_t wavelet_tree_fast::prev(uint64_t pos, TAlphabet c) const {
    assert(pos < size());
    assert(c < (1llu << logsigma()));

    return ::prev(*this, pos, c, MAX_ITER_WAVELET_TREE_STAT);
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
