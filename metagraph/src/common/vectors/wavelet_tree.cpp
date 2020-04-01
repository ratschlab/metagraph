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

bool wavelet_tree::operator==(const wavelet_tree &other) const {
    if (size() != other.size())
        return false;

    const wavelet_tree *bitmap_p[] = { this, &other };
    const sdsl::int_vector<> *bv_p[] = { nullptr, nullptr };

    for (int i : { 0, 1 }) {
        // TODO: write this properly for all possible templates
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

/* TODO: write this properly for all possible templates
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
*/
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
inline uint64_t next(const Vector &v, uint64_t i, TAlphabet c, size_t num_steps) {
    assert(i < v.size());

    if (v[i] == c)
        return i;

    for (size_t t = 1; t < num_steps; ++t) {
        if (i + t == v.size() || v[i + t] == c)
            return i + t;
    }

    uint64_t rk = v.rank(c, i) + 1;
    return rk <= v.count(c)
            ? v.select(c, rk)
            : v.size();
}

template <typename Vector>
inline uint64_t prev(const Vector &v, uint64_t i, TAlphabet c, size_t num_steps) {
    assert(i < v.size());

    for (size_t t = 0; t < num_steps; ++t, --i) {
        if (v[i] == c)
            return i;

        if (i == 0)
            return v.size();
    }

    uint64_t rk = v.rank(c, i);
    return rk ? v.select(c, rk)
              : v.size();
}

///////////////////////////////////////////////////////////////////
// wavelet_tree_sdsl_fast -- sdsl rank/select, int_vector access //
///////////////////////////////////////////////////////////////////

template <class t_wt_sdsl>
wavelet_tree_sdsl_fast<t_wt_sdsl>
::wavelet_tree_sdsl_fast(uint8_t logsigma, uint64_t size, TAlphabet c)
      : int_vector_(size, c, logsigma),
        wwt_(int_vector_),
        count_(1 << logsigma, 0) {
    count_[c] = size;
}

template <class t_wt_sdsl>
wavelet_tree_sdsl_fast<t_wt_sdsl>
::wavelet_tree_sdsl_fast(uint8_t logsigma, sdsl::int_vector<>&& vector) {
    if (vector.width() == logsigma) {
        int_vector_ = std::move(vector);
    } else {
        int_vector_ = pack_vector(vector, logsigma);
    }
    wwt_ = t_wt_sdsl(int_vector_),

    count_.resize(1 << this->logsigma());
    for (TAlphabet c = 0; c < count_.size(); ++c) {
        count_[c] = rank(c, size());
    }
}

template <class t_wt_sdsl>
wavelet_tree_sdsl_fast<t_wt_sdsl>
::wavelet_tree_sdsl_fast(uint8_t logsigma, t_wt_sdsl&& wwt)
      : int_vector_(pack_vector(wwt, logsigma)),
        wwt_(std::move(wwt)),
        count_(1 << logsigma) {
    for (TAlphabet c = 0; c < count_.size(); ++c) {
        count_[c] = rank(c, size());
    }
}

template <class t_wt_sdsl>
bool wavelet_tree_sdsl_fast<t_wt_sdsl>::load(std::istream &in) {
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
        std::cerr << "ERROR: Not enough memory to load wavelet_tree_sdsl_fast." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

template <class t_wt_sdsl>
void wavelet_tree_sdsl_fast<t_wt_sdsl>::serialize(std::ostream &out) const {
    int_vector_.serialize(out);
    wwt_.serialize(out);
}

template <class t_wt_sdsl>
uint64_t wavelet_tree_sdsl_fast<t_wt_sdsl>::rank(TAlphabet c, uint64_t i) const {
    assert(c < (1llu << logsigma()));
    return wwt_.rank(std::min(i + 1, size()), c);
}

template <class t_wt_sdsl>
uint64_t wavelet_tree_sdsl_fast<t_wt_sdsl>::select(TAlphabet c, uint64_t i) const {
    assert(i > 0 && size() > 0);
    assert(i <= rank(c, size() - 1));
    assert(c < (1llu << logsigma()));
    return wwt_.select(i, c);
}

template <class t_wt_sdsl>
TAlphabet wavelet_tree_sdsl_fast<t_wt_sdsl>::operator[](uint64_t i) const {
    assert(i < size());
    return int_vector_[i];
}

template <class t_wt_sdsl>
uint64_t wavelet_tree_sdsl_fast<t_wt_sdsl>::next(uint64_t i, TAlphabet c) const {
    assert(i < size());
    assert(c < (1llu << logsigma()));

    return ::next(*this, i, c, MAX_ITER_WAVELET_TREE_STAT);
}

template <class t_wt_sdsl>
uint64_t wavelet_tree_sdsl_fast<t_wt_sdsl>::prev(uint64_t i, TAlphabet c) const {
    assert(i < size());
    assert(c < (1llu << logsigma()));

    return ::prev(*this, i, c, MAX_ITER_WAVELET_TREE_STAT);
}

template <class t_wt_sdsl>
void wavelet_tree_sdsl_fast<t_wt_sdsl>::clear() {
    int_vector_ = sdsl::int_vector<>(0, 0, int_vector_.width());
    wwt_ = t_wt_sdsl();
}

template class wavelet_tree_sdsl_fast<sdsl::wt_huff<>>;


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

uint64_t wavelet_tree_dyn::next(uint64_t i, TAlphabet c) const {
    assert(i < size());
    assert(c < (1llu << logsigma()));

    return ::next(*this, i, c, MAX_ITER_WAVELET_TREE_DYN);
}

uint64_t wavelet_tree_dyn::prev(uint64_t i, TAlphabet c) const {
    assert(i < size());
    assert(c < (1llu << logsigma()));

    return ::prev(*this, i, c, MAX_ITER_WAVELET_TREE_DYN);
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
// wavelet_tree_sdsl: sdsl rank/select, wt access //
////////////////////////////////////////////////////

template <class t_wt_sdsl>
wavelet_tree_sdsl<t_wt_sdsl>
::wavelet_tree_sdsl(uint8_t logsigma, const t_wt_sdsl &wwt)
      : wwt_(wwt), logsigma_(logsigma), count_(1 << logsigma) {
    for (TAlphabet c = 0; c < count_.size(); ++c) {
        count_[c] = rank(c, size());
    }
}

template <class t_wt_sdsl>
wavelet_tree_sdsl<t_wt_sdsl>
::wavelet_tree_sdsl(uint8_t logsigma, t_wt_sdsl&& wwt)
      : wwt_(std::move(wwt)), logsigma_(logsigma), count_(1 << logsigma) {
    for (TAlphabet c = 0; c < count_.size(); ++c) {
        count_[c] = rank(c, size());
    }
}

template <class t_wt_sdsl>
uint64_t wavelet_tree_sdsl<t_wt_sdsl>::rank(TAlphabet c, uint64_t i) const {
    assert(c < (1llu << logsigma()));
    return wwt_.rank(std::min(i + 1, size()), c);
}

template <class t_wt_sdsl>
uint64_t wavelet_tree_sdsl<t_wt_sdsl>::select(TAlphabet c, uint64_t i) const {
    assert(i > 0 && size() > 0);
    assert(i <= rank(c, size() - 1));
    assert(c < (1llu << logsigma()));
    return wwt_.select(i, c);
}

template <class t_wt_sdsl>
TAlphabet wavelet_tree_sdsl<t_wt_sdsl>::operator[](uint64_t i) const {
    assert(i < size());
    return wwt_[i];
}

template <class t_wt_sdsl>
uint64_t wavelet_tree_sdsl<t_wt_sdsl>::next(uint64_t i, TAlphabet c) const {
    assert(i < size());
    assert(c < (1llu << logsigma()));

    return ::next(*this, i, c, MAX_ITER_WAVELET_TREE_SMALL);
}

template <class t_wt_sdsl>
uint64_t wavelet_tree_sdsl<t_wt_sdsl>::prev(uint64_t i, TAlphabet c) const {
    assert(i < size());
    assert(c < (1llu << logsigma()));

    return ::prev(*this, i, c, MAX_ITER_WAVELET_TREE_SMALL);
}

template <class t_wt_sdsl>
void wavelet_tree_sdsl<t_wt_sdsl>::serialize(std::ostream &out) const {
    wwt_.serialize(out);
}

template <class t_wt_sdsl>
bool wavelet_tree_sdsl<t_wt_sdsl>::load(std::istream &in) {
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
        std::cerr << "ERROR: Not enough memory to load wavelet_tree_sdsl" << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

template <class t_wt_sdsl>
sdsl::int_vector<> wavelet_tree_sdsl<t_wt_sdsl>::to_vector() const {
    sdsl::int_vector<> vector(wwt_.size(), 0, logsigma_);
    std::copy(wwt_.begin(), wwt_.end(), vector.begin());
    return vector;
}

template class wavelet_tree_sdsl<sdsl::wt_huff<>>;


////////////////////////////////////////////////////
// partite_vector based on bit_vectors, immutable //
////////////////////////////////////////////////////

template <class t_bv>
partite_vector<t_bv>::partite_vector(uint8_t logsigma, uint64_t size, TAlphabet c)
      : int_vector_(size, c, logsigma), bitmaps_(1 << logsigma) {
    bitmaps_.at(c) = t_bv(size, 1);
}

template <class t_bv>
partite_vector<t_bv>::partite_vector(uint8_t logsigma,
                                     sdsl::int_vector<>&& vector) {
    if (vector.width() == logsigma) {
        int_vector_ = std::move(vector);
    } else {
        int_vector_ = pack_vector(vector, logsigma);
    }

    std::vector<sdsl::bit_vector> bitmaps(1 << logsigma);

    for (size_t i = 0; i < int_vector_.size(); ++i) {

        TAlphabet c = int_vector_[i];

        assert(c < static_cast<uint64_t>(bitmaps.size()));

        if (!bitmaps[c].size())
            bitmaps[c] = sdsl::bit_vector(int_vector_.size(), 0);

        bitmaps[c][i] = 1;
    }

    for (auto &bv : bitmaps) {
        bitmaps_.emplace_back(std::move(bv));
    }
}

template <class t_bv>
uint64_t partite_vector<t_bv>::rank(TAlphabet c, uint64_t i) const {
    assert(c < bitmaps_.size());
    return bitmaps_[c].size() ? bitmaps_[c].rank1(i) : 0;
}

template <class t_bv>
uint64_t partite_vector<t_bv>::select(TAlphabet c, uint64_t i) const {
    assert(i > 0 && size() > 0);
    assert(c < bitmaps_.size());
    assert(i <= rank(c, size() - 1));
    return bitmaps_[c].select1(i);
}

template <class t_bv>
TAlphabet partite_vector<t_bv>::operator[](uint64_t i) const {
    assert(i < size());
    return int_vector_[i];
}

template <class t_bv>
uint64_t partite_vector<t_bv>::next(uint64_t i, TAlphabet c) const {
    assert(i < size());
    assert(c < (1llu << logsigma()));

    return bitmaps_[c].size() ? bitmaps_[c].next1(i) : size();
}

template <class t_bv>
uint64_t partite_vector<t_bv>::prev(uint64_t i, TAlphabet c) const {
    assert(i < size());
    assert(c < (1llu << logsigma()));

    return bitmaps_[c].size() ? bitmaps_[c].prev1(i) : size();
}

template <class t_bv>
void partite_vector<t_bv>::serialize(std::ostream &out) const {
    int_vector_.serialize(out);
    for (const auto &v : bitmaps_) {
        v.serialize(out);
    }
}

template <class t_bv>
bool partite_vector<t_bv>::load(std::istream &in) {
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
        std::cerr << "ERROR: Not enough memory to load partite_vector" << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

template <class t_bv>
void partite_vector<t_bv>::clear() {
    int_vector_.resize(0);
    bitmaps_.assign(bitmaps_.size(), t_bv());
}

template class partite_vector<bit_vector_stat>;
