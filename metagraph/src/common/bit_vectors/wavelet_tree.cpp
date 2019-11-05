#include "wavelet_tree.hpp"

#include "../utils.hpp"
#include "../serialization.hpp"

#include <cassert>

// TODO: run benchmarks and optimize these parameters
const size_t MAX_ITER_WAVELET_TREE_STAT = 1000;
const size_t MAX_ITER_WAVELET_TREE_DYN = 0;
const size_t MAX_ITER_WAVELET_TREE_SMALL = 10;

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

  const wavelet_tree *bitmap_p[] = {this, &other};
  const sdsl::int_vector<> *bv_p[] = {nullptr, nullptr};

  for (int i : {0, 1}) {
    if (dynamic_cast<const wavelet_tree_stat *>(bitmap_p[i])) {
      bv_p[i] = &dynamic_cast<const wavelet_tree_stat &>(*bitmap_p[i]).data();
    } else if (dynamic_cast<const wavelet_tree_fast *>(bitmap_p[i])) {
      bv_p[i] = &dynamic_cast<const wavelet_tree_fast &>(*bitmap_p[i]).data();
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

template <class WaveletTree> WaveletTree wavelet_tree::convert_to() {
  if (dynamic_cast<WaveletTree *>(this)) {
    return WaveletTree(dynamic_cast<WaveletTree &&>(*this));
  } else if (dynamic_cast<wavelet_tree_small *>(this) &&
             typeid(WaveletTree) == typeid(wavelet_tree_stat)) {
    return WaveletTree(
        logsigma(), std::move(dynamic_cast<wavelet_tree_small *>(this)->wwt_));
  } else if (dynamic_cast<wavelet_tree_stat *>(this) &&
             !dynamic_cast<wavelet_tree_stat *>(this)->requires_update_ &&
             typeid(WaveletTree) == typeid(wavelet_tree_small)) {
    return WaveletTree(
        logsigma(), std::move(dynamic_cast<wavelet_tree_stat *>(this)->wwt_));
  } else {
    uint8_t bits_per_number = logsigma();

    sdsl::int_vector<> vector;
    if (dynamic_cast<wavelet_tree_stat *>(this)) {
      vector = std::move(dynamic_cast<wavelet_tree_stat *>(this)->int_vector_);
      vector.resize(dynamic_cast<wavelet_tree_stat *>(this)->n_);
    } else if (dynamic_cast<wavelet_tree_fast *>(this)) {
      vector = std::move(dynamic_cast<wavelet_tree_fast *>(this)->int_vector_);
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
inline uint64_t next(const Vector &v, uint64_t pos, uint64_t value,
                     size_t num_steps) {
  assert(pos < v.size());

  if (v[pos] == value)
    return pos;

  for (size_t t = 1; t < num_steps; ++t) {
    if (pos + t == v.size() || v[pos + t] == value)
      return pos + t;
  }

  uint64_t rk = v.rank(value, pos) + 1;
  return rk <= v.rank(value, v.size() - 1) ? v.select(value, rk) : v.size();
}

template <typename Vector>
inline uint64_t prev(const Vector &v, uint64_t pos, uint64_t value,
                     size_t num_steps) {
  assert(pos < v.size());

  for (size_t t = 0; t < num_steps; ++t, --pos) {
    if (v[pos] == value)
      return pos;

    if (pos == 0)
      return v.size();
  }

  uint64_t rk = v.rank(value, pos);
  return rk ? v.select(value, rk) : v.size();
}

///////////////////////////////////////////////////////////
// wavelet_tree_stat sdsl rank/select, int_vector access //
///////////////////////////////////////////////////////////

wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma, uint64_t size,
                                     uint64_t value)
    : int_vector_(size, value, logsigma), n_(size) {}

template <class Vector>
wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma, const Vector &vector)
    : int_vector_(pack_vector(vector, logsigma)), n_(int_vector_.size()) {}

template wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma,
                                              const sdsl::int_vector<> &vector);
template wavelet_tree_stat::wavelet_tree_stat(
    uint8_t logsigma, const std::vector<uint8_t> &vector);
template wavelet_tree_stat::wavelet_tree_stat(
    uint8_t logsigma, const std::vector<uint64_t> &vector);
template wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma,
                                              const std::vector<int> &vector);

wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma,
                                     sdsl::int_vector<> &&vector) {
  if (vector.width() == logsigma) {
    int_vector_ = std::move(vector);
  } else {
    int_vector_ = pack_vector(vector, logsigma);
  }
  n_ = int_vector_.size();
}

wavelet_tree_stat::wavelet_tree_stat(uint8_t logsigma, sdsl::wt_huff<> &&wwt)
    : int_vector_(pack_vector(wwt, logsigma)), wwt_(std::move(wwt)),
      requires_update_(false), n_(wwt_.size()) {}

wavelet_tree_stat::wavelet_tree_stat(const wavelet_tree_stat &other) {
  *this = other;
}

wavelet_tree_stat::wavelet_tree_stat(wavelet_tree_stat &&other) noexcept {
  *this = std::move(other);
}

wavelet_tree_stat &wavelet_tree_stat::
operator=(const wavelet_tree_stat &other) {
  int_vector_ = other.int_vector_;
  n_ = other.n_;
  if (!other.requires_update_) {
    requires_update_ = false;
    wwt_ = other.wwt_;
  }
  return *this;
}

wavelet_tree_stat &wavelet_tree_stat::
operator=(wavelet_tree_stat &&other) noexcept {
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
    std::cerr << "ERROR: Not enough memory to load wavelet_tree_stat."
              << std::endl;
    return false;
  } catch (...) {
    return false;
  }
}

void wavelet_tree_stat::serialize(std::ostream &out) const {
  if (requires_update_)
    init_wt();

  int_vector_.serialize(out);
  wwt_.serialize(out);
}

void wavelet_tree_stat::set(uint64_t id, uint64_t val) {
  if (int_vector_[id] == val)
    return;

  if (!requires_update_) {
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
    std::copy_backward(int_vector_.begin() + id, int_vector_.begin() + n_ - 1,
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
    std::copy(int_vector_.begin() + id + 1, int_vector_.begin() + n_,
              int_vector_.begin() + id);
  n_--;
}

uint64_t wavelet_tree_stat::rank(uint64_t c, uint64_t i) const {
  if (requires_update_)
    init_wt();

  return wwt_.rank(std::min(i + 1, size()), c);
}

uint64_t wavelet_tree_stat::select(uint64_t c, uint64_t i) const {
  assert(i > 0 && size() > 0);

  if (requires_update_)
    init_wt();

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

void wavelet_tree_stat::init_wt() const {
  std::unique_lock<std::mutex> lock(mu_);

  if (!requires_update_)
    return;

  int_vector_.resize(n_);
  wwt_ = decltype(wwt_)(int_vector_);

  requires_update_ = false;
}

////////////////////////////////////////////////
// wavelet_tree_dyn xxsds/DYNAMIC rank/select //
////////////////////////////////////////////////

wavelet_tree_dyn::wavelet_tree_dyn(uint8_t logsigma) : dwt_(1ull << logsigma) {}

template <class Vector>
wavelet_tree_dyn::wavelet_tree_dyn(uint8_t logsigma, const Vector &vector)
    : dwt_(1ull << logsigma) {
  dwt_.push_many(1ull << logsigma, vector);
}

template wavelet_tree_dyn::wavelet_tree_dyn(uint8_t,
                                            const sdsl::int_vector<> &);
template wavelet_tree_dyn::wavelet_tree_dyn(uint8_t,
                                            const std::vector<uint8_t> &);
template wavelet_tree_dyn::wavelet_tree_dyn(uint8_t, const std::vector<int> &);
template wavelet_tree_dyn::wavelet_tree_dyn(uint8_t,
                                            const std::vector<uint64_t> &);

uint64_t wavelet_tree_dyn::rank(uint64_t c, uint64_t i) const {
  return size() > 0 ? dwt_.rank(std::min(i + 1, size()), c) : 0;
}

uint64_t wavelet_tree_dyn::select(uint64_t c, uint64_t i) const {
  assert(i > 0 && size() > 0 && i <= rank(c, size() - 1));
  return dwt_.select(i - 1, c);
}

uint64_t wavelet_tree_dyn::operator[](uint64_t id) const {
  assert(id < size());
  return dwt_.at(id);
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
  dwt_.insert(id, val);
}

void wavelet_tree_dyn::remove(uint64_t id) {
  assert(id < size());
  dwt_.remove(id);
}

uint8_t wavelet_tree_dyn::logsigma() const {
  return utils::code_length(const_cast<dwt_type &>(dwt_).alphabet_size()) - 1;
}

void wavelet_tree_dyn::serialize(std::ostream &out) const {
  dwt_.serialize(out);
}

bool wavelet_tree_dyn::load(std::istream &in) {
  if (!in.good())
    return false;

  // TODO: catch reading errors
  dwt_.load(in);
  return in.good();
}

void wavelet_tree_dyn::clear() { dwt_ = decltype(dwt_)(dwt_.alphabet_size()); }

sdsl::int_vector<> wavelet_tree_dyn::to_vector() const {
  sdsl::int_vector<> vector(dwt_.size(), 0, logsigma());
  for (size_t i = 0; i < vector.size(); ++i) {
    vector[i] = dwt_.at(i);
  }
  return vector;
}

////////////////////////////////////////////////////
// wavelet_tree_small sdsl rank/select, wt access //
////////////////////////////////////////////////////

template <class Vector>
wavelet_tree_small::wavelet_tree_small(uint8_t logsigma, const Vector &vector)
    : wwt_(pack_vector(vector, logsigma)), logsigma_(logsigma) {}

template wavelet_tree_small::wavelet_tree_small(
    uint8_t logsigma, const sdsl::int_vector<> &vector);
template wavelet_tree_small::wavelet_tree_small(
    uint8_t logsigma, const std::vector<uint8_t> &vector);
template wavelet_tree_small::wavelet_tree_small(
    uint8_t logsigma, const std::vector<uint64_t> &vector);
template wavelet_tree_small::wavelet_tree_small(uint8_t logsigma,
                                                const std::vector<int> &vector);

wavelet_tree_small::wavelet_tree_small(uint8_t logsigma,
                                       const sdsl::wt_huff<> &wwt)
    : wwt_(wwt), logsigma_(logsigma) {}

wavelet_tree_small::wavelet_tree_small(uint8_t logsigma, sdsl::wt_huff<> &&wwt)
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
  assert(false);
  throw std::runtime_error("Not supported");
}

void wavelet_tree_small::insert(uint64_t, uint64_t) {
  assert(false);
  throw std::runtime_error("Not supported");
}

void wavelet_tree_small::remove(uint64_t) {
  assert(false);
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
    std::cerr << "ERROR: Not enough memory to load wavelet_tree_small"
              << std::endl;
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

wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma, uint64_t size,
                                     uint64_t value)
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

template wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma,
                                              const sdsl::int_vector<> &vector);
template wavelet_tree_fast::wavelet_tree_fast(
    uint8_t logsigma, const std::vector<uint8_t> &vector);
template wavelet_tree_fast::wavelet_tree_fast(
    uint8_t logsigma, const std::vector<uint64_t> &vector);
template wavelet_tree_fast::wavelet_tree_fast(uint8_t logsigma,
                                              const std::vector<int> &vector);

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

  std::copy_backward(int_vector_.begin() + id, int_vector_.end() - 1,
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

  std::copy(int_vector_.begin() + id + 1, int_vector_.end(),
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
    std::cerr << "ERROR: Not enough memory to load wavelet_tree_fast"
              << std::endl;
    return false;
  } catch (...) {
    return false;
  }
}

void wavelet_tree_fast::clear() {
  int_vector_.resize(0);
  bitmaps_.assign(bitmaps_.size(), bit_vector_stat());
}
