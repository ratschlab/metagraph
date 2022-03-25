#ifndef __WAVELET_TREE_HPP__
#define __WAVELET_TREE_HPP__

#include <cstdint>
#include <cassert>

#include <sdsl/wt_huff.hpp>
#include <dynamic.hpp>

#include "bit_vector_sdsl.hpp"


class wavelet_tree {
  public:
    typedef uint64_t TAlphabet;

    virtual ~wavelet_tree() {}

    virtual bool operator==(const wavelet_tree &other) const final;
    virtual bool operator!=(const wavelet_tree &other) const final { return !(*this == other); }

    virtual uint64_t rank(TAlphabet c, uint64_t i) const = 0;
    virtual uint64_t select(TAlphabet c, uint64_t i) const = 0;
    virtual TAlphabet operator[](uint64_t i) const = 0;
    virtual std::pair<uint64_t, TAlphabet> inverse_select(uint64_t i) const = 0;

    // get the position of the next value |c| in subvector [i, ...]
    virtual uint64_t next(uint64_t i, TAlphabet c) const = 0;
    // get the position of the previous value |c| in subvector [..., i]
    // if doesn't exist, return size()
    virtual uint64_t prev(uint64_t i, TAlphabet c) const = 0;

    virtual uint64_t size() const = 0;
    virtual uint8_t logsigma() const = 0;
    virtual uint64_t count(TAlphabet c) const = 0;

    virtual bool load(std::istream &in) = 0;
    virtual void serialize(std::ostream &out) const = 0;

    virtual void clear() = 0;

    // FYI: This function invalidates the current object
    template <class WaveletTree>
    WaveletTree convert_to();

    virtual sdsl::int_vector<> to_vector() const = 0;
};


class wavelet_tree_sdsl_augmented : public wavelet_tree {
  public:
    virtual ~wavelet_tree_sdsl_augmented() {}
    virtual const sdsl::int_vector<>& data() const = 0;
};


template <class t_wt_sdsl = sdsl::wt_huff<>>
class wavelet_tree_sdsl_fast : public wavelet_tree_sdsl_augmented {
    friend wavelet_tree;

  public:
    typedef t_wt_sdsl wt_type;

    explicit wavelet_tree_sdsl_fast(uint8_t logsigma)
      : wavelet_tree_sdsl_fast(logsigma, sdsl::int_vector<>()) {}

    template <class Vector>
    wavelet_tree_sdsl_fast(uint8_t logsigma, const Vector &vector)
      : wavelet_tree_sdsl_fast(logsigma, pack_vector(vector, logsigma)) {}
    wavelet_tree_sdsl_fast(uint8_t logsigma, sdsl::int_vector<>&& vector);

    wavelet_tree_sdsl_fast(uint8_t logsigma, const t_wt_sdsl &wwt)
      : wavelet_tree_sdsl_fast(logsigma, t_wt_sdsl(wwt)) {}
    wavelet_tree_sdsl_fast(uint8_t logsigma, t_wt_sdsl&& wwt);

    uint64_t rank(TAlphabet c, uint64_t i) const;
    uint64_t select(TAlphabet c, uint64_t i) const;
    TAlphabet operator[](uint64_t i) const;
    std::pair<uint64_t, TAlphabet> inverse_select(uint64_t i) const;

    uint64_t next(uint64_t i, TAlphabet c) const;
    uint64_t prev(uint64_t i, TAlphabet c) const;

    uint64_t size() const { return int_vector_.size(); }
    uint8_t logsigma() const { return int_vector_.width(); }
    uint64_t count(TAlphabet c) const { return count_[c]; }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear();

    sdsl::int_vector<> to_vector() const { return int_vector_; }
    const sdsl::int_vector<>& data() const { return int_vector_; }

  private:
    sdsl::int_vector<> int_vector_;
    t_wt_sdsl wwt_;
    std::vector<uint64_t> count_;
};


// A straightforward and fast implementation of the wavelet tree interface
template <class t_bv = bit_vector_stat>
class partite_vector : public wavelet_tree_sdsl_augmented {
    friend wavelet_tree;

  public:
    typedef t_bv bv_type;

    explicit partite_vector(uint8_t logsigma)
      : partite_vector(logsigma, sdsl::int_vector<>()) {}

    template <class Vector>
    partite_vector(uint8_t logsigma, const Vector &vector)
      : partite_vector(logsigma, pack_vector(vector, logsigma)) {}

    partite_vector(uint8_t logsigma, sdsl::int_vector<>&& vector);

    uint64_t rank(TAlphabet c, uint64_t i) const;
    uint64_t select(TAlphabet c, uint64_t i) const;
    TAlphabet operator[](uint64_t i) const;
    std::pair<uint64_t, TAlphabet> inverse_select(uint64_t i) const;

    uint64_t next(uint64_t i, TAlphabet c) const;
    uint64_t prev(uint64_t i, TAlphabet c) const;

    uint64_t size() const { return int_vector_.size(); }
    uint8_t logsigma() const { return int_vector_.width(); }
    uint64_t count(TAlphabet c) const { return bitmaps_[c].num_set_bits(); }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear();

    sdsl::int_vector<> to_vector() const { return int_vector_; }
    const sdsl::int_vector<>& data() const { return int_vector_; }

  private:
    sdsl::int_vector<> int_vector_;
    std::vector<t_bv> bitmaps_;
};


class wavelet_tree_dyn : public wavelet_tree {
  public:
    explicit wavelet_tree_dyn(uint8_t logsigma);

    template <class Vector>
    wavelet_tree_dyn(uint8_t logsigma, const Vector &vector);

    uint64_t rank(TAlphabet c, uint64_t i) const;
    uint64_t select(TAlphabet c, uint64_t i) const;
    TAlphabet operator[](uint64_t i) const;
    std::pair<uint64_t, TAlphabet> inverse_select(uint64_t i) const;

    uint64_t next(uint64_t i, TAlphabet c) const;
    uint64_t prev(uint64_t i, TAlphabet c) const;

    void set(uint64_t i, TAlphabet c);
    void insert(uint64_t i, TAlphabet c);
    void remove(uint64_t i);

    uint64_t size() const { return dwt_.size(); }
    uint8_t logsigma() const;
    uint64_t count(TAlphabet c) const { return rank(c, size()); }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear();

    sdsl::int_vector<> to_vector() const;

  private:
    dyn::wt_str dwt_;
};


template <class t_wt_sdsl = sdsl::wt_huff<>>
class wavelet_tree_sdsl : public wavelet_tree {
    friend wavelet_tree;

  public:
    typedef t_wt_sdsl wt_type;

    explicit wavelet_tree_sdsl(uint8_t logsigma)
      : wavelet_tree_sdsl(logsigma, t_wt_sdsl()) {}

    template <class Vector>
    wavelet_tree_sdsl(uint8_t logsigma, const Vector &vector)
      : wavelet_tree_sdsl(logsigma, t_wt_sdsl(vector)) {}

    wavelet_tree_sdsl(uint8_t logsigma, const t_wt_sdsl &wwt)
      : wavelet_tree_sdsl(logsigma, t_wt_sdsl(wwt)) {}

    wavelet_tree_sdsl(uint8_t logsigma, t_wt_sdsl&& wwt);

    uint64_t rank(TAlphabet c, uint64_t i) const;
    uint64_t select(TAlphabet c, uint64_t i) const;
    TAlphabet operator[](uint64_t i) const;
    std::pair<uint64_t, TAlphabet> inverse_select(uint64_t i) const;

    uint64_t next(uint64_t i, TAlphabet c) const;
    uint64_t prev(uint64_t i, TAlphabet c) const;

    uint64_t size() const { return wwt_.size(); }
    uint8_t logsigma() const { return logsigma_; }
    uint64_t count(TAlphabet c) const { return count_[c]; }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear() { wwt_ = t_wt_sdsl(); count_.assign(count_.size(), 0); }

    sdsl::int_vector<> to_vector() const;

  private:
    t_wt_sdsl wwt_;
    uint8_t logsigma_;
    std::vector<uint64_t> count_;
};


typedef wavelet_tree_sdsl<> wavelet_tree_stat;

typedef wavelet_tree_sdsl<sdsl::wt_huff<sdsl::rrr_vector<63>>> wavelet_tree_small;

typedef partite_vector<> wavelet_tree_fast;

#endif // __WAVELET_TREE_HPP__
