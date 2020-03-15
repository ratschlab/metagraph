#ifndef __WAVELET_TREE_HPP__
#define __WAVELET_TREE_HPP__

#include <cstdint>
#include <mutex>
#include <atomic>

#include <sdsl/wt_huff.hpp>
#include <dynamic.hpp>

#include "bit_vector_stat.hpp"


class wavelet_tree {
  public:
    virtual ~wavelet_tree() {};

    virtual bool operator==(const wavelet_tree &other) const final;
    virtual bool operator!=(const wavelet_tree &other) const final { return !(*this == other); }

    virtual uint64_t rank(uint64_t c, uint64_t i) const = 0;
    virtual uint64_t select(uint64_t c, uint64_t i) const = 0;
    virtual uint64_t operator[](uint64_t id) const = 0;

    // get the position of the next value |val| in subvector [id, ...]
    virtual uint64_t next(uint64_t id, uint64_t val) const = 0;
    // get the position of the previous value |val| in subvector [..., id]
    // if doesn't exist, return size()
    virtual uint64_t prev(uint64_t id, uint64_t val) const = 0;

    virtual uint64_t size() const = 0;
    virtual uint8_t logsigma() const = 0;

    virtual bool load(std::istream &in) = 0;
    virtual void serialize(std::ostream &out) const = 0;

    virtual void clear() = 0;

    // FYI: This function invalidates the current object
    template <class WaveletTree>
    WaveletTree convert_to();

    virtual sdsl::int_vector<> to_vector() const = 0;
};


class wavelet_tree_stat : public wavelet_tree {
    friend wavelet_tree;

  public:
    explicit wavelet_tree_stat(uint8_t logsigma,
                               uint64_t size = 0, uint64_t value = 0);
    template <class Vector>
    wavelet_tree_stat(uint8_t logsigma, const Vector &vector);
    wavelet_tree_stat(uint8_t logsigma, sdsl::int_vector<>&& vector);
    wavelet_tree_stat(uint8_t logsigma, sdsl::wt_huff<>&& wwt);

    wavelet_tree_stat(const wavelet_tree_stat &other);
    wavelet_tree_stat(wavelet_tree_stat&& other) noexcept;
    wavelet_tree_stat& operator=(const wavelet_tree_stat &other);
    wavelet_tree_stat& operator=(wavelet_tree_stat&& other) noexcept;

    uint64_t rank(uint64_t c, uint64_t i) const;
    uint64_t select(uint64_t c, uint64_t i) const;
    uint64_t operator[](uint64_t id) const;

    uint64_t next(uint64_t id, uint64_t val) const;
    uint64_t prev(uint64_t id, uint64_t val) const;

    uint64_t size() const { return n_; }
    uint8_t logsigma() const { return int_vector_.width(); }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear();

    sdsl::int_vector<> to_vector() const;
    const sdsl::int_vector<>& data() const { return int_vector_; }

  private:
    void init_wt() const;

    mutable sdsl::int_vector<> int_vector_;
    mutable sdsl::wt_huff<> wwt_;
    mutable std::atomic_bool requires_update_ { true };
    mutable std::mutex mu_;
    uint64_t n_;
};


// FYI: this, in fact, isn't a wavelet tree
class wavelet_tree_fast : public wavelet_tree {
    friend wavelet_tree;

  public:
    explicit wavelet_tree_fast(uint8_t logsigma,
                               uint64_t size = 0, uint64_t value = 0);
    template <class Vector>
    wavelet_tree_fast(uint8_t logsigma, const Vector &vector);

    uint64_t rank(uint64_t c, uint64_t i) const;
    uint64_t select(uint64_t c, uint64_t i) const;
    uint64_t operator[](uint64_t id) const;

    uint64_t next(uint64_t id, uint64_t val) const;
    uint64_t prev(uint64_t id, uint64_t val) const;

    uint64_t size() const { return int_vector_.size(); }
    uint8_t logsigma() const { return int_vector_.width(); }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear();

    sdsl::int_vector<> to_vector() const { return int_vector_; }
    const sdsl::int_vector<>& data() const { return int_vector_; }

  private:
    sdsl::int_vector<> int_vector_;
    std::vector<bit_vector_stat> bitmaps_;
};


class wavelet_tree_dyn : public wavelet_tree {
  public:
    explicit wavelet_tree_dyn(uint8_t logsigma);

    template <class Vector>
    wavelet_tree_dyn(uint8_t logsigma, const Vector &vector);

    uint64_t rank(uint64_t c, uint64_t i) const;
    uint64_t select(uint64_t c, uint64_t i) const;
    uint64_t operator[](uint64_t id) const;

    uint64_t next(uint64_t id, uint64_t val) const;
    uint64_t prev(uint64_t id, uint64_t val) const;

    void set(uint64_t id, uint64_t val);
    void insert(uint64_t id, uint64_t val);
    void remove(uint64_t id);

    uint64_t size() const { return dwt_.size(); }
    uint8_t logsigma() const;

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear();

    sdsl::int_vector<> to_vector() const;

  private:
    using dwt_type = dyn::wt_str;
    dwt_type dwt_;
};


class wavelet_tree_small : public wavelet_tree {
    friend wavelet_tree;

  public:
    explicit wavelet_tree_small(uint8_t logsigma) : logsigma_(logsigma) {}

    template <class Vector>
    wavelet_tree_small(uint8_t logsigma, const Vector &vector);

    wavelet_tree_small(uint8_t logsigma, const sdsl::wt_huff<> &wwt);
    wavelet_tree_small(uint8_t logsigma, sdsl::wt_huff<>&& wwt);

    uint64_t rank(uint64_t c, uint64_t i) const;
    uint64_t select(uint64_t c, uint64_t i) const;
    uint64_t operator[](uint64_t id) const;

    uint64_t next(uint64_t id, uint64_t val) const;
    uint64_t prev(uint64_t id, uint64_t val) const;

    uint64_t size() const { return wwt_.size(); }
    uint8_t logsigma() const { return logsigma_; }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear() { wwt_ = sdsl::wt_huff<>(); }

    sdsl::int_vector<> to_vector() const;

  private:
    sdsl::wt_huff<> wwt_;
    uint8_t logsigma_;
};

#endif // __WAVELET_TREE_HPP__
