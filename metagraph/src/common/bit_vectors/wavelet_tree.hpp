#ifndef __WAVELET_TREE_HPP__
#define __WAVELET_TREE_HPP__

#include <cstdint>
#include <mutex>
#include <atomic>

#include <sdsl/wavelet_trees.hpp>
#include <libmaus2/wavelet/DynamicWaveletTree.hpp>

#include "bit_vector.hpp"


class wavelet_tree {
  public:
    virtual ~wavelet_tree() {};

    virtual uint64_t rank(uint64_t c, uint64_t i) const = 0;
    virtual uint64_t select(uint64_t c, uint64_t i) const = 0;
    virtual uint64_t operator[](uint64_t id) const = 0;

    // get the position of the next value |var| in subvector [id, ...]
    virtual uint64_t next(uint64_t id, uint64_t val) const = 0;
    // get the position of the previous value |var| in subvector [..., id]
    // if doesn't exist, return size()
    virtual uint64_t prev(uint64_t id, uint64_t val) const = 0;

    virtual void set(uint64_t id, uint64_t val) = 0;
    virtual void insert(uint64_t id, uint64_t val) = 0;
    virtual void remove(uint64_t id) = 0;
    virtual void clear() = 0;

    virtual uint64_t size() const = 0;
    virtual uint8_t logsigma() const = 0;

    virtual bool load(std::istream &in) = 0;
    virtual void serialize(std::ostream &out) const = 0;

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

    void set(uint64_t id, uint64_t val);
    void insert(uint64_t id, uint64_t val);
    void remove(uint64_t id);

    uint64_t size() const { return n_; }
    uint8_t logsigma() const { return int_vector_.width(); }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear();

    sdsl::int_vector<> to_vector() const;

  private:
    void init_wt();

    sdsl::int_vector<> int_vector_;
    sdsl::wt_huff<> wwt_;
    std::atomic_bool requires_update_ { true };
    std::mutex mu_;
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

    void set(uint64_t id, uint64_t val);
    void insert(uint64_t id, uint64_t val);
    void remove(uint64_t id);

    uint64_t size() const { return int_vector_.size(); }
    uint8_t logsigma() const { return int_vector_.width(); }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear();

    sdsl::int_vector<> to_vector() const { return int_vector_; }
    const sdsl::int_vector<>& get_vector() const { return int_vector_; }

  private:
    sdsl::int_vector<> int_vector_;
    std::vector<bit_vector_stat> bitmaps_;
};


class wavelet_tree_dyn : public wavelet_tree {
  public:
    explicit wavelet_tree_dyn(uint8_t logsigma);

    template <class Vector>
    wavelet_tree_dyn(uint8_t logsigma, const Vector &vector);

    wavelet_tree_dyn(const wavelet_tree_dyn &other);
    wavelet_tree_dyn(wavelet_tree_dyn&& other) noexcept;
    wavelet_tree_dyn& operator=(const wavelet_tree_dyn &other);
    wavelet_tree_dyn& operator=(wavelet_tree_dyn&& other) noexcept;

    uint64_t rank(uint64_t c, uint64_t i) const;
    uint64_t select(uint64_t c, uint64_t i) const;
    uint64_t operator[](uint64_t id) const;

    uint64_t next(uint64_t id, uint64_t val) const;
    uint64_t prev(uint64_t id, uint64_t val) const;

    void set(uint64_t id, uint64_t val);
    void insert(uint64_t id, uint64_t val);
    void remove(uint64_t id);

    uint64_t size() const { return wwt_->size(); }
    uint8_t logsigma() const { return wwt_->b; }

    bool load(std::istream &in);
    void serialize(std::ostream &out) const;

    void clear();

    sdsl::int_vector<> to_vector() const;

  private:
    std::unique_ptr<libmaus2::wavelet::DynamicWaveletTree<6, 64>> wwt_;
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

    void set(uint64_t id, uint64_t val);
    void insert(uint64_t id, uint64_t val);
    void remove(uint64_t id);

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