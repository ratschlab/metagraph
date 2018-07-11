#ifndef __WAVELET_TREE_HPP__
#define __WAVELET_TREE_HPP__

#include <cstdint>

#include <sdsl/wavelet_trees.hpp>
#include <libmaus2/wavelet/DynamicWaveletTree.hpp>


class wavelet_tree {
  public:
    virtual ~wavelet_tree() {};

    virtual uint64_t rank(uint64_t c, uint64_t i) const = 0;
    virtual uint64_t select(uint64_t c, uint64_t i) const = 0;
    virtual uint64_t size() const = 0;
    virtual void set(uint64_t id, uint64_t val) = 0;
    virtual uint64_t operator[](uint64_t id) const = 0;
    virtual void insert(uint64_t id, uint64_t val) = 0;
    virtual void remove(uint64_t id) = 0;
    virtual bool deserialise(std::istream &in) = 0;
    virtual void serialise(std::ostream &out) const = 0;
    virtual std::vector<uint64_t> to_vector() const = 0;

    // this invalidates the current object
    template <class WaveletTree>
    WaveletTree convert_to(uint64_t logsigma);
};

class wavelet_tree_small;

class wavelet_tree_stat : public wavelet_tree {
    friend class wavelet_tree_small;
    friend class wavelet_tree;

  public:
    typedef sdsl::wt_huff<> wt_type;

    explicit wavelet_tree_stat(uint64_t logsigma,
                               uint64_t size = 0,
                               uint64_t value = 0);

    template <class Vector>
    wavelet_tree_stat(uint64_t logsigma, const Vector &vector);
    wavelet_tree_stat(const wavelet_tree_stat &other);

    wavelet_tree_stat(uint64_t logsigma, wt_type&& wwt);
    wavelet_tree_stat(wavelet_tree_stat&& other);
    wavelet_tree_stat(sdsl::int_vector<>&& vector);

    wavelet_tree_stat& operator=(sdsl::int_vector<>&& vector);

    uint64_t size() const { return n_; }

    bool deserialise(std::istream &in);
    void serialise(std::ostream &out) const;

    void insert(uint64_t id, uint64_t val);
    void remove(uint64_t id);

    uint64_t rank(uint64_t c, uint64_t i) const;
    uint64_t select(uint64_t c, uint64_t i) const;

    uint64_t operator[](uint64_t id) const;
    std::vector<uint64_t> to_vector() const;

    void set(uint64_t id, uint64_t val);

  private:
    sdsl::int_vector<> int_vector_;
    wt_type wwt_;
    bool requires_update_ = true;
    uint64_t n_;

    void init_wt();
};


class wavelet_tree_small : public wavelet_tree {
    friend class wavelet_tree_stat;
    friend class wavelet_tree;

  public:
    typedef sdsl::wt_huff<> wt_type;

    explicit wavelet_tree_small(uint64_t logsigma,
                                uint64_t size = 0,
                                uint64_t value = 0);

    template <class Vector>
    wavelet_tree_small(uint64_t logsigma, const Vector &vector);

    wavelet_tree_small(const wavelet_tree_small &other);

    wavelet_tree_small(wavelet_tree_small&& other);
    wavelet_tree_small(wt_type&& wwt);
    wavelet_tree_small(sdsl::int_vector<>&& vector);

    wavelet_tree_small& operator=(sdsl::int_vector<>&& vector);

    uint64_t size() const { return wwt_.size(); }

    bool deserialise(std::istream &in);
    void serialise(std::ostream &out) const;

    void insert(uint64_t id, uint64_t val);
    void remove(uint64_t id);

    uint64_t rank(uint64_t c, uint64_t i) const;
    uint64_t select(uint64_t c, uint64_t i) const;

    uint64_t operator[](uint64_t id) const;
    std::vector<uint64_t> to_vector() const;

    void set(uint64_t id, uint64_t val);

  private:
    wt_type wwt_;
};


class wavelet_tree_dyn : public wavelet_tree {
  public:
    typedef libmaus2::wavelet::DynamicWaveletTree<6, 64> wt_type;

    explicit wavelet_tree_dyn(uint64_t logsigma) : wwt_(logsigma) {}

    template <class Vector>
    wavelet_tree_dyn(uint64_t logsigma, const Vector &W_stat);
    wavelet_tree_dyn(const wavelet_tree_dyn &other);

    wavelet_tree_dyn(wavelet_tree_dyn&& other);

    uint64_t size() const { return wwt_.size(); }

    bool deserialise(std::istream &in);
    void serialise(std::ostream &out) const;

    void insert(uint64_t id, uint64_t val);
    void remove(uint64_t id);

    uint64_t rank(uint64_t c, uint64_t i) const;
    uint64_t select(uint64_t c, uint64_t i) const;

    uint64_t operator[](uint64_t id) const;
    std::vector<uint64_t> to_vector() const;

    void set(uint64_t id, uint64_t val);

  private:
    wt_type wwt_;

    bool get_bit_raw(uint64_t id) const;
};


#endif // __WAVELET_TREE_HPP__
