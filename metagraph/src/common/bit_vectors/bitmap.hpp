#ifndef __BITMAP_HPP__
#define __BITMAP_HPP__

#include <functional>
#include <memory>
#include <set>

#include <sdsl/int_vector.hpp>

template <class T>
using VoidCall = std::function<void(T)>;


sdsl::bit_vector to_sdsl(const std::vector<bool> &vector);

void call_ones(const sdsl::bit_vector &vector,
               const VoidCall<uint64_t> &callback);

void call_ones(const sdsl::bit_vector &vector,
               uint64_t begin, uint64_t end,
               const VoidCall<uint64_t> &callback);

void call_zeros(const sdsl::bit_vector &vector,
                const VoidCall<uint64_t> &callback);

void call_zeros(const sdsl::bit_vector &vector,
                uint64_t begin, uint64_t end,
                const VoidCall<uint64_t> &callback);

uint64_t inner_prod(const sdsl::bit_vector &first,
                    const sdsl::bit_vector &second);


class bitmap {
  public:
    virtual ~bitmap() {}

    virtual void set(uint64_t id, bool val) = 0;

    virtual bool operator[](uint64_t id) const = 0;
    virtual uint64_t get_int(uint64_t id, uint32_t width) const = 0;

    virtual uint64_t size() const = 0;
    virtual uint64_t num_set_bits() const = 0;

    virtual void add_to(sdsl::bit_vector *other) const;
    virtual void call_ones(const VoidCall<uint64_t> &callback) const;
    virtual void call_ones_in_range(uint64_t begin, uint64_t end,
                                    const VoidCall<uint64_t> &callback) const = 0;
};


class bitmap_dyn : public bitmap {
  public:
    virtual ~bitmap_dyn() {}

    // indexes - positions of newly inserted elements in resulting vector
    virtual void insert_zeros(const std::vector<uint64_t> &indexes) = 0;

    virtual bitmap_dyn& operator|=(const bitmap &other);
};


class bitmap_set : public bitmap_dyn {
  public:
    explicit bitmap_set(uint64_t size = 0, bool value = 0);

    bitmap_set(uint64_t size, const std::set<uint64_t> &bits);
    bitmap_set(uint64_t size, std::set<uint64_t>&& bits) noexcept;
    bitmap_set(uint64_t size, std::initializer_list<uint64_t> init);

    virtual void set(uint64_t id, bool val) override;
    virtual void insert_zeros(const std::vector<uint64_t> &pos) override;

    virtual bool operator[](uint64_t id) const override { return bits_.count(id); }
    virtual uint64_t get_int(uint64_t id, uint32_t width) const override;

    virtual uint64_t size() const override { return size_; }
    virtual uint64_t num_set_bits() const override { return bits_.size(); }

    virtual void call_ones_in_range(uint64_t begin, uint64_t end,
                                    const VoidCall<uint64_t> &callback) const override;

    const std::set<uint64_t>& data() const { return bits_; }

  private:
    uint64_t size_;
    std::set<uint64_t> bits_;
};


class bitmap_vector : public bitmap_dyn {
  public:
    explicit bitmap_vector(uint64_t size = 0, bool value = 0);

    explicit bitmap_vector(const sdsl::bit_vector &vector);
    bitmap_vector(sdsl::bit_vector&& vector) noexcept;
    bitmap_vector(std::initializer_list<bool> bitmap);

    virtual void set(uint64_t id, bool val) override;
    virtual void insert_zeros(const std::vector<uint64_t> &pos) override;
    virtual bitmap_vector& operator|=(const bitmap &other) override;

    virtual bool operator[](uint64_t id) const override { return bit_vector_[id]; }
    virtual uint64_t get_int(uint64_t id, uint32_t width) const override {
        return bit_vector_.get_int(id, width);
    }

    virtual uint64_t size() const override { return bit_vector_.size(); }
    virtual uint64_t num_set_bits() const override { return num_set_bits_; }

    virtual void add_to(sdsl::bit_vector *other) const override;
    virtual void call_ones_in_range(uint64_t begin, uint64_t end,
                                    const VoidCall<uint64_t> &callback) const override;

    const sdsl::bit_vector& data() const { return bit_vector_; }

  private:
    uint64_t num_set_bits_ = 0;
    sdsl::bit_vector bit_vector_;
};


class bitmap_adaptive : public bitmap_dyn {
  public:
    explicit bitmap_adaptive(uint64_t size = 0, bool value = 0);

    explicit bitmap_adaptive(const sdsl::bit_vector &vector);
    bitmap_adaptive(sdsl::bit_vector&& vector) noexcept;
    bitmap_adaptive(std::initializer_list<bool> bitmap);

    bitmap_adaptive(uint64_t size, const std::set<uint64_t> &bits);
    bitmap_adaptive(uint64_t size, std::set<uint64_t>&& bits) noexcept;
    bitmap_adaptive(uint64_t size, std::initializer_list<uint64_t> bits);

    virtual void set(uint64_t id, bool val) override;
    virtual void insert_zeros(const std::vector<uint64_t> &pos) override;
    virtual bitmap_adaptive& operator|=(const bitmap &other) override;

    virtual bool operator[](uint64_t id) const override { return (*bitmap_)[id]; }
    virtual uint64_t get_int(uint64_t id, uint32_t width) const override {
        return bitmap_->get_int(id, width);
    }

    virtual uint64_t size() const override { return bitmap_->size(); }
    virtual uint64_t num_set_bits() const override { return bitmap_->num_set_bits(); }

    virtual void add_to(sdsl::bit_vector *other) const override { bitmap_->add_to(other); }
    virtual void call_ones_in_range(uint64_t begin, uint64_t end,
                                    const VoidCall<uint64_t> &callback) const override {
        bitmap_->call_ones_in_range(begin, end, callback);
    }

    const bitmap_dyn& data() const {
        assert(bitmap_.get());
        return *bitmap_.get();
    }

  private:
    void check_switch();
    void to_bit_vector();
    void to_set();

    std::unique_ptr<bitmap_dyn> bitmap_;

    static const size_t kRowCutoff;
    static const size_t kMaxNumIndicesLogRatio;
    static const size_t kNumIndicesMargin;
};

class bitmap_lazy : public bitmap {
  public:
    typedef std::function<bool(uint64_t)> BoolCallback;

    explicit bitmap_lazy(size_t size = 0, bool value = false);

    bitmap_lazy(BoolCallback callback,
                size_t size = -1,
                size_t num_set_bits = -1) noexcept;

    void set(uint64_t, bool) {
        throw std::runtime_error("Not implemented.");
    }

    bool operator[](uint64_t id) const { return in_bitmap_(id); }
    uint64_t get_int(uint64_t id, uint32_t width) const;

    uint64_t size() const;
    uint64_t num_set_bits() const;

    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const;

  private:
    BoolCallback in_bitmap_;
    size_t size_;
    size_t num_set_bits_;
};

#endif // __BITMAP_HPP__
