#ifndef __BITMAP_HPP__
#define __BITMAP_HPP__

#include <functional>
#include <memory>
#include <set>

#include <sdsl/int_vector.hpp>

#include "bitmap_builder.hpp"


/**
 * An abstract bitmap class: vector of a fixed size with elements 0 and 1.
 */
class bitmap {
  public:
    virtual ~bitmap() {}

    virtual bool operator==(const bitmap &other) const final;
    virtual bool operator!=(const bitmap &other) const final { return !(*this == other); }

    virtual bool operator[](uint64_t id) const = 0;
    virtual uint64_t get_int(uint64_t id, uint32_t width) const = 0;

    virtual uint64_t size() const = 0;
    virtual uint64_t num_set_bits() const = 0;

    virtual void add_to(sdsl::bit_vector *other) const = 0;
    template <class T>
    using VoidCall = std::function<void(T)>;
    virtual void call_ones(const VoidCall<uint64_t> &callback) const final;
    virtual void call_ones_in_range(uint64_t begin, uint64_t end,
                                    const VoidCall<uint64_t> &callback) const = 0;
};


class bitmap_dyn : public bitmap, public bitmap_builder {
  public:
    virtual ~bitmap_dyn() {}

    virtual void add_one(uint64_t id) { set(id, true); }
    virtual void set(uint64_t id, bool val) = 0;

    // indexes - positions of newly inserted elements in resulting vector
    virtual void insert_zeros(const std::vector<uint64_t> &indexes) = 0;

    virtual bitmap_dyn& operator|=(const bitmap &other);

    virtual InitializationData get_initialization_data() {
        return { size(), num_set_bits(),
                 [&](auto callback) { call_ones(callback); } };
    }
};


class bitmap_set : public bitmap_dyn {
  public:
    explicit bitmap_set(uint64_t size = 0, bool value = 0);

    bitmap_set(uint64_t size, const std::set<uint64_t> &bits);
    bitmap_set(uint64_t size, std::set<uint64_t>&& bits) noexcept;
    bitmap_set(uint64_t size, std::initializer_list<uint64_t> init);

    void set(uint64_t id, bool val) override;
    void insert_zeros(const std::vector<uint64_t> &pos) override;

    bool operator[](uint64_t id) const override { return bits_.count(id); }
    uint64_t get_int(uint64_t id, uint32_t width) const override;

    uint64_t size() const override { return size_; }
    uint64_t num_set_bits() const override { return bits_.size(); }

    void add_to(sdsl::bit_vector *other) const override;
    void call_ones_in_range(uint64_t begin, uint64_t end,
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

    void set(uint64_t id, bool val) override;
    void insert_zeros(const std::vector<uint64_t> &pos) override;
    bitmap_vector& operator|=(const bitmap &other) override;

    bool operator[](uint64_t id) const override { return bit_vector_[id]; }
    uint64_t get_int(uint64_t id, uint32_t width) const override {
        return bit_vector_.get_int(id, width);
    }

    uint64_t size() const override { return bit_vector_.size(); }
    uint64_t num_set_bits() const override { return num_set_bits_; }

    void add_to(sdsl::bit_vector *other) const override;
    void call_ones_in_range(uint64_t begin, uint64_t end,
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

    void set(uint64_t id, bool val) override;
    void insert_zeros(const std::vector<uint64_t> &pos) override;
    bitmap_adaptive& operator|=(const bitmap &other) override;

    bool operator[](uint64_t id) const override { return (*bitmap_)[id]; }
    uint64_t get_int(uint64_t id, uint32_t width) const override {
        return bitmap_->get_int(id, width);
    }

    uint64_t size() const override { return bitmap_->size(); }
    uint64_t num_set_bits() const override { return bitmap_->num_set_bits(); }

    void add_to(sdsl::bit_vector *other) const override { bitmap_->add_to(other); }
    void call_ones_in_range(uint64_t begin, uint64_t end,
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
    explicit bitmap_lazy(size_t size = 0, bool value = false);

    bitmap_lazy(std::function<bool(uint64_t)>&& callback,
                size_t size,
                size_t num_set_bits = -1) noexcept;

    bool operator[](uint64_t id) const { return is_bit_set_(id); }
    uint64_t get_int(uint64_t id, uint32_t width) const;

    uint64_t size() const { return size_; }
    uint64_t num_set_bits() const;

    void add_to(sdsl::bit_vector *other) const;
    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const;

  private:
    std::function<bool(uint64_t)> is_bit_set_;
    size_t size_;
    size_t num_set_bits_;
};

class bitmap_generator : public bitmap {
  public:
    explicit bitmap_generator(size_t size = 0, bool value = false);

    bitmap_generator(std::vector<uint64_t>&& set_bits, size_t size) noexcept
          : size_(size), num_set_bits_(set_bits.size()),
            generator_([s=std::move(set_bits)](const VoidCall<uint64_t> &callback) {
                           std::for_each(s.begin(), s.end(), callback);
                       }) {}

    bitmap_generator(std::function<void(const VoidCall<uint64_t>&)>&& generator,
                     size_t size,
                     size_t num_set_bits = -1) noexcept;

    bool operator[](uint64_t) const { throw std::runtime_error("Not implemented"); }
    uint64_t get_int(uint64_t, uint32_t) const { throw std::runtime_error("Not implemented"); }

    uint64_t size() const { return size_; }
    uint64_t num_set_bits() const { return num_set_bits_; }

    void add_to(sdsl::bit_vector *other) const;

    // This is inefficient when begin > 0 and end < size_. This will throw an
    // exception if it is called in that way.
    void call_ones_in_range(uint64_t begin, uint64_t end,
                            const VoidCall<uint64_t> &callback) const;

  private:
    size_t size_;
    size_t num_set_bits_;
    std::function<void(const VoidCall<uint64_t>&)> generator_;
};

#endif // __BITMAP_HPP__
