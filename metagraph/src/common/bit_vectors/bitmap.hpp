#ifndef __BITMAP_HPP__
#define __BITMAP_HPP__

#include <mutex>
#include <atomic>
#include <functional>
#include <memory>
#include <set>

#include <sdsl/int_vector.hpp>


class bitmap {
  public:
    virtual ~bitmap() {};

    virtual void set(uint64_t id, bool val) = 0;
    virtual bool operator[](uint64_t id) const = 0;
    virtual uint64_t get_int(uint64_t id, uint32_t width) const = 0;
    virtual uint64_t size() const = 0;
    virtual uint64_t num_set_bits() const = 0;

    virtual void call_ones(const std::function<void(uint64_t)> &callback) const = 0;
};


class bitmap_dyn : public bitmap {
  public:
    // indexes - positions of newly inserted elements in resulting vector
    virtual void insert_zeros(const std::vector<uint64_t> &indexes) = 0;
};


class bitmap_set : public bitmap_dyn {
  public:
    explicit bitmap_set(uint64_t size = 0, bool value = 0);

    bitmap_set(uint64_t size, const std::set<uint64_t> &bits);
    bitmap_set(uint64_t size, std::set<uint64_t>&& bits) noexcept;
    bitmap_set(uint64_t size, std::initializer_list<uint64_t> init);

    void set(uint64_t id, bool val);
    bool operator[](uint64_t id) const;
    uint64_t get_int(uint64_t id, uint32_t width) const;
    void insert_zeros(const std::vector<uint64_t> &pos);

    inline uint64_t size() const { return size_; }
    inline uint64_t num_set_bits() const { return bits_.size(); }

    void call_ones(const std::function<void(uint64_t)> &callback) const;

    inline const std::set<uint64_t>& data() const { return bits_; }
    inline std::set<uint64_t>& data() { return bits_; }

  private:
    size_t size_;
    std::set<uint64_t> bits_;
};


class bitmap_vector : public bitmap_dyn {
  public:
    explicit bitmap_vector(uint64_t size = 0, bool value = 0);

    explicit bitmap_vector(const sdsl::bit_vector &vector);
    bitmap_vector(sdsl::bit_vector&& vector) noexcept;
    bitmap_vector(std::initializer_list<bool> bitmap);

    void set(uint64_t id, bool val);
    inline bool operator[](uint64_t id) const { return bit_vector_[id]; }

    inline uint64_t get_int(uint64_t id, uint32_t width) const {
        return bit_vector_.get_int(id, width);
    }
    void insert_zeros(const std::vector<uint64_t> &pos);

    inline uint64_t size() const { return bit_vector_.size(); }
    inline uint64_t num_set_bits() const { return num_set_bits_; }

    void call_ones(const std::function<void(uint64_t)> &callback) const;

    inline const sdsl::bit_vector& data() const { return bit_vector_; }
    inline sdsl::bit_vector& data() { return bit_vector_; }

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

    void set(uint64_t id, bool val);
    inline bool operator[](uint64_t id) const { return (*bitmap_)[id]; }

    inline uint64_t get_int(uint64_t id, uint32_t width) const {
        return bitmap_->get_int(id, width);
    }

    void insert_zeros(const std::vector<uint64_t> &pos);

    inline uint64_t size() const { return bitmap_->size(); }
    inline uint64_t num_set_bits() const { return bitmap_->num_set_bits(); }

    inline void call_ones(const std::function<void(uint64_t)> &callback) const {
        bitmap_->call_ones(callback);
    }

    inline const bitmap_dyn& data() const {
        assert(bitmap_.get());
        return *bitmap_.get();
    }

    inline bitmap_dyn& data() {
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


void call_ones(const sdsl::bit_vector &vector,
               const std::function<void(uint64_t)> &callback);

void call_zeros(const sdsl::bit_vector &vector,
                const std::function<void(uint64_t)> &callback);

uint64_t inner_prod(const sdsl::bit_vector &first,
                    const sdsl::bit_vector &second);

#endif // __BITMAP_HPP__
