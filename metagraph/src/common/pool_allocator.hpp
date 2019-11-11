#ifndef __POOL_ALLOCATOR_HPP__
#define __POOL_ALLOCATOR_HPP__

#include <memory>
#include <vector>


/** A faster alternative to std::allocator<T>
 *
 * The code was copied and has been modified from:
 * https://probablydance.com/2014/11/09/plalloc-a-simple-stateful-allocator-for-node-based-containers/
 */
template <typename T>
class plalloc {
  public:
    typedef T value_type;

    plalloc() = default;
    template <typename U>
    plalloc(const plalloc<U>&) {}
    plalloc(const plalloc&) {}
    plalloc& operator=(const plalloc&) { return *this; }
    plalloc(plalloc&&) = default;
    plalloc& operator=(plalloc&&) = default;

    typedef std::true_type propagate_on_container_copy_assignment;
    typedef std::true_type propagate_on_container_move_assignment;
    typedef std::true_type propagate_on_container_swap;

    bool operator==(const plalloc &other) const { return this == &other; }
    bool operator!=(const plalloc &other) const { return !(*this == other); }

    T* allocate(size_t num_to_allocate) {
        if (num_to_allocate != 1)
            return static_cast<T*>(::operator new(sizeof(T) * num_to_allocate));

        if (available.size()) {
            T *result = available.back();
            available.pop_back();
            return result;
        }

        // first allocate 8, then double whenever
        // we run out of memory
        size_t to_allocate = 8 << memory.size();
        available.reserve(to_allocate);
        std::unique_ptr<value_holder[]> allocated(new value_holder[to_allocate]);
        value_holder *first_new = allocated.get();
        memory.emplace_back(std::move(allocated));
        size_t to_return = to_allocate - 1;
        for (size_t i = 0; i < to_return; ++i) {
            available.push_back(std::addressof(first_new[i].value));
        }
        return std::addressof(first_new[to_return].value);
    }
    void deallocate(T *ptr, size_t num_to_free) {
        if (num_to_free == 1) {
            available.push_back(ptr);
        } else {
            ::operator delete(ptr);
        }
    }

    // boilerplate that shouldn't be needed, except
    // libstdc++ doesn't use allocator_traits yet
    template<typename U>
    struct rebind {
        typedef plalloc<U> other;
    };

    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

    template<typename U, typename... Args>
    void construct(U *object, Args&&... args) {
        new (object) U(std::forward<Args>(args)...);
    }
    template<typename U, typename... Args>
    void construct(const U *object, Args &&... args) = delete;

    template<typename U>
    void destroy(U *object) { object->~U(); }

  private:
    union value_holder {
        value_holder() {}
        ~value_holder() {}
        T value;
    };

    std::vector<std::unique_ptr<value_holder[]>> memory;
    std::vector<T*> available;
};

#endif // __POOL_ALLOCATOR_HPP__
