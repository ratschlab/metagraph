#ifndef __SHARED_CACHE_HPP__
#define __SHARED_CACHE_HPP__

#include <algorithm>
#include <atomic>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <new>
#include <optional>
#include <type_traits>

#include <sys/mman.h>

#include "common/utils/atomic.hpp"


namespace mtg {
namespace common {

/**
 * A fixed-size cache that threads share without locks. Each key has one slot, so a new
 * key replaces the entry that was in its slot. A sequence number per slot lets a reader
 * see either a whole entry or a miss. Reads and writes never allocate or block; a write
 * that meets another write to the same slot is dropped.
 *
 * Keys and values must be trivially copyable. A capacity of 0 disables the cache.
 */
template <typename Key, typename Value, class Hash = std::hash<Key>>
class SharedCache {
    static_assert(std::is_trivially_copyable_v<Key>);
    static_assert(std::is_trivially_copyable_v<Value>);

    struct Entry {
        Key key;
        Value value;
    };
    static constexpr size_t kWords = (sizeof(Entry) + 7) / 8;

    // Plain words, accessed through atomic_ref. The sequence number is odd while the
    // slot is written and 0 while it is empty.
    struct Slot {
        uint64_t seq;
        uint64_t words[kWords];
    };

    // a table this large is mapped from the kernel rather than calloc'ed (see the constructor)
    static constexpr size_t kMapBytes = size_t(1) << 20;

    // |mapped|: the size of a mapped table, 0 for one from calloc
    struct Free {
        size_t mapped = 0;
        void operator()(Slot *slots) const {
            if (mapped) {
                munmap(slots, mapped);
            } else {
                std::free(slots);
            }
        }
    };

  public:
    // |capacity| is rounded up to a power of two (at least 2)
    explicit SharedCache(size_t capacity) {
        if (!capacity)
            return;
        capacity = std::bit_ceil(std::max(capacity, size_t(2)));
        // Zeroed pages are only touched when a slot is first used. A large table is mapped
        // from the kernel for this: calloc under an allocator that recycles freed blocks
        // (jemalloc, glibc) zeroes a recycled block in full and keeps its pages resident after
        // free, so that each short-lived cache (one per request) would write and hold all of
        // its pages. A mapping is unmapped, its pages returned, when the cache is destroyed.
        const size_t bytes = capacity * sizeof(Slot);
        if (bytes >= kMapBytes) {
            void *table = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                               MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
            if (table == MAP_FAILED)
                throw std::bad_alloc();
            slots_ = std::unique_ptr<Slot[], Free>(static_cast<Slot*>(table), Free { bytes });
        } else {
            slots_.reset(static_cast<Slot*>(std::calloc(capacity, sizeof(Slot))));
            if (!slots_)
                throw std::bad_alloc();
        }
        shift_ = 64 - std::countr_zero(capacity);
    }

    size_t capacity() const { return slots_ ? size_t(1) << (64 - shift_) : 0; }

    std::optional<Value> find(const Key &key) const {
        if (!slots_)
            return std::nullopt;

        Slot &slot = slot_of(key);
        uint64_t seq = utils::atomic_load(slot.seq, std::memory_order_acquire);
        if (!seq || seq % 2)
            return std::nullopt;

        uint64_t words[kWords];
        for (size_t i = 0; i < kWords; ++i) {
            words[i] = utils::atomic_load(slot.words[i], std::memory_order_relaxed);
        }
        std::atomic_thread_fence(std::memory_order_acquire);
        if (utils::atomic_load(slot.seq, std::memory_order_relaxed) != seq)
            return std::nullopt;

        Entry entry;
        std::memcpy(&entry, words, sizeof(Entry));
        if (!(entry.key == key))
            return std::nullopt;

        return entry.value;
    }

    void insert(const Key &key, const Value &value) {
        if (!slots_)
            return;

        Slot &slot = slot_of(key);
        std::atomic_ref<uint64_t> seq(slot.seq);
        uint64_t old = seq.load(std::memory_order_relaxed);
        if (old % 2 || !seq.compare_exchange_strong(old, old + 1, std::memory_order_acquire))
            return;

        std::atomic_thread_fence(std::memory_order_release);
        uint64_t words[kWords] = {};
        Entry entry { key, value };
        std::memcpy(words, &entry, sizeof(Entry));
        for (size_t i = 0; i < kWords; ++i) {
            std::atomic_ref<uint64_t>(slot.words[i]).store(words[i], std::memory_order_relaxed);
        }
        seq.store(old + 2, std::memory_order_release);
    }

  private:
    Slot& slot_of(const Key &key) const {
        // Fibonacci hashing spreads consecutive ids over the slots
        return slots_[(Hash()(key) * 0x9E3779B97F4A7C15ull) >> shift_];
    }

    std::unique_ptr<Slot[], Free> slots_;
    int shift_ = 64;
};

} // namespace common
} // namespace mtg

#endif // __SHARED_CACHE_HPP__
