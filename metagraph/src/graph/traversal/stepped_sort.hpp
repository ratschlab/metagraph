#ifndef __STEPPED_SORT_HPP__
#define __STEPPED_SORT_HPP__

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <type_traits>
#include <vector>


namespace mtg {
namespace graph {
namespace traversal {

/**
 * The sorts of a /resolve's loops after the work, taken in steps so that a reading of the
 * request's stop and deadline can fall inside a sort of millions of elements
 * (SPEC-labeled-traversal-core.md §4.5; lane T8-RL: before, each was one piece the deadline was
 * read around, up to about 0.5 s at the server's caps and more on a loaded host).
 *
 * |step(n)| is told of n more elements handled, once they reach kSortStepElements since it was
 * last told: an element of a partition's pass (by copies, told every kSortStepElements of them;
 * by Hoare's scheme, at its swaps), an element of a heap sort, and each element of a range of at
 * most kSortPieceElements once std::sort sorted it. |step| may throw, and the sort then ends
 * where it is (the range's elements permuted, none lost): its caller abandons what it sorts.
 *
 * The order is std::sort's (std::partial_sort's), by the same strict weak order |comp|.
 * Elements that |comp| leaves equivalent may come in another order than std::sort would leave
 * them in (an order std::sort does not specify either): a caller whose answer must not depend
 * on it sorts by a total order, or reads nothing of the order among equivalent elements.
 */
constexpr size_t kSortStepElements = 4096;
constexpr size_t kSortPieceElements = 16384;

namespace stepped_sort_detail {

// counts the elements handled and tells |step| of them every kSortStepElements
template <class Step>
class Counter {
  public:
    explicit Counter(const Step &step) : step_(step) {}
    void operator()() {
        if (++count_ == kSortStepElements) {
            count_ = 0;
            step_(kSortStepElements);
        }
    }
    // n more elements handled, told once the count passes kSortStepElements
    void add_counted(size_t n) {
        count_ += n;
        if (count_ >= kSortStepElements) {
            const size_t told = count_;
            count_ = 0;
            step_(told);
        }
    }
    // what is counted and not yet told
    void flush() {
        if (count_) {
            const size_t n = count_;
            count_ = 0;
            step_(n);
        }
    }

  private:
    const Step &step_;
    size_t count_ = 0;
};

// sifts the element at |root| down the max-heap of the |size| elements from |first|
template <class It, class Comp>
void sift_down(It first, ptrdiff_t root, ptrdiff_t size, Comp &comp) {
    while (true) {
        ptrdiff_t child = 2 * root + 1;
        if (child >= size)
            return;
        if (child + 1 < size && comp(first[child], first[child + 1]))
            ++child;
        if (!comp(first[root], first[child]))
            return;
        std::iter_swap(first + root, first + child);
        root = child;
    }
}

// [first, last) as a max-heap by |comp|, an element counted each sift
template <class It, class Comp, class Count>
void make_heap(It first, It last, Comp &comp, Count &count) {
    const ptrdiff_t n = last - first;
    for (ptrdiff_t root = n / 2; root-- > 0; ) {
        sift_down(first, root, n, comp);
        count();
    }
}

// the max-heap [first, last) sorted ascending, an element counted each sift
template <class It, class Comp, class Count>
void sort_heap(It first, It last, Comp &comp, Count &count) {
    for (ptrdiff_t end = last - first; end > 1; ) {
        --end;
        std::iter_swap(first, first + end);
        sift_down(first, 0, end, comp);
        count();
    }
}

// the median of the first, middle and last element moved to the front: *first is the pivot, and
// the last element is not less than it (so that the partition's forward scan stops before |last|)
template <class It, class Comp>
void median_to_front(It first, It last, Comp &comp) {
    It a = first, b = first + (last - first) / 2, c = last - 1;
    if (comp(*b, *a))
        std::iter_swap(a, b);
    if (comp(*c, *b)) {
        std::iter_swap(b, c);
        if (comp(*b, *a))
            std::iter_swap(a, b);
    }
    // *a <= *b <= *c: the median to the front
    std::iter_swap(a, b);
}

// [first, last) partitioned around *first (Hoare's scheme, the pivot left in place until the
// end): returns cut with every element of [first, cut) not greater than *cut, the pivot, and
// every element of (cut, last) not less than it. The elements the scans pass counted at each swap
template <class It, class Comp, class Count>
It partition_at_front(It first, It last, Comp &comp, Count &count) {
    It i = first;
    It j = last;
    while (true) {
        const It i0 = i;
        const It j0 = j;
        // (stops at last - 1 at the latest: median_to_front)
        while (comp(*++i, *first)) {}
        // (stops at first at the latest: the pivot is not less than itself)
        while (comp(*first, *--j)) {}
        count.add_counted((i - i0) + (j0 - j));
        if (!(i < j))
            break;
        std::iter_swap(i, j);
    }
    std::iter_swap(first, j);
    return j;
}

// |a|, |b|, |c| put in order
template <class It, class Comp>
void sort3(It a, It b, It c, Comp &comp) {
    if (comp(*b, *a))
        std::iter_swap(a, b);
    if (comp(*c, *b)) {
        std::iter_swap(b, c);
        if (comp(*b, *a))
            std::iter_swap(a, b);
    }
}

// the ninther of [first, last) (at least 8 elements) moved to the front: the median of the
// medians of three triples at its front, middle and back
template <class It, class Comp>
void ninther_to_front(It first, It last, Comp &comp) {
    const ptrdiff_t h = (last - first) / 2;
    sort3(first, first + h, last - 1, comp);
    sort3(first + 1, first + (h - 1), last - 2, comp);
    sort3(first + 2, first + (h + 1), last - 3, comp);
    sort3(first + (h - 1), first + h, first + (h + 1), comp);
    std::iter_swap(first, first + h);
}

// [first, last) of trivially copyable elements partitioned around (a copy of) *first, every
// element moved without a branch on the comparison (Lomuto's scheme): returns cut with [first,
// cut) less than the pivot, *cut the pivot, and (cut, last) not less than it. |tell|(n) told of
// the elements passed, kSortStepElements at a time
template <class It, class Comp, class Tell>
It partition_by_copies(It first, It last, Comp &comp, Tell &tell, bool equal) {
    using T = typename std::iterator_traits<It>::value_type;
    const T pivot = *first;
    It write = first + 1;
    for (It read = first + 1; read < last; ) {
        const It end = read + std::min<ptrdiff_t>(last - read, kSortStepElements);
        const ptrdiff_t n = end - read;
        for ( ; read < end; ++read) {
            const T x = *read;
            // (|equal|: not greater, for a pivot that nothing is less than: its equals)
            const bool front = equal ? !comp(pivot, x) : comp(x, pivot);
            *read = *write;
            *write = x;
            write += front;
        }
        tell(n);
    }
    if (equal)
        return write;
    std::iter_swap(first, write - 1);
    return write - 1;
}

// sort_in_steps with |max_depth| partitions of a range before it is heap-sorted (the tests reach
// the heap sort with a small one)
template <class It, class Comp, class Step>
void sort_ranges(It first, It last, Comp &comp, const Step &step, int max_depth) {
    if (last - first < 2)
        return;
    constexpr bool kByCopies
            = std::is_trivially_copyable<typename std::iterator_traits<It>::value_type>::value;
    Counter<Step> count(step);
    struct Range {
        It first;
        It last;
        int depth;
    };
    std::vector<Range> pending { { first, last, max_depth } };
    while (!pending.empty()) {
        Range r = pending.back();
        pending.pop_back();
        while (r.last - r.first > static_cast<ptrdiff_t>(kSortPieceElements)) {
            if (r.depth-- == 0) {
                make_heap(r.first, r.last, comp, count);
                sort_heap(r.first, r.last, comp, count);
                r.first = r.last;
                break;
            }
            // the range's sides: [r.first, low) and [high, r.last)
            It low;
            It high;
            if constexpr (kByCopies) {
                ninther_to_front(r.first, r.last, comp);
                auto tell = [&](ptrdiff_t n) { count.add_counted(n); };
                low = partition_by_copies(r.first, r.last, comp, tell, false);
                high = low + 1;
                // the pivot the least of the range: its equals moved to the front, sorted
                if (low == r.first)
                    high = partition_by_copies(r.first, r.last, comp, tell, true);
                const ptrdiff_t n = r.last - r.first;
                const ptrdiff_t ls = low - r.first;
                const ptrdiff_t hs = r.last - high;
                if (ls < n / 8 && ls >= 32) {
                    std::iter_swap(r.first, r.first + ls / 4);
                    std::iter_swap(low - 1, low - ls / 4);
                }
                if (hs < n / 8 && hs >= 32) {
                    std::iter_swap(high, high + hs / 4);
                    std::iter_swap(r.last - 1, r.last - hs / 4);
                }
            } else {
                median_to_front(r.first, r.last, comp);
                low = partition_at_front(r.first, r.last, comp, count);
                high = low + 1;
            }
            // the smaller side first, the larger on the stack: at most log2 n ranges there
            if (low - r.first < r.last - high) {
                pending.push_back({ high, r.last, r.depth });
                r.last = low;
            } else {
                pending.push_back({ r.first, low, r.depth });
                r.first = high;
            }
        }
        if (r.last - r.first > 1) {
            std::sort(r.first, r.last, comp);
            count.add_counted(r.last - r.first);
        }
    }
    count.flush();
}

} // namespace stepped_sort_detail

/**
 * std::sort(first, last, comp), in steps: a range of at most kSortPieceElements by std::sort; a
 * larger one partitioned, the smaller side sorted first (the larger kept on a stack of at most
 * log2 n ranges); a range partitioned more than 2 log2 n times heap-sorted (std::sort's own guard
 * against a quadratic sort). Elements trivially copyable (the runs, starts and ends of a
 * max_support round, a discovery's ranking) are partitioned without branches (Lomuto's scheme)
 * around the ninther of the range, the elements equal to a pivot that is the range's least moved
 * to its front in one more pass, a few elements of a side much smaller than an eighth of the
 * range swapped (pdqsort's guards against many equal elements and patterns): about as fast as
 * std::sort at 10 million elements (lane T8-RL, bench/sv); others (the grouping's candidates)
 * by Hoare's scheme around the median of three.
 */
template <class It, class Comp, class Step>
void sort_in_steps(It first, It last, Comp comp, const Step &step) {
    int max_depth = 0;
    for (ptrdiff_t m = last - first; m > 1; m >>= 1) {
        max_depth += 2;
    }
    stepped_sort_detail::sort_ranges(first, last, comp, step, max_depth);
}

/**
 * std::partial_sort(first, middle, last, comp), in steps: the |middle - first| least elements of
 * [first, last) sorted at [first, middle), by a heap of them that each later element is offered
 * to (an element counted each one offered and each one sifted); the order of [middle, last) is
 * unspecified, as std::partial_sort's.
 */
template <class It, class Comp, class Step>
void partial_sort_in_steps(It first, It middle, It last, Comp comp, const Step &step) {
    using namespace stepped_sort_detail;
    if (first == middle)
        return;
    Counter<Step> count(step);
    const ptrdiff_t m = middle - first;
    make_heap(first, middle, comp, count);
    for (It i = middle; i != last; ++i) {
        count();
        if (comp(*i, *first)) {
            std::iter_swap(i, first);
            sift_down(first, 0, m, comp);
        }
    }
    sort_heap(first, middle, comp, count);
    count.flush();
}

/**
 * The elements of |v| copied, a step each, into a buffer of |capacity| elements (at least
 * v.size()), which then replaces v's: a vector's growth (reserve) in steps. On a throwing step
 * |v| is unchanged.
 */
template <class T, class Step>
void reserve_in_steps(std::vector<T> &v, size_t capacity, const Step &step) {
    if (capacity <= v.capacity())
        return;
    std::vector<T> grown;
    grown.reserve(capacity);
    for (size_t begin = 0; begin < v.size(); begin += kSortStepElements) {
        const size_t end = std::min(v.size(), begin + kSortStepElements);
        grown.insert(grown.end(), v.begin() + begin, v.begin() + end);
        step(end - begin);
    }
    v.swap(grown);
}

/**
 * std::unique on |v| (adjacent equal elements, by ==, kept once) and the vector cut to what
 * stays, a step each element.
 */
template <class T, class Step>
void unique_in_steps(std::vector<T> &v, const Step &step) {
    size_t kept = 0;
    for (size_t begin = 0; begin < v.size(); begin += kSortStepElements) {
        const size_t end = std::min(v.size(), begin + kSortStepElements);
        for (size_t i = begin; i < end; ++i) {
            if (!kept || !(v[kept - 1] == v[i]))
                v[kept++] = v[i];
        }
        step(end - begin);
    }
    v.resize(kept);
}

/**
 * |v| resized to |n| value-initialized elements (n not less than v.size()), its capacity
 * reserved first and the elements added kSortStepElements at a time, a step each.
 */
template <class T, class Step>
void assign_zeros_in_steps(std::vector<T> &v, size_t n, const Step &step) {
    v.clear();
    v.reserve(n);
    while (v.size() < n) {
        const size_t add = std::min(n - v.size(), kSortStepElements);
        v.resize(v.size() + add);
        step(add);
    }
}

} // namespace traversal
} // namespace graph
} // namespace mtg

#endif // __STEPPED_SORT_HPP__
