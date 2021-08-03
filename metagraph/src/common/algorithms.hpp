#ifndef __ALGORITHMS_HPP__
#define __ALGORITHMS_HPP__

#include <vector>
#include <numeric>
#include <functional>
#include <utility>
#include <cassert>
#include <random>
#include <set>


namespace utils {

    /**
     *  This function checks whether two given strings are identical.
     */
    template <class String>
    bool seq_equal(const String &s1, const String &s2, size_t start = 0) {
        if (s1.size() != s2.size())
            return false;

        for (size_t i = start; i < s1.size(); ++i) {
            if (s1.at(i) != s2.at(i))
                return false;
        }
        return true;
    }

    /**
     *  This function checks whether string s1 is co-lexicographically
     *  greater than s2.
     */
    template <class String>
    bool colexicographically_greater(const String &s1, const String &s2) {
        size_t ss1 = s1.size();
        size_t ss2 = s2.size();
        for (size_t i = 1; i <= std::min(ss1, ss2); ++i) {
            if (s1.at(ss1 - i) != s2.at(ss2 - i))
                return (s1.at(ss1 - i) > s2.at(ss2 - i));
        }
        return ss1 > ss2;
    }

    // For each occurrence of |label| in |array|, mark segment
    // of length |segment_length| starting at that position.
    // Example: [X]***[X]****[X]***
    //    --->  [111]0[111]00[111]0
    template <class Vector>
    inline std::vector<bool> drag_and_mark_segments(const Vector &array,
                                                    typename Vector::value_type label,
                                                    size_t segment_length) {
        std::vector<bool> mask(array.size(), false);
        size_t last_occurrence
            = std::find(array.data(), array.data() + array.size(), label)
                - array.data();

        for (size_t i = last_occurrence; i < array.size(); ++i) {
            if (array[i] == label)
                last_occurrence = i;

            if (i - last_occurrence < segment_length)
                mask[i] = true;
        }

        return mask;
    }

    template <class AIt, class BIt>
    uint64_t count_intersection(AIt first_begin, AIt first_end,
                                BIt second_begin, BIt second_end) {
        assert(std::is_sorted(first_begin, first_end));
        assert(std::is_sorted(second_begin, second_end));
        assert(std::adjacent_find(first_begin, first_end) == first_end);
        assert(std::adjacent_find(second_begin, second_end) == second_end);

        uint64_t count = 0;

        while (first_begin != first_end && second_begin != second_end) {
            first_begin = std::lower_bound(first_begin, first_end, *second_begin);

            if (first_begin == first_end)
                break;

            second_begin = std::lower_bound(second_begin, second_end, *first_begin);

            if (second_begin == second_end)
                break;

            if (*first_begin == *second_begin) {
                ++count;
                ++first_begin;
                ++second_begin;
            }
        }

        return count;
    }

    // Return true if the two sorted ranges share a common element
    template <class InputIt1, class InputIt2>
    constexpr bool share_element(InputIt1 first1,
                                 InputIt1 last1,
                                 InputIt2 first2,
                                 InputIt2 last2) {
        assert(std::is_sorted(first1, last1));
        assert(std::is_sorted(first2, last2));

        while (first1 != last1 && first2 != last2) {
            if (*first1 < *first2) {
                first1 = std::lower_bound(first1, last1, *first2);
            } else if (*first2 < *first1) {
                first2 = std::lower_bound(first2, last2, *first1);
            } else {
                return true;
            }
        }

        return false;
    }

    // Intersect the sorted ranges a1 and b1 with corresponding sorted ranges of
    // sorted ranges a2 and b2 (of equal length).
    // i.e., For each shared element between a1 and b1, intersect the corresponding
    // ranges in a2 and b2.
    template <class OutType, class SetOp,
              class InputIt1, class InputIt2, class InputIt3, class InputIt4,
              class OutputIt1, class OutputIt2>
    constexpr void indexed_set_op(InputIt1 a1_begin,
                                  InputIt1 a1_end,
                                  InputIt2 a2_begin,
                                  InputIt3 b1_begin,
                                  InputIt3 b1_end,
                                  InputIt4 b2_begin,
                                  OutputIt1 out1,
                                  OutputIt2 out2) {
        SetOp set_op;

        while (a1_begin != a1_end && b1_begin != b1_end) {
            if (*a1_begin < *b1_begin) {
                auto a1_begin_next = std::lower_bound(a1_begin, a1_end, *b1_begin);
                a2_begin += a1_begin_next - a1_begin;
                a1_begin = a1_begin_next;
            } else if (*b1_begin < *a1_begin) {
                auto b1_begin_next = std::lower_bound(b1_begin, b1_end, *a1_begin);
                b2_begin += b1_begin_next - b1_begin;
                b1_begin = b1_begin_next;
            } else {
                OutType merged;
                set_op(a2_begin->begin(), a2_begin->end(),
                       b2_begin->begin(), b2_begin->end(),
                       std::back_inserter(merged));
                if (merged.size()) {
                    *out1 = *a1_begin;
                    ++out1;
                    *out2 = std::move(merged);
                    ++out2;
                }
                ++a1_begin;
                ++b1_begin;
                ++a2_begin;
                ++b2_begin;
            }
        }
    }

    // Bitmap |new_indexes| marks positions of inserted values in the final vector
    template <class Vector, class Bitmap>
    void insert(Vector *vector,
                const Bitmap &new_indexes,
                const typename Vector::value_type &value) {
        assert(vector);
        assert(new_indexes.size() == vector->size() + new_indexes.num_set_bits());

        if (!new_indexes.num_set_bits())
            return;

        vector->resize(vector->size() + new_indexes.num_set_bits());

        // need to move everything to the right as we call ones from left to right
        std::move_backward(vector->begin(),
                           vector->end() - new_indexes.num_set_bits(),
                           vector->end());

        size_t i = new_indexes.num_set_bits();
        size_t curpos = 0;

        // call ones from left to right
        new_indexes.call_ones([&](auto new_index) {
            assert(new_index >= curpos);

            // move block
            while (curpos < new_index) {
                (*vector)[curpos++] = std::move((*vector)[i++]);
            }
            // insert new value
            (*vector)[curpos++] = value;
        });
    }

    // new_indexes - positions of inserted values in the final vector
    template <class Vector>
    void insert(Vector *vector,
                const std::vector<uint64_t> &new_indexes,
                const typename Vector::value_type &value) {
        assert(vector);
        assert(std::is_sorted(new_indexes.begin(), new_indexes.end()));
        assert(!new_indexes.size()
                    || new_indexes.back() < vector->size() + new_indexes.size());

        if (!new_indexes.size())
            return;

        vector->resize(vector->size() + new_indexes.size());

        uint64_t i = vector->size() - 1;
        uint64_t shift = new_indexes.size();

        for (auto it = new_indexes.rbegin(); it != new_indexes.rend(); ++it) {
            while (i > *it) {
                assert(i >= shift && "Invalid indexes for insertion");
                (*vector)[i] = std::move((*vector)[i - shift]);
                i--;
            }
            // insert new value
            shift--;
            (*vector)[i--] = value;
        }
    }

    template <class Array, class Mask>
    void erase(Array *vector, const Mask &erase_mask) {
        assert(vector);
        assert(vector->size() == erase_mask.size());

        size_t j = 0;
        for (size_t i = 0; i < erase_mask.size(); ++i) {
            if (!erase_mask[i])
                (*vector)[j++] = (*vector)[i];
        }
        vector->resize(j);
    }

    template <typename T, class Vector = std::vector<T>>
    Vector arange(T first, size_t size) {
        static_assert(std::is_same_v<typename Vector::value_type, T>);
        Vector result(size);
        std::iota(result.begin(), result.end(), first);
        return result;
    }

    std::vector<uint64_t> sample_indexes(uint64_t universe_size,
                                         uint64_t sample_size,
                                         std::mt19937 &gen);

    template <typename T>
    T get_quantile(const std::vector<std::pair<T, uint64_t>> &count_hist, double q) {
        assert(q >= 0.0);
        assert(q <= 1.0);

        assert(std::is_sorted(count_hist.begin(), count_hist.end(),
                              [](const auto &first, const auto &second) {
                                  return first.first < second.first;
                              }));

        uint64_t sum = 0;
        for (const auto &pair : count_hist) {
            sum += pair.second;
        }

        const double threshold = q * sum;
        uint64_t partial_sum = 0;

        for (const auto &pair : count_hist) {
            partial_sum += pair.second;
            if (partial_sum >= threshold)
                return pair.first;
        }

        assert(false);
        return count_hist.back().first;
    }

    /**
     * Smoothing of an array with a sliding window.
     *
     * If |window_size| is even, the left side of the window is longer than
     * the right one by 1. For example, for window size 6: .....?...
     *                                                       ^--*-^
     * On the sides of the array, the window is truncated:     .?.......
     *                                                       XX^*-^
     */
    template <typename T>
    void smooth_vector(size_t window_size, std::vector<T> *v_p) {
        auto &v = *v_p;

        if (window_size <= 1)
            return;

        size_t left = window_size / 2;
        size_t right = window_size - left - 1;

        assert(left >= right);

        if (right + 1 >= v.size()) {
            T sum = std::accumulate(v.begin(), v.end(), (T)0);
            std::fill(v.begin(), v.end(), std::round((double)sum / v.size()));

        } else {
            std::vector<double> psum(v.size());
            std::partial_sum(v.begin(), v.end(), psum.begin());
            // Now we transform it to the following:
            // [1] [2] [3] [4] [4]
            //      X   ^---^---*
            //  X   ^---^---*---^
            //  ^---^---*---^
            //  ^---*---^
            //  *---^
            size_t i = 0;
            // truncated window
            //     ..........
            // ^---*---^
            for ( ; i <= left && i + right < v.size(); ++i) {
                v[i] = std::round(psum[i + right] / (i + right + 1));
            }
            // full window
            //     ..............
            //     ^---*---^
            for ( ; i + right < v.size(); ++i) {
                assert(i > left);
                v[i] = std::round((psum[i + right] - psum[i - left - 1]) / window_size);
            }
            // window fully covers the sequence
            //     ......
            //   ^---*---^
            for ( ; i <= left; ++i) {
                assert(i < v.size());
                v[i] = std::round(psum.back() / v.size());
            }
            // truncated window
            //     ..........
            //       ^---*---^
            for ( ; i < v.size(); ++i) {
                assert(i > left);
                v[i] = std::round((psum.back() - psum[i - left - 1]) / (left + v.size() - i));
            }
        }
    }

} // namespace utils

#endif // __ALGORITHMS_HPP__
