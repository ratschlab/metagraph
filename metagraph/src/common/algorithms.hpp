#ifndef __ALGORITHMS_HPP__
#define __ALGORITHMS_HPP__

#include <algorithm>
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
            = std::find(array.begin(), array.end(), label) - array.begin();

        for (size_t i = last_occurrence; i < array.size(); ++i) {
            if (array[i] == label)
                last_occurrence = i;

            if (i - last_occurrence < segment_length)
                mask[i] = true;
        }

        return mask;
    }

    template <class It1, class It2>
    constexpr uint64_t count_intersection(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end) {
        assert(std::is_sorted(a_begin, a_end));
        assert(std::is_sorted(b_begin, b_end));
        assert(std::adjacent_find(a_begin, a_end) == a_end);
        assert(std::adjacent_find(b_begin, b_end) == b_end);

        uint64_t count = 0;

        while (a_begin != a_end && b_begin != b_end) {
            if (*a_begin < *b_begin) {
                ++a_begin;
            } else {
                if (*a_begin == *b_begin) {
                    ++count;
                    ++a_begin;
                }
                ++b_begin;
            }
        }

        return count;
    }

    // Return true if the two sorted ranges share a common element
    template <class It1, class It2>
    constexpr bool share_element(It1 a_begin, It1 a_end, It2 b_begin, It2 b_end) {
        assert(std::is_sorted(a_begin, a_end));
        assert(std::is_sorted(b_begin, b_end));

        while (a_begin != a_end && b_begin != b_end) {
            if (*a_begin < *b_begin) {
                ++a_begin;
            } else {
                if (*a_begin == *b_begin)
                    return true;
                ++b_begin;
            }
        }

        return false;
    }

    // Intersect sorted ranges index1 and index2 with their corresponding values
    // stored in value1 and value2, respectively.
    // For each element shared between index1 and index2, invoke callback
    // for that element and its corresponding values.
    // For each element in index1 not in index2, invoke callback_diff1 for that
    // element and its corresponding values.
    // For each element in index2 not in index1, invoke callback_diff2 for that
    // element and its corresponding values.
    template <class InIt1, class InIt2, class InIt3, class InIt4,
              class Callback, class CallbackDiff1, class CallbackDiff2>
    constexpr void match_indexed_values(InIt1 index1_begin, InIt1 index1_end,
                                        InIt2 value1_begin,
                                        InIt3 index2_begin, InIt3 index2_end,
                                        InIt4 value2_begin,
                                        const Callback &callback,
                                        const CallbackDiff1 &callback_diff1,
                                        const CallbackDiff2 &callback_diff2) {
        while (index1_begin != index1_end || index2_begin != index2_end) {
            if (index1_begin == index1_end) {
                callback_diff2(*index2_begin, *value2_begin);
                ++index2_begin;
                ++value2_begin;
            } else if (index2_begin == index2_end || *index1_begin < *index2_begin) {
                callback_diff1(*index1_begin, *value1_begin);
                ++index1_begin;
                ++value1_begin;
            } else {
                if (*index1_begin == *index2_begin) {
                    callback(*index1_begin, *value1_begin, *value2_begin);
                    ++index1_begin;
                    ++value1_begin;
                } else {
                    callback_diff2(*index2_begin, *value2_begin);
                }

                ++index2_begin;
                ++value2_begin;
            }
        }
    }

    // Output the intersection of the sorted ranges a and b to intersection_out,
    // and the elements from a not in b to diff_out (i.e., the set difference of a and b)
    template <class AIt, class BIt, class OutIt, class OutIt2>
    constexpr void set_intersection_difference(AIt a_begin, AIt a_end,
                                               BIt b_begin, BIt b_end,
                                               OutIt intersection_out, OutIt2 diff_out) {
        while (a_begin != a_end) {
            if (b_begin == b_end || *a_begin < *b_begin) {
                *diff_out = *a_begin;
                ++diff_out;
                ++a_begin;
            } else if (*a_begin > *b_begin) {
                ++b_begin;
            } else {
                *intersection_out = *a_begin;
                ++intersection_out;
                ++a_begin;
                ++b_begin;
            }
        }
    }

    // Intersect sorted ranges index1 and index2 with their corresponding values
    // stored in value1 and value2, respectively.
    // For each element shared between index1 and index2, invoke the callback
    // for that element and its corresponding values.
    template <class InIt1, class InIt2, class InIt3, class InIt4, class Callback>
    constexpr void match_indexed_values(InIt1 index1_begin, InIt1 index1_end,
                                        InIt2 value1_begin,
                                        InIt3 index2_begin, InIt3 index2_end,
                                        InIt4 value2_begin,
                                        const Callback &callback) {
        while (index1_begin != index1_end && index2_begin != index2_end) {
            if (*index1_begin < *index2_begin) {
                ++index1_begin;
                ++value1_begin;
            } else {
                if (*index1_begin == *index2_begin) {
                    callback(*index1_begin, *value1_begin, *value2_begin);
                    ++index1_begin;
                    ++value1_begin;
                }
                ++index2_begin;
                ++value2_begin;
            }
        }
    }

    template <typename It1, typename It2, typename Out>
    constexpr void set_intersection(It1 a_begin, It1 a_end,
                                    It2 b_begin, It2 b_end,
                                    Out out,
                                    int64_t delta = 0) {
        while (a_begin != a_end && b_begin != b_end) {
            if (*a_begin + delta < *b_begin) {
                ++a_begin;
            } else if (*a_begin + delta > *b_begin) {
                ++b_begin;
            } else {
                *out = *a_begin;
                ++a_begin;
                ++b_begin;
                ++out;
            }
        }
    }

    template <typename It1, typename It2, typename Out>
    constexpr void set_difference(It1 a_begin, It1 a_end,
                                  It2 b_begin, It2 b_end,
                                  Out out,
                                  int64_t delta = 0) {
        while (a_begin != a_end) {
            if (b_begin == b_end || *a_begin + delta < *b_begin) {
                *out = *a_begin;
                ++a_begin;
                ++out;
            } else if (*a_begin + delta > *b_begin) {
                ++b_begin;
            } else {
                ++a_begin;
                ++b_begin;
            }
        }
    }

    template <typename It1, typename It2, typename Out>
    constexpr void set_union(It1 a_begin, It1 a_end,
                             It2 b_begin, It2 b_end,
                             Out out,
                             int64_t delta = 0) {
        while (a_begin != a_end || b_begin != b_end) {
            if (b_begin == b_end) {
                *out = *a_begin;
                ++a_begin;
                ++out;
            } else if (a_begin == a_end || *a_begin + delta > *b_begin) {
                *out = *b_begin - delta;
                ++b_begin;
                ++out;
            } else {
                if (*a_begin + delta == *b_begin)
                    ++b_begin;

                *out = *a_begin;
                ++a_begin;
                ++out;
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
