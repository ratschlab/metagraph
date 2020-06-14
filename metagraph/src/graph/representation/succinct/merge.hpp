#pragma once
#include <algorithm>
#include <cstddef>
#include <limits>
#include <random>
#include <iostream>

namespace mtg {
namespace common {
template <typename T>
size_t find_next_X_block(T *arr,
                         size_t beg,
                         size_t z,
                         size_t y,
                         size_t k,
                         size_t f,
                         size_t buf1,
                         size_t buf2) {
    T min1 = std::numeric_limits<T>::max();
    T min2 = std::numeric_limits<T>::max();
    size_t m = (static_cast<int64_t>(z - beg - f) / static_cast<int64_t>(k)) * k + f + beg;
    if (m <= z)
        m += k;
    size_t i = m; // find from m; the start of the block adjacent to the right of z
    size_t result = 0;
    while (i + k <= y) {
        if (i != buf1 && i != buf2) {
            size_t j = (i < buf1 && buf1 < i + k) ? m - 1 : i + k - 1;
            if (arr[i] <= min1 && arr[j] <= min2) {
                result = i;
                min1 = arr[i];
                min2 = arr[j];
            }
        }
        i += k;
    }
    return result;
}

/**
 * Merges the  remaining Y-elements A[y..end) with the buffer elementsA[z..y−1]. Note that
 * A[y..end) is  already  sorted/
 */
template <typename T>
void merge_b_and_y(T *arr, size_t z, size_t y, size_t end) {
    while (z < y && y < end) {
        T *Aj = std::min_element(arr + z, arr + y);
        if (*Aj <= arr[y]) {
            std::swap(arr[z], *Aj);
        } else {
            std::swap(arr[z], arr[y]);
            y++;
        }
        z++;
    }
    if (z < y)
        std::sort(arr + z, arr + end);
}

/**
 * In place merge-sort implementation, based on "A simple algorithm for in-place merging"
 * by Jing-Chao Chen,
 * @param T the elements to be merged in-place. The two sorted arrays are beg..mid-1 and
 * mid..end-1.
 * @param beg index of the first element of the first array
 * @param mid index the first element of the second array
 * @param end total number of elements
 * @param k the merge block size
 */
template <typename T>
void merge(T *arr, size_t beg, size_t mid, size_t end, size_t k) {
    assert(beg <= mid && mid < end);
    static constexpr size_t MAX = std::numeric_limits<size_t>::max();
    size_t f = (mid - beg) % k;
    size_t x = (f == 0) ? mid - 2 * k : mid - k - f; // current element in 1st half
    if (mid - beg < 2 * k)
        x = beg;
    T tmp = arr[x];
    arr[x] = arr[beg];
    size_t z = beg; // z is the "whole" position; everything to its left is sorted
    size_t y = mid; // current element in 2nd half
    size_t buf1 = x + 1; // buf1 and buf2 are the 2 buffers for temporary data
    size_t buf2 = mid - k;
    while (y - z > 2 * k) { // test if X-elements (excluding the buffer elements) are exhausted
        if (y >= end || arr[x] <= arr[y]) {
            arr[z] = arr[x];
            arr[x] = arr[buf1];
            x++;
            if ((x - beg) % k == f) { // x reached the start position of the new block
                if (z < x - k)
                    buf2 = x - k;
                x = find_next_X_block(arr, beg, z, y, k, f, buf1, buf2);
            }
        } else {
            arr[z] = arr[y];
            arr[y] = arr[buf1];
            y++;
            if ((y - mid) % k == 0)
                buf2 = y - k;
        }
        z++;
        arr[buf1] = arr[z];
        if (z == x)
            x = buf1;
        if (z == buf2)
            buf2 = MAX;
        buf1++;
        if ((buf1 - beg) % k == f)
            buf1 = buf2 == MAX ? buf1 - k : buf2;
    }
    arr[z] = tmp;
    merge_b_and_y(arr, z, y, end);
}
} // namespace common
} // namespace mtg
/*
int main() {
    std::mt19937 rng(123456);
    std::uniform_int_distribution<int> uniform_dist(1, 10);
    for (uint32_t size = 10; size < 1000; ++size) {
        for (uint32_t k = sqrt(size); k < 2 * sqrt(size); ++k) {
            for (int mid_diff = -10; mid_diff < 10; ++mid_diff) {
                std::vector<int> a(size);
                size_t mid = std::max(std::min(1, static_cast<int>(a.size()) / 2 + mid_diff),
                                      static_cast<int>(size - 2));
                a[0] = uniform_dist(rng);
                for (uint32_t i = 1; i < mid; ++i) {
                    a[i] = a[i - 1] + uniform_dist(rng);
                }
                a[mid] = uniform_dist(rng);
                for (uint32_t i = mid + 1; i < a.size(); ++i) {
                    a[i] = a[i - 1] + uniform_dist(rng);
                }
                mtg::common::merge(a.data(), 0, mid, a.size(), k);
                bool broken = false;
                for (uint32_t i = 1; i < a.size(); ++i) {
                    if (a[i] < a[i - 1]) {
                        broken = true;
                        break;
                    }
                }
                if (broken) {
                    std::cout << a[0] << " ";
                    for (uint32_t i = 1; i < a.size(); ++i) {
                        std::cout << a[i] << " ";
                        if (a[i] < a[i - 1]) {
                            std::cout << "\nBroken at " << i << " " << size << " " << k
                                      << std::endl;
                        }
                    }
                    std::exit(0);
                }
            }
        }
    }
}
*/
