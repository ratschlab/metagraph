#pragma once
#include <algorithm>
#include <cstddef>
#include <limits>
#include <random>
namespace mtg {
namespace common {
template <typename T>
size_t
find_next_X_block(T *A, size_t x0, size_t z, size_t y, size_t k, size_t f, size_t b1, size_t b2) {
    T min1 = std::numeric_limits<T>::max();
    T min2 = std::numeric_limits<T>::max();
    size_t m = ((z - x0 - f) / k) * k + f + x0;
    if (m <= z)
        m += k;
    size_t i = m; // find from m; the start of the block adjacent to the right of z
    size_t result = 0;
    while (i + k <= y) {
        if (i != b1 && i != b2) {
            size_t j = (i < b1 && b1 < i + k) ? m - 1 : i + k - 1;
            if (A[i] <= min1 && A[j] <= min2) {
                result = i;
                min1 = A[i];
                min2 = A[j];
            }
        }
        i += k;
    }
    return result;
}

template <typename T>
void mergeBandY(T *A, size_t z, size_t y, size_t yn) {
    while (z < y && y < yn) {
        T *Aj = std::min_element(A + z, A + y);
        if (*Aj <= A[y]) {
            std::swap(A[z], *Aj);
        } else {
            std::swap(A[z], A[y]);
            y++;
        }
        z++;
    }
    if (z < y)
        std::sort(A + z, A + yn);
}

template <typename T>
void merge(T *A, size_t beg, size_t mid, size_t end, size_t k) {
    size_t f = (mid - beg) % k;
    size_t x = (f == 0) ? mid - 2 * k : mid - k - f;
    T tmp = A[x];
    A[x] = A[beg];
    size_t z = beg;
    size_t y = mid;
    size_t buf1 = x + 1;
    size_t buf2 = mid - k;
    while (y - z > 2 * k) { // test if X-elements (excluding the buffer elements) are exhausted
        if (A[x] <= A[y] || y >= end) {
            A[z] = A[x];
            A[x] = A[buf1];
            x++;
            if ((x - beg) % k == f) { // x reached the start position of the new block
                if (z < x - k)
                    buf2 = x - k;
                x = find_next_X_block(A, beg, z, y, k, f, buf1, buf2);
            }
        } else {
            A[z] = A[y];
            A[y] = A[buf1];
            y++;
            if ((y - mid) % k == 0)
                buf2 = y - k;
        }
        z++;
        A[buf1] = A[z];
        if (z == x)
            x = buf1;
        if (z == buf2)
            buf2 = -1;
        buf1 = buf1 + 1;
        if ((buf1 - beg) % k == f)
            buf1 = buf2 == -1 ? buf1 - k : buf2;
    }
    A[z] = tmp;
    mergeBandY(A, z, y, end);
}
} // namespace common
} // namespace mtg
/*
int main() {
    std::mt19937 rng(123456);
    std::uniform_int_distribution<int> uniform_dist(1, 10);
    std::vector<int> a(100);
    size_t mid = a.size()/2;
    a[0] = uniform_dist(rng);
    for (uint32_t i = 1; i < mid; ++i) {
        a[i] = a[i - 1] + uniform_dist(rng);
    }
    a[mid] = uniform_dist(rng);
    for (uint32_t i = mid+1; i < a.size(); ++i) {
        a[i] = a[i - 1] + uniform_dist(rng);
    }
    merge(a.data(), 0, mid, a.size(), 5);
    for (const auto v : a) {
        std::cout << v << " ";
    }
}
*/
