
#include <iostream>
#include <vector>
#include "wavelet_tree.hpp"
#include "unix_tools.hpp"
#include "vector_with_rank_support.hpp"

int x;
void populate_wt(int size) {
    Timer t;
    wavelet_tree_dyn wt(3);
    for (int i = 0; i < size; i++) {
        int pos = rand() % (wt.size() + 1);
        int val = rand() % 5;
        wt.insert(pos, val);
        auto xy = wt.rank(val, pos);
        x += wt.select(val, xy - 1);
    }
    cout << "D" << size << ":" << t.elapsed() / size << " sec/e" << endl;
}

void populate_wt_ns(int size) {
    Timer t;
    wavelet_tree_dyn wt(3);
    for (int i = 0; i < size; i++) {
        int pos = rand() % (wt.size() + 1);
        int val = rand() % 5;
        wt.insert(pos, val);
        x += wt.rank(val, pos);
    }
    cout << "R" << size << ":" << t.elapsed() / size << " sec/e" << endl;
}

void populate_vec(int size) {
    Timer t;
    VectorWithRankSupport<char> wt;
    for (int i = 0; i < size; i++) {
        int pos = rand() % (wt.size() + 1);
        int val = rand() % 5;
        wt.insert(pos, val);
        x += wt.rank(val, pos);
    }
    cout << "V" << size << ":" << t.elapsed() / size << " sec/e" << endl;
}

using namespace std;
int main(int argc, char **argv) {
    for (int size : { 100, 1'000, 10'000, 100'000, 1'000'000, 10'000'000 }) {
        populate_wt(size);
        populate_wt_ns(size);
        populate_vec(size);
    }

    return 0;
}
