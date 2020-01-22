
#include <iostream>
#include <vector>
#include "wavelet_tree.hpp"

using namespace std;
int main(int argc, char **argv) {
    int number_of_wt = 50;
    int initial_size = 1908;
    vector<int> initial_content(initial_size, 6);
    vector<wavelet_tree_dyn> size_test(number_of_wt, wavelet_tree_dyn(3, initial_content));
#pragma omp parallel for num_threads(5)
    for (int i = 0; i < number_of_wt; i++) {
        for (int j = 0; j < 2'000'00; j++) {
            // int pick = rand() % number_of_wt;
            int pick = i;
            size_test[pick].insert(rand() % size_test[pick].size(), rand() % 6);
        }
    }
    for (int i = 0; i < number_of_wt; i++) {
        for (int j = 0; j < initial_size; j++) {
            cout << size_test[i].select(6, j + 1) << ",";
        }
        cout << endl;
        assert(size_test[i].rank(6, size_test[i].size() - 1) == initial_size);
    }

    return 0;
}
//
// Starting tracking the heap
// Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0001.heap (128 MB currently
// in use) Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0002.heap (257 MB
// currently in use) Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0003.heap
// (386 MB currently in use) Dumping heap profile to
// /media/studenyj/Linux/wavelet.hprof.0004.heap (515 MB currently in use) Dumping heap
// profile to /media/studenyj/Linux/wavelet.hprof.0005.heap (644 MB currently in use) Dumping
// heap profile to /media/studenyj/Linux/wavelet.hprof.0006.heap (774 MB currently in use)
// Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0007.heap (903 MB currently
// in use) Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0008.heap (1032 MB
// currently in use) Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0009.heap
// (1161 MB currently in use) Dumping heap profile to
// /media/studenyj/Linux/wavelet.hprof.0010.heap (1290 MB currently in use) Dumping heap
// profile to /media/studenyj/Linux/wavelet.hprof.0011.heap (1419 MB currently in use) Dumping
// heap profile to /media/studenyj/Linux/wavelet.hprof.0012.heap (1548 MB currently in use)
// Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0013.heap (1677 MB
// currently in use) Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0014.heap
// (1806 MB currently in use) Dumping heap profile to
// /media/studenyj/Linux/wavelet.hprof.0015.heap (1935 MB currently in use) Dumping heap
// profile to /media/studenyj/Linux/wavelet.hprof.0016.heap (2064 MB currently in use) Dumping
// heap profile to /media/studenyj/Linux/wavelet.hprof.0017.heap (2193 MB currently in use)
// Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0018.heap (2322 MB
// currently in use) Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0019.heap
// (2451 MB currently in use) Dumping heap profile to
// /media/studenyj/Linux/wavelet.hprof.0020.heap (2580 MB currently in use) Dumping heap
// profile to /media/studenyj/Linux/wavelet.hprof.0021.heap (2709 MB currently in use) Dumping
// heap profile to /media/studenyj/Linux/wavelet.hprof.0022.heap (2838 MB currently in use)
// Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0023.heap (2967 MB
// currently in use) Dumping heap profile to /media/studenyj/Linux/wavelet.hprof.0024.heap
// (3096 MB currently in use) Dumping heap profile to
// /media/studenyj/Linux/wavelet.hprof.0025.heap (3225 MB currently in use) Dumping heap
// profile to /media/studenyj/Linux/wavelet.hprof.0026.heap (Exiting, 72 bytes in use)