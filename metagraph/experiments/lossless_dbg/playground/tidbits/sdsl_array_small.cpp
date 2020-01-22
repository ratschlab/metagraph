#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/construct.hpp>

using namespace std;
using namespace sdsl;

void int_vector_compile_time_width_8() {
    sdsl::int_vector<8> test = { 1, 2, 3, 4, 5, 6, 7, 8 };
    wt_rlmn<> wt;
    construct_im(wt, test, 0);
    for (auto e : wt) {
        cout << int(e) << " ";
    }
    cout << endl;
}

void int_vector_compile_time_width_5() {
    sdsl::int_vector<5> test = { 1, 2, 3, 4, 5, 6, 7, 8 };
    wt_rlmn<> wt;
    construct_im(wt, test, 0);
    for (auto e : wt) {
        cout << int(e) << " ";
    }
    cout << endl;
}

void int_vector_auto_width() {
    sdsl::int_vector<> test = { 1, 2, 3, 4, 5, 6, 7, 8 };
    wt_rlmn<> wt;
    construct_im(wt, test, 0);
    for (auto e : wt) {
        cout << int(e) << " ";
    }
    cout << endl;
}

void int_vector_adj_width() {
    sdsl::int_vector<> test = { 1, 2, 3, 4, 5, 6, 7, 8 };
    test.width(8);
    wt_rlmn<> wt;
    construct_im(wt, test, 0);
    for (auto e : wt) {
        cout << int(e) << " ";
    }
    cout << endl;
}

int main(int argc, char **argv) {
    int_vector_compile_time_width_8(); // 1 2 3 4 5 6 7 8

    int_vector_compile_time_width_5(); // Warning: Width of int_vector<8> was specified as
                                       // 5 Length is 40 bits 65 12 82 204 65

    int_vector_auto_width(); // Warning: Width of int_vector<8> was specified as 64
                             // Length is 512 bits
                             // 1 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 4 0 0 0 0
                             // 0 0 0 5 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0 8 0 0 0 0 0 0

    int_vector_adj_width(); // 1 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 4 0 0 0 0 0
                            // 0 0 5 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 7 0 0 0 0 0 0 0 8 0 0 0 0 0 0
    return 0;
}