#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/construct.hpp>

using namespace std;
using namespace sdsl;


void int_vector_adj_width() {
    sdsl::int_vector<> test = {1,2,3,4,5,6,7,8};
    test.width(8);
    wt_rlmn<> wt;
    construct_im(wt,test,0);
    for(auto e : wt) {
        cout << int(e) << " ";
    }
    cout << endl;
}

void int_vector_compile_time_width_5() {
    sdsl::int_vector<5> test = {1,2,3,4,5,6,7,8};
    wt_rlmn<> wt;
    construct_im(wt,test,0);
    for(auto e : wt) {
        cout << int(e) << " ";
    }
    cout << endl;
}

void int_vector_auto_width() {
    sdsl::int_vector<> test = {1,2,3,4,5,6,7,8};
    wt_rlmn<> wt;
    construct_im(wt,test,0);
    for(auto e : wt) {
        cout << int(e) << " ";
    }
    cout << endl;
}



void int_vector_width_during_construction() {
    sdsl::int_vector<> test(10, 0, 5);
    test[0] = 0;
    test[1] = 1;
    test[2] = 2;
    test[3] = 3;
    test[4] = 4;
    test[5] = 5;
    test[6] = 6;
    test[7] = 7;
    test[8] = 8;
    wt_rlmn<> wt;
    construct_im(wt,test,0);
    cout << int(wt[4]) << endl;
}



void t6() {
    char test[8] = {1,2,3,4,5,6,7,8};
    wt_rlmn<> wt;
    construct_im(wt,test,1);
    cout << int(wt[4]) << endl;
}

void t7() {
    int_vector_buffer<8> test("test",
                              std::ios::out,
                              1024*1024,
                              1,
                              false);
    test.push_back(0);
    test.push_back(1);
    test.push_back(2);
    test.push_back(3);
    test.push_back(4);
    test.push_back(5);
    test.push_back(6);
    test.push_back(7);
    wt_rlmn<> wt(test,test.size());
    cout << int(wt[3]) << endl;
}



int main(int argc,char** argv) {
    int_vector_adj_width();
    int_vector_compile_time_width_5();
    int_vector_auto_width();
    //t6(); // wrong, assertion failed: (i < size()), function operator[], file ...../sdsl-lite/include/sdsl/wt_rlmn.hpp
    //t7(); // works, not using construct_im
    return 0;
}