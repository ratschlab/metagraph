
#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/construct.hpp>

using namespace std;
using namespace sdsl;
void t1() {
    vector<char> test = {1,2,3,4,5,6,7,8,0};
    wt_rlmn<> wt;
    construct_im(wt,test,1);
    cout << int(wt[3]) << endl;
}

void t2() {
    sdsl::int_vector<8> test = {1,2,3,4,5,6,7,8};
    wt_rlmn<> wt;
    construct_im(wt,test,1);
    cout << int(wt[3]) << endl;
}

void t3() {
    sdsl::int_vector<> test = {1,2,3,4,5,6,7,8};
    test.width(8);
    wt_rlmn<> wt;
    construct_im(wt,test,0);
    cout << int(wt[3]) << endl;
}

void t4() {
    // only this works
    sdsl::int_vector<> test = {1,2,3,4,5,6,7,8};
    wt_rlmn<> wt;
    construct_im(wt,test,0);
    cout << int(wt[3]) << endl;
}

void t5() {
    int_vector_buffer<> ivb("test",
                            std::ios::out, // we do not want to open an existing file, but create a new one
                            1024*1024,     // use a buffer of about 1MB
                            1,            // use 64bit for each integer
                            false);
    ivb.push_back(0);
    ivb.push_back(1);
    ivb.push_back(2);
    ivb.push_back(250);
    ivb.push_back(4);
    ivb.push_back(5);
    ivb.push_back(6);
    ivb.push_back(7);
    wt_rlmn<> wt(ivb,ivb.size());
    cout << int(wt[3]) << endl;
}

int main(int argc,char** argv) {
    sdsl::int_vector<> test(10, 0, 3);
    test[3] = 5;
    test[4] = 4;
    wt_rlmn<> wt;
    construct_im(wt,test,0);
    cout << int(wt[3]) << endl;
    cout << int(wt[4]) << endl;
    t1();
    t2();
    t3();
    t4();
    return 0;
}