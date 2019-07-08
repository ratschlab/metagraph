
#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/construct.hpp>

using namespace std;
using namespace sdsl;
void t1() {
    vector<int16_t> test = {0,1,2,3,4,5,6,7};
    wt_rlmn<> wt;
    construct_im(wt,test,2);
    cout << int(wt[3]) << endl;
}
void t2() {
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
    sdsl::int_vector<8> test = {0,1,2,3,4,5,6,7};
    wt_rlmn<> wt;
    construct_im(wt,test,1);
    cout << int(wt[3]) << endl;
    t1();
    t2();
    return 0;
}