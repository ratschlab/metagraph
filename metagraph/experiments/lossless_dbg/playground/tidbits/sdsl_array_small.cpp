
#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/construct.hpp>

using namespace std;
using namespace sdsl;
void t1() {
    vector<char> test = {0,1,2,3,4,5,6,7};
    wt_rlmn<> wt;
    construct_im(wt,test,1);
    cout << int(wt[3]) << endl;
}

int main(int argc,char** argv) {
    sdsl::int_vector<8> test = {0,1,2,3,4,5,6,7};
    wt_rlmn<> wt;
    construct_im(wt,test,1);
    cout << int(wt[3]) << endl;
    t1();
    return 0;
}