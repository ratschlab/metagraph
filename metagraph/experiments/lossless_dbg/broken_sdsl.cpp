//
// Created by Jan Studený on 2019-04-02.
//



#include <iostream>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/construct.hpp>

using namespace std;
using namespace sdsl;

int main() {
    wt_rlmn<> wt;
    construct_im(wt,"ABCD",1);

    return 0;
}