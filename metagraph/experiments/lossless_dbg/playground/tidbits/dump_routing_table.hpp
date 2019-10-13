//
// Created by Jan Studený on 2019-07-06.
//

#ifndef __DUMP_ROUTING_TABLE_HPP__
#define __DUMP_ROUTING_TABLE_HPP__
#include <iostream>
#include <sdsl/wt_rlmn.hpp>
#include "utilities.hpp"
int main_dumprt(int argc, char *argv[]) {
    sdsl::wt_rlmn<> wt;
    wt.load(cin);
    for(int64_t i=0;i<wt.size();i++) {
        cout << decode(wt[i]);
    }
    return 0;
}

#endif // __DUMP_ROUTING_TABLE_HPP__
