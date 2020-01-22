//
// Created by Jan Studený on 2019-07-08.
//


#include <iostream>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/construct.hpp>

#include "configuration.hpp"
using namespace sdsl;
using namespace std;
int main(int argc, char *argv[]) {
    WaveletTreeRLMN wt;
    wt.load(cin);
    int_vector<> dec(wt.size());
    for (uint64_t i = 0; i < wt.size(); i++) {
        dec[i] = wt[i];
    }
    using t_bitvector = hyb_vector<>;
    sdsl::wt_rlmn<t_bitvector, t_bitvector::rank_1_type, t_bitvector::select_1_type,
                  wt_huff<t_bitvector, t_bitvector::rank_1_type, t_bitvector::select_1_type,
                          t_bitvector::select_0_type, byte_tree<>>>
            wt_b;
    construct_im(wt_b, dec, 0);
    sdsl::write_structure<HTML_FORMAT>(wt_b, cout);
    // wt_b.serialize(cout);
    return 0;
}