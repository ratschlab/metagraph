//
// Created by Jan Studený on 2019-07-08.
//

#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/construct.hpp>

using namespace std;
using namespace sdsl;

int main(int argc,char** argv) {
    enc_vector<> enc;
    enc.load(cin);
    int_vector<> vec(enc.size());
    for(int i=0; i < enc.size();i++) {
        vec[i] = enc[i];
    }
    util::bit_compress(vec);
    vec.serialize(cout);
    return 0;
}