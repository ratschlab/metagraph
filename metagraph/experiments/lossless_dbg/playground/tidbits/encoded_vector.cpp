#include <stdint.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/enc_vector.hpp>
#include <sdsl/util.hpp>
#include <algorithm> // for sort

using namespace sdsl;
using namespace sdsl::util;
using namespace std;

template <typename T>
void print_size(T& vec,string name="") {
    std::cout << name << " uses " << size_in_mega_bytes(vec) << " MB." << std::endl;
}
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double ran_expo(double lambda){
    double u;
    u = rand() / (RAND_MAX + 1.0);
    return -log(1- u) / lambda;
}

int main(int argc,char** argv) {
    uint64_t size = 1'000'000;

    // create a int_vector
    int_vector<> vec(size);

    //  fill with random numbers that use <=12 bits
    long mx = 0;
    for(int_vector<32>::size_type i=0;i<vec.size();i++) {
        long val = rand() % 10;
        //long val = ran_expo(0.01);
        vec[i] = val;
        mx = max(val,mx);
    }
    cout << "MX" << mx << endl;
    enc_vector<> enc0(vec);
    print_size(enc0,"<>");
    enc_vector<coder::elias_delta,128999, 12> enc(vec);
    print_size(enc,"12");
    enc_vector<coder::elias_delta,128, 12> enc2(vec);
    print_size(enc2,"12");
    enc_vector<coder::elias_delta,128, 8> enc3(vec);
    print_size(enc2,"8");
    enc_vector<coder::fibonacci,128, 12> enc4(vec);
    print_size(enc4,"fib");
    print_size(vec,"orig");
    util::bit_compress(vec);
    print_size(vec,"comp");
//Output:
//MX255
//<> uses 5.19851 MB.
//12 uses 5.19095 MB.
//12 uses 5.19851 MB.
//8 uses 5.19851 MB.
//fib uses 6.14244 MB.
//orig uses 7.6294 MB.
//comp uses 0.953683 MB.
    return EXIT_SUCCESS;
}
