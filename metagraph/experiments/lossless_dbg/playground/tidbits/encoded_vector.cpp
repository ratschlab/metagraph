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
        //vec[i] = rand() % 4096;
        long val = ran_expo(0.01);
        vec[i] = val;
        mx = max(val,mx);
    }
    cout << "MX" << mx << endl;
    enc_vector<coder::elias_delta,128999, 12> enc(vec);
    print_size(enc,"12");
    enc_vector<coder::elias_delta,128, 12> enc2(vec);
    print_size(enc2,"12");
    enc_vector<coder::elias_delta,128, 8> enc3(vec);
    print_size(enc2,"8");
    enc_vector<coder::fibonacci,128, 12> enc4(vec);
    print_size(enc4,"fib");



    // bit compress
    util::bit_compress(vec);
    std::cout << "vec uses " << size_in_mega_bytes(vec) << " MB." << std::endl;

    // create a second int vector with specified int width and default value 5
    int_vector<17> vec_17(size,5);
    std::cout << "vec_17 uses " << size_in_mega_bytes(vec_17) << " MB." << std::endl;

    // create default vector and resize to 17 bits
    int_vector<> vec_res(size);
    std::cout << "vec_res uses " << size_in_mega_bytes(vec_res) << " MB." << std::endl;
    // data is lost here
    vec_res.width(17);
    vec_res.resize(size);
    std::cout << "vec_res uses " << size_in_mega_bytes(vec_res) << " MB." << std::endl;

    // use the stl sort function to sort an int_vector
    std::sort(vec.begin(), vec.end());

    // print out content:
    int_vector<>::iterator it;
    std::cout << "vec contains:";
    for (it=vec.begin(); it!=vec.end(); ++it) std::cout << " " << *it;

    // clear up memory used by bv
    clear(vec);
    clear(vec_17);
    clear(vec_res);

    return EXIT_SUCCESS;
}
