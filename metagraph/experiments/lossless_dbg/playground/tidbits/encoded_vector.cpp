#include <stdint.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <algorithm> // for sort

using namespace sdsl;

int main(int argc,char** argv) {
    uint64_t size = 10000;

    // create a int_vector
    int_vector<> vec(size);

    //  fill with random numbers that use <=12 bits
    for(int_vector<32>::size_type i=0;i<vec.size();i++) {
        vec[i] = rand() % 4096;
    }
    std::cout << "vec uses " << util::get_size_in_mega_bytes(vec) << " MB." << std::endl;

    // bit compress
    util::bit_compress(vec);
    std::cout << "vec uses " << util::get_size_in_mega_bytes(vec) << " MB." << std::endl;

    // create a second int vector with specified int width and default value 5
    int_vector<17> vec_17(size,5);
    std::cout << "vec_17 uses " << util::get_size_in_mega_bytes(vec_17) << " MB." << std::endl;

    // create default vector and resize to 17 bits
    int_vector<> vec_res(size);
    std::cout << "vec_res uses " << util::get_size_in_mega_bytes(vec_res) << " MB." << std::endl;
    // data is lost here
    vec_res.set_int_width(17);
    vec_res.resize(size);
    std::cout << "vec_res uses " << util::get_size_in_mega_bytes(vec_res) << " MB." << std::endl;

    // use the stl sort function to sort an int_vector
    std::sort(vec.begin(), vec.end());

    // print out content:
    int_vector<>::iterator it;
    std::cout << "vec contains:";
    for (it=vec.begin(); it!=vec.end(); ++it) std::cout << " " << *it;

    // clear up memory used by bv
    sdsl::util::clear(vec);
    sdsl::util::clear(vec_17);
    sdsl::util::clear(vec_res);

    return EXIT_SUCCESS;
}
