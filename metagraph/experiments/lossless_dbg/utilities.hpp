//
//  utils.hpp
//  lossless_dbg
//
//  Created by Jan Studený on 11/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#ifndef utils_h
#define utils_h

#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <nlohmann/json.hpp>
#include "sequence_io.hpp"

using namespace std;

//namespace nlohmann {
//    template <typename T>
//    struct adl_serializer<std::map<int, T>> {
//    static void to_json(json& j, std::unordered_map<int, T> const& map) {
//        for (auto& [key, value] : map) {
//            // map int to string
//            j[std::to_string(key)] = value;
//        }
//    }
//};
//}

#define x first
#define y second
#define all(x) begin(x),end(x)

namespace fs = std::filesystem;

#define local_file(filename) (fs::path(__FILE__).parent_path() / (filename))

void save_string(const string &to_save,const string &filename);

//template<typename POD>
//std::ostream& serialize(std::ostream& os, std::vector<POD> const& v);
//
//template<typename POD>
//std::istream& deserialize(std::istream& is, std::vector<POD>& v);

struct d_t {
    template<typename T>
    d_t operator,(const T &first) {
        std::cerr << ' ' <<  x;
        return *this;
    }
};



#define D(args ...) { d_t(), "|", __LINE__, "|", #args, ":", args, "\n"; }

void transform_to_fasta(const string &filename,const vector<string>& reads);

void write_reads_to_fasta(const vector<string>& reads,const string &filename);

vector<string> read_reads_from_fasta(const string &filename);

// openmp reductions

void reduce_maps(\
                 std::map<int, int>& output, \
                 std::map<int, int>& input);

//#pragma omp declare reduction(map_reduction : \
//std::map<int, int> : \
//reduce_maps(omp_out, omp_in)) \
//initializer(omp_priv(omp_orig))


template<typename POD>
inline std::istream &deserialize(std::istream &is, vector<POD> &v) {
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
                  "Can only deserialize POD types with this function");

    decltype(v.size()) size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    v.resize(size);
    is.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(POD));
    return is;
}

template<typename POD>
inline std::ostream &serialize(std::ostream &os, const vector<POD> &v) {
    // this only works on built in data types (PODs)
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
                  "Can only serialize POD types with this function");

    auto size = v.size();
    os.write(reinterpret_cast<char const*>(&size), sizeof(size));
    os.write(reinterpret_cast<char const*>(v.data()), v.size() * sizeof(POD));
    return os;
}

#endif /* utils_h */
