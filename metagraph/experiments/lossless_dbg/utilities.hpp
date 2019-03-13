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

void save_string(const string &to_save,const string &filename) {
    ofstream myfile;
    myfile.open (filename);
    myfile << to_save;
    myfile.close();
}

struct d_t {
    template<typename T> d_t & operator,(const T & x) {
        std::cerr << ' ' <<  x;
        return *this;
    }
};

#define D(args ...) { d_t, "|", __LINE__, "|", #args, ":", args, "\n"; }

void transform_to_fasta(const string &filename,const vector<string>& reads) {
    ofstream myfile;
    myfile.open (filename);
    for(auto& read : reads) {
        myfile << ">" << endl;
        myfile << read << endl;
    }
    myfile.close();
}

void write_reads_to_fasta(const vector<string>& reads,const string &filename) {
    transform_to_fasta(filename,reads);
}

vector<string> read_reads_from_fasta(const string &filename) {
    vector<string> result;
    read_fasta_file_critical(
                             filename,
                             [&](kseq_t* read) {
                                 result.push_back(read->seq.s);
                             });
    return result;
}

// openmp reductions

void reduce_maps(\
                 std::map<int, int>& output, \
                 std::map<int, int>& input)
{
    for (auto& X : input) {
        output[X.first] += X.second;
    }
}

#pragma omp declare reduction(map_reduction : \
std::map<int, int> : \
reduce_maps(omp_out, omp_in)) \
initializer(omp_priv(omp_orig))

#endif /* utils_h */
