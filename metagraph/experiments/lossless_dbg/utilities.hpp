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



#define PRINT_VAR(args ...) { d_t(), "|", __LINE__, "|", #args, ":", args, "\n"; }

void transform_to_fasta(const string &filename,const vector<string>& reads);

void write_reads_to_fasta(const vector<string>& reads,const string &filename);

vector<string> read_reads_from_fasta(const string &filename);

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

inline string decode_from_input(const vector<string>& input,const pair<string,int64_t>& path_id,int64_t kmer_length) {
    int64_t current_relative_id = 0;
    for(auto& sequence : input) {
        if (sequence.substr(0,kmer_length) == path_id.first) {
            if (current_relative_id == path_id.second) {
                return sequence;
            }
            else {
                current_relative_id++;
            }
        }
    }
    return "";
}

#endif /* utils_h */
