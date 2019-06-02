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
#include "kmer_extractor.hpp"
#include "sequence_io.hpp"

using namespace std;
using ll = long long;

#define x first
#define y second
#define all(x) begin(x),end(x)

namespace fs = std::filesystem;

#define local_file(filename) (fs::path(__FILE__).parent_path() / (filename))

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

inline void save_string(const string &to_save, const string &filename) {
    ofstream myfile;
    myfile.open (filename);
    myfile << to_save;
    myfile.close();
}

inline void transform_to_fasta(const string &filename, const vector <string> &reads) {
    ofstream myfile;
    myfile.open (filename);
    for(auto& read : reads) {
        myfile << ">" << endl;
        myfile << read << endl;
    }
    myfile.close();
}


inline void write_reads_to_fasta(const vector <string> &reads, const string &filename) {
    transform_to_fasta(filename,reads);
}

inline void reduce_maps(std::map<int, int> &output, std::map<int, int> &input) {
    for (auto& X : input) {
        output[X.first] += X.second;
    }
}

inline int8_t encode(char c) {
    if (c == '#') return 6;//alphabet_decoder.alph_size;
    if (c == '$') return 0;
    return KmerExtractor::encode(c);
}
inline char decode(int8_t c) {
    if (c == 6) return '#';
    if (c == 0) return '$';
    return KmerExtractor::encode(c);
}

inline string& clamp_alphabet(string& text,const string& alphabet="$ACGTN",char replacement='N') {
    for(auto& c : text) {
        if (alphabet.find(c) == string::npos) {
            c = replacement;
        }
    }
    return text;
}

inline vector <string> read_reads_from_fasta(const string &filename) {
    vector<string> result;
    read_fasta_file_critical(
            filename,
            [&](kseq_t* read) {
                string read_seq = read->seq.s;
                clamp_alphabet(read_seq);
                result.push_back(read_seq);
            });
    return result;
}

template <typename... Args>
inline void doPrint(std::ostream& out, Args&&... args)
{
    ((out << ',' << std::forward<Args>(args)), ...);
}
template <typename Test,typename Reference>
class IdentityComparator : public Reference {
public:
    using Reference::Reference;
    template<typename ...Args>
    int64_t rank(Args... args) const {
        auto target = Reference::rank(args...);
        auto value = t.rank(args...);
        if (target != value) {
            doPrint(cout,args...);
        }
        assert(target==value);
        return target;
    }

    template<typename ...Args>
    int64_t size(Args... args) const {
        auto target = Reference::size(args...);
        auto value = t.size(args...);
        if (target != value) {
            doPrint(cout,args...);
        }
        assert(target==value);
        return target;
    }

    template<typename ...Args>
    int64_t select(Args... args) const {
        auto target = Reference::select(args...);
        auto value = t.select(args...);
        if (target != value) {
            doPrint(cout,args...);
        }
        assert(target==value);
        return target;
    }

    template<typename ...Args>
    int64_t operator[](Args... args) const {
        auto target = Reference::operator[](args...);
        auto value = t.operator[](args...);
        if (target != value) {
            doPrint(cout,args...);
        }
        assert(target==value);
        return target;
    }

    template<typename ...Args>
    void insert(Args... args) {
        Reference::insert(args...);
        t.insert(args...);
    }

    Test t;
};

#endif /* utils_h */
