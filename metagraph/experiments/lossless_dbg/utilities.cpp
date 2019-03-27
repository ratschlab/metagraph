//
// Created by Jan Studený on 2019-03-21.
//

#include "utilities.hpp"

void save_string(const string &to_save, const string &filename) {
    ofstream myfile;
    myfile.open (filename);
    myfile << to_save;
    myfile.close();
}

void transform_to_fasta(const string &filename, const vector <string> &reads) {
    ofstream myfile;
    myfile.open (filename);
    for(auto& read : reads) {
        myfile << ">" << endl;
        myfile << read << endl;
    }
    myfile.close();
}

template<typename T>
d_t &d_t::operator,(const T &first) {
    std::cerr << ' ' <<  x;
    return *this;
}

void write_reads_to_fasta(const vector <string> &reads, const string &filename) {
    transform_to_fasta(filename,reads);
}

void reduce_maps(std::map<int, int> &output, std::map<int, int> &input) {
    for (auto& X : input) {
        output[X.first] += X.second;
    }
}

vector <string> read_reads_from_fasta(const string &filename) {
    vector<string> result;
    read_fasta_file_critical(
            filename,
            [&](kseq_t* read) {
                result.push_back(read->seq.s);
            });
    return result;
}

template<typename POD>
std::istream &deserialize(std::istream &is, vector<POD> &v) {
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
                  "Can only deserialize POD types with this function");

    decltype(v.size()) size;
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    v.resize(size);
    is.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(POD));
    return is;
}

template<typename POD>
std::ostream &serialize(std::ostream &os, const vector<POD> &v) {
    // this only works on built in data types (PODs)
    static_assert(std::is_trivial<POD>::value && std::is_standard_layout<POD>::value,
                  "Can only serialize POD types with this function");

    auto size = v.size();
    os.write(reinterpret_cast<char const*>(&size), sizeof(size));
    os.write(reinterpret_cast<char const*>(v.data()), v.size() * sizeof(POD));
    return os;
}
