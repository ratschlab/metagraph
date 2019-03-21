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
