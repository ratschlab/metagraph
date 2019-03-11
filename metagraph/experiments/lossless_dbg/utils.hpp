//
//  utils.hpp
//  lossless_dbg
//
//  Created by Jan Studený on 11/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#ifndef utils_h
#define utils_h

#include "sequence_io.hpp"


void transform_to_fasta(const string &filename,vector<string> reads) {
    ofstream myfile;
    myfile.open (filename);
    for(auto& read : reads) {
        myfile << ">" << endl;
        myfile << read << endl;
    }
    myfile.close();
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
