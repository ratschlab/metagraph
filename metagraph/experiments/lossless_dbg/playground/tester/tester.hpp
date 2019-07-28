//
// Created by Jan Studený on 2019-07-04.
//

#ifndef METAGRAPH_TESTER_HPP
#define METAGRAPH_TESTER_HPP

//
// Created by Jan Studený on 2019-03-11.
//
#include <iostream>
#include <map>
#include <tclap/CmdLine.h>

using TCLAP::ValueArg;
using TCLAP::MultiArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::UnlabeledMultiArg;
using TCLAP::ValuesConstraint;

#include "utilities.hpp"

using namespace std;
using namespace std::string_literals;


using node_index = SequenceGraph::node_index;



int main_tester(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Compress reads",' ', "0.1");
    TCLAP::ValueArg<std::string> leftArg("l",
                                         "left_hand_side",
                                         "FASTA/Q file that should be compressed",
                                         true,
                                         "",
                                         "string");
    TCLAP::ValueArg<std::string> rightArg("r",
                                          "right_hand_side",
                                          "FASTA/Q file that should be compressed",
                                          true,
                                          "",
                                          "string");
    cmd.add(leftArg);
    cmd.add(rightArg);
    cmd.parse(argc, argv);
    auto left = leftArg.getValue();
    auto right = rightArg.getValue();
    auto left_reads_ordered = read_reads_from_fasta(left);
    auto right_reads_ordered = read_reads_from_fasta(right);
    auto left_reads = multiset<string>(all(left_reads_ordered));
    auto right_reads = multiset<string>(all(right_reads_ordered));

    if (left_reads == right_reads) {
        cout << "Reads are identical up to ordering." << endl;
        return 0;
    }
    else {
        cout << "Files differ!!!" << endl;
        cerr << "Files differ!!!" << endl;
        cout << "l: " << left_reads.size() << endl;
        cout << "r: " << right_reads.size() << endl;
        std::multiset<string> only_in_left;
        std::set_difference(all(left_reads),
                            all(right_reads),
                            std::inserter(only_in_left, only_in_left.begin()));
        std::multiset<string> only_in_right;
        std::set_difference(all(right_reads),
                            all(left_reads),
                            std::inserter(only_in_right, only_in_right.begin()));
        cout << "----- Only in left -----" << endl;
        for(auto &e : only_in_left) {
            cout << e << endl;
        }
        cout << "----- Only in right -----" << endl;
        for(auto &e : only_in_right) {
            cout << e << endl;
        }
        return -1;
    }
}


#endif //METAGRAPH_TESTER_HPP
