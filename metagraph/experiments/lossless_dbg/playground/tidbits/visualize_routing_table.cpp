//
// Created by Jan Studený on 2019-03-11.
//
#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <tclap/CmdLine.h>
#include <random>
#include <nlohmann/json.hpp>
#include <sdsl/io.hpp>
#include <sdsl/wt_rlmn.hpp>

using TCLAP::ValueArg;
using TCLAP::MultiArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::UnlabeledMultiArg;
using TCLAP::ValuesConstraint;


using namespace std;
using namespace sdsl;
using namespace nlohmann;
using namespace std::string_literals;

int main(int argc, char *argv[]) {
    // TCLAP::CmdLine cmd("Compress reads",' ', "0.1");
    // TCLAP::ValueArg<std::string> inputArg("i",
    //                                       "input",
    //                                       "Graph to use as a reference in compression",
    //                                       true,
    //                                       "",
    //                                       "string",cmd);
    // cmd.parse(argc, argv);
    // ifstream joins_file(inputArg.getValue());
    sdsl::wt_rlmn<> wt;
    wt.load(cin);
    sdsl::write_structure<JSON_FORMAT>(wt, cout);
    return 0;
}

