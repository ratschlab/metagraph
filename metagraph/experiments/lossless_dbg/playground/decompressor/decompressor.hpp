//
// Created by Jan Studený on 2019-03-11.
//
#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <tclap/CmdLine.h>
#include <filesystem>

using TCLAP::ValueArg;
using TCLAP::MultiArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::UnlabeledMultiArg;
using TCLAP::ValuesConstraint;

using json = nlohmann::json;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#pragma clang diagnostic ignored "-Wcomma"

#include "samplers.hpp"
#include "utilities.hpp"

#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"



using namespace std;
using namespace std::string_literals;


using node_index = SequenceGraph::node_index;


int main_decompressor(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Decompress reads",' ', "0.1");
    TCLAP::ValueArg<std::string> inputArg("i",
                                         "input",
                                         "Folder where to store the compressed files.",
                                         true,
                                         "",
                                         "string",cmd);
    TCLAP::ValueArg<std::string> outputArg("o",
                                          "output",
                                          "FASTA/Q file that should be compressed",
                                          true,
                                          "",
                                          "string",cmd);
    TCLAP::ValueArg<std::string> statisticsArg("s",
                                               "statistics",
                                               "Filename of json file that will output statistics about compressed file.",
                                               false,
                                               "statistics.json",
                                               "filename",cmd);
    TCLAP::ValueArg<int> numThreadsArg("p",
                                       "threads",
                                       "Number of threads to use for parallel computation.",
                                       false,
                                       1,
                                       "int",cmd);
    TCLAP::ValueArg<std::string> graphArg("g",
                                          "graph",
                                          "Graph to use as a reference in decompression",
                                          false,
                                          "",
                                          "string",cmd);
    cmd.parse(argc, argv);
    auto input_folder = inputArg.getValue();
    auto output_filename = outputArg.getValue();
    omp_set_num_threads(numThreadsArg.getValue());
    set_num_threads(numThreadsArg.getValue());
    auto statistics_filename = statisticsArg.getValue();
    auto db = PathDatabaseWavelet<>::deserialize(input_folder);
    auto reads = db.decode_all_reads();
    auto statistics = db.get_statistics(0);
    save_string(statistics.dump(4),statistics_filename);
    write_reads_to_fasta(reads,output_filename);

    return 0;
}
