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

#include "dbg_succinct.hpp"
#include "sequence_graph.hpp"
#include "dbg_succinct_construct.hpp"
#include "dbg_hash.hpp"


#include "path_database_baseline_wavelet.hpp"
#include "samplers.hpp"
#include "utilities.hpp"

#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"



using namespace std;
using namespace std::string_literals;


using node_index = SequenceGraph::node_index;



int main(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Compress reads",' ', "0.1");
    TCLAP::ValueArg<std::string> inputArg("i",
                                         "input",
                                         "FASTA/Q file that should be compressed",
                                         true,
                                         "",
                                         "string");
    TCLAP::ValueArg<std::string> statisticsArg("s",
                                          "statistics",
                                          "Filename of json file that will output statistics about compressed file.",
                                          false,
                                          "statistics.json",
                                          "string");
    TCLAP::ValueArg<fs::path> compressedArg("o",
                                          "output",
                                          "Folder where to store the compressed files.",
                                          false,
                                          "",
                                          "string");
    cmd.add(inputArg);
    cmd.add(statisticsArg);
    cmd.add(compressedArg);
    cmd.parse(argc, argv);
    auto input_filename = inputArg.getValue();
    auto statistics_filename = statisticsArg.getValue();
    auto reads = read_reads_from_fasta(input_filename);
    auto db = PathDatabaseBaselineWavelet(reads);
    db.encode(reads);
    if (statisticsArg.isSet()) {
        //auto statistics = db.get_statistics();
        throw std::runtime_error("Not supported yet");
        //save_string(statistics.dump(4),statistics_filename);
    }
    if (compressedArg.isSet()) {
        fs::path compress_folder = compressedArg.getValue();
        db.serialize(compress_folder);
    }

    return 0;
}
