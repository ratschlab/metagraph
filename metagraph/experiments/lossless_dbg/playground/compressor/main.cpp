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

#define MEMOIZE
#define ALL_EDGES_COVERED

#include "path_database_baseline_wavelet_deprecated.hpp"
#include "path_database_baseline_wavelet.hpp"
#include "samplers.hpp"
#include "utilities.hpp"


#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"



using namespace std;
using namespace std::string_literals;


using node_index = SequenceGraph::node_index;


template<class Database>
Database compressReads(ValueArg<string> &compressedArg, const vector<string> &reads) {
    auto db = Database(reads);
    db.encode(reads);

    if (compressedArg.isSet()) {
        fs::path compress_folder = compressedArg.getValue();
        db.serialize(compress_folder);
    }
    return db;
}

int main(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Compress reads",' ', "0.1");
    TCLAP::ValueArg<std::string> inputArg("i",
                                         "input",
                                         "FASTA/Q file that should be compressed",
                                         true,
                                         "",
                                         "string",cmd);
    TCLAP::ValueArg<std::string> statisticsArg("s",
                                          "statistics",
                                          "Filename of json file that will output statistics about compressed file.",
                                          false,
                                          "statistics.json",
                                          "string",cmd);
    TCLAP::ValueArg<std::string> compressedArg("o",
                                          "output",
                                          "Folder where to store the compressed files.",
                                          false,
                                          "",
                                          "string",cmd);
    std::vector<std::string> regimes {
            "wavelet",
            "wavelet_old",
    };
    ValuesConstraint<std::string> regime_constraint(regimes);
    ValueArg<std::string> compressor_type("c",
                                          "compressor-type",
                                          "Which compressor to choose",
                                          false,
                                          "wavelet",
                                          &regime_constraint, cmd);
    cmd.parse(argc, argv);
    auto input_filename = inputArg.getValue();
    auto statistics_filename = statisticsArg.getValue();
    auto reads = read_reads_from_fasta(input_filename);
    auto compressor = compressor_type.getValue();
    if (compressor == "wavelet") {
        auto db = compressReads<PathDatabaseBaselineWavelet<>>(compressedArg, reads);
        if (statisticsArg.isSet()) {
            auto statistics = db.get_statistics();
            save_string(statistics.dump(4),statistics_filename);
        }
    }
    else if (compressor == "wavelet_old") {
        compressReads<PathDatabaseBaselineWaveletDeprecated>(compressedArg, reads);
    }

    return 0;
}

