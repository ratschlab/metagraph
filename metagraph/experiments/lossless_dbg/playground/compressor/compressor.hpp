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

#include "path_database_wavelet.hpp"
#include "samplers.hpp"
#include "utilities.hpp"
#include "threading.hpp"


#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"



using namespace std;
using namespace std::string_literals;


using node_index = SequenceGraph::node_index;


template<class Database>
Database compressReads(ValueArg<string> &compressedArg,
                       const vector<string> &filenames,
                       const vector<string> &reads) {
    auto db = Database(filenames);
    db.encode(reads);

    if (compressedArg.isSet()) {
        fs::path compress_folder = compressedArg.getValue();
        db.serialize(compress_folder);
    }
    return db;
}

template<class Database>
Database compressReadsDeprecated(ValueArg<string> &compressedArg,
                                 const vector<string> &reads,int kmer_length) {
    auto db = Database(reads,kmer_length);
    db.encode(reads);

    if (compressedArg.isSet()) {
        fs::path compress_folder = compressedArg.getValue();
        db.serialize(compress_folder);
    }
    return db;
}

int main_compressor(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Compress reads",' ', "0.1");

    //TODO:
    //TCLAP::ValueArg<std::string> inputArg("g",
    //                                      "graph",
    //                                      "DeBruijnGraph to use as a reference in compression",
    //                                      true,
    //                                      "",
    //                                      "string",cmd);

    TCLAP::ValueArg<std::string> inputArg("i",
                                         "input",
                                         "FASTA/Q file that should be compressed",
                                         true,
                                         "",
                                         "string",cmd);
    TCLAP::ValueArg<int> kmerLengthArg("k",
                                       "kmer-length",
                                       "Length of the kmers for graph construction",
                                       false,
                                       1,
                                       "int",cmd);

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
    TCLAP::ValueArg<int> numThreadsArg("p",
                                       "threads",
                                       "Number of threads to use for parallel computation.",
                                       false,
                                       1,
                                       "int",cmd);


    std::vector<std::string> regimes {
            "wavelet",
    };
    ValuesConstraint<std::string> regime_constraint(regimes);
    ValueArg<std::string> compressor_type("c",
                                          "compressor-type",
                                          "Which compressor to choose",
                                          false,
                                          "wavelet",
                                          &regime_constraint, cmd);
    cmd.parse(argc, argv);
    set_num_threads(numThreadsArg.getValue());
    auto input_filename = inputArg.getValue();
    auto statistics_filename = statisticsArg.getValue();
    // TODO: Don't read all reads to memory
    auto reads = read_reads_from_fasta(input_filename);
    auto kmer_length = kmerLengthArg.getValue();
    auto compressor = compressor_type.getValue();
    if (compressor == "wavelet") {
        auto db = compressReadsDeprecated<PathDatabaseWavelet<>>(compressedArg,
                                                       reads,kmer_length);
        if (statisticsArg.isSet()) {
            auto statistics = db.get_statistics();
            save_string(statistics.dump(4),statistics_filename);
        }
    }
    return 0;
}

