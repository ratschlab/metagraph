//
// Created by Jan Studený on 2019-03-11.
//
#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <tclap/CmdLine.h>
#include <omp.h>

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
//#define ALL_EDGES_COVERED

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
		const vector<string> &reads,int64_t kmer_length) {
	auto db = new Database(reads,kmer_length);
	db.encode(reads);

	if (compressedArg.isSet()) {
		fs::path compress_folder = compressedArg.getValue();
		db.serialize(compress_folder);
	}
	return db;
}

template<class DatabaseT>
void compress_reads(ValueArg<std::string> &graphArg, const ValueArg<std::string> &statisticsArg,
                    ValueArg<std::string> &compressedArg, const string &statistics_filename,
                    const vector<string> &reads, int64_t kmer_length, shared_ptr<BetterDBGSuccinct> &graph) {
    unique_ptr<DatabaseT> pd;
    if (graphArg.isSet()) {
        Timer timer;
        cerr << "Started loading the graph" << endl;
        graph->load(graphArg.getValue());
        cerr << "Finished loading the graph in " << timer.elapsed() << " sec." << endl;
        pd.reset(new DatabaseT(graph));
    } else {
        pd.reset(new DatabaseT(reads,kmer_length));
    }
    pd->encode(reads);
    if (compressedArg.isSet()) {
        fs::path compress_folder = compressedArg.getValue();
        pd->serialize(compress_folder);
    }
    if (statisticsArg.isSet()) {
        auto statistics = pd->get_statistics();
        save_string(statistics.dump(4),statistics_filename);
    }
}

int main_compressor(int argc, char *argv[]) {
	TCLAP::CmdLine cmd("Compress reads",' ', "0.1");

	TCLAP::ValueArg<std::string> graphArg("g",
			"graph",
			"Graph to use as a reference in compression",
			false,
			"",
			"string",cmd);
    std::vector<std::string> rerouting {
            "yes","no"
    };
    ValuesConstraint<std::string> rerouting_constraint(rerouting);
    TCLAP::ValueArg<std::string> pathReroutingArg("r",
                                          "path-rerouting",
                                          "Enable rerouting to improve space efficiency.",
                                          false,
                                          "yes",
                                          &rerouting_constraint,
                                          cmd);

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
			21,
			"int64_t",cmd);
	TCLAP::ValueArg<std::string> statisticsArg("s",
			"statistics",
			"Filename of json file that will output statistics about compressed file.",
			false,
			"statistics.json",
			"string",cmd);
	TCLAP::ValueArg<std::string> compressedArg("o",
			"output",
			"Folder where to store the compressed files.",
			true,
			"",
			"string",cmd);
	TCLAP::ValueArg<int> numThreadsArg("p",
			"threads",
			"Number of threads to use for parallel computation.",
			false,
			1,
			"int64_t",cmd);


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
    omp_set_num_threads(numThreadsArg.getValue());
	set_num_threads(numThreadsArg.getValue());
	auto input_filename = inputArg.getValue();
	auto statistics_filename = statisticsArg.getValue();
	bool use_transformations = pathReroutingArg.getValue() == "yes";
	// TODO: Don't read all reads to memory
    Timer read_timer;
    cerr << "Started loading the reads." << endl;
	auto reads = read_reads_from_fasta(input_filename);
    cerr << "Finished loading the reads in " << read_timer.elapsed() << " sec." << endl;
	auto kmer_length = kmerLengthArg.getValue();
	auto compressor = compressor_type.getValue();
    auto graph = std::make_shared<DBGSuccinct>(21);
	if (compressor == "wavelet") {
	    if (use_transformations) {
            compress_reads<PathDatabaseWavelet<>>(graphArg, statisticsArg, compressedArg, statistics_filename, reads, kmer_length, graph);
        }
	    else {
            compress_reads<PathDatabaseWaveletWithtoutTransformation<>>(graphArg, statisticsArg, compressedArg, statistics_filename, reads, kmer_length, graph);
	    }
    }
	return 0;
}

