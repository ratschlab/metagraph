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


//#define MASK_DUMMY_KMERS
//#define ALL_EDGES_COVERED
//#define DISABLE_PARALELIZATION

#include "path_database_wavelet.hpp"
#include "dynamic_path_database.hpp"
#include "reference_dynamic_incoming_table.hpp"
#include "reference_dynamic_routing_table.hpp"
#include "samplers.hpp"
#include "utilities.hpp"
#include "threading.hpp"
#include "configuration.hpp"



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
void compress_store_reads(ValueArg<std::string> &graphArg, const ValueArg<std::string> &statisticsArg,
                          ValueArg<std::string> &compressedArg, const string &statistics_filename,
                          const vector<string> &reads, shared_ptr<DBGSuccinct> &graph, int chunks,
                          uint64_t stat_verbosity) {
    unique_ptr<DatabaseT> pd;
        VerboseTimer graph_loading_timer("loading the graph");
        if (not graph->load(graphArg.getValue())) {
            throw std::runtime_error("Error loading graph");
        }
#ifdef MASK_DUMMY_KMERS
    mg::common::logger->info("Masking dummy kmers");
    graph->mask_dummy_kmers(1, false);
#else
    alt_assert(graph->get_node_sequence(1) == std::string(graph->get_k(), '$'));
#endif
    graph_loading_timer.finished();
        pd.reset(new DatabaseT(graph,chunks));
    pd->encode(reads);
    if (compressedArg.isSet()) {
        fs::path compress_folder = compressedArg.getValue();
        pd->serialize(compress_folder);
    }
    if (statisticsArg.isSet()) {
        auto statistics = pd->get_statistics(stat_verbosity);
        save_string(statistics.dump(4),statistics_filename);
    }
}

template<class DatabaseT>
void compress_reads(ValueArg<std::string> &graphArg, const ValueArg<std::string> &statisticsArg,
                    const string &statistics_filename, const vector<string> &reads, shared_ptr<DBGSuccinct> &graph,
                    uint64_t stat_verbosity) {
	unique_ptr<DatabaseT> pd;
		VerboseTimer graph_loading_timer("loading the graph");
        if (not graph->load(graphArg.getValue())) {
            throw std::runtime_error("Error loading graph");
        }
#ifdef MASK_DUMMY_KMERS
    graph->mask_dummy_kmers(1, false);
#else
    alt_assert(graph->get_node_sequence(1) == std::string(graph->get_k(), '$'));
#endif
    graph_loading_timer.finished();
		pd.reset(new DatabaseT(graph));
	pd->encode(reads);
	if (statisticsArg.isSet()) {
		auto statistics = pd->get_statistics(stat_verbosity);
		save_string(statistics.dump(4),statistics_filename);
	}
}

int main_compressor(int argc, char *argv[]) {
	TCLAP::CmdLine cmd("Compress reads",' ', "0.1");

	TCLAP::ValueArg<std::string> graphArg("g",
			"graph",
			"Graph to use as a reference in compression",
			true,
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
			"filename",cmd);
	TCLAP::ValueArg<std::string> statisticsArg("s",
			"statistics",
			"Filename of json file that will output statistics about compressed file.",
			false,
			"statistics.json",
        "filename",cmd);
    TCLAP::ValueArg<int> verbosityArg("x",
                                      "statistics-verbosity",
                                      "Level of detail of the statistics",
                                      false,
                                      0u,
                                      "int64_t",cmd);
	TCLAP::ValueArg<std::string> compressedArg("o",
			"output",
			"Folder where to store the compressed files.",
			false,
			"",
			"filename",cmd);
    TCLAP::ValueArg<int> reducedCoverageArg("v",
                                               "reduced-coverage",
                                               "Enable / disable reduced coverage technique",
                                               false,
                                               1,
                                               "filename",cmd);
	TCLAP::ValueArg<int> numThreadsArg("p",
			"threads",
			"Number of threads to use for parallel computation.",
			false,
			1,
			"int",cmd);
	TCLAP::ValueArg<int> chunksArg("u",
									   "chunks",
									   "Number of chunks for routing and incoming table (less chunks decreases memory consumption but increases compression time) -1 = inf disables chunking.",
									   false,
									   DefaultChunks,
									   "int",cmd);


	std::vector<std::string> regimes {
		"wavelet",
		"wavelet_debug",
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
    // TODO: Don't read all reads to memory #noted
    VerboseTimer read_timer("loading the reads");
	auto reads = read_reads_from_fasta(input_filename);
    read_timer.finished();
	auto compressor = compressor_type.getValue();
	auto chunks = chunksArg.getValue();
	auto stat_verbosity = verbosityArg.getValue();
    auto graph = std::make_shared<DBGSuccinct>(21);
    fs::create_directories(fs::path(compressedArg.getValue()) / "path_encoder.flag");
    if (compressor == "wavelet_debug") {
        compress_store_reads<DynamicPathDatabaseCore<
                IdentityComparator<DynamicRoutingTable<>, ReferenceDynamicRoutingTable<>>,
                IdentityComparator<DynamicIncomingTable<>, ReferenceDynamicIncomingTable<>>,
                IdentityComparator<ExitBarrier<>, ReferenceExitBarrier<>>
        >>(graphArg, statisticsArg, compressedArg, statistics_filename, reads, graph, chunks, stat_verbosity);
    }
	if (compressor == "wavelet") {
	    if (use_transformations) {
	    	if (chunks > 0) {
                compress_store_reads<PathDatabaseWavelet<>>(graphArg, statisticsArg, compressedArg, statistics_filename,
                                                            reads, graph, chunks, stat_verbosity);
	    	}
            else {
            	throw std::runtime_error("Using transformation is not implemented yet");
            }
        }
	    else {
	    	if (chunks > 0) {
	    	    if (reducedCoverageArg.getValue()) {
                    compress_store_reads<PathDatabaseWaveletWithtoutTransformation<true>>(graphArg, statisticsArg,
                                                                                          compressedArg,
                                                                                          statistics_filename, reads,
                                                                                          graph, chunks,
                                                                                          stat_verbosity);
	    	    }
	    	    else {
                    compress_store_reads<PathDatabaseWaveletWithtoutTransformation<false>>(graphArg, statisticsArg,
                                                                                           compressedArg,
                                                                                           statistics_filename, reads,
                                                                                           graph, chunks,
                                                                                           stat_verbosity);
	    	    }
	    	}
			else {
			    if (not reducedCoverageArg.getValue()) { throw std::runtime_error("Not implemented"); }
                compress_store_reads<PathDatabaseWaveletWithtoutTransformation<true, ReferenceDynamicRoutingTable<>, ReferenceDynamicIncomingTable<>>>(
                        graphArg, statisticsArg, compressedArg, statistics_filename, reads, graph, chunks,
                        stat_verbosity);
			}
	    }
    }
	return 0;
}

