#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/map.h>
#include <seqan/index.h>

#include <vector>
#include <string>

#include <btree_map.h>
#include <dbg_seqan.hpp>
#include <dbg_succinct_libmaus.hpp>
//#include <dmm_tree.hpp>
#include <unix_tools.hpp>

//#include <lemon/list_graph.h>

using namespace seqan;

struct CFG 
{
    //String<char> > fname;
    std::vector<std::string> fname;
    bool verbose;
    unsigned int k;

    CFG() :
        verbose(false), k(31)
    {}
};

// parse command line arguments and options
ArgumentParser::ParseResult
parseCommandLine(CFG & config, int argc, char const ** argv) 
{
    ArgumentParser parser("metagraph");

    // add program meta information
    setShortDescription(parser, "A comprehensive graph represenatation of metagenome information");
    setVersion(parser, "0.1");

    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] <\\fIFasta filename\\fP>");
    addDescription(parser,
                  "This program is the first implementation of a "
                  "meta-metagenome graph for identifciation and annotation "
                  "purposes.");

    // add arguments
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "TEXT", 1));

    // add options
    addOption(parser, ArgParseOption("v", "verbose", "Increase verbosity level of output"));
    addOption(parser, ArgParseOption("k", "kmer_length", "Length of the k-mer to use", ArgParseArgument::INTEGER, "INT"));

    // parse command line
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // everything correct?
    if (res != ArgumentParser::PARSE_OK)
        return res;
    
    // fill config object 
    getOptionValue(config.verbose, parser, "verbose");
    getOptionValue(config.k, parser, "k");
    config.fname = getArgumentValues(parser, 0);

    // do any logic / error checking here (e.g., mutually exclusive options)

    // return 
    return ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv) {
    
    // command line parsing
    CFG config;
    ArgumentParser::ParseResult res = parseCommandLine(config, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    if (config.verbose)
        std::cout << CharString("Welcome to MetaGraph") << std::endl;

    // build the graph from the k-mers in the string
    //typedef unsigned int TCargo;
    //typedef Graph<Directed<TCargo> > TGraph; // --> we can use this when we want to have edge weights

    // create graph object
    //DBG_seqan* graph = new DBG_seqan(config.k);
    DBG_succ* graph = new DBG_succ(config.k);

    // read from fasta stream 
    CharString id;
    String<Dna5> seq;

    // iterate over input files
    for (unsigned int f = 0; f < length(config.fname); ++f) {

        if (config.verbose) {
            std::cout << std::endl << "Parsing " << config.fname[f] << std::endl;
        }

        // open stream to fasta file
        SequenceStream stream;
        open(stream, toCString(config.fname.at(f)), SequenceStream::READ, SequenceStream::FASTA);
        if (!isGood(stream)) {
            std::cerr << "ERROR while opening input file " << config.fname.at(f) << std::endl;
            return 1;
        }

        while (!atEnd(stream)) {
            if (readRecord(id, seq, stream) != 0) {
                std::cerr << "ERROR while reading from " << config.fname.at(f) << std::endl;
                return 1;
            }
            // add all k-mers of seq to the graph
            graph->add_seq(seq);
        }
        //graph->update_counters();
        //graph->print_stats();
        fprintf(stdout, "current mem usage: %lu MB\n", get_curr_mem() / (1<<20));
    }

    delete graph;

    return 0;
}
