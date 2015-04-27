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

using namespace seqan;

struct CFG 
{
    String<char> fname;
    bool verbose;

    CFG() :
        verbose(false)
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
    addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "TEXT"));

    // add options
    addOption(parser, ArgParseOption("v", "verbose", "Increase verbosity level of output"));

    // parse command line
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // everything correct?
    if (res != ArgumentParser::PARSE_OK)
        return res;
    
    // fill config object 
    getOptionValue(config.verbose, parser, "verbose");
    getArgumentValue(config.fname, parser, 0); 

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

    // open stream to fasta file
    SequenceStream stream;
    open(stream, toCString(config.fname), SequenceStream::READ, SequenceStream::FASTA);
    if (!isGood(stream)) {
        std::cerr << "ERROR while opening input file " << config.fname << std::endl;
        return 1;
    }

    // build the graph from the k-mers in the string
    typedef unsigned int TCargo;
    typedef Graph<Directed<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    // read from fasta stream 
    CharString id;
    String<Dna5> seq;
    Map<Pair<String<Dna5>, TVertexDescriptor> > kmer_map;
    TGraph g;
    while (!atEnd(stream)) {
        if (readRecord(id, seq, stream) != 0) {
            std::cerr << "ERROR while reading from " << config.fname << std::endl;
            return 1;
        }

        // iterate over k-mers of current sequence
        Shape<Dna5, UngappedShape<31> > kmer_hash;
        hashInit(kmer_hash, begin(seq));
        TVertexDescriptor last_node;
        for (unsigned i = 1; i < length(seq) - length(kmer_hash) + 1; ++i) {
            if (i % 1000000 == 0)
                fprintf(stdout, "%i\n", i);
            TVertexDescriptor current_node = addVertex(g);
            add(kmer_map, hashNext(kmer_hash, begin(seq) + i), current_node);
            if (i > 1)
                addEdge(g, last_node, current_node);
            last_node = current_node;
        }

        //std::cout << id << ":" << seq << std::endl;
    }



    return 0;
}
