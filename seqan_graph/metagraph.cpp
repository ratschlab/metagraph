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
    Shape<Dna5, UngappedShape<11> > hash_func;
    typedef Value<Shape<Dna5, UngappedShape<11> > >::Type THash;
    Map<Pair<THash, TVertexDescriptor> > kmer_map;
    TGraph g;

    unsigned int cnt_first = 0;
    unsigned int cnt_recurr = 0;

    while (!atEnd(stream)) {
        if (readRecord(id, seq, stream) != 0) {
            std::cerr << "ERROR while reading from " << config.fname << std::endl;
            return 1;
        }

        // iterate over k-mers of current sequence
        THash hash;
        hashInit(hash_func, begin(seq));
        TVertexDescriptor last_node;
        TVertexDescriptor current_node;
        for (unsigned i = 0; i < length(seq) - length(hash_func) + 1; ++i) {
            //if (i > 0 && i % 1000 == 0)
            //    return 0;
            if (i > 0 && i % 100000 == 0) {
                std::cout << "." << std::flush;
                //fprintf(stdout, ".");
                if (i % 1000000 == 0)
                    fprintf(stdout, "%i\n", i);
            }
            hash = hashNext(hash_func, begin(seq) + i);
            //std::cout << hash << " " << infix(seq, i, i+4) << " ";
            if (!hasKey(kmer_map, hash)) {
                current_node = addVertex(g);
                add(kmer_map, hash, current_node);
                ++cnt_first;
                //std::cout << "new" << std::endl;
            } else {
                current_node = kmer_map[hash];
                ++cnt_recurr;
                //std::cout << "recurr" << std::endl;
            }
            if (i > 1)
                addEdge(g, last_node, current_node);
            last_node = current_node;
        }
        std::cout << std::endl << "k-mers opening new nodes: " << cnt_first << std::endl;
        std::cout << "k-mers already present in nodes: " << cnt_recurr << std::endl;
        //std::cout << id << ":" << seq << std::endl;
    }



    return 0;
}
