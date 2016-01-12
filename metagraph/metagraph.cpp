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

#include <datatypes.hpp>
#include <dbg_succinct_libmaus.hpp>
#include <config.hpp>

// parse command line arguments and options
seqan::ArgumentParser::ParseResult
parseCommandLine(CFG & config, int argc, char const ** argv) 
{
    seqan::ArgumentParser parser("metagraph");

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
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT", true /*can be given multiple time*/, 1 /* min number of times*/));

    // add options
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Increase verbosity level of output"));
    addOption(parser, seqan::ArgParseOption("i", "integrate", "integrates fasta into given (-I) graph"));
    addOption(parser, seqan::ArgParseOption("k", "kmer_length", "Length of the k-mer to use", seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("d", "distance", "Max allowed alignment distance", seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption("O", "outfile_base", "basename for graph output files", seqan::ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("S", "sql_base", "basename for SQL output files", seqan::ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("I", "infile_base", "basename for loading graph input files", seqan::ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("m", "merge", "list of graph file basenames used as input for merging", seqan::ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("c", "compare", "list of graph file basenames used as input for comparison", seqan::ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("a", "align", "TODO", seqan::ArgParseArgument::STRING, "TEXT"));

    // set defaults
    setDefaultValue(parser, "O", "");
    setDefaultValue(parser, "I", "");
    setDefaultValue(parser, "S", "");
    setDefaultValue(parser, "m", "");
    setDefaultValue(parser, "c", "");
    setDefaultValue(parser, "a", "");
    setDefaultValue(parser, "d", 0);

    // parse command line
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // everything correct?
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    
    // fill config object 
    config.verbose = isSet(parser, "verbose");
    config.integrate = isSet(parser, "integrate");
    getOptionValue(config.k, parser, "k");
    getOptionValue(config.outfbase, parser, "O");
    getOptionValue(config.infbase, parser, "I");
    getOptionValue(config.sqlfbase, parser, "S");
    getOptionValue(config.merge, parser, "m");
    getOptionValue(config.compare, parser, "c");
    getOptionValue(config.align, parser, "a");
    getOptionValue(config.distance, parser, "d");
    config.fname = seqan::getArgumentValues(parser, 0);

    // do any logic / error checking here (e.g., mutually exclusive options)

    // return 
    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv) {

    typedef seqan::ModifiedAlphabet<seqan::Dna5, seqan::ModExpand<'X'> > Dna5F; 

    CFG config;

    seqan::ArgumentParser::ParseResult res = parseCommandLine(config, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    if (config.verbose)
        std::cout << seqan::CharString("Welcome to MetaGraph") << std::endl;

    // create graph object
    DBG_succ* graph = NULL;

    if (! config.compare.empty()) {
        int cnt = 0;
        std::string token;
        std::stringstream comparelist(config.compare);
        while (std::getline(comparelist, token, ',')) {
            if (cnt == 0) {
                std::cout << "Opening file " << token << std::endl;
                graph = new DBG_succ(token, config);
            } else {
                std::cout << "Opening file for comparison ..." << token << std::endl;
                DBG_succ* graph_ = new DBG_succ(token, config);
                bool identical = graph->compare(graph_);
                if (identical) {
                    std::cout << "Graphs are identical" << std::endl;
                } else {
                    std::cout << "Graphs are not identical" << std::endl;
                }
                delete graph_;
            }
            cnt++;
        }
    } else if (! config.merge.empty()) {
        int cnt = 0;
        std::string token;
        std::stringstream mergelist(config.merge);
        while (std::getline(mergelist, token, ',')) {
            if (cnt == 0) {
                std::cout << "Opening file " << token << std::endl;
                graph = new DBG_succ(token, config);
                //graph->print_seq();
            } else {
                std::cout << "Opening file for merging ..." << token << std::endl;

                DBG_succ* graph_ = new DBG_succ(token, config);
                //graph->merge(graph_);
                
                DBG_succ* graph__ = new DBG_succ(graph_->get_k(), config, false);
                graph__->merge(graph, graph_);
                delete graph;
                graph = graph__;

                delete graph_;
                std::cerr << "... done merging." << std::endl;
            }
            cnt++;
        }
        //graph->print_seq();
    }
    else if (! config.align.empty()) {

      // need these two inputs
      if (config.align.empty() || config.infbase.empty()) {
        std::cerr << "Requires <fasta file> and a <de bruijn graph> to align reads." << config.align << std::endl;
        exit(1);
      }

      DBG_succ* graph = new DBG_succ(config.infbase, config);

      seqan::CharString id;
      seqan::String<seqan::Dna5> seq;
      seqan::SequenceStream stream;

      // open stream to unput fasta
      open(stream, seqan::toCString(config.align), seqan::SequenceStream::READ, seqan::SequenceStream::FASTA);

      if (!seqan::isGood(stream)) {
        std::cerr << "ERROR while opening input file " << config.align << std::endl;
        exit(1);
      }
           
      while (!seqan::atEnd(stream)) {
          if (seqan::readRecord(id, seq, stream) != 0) {
              std::cerr << "!!!ERROR while reading from " << config.align << std::endl;
              exit(1);
          }

          graph->print_seq();

          uint64_t aln_len = seqan::length(seq);

          if (config.distance > 0) {
              std::vector<std::vector<HitInfo> > graphindices = graph->align_fuzzy(seq, aln_len, config.distance);

              for (size_t i = 0; i < graphindices.size(); ++i) {
                 //#std::cout << seqan::infix(seq, i, i + graph->get_k()) << ": ";
                  std::cout << seqan::infix(seq, i, i + aln_len) << ": ";
                  for (size_t l = 0; l < graphindices.at(i).size(); ++l) {
                      HitInfo curr_hit(graphindices.at(i).at(l));
                      for (size_t m = 0; m < curr_hit.path.size(); ++m) {
                          std::cout << curr_hit.path.at(m) << ':';
                      }
                      for (size_t m = curr_hit.rl; m <= curr_hit.ru; ++m) {
                          std::cout << m << " ";
                      }
                      std::cout << "[" << curr_hit.cigar << "] ";
                  }
                  std::cout << std::endl;
              }

          } else {
              std::vector<uint64_t> graphindices = graph->align(seq);

              for (size_t i = 0; i < graphindices.size(); ++i) {
                  for (uint64_t j = 0; j < graph->get_k(); ++j) {
                      std::cout << seq[i+j];
                  }
                  std::cout << ": " << graphindices[i] << std::endl;
              }
          }

      }

      // close stream
      seqan::close(stream);
    }
    else {
      if (! config.infbase.empty()) {
        graph = new DBG_succ(config.infbase, config);
      } else {
        graph = new DBG_succ(config.k, config);
      }

        if (config.infbase.empty() || config.integrate) {
            // read from fasta stream 
            seqan::CharString id;
            seqan::String<seqan::Dna5> seq;
            // iterate over input files
            for (unsigned int f = 0; f < seqan::length(config.fname); ++f) {

                if (config.verbose) {
                    std::cout << std::endl << "Parsing " << config.fname[f] << std::endl;
                }

                // open stream to fasta file
                seqan::SequenceStream stream;
                open(stream, seqan::toCString(config.fname.at(f)), seqan::SequenceStream::READ, seqan::SequenceStream::FASTA);
                if (!seqan::isGood(stream)) {
                    std::cerr << "ERROR while opening input file " << config.fname.at(f) << std::endl;
                    exit(1);
                }

                while (!seqan::atEnd(stream)) {
                    if (seqan::readRecord(id, seq, stream) != 0) {
                        std::cerr << "!!!ERROR while reading from " << config.fname.at(f) << std::endl;
                        exit(1);
                    }
                    // add all k-mers of seq to the graph
                    graph->add_seq(seq);
                }
                //graph->update_counters();
                //graph->print_stats();
                
                //fprintf(stdout, "current mem usage: %lu MB\n", get_curr_mem() / (1<<20));
            }
            //graph->print_seq();
            //seqan::String<Dna5F> test("ACC");
            //std::cout << "ACC: " << graph->index(test) << std::endl;
        }
    }

    // graph output
    if (!config.sqlfbase.empty())
        graph->toSQL();
    if (!config.outfbase.empty())
        graph->toFile();

    delete graph;

    return 0;
}
