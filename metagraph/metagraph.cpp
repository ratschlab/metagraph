#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <zlib.h>

#include "datatypes.hpp"
#include "dbg_succinct_libmaus.hpp"
#include "config.hpp"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

Config* config;

int main(int argc, char const ** argv) {

    // parse command line arguments and options
    config = new Config(argc, argv);

    if (config->verbose)
        std::cout << "Welcome to MetaGraph" << std::endl;

    // create graph pointer
    DBG_succ* graph = NULL;

    switch (config->identity) {

        case Config::compare: {
            for (unsigned int f = 0; f < config->fname.size(); ++f) {
                if (f == 0) {
                    std::cout << "Opening file " << config->fname.at(f) << std::endl;
                    graph = new DBG_succ(config->fname.at(f), config);
                } else {
                    std::cout << "Opening file for comparison ..." << config->fname.at(f) << std::endl;
                    DBG_succ* graph_ = new DBG_succ(config->fname.at(f), config);
                    bool identical = graph->compare(graph_);
                    if (identical) {
                        std::cout << "Graphs are identical" << std::endl;
                    } else {
                        std::cout << "Graphs are not identical" << std::endl;
                    }
                    delete graph_;
                }
            }
        } break;

        case Config::merge: {
            for (unsigned int f = 0; f < config->fname.size(); ++f) {
                if (f == 0) {
                    std::cout << "Opening file " << config->fname.at(f) << std::endl;
                    graph = new DBG_succ(config->fname.at(f), config);
                } else {
                    std::cout << "Opening file for merging ..." << config->fname.at(f) << std::endl;
                    DBG_succ* graph_ = new DBG_succ(config->fname.at(f), config);
                    
                    DBG_succ* graph__ = new DBG_succ(graph_->get_k(), config, false);
                    graph__->merge(graph, graph_);
                    delete graph;
                    graph = graph__;

                    delete graph_;
                    std::cerr << "... done merging." << std::endl;
                }
            }
        } break;

        case Config::stats: {
            std::ofstream outstream;
            if (!config->outfbase.empty()) {
                outstream.open((config->outfbase + ".stats.dbg").c_str());
                outstream << "file\tnodes\tedges\tk" << std::endl;
            }
            for (unsigned int f = 0; f < config->fname.size(); ++f) {
                DBG_succ* graph_ = new DBG_succ(config->fname.at(f), config);
                if (!config->quiet) {
                    std::cout << "Statistics for file " << config->fname.at(f) << std::endl;
                    std::cout << "nodes: " << graph_->get_node_count() << std::endl;
                    std::cout << "edges: " << graph_->get_edge_count() << std::endl;
                    std::cout << "k: " << graph_->get_k() << std::endl; 
                }
                if (!config->outfbase.empty()) {
                    outstream << config->fname.at(f) << "\t" 
                              << graph_->get_node_count() << "\t" 
                              << graph_->get_edge_count() << "\t"
                              << graph_->get_k() << std::endl;
                }
                if (config->print_graph)
                    graph_->print_seq();
                delete graph_;
            }
            if (!config->outfbase.empty())
                outstream.close();
        } break;


        case Config::align: {
            // check required inputs
            if (config->infbase.empty()) {
              std::cerr << "Requires input <de bruijn graph> to align reads." << config->align << std::endl;
              exit(1);
            }

            DBG_succ* graph = new DBG_succ(config->infbase, config);

            for (unsigned int f = 0; f < config->fname.size(); ++f) {
                std::cout << "Opening file for alignment ..." << config->fname.at(f) << std::endl;

                // open stream to input fasta
                gzFile input_p = gzopen(config->fname.at(f).c_str(), "r");
                kseq_t *read_stream = kseq_init(input_p);
                if (read_stream == NULL) {
                  std::cerr << "ERROR while opening input file " << config->align << std::endl;
                  exit(1);
                }

                while (kseq_read(read_stream) >= 0) {

                    //graph->print_seq();
                    uint64_t aln_len = read_stream->seq.l;

                    if (config->distance > 0) {
                        std::vector<std::vector<HitInfo> > graphindices = graph->align_fuzzy(read_stream->seq, aln_len, config->distance);

                        for (size_t i = 0; i < graphindices.size(); ++i) {
                            int print_len = (i + aln_len < read_stream->seq.l) ? aln_len : (read_stream->seq.l - i);
                            printf("%.*s: ", print_len, read_stream->seq.s + i);

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
                        std::vector<uint64_t> graphindices = graph->align(read_stream->seq);

                        for (size_t i = 0; i < graphindices.size(); ++i) {
                            for (uint64_t j = 0; j < graph->get_k(); ++j) {
                                std::cout << read_stream->seq.s[i+j];
                            }
                            std::cout << ": " << graphindices[i] << std::endl;
                        }
                    }
                }

                // close stream
                kseq_destroy(read_stream);
                gzclose(input_p);
            }
        } break;

        case Config::build: {
            if (!config->infbase.empty()) {
                graph = new DBG_succ(config->infbase, config);
            } else {
                graph = new DBG_succ(config->k, config);
            }

            // iterate over input files
            for (unsigned int f = 0; f < config->fname.size(); ++f) {

                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << config->fname[f] << std::endl;
                }

                // open stream to fasta file
                gzFile input_p = gzopen(config->fname[f].c_str(), "r");
                kseq_t *read_stream = kseq_init(input_p);

                if (read_stream == NULL) {
                    std::cerr << "ERROR while opening input file " << config->fname.at(f) << std::endl;
                    exit(1);
                }

                while (kseq_read(read_stream) >= 0) {
                    // add all k-mers of seq to the graph
                    graph->add_seq(read_stream->seq);
                }
                kseq_destroy(read_stream);
                gzclose(input_p);

                //graph->update_counters();
                //graph->print_stats();
                
                //fprintf(stdout, "current mem usage: %lu MB\n", get_curr_mem() / (1<<20));
            }
            //graph->print_seq();
        } break;
    }

    // output and cleanup
    if (graph) {
        // graph output
        if (config->print_graph)
            graph->print_seq();
        if (!config->sqlfbase.empty())
            graph->toSQL();
        if (!config->outfbase.empty())
            graph->toFile();

        delete graph;
    }
    delete config;

    return 0;
}
