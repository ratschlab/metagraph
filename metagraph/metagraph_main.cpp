#include <iostream>
#include <fstream>
#include <ctime>

#include <vector>
#include <set>
#include <deque>
#include <string>
#include <zlib.h>
#include <pthread.h>
#include <map>

#include "datatypes.hpp"
#include "dbg_succinct_libmaus.hpp"
#include "config.hpp"
#include "helpers.hpp"
#include "kseq.h"
#include "utils.hpp"
#include "vcfparse.h"
#include "compare.hpp"
#include "traverse.hpp"
#include "merge.hpp"
#include "annotate.hpp"
#include "construct.hpp"

KSEQ_INIT(gzFile, gzread)

Config* config;
ParallelMergeContainer* merge_data;
ParallelAnnotateContainer* anno_data = new ParallelAnnotateContainer();

pthread_mutex_t mutex_merge_result = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_bin_idx = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_annotate = PTHREAD_MUTEX_INITIALIZER;
pthread_attr_t attr;  

void parallel_merge_collect(DBG_succ* result) {

    for (size_t i = 0; i < merge_data->result.size(); ++i) {
        if (merge_data->result.at(i)) {
            construct::append_graph(result, merge_data->result.at(i));
            delete merge_data->result.at(i);
        }
    }
    merge_data->result.clear();
    merge_data->bins_done = 0;
    result->p = result->succ_W(1, 0);
}

bool operator==(const AnnotationSet& lhs, const AnnotationSet& rhs) {
    return (lhs.annotation == rhs.annotation); 
}

void get_RAM() {
    //output total RAM usage
    FILE* sfile = fopen("/proc/self/status","r");
    char line[128];
    while (fgets(line, 128, sfile) != NULL) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            break;
        }
    }
    std::cout << line;
    fclose(sfile);
}


/*
 * Distribute the annotation of all k-mers in a sequence over 
 * a number of parallel bins.
 */
void *parallel_annotate_wrapper(void *arg) {

    uint64_t curr_idx, start, end;

    while (true) {
        pthread_mutex_lock (&mutex_bin_idx);
        if (anno_data->idx == anno_data->total_bins) {
            pthread_mutex_unlock (&mutex_bin_idx);
            break;
        } else {
            curr_idx = anno_data->idx;
            anno_data->idx++;
            pthread_mutex_unlock (&mutex_bin_idx);
            
            start = curr_idx * anno_data->binsize;
            end = std::min(((curr_idx + 1) * anno_data->binsize) + anno_data->graph->k - 1, anno_data->seq->l);
            //std::cerr << "start " << start << " end " << end << std::endl;
            annotate::annotate_seq(anno_data->graph, *(anno_data->seq), *(anno_data->label), start, end, &mutex_annotate);
        }
    }
    pthread_exit((void*) 0);
}


/*
 * Distribute the merging of a set of graph structures over
 * bins, such that n parallel threads are used. 
 */
void *parallel_merge_wrapper(void *arg) {

    unsigned int curr_idx;
    DBG_succ* graph;

    while (true) {
        pthread_mutex_lock (&mutex_bin_idx);
        if (merge_data->idx == merge_data->ref_bins.size()) {
            pthread_mutex_unlock (&mutex_bin_idx);
            break;
        } else {
            curr_idx = merge_data->idx;
            merge_data->idx++;
            pthread_mutex_unlock (&mutex_bin_idx);
            
            graph = new DBG_succ(merge_data->k, config, false);
            std::vector<uint64_t> kv;
            std::vector<uint64_t> nv;                                                                                                            
            for (size_t i = 0; i < merge_data->graphs.size(); i++) {
                kv.push_back(merge_data->bins.at(i).at(curr_idx).first);
                nv.push_back(merge_data->bins.at(i).at(curr_idx).second + 1);
            }
            merge::merge(graph, merge_data->graphs, kv, nv, true);

            pthread_mutex_lock (&mutex_merge_result);
            merge_data->result.at(curr_idx) = graph;
            merge_data->bins_done++;
            if (config->verbose) 
                std::cout << "finished bin " << curr_idx + 1 << " (" << merge_data->bins_done << "/" << merge_data->ref_bins.size() << ")" << std::endl;
            pthread_mutex_unlock (&mutex_merge_result);

        }
    }
    pthread_exit((void*) 0);
}
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
                    bool identical = compare::compare(graph, graph_);
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

            // collect results on an external merge
            if (config->collect > 1) {
                std::string fname;
                for (uint64_t f = 0; f < config->collect; f++) {
                    fname = config->outfbase + "." + std::to_string(f) + "_" + std::to_string(config->collect);
                    std::cout << "Opening file " << fname << std::endl;
                    if (config->fast) {
                        if (f == 0) {
                            graph = new DBG_succ(utils::kFromFile(fname), config, false);
                            graph->last_stat.push_back(0);
                            graph->W_stat.push_back(0);
                        }
                        DBG_succ* graph_to_append = new DBG_succ(fname, config);
                        construct::append_graph_static(graph, graph_to_append);
                        delete graph_to_append;
                    } else {
                        if (f == 0) {
                            graph = new DBG_succ(fname, config);
                        } else {
                            DBG_succ* graph_to_append = new DBG_succ(fname, config);
                            construct::append_graph(graph, graph_to_append);
                            delete graph_to_append;
                        }
                    }
                }
                if (config->fast) {
                    graph->toDynamic();
                }
                graph->p = graph->succ_W(1, 0);

            // run normal merge procedure
            } else {

                // some preliminaries to make command line options consistent
               // if ((config->parts_total > 1) && (config->parts_total > (config->parallel * config->bins_per_thread)))
                //    config->bins_per_thread = config->parts_total / config->parallel;
                            
                std::vector<DBG_succ*> graphs;
                std::vector<uint64_t> kv;
                std::vector<uint64_t> nv;
                for (unsigned int f = 0; f < config->fname.size(); ++f) {
                        std::cout << "Opening file " << config->fname.at(f) << std::endl;
                        graph = new DBG_succ(config->fname.at(f), config);
                        graphs.push_back(graph);
                        kv.push_back(1);
                        nv.push_back(graph->get_size());
                }
                DBG_succ* target_graph = new DBG_succ(graphs.front()->get_k(), config, false);

                if ((config->parallel > 1) || (config->parts_total > 1)) {
                    pthread_t* threads = NULL; 
                    merge_data = new ParallelMergeContainer();

                    // get bins in graphs according to required threads
                    if (config->verbose)
                        std::cout << "Collecting reference bins" << std::endl;
                    std::cerr << "parallel " << config->parallel << " per thread " << config->bins_per_thread << " parts total " << config->parts_total << std::endl;
                    merge_data->ref_bins = merge::get_bins(graphs.front(), config->parallel * config->bins_per_thread * config->parts_total);

                    // only work on subset of the bins when requested
                    if (config->parts_total > 1) {
                        merge_data->subset_bins(config->part_idx, config->parts_total, config->parallel * config->bins_per_thread);
                    }
                    merge_data->bins.push_back(merge_data->ref_bins);

                    if (config->verbose)
                        std::cout << "Collecting relative bins" << std::endl;
                    for (size_t i = 1; i < graphs.size(); i++) {
                        std::cerr << "getting bins for " << i << ": " << config->fname.at(i) << std::endl; 
                        merge_data->bins.push_back(merge::get_bins_relative(graphs.at(i), graphs.front(), merge_data->ref_bins, merge_data->first, merge_data->last));
                    }
                    for (size_t i = 0; i < graphs.size(); i++) {
                        for (size_t ii = 0; ii < merge_data->bins.at(i).size(); ii++) {
                            if (merge_data->bins.at(i).at(ii).first > merge_data->bins.at(i).at(ii).second) {
                               merge_data->bins.at(i).at(ii) = std::make_pair(graphs.at(i)->get_size(), graphs.at(i)->get_size() - 1); 
                            }
                        }
                    }

                    // print bin stats
                    if (config->verbose) {
                        merge_data->get_bin_stats();
                    }

                    // prepare data shared by threads
                    merge_data->idx = 0;
                    merge_data->k = graph->get_k();
                    merge_data->graphs = graphs;
                    for (size_t i = 0; i < merge_data->ref_bins.size(); i++)
                        merge_data->result.push_back(NULL);
                    merge_data->bins_done = 0;

                    // create threads
                    threads = new pthread_t[config->parallel]; 
                    pthread_attr_init(&attr);
                    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

                    // do the work
                    for (size_t tid = 0; tid < config->parallel; tid++) {
                       pthread_create(&threads[tid], &attr, parallel_merge_wrapper, (void *) tid); 
                       std::cerr << "starting thread " << tid << std::endl;
                    }

                    // join threads
                    if (config->verbose)
                        std::cout << "Waiting for threads to join" << std::endl;
                    for (size_t tid = 0; tid < config->parallel; tid++) {
                        pthread_join(threads[tid], NULL);
                    }
                    delete[] threads;

                    // collect results
                    std::cerr << "Collecting results" << std::endl;
                    parallel_merge_collect(target_graph);

                    graph = target_graph;

                    delete merge_data;

                } else {
                    merge::merge(target_graph, graphs, kv, nv);
                    graph = target_graph;
                }
                for (size_t f = 0; f < graphs.size(); f++)
                    delete graphs.at(f);
                std::cerr << "... done merging." << std::endl;
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
                /*graph_->W = new libmaus2::wavelet::DynamicWaveletTree<6, 64> (3);
                graph_->W->insert(1ull, 0);
                graph_->W->insert(7ull, 1);
                graph_->W->insert(4ull, 2);
                graph_->W->insert(3ull, 3);
                graph_->W->insert(2ull, 4);
                graph_->W->insert(5ull, 5);

                exit(1);
                */
                if (!config->quiet) {
                    std::cout << "Statistics for file " << config->fname.at(f) << std::endl;
                    std::cout << "nodes: " << graph_->get_node_count() << std::endl;
                    std::cout << "edges: " << graph_->get_edge_count() << std::endl;
                    std::cout << "k: " << graph_->get_k() << std::endl; 
                    //graph_->traversalHash();
                    /*std::deque<uint64_t> tmp;
                    tmp.push_back(4);
                    tmp.push_back(6);
                    tmp.push_back(6);
                    tmp.push_back(6);
                    tmp.push_back(6);
                    std::cout << graph_->index_predecessor(tmp) << std::endl; */
                }
                if (!config->outfbase.empty()) {
                    outstream << config->fname.at(f) << "\t" 
                              << graph_->get_node_count() << "\t" 
                              << graph_->get_edge_count() << "\t"
                              << graph_->get_k() << std::endl;
                }
                if (config->print_graph)
                    graph_->print_seq();
                if (config->print_graph_succ)
                    graph_->print_state_str();

                /*DBG_succ* graph_tut = new DBG_succ(config->k, config);
                std::cerr << "inserting step by step" << std::endl;
                for (size_t i = 0; i < 10000000; ++i)
                    graph_tut->last->insertBit(i % 2, i);
                std::cerr << "done" << std::endl;

                std::cerr << "construct anew" << std::endl;
                graph_tut->last = new libmaus2::bitbtree::BitBTree<6, 64>(10000000ull, false); 
                for (size_t i = 0; i < 10000000; ++i)
                    graph_tut->last->setBitQuick(i, i % 2);
                std::cerr << "done" << std::endl;
                */

                /*std::cerr << graph_->get_edge_count() << std::endl;
                std::vector<uint64_t> result = graph_->split_range(1, graph_->get_edge_count(), 0);
                std::cout << "Results for d = 0" << std::endl;
                for (size_t i = 0; i < result.size(); ++i) {
                    std::cout << "  " << result.at(i) << std::endl;
                }
                result = graph_->split_range(2, 30, 1);
                std::cout << "Results for d = 1" << std::endl;
                for (size_t i = 0; i < result.size(); ++i) {
                    std::cout << "  " << result.at(i) << std::endl;
                }
                result = graph_->split_range(12, 20, 2);
                std::cout << "Results for d = 2" << std::endl;
                for (size_t i = 0; i < result.size(); ++i) {
                    std::cout << "  " << result.at(i) << std::endl;
                }*/

                delete graph_;
            }
            if (!config->outfbase.empty())
                outstream.close();

        } break;

        case Config::dump: {
            //for (unsigned int f = 0; f < config->fname.size(); ++f) {
                //DBG_succ* graph_ = new DBG_succ(config->fname.at(f), config);
                DBG_succ* graph_ = new DBG_succ(config->infbase, config);
                // checks whether annotation exists and creates an empty one if not
                graph_->annotationFromFile();
                graph_->print_adj_list();
                delete graph_;
            //}
        } break;


        case Config::align: {

            // load graph 
            if (config->infbase.empty()) {
              std::cerr << "Requires input <de bruijn graph> to align reads." << config->align << std::endl;
              exit(1);
            }
            graph = new DBG_succ(config->infbase, config);

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

            if (config->verbose)
                std::cerr << "k is " << config->k << std::endl;

            //enumerate all suffices
            unsigned int suffix_len = (unsigned int)ceil(log2(config->nsplits)/log2(graph->alph_size-1));
            //should be set to at most k-1 so that W calculation is correct
            suffix_len = std::min(suffix_len, (unsigned int)graph->k-1);
            std::deque<std::string> suffices = {""};
            for (size_t i=0;i<suffix_len;++i) {
                while (suffices[0].length() < suffix_len) {
                    for (size_t j=0;j<graph->alph_size;++j) {
                        suffices.push_back(graph->alphabet[j] + suffices[0]);
                    }
                    suffices.pop_front();
                }
            }
            assert(suffices.size() == pow(graph->alph_size, suffix_len));

            //generate sink node
            std::string starts = std::string(graph->k-1, graph->alphabet[graph->alph_size-1]) + graph->alphabet[0] + graph->alphabet[0];
            std::string sinks = graph->alphabet.substr(0,1);
            kstring_t graphsink = {1, 1, &sinks[0u]};
            kstring_t start = {graph->k+1, graph->k+1, &starts[0u]};

            clock_t tstart, timelast;
            
            //one pass per suffix
            for (size_t j=0;j<suffices.size();++j) {
                std::cout << "Suffix: " << suffices[j] << "\n";
                //add sink nodes
                construct::add_seq_alt(graph, start, false, config->parallel, suffices[j]);
                construct::add_seq_alt(graph, graphsink, true, config->parallel, suffices[j]);
                if (suffices[j].find("$") == std::string::npos) {
                    // iterate over input files
                    for (unsigned int f = 0; f < config->fname.size(); ++f) {
                        if (config->verbose) {
                            std::cout << std::endl << "Parsing " << config->fname[f] << std::endl;
                        }
                        // open stream
                        gzFile input_p = gzopen(config->fname[f].c_str(), "r");

                        //Get extension of file
                        int dotind = config->fname.at(f).rfind(".");
                        std::string ext = "";
                        if (dotind >= 0) {
                            if (config->fname.at(f).substr(dotind) == ".gz") {
                                int nextind = config->fname.at(f).substr(0,dotind-1).rfind(".");
                                ext = config->fname.at(f).substr(nextind, dotind-nextind);
                            } else {
                                ext = config->fname.at(f).substr(dotind);
                            }
                        }
                        if (dotind >= 0 && ext == ".vcf") {
                            //READ FROM VCF
                            uint64_t nbp = 0;
                            uint64_t nbplast = 0;
                            tstart = clock();
                            timelast = clock();
                            vcfparse *vcf = vcf_init(config->refpath.c_str(), config->fname.at(f).c_str(), graph->k);
                            if (!vcf) {
                                std::cerr << "ERROR reading VCF " << config->fname.at(f) << std::endl;
                                exit(1);
                            }
                            std::cerr << "Loading VCF with " << config->parallel << " threads per line\n";
                            for (size_t i=1; vcf_get_seq(vcf);++i) {
                                if (i % 10000 == 0) {
                                    std::cout << "." << std::flush;
                                    if (i % 100000 == 0) {
                                        fprintf(stdout, "%lu - bp %lu / runtime %lu / BPph %lu\n", i, nbp, (clock() - tstart)/CLOCKS_PER_SEC, 60*60*CLOCKS_PER_SEC*(nbp-nbplast)/(clock()-timelast));
                                        nbplast = nbp;
                                        timelast = clock();
                                    }
                                }
                                nbp += vcf->seq.l;
                                construct::add_seq_alt(graph, vcf->seq, false, config->parallel, suffices[j]);
                            }
                            vcf_destroy(vcf);
                        } else {
                            //READ FROM FASTA
                            //TODO: handle read_stream->qual
                            kseq_t *read_stream = kseq_init(input_p);
                            if (read_stream == NULL) {
                                std::cerr << "ERROR while opening input file " << config->fname.at(f) << std::endl;
                                exit(1);
                            }
                            for (size_t i=1;kseq_read(read_stream) >= 0; ++i) {
                            //while (kseq_read(read_stream) >= 0) {
                                // possibly reverse k-mers
                                if (config->reverse)
                                    reverse_complement(read_stream->seq);                    
                                // add all k-mers of seq to the graph
                                //if CID==
                                //0: normal
                                //1: bridge (from reference)
                                //even>1: normal from read
                                //odd>1: bridge from read
                                construct::add_seq_alt(graph, read_stream->seq, true, config->parallel, suffices[j], read_stream->qual.l ? i*2 : 0);
                            }
                            kseq_destroy(read_stream);
                        }
                        gzclose(input_p);
                        //graph->update_counters();
                        //graph->print_stats();
                        //fprintf(stdout, "current mem usage: %lu MB\n", get_curr_mem() / (1<<20));
                    }
                }
                get_RAM();
                //append to succinct representation and clear kmer list
                tstart = clock();
                std::cout << "Sorting kmers and appending succinct representation from current bin\t";
                construct::construct_succ(graph, config->parallel);
                std::cout << (clock()-tstart)/CLOCKS_PER_SEC << "\n";
            }
            //TODO: cleanup
            tstart = clock();
            std::cerr << "Converting static graph to dynamic\t";
            graph->toDynamic();
            construct::clean_bridges(graph);
            std::cout << (clock()-tstart)/CLOCKS_PER_SEC << "\n";
            //graph->print_seq();
        } break;

        case Config::annotate: {

            // load graph 
            if (config->infbase.empty()) {
              std::cerr << "Requires input <de bruijn graph> for annotation. Use option -I. " << std::endl;
              exit(1);
            }
            graph = new DBG_succ(config->infbase, config);

            // load annotatioun (if file does not exist, empty annotation is created)
            graph->annotationFromFile();

            // set up rocksdb
           // rocksdb::Options options;
           // options.create_if_missing = true;
           // rocksdb::DB* db;
           // rocksdb::Status status = rocksdb::DB::Open(options, config->dbpath, &db);
            
            uint64_t total_seqs = 0;
            
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

                    if (config->reverse)
                        reverse_complement(read_stream->seq);                    

                    if (config->parallel > 1) {
                        pthread_t* threads = NULL; 

                        anno_data->seq = &(read_stream->seq);
                        anno_data->label = &(read_stream->name);
                        anno_data->graph = graph;
                        anno_data->idx = 0;
                        anno_data->binsize = (read_stream->seq.l + 1) / config->parallel * config->bins_per_thread;
                        anno_data->total_bins = ((read_stream->seq.l + anno_data->binsize - 1) / anno_data->binsize); 

                        // create threads
                        threads = new pthread_t[config->parallel]; 
                        pthread_attr_init(&attr);
                        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

                        // do the work
                        for (size_t tid = 0; tid < config->parallel; tid++) {
                           pthread_create(&threads[tid], &attr, parallel_annotate_wrapper, (void *) tid); 
                           //std::cerr << "starting thread " << tid << std::endl;
                        }

                        // join threads
                        //if (config->verbose)
                        //    std::cout << "Waiting for threads to join" << std::endl;
                        for (size_t tid = 0; tid < config->parallel; tid++) {
                            pthread_join(threads[tid], NULL);
                        }
                        delete[] threads;

                        total_seqs += 1;

                        //if (config->verbose)
                        //    std::cout << "added labels for " << total_seqs << " sequences, last was " << std::string(read_stream->name.s) << std::endl;
                    } else {
                        annotate::annotate_seq(graph, read_stream->seq, read_stream->name);
                    }
                    if (config->verbose) {
                        std::cout << "entries in annotation map: " << graph->combination_count << std::endl << "length of combination vector: " << graph->combination_vector.size() << std::endl;
                        std::cout << "added labels for " << total_seqs << " sequences, last was " << std::string(read_stream->name.s) << std::endl;
                    }
                    /*for (std::unordered_map<uint32_t, std::set<uint32_t> >::iterator ittt = graph->annotation_map.begin(); ittt != graph->annotation_map.end(); ++ittt) {
                        std::cerr << "map : " << ittt->first << ":";
                        for (std::set<uint32_t>::iterator it4 = ittt->second.begin(); it4 != ittt->second.end(); ++it4) {
                            std::cerr << " " << *it4;
                        }
                        std::cerr << std::endl;
                    }*/
                }
                kseq_destroy(read_stream);
                gzclose(input_p);
            }

            if (config->print_graph)
                annotate::annotationToScreen(graph);

            graph->annotationToFile();

        } break;

        case Config::classify: {

            // load graph 
            if (config->infbase.empty()) {
              std::cerr << "Requires input <de bruijn graph> for annotation. Use option -I. " << std::endl;
              exit(1);
            }
            graph = new DBG_succ(config->infbase, config);

            // load annotatioun (if file does not exist, empty annotation is created)
            graph->annotationFromFile();

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
                
                std::set<uint32_t> labels_fwd;
                std::set<uint32_t> labels_rev;
                while (kseq_read(read_stream) >= 0) {

                    std::cout << std::string(read_stream->name.s) << "\t";
                    labels_fwd = annotate::classify_read(graph, read_stream->seq, config->distance);
                    for (std::set<uint32_t>::iterator it = labels_fwd.begin(); it != labels_fwd.end(); it++)
                        std::cout << graph->id_to_label.at(*it) << ":";
                    std::cout << "\t";

                    reverse_complement(read_stream->seq);                    
                    labels_rev = annotate::classify_read(graph, read_stream->seq, config->distance);
                    for (std::set<uint32_t>::iterator it = labels_rev.begin(); it != labels_rev.end(); it++)
                        std::cout << graph->id_to_label.at(*it) << ":";
                    std::cout << std::endl;
                }
                kseq_destroy(read_stream);
                gzclose(input_p);
            }
 


        } break;
    }

    // output and cleanup
    if (graph) {
        // graph output
        if (config->print_graph)
            graph->print_seq();
        if (!config->sqlfbase.empty())
            traverse::toSQL(graph);
        if (!config->outfbase.empty())
            graph->toFile(config->parts_total, config->part_idx);

        delete graph;
    }
    delete config;

    return 0;
}



