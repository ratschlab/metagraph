#include <fstream>
#include <ctime>

#include "dbg_succinct.hpp"
#include "config.hpp"
#include "helpers.hpp"
#include "utils.hpp"
#include "vcf_parser.hpp"
#include "traverse.hpp"
#include "dbg_succinct_merge.hpp"
#include "annotate.hpp"
#include "unix_tools.hpp"


KSEQ_INIT(gzFile, gzread)


/*
 * Helper function to determine the bin boundaries, given
 * a number of bins.
 */
std::vector<std::pair<uint64_t, uint64_t>> get_bins(DBG_succ *G, uint64_t bins) {

    uint64_t nodes = G->rank_last(G->get_W().size() - 1);
    uint64_t orig_bins = bins;
    std::cerr << "working with " << orig_bins << " orig bins; "
                                     << nodes << " nodes" <<  std::endl;
    if (bins > nodes) {
        std::cerr << "[WARNING] There are max "
                  << nodes << " slots available for binning. Your current choice is "
                  << bins << " which will create "
                  << bins - nodes << " empty slots." << std::endl;
        bins = nodes;
    }

    std::vector<std::pair<uint64_t, uint64_t>> result;
    uint64_t binsize = (nodes + bins - 1) / bins;
    uint64_t thresh = (nodes - (bins * (nodes / bins))) * binsize;
    uint64_t pos = 1;
    for (uint64_t i = 0; i < nodes;) {
        if (i >= thresh) {
            binsize = nodes / bins;
        }
        result.push_back(
            std::make_pair(pos, G->select_last(std::min(nodes, i + binsize)))
        );
        pos = result.back().second + 1;
        i += binsize;
    }

    for (uint64_t i = bins; i < orig_bins; i++) {
        //result.push_back(std::make_pair(pos, pos));
        result.push_back(std::make_pair(1, 0));
    }

    std::cerr << "created " << result.size() << " bins" << std::endl;
    return result;
}


std::vector<std::pair<uint64_t, uint64_t>> get_bins_relative(
                                                DBG_succ* G_from,
                                                DBG_succ* G_to,
                                                std::vector<std::pair<uint64_t, uint64_t>> ref_bins,
                                                uint64_t first_pos,
                                                uint64_t last_pos
                                            ) {

    std::vector<std::pair<uint64_t, uint64_t>> result;
    uint64_t pos = (first_pos == 0) ? 1 : G_from->colex_upper_bound(G_to->get_node_seq(first_pos)) + 1;
    uint64_t upper;
    for (size_t i = 0; i < ref_bins.size(); i++) {
        if (ref_bins.at(i).second == 0) { // this happens if we have more bins than nodes
            result.push_back(std::make_pair(0, 0));
        } else {
            upper = G_from->colex_upper_bound(G_to->get_node_seq(ref_bins.at(i).second));
            std::cerr << "ref bin " << ref_bins.at(i).second << " rel upper " << upper << std::endl;
            result.push_back(std::make_pair(pos, upper));
            pos = upper + 1;
        }
    }
    result.back().second = (last_pos == 0) ? G_from->get_W().size() - 1 : result.back().second;
    return result;
}


struct ParallelMergeContainer {
    std::vector<std::pair<uint64_t, uint64_t> > ref_bins;
    std::vector<std::vector<std::pair<uint64_t, uint64_t> > > bins;
    std::vector<DBG_succ*> result;
    std::vector<DBG_succ*> graphs;
    unsigned int idx;
    unsigned int k;
    unsigned int bins_done;
    unsigned int first = 0;
    unsigned int last = 0;

    /* Helper function to subset the bins to the chunk
     * computed in the current distributed compute.
     */
    void subset_bins(unsigned int idx, unsigned int total, unsigned int bins_per_part){

        //std::cerr << "ref bins " << ref_bins.size() << " total " << total << " per part " << bins_per_part << std::endl;
        assert(ref_bins.size() == (total * bins_per_part));

        std::vector< std::pair<uint64_t, uint64_t> > new_ref_bins;
        //std::cerr << "min: " << binsize_min << " max: " << binsize_max << " thresh: " << threshold << " total: " << total << std::endl;

        /*size_t start, end;
        if (idx < threshold) {
            start = binsize_max * idx;
            end = (idx == (total - 1)) ? ref_bins.size() : binsize_max * (idx + 1);
        } else {
            start = (threshold * binsize_max) + ((idx - threshold) * binsize_min);
            end = (idx == (total - 1)) ? ref_bins.size() : (threshold * binsize_max) + ((idx - threshold + 1) * binsize_min);
        }*/
        size_t start = idx * bins_per_part;
        size_t end = (idx + 1) * bins_per_part;

        if (start > 0)
            first = ref_bins.at(start - 1).second;
        if (end < ref_bins.size())
            last = ref_bins.at(end).second;

        for (size_t i = start; i < end; i++) {
            new_ref_bins.push_back(ref_bins.at(i));
        }

        ref_bins = new_ref_bins;
    }

    void print_bins() {
        for (size_t ii = 0; ii < bins.size(); ii++) {
            std::cerr << "graph " << ii + 1 << std::endl;
            for (size_t i = 0; i < bins.at(ii).size(); i++)
                std::cerr << bins.at(ii).at(i).first << " - " << bins.at(ii).at(i).second << std::endl;
        }
    }

    /* Show an overview of the distribution of merging bin
     * sizes.
     */
    void get_bin_stats() {
        size_t min_bin = 0, max_bin = 0, total_bin = 0;

        size_t cum_size;
        for (size_t i = 0; i < ref_bins.size(); ++i) {
            cum_size = 0;
            for (size_t ii = 0; ii < bins.size(); ii++) {
                cum_size += (bins.at(ii).at(i).first == 0) ? 0 : bins.at(ii).at(i).second - bins.at(ii).at(i).first + 1;
            }
            if (cum_size > 0) {
                min_bin = (min_bin == 0) ? cum_size : std::min(min_bin, cum_size);
                max_bin = (max_bin == 0) ? cum_size : std::max(max_bin, cum_size);
            }
            total_bin += cum_size;
        }

        std::cout << std::endl;
        std::cout << "Total number of bins: " << ref_bins.size() << std::endl;
        std::cout << "Total size: " << total_bin << std::endl;
        std::cout << "Smallest bin: " << min_bin << std::endl;
        std::cout << "Largest bin: " << max_bin << std::endl;
        std::cout << "Average bin size: " << total_bin / ref_bins.size() << std::endl << std::endl;
    }
};


struct ParallelAnnotateContainer {
    kstring_t *seq;
    kstring_t *label;
    DBG_succ *graph;
    Config *config;
    uint64_t idx;
    uint64_t binsize;
    uint64_t total_bins;
    //pthread_mutex_t* anno_mutex;
};


ParallelMergeContainer *merge_data;
// ParallelAnnotateContainer *anno_data = new ParallelAnnotateContainer();

pthread_mutex_t mutex_merge_result = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_bin_idx = PTHREAD_MUTEX_INITIALIZER;
// pthread_mutex_t mutex_annotate = PTHREAD_MUTEX_INITIALIZER;
pthread_attr_t attr;

const std::vector<std::string> annots = {
  "AC_AFR", "AC_EAS", "AC_AMR", "AC_ASJ",
  "AC_FIN", "AC_NFE", "AC_SAS", "AC_OTH"
};

void parallel_merge_collect(DBG_succ *result) {

    for (size_t i = 0; i < merge_data->result.size(); ++i) {
        if (merge_data->result.at(i)) {
            result->append_graph(merge_data->result.at(i));
            delete merge_data->result.at(i);
        }
    }
    merge_data->result.clear();
    merge_data->bins_done = 0;
    result->p_ = result->succ_W(1, 0);
}

/*
 * Distribute the annotation of all k-mers in a sequence over
 * a number of parallel bins.
 */
// void* parallel_annotate_wrapper(void *) {

//     uint64_t curr_idx, start, end;

//     while (true) {
//         pthread_mutex_lock (&mutex_bin_idx);
//         if (anno_data->idx == anno_data->total_bins) {
//             pthread_mutex_unlock (&mutex_bin_idx);
//             break;
//         } else {
//             curr_idx = anno_data->idx;
//             anno_data->idx++;
//             pthread_mutex_unlock (&mutex_bin_idx);

//             start = curr_idx * anno_data->binsize;
//             //end = std::min(((curr_idx + 1) * anno_data->binsize) + anno_data->graph->k - 1, anno_data->seq->l);
//             end = std::min((curr_idx + 1) * anno_data->binsize,
//                            static_cast<uint64_t>(anno_data->seq->l));
//             //std::cerr << "start " << start << " end " << end << std::endl;
//             annotate::annotate_seq(anno_data->graph,
//                                    anno_data->config,
//                                    *(anno_data->seq),
//                                    *(anno_data->label),
//                                    start, end, &mutex_annotate);
//         }
//     }
//     pthread_exit((void*) 0);
// }


/*
 * Distribute the merging of a set of graph structures over
 * bins, such that n parallel threads are used.
 */
void* parallel_merge_wrapper(void *config_) {
    Config *config = static_cast<Config *>(config_);

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

            graph = new DBG_succ(merge_data->k, false);
            std::vector<uint64_t> kv;
            std::vector<uint64_t> nv;
            for (size_t i = 0; i < merge_data->graphs.size(); i++) {
                kv.push_back(merge_data->bins.at(i).at(curr_idx).first);
                nv.push_back(merge_data->bins.at(i).at(curr_idx).second + 1);
            }
            merge::merge(graph, merge_data->graphs, kv, nv);

            pthread_mutex_lock (&mutex_merge_result);
            merge_data->result.at(curr_idx) = graph;
            merge_data->bins_done++;
            if (config->verbose)
                std::cout << "finished bin " << curr_idx + 1
                          << " (" << merge_data->bins_done
                                  << "/" << merge_data->ref_bins.size()
                          << ")" << std::endl;
            pthread_mutex_unlock (&mutex_merge_result);
        }
    }
    pthread_exit((void*) 0);
}


int main(int argc, const char *argv[]) {

    // parse command line arguments and options
    Config *config = new Config(argc, argv);

    if (config->verbose)
        std::cout << "Welcome to MetaGraph" << std::endl;

    // create graph pointer
    DBG_succ *graph = NULL;

    const auto &files = config->fname;

    switch (config->identity) {
        //TODO: allow for building by appending/adding to an existing graph
        case Config::BUILD: {
            if (config->infbase.size()) {
                graph = new DBG_succ();
                graph->load(config->infbase);
                // graph->annotationFromFile(config->infbase + ".anno.dbg");
            } else {
                graph = new DBG_succ(config->k);
            }

            if (config->verbose)
                std::cerr << "k is " << graph->k_ << std::endl;

            if (config->fast) {
                graph->switch_state(Config::CSTR);

                //enumerate all suffices
                std::deque<std::string> suffices = graph->generate_suffices(config->nsplits);

                clock_t tstart, timelast;

                //one pass per suffix
                for (size_t j = 0; j < suffices.size(); ++j) {
                    std::cout << "Suffix: " << suffices[j] << "\n";
                    //add sink nodes
                    graph->add_sink(config->parallel, suffices[j]);

                    if (suffices[j].find("$") == std::string::npos) {
                        // iterate over input files
                        for (unsigned int f = 0; f < files.size(); ++f) {
                            if (config->verbose) {
                                std::cout << std::endl << "Parsing " << files[f] << std::endl;
                            }
                            // open stream
                            gzFile input_p = gzopen(files[f].c_str(), "r");

                            if (utils::get_filetype(files[f]) == "VCF") {
                                //READ FROM VCF
                                uint64_t nbp = 0;
                                uint64_t nbplast = 0;
                                tstart = clock();
                                timelast = clock();
                                vcf_parser vcf;
                                if (!vcf.init(config->refpath, files[f], graph->k_)) {
                                    std::cerr << "ERROR reading VCF " << files[f] << std::endl;
                                    exit(1);
                                }
                                std::cerr << "Loading VCF with " << config->parallel << " threads per line\n";
                                std::string sequence;
                                std::string annotation;
                                for (size_t i = 1; vcf.get_seq(annots, &sequence, &annotation); ++i) {
                                    if (i % 10'000 == 0) {
                                        std::cout << "." << std::flush;
                                        if (i % 100'000 == 0) {
                                            fprintf(stdout, "%zu - bp %" PRIu64 " / runtime %lu / BPph %" PRIu64 "\n",
                                                            i,
                                                            nbp,
                                                            (clock() - tstart) / CLOCKS_PER_SEC,
                                                            uint64_t(60) * uint64_t(60)
                                                                * CLOCKS_PER_SEC * (nbp - nbplast)
                                                                / (clock() - timelast));
                                            nbplast = nbp;
                                            timelast = clock();
                                        }
                                    }
                                    annotation = "VCF:" + annotation;
                                    nbp += sequence.length();
                                    graph->add_sequence_fast(sequence,
                                                             false, config->parallel,
                                                             suffices[j]);
                                }
                            } else {
                                //READ FROM FASTA
                                //TODO: handle read_stream->qual
                                kseq_t *read_stream = kseq_init(input_p);
                                if (read_stream == NULL) {
                                    std::cerr << "ERROR while opening input file " << files[f] << std::endl;
                                    exit(1);
                                }
                                for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
                                //while (kseq_read(read_stream) >= 0) {
                                    // possibly reverse k-mers
                                    if (config->reverse)
                                        reverse_complement(read_stream->seq);
                                    // add all k-mers of seq to the graph
                                    graph->add_sequence_fast(std::string(read_stream->seq.s, read_stream->seq.l),
                                                             true, config->parallel, suffices[j]);
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
                    graph->construct_succ(config->parallel);
                    std::cout << (clock()-tstart)/CLOCKS_PER_SEC << "\n\n";
                }
                //TODO: cleanup
                tstart = clock();
                std::cerr << "Converting static graph to dynamic\t";
                graph->switch_state(Config::DYN);
                std::cout << (clock()-tstart)/CLOCKS_PER_SEC << "\n";
            } else {
                //slower method
                //TODO: merge in optimizations from seqmerge branch
                //TODO: add VCF parsing (needs previous merge)
                for (unsigned int f = 0; f < files.size(); ++f) {
                    if (config->verbose) {
                        std::cout << std::endl << "Parsing " << files[f] << std::endl;
                    }
                    // open stream
                    gzFile input_p = gzopen(files[f].c_str(), "r");
                    if (utils::get_filetype(files[f]) == "VCF") {
                        std::cerr << "ERROR: this method of reading VCFs not yet implemented" << std::endl;
                        exit(1);
                    } else {
                        kseq_t *read_stream = kseq_init(input_p);
                        if (read_stream == NULL) {
                            std::cerr << "ERROR while opening input file " << files[f] << std::endl;
                            exit(1);
                        }
                        for (size_t i = 1; kseq_read(read_stream) >= 0; ++i) {
                            if (config->reverse)
                                reverse_complement(read_stream->seq);
                            graph->add_sequence(std::string(read_stream->seq.s, read_stream->seq.l));
                        }
                        kseq_destroy(read_stream);
                    }
                    gzclose(input_p);
                }
            }
            graph->switch_state(config->state);
            config->infbase = config->outfbase;
            // graph->annotationToFile(config->infbase + ".anno.dbg");
            //graph->print_seq();
            break;
        }
        case Config::ANNOTATE: {
           //  // load graph
           //  if (config->infbase.empty()) {
           //    std::cerr << "Requires input <de bruijn graph> for annotation. Use option -I. " << std::endl;
           //    exit(1);
           //  }
           //  graph = new DBG_succ(config->infbase);

           //  // load annotation (if file does not exist, empty annotation is created)
           //  // graph->annotationFromFile(config->infbase + ".anno.dbg");

           //  // set up rocksdb
           // // rocksdb::Options options;
           // // options.create_if_missing = true;
           // // rocksdb::DB* db;
           // // rocksdb::Status status = rocksdb::DB::Open(options, config->dbpath, &db);

           //  uint64_t total_seqs = 0;

           //  // iterate over input files
           //  for (unsigned int f = 0; f < files.size(); ++f) {

           //      if (config->verbose) {
           //          std::cout << std::endl << "Parsing " << files[f] << std::endl;
           //      }

           //      // open stream to fasta file
           //      gzFile input_p = gzopen(files[f].c_str(), "r");
           //      kseq_t *read_stream = kseq_init(input_p);

           //      if (read_stream == NULL) {
           //          std::cerr << "ERROR while opening input file " << files[f] << std::endl;
           //          exit(1);
           //      }

           //      while (kseq_read(read_stream) >= 0) {

           //          if (config->reverse)
           //              reverse_complement(read_stream->seq);

           //          if (config->parallel > 1) {
           //              pthread_t* threads = NULL;

           //              anno_data->seq = &(read_stream->seq);
           //              anno_data->label = &(read_stream->name);
           //              anno_data->graph = graph;
           //              anno_data->config = config;
           //              anno_data->idx = 0;
           //              anno_data->binsize = (read_stream->seq.l + 1) / (config->parallel * config->bins_per_thread);
           //              anno_data->total_bins = ((read_stream->seq.l + anno_data->binsize - 1) / anno_data->binsize);

           //              // create threads
           //              threads = new pthread_t[config->parallel];
           //              pthread_attr_init(&attr);
           //              pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

           //              // do the work
           //              for (size_t tid = 0; tid < config->parallel; tid++) {
           //                 pthread_create(&threads[tid], &attr, parallel_annotate_wrapper, (void *) tid);
           //                 //std::cerr << "starting thread " << tid << std::endl;
           //              }

           //              // join threads
           //              //if (config->verbose)
           //              //    std::cout << "Waiting for threads to join" << std::endl;
           //              for (size_t tid = 0; tid < config->parallel; tid++) {
           //                  pthread_join(threads[tid], NULL);
           //              }
           //              delete[] threads;

           //              total_seqs += 1;

           //              //if (config->verbose)
           //              //    std::cout << "added labels for " << total_seqs
           //              //              << " sequences, last was " << std::string(read_stream->name.s) << std::endl;
           //          } else {
           //              annotate::annotate_seq(graph, config, read_stream->seq, read_stream->name);
           //          }
           //          if (config->verbose) {
           //              //std::cout << "entries in annotation map: " << graph->combination_count << std::endl
           //              //          << "length of combination vector: " << graph->combination_vector.size() << std::endl;
           //              //std::cout << "added labels for " << total_seqs << " sequences, last was " << std::string(read_stream->name.s) << std::endl;
           //          }
           //      }
           //      kseq_destroy(read_stream);
           //      gzclose(input_p);
           //  }

           //  if (config->print_graph)
           //      annotate::annotationToScreen(graph);

           //  graph->annotationToFile(config->infbase + ".anno.dbg");
            break;
        }
        case Config::CLASSIFY: {
            // load graph
            if (config->infbase.empty()) {
              std::cerr << "Requires input <de bruijn graph> for annotation. Use option -I. " << std::endl;
              exit(1);
            }
            graph = new DBG_succ();
            graph->load(config->infbase);

            // load annotatioun (if file does not exist, empty annotation is created)
            // graph->annotationFromFile(config->infbase + ".anno.dbg");

            // iterate over input files
            for (unsigned int f = 0; f < files.size(); ++f) {

                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << files[f] << std::endl;
                }

                // open stream to fasta file
                gzFile input_p = gzopen(files[f].c_str(), "r");
                kseq_t *read_stream = kseq_init(input_p);

                if (read_stream == NULL) {
                    std::cerr << "ERROR while opening input file " << files[f] << std::endl;
                    exit(1);
                }

                std::set<uint32_t> labels_fwd;
                std::set<uint32_t> labels_rev;
                while (kseq_read(read_stream) >= 0) {

                    std::cout << std::string(read_stream->name.s) << ": ";
                    // labels_fwd = annotate::classify_read(graph, read_stream->seq, config->distance);
                    // for (auto it = labels_fwd.begin(); it != labels_fwd.end(); it++)
                    //     std::cout << graph->id_to_label.at(*it + 1) << ":";
                    // std::cout << "\t";

                    reverse_complement(read_stream->seq);
                    // labels_rev = annotate::classify_read(graph, read_stream->seq, config->distance);
                    // for (auto it = labels_rev.begin(); it != labels_rev.end(); it++)
                    //     std::cout << graph->id_to_label.at(*it + 1) << ":";
                    // std::cout << std::endl;
                }
                kseq_destroy(read_stream);
                gzclose(input_p);
            }
            break;
        }
        case Config::COMPARE: {
            for (size_t f = 0; f < files.size(); ++f) {
                if (f == 0) {
                    std::cout << "Opening file " << files[f] << std::endl;
                    graph = new DBG_succ();
                    graph->load(files[f]);
                } else {
                    std::cout << "Opening file for comparison ..." << files[f] << std::endl;
                    DBG_succ *graph_ = new DBG_succ();
                    graph_->load(files[f]);
                    if (*graph == *graph_) {
                        std::cout << "Graphs are identical" << std::endl;
                    } else {
                        std::cout << "Graphs are not identical" << std::endl;
                    }
                    delete graph_;
                }
            }
            break;
        }
        case Config::MERGE: {
            // collect results on an external merge
            if (config->collect > 1) {
                std::string fname;
                for (uint64_t f = 0; f < config->collect; f++) {
                    fname = config->outfbase + "." + std::to_string(f) + "_" + std::to_string(config->collect);
                    std::cout << "Opening file " << fname << std::endl;
                    if (config->fast) {
                        if (f == 0) {
                            graph = new DBG_succ(utils::kFromFile(fname), false);
                            graph->last_stat.push_back(0);
                            graph->W_stat.push_back(0);
                        }
                        DBG_succ* graph_to_append = new DBG_succ();
                        graph_to_append->load(fname);
                        graph->append_graph_static(graph_to_append);
                        delete graph_to_append;
                    } else {
                        if (f == 0) {
                            graph = new DBG_succ();
                            graph->load(fname);
                        } else {
                            DBG_succ* graph_to_append = new DBG_succ();
                            graph_to_append->load(fname);
                            graph->append_graph(graph_to_append);
                            delete graph_to_append;
                        }
                    }
                }
                if (config->fast) {
                    //graph->toDynamic();
                    graph->switch_state(Config::DYN);
                }
                graph->p_ = graph->succ_W(1, 0);

            // run normal merge procedure
            } else {

                // some preliminaries to make command line options consistent
                // if ((config->parts_total > 1) && (config->parts_total > (config->parallel * config->bins_per_thread)))
                //    config->bins_per_thread = config->parts_total / config->parallel;

                std::vector<DBG_succ*> graphs;
                std::vector<uint64_t> kv;
                std::vector<uint64_t> nv;
                for (unsigned int f = 0; f < files.size(); ++f) {
                        std::cout << "Opening file " << files[f] << std::endl;
                        graph = new DBG_succ();
                        graph->load(files[f]);
                        graphs.push_back(graph);
                        kv.push_back(1);
                        nv.push_back(graph->get_W().size());
                }
                DBG_succ* target_graph = new DBG_succ(graphs.front()->get_k(), false);

                if ((config->parallel > 1) || (config->parts_total > 1)) {
                    pthread_t *threads = NULL;
                    merge_data = new ParallelMergeContainer();

                    // get bins in graphs according to required threads
                    if (config->verbose)
                        std::cout << "Collecting reference bins" << std::endl;
                    std::cerr << "parallel " << config->parallel
                              << " per thread " << config->bins_per_thread
                              << " parts total " << config->parts_total << std::endl;
                    merge_data->ref_bins = get_bins(graphs.front(), config->parallel * config->bins_per_thread * config->parts_total);

                    // only work on subset of the bins when requested
                    if (config->parts_total > 1) {
                        merge_data->subset_bins(config->part_idx, config->parts_total, config->parallel * config->bins_per_thread);
                    }
                    merge_data->bins.push_back(merge_data->ref_bins);

                    if (config->verbose)
                        std::cout << "Collecting relative bins" << std::endl;
                    for (size_t i = 1; i < graphs.size(); i++) {
                        std::cerr << "getting bins for " << i << ": " << files[i] << std::endl;
                        merge_data->bins.push_back(get_bins_relative(graphs.at(i), graphs.front(), merge_data->ref_bins,
                                                                     merge_data->first, merge_data->last));
                    }
                    for (size_t i = 0; i < graphs.size(); i++) {
                        for (size_t ii = 0; ii < merge_data->bins.at(i).size(); ii++) {
                            if (merge_data->bins.at(i).at(ii).first > merge_data->bins.at(i).at(ii).second) {
                               merge_data->bins.at(i).at(ii) = std::make_pair(graphs.at(i)->get_W().size(),
                                                                              graphs.at(i)->get_W().size() - 1);
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
                       pthread_create(&threads[tid], &attr, parallel_merge_wrapper, (void *) config);
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
            break;
        }
        case Config::STATS: {
            std::ofstream outstream;
            if (!config->outfbase.empty()) {
                outstream.open((config->outfbase + ".stats.dbg").c_str());
                outstream << "file\tnodes\tedges\tk" << std::endl;
            }
            for (unsigned int f = 0; f < files.size(); ++f) {
                DBG_succ* graph_ = new DBG_succ();
                graph_->load(files[f]);
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
                    std::cout << "Statistics for file " << files[f] << std::endl;
                    std::cout << "nodes: " << graph_->num_nodes() << std::endl;
                    std::cout << "edges: " << graph_->num_edges() << std::endl;
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
                    outstream << files[f] << "\t"
                              << graph_->num_nodes() << "\t"
                              << graph_->num_edges() << "\t"
                              << graph_->get_k() << std::endl;
                }
                if (config->print_graph)
                    graph_->print_seq();
                if (config->print_graph_succ)
                    graph_->print_state();

                std::ifstream instream((files[f] + ".anno.dbg").c_str());
                if (instream.good()) {
                    size_t anno_size = libmaus2::util::NumberSerialisation::deserialiseNumber(instream);
                    std::cout << "annot: " << anno_size << std::endl;
                }
                instream.close();

                /*DBG_succ* graph_tut = new DBG_succ(config->k);
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

                /*std::cerr << graph_->num_edges() << std::endl;
                std::vector<uint64_t> result = graph_->split_range(1, graph_->num_edges(), 0);
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
            break;
        }
        case Config::DUMP: {
            //for (unsigned int f = 0; f < files.size(); ++f) {
                //DBG_succ* graph_ = new DBG_succ(files[f]);
                DBG_succ *graph_ = new DBG_succ();
                graph_->load(config->infbase);
                // checks whether annotation exists and creates an empty one if not
                // graph_->annotationFromFile(config->infbase + ".anno.dbg");
                graph_->print_adj_list();
                //graph_->print_adj_list(config->outfbase);
                delete graph_;
            //}
            break;
        }
        case Config::ALIGN: {
            // load graph
            if (config->infbase.empty()) {
              std::cerr << "Requires input <de bruijn graph> to align reads." << config->ALIGN << std::endl;
              exit(1);
            }
            graph = new DBG_succ();
            graph->load(config->infbase);

            for (unsigned int f = 0; f < files.size(); ++f) {
                std::cout << "Opening file for alignment ..." << files[f] << std::endl;

                // open stream to input fasta
                gzFile input_p = gzopen(files[f].c_str(), "r");
                kseq_t *read_stream = kseq_init(input_p);
                if (read_stream == NULL) {
                  std::cerr << "ERROR while opening input file " << config->ALIGN << std::endl;
                  exit(1);
                }

                while (kseq_read(read_stream) >= 0) {

                    //graph->print_seq();
                    uint64_t aln_len = read_stream->seq.l;
                    //bool reverse = false;

                    if (config->distance > 0) {
                        // since we set aln_len = read_stream->seq.l, we only get a single hit vector
                        auto graphindices = graph->align_fuzzy(
                            std::string(read_stream->seq.s, read_stream->seq.l),
                            aln_len,
                            config->distance
                        );
                        //for (size_t i = 0; i < graphindices.size(); ++i) {
                        size_t i = 0;
                        int print_len = (i + aln_len < read_stream->seq.l) ? aln_len : (read_stream->seq.l - i);
                        printf("%.*s: ", print_len, read_stream->seq.s + i);

                        for (size_t l = 0;  l < graphindices.at(i).size(); ++l) {
                            HitInfo curr_hit(graphindices.at(i).at(l));
                            for (size_t m = 0; m < curr_hit.path.size(); ++m) {
                                std::cout << curr_hit.path.at(m) << ':';
                            }
                            for (size_t m = curr_hit.rl; m <= curr_hit.ru; ++m) {
                                std::cout << m << " ";
                            }
                            std::cout << "[" << curr_hit.cigar << "] ";
                        }
                        //}

                        // try reverse
                        if (graphindices.at(i).size() == 0) {
                            reverse_complement(read_stream->seq);
                            graphindices = graph->align_fuzzy(
                                std::string(read_stream->seq.s, read_stream->seq.l),
                                aln_len, config->distance
                            );
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
                        }
                        std::cout << std::endl;
                    } else {
                        std::vector<uint64_t> graphindices = graph->align(
                            std::string(read_stream->seq.s, read_stream->seq.l)
                        );

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
            break;
        }
        case Config::NO_IDENTITY: {
            assert(false);
            break;
        }
    }

    // output and cleanup
    if (graph) {
        // graph output
        if (config->print_graph)
            graph->print_seq();
        if (!config->sqlfbase.empty())
            traverse::toSQL(graph, config->fname, config->sqlfbase);
        if (!config->outfbase.empty())
            graph->serialize(config->outfbase
                                + "." + std::to_string(config->part_idx)
                                + "_" + std::to_string(config->parts_total));
        delete graph;
    }
    delete config;

    return 0;
}
