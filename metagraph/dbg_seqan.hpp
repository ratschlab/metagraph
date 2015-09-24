#ifndef __DBG_FULL_HPP__
#define __DBG_FULL_HPP__

/**
 * This class contains the uncompressed implementation of
 * the de bruijn graph using the Seqan graph class.
 */

#include <iostream>

#include <seqan/basic.h>
#include <seqan/graph_types.h>
#include <seqan/sequence.h>

using namespace seqan;

class DBG_seqan {

    public:
    /**
     * set up different seqan types that we will need to build the graph
     */
    typedef Value<Shape<Dna5, SimpleShape> >::Type THash;
    typedef Graph<Directed<void> > TGraph; // --> use this, when we don't need edge weights
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    Map<Pair<THash, TVertexDescriptor> > kmer_map;
    //typedef btree::btree_map<THash, TVertexDescriptor> TMap;
    //TMap *kmer_map = new TMap();
    Shape<Dna5, SimpleShape> hash_func;
    TGraph g;

    int k;
    unsigned int cnt_first_total;
    unsigned int cnt_recurr_total;
    unsigned int cnt_first;
    unsigned int cnt_recurr;

    /** 
     * construct graph object, based on given k-mer size
     */
    DBG_seqan(
        const int k
    ) : k(k) {
        resize(hash_func, k);
        cnt_first_total = 0;
        cnt_recurr_total = 0;
        cnt_first = 0;
        cnt_recurr = 0;
    }

    /**
     * Destructor
     */
    ~DBG_seqan() {};

    /**
     * This function takes a sequence string as input and adds all k-mers
     * of the sequence to the graph
     */
    void
    add_seq (
        String<Dna5> seq
    ) {

        // init hash function
        THash hash;
        hashInit(hash_func, begin(seq));

        TVertexDescriptor last_node;
        TVertexDescriptor current_node;
        //TVertexDescriptor *current_node = &addVertex(g);

        for (unsigned i = 0; i < length(seq) - length(hash_func) + 1; ++i) {
            if (i > 0 && i % 100000 == 0) {
                std::cout << "." << std::flush;
                if (i % 1000000 == 0)
                    fprintf(stdout, "%i\n", i);
            }
            hash = hashNext(hash_func, begin(seq) + i);
            //std::cout << hash << " " << infix(seq, i, i+4) << " ";
            if (!hasKey(kmer_map, hash)) {
            //if (kmer_map->find(hash) == kmer_map->end()) {
                current_node = addVertex(g);
                add(kmer_map, hash, current_node);
                //kmer_map->insert(std::make_pair(hash, current_node));
                ++cnt_first;
            } else {
                current_node = kmer_map[hash];
                //current_node = *(kmer_map->find(hash));
                ++cnt_recurr;
            }
            if (i > 1)
                addEdge(g, last_node, current_node);
            last_node = current_node;
        }
    };

    void update_counters()
    {
        cnt_first_total += cnt_first;
        cnt_recurr_total += cnt_recurr;
    };

    void print_stats() 
    {
        std::cout << std::endl << "k-mers opened new nodes: " << cnt_first << " / " << cnt_first_total << " (last / total)" << std::endl;
        std::cout << "k-mers already present in nodes: " << cnt_recurr << " / " << cnt_recurr_total << " (last / total)" << std::endl;
    };
};
#endif
