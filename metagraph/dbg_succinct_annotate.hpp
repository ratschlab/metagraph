#ifndef __DBG_BLOOM_ANNOTATOR_HPP__
#define __DBG_BLOOM_ANNOTATOR_HPP__

#include <unordered_set>

#include "dbg_succinct.hpp"
#include "annotate.hpp"
#include "kmer.hpp"


//TODO
class DeBruijnGraphWrapper {
  public:
    DeBruijnGraphWrapper();


};

class Annotator {
  public:
    Annotator(DBG_succ *graph, double bloom_size_factor)
        : graph_(graph),
          bloom_size_factor_(bloom_size_factor),
          total_traversed_(0) {}

    void init_exact_hasher() {
        assert(!annotation_exact);
        annotation_exact = new annotate::HashAnnotation<annotate::ExactFilter>();
    }

    ~Annotator() {
        if (annotation_exact)
            delete annotation_exact;
    }

    void set_graph(DBG_succ *graph) {
        graph_ = graph;
    }

    bool exact_enabled() {
        return annotation_exact != NULL;
    }

    void add_sequences(const std::vector<KMer>::iterator &begin,
                       const std::vector<KMer>::iterator &end,
                       size_t column) {
        if (column >= annotation.size()) {
            annotation.resize(column);
            sizes_v.emplace_back(end - begin);
            annotation.append_bit((size_t)((double)sizes_v.back() * bloom_size_factor_));
        }
        if (annotation[column].size() == 0) {
            annotation[column].resize((size_t)((double)(end - begin) * bloom_size_factor_));
        }
        for (auto it = begin; it != end; ++it) {
            annotation.insert(it->begin(), it->end(), column);
            if (annotation_exact)
                annotation_exact->insert(it->begin(), it->end(), column);
        }
    }

    template <typename T>
    void add_sequences(const std::vector<KMer>::iterator &begin,
                       const std::vector<KMer>::iterator &end,
                       T *cat_begin, T *cat_end) {
        for (auto it = begin; it != end; ++it) {
            annotation.insert(it->begin(), it->end(), cat_begin, cat_end);
            if (annotation_exact)
                annotation_exact->insert(it->begin(), it->end(), cat_begin, cat_end);
        }
    }

    void add_column(const std::vector<KMer>::iterator &begin,
                      const std::vector<KMer>::iterator &end) {
        sizes_v.emplace_back(end - begin);
        annotation.append_bit((size_t)((double)sizes_v.back() * bloom_size_factor_));
        if (annotation_exact)
            annotation_exact->append_bit();
        add_sequences(begin, end, annotation.size() - 1);
    }

    void set_junction(int size = -1) {
        if (size == -1) {
            //size = *(std::max_element(sizes_v.begin(), sizes_v.end()));
            size = std::accumulate(sizes_v.begin(), sizes_v.end(), 0);
        }
        annotation.set_cont(size);

        //compute junction bits
        for (size_t i = 1; i < graph_->get_W().size(); ++i) {
            //outdegree > 1
            if (!graph_->get_last(i) && graph_->get_last(i - 1)) {
                auto kmer = graph_->get_node_seq(i);
                kmer.emplace_back(0);
                //for every edge with last == 0
                auto j = i;
                for (; !graph_->get_last(j); ++j) {
                    kmer.back() = graph_->get_W(j);
                    KMer int_kmer(KMer::pack_kmer(kmer.begin(), kmer.size()));
                    annotation.set_junction(int_kmer.begin(), int_kmer.end());
                }
                //for last outgoing edge
                kmer.back() = graph_->get_W(j);
                KMer int_kmer(KMer::pack_kmer(kmer.begin(), kmer.size()));
                annotation.set_junction(int_kmer.begin(), int_kmer.end());
            }
            //indegree > 1
            if (graph_->get_W(i) >= DBG_succ::alph_size) {
                size_t j = graph_->outgoing(i, graph_->get_W(i));
                auto kmer = graph_->get_node_seq(j);
                kmer.emplace_back(0);
                for (; !graph_->get_last(j); ++j) {
                    kmer.back() = graph_->get_W(j);
                    KMer int_kmer(KMer::pack_kmer(kmer.begin(), kmer.size()));
                    annotation.set_junction(int_kmer.begin(), int_kmer.end());
                }
                kmer.back() = graph_->get_W(j);
                KMer int_kmer(KMer::pack_kmer(kmer.begin(), kmer.size()));
                annotation.set_junction(int_kmer.begin(), int_kmer.end());
            }
        }
    }

    template <typename KmerVector>
    std::vector<size_t> annotation_from_seq(KmerVector kvector) {
        KMer kmer(KMer::pack_kmer(kvector, kvector.size()));
        return annotation.find(kmer.begin(), kmer.end());
    }

    std::vector<size_t> annotation_from_index(size_t i) {
        auto seq = graph_->get_node_seq(i);
        return annotation_from_seq(seq);
    }

    std::vector<size_t> annotation_from_kmer(KMer &kmer) {
        return annotation.find(kmer.begin(), kmer.end());
    }

    std::vector<size_t> annotation_exact_from_kmer(KMer &kmer) {
        if (!annotation_exact) {
            std::cerr << "ERROR: exact hasher not initialized.\n";
            exit(1);
        }
        return annotation_exact->find(kmer.begin(), kmer.end());
    }

    KMer kmer_from_index(const DBG_succ::node_iterator &i) {
        auto kmer = graph_->get_node_seq(i);
        kmer.emplace_back(graph_->get_W(i) % DBG_succ::alph_size);
        return KMer(KMer::pack_kmer(kmer, kmer.size()));
    }

    std::vector<size_t> get_annotation(const DBG_succ::node_iterator &i) {
        KMer int_kmer = kmer_from_index(i);
        return annotation_from_kmer(int_kmer);
    }


    std::vector<size_t> get_annotation_corrected(const DBG_succ::node_iterator &i,
                                                 size_t path_cutoff = 50) {
        //initial raw annotation
        auto orig_kmer = graph_->get_node_seq(i);
        orig_kmer.emplace_back(graph_->get_W(i) % DBG_succ::alph_size);
        auto ckmer = orig_kmer;

        auto curannot = annotation_from_seq(ckmer);

        size_t pcount_old = annotate::HashAnnotation<>::bigint_popcount(curannot);
        size_t pcount_new = 0;
        auto nextannot = curannot;
        size_t path = 0;
        DBG_succ::node_iterator j = i;

        //TODO: checking for ckmer.back() and ckmer.front() means that no
        // dummy edges can be annotation. Either find a better heuristic, or find
        // a way to label them.

        while (ckmer.back() && pcount_old && path++ < path_cutoff) {
            //move forward
            j = graph_->outgoing(j, ckmer.back());

            ckmer.emplace_back(graph_->get_W(j) % DBG_succ::alph_size);
            ckmer.pop_front();

            total_traversed_++;
            //update_kmer(graph_->get_k(), nextW, curW, &ckmer);

            //check outdegree
            if (!graph_->get_last(j) || !graph_->get_last(j - 1))
                break;

            //curW = nextW;

            //bitwise AND annotations
            nextannot = annotate::HashAnnotation<>::merge_and(
                    curannot, 
                    //annotation.find(ckmer.begin(), ckmer.end())
                    annotation_from_seq(ckmer)
            );

            //check popcounts
            pcount_new = annotate::HashAnnotation<>::bigint_popcount(nextannot);
            assert(pcount_new <= pcount_old);

            //if zero, then start of new sequence
            if (pcount_new == 0) {
                break;
            }

            //the new annotatation is "correct", so accept it
            curannot = nextannot;

            //path length stopping conditions
            if (pcount_new < pcount_old) {
                path = 0;
                pcount_old = pcount_new;
            }
            //std::cout << pcount_old << " " << path << "\n";
        }

        //backward correction
        ckmer = orig_kmer;
        j = i;
        path = 0;
        while (ckmer.front()
                && graph_->indegree(j) == 1
                && pcount_old
                && path++ < path_cutoff) {
            j = graph_->get_minus_k_value(j, 0).second;

            ckmer.emplace_front(
                graph_->get_minus_k_value(j, graph_->get_k() - 1).first
            );
            ckmer.pop_back();

            total_traversed_++;
            nextannot = annotate::HashAnnotation<>::merge_and(
                    curannot, 
                    //annotation.find(ckmer.begin(), ckmer.end())
                    annotation_from_seq(ckmer)
            );

            pcount_new = annotate::HashAnnotation<>::bigint_popcount(nextannot);
            assert(pcount_new <= pcount_old);

            if (pcount_new == 0) {
                break;
            }
            curannot = nextannot;

            if (pcount_new < pcount_old) {
                path = 0;
                pcount_old = pcount_new;
            }
        }


		return curannot;
    }

    std::vector<uint8_t> test_fp(const DBG_succ::node_iterator &i) {
        if (!annotation_exact) {
            std::cerr << "ERROR: exact hasher not initialized.\n";
            exit(1);
        }
        KMer int_kmer = kmer_from_index(i);
        auto test = annotation_from_kmer(int_kmer);
        auto test_exact = annotation_exact_from_kmer(int_kmer);
        auto curannot = get_annotation_corrected(i);
        auto jt = test.begin();
        auto kt = test_exact.begin();
        auto lt = curannot.begin();
        std::vector<uint8_t> stats(3);
        for (; jt != test.end(); ++jt, ++kt, ++lt) {
            //check for false negatives
            if ((*jt | *kt) != *jt) {
                std::cerr << "Encoding " << i << " failed\n";
                exit(1);
            }
            //correction introduced extra bits
            if ((*lt | *jt) != *jt) {
                std::cerr << "False bits added\n";
                exit(1);
            }
            //false positives before correction
            if (!stats[0] && (*jt | *kt) != *kt) {
                stats[0] = 1;
            }
            //false positives after correction
            if (!stats[1] && (*lt | *kt) != *kt) {
                stats[1] = 1;
            }
            //false negatives after correction
            if (!stats[2] && (*lt | *kt) != *lt) {
                stats[2] = 1;
            }
            if (stats[0] && stats[1] && stats[2]) {
                break;
            }
        }
        return stats;
    }

    void test_fp_all(size_t step = 1) {
        size_t fp = 0;
        size_t fp_pre = 0;
        size_t fn = 0;
        //size_t total = graph_->get_W().size() / step;
        size_t total = 0;
        for (DBG_succ::node_iterator i = 1; i < graph_->get_W().size(); i += step) {
            total++;
            auto stats = test_fp(i);
            fp_pre += stats[0];
            fp += stats[1];
            fn += stats[2];
        }
        std::cout << "\n";
        std::cout << "Total:\t" << total << "\n";
        std::cout << "Post:\t"
                  << "FP:\t" << fp << " "
                  << (double)fp / (double)total << "\t"
                  << "FN:\t" << fn
                  << "\n";
        std::cout << "Pre:\t"
                  << "FP:\t" << fp_pre << " "
                  << (double)fp_pre / (double)total << "\t"
                  << "\n";
        std::cout << "Total traversed: " << total_traversed_ << "\n";
    }

    void serialize(std::ostream &out) {
        annotation.serialize(out);
    }

    void serialize(const std::string &filename) {
        std::ofstream out(filename + ".annot.dbg");
        serialize(out);
        out.close();
    }

  private:
    DBG_succ *graph_;
    double bloom_size_factor_;

    //TODO: get rid of this if not using degree Bloom filter
    std::vector<size_t> sizes_v;

    size_t total_traversed_ = 0;

    //TODO: set this value
    annotate::HashAnnotation<annotate::BloomFilter<1>> annotation;
    annotate::HashAnnotation<annotate::ExactFilter> *annotation_exact = NULL;
};

#endif // __DBG_BLOOM_ANNOTATOR_HPP__
