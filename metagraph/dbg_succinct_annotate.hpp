#include "dbg_succinct.hpp"
#include "annotate.hpp"
#include "kmer.hpp"
#include <unordered_set>

class Annotator {
    public:

        Annotator(DBG_succ *graph, double&& bloom_size_factor) 
            : graph(graph), bloom_size_factor(bloom_size_factor) { }
        Annotator(DBG_succ *graph, const double &bloom_size_factor) 
            : graph(graph), bloom_size_factor(bloom_size_factor) { }

        void set_graph(DBG_succ *graph) {
            this->graph = graph;
        }

        void add_sequences(const std::vector<KMer>::iterator &begin,
                           const std::vector<KMer>::iterator &end,
                           const size_t &category) {
            for (auto it = begin; it != end; ++it) {
                annotation.insert(it->begin(), it->end(), category);
                annotation_exact.insert(it->begin(), it->end(), category);
            }
        }

        template <typename T>
        void add_sequences(const std::vector<KMer>::iterator &begin,
                           const std::vector<KMer>::iterator &end,
                           T *cat_begin, T *cat_end) {
            for (auto it = begin; it != end; ++it) {
                annotation.insert(it->begin(), it->end(), cat_begin, cat_end);
                annotation_exact.insert(it->begin(), it->end(), cat_begin, cat_end);
            }
        }

        void add_category(const std::vector<KMer>::iterator &begin,
                          const std::vector<KMer>::iterator &end) {
            sizes_v.emplace_back(end - begin);
            annotation.append_bit((size_t)((double)sizes_v.back() * bloom_size_factor));
            annotation_exact.append_bit();
            add_sequences(begin, end, annotation.size() - 1);
        }

        void set_junction(int size = -1) {
            if (size == -1) {
                size = *(std::max_element(sizes_v.begin(), sizes_v.end()));
            }
            annotation.set_cont(size);
            //compute using succ0_last
            for (size_t i = 1; i < graph->get_W().size(); i = graph->succ0_last(i)) {
                auto kmer = graph->get_node_seq(i);
                kmer.emplace_back(0);
                size_t next = graph->succ_last(i);
                for (; i <= next; ++i) {
                    kmer.back() = graph->get_W(i) % DBG_succ::alph_size;
                    KMer int_kmer(KMer::pack_kmer(kmer.begin(), kmer.size()));
                    annotation.set_junction(int_kmer.begin(), int_kmer.end());
                }
            }
            //compute using succ_last
            /*
            size_t last = 1;
            size_t cur;
            while (last < graph->get_W().size()) {
                cur = graph->succ_last(last + 1);
                if (cur - last > 1) {
                    auto kmer = graph->get_node_seq(cur);
                    kmer.push_back(0);
                    for (auto it = last + 1; it <= cur; ++it) {
                        kmer.back() = graph->get_W(it) % DBG_succ::alph_size;
                        KMer int_kmer(KMer::pack_kmer(kmer.begin(), kmer.size()));
                        annotation.set_junction(int_kmer.begin(), int_kmer.end());
                    }
                } else {
                    if (graph->indegree(cur) > 1) {
                        auto kmer = graph->get_node_seq(cur);
                        kmer.emplace_back(graph->get_W(cur) % DBG_succ::alph_size);
                        KMer int_kmer(KMer::pack_kmer(kmer.begin(), kmer.size()));
                        annotation.set_junction(int_kmer.begin(), int_kmer.end());
                    }
                }
                last = cur;
            }
            */
        }

        std::vector<size_t> annotation_from_kmer(KMer &kmer) {
            return annotation.find(kmer.begin(), kmer.end());
        }

        std::vector<size_t> annotation_exact_from_kmer(KMer &kmer) {
            return annotation_exact.find(kmer.begin(), kmer.end());
        }

        KMer kmer_from_index(const DBG_succ::node_iterator &i) {
            auto kmer = graph->get_node_seq(i);
            kmer.emplace_back(graph->get_W(i) % DBG_succ::alph_size);
            return KMer(KMer::pack_kmer(kmer.begin(), kmer.size()));
        }

        std::vector<size_t> get_annotation(const DBG_succ::node_iterator &i) {
            KMer int_kmer = kmer_from_index(i);
            return annotation_from_kmer(int_kmer);
        }


        std::vector<size_t> get_annotation_corrected(const DBG_succ::node_iterator &i) {
            std::unordered_set<size_t> visited;

            KMer int_kmer = kmer_from_index(i);
            auto curannot = annotation_from_kmer(int_kmer);

            KMer ckmer = int_kmer;
            DBG_succ::node_iterator j = i;
            TAlphabet curW = graph->get_W(i) % DBG_succ::alph_size;
			TAlphabet newW;
            visited.insert(i);
            size_t pcount_old = annotate::HashAnnotation<>::bigint_popcount(curannot);
            size_t pcount_new = 0;
            auto nextannot = curannot;
            size_t path = 0;
            std::string str_kmer = ckmer.to_string(DBG_succ::alphabet).substr(0,2);
            if (str_kmer[0] != '$' && str_kmer[1] != '$' && pcount_old > 1) {
                while (true) {
                    j = graph->outgoing(j, graph->get_W(j));
                    newW = graph->get_W(j) % DBG_succ::alph_size;
                    if (!newW) {
                        break;
                    }   
                    update_kmer(graph->get_k(), newW, curW, ckmer.begin_256());
                    if (graph->indegree(j) > 1 
                     || annotation.is_junction(ckmer.begin(), ckmer.end()))
                        break;
                    if (!visited.insert(j).second) {
                        break;
                    }   
    
                    curW = newW;
                    nextannot = annotation.find(ckmer.begin(), ckmer.end());
                    annotate::HashAnnotation<>::merge_and(curannot, nextannot);
                    pcount_new = annotate::HashAnnotation<>::bigint_popcount(curannot);
                    if (pcount_new < pcount_old) {
                        path = 0;
                    } else {
                        path++;
                    }   
                    if (pcount_new == 1 || path > 5) {
                        break;
                    }   
                }   
            }   
            //TODO:move backward


			return curannot;
        }

        std::vector<uint8_t> test_fp(const DBG_succ::node_iterator &i) {
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
            size_t total = graph->get_W().size() / step;
            for (DBG_succ::node_iterator i = 1; i < graph->get_W().size(); i += step) {
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
        DBG_succ *graph = NULL;
        double bloom_size_factor = 0;
        std::vector<size_t> sizes_v;

        //TODO: set this value
        annotate::HashAnnotation<annotate::BloomFilter<5>> annotation;
        annotate::HashAnnotation<annotate::ExactFilter> annotation_exact;
};
