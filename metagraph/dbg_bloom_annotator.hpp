#ifndef __DBG_BLOOM_ANNOTATOR_HPP__
#define __DBG_BLOOM_ANNOTATOR_HPP__

#include <unordered_set>

#include "dbg_succinct.hpp"
#include "annotate.hpp"
#include "kmer.hpp"
#include "hashers.hpp"


//TODO
class DeBruijnGraphWrapper {
  public:
    typedef uint64_t edge_index;
};


// class DBGSuccAnnotWrapper : public DeBruijnGraphWrapper {
//   public:
//     DBGSuccAnnotWrapper(const DBG_succ *graph) : graph_(graph) {}

//     size_t get_k() const { return graph_->get_k(); }

//     edge_index first_edge() const { return 1; }
//     edge_index last_edge() const { return graph_->get_W()->size() - 1; }

//     std::string encode_sequence(const std::string &sequence) const {
//         auto result = sequence;
//         for (char &c : result) {
//             if (graph_.alphabet.find(c) == std::string::npos)
//                 c = 'N';
//         }
//         return result;
//     }

//     std::string get_node_kmer(edge_index i) const {
//         return 'C';
//     }

//     char get_edge_label(edge_index i) const {
//         return 'C';
//     }

//     bool is_the_only_outgoing(edge_index i) const {
//         return true;
//     }

//     bool is_node_indegree_one(edge_index i) const {
//         return true;
//     }

//   private:
//     const DBG_succ *graph_;
// };


class PreciseAnnotator {
  public:
    PreciseAnnotator(size_t k) : k_(k) {}

    void add_sequence(const std::string &sequence, size_t column) {
        assert(sequence.size() >= k_ + 1);

        if (column >= annotation_exact.size())
            annotation_exact.resize(column + 1);

        for (size_t i = 0; i + k_ + 1 < sequence.size(); ++i) {
            annotation_exact.insert(&sequence[i], &sequence[i] + k_ + 1, column);
        }
    }

    void add_column(const std::string &sequence) {
        add_sequence(sequence, annotation_exact.size());
    }

    std::vector<size_t> annotation_from_kmer(const std::string &kmer) const {
        return annotation_exact.find(kmer.data(), kmer.data() + kmer.size());
    }

  private:
    size_t k_;

    annotate::HashAnnotation<annotate::ExactFilter> annotation_exact;
};


class BloomAnnotator {
  public:
    BloomAnnotator(size_t num_hash_functions,
                   DBG_succ *graph, double bloom_size_factor)
        : annotation(num_hash_functions),
          graph_(graph),
          bloom_size_factor_(bloom_size_factor),
          total_traversed_(0) {}

    void add_sequence(const std::string &sequence, size_t column) {
        assert(graph_);
        assert(sequence.size() >= graph_->get_k() + 1);

        if (column >= annotation.size())
            annotation.resize(column + 1);

        if (annotation[column].size() == 0) {
            annotation[column].resize(
                bloom_size_factor_ * (sequence.size() - graph_->get_k()) + 1
            );
        }
        for (size_t i = 0; i + graph_->get_k() + 1 < sequence.size(); ++i) {
            annotation.insert(&sequence[i],
                              &sequence[i] + graph_->get_k() + 1,
                              column);
        }
    }

    void add_column(const std::string &sequence) {
        add_sequence(sequence, annotation.size());
    }

    std::vector<size_t> annotation_from_kmer(const std::string &kmer) const {
        return annotation.find(kmer.data(), kmer.data() + kmer.size());
    }

    std::string kmer_from_index(DBG_succ::node_iterator index) {
        std::string kmer = graph_->get_node_str(index);
        return kmer + graph_->decode(graph_->get_W(index));
    }

    std::vector<size_t> get_annotation(DBG_succ::node_iterator index) {
        return annotation_from_kmer(kmer_from_index(index));
    }

    std::vector<size_t> get_annotation_corrected(DBG_succ::node_iterator i,
                                                 size_t path_cutoff = 50) {
        //initial raw annotation
        std::string orig_kmer = kmer_from_index(i);

        auto curannot = annotation_from_kmer(orig_kmer);

        // Dummy edges are not supposed to be annotated
        if (orig_kmer.find('$') != std::string::npos) {
            curannot.assign(curannot.size(), 0);
            return curannot;
        }

        size_t pcount_old = annotate::HashAnnotation<>::bigint_popcount(curannot);
        if (!pcount_old)
            return curannot;

        //TODO: checking for cur_kmer.back() and cur_kmer.front() means that no
        // dummy edges can be annotation. Either find a better heuristic, or find
        // a way to label them.

// std::cout << "\n\n\n\n" << orig_kmer << "\t";
// for (size_t i = 0; i < curannot.size(); ++i) {
//     if (curannot[i])
//         std::cout << " " << i;
// }
// std::cout << std::endl;

        std::string cur_kmer = orig_kmer;
        auto j = i;
        size_t path = 0;
        while (path++ < path_cutoff) {
            total_traversed_++;

            //move forward
            j = graph_->traverse(j, cur_kmer.back());

            cur_kmer = cur_kmer.substr(1) + graph_->decode(graph_->get_W(j));

            //check outdegree
            if (cur_kmer.back() == '$'
                    || !graph_->get_last(j) || !graph_->get_last(j - 1))
                break;

            //bitwise AND annotations
            auto nextannot = annotate::HashAnnotation<>::merge_and(
                curannot,
                annotation_from_kmer(cur_kmer)
            );

// std::cout << cur_kmer << "\t";
// for (size_t i = 0; i < nextannot.size(); ++i) {
//     if (nextannot[i])
//         std::cout << " " << i;
// }
// std::cout << std::endl;

            //check popcounts
            size_t pcount_new = annotate::HashAnnotation<>::bigint_popcount(nextannot);
            assert(pcount_new <= pcount_old);

            //if zero, then start of new sequence
            if (pcount_new == 0)
                break;

            //path length stopping conditions
            if (pcount_new < pcount_old) {
                curannot = nextannot;
                path = 0;
                pcount_old = pcount_new;
            }
        }

        //backward correction
        cur_kmer = orig_kmer;
        j = i;
        path = 0;
        while (graph_->indegree(j) == 1 && path++ < path_cutoff) {
            total_traversed_++;

            j = graph_->get_minus_k_value(j, 0).second;

            cur_kmer.pop_back();
            cur_kmer.insert(
                cur_kmer.begin(),
                graph_->decode(graph_->get_minus_k_value(j, graph_->get_k() - 1).first)
            );
            if (cur_kmer.front() == '$')
                break;

            auto nextannot = annotate::HashAnnotation<>::merge_and(
                curannot,
                annotation_from_kmer(cur_kmer)
            );

// std::cout << cur_kmer << "\t";
// for (size_t i = 0; i < nextannot.size(); ++i) {
//     if (nextannot[i])
//         std::cout << " " << i;
// }
// std::cout << std::endl;

            auto pcount_new = annotate::HashAnnotation<>::bigint_popcount(nextannot);
            assert(pcount_new <= pcount_old);

            if (pcount_new == 0)
                break;

            if (pcount_new < pcount_old) {
                curannot = nextannot;
                path = 0;
                pcount_old = pcount_new;
            }
        }

		return curannot;
    }

    std::vector<uint8_t> test_fp(DBG_succ::node_iterator i,
                                 const PreciseAnnotator &annotation_exact) {

        auto int_kmer = kmer_from_index(i);

        auto test = annotation_from_kmer(int_kmer);
        auto test_exact = annotation_exact.annotation_from_kmer(int_kmer);

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
                std::cout << "FN: " << int_kmer << std::endl;
            }
            if (stats[0] && stats[1] && stats[2]) {
                break;
            }
        }
        return stats;
    }

    void test_fp_all(const PreciseAnnotator &annotation_exact, size_t step = 1) {
        size_t fp = 0;
        size_t fp_pre = 0;
        size_t fn = 0;
        //size_t total = graph_->get_W().size() / step;
        size_t total = 0;
        for (DBG_succ::node_iterator i = 1; i < graph_->get_W().size(); i += step) {
            total++;
            auto stats = test_fp(i, annotation_exact);
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
    annotate::HashAnnotation<annotate::BloomFilter> annotation;

    DBG_succ *graph_;
    double bloom_size_factor_;

    //TODO: get rid of this if not using degree Bloom filter
    std::vector<size_t> sizes_v;

    size_t total_traversed_;
};

#endif // __DBG_BLOOM_ANNOTATOR_HPP__
