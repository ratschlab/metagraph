#include "dbg_bloom_annotator.hpp"

#include <fstream>


namespace annotate {

void PreciseAnnotator::add_sequence(const std::string &sequence, size_t column) {
    std::string preprocessed_seq = graph_.encode_sequence(sequence);

    // Don't annotate short sequences
    if (preprocessed_seq.size() < graph_.get_k() + 1)
        return;

    if (column >= annotation_exact.size())
        annotation_exact.resize(column + 1);

    for (size_t i = 0; i + graph_.get_k() < preprocessed_seq.size(); ++i) {
        annotation_exact.insert(
            &preprocessed_seq[i],
            &preprocessed_seq[i] + graph_.get_k() + 1, column
        );
    }
}

void PreciseAnnotator::add_column(const std::string &sequence) {
    add_sequence(sequence, annotation_exact.size());
}

std::vector<size_t>
PreciseAnnotator::annotation_from_kmer(const std::string &kmer) const {
    assert(kmer.length() == graph_.get_k() + 1);
    return annotation_exact.find(kmer.data(), kmer.data() + kmer.size());
}

BloomAnnotator::BloomAnnotator(size_t num_hash_functions,
                               const DeBruijnGraphWrapper &graph,
                               double bloom_size_factor,
                               bool verbose)
    : annotation(num_hash_functions),
      graph_(graph),
      bloom_size_factor_(bloom_size_factor),
      total_traversed_(0),
      verbose_(verbose) {}

void BloomAnnotator::add_sequence(const std::string &sequence, size_t column) {
    std::string preprocessed_seq = graph_.encode_sequence(sequence);

    // Don't annotate short sequences
    if (preprocessed_seq.size() < graph_.get_k() + 1)
        return;

    if (column >= annotation.size())
        annotation.resize(column + 1);

    if (annotation[column].size() == 0) {
        annotation[column].resize(
            bloom_size_factor_
                * (preprocessed_seq.size() - graph_.get_k()) + 1
        );
    }

    for (size_t i = 0; i + graph_.get_k() < preprocessed_seq.size(); ++i) {
        annotation.insert(
            &preprocessed_seq[i],
            &preprocessed_seq[i] + graph_.get_k() + 1,
            column
        );
    }
}

void BloomAnnotator::add_column(const std::string &sequence) {
    add_sequence(sequence, annotation.size());
}

std::vector<size_t>
BloomAnnotator::annotation_from_kmer(const std::string &kmer) const {
    assert(kmer.length() == graph_.get_k() + 1);
    return annotation.find(kmer.data(), kmer.data() + kmer.size());
}

std::vector<size_t>
BloomAnnotator::get_annotation(DeBruijnGraphWrapper::edge_index i) const {
    return annotation_from_kmer(kmer_from_index(i));
}

std::vector<size_t>
BloomAnnotator::get_annotation_corrected(DeBruijnGraphWrapper::edge_index i,
                                         size_t path_cutoff) const {
    //initial raw annotation
    std::string orig_kmer = kmer_from_index(i);

    auto curannot = annotation_from_kmer(orig_kmer);

    // Dummy edges are not supposed to be annotated
    if (graph_.is_dummy_edge(orig_kmer)) {
        curannot.assign(curannot.size(), 0);
        return curannot;
    }

    size_t pcount_old
        = annotate::popcount(curannot);

    if (!pcount_old)
        return curannot;

    std::string cur_kmer = orig_kmer;
    auto j = i;
    size_t path = 0;
    while (path++ < path_cutoff) {
        const_cast<size_t&>(total_traversed_)++;

        //traverse forward
        j = graph_.next_edge(j, cur_kmer.back());

        cur_kmer = kmer_from_index(j);

        //check outdegree
        if (graph_.is_dummy_edge(cur_kmer)
                || !graph_.has_the_only_outgoing_edge(j))
            break;

        //bitwise AND annotations
        auto nextannot = annotate::merge_and(
            curannot,
            annotation_from_kmer(cur_kmer)
        );

        //check popcounts
        size_t pcount_new
            = annotate::popcount(nextannot);

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
    while (graph_.has_the_only_incoming_edge(j)
            && path++ < path_cutoff) {
        const_cast<size_t&>(total_traversed_)++;

        j = graph_.prev_edge(j);

        cur_kmer = kmer_from_index(j);
        if (graph_.is_dummy_edge(cur_kmer))
            break;

        auto nextannot = annotate::merge_and(
            curannot,
            annotation_from_kmer(cur_kmer)
        );

        auto pcount_new
            = annotate::popcount(nextannot);

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

void BloomAnnotator::test_fp_all(const PreciseAnnotator &annotation_exact,
                                 size_t step) const {
    size_t fp = 0;
    size_t fp_pre = 0;
    size_t fn = 0;
    size_t total = 0;
    for (auto i = graph_.first_edge(); i <= graph_.last_edge(); i += step) {
        if (graph_.is_dummy_edge(kmer_from_index(i)))
            continue;
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

void BloomAnnotator::serialize(std::ostream &out) const {
    annotation.serialize(out);
}

void BloomAnnotator::serialize(const std::string &filename) const {
    std::ofstream out(filename + ".annot.dbg");
    serialize(out);
    out.close();
}

std::vector<size_t>
BloomAnnotator::unpack(const std::vector<size_t> &packed) {
    std::vector<size_t> labels;
    for (size_t i = 0; i < packed.size() * 64; ++i) {
        if (packed[i / 64] & (1llu << (i % 64)))
            labels.push_back(i);
    }
    return labels;
}

std::string
BloomAnnotator::kmer_from_index(DeBruijnGraphWrapper::edge_index index) const {
    return graph_.get_node_kmer(index) + graph_.get_edge_label(index);
}

std::vector<uint8_t>
BloomAnnotator::test_fp(DeBruijnGraphWrapper::edge_index i,
                        const PreciseAnnotator &annotation_exact) const {

    auto int_kmer = kmer_from_index(i);

    auto test = annotation_from_kmer(int_kmer);
    auto test_exact = annotation_exact.annotation_from_kmer(int_kmer);

    auto curannot = get_annotation_corrected(i);

    auto jt = test.begin();
    auto kt = test_exact.begin();
    auto lt = curannot.begin();

    std::vector<uint8_t> stats(3, 0);

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
            if (verbose_)
                std::cout << "FP: " << int_kmer << std::endl;
        }
        //false negatives after correction
        if (!stats[2] && (*lt | *kt) != *lt) {
            stats[2] = 1;
            if (verbose_) {
                std::cout << "FN: " << int_kmer << std::endl;
                auto unpacked_labels = unpack(test_exact);
                std::cout << "True annotation:\t";
                for (size_t value : unpacked_labels) {
                    std::cout << value << " ";
                }
                std::cout << std::endl;
                unpacked_labels = unpack(curannot);
                std::cout << "Corrected annotation:\t";
                for (size_t value : unpacked_labels) {
                    std::cout << value << " ";
                }
                std::cout << std::endl;
            }
        }
        if (stats[0] && stats[1] && stats[2]) {
            break;
        }
    }
    return stats;
}

} // namespace annotate
