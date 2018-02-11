#include "dbg_hash.hpp"

const std::string kAlphabet = "ACGTN$";


std::string DBGHash::encode_sequence(const std::string &sequence) const {
    std::string result = sequence;

    for (char &c : result) {
        if (kAlphabet.find(c) == std::string::npos)
            c = 'N';
    }
    return result;
}

std::string DBGHash::get_node_kmer(edge_index i) const {
    std::string kmer = kmers_[i];
    kmer.pop_back();
    return kmer;
}

bool DBGHash::has_the_only_outgoing_edge(edge_index i) const {
    std::string kmer = kmers_[i];

    size_t num_edges = 0;
    for (char c : kAlphabet) {
        kmer.back() = c;
        num_edges += (indices_.find(kmer) != indices_.end());
    }

    return num_edges == 1;
}

bool DBGHash::has_the_only_incoming_edge(edge_index i) const {
    std::string kmer("$");
    kmer += kmers_[i];
    kmer.pop_back();

    size_t num_edges = 0;
    for (char c : kAlphabet) {
        kmer.front() = c;
        num_edges += (indices_.find(kmer) != indices_.end());
    }

    return num_edges == 1;
}

bool DBGHash::is_dummy_edge(const std::string &kmer) const {
    assert(kmer.length() == k_ + 1);
    for (char c : kmer) {
        if (is_dummy_label(c))
            return true;
    }
    return false;
}

bool DBGHash::is_dummy_label(char c) const {
    return c == '$';
}

DBGHash::edge_index DBGHash::next_edge(edge_index i, char edge_label) const {
    std::string kmer = kmers_[i];
    kmer.back() = edge_label;

    kmer.erase(kmer.begin());
    kmer.push_back('$');

    for (char c : kAlphabet) {
        kmer.back() = c;
        auto it = indices_.find(kmer);
        if (it != indices_.end())
            return it->second;
    }
    assert(false && "Can't traverse if there are no outgoing edges");
    return 0;
}

DBGHash::edge_index DBGHash::prev_edge(edge_index i) const {
    std::string kmer("$");
    kmer += kmers_[i];
    kmer.pop_back();

    for (char c : kAlphabet) {
        kmer.front() = c;
        auto it = indices_.find(kmer);
        if (it != indices_.end())
            return it->second;
    }
    assert(false && "Can't traverse if there are not incoming edges");
    return 0;
}

void DBGHash::add_sequence(const std::string &sequence, bool rooted) {
    // Don't annotate short sequences
    if (sequence.size() < k_ + 1)
        return;

    std::string transformed_seq = (!rooted ? std::string(k_ + 1, '$') : std::string())
                                    + encode_sequence(sequence)
                                    + (!rooted ? std::string(1, '$') : std::string());

    for (size_t i = 0; i + k_ < transformed_seq.size(); ++i) {
        std::string kmer = transformed_seq.substr(i, k_ + 1);
        if (indices_.find(kmer) == indices_.end()) {
            indices_[kmer] = kmers_.size();
            kmers_.push_back(kmer);
        }
    }
}

size_t DBGHash::get_num_edges() const {
    return kmers_.size();
}
