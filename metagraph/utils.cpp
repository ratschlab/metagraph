#include "utils.hpp"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>


namespace utils {

typedef uint64_t TAlphabet;

uint64_t kFromFile(const std::string &infbase) {
    uint64_t k = 0;
    std::ifstream instream((infbase + ".F.dbg").c_str()); 
    std::string line;
    size_t mode = 0;
    while (std::getline(instream, line)) {
        if (strcmp(line.c_str(), ">F") == 0 || strcmp(line.c_str(), ">p") == 0) {
            mode = 1;
        } else if (strcmp(line.c_str(), ">k") == 0) {
            mode = 2;
        } else {
            if (mode == 2) {
                k = strtoul(line.c_str(), NULL, 10);
            } else if (mode != 1) {
                fprintf(stderr, "ERROR: input file corrupted\n");
                exit(1);
            }
        }
    }
    instream.close();
    return k;
}

/**
* This function takes a pointer to a graph structure G1 and a corresponding node index k1_node
* as well as a pointer to a second graph structure G2 and a corresponding node index k2_node. It
* returns a pair of bool with the first value set to true if G1(k1_node) < G2(k2_node) and the
* second value set to true if G2(k2_node) < G1(k1_node).
*/
std::pair<bool, bool> compare_nodes(const DBG_succ *G1, uint64_t k1_node,
                                    const DBG_succ *G2, uint64_t k2_node) {
    assert(G1->get_k() == G2->get_k());

    std::pair<TAlphabet, uint64_t> k1_val;
    std::pair<TAlphabet, uint64_t> k2_val;

    for (uint64_t curr_k = 0; curr_k < G1->get_k(); ++curr_k) {
        k1_val = G1->get_minus_k_value(k1_node, 0);
        k2_val = G2->get_minus_k_value(k2_node, 0);
        if (k1_val.first != k2_val.first) {
            break;
        }
        k1_node = k1_val.second;
        k2_node = k2_val.second;
    }

    return std::make_pair(k1_val.first < k2_val.first, k2_val.first < k1_val.first);
}


std::vector<bool> smallest_nodes(const std::vector<const DBG_succ*> &G,
                                 const std::vector<uint64_t> &k,
                                 const std::vector<uint64_t> &n) {
    std::vector<bool> ignore(G.size());
    for (size_t i = 0; i < G.size(); i++) {
        ignore[i] = (k[i] >= n[i]);
    }
    std::vector<TAlphabet> k_val(G.size());
    std::vector<uint64_t> k_tmp = k;

    for (uint64_t curr_k = 0; curr_k < G.at(0)->get_k(); ++curr_k) {
        TAlphabet min = DBG_succ::alph_size;
        for (size_t i = 0; i < G.size(); i++) {
            if (!ignore[i]) {
                auto ret = G[i]->get_minus_k_value(k_tmp[i], 0);
                k_val[i] = ret.first;
                k_tmp[i] = ret.second;
                min = std::min(min, ret.first);
            }
        }

        size_t cnt = 0;
        for (size_t i = 0; i < G.size(); i++) {
            if (ignore[i])
                continue;

            if (k_val[i] == min) {
                cnt += 1;
            } else {
                ignore[i] = true;
            }
        }

        if (cnt == 0)
            return std::vector<bool>(G.size(), false);

        if (cnt == 1)
            break;
    }

    uint64_t min_edge = DBG_succ::alph_size;    
    // get minimal outgoing edge
    std::vector<TAlphabet> edge_label(G.size());
    for (size_t i = 0; i < G.size(); i++) {
        if (!ignore[i]) {
            edge_label[i] = G[i]->get_W(k[i]) % DBG_succ::alph_size;
            min_edge = std::min(min_edge, edge_label[i]);
        }
    }

    std::vector<bool> smallest(G.size(), false);
    for (size_t i = 0; i < G.size(); i++) {
        if (!ignore[i] && edge_label[i] == min_edge) {
            smallest[i] = true;
        }
    }
    return smallest;
}


/**
 *  Returns the input file type, given a filename
 */
std::string get_filetype(const std::string &fname) {
    size_t dotind = fname.rfind(".");
    if (dotind == std::string::npos)
        return "";

    std::string ext = fname.substr(dotind);

    if (ext == ".gz") {
        size_t nextind = fname.substr(0, dotind - 1).rfind(".");
        if (nextind == std::string::npos)
            return "";

        ext = fname.substr(nextind, dotind - nextind);
    }

    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    if (ext == ".vcf") {
        return "VCF";
    } else if (ext == ".fq") {
        return "FASTQ";
    } else {
        return "FASTA";
    }
}

/**
 * Given a minimum number of splits,
 * generate a list of suffices from the alphabet.
 */
std::deque<std::string> generate_strings(const std::string &alphabet,
                                         size_t length) {

    std::deque<std::string> suffices = { "" };
    while (suffices[0].length() < length) {
        for (const char c : alphabet) {
            suffices.push_back(c + suffices[0]);
        }
        suffices.pop_front();
    }
    assert(suffices.size() == std::pow(alphabet.size(), length));
    return suffices;
}


} // namespace utils
