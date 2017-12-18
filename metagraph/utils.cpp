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
    //std::cerr << "k1_val: " << k1_val.first << " k2_val: " << k2_val.first << " curr_k:" << curr_k << std::endl;
    return std::make_pair(k1_val.first < k2_val.first, k2_val.first < k1_val.first);
}


std::pair<std::vector<bool>, uint64_t> compare_nodes(std::vector<const DBG_succ*> G,
                                                     std::vector<uint64_t> k,
                                                     std::vector<uint64_t> n,
                                                     size_t &cnt) {

    std::vector<bool> result (G.size(), false);
    std::vector<bool> ignore (G.size(), false);
    std::vector<std::pair<TAlphabet, uint64_t> > k_val;
    std::pair<TAlphabet, uint64_t> min;
    std::vector<uint64_t> k_tmp (k);
    size_t s = G.size();

    uint64_t curr_k = 0;
    while (curr_k < G.at(0)->get_k()) {
        k_val.clear();
        for (size_t i = 0; i < s; i++) {
            //std::cerr << "curr_k: " << curr_k << " - i: " << i; 
            if ((k.at(i) < n.at(i)) && !ignore.at(i)) {
                k_val.push_back(G.at(i)->get_minus_k_value(k_tmp.at(i), 0));
            } else {
                k_val.push_back(std::make_pair(G.at(0)->alph_size, 0));
                ignore.at(i) = true;
            }
            //std::cerr << " k_val.first: " << k_val.back().first << " k_val.second: " << k_val.back().second << std::endl;
        }
        
        min = *std::min_element(k_val.begin(), k_val.end(),
                [](const std::pair<TAlphabet, uint64_t> &a,
                   const std::pair<TAlphabet, uint64_t> &b) {
                     return a.first < b.first;
                }
        );
        cnt = 0;
        for (size_t i = 0; i < s; i++)
            if ((k_val.at(i).first == min.first) && !ignore.at(i)) {
                cnt++;
            } else {
                ignore.at(i) = true;
            }
        if (cnt == 1)
            break;
        ++curr_k;
        for (size_t i = 0; i < s; i++) {
            if (!ignore.at(i))
                k_tmp.at(i) = k_val.at(i).second;
        }
    }

    //std::cerr << "cnt: " << cnt << " s: " << s << std::endl;
    if (cnt == 0) {
        return std::make_pair(result, 0L);
    } else {
        uint64_t min_edge = G.at(0)->alph_size;    
        uint64_t max_edge_val = 0; 
        // get minimal outgoing edge
        for (size_t i = 0; i < s; i++) {
            if ((k_val.at(i).first == min.first) && !ignore.at(i))
                min_edge = std::min(min_edge, G.at(i)->get_W(k.at(i)) % G.at(0)->alph_size);
        }

        for (size_t i = 0; i < s; i++) {
            //std::cerr << "i: " << i << " k: " << k.at(i) << " W: " << G.at(i)->get_W(k.at(i)) << " min.f: " << min.first <<  " result: ";
            result.at(i) = ((k_val.at(i).first == min.first) && (G.at(i)->get_W(k.at(i)) % G.at(i)->alph_size <= min_edge) && !ignore.at(i));
            //std::cerr << result.at(i) << std::endl;
            if (result.at(i))
                max_edge_val = std::max(max_edge_val, G.at(i)->get_W(k.at(i)));
        }
        return std::make_pair(result, max_edge_val);
    }
}


/**
 *  Returns the input file type, given a filename
 */
std::string get_filetype(const std::string &fname) {
    int dotind = fname.rfind(".");
    std::string ext;
    if (dotind < 0) {
        std::cerr << "ERROR: Filetype unknown\n";
        exit(1);
    }
    if (fname.substr(dotind) == ".gz") {
        int nextind = fname.substr(0, dotind - 1).rfind(".");
        ext = fname.substr(nextind, dotind - nextind);
    } else {
        ext = fname.substr(dotind);
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
