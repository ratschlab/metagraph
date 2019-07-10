#include <iostream>
#include <regex>
#include <cstdlib>

using namespace std;


int main(){
    enum op {M, X, D, I, S};
    string cigar;
    uint64_t qlen;
    uint64_t qlen_from_cigar = 0;
    regex match_regex("[0-9]+M");
    regex del_regex("[0-9]+D");
    regex in_regex("[0-9]+I");
    regex mis_match_regex("[0-9]+X");
    regex clipped_regex("[0-9]+S");
    pair<regex, int8_t> regex_array[] = {{match_regex, op::M}, {del_regex, op::D},
                                         {in_regex, op::I}, {mis_match_regex, op::X},
                                         {clipped_regex, op::S}};
    smatch match;
    while (cin >> cigar >> qlen) {
        int64_t score = 0;
        qlen_from_cigar = 0;
        for (auto it = std::begin(regex_array); it != std::end(regex_array); ++ it) {
            auto regex_pair = *it;
            string cigar_cpy = cigar;
            while (regex_search(cigar_cpy, match, regex_pair.first)) {
                string score_str = match.str();
                score_str.pop_back();
                auto single_score = atoi(score_str.c_str());
                if (regex_pair.second == op::M) {
                    score += single_score * 2;
                    qlen_from_cigar += single_score;
                } else if (regex_pair.second == op::I) {
                    score -= (single_score - 1) + 3 ;
                    qlen_from_cigar += single_score;
                } else if (regex_pair.second == op::D) {
                    score -= (single_score - 1) + 3 ;
                } else if (regex_pair.second == op::X) {
                    score -= single_score * 2 ;
                    qlen_from_cigar += single_score;
                } else if (regex_pair.second == op::S) {
                    qlen_from_cigar += single_score;
                }
                cigar_cpy = match.suffix();
            }
        }
        if (qlen == qlen_from_cigar)
            cout << score << " " << qlen << " " << qlen_from_cigar << endl;
    }
    return 0;
}
