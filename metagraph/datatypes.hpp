#ifndef __DATATYPES_HPP__
#define __DATATYPES_HPP__

#include <functional>
#include <set>
#include <unordered_set>

class DBG_succ;

struct HitInfo {
    uint64_t rl;
    uint64_t ru;
    uint64_t str_pos;
    uint64_t graph_pos;
    uint64_t distance;
    std::string cigar;
    std::vector<uint64_t> path;

    HitInfo(uint64_t rl_, uint64_t ru_, uint64_t str_pos_, uint64_t graph_pos_, uint64_t distance_, std::string cigar_, std::vector<uint64_t> path_):
        rl(rl_),
        ru(ru_),
        str_pos(str_pos_),
        graph_pos(graph_pos_),
        distance(distance_),
        cigar(cigar_),
        path(path_) {}

    HitInfo(const HitInfo& other):
        rl(other.rl),
        ru(other.ru),
        str_pos(other.str_pos),
        graph_pos(other.graph_pos),
        distance(other.distance),
        cigar(other.cigar),
        path(other.path) {}
};

class HitInfoCompare {
    bool is_reverse;
public:
    HitInfoCompare(const bool& is_reverse_ = false) {
        is_reverse = is_reverse_;
    }

    bool operator() (const HitInfo& lhs, const HitInfo& rhs) const {
        if (is_reverse) return (lhs.distance < rhs.distance);
        else return (lhs.distance > rhs.distance);
    }
};

struct AnnotationSet {
    std::set<std::string> annotation;
};

struct AnnotationHash {
    std::uint16_t operator()(const std::set<std::string> &a) const {
        std::set<std::string>::iterator it = a.begin();
        std::uint16_t h1 = (uint16_t) std::hash<std::string>{}(*it);
        it++;
        if (a.size() > 1) {
            for (; it != a.end(); ++it) {
                std::uint16_t h2 = (uint16_t) std::hash<std::string>{}(*it);
                h1 = h1 ^ (h2 << 1);
            }
        }
        return h1;
    }
};

struct ParallelMergeContainer {
    std::vector<std::pair<uint64_t, uint64_t> > bins_g1;
    std::vector<std::pair<uint64_t, uint64_t> > bins_g2;
    std::vector<DBG_succ*> result;
    DBG_succ* graph1;
    DBG_succ* graph2;
    unsigned int idx;
    unsigned int k;
    unsigned int bins_done;
};

#endif
