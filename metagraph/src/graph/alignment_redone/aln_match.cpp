#include "aln_match.hpp"

namespace mtg::graph::align_redone {

using mtg::common::logger;

std::string spell_path(const DeBruijnGraph &graph,
                       const std::vector<DeBruijnGraph::node_index> &path) {
    std::string seq;

    if (path.empty())
        return seq;

    seq.reserve(path.size() + graph.get_k() - 1);

    size_t num_dummy = 0;
    size_t num_unknown = 0;
    if (path.front()) {
        seq += graph.get_node_sequence(path.front());
    } else {
        num_unknown = graph.get_k();
        seq += std::string(num_unknown, '$');
        num_dummy = 1;
    }

    for (size_t i = 1; i < path.size(); ++i) {
        if (path[i]) {
            if (num_dummy) {
                seq += '$';
                ++num_unknown;
                std::string next_seq = graph.get_node_sequence(path[i]);
                auto it = seq.end() - next_seq.size();
                for (char c : next_seq) {
                    if (*it == '$' && c != '$') {
                        --num_unknown;
                        *it = c;
                    }

                    ++it;
                }
                num_dummy = 0;
            } else {
                char next = '\0';
                graph.call_outgoing_kmers(path[i - 1], [&](auto next_node, char c) {
                    if (next_node == path[i])
                        next = c;
                });

                if (!next) {
                    logger->error("Invalid edge {} {}\t{} {}",
                                  path[i - 1], path[i],
                                  graph.get_node_sequence(path[i - 1]),
                                  graph.get_node_sequence(path[i]));
                    throw std::runtime_error("");
                }

                seq += next;
                if (num_unknown) {
                    std::string next_seq = graph.get_node_sequence(path[i]);
                    auto it = seq.end() - next_seq.size();
                    for (char c : next_seq) {
                        if (*it == '$' && c != '$') {
                            --num_unknown;
                            *it = c;
                        }

                        ++it;
                    }
                }
            }
        } else {
            seq += '$';
            ++num_dummy;
            ++num_unknown;
        }
    }

    assert(seq.size() == path.size() + graph.get_k() - 1);
    return seq;
}

bool operator==(const Match &a, const Match &b) {
    return a.get_orientation() == b.get_orientation()
        && a.get_query().begin() == b.get_query().begin()
        && a.get_query().end() == b.get_query().end()
        && a.get_clipping() == b.get_clipping()
        && a.get_end_clipping() == b.get_end_clipping()
        && a.get_end_trim() == b.get_end_trim()
        && a.get_path() == b.get_path();
}

bool operator==(const Anchor &a, const Anchor &b) {
    return static_cast<const Match&>(a) == static_cast<const Match&>(b)
        && a.get_label_class() == b.get_label_class();
}

bool operator==(const Alignment &a, const Alignment &b) {
    return static_cast<const Match&>(a) == static_cast<const Match&>(b)
        && a.get_cigar() == b.get_cigar();
}

std::ostream& operator<<(std::ostream &out, const Match &a) {
    Cigar cigar = a.generate_cigar();
    out << fmt::format("{}\t{}\t*\t{}\t{}\t{}",
                       !a.get_orientation() ? "+" : "-",
                       a.get_spelling(),
                       cigar.get_num_matches(),
                       cigar.to_string(),
                       a.get_end_trim());
    return out;
}

size_t Anchor::trim_end() {
    assert(get_spelling().size() == get_path_spelling().size() - get_end_trim());
    if (suffix_.empty())
        return 0;

    size_t num_trimmed = 0;

#ifndef NDEBUG
    std::string old_spelling(get_spelling());
#endif

    while (suffix_.size() && path_.size() > 1) {
        suffix_.pop_back();
        path_.pop_back();
        ++num_trimmed;
        assert(get_spelling().size() == get_path_spelling().size() - get_end_trim());
        assert(get_spelling() == old_spelling);
    }

    assert(get_spelling() == old_spelling);

    return num_trimmed;
}

size_t Alignment::trim_end() {
    assert(get_spelling().size() == get_path_spelling().size() - get_end_trim());
    if (end_trim_ == 0)
        return 0;

#ifndef NDEBUG
    std::string old_spelling(get_spelling());
#endif

    size_t num_trimmed = 0;
    while (end_trim_ && path_.size() > 1) {
        assert(path_spelling_.size());
        path_spelling_.pop_back();
        path_.pop_back();
        --end_trim_;
        ++num_trimmed;
        assert(get_spelling().size() == get_path_spelling().size() - get_end_trim());
        assert(get_spelling() == old_spelling);
    }

    assert(get_spelling() == old_spelling);

    return num_trimmed;
}

void Anchor::append(const Anchor &other, const DeBruijnGraph *graph) {
    assert(query_.data() == other.query_.data());
    assert(query_.data() + query_.size() == other.query_.data() + other.query_.size());
    assert(label_class_ == other.label_class_);

    if (empty()) {
        *this = other;
        return;
    }

    assert(other.seed_.data() >= seed_.data());

    if (seed_.data() + seed_.size() >= other.seed_.data() + other.seed_.size())
        return;

    if (seed_.data() == other.seed_.data()) {
        assert(coord_ == other.coord_);
        *this = other;
        return;
    }

    auto next = path_.begin() + (other.seed_.data() - seed_.data());
    assert(next <= path_.end());

    seed_ = std::string_view(
        seed_.data(),
        other.seed_.data() + other.seed_.size() - seed_.data()
    );

    path_.erase(next, path_.end());
    path_.insert(path_.end(), other.path_.begin(), other.path_.end());
    suffix_ = other.suffix_;

    assert(!graph || is_spelling_valid(*graph));
}

} // namespace mtg::graph::align

namespace std {

std::size_t hash<mtg::graph::align_redone::Match>::operator()(const mtg::graph::align_redone::Match &a) const {
    return (a.get_path().size() ? a.get_path()[0] : 0) + a.get_clipping() + a.get_end_clipping();
}

std::size_t hash<mtg::graph::align_redone::Anchor>::operator()(const mtg::graph::align_redone::Anchor &a) const {
    return std::hash<mtg::graph::align_redone::Match>()(a);
}

std::size_t hash<mtg::graph::align_redone::Alignment>::operator()(const mtg::graph::align_redone::Alignment &a) const {
    return std::hash<std::string>()(a.get_cigar().to_string());
}

} // namespace std
