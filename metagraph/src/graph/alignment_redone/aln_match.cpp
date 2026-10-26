#include "aln_match.hpp"

#include "graph/alignment/alignment.hpp"

namespace mtg::graph::align_redone {

std::string spell_path(const DeBruijnGraph &graph,
                       const std::vector<DeBruijnGraph::node_index> &path) {
    if (path.empty())
        return {};

    return align::spell_path(graph, path, 0);
}

bool operator==(const Match &a, const Match &b) {
    return a.get_orientation() == b.get_orientation()
        && a.get_query().begin() == b.get_query().begin()
        && a.get_query().end() == b.get_query().end()
        && a.get_clipping() == b.get_clipping()
        && a.get_end_clipping() == b.get_end_clipping()
        && a.get_offset() == b.get_offset()
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

size_t Alignment::trim_end() {
    assert(get_offset() + get_spelling().size() + get_end_trim() == path_spelling_.size());
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
        label_classes_.pop_back();
        --end_trim_;
        ++num_trimmed;
        assert(get_offset() + get_spelling().size() + get_end_trim() == path_spelling_.size());
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
