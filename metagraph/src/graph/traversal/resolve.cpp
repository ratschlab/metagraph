#include "resolve.hpp"

#include <algorithm>
#include <cctype>

#include <tsl/hopscotch_map.h>


namespace mtg {
namespace graph {
namespace traversal {

using annot::matrix::BinaryMatrix;
using annot::matrix::MultiIntMatrix;

std::string in_graph_case(std::string_view sequence) {
    std::string seq(sequence);
#if ! _DNA_CASE_SENSITIVE_GRAPH
    for (char &c : seq) {
        c = std::toupper(static_cast<unsigned char>(c));
    }
#endif
    return seq;
}

namespace {

// The rows a profile of explicit labels decodes at once
constexpr size_t kResolveBatchRows = 4096;

// The support of one label, built k-mer by k-mer in ascending order (each k-mer at most once)
struct LabelSupport {
    uint64_t kmers = 0;
    std::vector<KmerInterval> runs;

    void add(uint64_t i) {
        // present at the k-mer once, however many hits name it
        if (!runs.empty() && runs.back().end == i + 1)
            return;
        kmers++;
        if (!runs.empty() && runs.back().end == i) {
            runs.back().end = i + 1;
        } else {
            runs.push_back({ i, i + 1 });
        }
    }
};

} // namespace

SupportProfile resolve_support(const LabelOracle &oracle,
                               std::string_view query,
                               const std::vector<std::string> &labels) {
    SupportProfile profile;
    profile.k = oracle.get_k();
    profile.regime = oracle.regime();

    if (labels.empty())
        throw std::invalid_argument("Specify the labels to report on");
    if (query.size() < profile.k)
        throw std::invalid_argument("Query shorter than k");

    profile.num_kmers = query.size() - profile.k + 1;
    const std::vector<node_index> keys = map_to_nodes(oracle.graph(), query);
    assert(keys.size() == profile.num_kmers);

    profile.graph_runs = runs_of(profile.num_kmers, [&](uint64_t i) { return keys[i] != npos; });

    std::vector<LabelRef> refs;
    for (const auto &name : labels) {
        refs.push_back(oracle.resolve_label(name));
    }
    LabelQuery query_labels(oracle, refs);
    std::vector<LabelSupport> support(refs.size());
    // decode the rows kResolveBatchRows at a time
    for (size_t begin = 0; begin < keys.size(); begin += kResolveBatchRows) {
        const size_t end = std::min(keys.size(), begin + kResolveBatchRows);
        const auto hits = query_labels.fetch(
                std::vector<node_index>(keys.begin() + begin, keys.begin() + end));
        for (size_t i = 0; i < hits.size(); ++i) {
            for (LabelId l : hits[i]) {
                support[l].add(begin + i);
            }
        }
    }

    for (LabelId l = 0; l < refs.size(); ++l) {
        LabelProfile lp;
        lp.label = refs[l];
        lp.kmers_supported = support[l].kmers;
        lp.runs = std::move(support[l].runs);
        profile.labels.push_back(std::move(lp));
    }
    return profile;
}

} // namespace traversal
} // namespace graph
} // namespace mtg
