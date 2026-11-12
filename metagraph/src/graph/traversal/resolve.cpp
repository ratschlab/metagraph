#include "resolve.hpp"

#include <algorithm>
#include <cctype>

#include <tsl/hopscotch_map.h>

#include "graph/annotated_dbg.hpp"
#include "annotation/coord_to_header.hpp"


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
                               const ResolveOptions &options) {
    SupportProfile profile;
    profile.k = oracle.get_k();
    profile.regime = oracle.regime();

    if (options.labels.empty() == !options.discover)
        throw std::invalid_argument("Specify either an explicit label list or discovery, not both");
    if (query.size() < profile.k)
        throw std::invalid_argument("Query shorter than k");
    if (options.discover && !options.discover_max_labels)
        throw std::invalid_argument("Discovery needs discover_max_labels of at least 1");

    profile.num_kmers = query.size() - profile.k + 1;
    const std::vector<node_index> keys = map_to_nodes(oracle.graph(), query);
    assert(keys.size() == profile.num_kmers);

    profile.graph_runs = runs_of(profile.num_kmers, [&](uint64_t i) { return keys[i] != npos; });

    std::vector<LabelRef> refs;
    std::vector<LabelSupport> support;
    if (!options.labels.empty()) {
        for (const auto &name : options.labels) {
            refs.push_back(oracle.resolve_label(name));
        }
        LabelQuery query_labels(oracle, refs);
        support.resize(refs.size());
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
    } else {
        if (options.discover_kind == LabelKind::HEADER && !oracle.coord_to_header())
            throw std::invalid_argument("Header discovery requires a CoordToHeader index");
        // the ranking of /search with signatures, so that both report the same labels,
        // counts and runs
        const bool headers = options.discover_kind == LabelKind::HEADER;
        size_t total = 0;
        auto top = oracle.anno_graph().get_top_label_signatures(
                query, options.discover_max_labels, 0.0, 0.0, &total, headers);
        for (const auto &[name, count, presence] : top) {
            // a name /traverse would refuse as ambiguous is left out, so that every label
            // can be resubmitted
            if (oracle.is_ambiguous(name))
                continue;
            LabelRef ref = oracle.resolve_label(name);
            assert(ref.kind == options.discover_kind);
            LabelSupport found;
            for (uint64_t i = 0; i < presence.size(); ++i) {
                if (presence[i])
                    found.add(i);
            }
            refs.push_back(std::move(ref));
            support.push_back(std::move(found));
        }
        // kept: the labels returned, after the ambiguous ones are left out
        if (total > options.discover_max_labels)
            profile.labels_truncated = LabelTruncation { refs.size(), total };
    }

    for (LabelId l = 0; l < refs.size(); ++l) {
        LabelProfile lp;
        lp.label = refs[l];
        lp.kmers_supported = support[l].kmers;
        lp.runs = std::move(support[l].runs);
        profile.labels.push_back(std::move(lp));
    }
    // discovered labels, most supported first: the ranking only sorts when it cuts the list
    if (options.discover) {
        std::sort(profile.labels.begin(), profile.labels.end(),
                  [](const LabelProfile &a, const LabelProfile &b) {
                      return std::make_tuple(b.kmers_supported, a.label.column, a.label.seq_id)
                           < std::make_tuple(a.kmers_supported, b.label.column, b.label.seq_id);
                  });
    }
    return profile;
}

} // namespace traversal
} // namespace graph
} // namespace mtg
