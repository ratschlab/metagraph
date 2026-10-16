#include "context_cursor.hpp"

#include "graph/graph_extensions/node_first_cache.hpp"
#include "common/vector_map.hpp"

#include <algorithm>
#include <cctype>
#include <map>
#include <stdexcept>
#include <unordered_set>

namespace mtg::graph::traversal {
namespace {
using Matrix = annot::matrix::MultiIntMatrix;
using Row = Matrix::RowTuples;

size_t base_index(char base) {
    auto index = std::string_view("ACGT").find(base);
    if (index == std::string_view::npos)
        throw std::invalid_argument("ContextCursor requires ACGT sequences");
    return index;
}

std::vector<Row> read_rows(const Matrix &matrix,
                           const std::vector<annot::matrix::BinaryMatrix::Row> &ids) {
    auto rows = matrix.get_row_tuples(ids);
    for (auto &row : rows) {
        // BRWT returns columns in tree order, which need not be label order.
        if (!std::ranges::is_sorted(row, {}, &Row::value_type::first))
            std::ranges::sort(row, {}, &Row::value_type::first);
        for (auto &[column, coords] : row)
            if (!std::ranges::is_sorted(coords))
                std::ranges::sort(coords);
    }
    return rows;
}

bool contains(const Row &row, uint64_t column, uint64_t coordinate) {
    auto it = std::ranges::lower_bound(row, column, {}, &Row::value_type::first);
    return it != row.end() && it->first == column
        && std::ranges::binary_search(it->second, coordinate);
}
} // namespace

ContextCursor::ContextCursor(const AnnotatedDBG &annotated, std::string seed,
                             const std::vector<std::string> &columns,
                             size_t max_occurrences)
      : annotated_(annotated), sequence_(std::move(seed)) {
    const auto &graph = annotated_.get_graph();
    const auto *mapping = annotated_.get_coord_to_header();
    const auto *matrix = dynamic_cast<const Matrix *>(&annotated_.get_annotator().get_matrix());
    if (graph.get_mode() != DeBruijnGraph::BASIC || graph.get_k() < 2
            || graph.alphabet().find_first_not_of("$ACGT") != std::string::npos
            || !matrix || !mapping)
        throw std::invalid_argument("ContextCursor requires BASIC DNA graph (k >= 2), coordinates and .seqs");
    if (!annotated_.check_compatibility()
            || mapping->num_columns() != matrix->get_binary_matrix().num_columns())
        throw std::invalid_argument("Incompatible graph, coordinates and record mapping");
    if (sequence_.size() < graph.get_k())
        throw std::invalid_argument("Context seed must contain at least k bases");
    if (!max_occurrences)
        throw std::invalid_argument("Context occurrence limit must be positive");
    for (char &base : sequence_) {
        base = std::toupper(static_cast<unsigned char>(base));
        base_index(base);
    }
    std::unordered_set<uint64_t> selected;
    for (const auto &column : columns)
        selected.insert(annotated_.get_annotator().get_label_encoder().encode(column));

    std::vector<Node> nodes;
    graph.map_to_nodes_sequentially(sequence_, [&](Node node) { nodes.push_back(node); });
    if (std::ranges::find(nodes, 0) != nodes.end())
        return;
    first_ = nodes.front();
    last_ = nodes.back();
    std::vector<annot::matrix::BinaryMatrix::Row> rows;
    std::vector<size_t> positions;
    VectorMap<Node, size_t> row_positions;
    row_positions.reserve(nodes.size());
    for (auto node : nodes) {
        auto [it, inserted] = row_positions.try_emplace(node, rows.size());
        if (inserted)
            rows.push_back(AnnotatedDBG::graph_to_anno_index(node));
        positions.push_back(it->second);
    }
    const auto tuples = read_rows(*matrix, rows);

    // Decode each distinct k-mer once. Start intersections at the least
    // frequent query position to limit candidate occurrences.
    std::vector<size_t> frequencies(tuples.size());
    for (size_t i = 0; i < tuples.size(); ++i) {
        for (const auto &[column, coords] : tuples[i])
            if (selected.empty() || selected.contains(column))
                frequencies[i] += coords.size();
    }
    size_t anchor = std::ranges::min_element(positions, {},
        [&](size_t row) { return frequencies[row]; }) - positions.begin();
    for (const auto &[column, coords] : tuples[positions[anchor]]) {
        if (!selected.empty() && !selected.contains(column))
            continue;
        for (auto coordinate : coords) {
            if (coordinate >= mapping->num_kmers(column))
                throw std::runtime_error("Coordinate outside indexed record mapping");
            auto [record, offset] = mapping->map_single_coord(column, coordinate);
            auto length = mapping->num_kmers_in_sequence(column, record);
            if (offset < anchor || nodes.size() > length - (offset - anchor))
                continue; // includes apparent matches across adjacent records
            auto start = coordinate - anchor;
            bool match = true;
            for (size_t i = 0; i < positions.size(); ++i) {
                if (!contains(tuples[positions[i]], column, start + i)) { match = false; break; }
            }
            if (match) {
                if (hits_.size() == max_occurrences)
                    throw std::length_error("Context occurrence limit exceeded");
                hits_.push_back({{column, record, offset - anchor}, start, length});
            }
        }
    }
}

ContextCursor::ContextCursor(const ContextCursor &other)
      : annotated_(other.annotated_), sequence_(other.sequence_), first_(other.first_),
        last_(other.last_), hits_(other.hits_), boundaries_ready_(other.boundaries_ready_),
        boundaries_(other.boundaries_), neighbors_(other.neighbors_) {
    // Keep the mutable navigation cache private to each independently used state.
}

ContextCursor::ContextCursor(ContextCursor &&other) noexcept = default;
ContextCursor::~ContextCursor() = default;

std::vector<ContextCursor::Occurrence> ContextCursor::occurrences() const {
    std::vector<Occurrence> result;
    result.reserve(hits_.size());
    for (const auto &hit : hits_)
        result.push_back(hit.occurrence);
    return result;
}

void ContextCursor::find_boundaries(std::optional<Side> side) {
    const bool left = !boundaries_ready_[0] && (!side || *side == Side::LEFT);
    const bool right = !boundaries_ready_[1] && (!side || *side == Side::RIGHT);
    if (!left && !right)
        return;
    boundaries_.resize(hits_.size());
    if (!hits_.empty()) {
        const auto &graph = annotated_.get_graph();
        auto &matrix = dynamic_cast<const Matrix &>(annotated_.get_annotator().get_matrix());
        std::vector<annot::matrix::BinaryMatrix::Row> rows;
        std::vector<std::pair<size_t, char>> edges;
        auto add = [&](size_t side, Node node, char base) {
            if (base == '$')
                return; // unmasked succinct graphs expose dummy edges
            neighbors_[side][base_index(base)] = node;
            rows.push_back(AnnotatedDBG::graph_to_anno_index(node));
            edges.emplace_back(side, base);
        };
        if (left) {
            auto emit = [&](auto n, char b) { add(0, n, b); };
            if (!first_cache_) {
                if (const auto *succinct = dynamic_cast<const DBGSuccinct *>(&graph))
                    first_cache_ = std::make_unique<NodeFirstCache>(*succinct, 4096);
            }
            if (first_cache_)
                first_cache_->call_incoming_kmers(first_, emit);
            else
                graph.call_incoming_kmers(first_, emit);
        }
        if (right)
            graph.call_outgoing_kmers(last_, [&](auto n, char b) { add(1, n, b); });
        auto tuples = read_rows(matrix, rows);
        auto kmers = sequence_.size() - graph.get_k() + 1;
        for (size_t i = 0; i < hits_.size(); ++i) {
            const auto &hit = hits_[i];
            for (size_t e = 0; e < edges.size(); ++e) {
                auto [side, base] = edges[e];
                if ((!side && !hit.occurrence.offset)
                        || (side && kmers == hit.record_kmers - hit.occurrence.offset))
                    continue;
                auto coordinate = side ? hit.global + kmers : hit.global - 1;
                if (contains(tuples[e], hit.occurrence.column, coordinate)) {
                    if (boundaries_[i][side])
                        throw std::runtime_error("Multiple bases at one source coordinate");
                    boundaries_[i][side] = base;
                }
            }
            // An ambiguous/unindexed source k-mer can stop an otherwise valid
            // record. Do not misreport this as a record end or invent a base.
            if (left && hit.occurrence.offset && !boundaries_[i][0])
                boundaries_[i][0] = '?';
            if (right && kmers < hit.record_kmers - hit.occurrence.offset && !boundaries_[i][1])
                boundaries_[i][1] = '?';
        }
    }
    boundaries_ready_[0] |= left;
    boundaries_ready_[1] |= right;
}

std::vector<ContextCursor::Choice> ContextCursor::choices() {
    find_boundaries();
    std::map<std::pair<char, char>, size_t> counts;
    for (const auto &boundary : boundaries_)
        ++counts[{boundary[0], boundary[1]}];
    std::vector<Choice> result;
    for (const auto &[pair, count] : counts)
        result.push_back({pair.first, pair.second, count});
    return result;
}

bool ContextCursor::extend(Side side, char base) {
    base = std::toupper(static_cast<unsigned char>(base));
    auto index = base_index(base);
    find_boundaries(side);
    auto direction = side == Side::LEFT ? 0 : 1;
    if (std::ranges::none_of(boundaries_, [&](const auto &b) { return b[direction] == base; }))
        return false;

    // Grow the string before changing hits, so allocation failure leaves the
    // occurrence set intact. Compact in place: closure often keeps every hit.
    if (!direction)
        sequence_.insert(sequence_.begin(), base);
    else
        sequence_.push_back(base);
    size_t retained = 0;
    for (size_t i = 0; i < hits_.size(); ++i) {
        if (boundaries_[i][direction] == base) {
            auto hit = hits_[i];
            if (!direction) { --hit.global; --hit.occurrence.offset; }
            hits_[retained] = hit;
            auto boundary = boundaries_[i];
            boundary[direction] = 0;
            boundaries_[retained++] = boundary;
        }
    }
    hits_.resize(retained);
    boundaries_.resize(retained);
    if (!direction)
        first_ = neighbors_[0][index];
    else
        last_ = neighbors_[1][index];
    boundaries_ready_[direction] = false;
    neighbors_[direction] = {};
    return true;
}

ContextCursor::Closure ContextCursor::close(size_t max_bases_per_side) {
    std::array<size_t, 2> added{};
    std::array<StopReason, 2> stops{};
    for (size_t side = 0; side < 2; ++side) {
        while (true) {
            if (hits_.empty()) { stops[side] = StopReason::NO_OCCURRENCES; break; }
            if (added[side] == max_bases_per_side) { stops[side] = StopReason::RADIUS_LIMIT; break; }
            find_boundaries(side ? Side::RIGHT : Side::LEFT);
            auto boundary = [side](const auto &b) { return b[side]; };
            if (std::ranges::find(boundaries_, '\0', boundary) != boundaries_.end()) {
                stops[side] = StopReason::RECORD_END;
                break;
            }
            if (std::ranges::find(boundaries_, '?', boundary) != boundaries_.end()) {
                stops[side] = StopReason::UNINDEXED;
                break;
            }
            char base = boundaries_.front()[side];
            if (!std::ranges::all_of(boundaries_, [&](const auto &b) { return b[side] == base; })) {
                stops[side] = StopReason::BRANCH;
                break;
            }
            extend(side ? Side::RIGHT : Side::LEFT, base);
            ++added[side];
        }
    }
    return {added[0], added[1], stops[0], stops[1]};
}

} // namespace mtg::graph::traversal
