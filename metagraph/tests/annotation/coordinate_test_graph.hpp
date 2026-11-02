#ifndef __COORDINATE_TEST_GRAPH_HPP__
#define __COORDINATE_TEST_GRAPH_HPP__

#include <numeric>
#include <stdexcept>

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/binary_matrix/multi_brwt/brwt_builders.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

namespace mtg::test {
struct CoordinateRecord { std::string sequence, column, header; };

// Small real coordinate-matrix fixture, with adjacent records in each column.
inline std::unique_ptr<graph::AnnotatedDBG>
coordinate_test_graph(size_t k, const std::vector<CoordinateRecord> &records,
                      bool row_diff = false, bool brwt = false) {
    using namespace annot::matrix;
    auto graph = std::make_shared<graph::DBGSuccinct>(k);
    annot::LabelEncoder<> encoder;
    for (const auto &record : records) {
        graph->add_sequence(record.sequence);
        encoder.insert_and_encode(record.column);
    }
    graph->mask_dummy_kmers(1, false);
    graph->switch_state(graph::boss::BOSS::STAT);
    std::vector<std::vector<std::vector<uint64_t>>> coords(
        encoder.size(), std::vector<std::vector<uint64_t>>(graph->max_index()));
    std::vector<std::vector<std::string>> headers(encoder.size());
    std::vector<std::vector<uint64_t>> counts(encoder.size());
    std::vector<uint64_t> offsets(encoder.size());
    for (const auto &record : records) {
        auto col = encoder.encode(record.column);
        size_t pos = 0;
        graph->map_to_nodes_sequentially(record.sequence, [&](auto node) {
            if (node)
                coords[col][node - 1].push_back(offsets[col] + pos);
            ++pos;
        });
        offsets[col] += pos;
        headers[col].push_back(record.header);
        counts[col].push_back(pos);
    }
    sdsl::bit_vector anchors(graph->max_index(), 1);
    std::unique_ptr<decltype(coords)> expected;
    if (row_diff) {
        expected = std::make_unique<decltype(coords)>(coords);
        const auto &original = *expected;
        graph->call_nodes([&](auto node) {
            if (!graph->has_single_outgoing(node))
                return;
            uint64_t next = 0;
            graph->adjacent_outgoing_nodes(node, [&](auto n) { next = n; });
            if (next <= node)
                return; // strictly increasing paths cannot form anchor-free cycles
            anchors[node - 1] = 0;
            for (size_t col = 0; col < coords.size(); ++col) {
                auto current = original[col][node - 1];
                for (auto &coordinate : current)
                    ++coordinate;
                auto &diff = coords[col][node - 1];
                diff.clear();
                const auto &successor = original[col][next - 1];
                std::set_symmetric_difference(current.begin(), current.end(),
                                              successor.begin(), successor.end(),
                                              std::back_inserter(diff));
            }
        });
    }
    std::vector<std::unique_ptr<bit_vector>> bits;
    std::vector<bit_vector_smart> delimiters;
    std::vector<sdsl::int_vector<>> values;
    for (const auto &column : coords) {
        sdsl::bit_vector membership(column.size(), 0);
        std::vector<uint64_t> flat;
        std::vector<bool> marks{true};
        for (size_t row = 0; row < column.size(); ++row) {
            if (column[row].empty())
                continue;
            membership[row] = 1;
            for (auto coord : column[row]) {
                flat.push_back(coord);
                marks.push_back(false);
            }
            marks.push_back(true);
        }
        sdsl::bit_vector delim(marks.size(), 0);
        for (size_t i = 0; i < marks.size(); ++i)
            delim[i] = marks[i];
        sdsl::int_vector<> vals(flat.size(), 0, 64);
        std::copy(flat.begin(), flat.end(), vals.begin());
        bits.push_back(std::make_unique<bit_vector_sd>(membership));
        delimiters.emplace_back(delim);
        values.push_back(std::move(vals));
    }
    using Annotator = graph::AnnotatedDBG::Annotator;
    auto wrap = [&](auto matrix) -> std::unique_ptr<Annotator> {
        using Base = std::remove_reference_t<decltype(*matrix)>;
        if (row_diff) {
            struct MutableRowDiff : TupleRowDiff<Base> {
                using TupleRowDiff<Base>::TupleRowDiff;
                void set_anchors(const sdsl::bit_vector &anchors) {
                    this->anchor_ = bit_vector_small(anchors);
                }
            };
            auto diffs = std::make_unique<MutableRowDiff>(graph.get(), std::move(*matrix));
            diffs->set_anchors(anchors);
            return std::make_unique<annot::StaticBinRelAnnotator<TupleRowDiff<Base>>>(
                std::move(diffs), encoder);
        }
        return std::make_unique<annot::StaticBinRelAnnotator<Base>>(std::move(matrix), encoder);
    };
    std::unique_ptr<Annotator> annotation;
    if (brwt) {
        auto reversed = [](const BRWTBottomUpBuilder::VectorPtrs &columns) {
            BRWTBottomUpBuilder::Partition partition(1);
            for (size_t i = columns.size(); i; --i)
                partition.front().push_back(i - 1);
            return partition;
        };
        annotation = wrap(std::make_unique<TupleCSCMatrix<BRWT>>(
            BRWTBottomUpBuilder::build(std::move(bits), reversed),
            std::move(delimiters), std::move(values)));
    } else {
        annotation = wrap(std::make_unique<TupleCSCMatrix<ColumnMajor>>(
            ColumnMajor(std::move(bits)), std::move(delimiters), std::move(values)));
    }
    if (expected) {
        std::vector<BinaryMatrix::Row> ids;
        graph->call_nodes([&](auto node) { ids.push_back(node - 1); });
        auto decoded = dynamic_cast<const MultiIntMatrix &>(annotation->get_matrix()).get_row_tuples(ids);
        for (size_t index = 0; index < decoded.size(); ++index) {
            const auto row = ids[index];
            size_t entry = 0;
            for (size_t col = 0; col < expected->size(); ++col) {
                const auto &coords = (*expected)[col][row];
                if (coords.empty())
                    continue;
                if (entry >= decoded[index].size() || decoded[index][entry].first != col
                        || std::vector<uint64_t>(decoded[index][entry].second.begin(),
                                                  decoded[index][entry].second.end()) != coords)
                    throw std::logic_error("Coordinate fixture row-diff roundtrip failed");
                ++entry;
            }
            if (entry != decoded[index].size())
                throw std::logic_error("Unexpected coordinate in row-diff fixture");
        }
    }
    return std::make_unique<graph::AnnotatedDBG>(
        graph, std::move(annotation), false,
        std::make_unique<annot::CoordToHeader>(std::move(headers), std::move(counts)));
}
} // namespace mtg::test
#endif
