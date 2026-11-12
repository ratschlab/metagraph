#include <algorithm>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include <spdlog/fmt/fmt.h>
#include "gtest/gtest.h"

#include "coordinate_test_graph.hpp"
#include "dense_search_values.hpp"
#include "cli/query.hpp"
#include "graph/alignment/score_kmer_presence_mask.hpp"


// The values per query position of a search's rows, held where a row matches (PresenceRuns,
// CountRuns, PositionCoords; train 5a X2): they stand for the vectors of the query's k-mers the
// queries returned before, and every output written from them (the x/o mask, the counts' runs,
// the coordinate ranges, the BIGSI score, the full strings of --verbose-output) is the one the
// vectors gave, byte for byte.

namespace mtg::test {
namespace {

using graph::PresenceRuns;
using graph::CountRuns;
using graph::PositionCoords;

// The writers of the vectors as the queries returned them (cli/query.cpp before the runs):
// the counts' runs "i=v" / "i-j=v" of one nonzero value
std::vector<std::string> dense_count_runs(const std::vector<size_t> &abundances) {
    std::vector<std::string> run_strings;
    std::pair<size_t, size_t> last(0, abundances.at(0));
    for (size_t i = 1; i <= abundances.size(); ++i) {
        if (i < abundances.size() && abundances[i] == last.second)
            continue; // extend run

        // end of run
        if (last.second) {
            if (i == last.first + 1) {
                run_strings.push_back(fmt::format("{}={}", last.first, last.second));
            } else {
                run_strings.push_back(fmt::format("{}-{}={}", last.first, i - 1, last.second));
            }
        }

        // start new run
        if (i < abundances.size()) {
            last.first = i;
            last.second = abundances[i];
        }
    }
    return run_strings;
}

// Masks of every kind the runs must reproduce: random bits of a density, runs of random
// lengths (long ones covering whole 64-bit words of ones or zeros), words of one bit
sdsl::bit_vector random_mask(std::mt19937 &rng, size_t size, int kind) {
    sdsl::bit_vector mask(size, false);
    if (kind < 6) {
        const double density[] = { 0.0, 0.02, 0.3, 0.5, 0.98, 1.0 };
        std::bernoulli_distribution bit(density[kind]);
        for (size_t i = 0; i < size; ++i) {
            mask[i] = bit(rng);
        }
    } else if (kind == 6) {
        std::geometric_distribution<size_t> length(0.02);
        bool bit = std::bernoulli_distribution(0.5)(rng);
        for (size_t i = 0; i < size; bit = !bit) {
            size_t end = std::min(size, i + 1 + length(rng));
            for ( ; i < end; ++i) {
                mask[i] = bit;
            }
        }
    } else {
        // whole words, with single flips between them
        for (size_t i = 0; i < size; ++i) {
            mask[i] = (i / 64) % 2 ? (i % 64 != 5) : (i % 64 == 63);
        }
    }
    return mask;
}

TEST(AnnotatedSearchPositions, PresenceRunsStandForTheirMask) {
    std::mt19937 rng(42);
    std::vector<size_t> sizes;
    for (size_t size = 0; size <= 300; ++size) {
        sizes.push_back(size);
    }
    for (size_t size : { 511, 512, 513, 1000, 4096, 5000, 20000 }) {
        sizes.push_back(size);
    }
    for (size_t size : sizes) {
        for (int kind = 0; kind < 8; ++kind) {
            const sdsl::bit_vector mask = random_mask(rng, size, kind);
            const PresenceRuns runs(mask);
            SCOPED_TRACE(::testing::Message() << size << " " << kind);
            ASSERT_EQ(size, runs.size);
            ASSERT_EQ(mask, runs.to_dense());
            ASSERT_EQ(sdsl::util::cnt_one_bits(mask), runs.count());
            // canonical: none empty, ascending, no two adjacent
            for (size_t r = 0; r < runs.runs.size(); ++r) {
                ASSERT_LT(runs.runs[r].first, runs.runs[r].second);
                if (r) {
                    ASSERT_LT(runs.runs[r - 1].second, runs.runs[r].first);
                }
            }
            // the outputs: the x/o mask, the full string, the score for every k
            ASSERT_EQ(cli::encode_presence_mask(mask), cli::encode_presence_mask(runs));
            ASSERT_EQ(sdsl::util::to_string(mask), sdsl::util::to_string(runs.to_dense()));
            for (size_t k : { 1, 2, 3, 4, 5, 11, 31, 63 }) {
                ASSERT_EQ(graph::align::score_kmer_presence_mask(k, mask),
                          graph::align::score_kmer_presence_mask(k, runs.size, runs.runs))
                    << "k " << k;
            }
        }
    }
}

TEST(AnnotatedSearchPositions, PresenceRunsAddedInOrder) {
    PresenceRuns runs(10);
    runs.add(0);
    runs.add(1);
    runs.add(3);
    runs.add_run(4, 6);
    runs.add_run(6, 6);
    runs.add(9);
    EXPECT_EQ((std::vector<std::pair<size_t, size_t>>{ { 0, 2 }, { 3, 6 }, { 9, 10 } }),
              runs.runs);
    EXPECT_EQ(6u, runs.count());
    EXPECT_EQ("x2o1x3o3x1", cli::encode_presence_mask(runs));
    EXPECT_EQ("o7", cli::encode_presence_mask(PresenceRuns(7)));
    EXPECT_EQ("", cli::encode_presence_mask(PresenceRuns()));
    EXPECT_EQ(0, graph::align::score_kmer_presence_mask(31, 0, {}));
}

TEST(AnnotatedSearchPositions, CountRunsStandForTheirValues) {
    std::mt19937 rng(43);
    for (size_t size = 1; size <= 400; ++size) {
        for (int kind = 0; kind < 4; ++kind) {
            std::vector<size_t> values(size, 0);
            if (kind == 0) {
                // a few values, each position on its own
                std::uniform_int_distribution<size_t> value(0, 3);
                for (auto &v : values) {
                    v = value(rng);
                }
            } else {
                // runs of one value, zeros between some
                std::geometric_distribution<size_t> length(kind == 1 ? 0.3 : 0.03);
                std::uniform_int_distribution<size_t> value(0, kind == 3 ? 1'000'000 : 2);
                for (size_t i = 0; i < size; ) {
                    const size_t v = value(rng);
                    for (size_t end = std::min(size, i + 1 + length(rng)); i < end; ++i) {
                        values[i] = v;
                    }
                }
            }
            const CountRuns runs(values);
            SCOPED_TRACE(::testing::Message() << size << " " << kind);
            ASSERT_EQ(values, runs.to_dense());
            ASSERT_EQ(static_cast<size_t>(std::count_if(values.begin(), values.end(),
                                                        [](size_t v) { return v > 0; })),
                      runs.count());
            for (const auto &run : runs.runs) {
                ASSERT_LT(run.begin, run.end);
                ASSERT_LT(0u, run.value);
            }
            ASSERT_EQ(dense_count_runs(values), cli::count_run_strings(runs));
        }
    }
    // runs added position by position merge where adjacent and of one value; a zero is no value
    CountRuns runs(8);
    runs.add(0, 2);
    runs.add(1, 2);
    runs.add(2, 0);
    runs.add(3, 2);
    runs.add(4, 5);
    runs.add(5, 5);
    runs.add(7, 1);
    EXPECT_EQ((std::vector<CountRuns::Run>{ { 0, 2, 2 }, { 3, 4, 2 }, { 4, 6, 5 }, { 7, 8, 1 } }),
              runs.runs);
    EXPECT_EQ((std::vector<std::string>{ "0-1=2", "3=2", "4-5=5", "7=1" }),
              cli::count_run_strings(runs));
}

TEST(AnnotatedSearchPositions, PositionCoordsStandForTheirCoordinates) {
    std::mt19937 rng(44);
    for (size_t size = 0; size <= 200; ++size) {
        for (int kind = 0; kind < 3; ++kind) {
            // coordinates that continue from one position to the next (ranges), with gaps,
            // new starts and several per position
            std::vector<SmallVector<uint64_t>> values(size);
            std::bernoulli_distribution present(kind == 0 ? 0.2 : (kind == 1 ? 0.7 : 1.0));
            std::bernoulli_distribution restart(0.1);
            std::uniform_int_distribution<size_t> how_many(1, 3);
            std::uniform_int_distribution<uint64_t> coord(0, 50);
            std::vector<uint64_t> last;
            for (size_t i = 0; i < size; ++i) {
                if (!present(rng)) {
                    last.clear();
                    continue;
                }
                std::vector<uint64_t> next;
                if (last.empty() || restart(rng)) {
                    for (size_t n = how_many(rng); n; --n) {
                        next.push_back(coord(rng));
                    }
                } else {
                    for (uint64_t c : last) {
                        next.push_back(c + 1);
                    }
                }
                std::sort(next.begin(), next.end());
                next.erase(std::unique(next.begin(), next.end()), next.end());
                values[i].assign(next.begin(), next.end());
                last = next;
            }
            const PositionCoords coords(values);
            SCOPED_TRACE(::testing::Message() << size << " " << kind);
            ASSERT_EQ(values, coords.to_dense());
            ASSERT_EQ(static_cast<size_t>(std::count_if(values.begin(), values.end(),
                                                        [](const auto &v) { return !v.empty(); })),
                      coords.count());
            ASSERT_EQ(cli::collapse_coord_ranges(values), cli::collapse_coord_ranges(coords));
        }
    }
    // coordinates added one by one: a position's in their order, a new position after it
    PositionCoords coords(6);
    coords.add(1, 7);
    coords.add(1, 9);
    coords.add(2, 8);
    coords.add(2, 10);
    coords.add(4, 3);
    EXPECT_EQ((std::vector<std::pair<size_t, size_t>>{ { 1, 3 }, { 4, 5 } }), coords.runs);
    EXPECT_EQ((std::vector<SmallVector<uint64_t>>{ {}, { 7, 9 }, { 8, 10 }, {}, { 3 }, {} }),
              coords.to_dense());
    EXPECT_EQ((std::vector<std::string>{ "1-7-8", "1-9-10", "4-3" }),
              cli::collapse_coord_ranges(coords));
}

// The records and query positions of a coordinate index, by brute force: the occurrences of
// each query k-mer in each record (a k-mer with a character outside ACGT is in no record)
struct Occurrences {
    std::string label;
    // per query position, the k-mer's positions in the record (or in the column's records
    // concatenated, without a record mapping)
    std::vector<std::vector<uint64_t>> coords;
    size_t kmers_in_target;

    size_t matched() const {
        return std::count_if(coords.begin(), coords.end(),
                             [](const auto &c) { return !c.empty(); });
    }
};

std::vector<Occurrences> occurrences(size_t k, const std::vector<CoordinateRecord> &records,
                                     const std::string &query, bool by_record) {
    const size_t num_kmers = query.size() < k ? 0 : query.size() - k + 1;
    std::vector<Occurrences> result;
    std::vector<std::string> columns;
    std::vector<uint64_t> offsets;
    for (const auto &record : records) {
        size_t c = std::find(columns.begin(), columns.end(), record.column) - columns.begin();
        if (c == columns.size()) {
            columns.push_back(record.column);
            offsets.push_back(0);
            if (!by_record)
                result.push_back({ record.column, std::vector<std::vector<uint64_t>>(num_kmers), 0 });
        }
        const size_t record_kmers = record.sequence.size() - k + 1;
        if (by_record)
            result.push_back({ record.header, std::vector<std::vector<uint64_t>>(num_kmers), record_kmers });
        auto &entry = by_record ? result.back() : result[c];
        const uint64_t offset = by_record ? 0 : offsets[c];
        for (size_t p = 0; p < num_kmers; ++p) {
            const std::string_view kmer(query.data() + p, k);
            if (kmer.find_first_not_of("ACGT") != std::string_view::npos)
                continue;
            for (size_t j = 0; j < record_kmers; ++j) {
                if (std::string_view(record.sequence.data() + j, k) == kmer)
                    entry.coords[p].push_back(offset + j);
            }
        }
        offsets[c] += record_kmers;
    }
    result.erase(std::remove_if(result.begin(), result.end(),
                                [](const auto &o) { return !o.matched(); }),
                 result.end());
    return result;
}

std::string random_sequence(std::mt19937 &rng, size_t length) {
    std::string s(length, 'A');
    std::uniform_int_distribution<int> base(0, 3);
    for (char &c : s) {
        c = "ACGT"[base(rng)];
    }
    return s;
}

// Every row of the three per-position queries, at every cut, holds the values of its record
// (its column without a record mapping) at every query position: the matched positions, the
// occurrences and the coordinates, by brute force; at discovery_fraction 1 every position
void check_positions(const graph::AnnotatedDBG &index, size_t k,
                     const std::vector<CoordinateRecord> &records, const std::string &query,
                     bool by_record) {
    const auto expected = occurrences(k, records, query, by_record);
    const size_t num_kmers = query.size() < k ? 0 : query.size() - k + 1;
    // a row returned is one of the expected rows of its label
    auto find = [&](const std::string &label, const auto &matches) {
        size_t found = 0;
        for (const auto &e : expected) {
            if (e.label == label && matches(e))
                ++found;
        }
        return found;
    };
    for (size_t top : { size_t(1), size_t(2), std::numeric_limits<size_t>::max() }) {
        SCOPED_TRACE(::testing::Message() << query << " top " << top);
        const auto coords = index.get_kmer_coordinates(query, top, 0, 0);
        const auto counts = index.get_kmer_counts(query, top, 0, 0);
        const auto masks = index.get_top_label_signatures(query, top, 0, 0);
        ASSERT_EQ(std::min(top, expected.size()), coords.size());
        ASSERT_EQ(coords.size(), counts.size());
        ASSERT_EQ(coords.size(), masks.size());
        for (size_t r = 0; r < coords.size(); ++r) {
            const auto &[label, count, positions, length] = coords[r];
            ASSERT_EQ(num_kmers, positions.size);
            auto dense_coords = positions.to_dense();
            for (auto &c : dense_coords) {
                std::sort(c.begin(), c.end());
            }
            EXPECT_LE(1u, find(label, [&](const Occurrences &e) {
                if (e.matched() != count || e.kmers_in_target != length)
                    return false;
                for (size_t p = 0; p < num_kmers; ++p) {
                    if (std::vector<uint64_t>(dense_coords[p].begin(), dense_coords[p].end())
                            != e.coords[p])
                        return false;
                }
                return true;
            })) << label;
        }
        for (const auto &[label, count, values] : counts) {
            ASSERT_EQ(num_kmers, values.size);
            const auto dense_counts = values.to_dense();
            EXPECT_LE(1u, find(label, [&](const Occurrences &e) {
                if (e.matched() != count)
                    return false;
                for (size_t p = 0; p < num_kmers; ++p) {
                    if (dense_counts[p] != e.coords[p].size())
                        return false;
                }
                return true;
            })) << label;
        }
        for (const auto &[label, count, mask] : masks) {
            ASSERT_EQ(num_kmers, mask.size);
            const auto dense_mask = mask.to_dense();
            EXPECT_LE(1u, find(label, [&](const Occurrences &e) {
                if (e.matched() != count)
                    return false;
                for (size_t p = 0; p < num_kmers; ++p) {
                    if (dense_mask[p] != !e.coords[p].empty())
                        return false;
                }
                return true;
            })) << label;
        }
    }
    // discovery_fraction 1: the rows matching every position, each marked everywhere
    const auto all = index.get_top_label_signatures(query, std::numeric_limits<size_t>::max(), 1.0);
    for (const auto &[label, count, mask] : all) {
        EXPECT_EQ(num_kmers, count);
        EXPECT_EQ(PresenceRuns(sdsl::bit_vector(num_kmers, true)), mask);
    }
}

TEST(AnnotatedSearchPositions, RowsHoldTheirRecordsValuesAtEveryPosition) {
    std::mt19937 rng(45);
    for (size_t round = 0; round < 12; ++round) {
        const size_t k = 3 + round % 3;
        // short records on few columns, so that k-mers repeat within and across records
        std::vector<CoordinateRecord> records;
        std::uniform_int_distribution<size_t> length(k, 40);
        const std::vector<std::string> column_names { "F0", "F1", "F2" };
        for (size_t i = 0; i < 8; ++i) {
            records.push_back({ random_sequence(rng, length(rng)), column_names[i % 3],
                                // one header in two columns
                                i == 5 ? std::string("r1") : "r" + std::to_string(i) });
        }
        // queries: a record, pieces of records joined, with a foreign k-mer and an N, a
        // repeated record (positions sharing a row), a random sequence
        std::vector<std::string> queries {
            records[0].sequence,
            records[1].sequence.substr(0, k + 2) + "N" + records[2].sequence + records[1].sequence,
            records[3].sequence + records[3].sequence + records[3].sequence,
            records[4].sequence + random_sequence(rng, 2 * k) + records[6].sequence,
            random_sequence(rng, 60),
            records[7].sequence.substr(0, k)
        };
        for (bool row_diff : { false, true }) {
            for (bool brwt : { false, true }) {
                for (bool by_record : { true, false }) {
                    SCOPED_TRACE(::testing::Message() << "round " << round << " row_diff "
                                 << row_diff << " brwt " << brwt << " by_record " << by_record);
                    auto index = coordinate_test_graph(k, records, row_diff, brwt, false, by_record);
                    for (const auto &query : queries) {
                        check_positions(*index, k, records, query, by_record);
                    }
                }
            }
        }
    }
}

// A query of many records' first k-mers, each record matched at one position: the rows
// returned times the query's k-mers is large, the positions matched are not (one per row).
// The values are held at the matched positions only: one run per row. (A record's first base
// differs from the base after the previous record's first k-mer, so that the k-mers across
// the junction do not continue the previous record.)
TEST(AnnotatedSearchPositions, ManyRecordsOneMatchEach) {
    std::mt19937 rng(46);
    const size_t k = 21;
    std::vector<CoordinateRecord> records;
    std::string query;
    for (size_t i = 0; i < 400; ++i) {
        std::string sequence = random_sequence(rng, 40);
        if (i && sequence[0] == records.back().sequence[k])
            sequence[0] = sequence[0] == 'A' ? 'C' : 'A';
        records.push_back({ sequence, "F" + std::to_string(i % 4), "r" + std::to_string(i) });
        query += records.back().sequence.substr(0, k);
    }
    auto index = coordinate_test_graph(k, records);
    const size_t num_kmers = query.size() - k + 1;
    const size_t top = std::numeric_limits<size_t>::max();
    const auto masks = index->get_top_label_signatures(query, top, 0, 0);
    const auto counts = index->get_kmer_counts(query, top, 0, 0);
    const auto coords = index->get_kmer_coordinates(query, top, 0, 0);
    ASSERT_EQ(records.size(), masks.size());
    ASSERT_EQ(records.size(), counts.size());
    ASSERT_EQ(records.size(), coords.size());
    for (size_t r = 0; r < records.size(); ++r) {
        const auto &[label, count, mask] = masks[r];
        const size_t i = std::stoul(label.substr(1));
        EXPECT_EQ(1u, count);
        EXPECT_EQ(num_kmers, mask.size);
        ASSERT_EQ(1u, mask.runs.size());
        EXPECT_EQ(std::make_pair(i * k, i * k + 1), mask.runs[0]);
        const std::string before = i ? fmt::format("o{}", i * k) : "";
        const size_t after = num_kmers - i * k - 1;
        EXPECT_EQ(before + "x1" + (after ? fmt::format("o{}", after) : ""),
                  cli::encode_presence_mask(mask));
        EXPECT_EQ(cli::encode_presence_mask(mask.to_dense()), cli::encode_presence_mask(mask));
        EXPECT_EQ(graph::align::score_kmer_presence_mask(k, mask.to_dense()),
                  graph::align::score_kmer_presence_mask(k, mask.size, mask.runs));

        const auto &[count_label, matched, values] = counts[r];
        EXPECT_EQ(label, count_label);
        ASSERT_EQ(1u, values.runs.size());
        EXPECT_EQ((CountRuns::Run{ i * k, i * k + 1, 1 }), values.runs[0]);
        EXPECT_EQ(std::vector<std::string>{ fmt::format("{}=1", i * k) },
                  cli::count_run_strings(values));

        const auto &[coord_label, coord_count, positions, length] = coords[r];
        EXPECT_EQ(label, coord_label);
        EXPECT_EQ(40 - k + 1, length);
        ASSERT_EQ(1u, positions.count());
        EXPECT_EQ(std::vector<std::string>{ fmt::format("{}-0", i * k) },
                  cli::collapse_coord_ranges(positions));
    }
}

// The rows a SeqSearchResult writes (the server's /search body, the CLI's JSON and text) from
// the values held at the matched positions are the ones the vectors gave: the writers of the
// vectors (encode_presence_mask, the score, the counts' runs, collapse_coord_ranges of the
// vectors and the full strings of --verbose-output) applied to to_dense()
TEST(AnnotatedSearchPositions, WrittenAsTheirVectors) {
    std::mt19937 rng(47);
    const size_t k = 31;
    for (size_t size : { 1, 2, 3, 63, 64, 65, 200, 1000 }) {
        cli::SeqSearchResult::LabelSigVec signatures;
        cli::SeqSearchResult::LabelCountAbundancesVec counts;
        cli::SeqSearchResult::LabelCountCoordsVec coords;
        for (int kind = 0; kind < 8; ++kind) {
            const std::string label = "label" + std::to_string(kind);
            const sdsl::bit_vector mask = random_mask(rng, size, kind);
            if (!sdsl::util::cnt_one_bits(mask))
                continue;
            std::vector<size_t> values(size, 0);
            std::vector<SmallVector<uint64_t>> positions(size);
            std::uniform_int_distribution<size_t> value(1, 3);
            for (size_t i = 0; i < size; ++i) {
                if (mask[i]) {
                    values[i] = value(rng);
                    positions[i] = { 10 + i, 500 + 2 * i };
                }
            }
            const size_t matched = sdsl::util::cnt_one_bits(mask);
            signatures.emplace_back(label, matched, PresenceRuns(mask));
            counts.emplace_back(label, matched, CountRuns(values));
            coords.emplace_back(label, matched, PositionCoords(positions), kind % 2 ? 0 : 5000);
        }
        auto write = [&](auto &&rows) {
            return cli::SeqSearchResult(cli::QuerySequence { 0, "q", "" }, std::move(rows));
        };
        const auto sig_result = write(cli::SeqSearchResult::LabelSigVec(signatures));
        const auto count_result = write(cli::SeqSearchResult::LabelCountAbundancesVec(counts));
        const auto coord_result = write(cli::SeqSearchResult::LabelCountCoordsVec(coords));
        for (bool verbose : { false, true }) {
            SCOPED_TRACE(::testing::Message() << size << " verbose " << verbose);
            const Json::Value sig_json = sig_result.to_json(verbose, k);
            const Json::Value count_json = count_result.to_json(verbose, k);
            const Json::Value coord_json = coord_result.to_json(verbose, k);
            std::string sig_text = "0\tq";
            std::string count_text = "0\tq";
            std::string coord_text = "0\tq";
            for (Json::ArrayIndex r = 0; r < signatures.size(); ++r) {
                const auto &[label, count, mask_runs] = signatures[r];
                const sdsl::bit_vector mask = mask_runs.to_dense();
                const std::string mask_string = verbose ? sdsl::util::to_string(mask)
                                                        : cli::encode_presence_mask(mask);
                const int32_t score = graph::align::score_kmer_presence_mask(k, mask);
                EXPECT_EQ(mask_string, sig_json["results"][r]["signature"]["presence_mask"].asString());
                EXPECT_EQ(score, sig_json["results"][r]["signature"]["score"].asInt());
                EXPECT_EQ(count, sig_json["results"][r]["kmer_count"].asUInt64());
                sig_text += fmt::format("\t<{}>:{}:{}:{}", label, count, mask_string, score);

                const std::vector<size_t> values = std::get<2>(counts[r]).to_dense();
                std::string abundances;
                if (verbose) {
                    for (size_t i = 0; i < values.size(); ++i) {
                        abundances += (i ? ":" : "") + std::to_string(values[i]);
                    }
                } else {
                    const auto runs = dense_count_runs(values);
                    for (size_t i = 0; i < runs.size(); ++i) {
                        abundances += (i ? ":" : "") + runs[i];
                    }
                }
                EXPECT_EQ(abundances, count_json["results"][r]["kmer_abundances"].asString());
                count_text += "\t<" + label + ">:" + abundances;

                const auto &[coord_label, coord_count, coord_values, length] = coords[r];
                const auto tuples = coord_values.to_dense();
                std::string written;
                if (verbose) {
                    for (size_t i = 0; i < tuples.size(); ++i) {
                        written += i ? ":" : "";
                        for (size_t j = 0; j < tuples[i].size(); ++j) {
                            written += (j ? "," : "") + std::to_string(tuples[i][j]);
                        }
                    }
                } else {
                    const auto ranges = cli::collapse_coord_ranges(tuples);
                    for (size_t i = 0; i < ranges.size(); ++i) {
                        written += (i ? ":" : "") + ranges[i];
                    }
                }
                EXPECT_EQ(written, coord_json["results"][r]["kmer_coords"].asString());
                EXPECT_EQ(length != 0, coord_json["results"][r].isMember("kmers_in_target"));
                coord_text += "\t<" + coord_label + ">"
                        + (length ? "/" + std::to_string(length) : "") + ":" + written;
            }
            EXPECT_EQ(sig_text, sig_result.to_string(";", false, verbose, k));
            EXPECT_EQ(count_text, count_result.to_string(";", false, verbose, k));
            EXPECT_EQ(coord_text, coord_result.to_string(";", false, verbose, k));
        }
    }
}

} // namespace
} // namespace mtg::test
