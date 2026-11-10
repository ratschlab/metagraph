#include "annotation_remap.hpp"

#include <atomic>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <stdexcept>

#include <sdsl/int_vector_buffer.hpp>

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/logger.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace annot {

using graph::DeBruijnGraph;
using node_index = DeBruijnGraph::node_index;
using mtg::common::logger;

static uint64_t width_for(uint64_t max_value) {
    return max_value ? sdsl::bits::hi(max_value) + 1 : 1;
}

// Walk `source`, look every k-mer up in `target`, and pass each matched pair to
// store(source_row, target_node). Throws if the graphs have different k-mers.
template <class Store>
static void match_rows(const DeBruijnGraph &source,
                       const DeBruijnGraph &target,
                       size_t num_threads,
                       const Store &store) {
    if (source.get_mode() != DeBruijnGraph::BASIC
            || target.get_mode() != DeBruijnGraph::BASIC)
        throw std::invalid_argument("Row remapping supports only graphs in basic mode");

    if (source.get_k() != target.get_k())
        throw std::invalid_argument("Row remapping needs graphs with the same k");

    sdsl::bit_vector found(target.max_index(), false);

    std::atomic<uint64_t> num_mapped(0);
    std::atomic<uint64_t> num_missing(0);
    std::atomic<uint64_t> num_repeated(0);

    source.call_sequences([&](const std::string &sequence, const std::vector<node_index> &path) {
        uint64_t mapped = 0;
        uint64_t missing = 0;
        uint64_t repeated = 0;
        size_t i = 0;
        target.map_to_nodes_sequentially(sequence, [&](node_index target_node) {
            node_index source_node = path.at(i++);
            if (!target_node) {
                ++missing;
                return;
            }
            if (fetch_and_set_bit(found.data(), target_node - 1, true,
                                  std::memory_order_relaxed)) {
                ++repeated;
                return;
            }
            store(source_node - 1, target_node);
            ++mapped;
        });
        num_mapped += mapped;
        num_missing += missing;
        num_repeated += repeated;
    }, num_threads);

    if (num_missing || num_repeated || num_mapped != target.num_nodes()) {
        throw std::runtime_error(fmt::format(
            "The graphs have different k-mers: {} source k-mers are missing from the"
            " target graph, {} target nodes were matched more than once, and {} of"
            " {} target nodes were matched",
            num_missing.load(), num_repeated.load(), num_mapped.load(), target.num_nodes()));
    }

    logger->trace("Matched {} of {} source rows to target rows",
                  num_mapped.load(), source.max_index());
}

sdsl::int_vector<> map_rows_to_graph(const DeBruijnGraph &source,
                                     const DeBruijnGraph &target,
                                     size_t num_threads) {
    sdsl::int_vector<> row_map(source.max_index(), 0, width_for(target.max_index()));

    // Entries of row_map are packed, so concurrent writes are guarded by a lock
    // per block. A block of 2^12 entries spans a whole number of 64-bit words.
    constexpr uint64_t kBlockBits = 12;
    constexpr uint64_t kNumLocks = 1 << 12;
    std::unique_ptr<std::mutex[]> locks(new std::mutex[kNumLocks]);

    match_rows(source, target, num_threads, [&](uint64_t row, node_index target_node) {
        std::lock_guard<std::mutex> lock(locks[(row >> kBlockBits) % kNumLocks]);
        row_map[row] = target_node;
    });

    return row_map;
}

void map_rows_to_graph(const DeBruijnGraph &source,
                       const DeBruijnGraph &target,
                       const std::string &filename,
                       const std::filesystem::path &tmp_dir,
                       size_t num_threads) {
    const uint64_t num_rows = source.max_index();
    const uint64_t width = width_for(target.max_index());

    // The walk visits source rows in no particular order. Each pair is appended to
    // the file of its bucket, a range of 2^bucket_bits source rows, so that every
    // bucket can then be placed in memory on its own. A bucket spans a whole
    // number of 64-bit words of the row map.
    constexpr uint64_t kMaxBuckets = 256;
    constexpr size_t kBufferSize = 1 << 16;
    uint64_t bucket_bits = 6;
    while (((num_rows + (1ull << bucket_bits) - 1) >> bucket_bits) > kMaxBuckets) {
        ++bucket_bits;
    }
    const uint64_t num_buckets = (num_rows + (1ull << bucket_bits) - 1) >> bucket_bits;
    const size_t record_size = (bucket_bits + width + 7) / 8;

    struct Bucket {
        std::mutex mutex;
        std::string filename;
        std::ofstream out;
        std::vector<char> buffer;
    };
    std::vector<Bucket> buckets(num_buckets);

    auto remove_buckets = [&]() {
        for (auto &bucket : buckets) {
            bucket.out.close();
            std::error_code ec;
            std::filesystem::remove(bucket.filename, ec);
        }
    };

    try {
        for (uint64_t b = 0; b < num_buckets; ++b) {
            buckets[b].filename = (tmp_dir / fmt::format("bucket.{}", b)).string();
            buckets[b].out.open(buckets[b].filename, std::ios::binary);
            buckets[b].buffer.reserve(kBufferSize + record_size);
        }

        match_rows(source, target, num_threads, [&](uint64_t row, node_index target_node) {
            unsigned __int128 record = row & ((1ull << bucket_bits) - 1);
            record = (record << width) | target_node;
            Bucket &bucket = buckets[row >> bucket_bits];
            std::lock_guard<std::mutex> lock(bucket.mutex);
            const char *bytes = reinterpret_cast<const char *>(&record);
            bucket.buffer.insert(bucket.buffer.end(), bytes, bytes + record_size);
            if (bucket.buffer.size() >= kBufferSize) {
                bucket.out.write(bucket.buffer.data(), bucket.buffer.size());
                bucket.buffer.clear();
            }
        });

        for (auto &bucket : buckets) {
            bucket.out.write(bucket.buffer.data(), bucket.buffer.size());
            bucket.buffer = {};
            bucket.out.close();
            if (!bucket.out)
                throw std::runtime_error("Cannot write to " + bucket.filename);
        }

        {
            std::ofstream out(filename, std::ios::binary);
            sdsl::int_vector<>::write_header(num_rows * width, width, out);
            if (!out)
                throw std::runtime_error("Cannot write to " + filename);
        }
        const uint64_t header_size = std::filesystem::file_size(filename);
        std::filesystem::resize_file(filename, header_size + ((num_rows * width + 63) / 64) * 8);

        // Each thread holds one bucket of the row map in memory.
        std::atomic<bool> failed(false);
        #pragma omp parallel for num_threads(std::min<size_t>(std::max<size_t>(1, num_threads), 16)) schedule(dynamic)
        for (uint64_t b = 0; b < num_buckets; ++b) {
            const uint64_t begin = b << bucket_bits;
            sdsl::int_vector<> rows(std::min<uint64_t>(num_rows - begin, 1ull << bucket_bits), 0, width);

            std::ifstream in(buckets[b].filename, std::ios::binary);
            std::vector<char> buffer(record_size * (kBufferSize / record_size));
            while (in.read(buffer.data(), buffer.size()) || in.gcount()) {
                for (std::streamsize i = 0; i + record_size <= static_cast<size_t>(in.gcount()); i += record_size) {
                    unsigned __int128 record = 0;
                    std::memcpy(&record, buffer.data() + i, record_size);
                    rows[static_cast<uint64_t>(record >> width)]
                        = static_cast<uint64_t>(record) & sdsl::bits::lo_set[width];
                }
            }

            std::fstream out(filename, std::ios::binary | std::ios::in | std::ios::out);
            out.seekp(header_size + (begin * width) / 8);
            out.write(reinterpret_cast<const char *>(rows.data()), ((rows.size() * width + 63) / 64) * 8);
            out.close();
            if (!out || in.bad())
                failed = true;
        }
        if (failed)
            throw std::runtime_error("Cannot write row map to " + filename);
    } catch (...) {
        remove_buckets();
        throw;
    }
    remove_buckets();
}

// A row map held in memory. Every thread reads the same vector.
struct RowMapInMemory {
    const sdsl::int_vector<> &row_map;

    uint64_t size() const { return row_map.size(); }
    const sdsl::int_vector<>& open() const { return row_map; }
};

// A row map in a file. Every thread opens its own buffered reader, so the map
// is never held in memory. Reading it through mmap is slow on network file systems.
struct RowMapInFile {
    explicit RowMapInFile(const std::string &filename)
          : filename(filename), num_rows(open().size()) {}

    uint64_t size() const { return num_rows; }
    sdsl::int_vector_buffer<> open() const { return sdsl::int_vector_buffer<>(filename); }

    std::string filename;
    uint64_t num_rows;
};

template <class RowMap>
static bool check_row_map(const DeBruijnGraph &source,
                          const DeBruijnGraph &target,
                          const RowMap &row_map,
                          uint64_t num_samples) {
    if (row_map.size() != source.max_index())
        return false;

    const uint64_t step = std::max<uint64_t>(1, row_map.size() / std::max<uint64_t>(1, num_samples));
    uint64_t num_kept = 0;
    uint64_t max_target_row = 0;
    bool matches = true;
    #pragma omp parallel num_threads(get_num_threads()) reduction(+:num_kept) reduction(max:max_target_row) reduction(&&:matches)
    {
        auto &&rows = row_map.open();
        #pragma omp for schedule(static, 1 << 20)
        for (uint64_t row = 0; row < row_map.size(); ++row) {
            uint64_t target_row = rows[row];
            if (!target_row)
                continue;

            ++num_kept;
            max_target_row = std::max(max_target_row, target_row);
            if (row % step == 0
                    && target.kmer_to_node(source.get_node_sequence(row + 1)) != target_row)
                matches = false;
        }
    }
    return matches && num_kept == target.num_nodes() && max_target_row <= target.max_index();
}

bool check_row_map(const DeBruijnGraph &source,
                   const DeBruijnGraph &target,
                   const sdsl::int_vector<> &row_map,
                   uint64_t num_samples) {
    return check_row_map(source, target, RowMapInMemory { row_map }, num_samples);
}

bool check_row_map(const DeBruijnGraph &source,
                   const DeBruijnGraph &target,
                   const std::string &row_map_filename,
                   uint64_t num_samples) {
    return check_row_map(source, target, RowMapInFile(row_map_filename), num_samples);
}

template <class RowMap>
static std::unique_ptr<MultiLabelAnnotation<std::string>>
remap_annotation_rows(const MultiLabelAnnotation<std::string> &annotation,
                      const DeBruijnGraph &source,
                      const RowMap &row_map,
                      uint64_t target_num_rows,
                      size_t num_threads) {
    using matrix::BRWT;
    using matrix::RowDiff;

    num_threads = std::max<size_t>(1, num_threads);

    if (annotation.get_matrix().num_rows() != row_map.size())
        throw std::invalid_argument("The row map does not match the annotation's rows");

    uint64_t max_target_row = 0;
    #pragma omp parallel num_threads(num_threads) reduction(max:max_target_row)
    {
        auto &&rows = row_map.open();
        #pragma omp for schedule(static, 1 << 20)
        for (uint64_t row = 0; row < row_map.size(); ++row) {
            max_target_row = std::max<uint64_t>(max_target_row, rows[row]);
        }
    }
    if (max_target_row > target_num_rows)
        throw std::invalid_argument("The row map points past the end of the target graph");

    const auto &label_encoder = annotation.get_label_encoder();

    if (const auto *brwt = dynamic_cast<const MultiBRWTAnnotator *>(&annotation)) {
        auto &&rows = row_map.open();
        auto matrix = brwt->get_matrix().remap_rows(rows, target_num_rows, num_threads);
        return std::make_unique<MultiBRWTAnnotator>(std::move(matrix), label_encoder);
    }

    const auto *row_diff = dynamic_cast<const RowDiffBRWTAnnotator *>(&annotation);
    if (!row_diff)
        throw std::invalid_argument("Row remapping supports only brwt and row_diff_brwt annotations");

    const RowDiff<BRWT> &matrix = row_diff->get_matrix();
    const auto &anchor = matrix.anchor();
    const auto &fork_succ = matrix.fork_succ();

    if (anchor.size() != row_map.size())
        throw std::invalid_argument("The annotation has no anchors for every row");

    if (!fork_succ.size() && !dynamic_cast<const graph::DBGSuccinct *>(&source))
        throw std::invalid_argument("The annotation has no successor bitmap, which is"
                                    " needed for graphs other than succinct");

    // Keep each row's successor from the source graph. Successors are marked by
    // node index, so the bitmap has one extra entry.
    sdsl::bit_vector new_anchor(target_num_rows, false);
    sdsl::bit_vector new_succ(target_num_rows + 1, false);
    // Mark the successors by source row first, so that the row map is read only
    // in row order.
    sdsl::bit_vector is_succ(row_map.size(), false);
    std::atomic<uint64_t> num_lost_successors(0);
    #pragma omp parallel num_threads(num_threads)
    {
        auto &&rows = row_map.open();
        #pragma omp for schedule(static, 1 << 20)
        for (uint64_t row = 0; row < row_map.size(); ++row) {
            uint64_t target_row = rows[row];
            if (!target_row)
                continue;

            if (anchor[row]) {
                set_bit(new_anchor.data(), target_row - 1, true, std::memory_order_relaxed);
                continue;
            }

            node_index succ = row_diff_successor(source, row + 1, fork_succ);
            if (!succ || succ > row_map.size()) {
                ++num_lost_successors;
                continue;
            }
            set_bit(is_succ.data(), succ - 1, true, std::memory_order_relaxed);
        }

        #pragma omp for schedule(static, 1 << 20)
        for (uint64_t row = 0; row < row_map.size(); ++row) {
            if (!is_succ[row])
                continue;

            if (uint64_t target_row = rows[row]) {
                set_bit(new_succ.data(), target_row, true, std::memory_order_relaxed);
            } else {
                ++num_lost_successors;
            }
        }
    }
    is_succ = sdsl::bit_vector();

    if (num_lost_successors) {
        throw std::runtime_error(fmt::format(
            "{} row-diff successors of kept rows are dropped",
            num_lost_successors.load()));
    }

    logger->trace("Remapped anchors and successors");

    auto &&rows = row_map.open();
    auto diffs = matrix.diffs().remap_rows(rows, target_num_rows, num_threads);
    logger->trace("Remapped row-diff matrix: {} of {} relations kept",
                  diffs->num_relations(), matrix.diffs().num_relations());

    auto result = std::make_unique<RowDiff<BRWT>>(nullptr, std::move(*diffs));
    result->set_anchor(matrix::IRowDiff::anchor_bv_type(std::move(new_anchor)));
    result->set_fork_succ(matrix::IRowDiff::fork_succ_bv_type(std::move(new_succ)));

    return std::make_unique<RowDiffBRWTAnnotator>(std::move(result), label_encoder);
}

std::unique_ptr<MultiLabelAnnotation<std::string>>
remap_annotation_rows(const MultiLabelAnnotation<std::string> &annotation,
                      const DeBruijnGraph &source,
                      const sdsl::int_vector<> &row_map,
                      uint64_t target_num_rows,
                      size_t num_threads) {
    return remap_annotation_rows(annotation, source, RowMapInMemory { row_map },
                                 target_num_rows, num_threads);
}

std::unique_ptr<MultiLabelAnnotation<std::string>>
remap_annotation_rows(const MultiLabelAnnotation<std::string> &annotation,
                      const DeBruijnGraph &source,
                      const std::string &row_map_filename,
                      uint64_t target_num_rows,
                      size_t num_threads) {
    return remap_annotation_rows(annotation, source, RowMapInFile(row_map_filename),
                                 target_num_rows, num_threads);
}

} // namespace annot
} // namespace mtg
