#include "annotation/coord_to_header.hpp"

#include <algorithm>

#include <tsl/hopscotch_map.h>
#include <tsl/hopscotch_set.h>

#include "common/logger.hpp"
#include "common/serialization.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/string_utils.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/base/sequence_graph.hpp"

namespace mtg {
namespace annot {

using Tuple = CoordToHeader::Tuple;
using mtg::common::logger;

namespace {

// a record index hashed and compared by the header it names (two types: the set derives from
// both)
struct HeaderOfRecordHash {
    const std::vector<std::string> *headers;
    size_t operator()(size_t i) const {
        return std::hash<std::string_view>()((*headers)[i]);
    }
};

struct HeaderOfRecordEqual {
    const std::vector<std::string> *headers;
    bool operator()(size_t a, size_t b) const {
        return (*headers)[a] == (*headers)[b];
    }
};

} // namespace

// The records in column order get global ids. The index keeps the ids sorted by the hash of
// their header, then by id: 8 bytes per record. A lookup binary-searches the hash and compares
// the headers of the records that share it.
struct CoordToHeader::HeaderIndex {
    std::vector<uint64_t> offsets;  // the global id of each column's first record, and the total
    std::vector<uint64_t> by_hash;  // global record ids, by (header hash, id)
};

namespace {

uint64_t header_hash(std::string_view header) {
    return std::hash<std::string_view>()(header);
}

} // namespace

CoordToHeader::CoordToHeader() {}
CoordToHeader::~CoordToHeader() {}

const CoordToHeader::HeaderIndex& CoordToHeader::header_index() const {
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    if (!header_index_) {
        auto index = std::make_unique<HeaderIndex>();
        index->offsets.assign(headers_.size() + 1, 0);
        for (Column col = 0; col < headers_.size(); ++col) {
            index->offsets[col + 1] = index->offsets[col] + headers_[col].size();
        }
        std::vector<std::pair<uint64_t, uint64_t>> hashed;
        hashed.reserve(index->offsets.back());
        for (Column col = 0; col < headers_.size(); ++col) {
            for (size_t s = 0; s < headers_[col].size(); ++s) {
                hashed.emplace_back(header_hash(headers_[col][s]), index->offsets[col] + s);
            }
        }
        std::sort(hashed.begin(), hashed.end());
        index->by_hash.reserve(hashed.size());
        for (const auto &[hash, id] : hashed) {
            index->by_hash.push_back(id);
        }
        header_index_ = std::move(index);
    }
    // the index is immutable once built, so lookups run outside the lock
    return *header_index_;
}

void CoordToHeader::index_headers() const {
    header_index();
}

// The records holding |header|, at most |limit| of them, in column order
std::vector<std::pair<CoordToHeader::Column, size_t>>
CoordToHeader::records_of(std::string_view header, size_t limit) const {
    const HeaderIndex &index = header_index();
    auto record = [&](uint64_t id) {
        const Column col = std::upper_bound(index.offsets.begin(), index.offsets.end(), id)
                            - index.offsets.begin() - 1;
        return std::make_pair(col, static_cast<size_t>(id - index.offsets[col]));
    };
    auto name = [&](uint64_t id) -> std::string_view {
        const auto [col, s] = record(id);
        return headers_[col][s];
    };
    const uint64_t hash = header_hash(header);
    auto it = std::partition_point(index.by_hash.begin(), index.by_hash.end(),
                                   [&](uint64_t id) { return header_hash(name(id)) < hash; });
    std::vector<std::pair<Column, size_t>> records;
    for (; it != index.by_hash.end() && records.size() < limit
            && header_hash(name(*it)) == hash; ++it) {
        if (name(*it) == header)
            records.push_back(record(*it));
    }
    return records;
}

std::optional<std::pair<CoordToHeader::Column, size_t>>
CoordToHeader::find_header(std::string_view header) const {
    // ids of equal hash are in id order, which is column order
    const auto records = records_of(header, 1);
    if (records.empty())
        return std::nullopt;
    return records.front();
}

bool CoordToHeader::is_repeated_header(std::string_view header) const {
    return records_of(header, 2).size() == 2;
}

void CoordToHeader::find_duplicate_headers(const std::vector<std::string> &headers,
                                           const std::function<bool(size_t, size_t)> &callback) {
    tsl::hopscotch_set<size_t, HeaderOfRecordHash, HeaderOfRecordEqual>
            first(std::max<size_t>(headers.size(), 1),
                  HeaderOfRecordHash { &headers }, HeaderOfRecordEqual { &headers });
    for (size_t i = 0; i < headers.size(); ++i) {
        auto [it, inserted] = first.insert(i);
        if (!inserted && !callback(*it, i))
            return;
    }
}

CoordToHeader::CoordToHeader(std::vector<std::vector<std::string>> &&headers,
                             std::vector<std::vector<uint64_t>> &&num_kmers)
      : headers_(std::move(headers)), coord_offsets_(num_kmers.size()) {
    assert(headers_.size() == num_kmers.size());
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t j = 0; j < num_kmers.size(); ++j) {
        assert(num_kmers[j].size() == num_sequences(j));
        if (num_kmers[j].empty())
            continue;
        auto &offsets = num_kmers[j];
        // Ensure no zero k-mer counts (each sequence must have at least one k-mer)
        assert(std::find(offsets.begin(), offsets.end(), 0) == offsets.end());
        std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
        coord_offsets_[j] = bit_vector_sd([&](const auto &callback) {
            for (uint64_t cur_coord : offsets) {
                callback(cur_coord - 1);
            }
        }, offsets.back(), offsets.size());
    }
    reset_sequence_starts();
}

void CoordToHeader::reset_sequence_starts() {
    sequence_starts_once_ = std::make_unique<std::once_flag[]>(num_columns());
    sequence_starts_.assign(num_columns(), {});
}

const std::vector<uint64_t>& CoordToHeader::sequence_starts(Column col) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    std::call_once(sequence_starts_once_[col], [&]() {
        const auto &offsets = coord_offsets_[col];
        std::vector<uint64_t> starts;
        starts.reserve(num_sequences(col) + 1);
        starts.push_back(0);
        // coord_offsets_ has a set bit at the last coordinate of each sequence
        for (size_t i = 1; i <= num_sequences(col); ++i) {
            starts.push_back(offsets.select1(i) + 1);
        }
        assert(starts.back() == offsets.size());
        sequence_starts_[col] = std::move(starts);
    });
    return sequence_starts_[col];
}

void CoordToHeader::throw_coord_out_of_range(Column col, uint64_t coord) const {
    throw std::out_of_range(fmt::format("Coordinate {} for column {} out of range "
            "(CoordToHeader has {} coordinates for that column)",
            coord, col, coord_offsets_[col].size()));
}

std::pair<size_t, uint64_t>
CoordToHeader::map_single_coord(Column col, uint64_t coord) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    const auto &offsets = coord_offsets_[col];
    if (coord >= offsets.size()) {
        throw std::out_of_range(fmt::format("Coordinate {} for column {} out of range "
                "(CoordToHeader has {} coordinates for that column)", coord, col, offsets.size()));
    }
    size_t header = coord ? offsets.rank1(coord - 1) : 0;
    uint64_t local_coord = !header ? coord : coord - offsets.select1(header) - 1;
    return { header, local_coord };
}

std::pair<uint64_t, uint64_t>
CoordToHeader::sequence_coords(Column col, size_t seq_id) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    const auto &offsets = coord_offsets_[col];
    if (seq_id >= num_sequences(col)) {
        throw std::out_of_range(fmt::format("Sequence id {} out of range for column {} "
                "({} sequences)", seq_id, col, num_sequences(col)));
    }
    // coord_offsets_ has a set bit at the last coordinate of each sequence
    return { seq_id ? offsets.select1(seq_id) + 1 : 0, offsets.select1(seq_id + 1) };
}

uint64_t CoordToHeader::num_kmers_in_sequence(Column col, size_t seq_id) const {
    const auto [first, last] = sequence_coords(col, seq_id);
    return last - first + 1;
}

bool CoordToHeader::load(const std::string &filename_base) {
    const std::string path = utils::make_suffix(filename_base, kExtension);
    std::unique_ptr<std::ifstream> in = utils::open_ifstream(path);
    if (!in->good()) {
        logger->error("Cannot open CoordToHeader file '{}': {}", path,
                      utils::file_read_failure_detail(path));
        return false;
    }
    try {
        // read into new vectors, so that a failed load leaves this object as it was
        uint64_t num_columns = load_number(*in);
        std::vector<std::vector<std::string>> headers(num_columns);
        std::vector<bit_vector_sd> coord_offsets(num_columns);

        for (uint64_t i = 0; i < num_columns; ++i) {
            load_string_vector(*in, &headers[i]);
            coord_offsets[i].load(*in);
        }
        // the index views the headers being replaced: drop it with them
        std::lock_guard<std::mutex> lock(header_index_mutex_);
        header_index_.reset();
        headers_.swap(headers);
        coord_offsets_.swap(coord_offsets);
        reset_sequence_starts();
        return true;

    } catch (const std::exception &e) {
        logger->error("Cannot load CoordToHeader from '{}': {} (caught: {})", path,
                      utils::file_read_failure_detail(path), e.what());
        return false;
    } catch (...) {
        logger->error("Cannot load CoordToHeader from '{}': {} (caught unknown exception)",
                      path, utils::file_read_failure_detail(path));
        return false;
    }
}

void CoordToHeader::serialize(const std::string &filename_base) const {
    auto fname = utils::make_suffix(filename_base, kExtension);
    std::ofstream out = utils::open_new_ofstream(fname);
    if (!out)
        throw std::ios_base::failure("Couldn't open file " + fname + " for writing");
    serialize_number(out, headers_.size());
    for (size_t i = 0; i < headers_.size(); ++i) {
        serialize_string_vector(out, headers_[i]);
        coord_offsets_[i].serialize(out);
    }
}

} // namespace annot
} // namespace mtg
