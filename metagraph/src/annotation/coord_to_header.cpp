#include "annotation/coord_to_header.hpp"

#include <tsl/hopscotch_map.h>

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

// header -> (column, seq_id) of its first occurrence in column order
struct CoordToHeader::HeaderIndex {
    tsl::hopscotch_map<std::string_view, std::pair<Column, size_t>> map;
};

CoordToHeader::CoordToHeader() {}
CoordToHeader::~CoordToHeader() {}

// The index is never shared or carried over: its keys view the headers of the object that
// built it, so a copy builds its own on first use
CoordToHeader::CoordToHeader(const CoordToHeader &other)
      : headers_(other.headers_), coord_offsets_(other.coord_offsets_) {}

size_t CoordToHeader::build_header_index() const {
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    build_header_index_locked();
    return header_index_->map.size();
}

void CoordToHeader::build_header_index_locked() const {
    if (!header_index_) {
        auto index = std::make_unique<HeaderIndex>();
        size_t total = 0;
        for (const auto &column : headers_) {
            total += column.size();
        }
        index->map.reserve(total);
        for (Column col = 0; col < headers_.size(); ++col) {
            for (size_t s = 0; s < headers_[col].size(); ++s) {
                index->map.try_emplace(std::string_view(headers_[col][s]), col, s);
            }
        }
        logger->trace("Built header index with {} sequences in {} columns", total,
                      headers_.size());
        header_index_ = std::move(index);
        header_index_builds_++;
    }
}

std::optional<std::pair<CoordToHeader::Column, size_t>>
CoordToHeader::find_header(std::string_view header) const {
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    build_header_index_locked();
    auto it = header_index_->map.find(header);
    if (it == header_index_->map.end())
        return std::nullopt;
    return it->second;
}

size_t CoordToHeader::num_header_index_builds() const {
    std::lock_guard<std::mutex> lock(header_index_mutex_);
    return header_index_builds_;
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
}

void CoordToHeader::map_to_local_coords(std::vector<RowTuples> *rows) const {
    assert(rows);
    for (auto &row : *rows) {
        for (auto &[col, coords] : row) {
            const size_t n = num_sequences(col);
            for (uint64_t &coord : coords) {
                auto [seq_id, local_coord] = map_single_coord(col, coord);
                assert(n);
                if (local_coord > std::numeric_limits<uint64_t>::max() / n) {
                    throw std::runtime_error(fmt::format("Local coordinate {} is too large to "
                            "pack with a seq_id into a single 64-bit integer "
                            "({} sequences in column {})", local_coord, n, col));
                }
                coord = local_coord * n + seq_id;
            }
        }
    }
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

CoordToHeader::SequenceRange
CoordToHeader::sequence_range(Column col, uint64_t coord) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    const auto &offsets = coord_offsets_[col];
    if (coord >= offsets.size()) {
        throw std::out_of_range(fmt::format("Coordinate {} for column {} out of range "
                "(CoordToHeader has {} coordinates for that column)", coord, col, offsets.size()));
    }
    // as map_single_coord: a set bit marks the last coordinate of each sequence
    const size_t header = coord ? offsets.rank1(coord - 1) : 0;
    const uint64_t first = header ? offsets.select1(header) + 1 : 0;
    return { header, first, offsets.select1(header + 1) };
}

uint64_t CoordToHeader::last_coord(Column col, size_t seq_id) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    if (seq_id >= num_sequences(col)) {
        throw std::out_of_range(fmt::format("Sequence id {} out of range for column {} "
                "({} sequences)", seq_id, col, num_sequences(col)));
    }
    return coord_offsets_[col].select1(seq_id + 1);
}

uint64_t CoordToHeader::num_kmers_in_sequence(Column col, size_t seq_id) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    const auto &offsets = coord_offsets_[col];
    if (seq_id >= num_sequences(col)) {
        throw std::out_of_range(fmt::format("Sequence id {} out of range for column {} "
                "({} sequences)", seq_id, col, num_sequences(col)));
    }
    // coord_offsets_ has a set bit at the partial-sum boundary of each
    // sequence, so the k-mer count for sequence s is
    //     select1(s+1) - (s == 0 ? -1 : select1(s)).
    uint64_t end = offsets.select1(seq_id + 1);
    uint64_t start = seq_id ? offsets.select1(seq_id) + 1 : 0;
    return end - start + 1;
}

bool CoordToHeader::load(const std::string &filename_base) {
    const std::string path = utils::make_suffix(filename_base, kExtension);
    std::unique_ptr<std::ifstream> in = utils::open_ifstream(path);
    if (!in->good()) {
        logger->error("Cannot open CoordToHeader file '{}': {}", path,
                      utils::file_read_failure_detail(path));
        return false;
    }
    {
        // the index views the headers about to be replaced
        std::lock_guard<std::mutex> lock(header_index_mutex_);
        header_index_.reset();
    }

    try {
        uint64_t num_columns = load_number(*in);
        headers_.resize(num_columns);
        coord_offsets_.resize(num_columns);

        for (uint64_t i = 0; i < num_columns; ++i) {
            load_string_vector(*in, &headers_[i]);
            coord_offsets_[i].load(*in);
        }
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
