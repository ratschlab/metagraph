#include "annotation/coord_to_header.hpp"

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

CoordToHeader::CoordToHeader(std::vector<std::vector<std::string>> &&headers,
                             std::vector<std::vector<uint64_t>> &&num_kmers)
      : headers_(std::move(headers)), coord_offsets_(num_kmers.size()) {
    assert(headers_.size() == num_kmers.size());
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t j = 0; j < num_kmers.size(); ++j) {
        assert(num_kmers[j].size() == num_headers(j));
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
            // n > 0 is guaranteed: map_single_coord throws below if col is
            // out of range or coords refers to an empty column.
            const size_t n = num_headers(col);
            for (uint64_t &coord : coords) {
                auto [seq_id, local_coord] = map_single_coord(col, coord);
                if (local_coord > std::numeric_limits<uint64_t>::max() / n) {
                    throw std::runtime_error(fmt::format("Local coordinate {} is too large to "
                            "pack with a seq_id into a single 64-bit integer "
                            "({} sequences in column {})",
                            local_coord, n, col));
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
                "(CoordToHeader has {} coordinates for that column)",
                coord, col, offsets.size()));
    }
    size_t header = coord ? offsets.rank1(coord - 1) : 0;
    uint64_t local_coord = !header ? coord : coord - offsets.select1(header) - 1;
    return { header, local_coord };
}

uint64_t CoordToHeader::num_kmers_in_sequence(Column col, size_t seq_id) const {
    if (col >= num_columns()) {
        throw std::out_of_range(fmt::format("Column {} out of range "
                "(CoordToHeader has {} columns)", col, num_columns()));
    }
    const auto &offsets = coord_offsets_[col];
    if (seq_id >= num_headers(col)) {
        throw std::out_of_range(fmt::format("Sequence id {} out of range for column {} "
                "({} sequences)", seq_id, col, num_headers(col)));
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
