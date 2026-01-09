#include "annotation/coord_to_header.hpp"

#include "common/serialization.hpp"
#include "common/utils/file_utils.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/base/sequence_graph.hpp"

namespace mtg {
namespace annot {

using Tuple = CoordToHeader::Tuple;

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
            const auto &offsets = coord_offsets_.at(col);
            for (uint64_t &coord : coords) {
                if (coord >= offsets.size()) {
                    throw std::runtime_error("Querying coordinate " + std::to_string(coord) + " for"
                            + " column " + std::to_string(col) + " while CoordToHeader has only "
                            + std::to_string(offsets.size()) + " coordinates for that column");
                }
                size_t header = coord ? offsets.rank1(coord - 1) : 0;
                // Convert global coordinate to local (sequence-based) coordinate
                uint64_t local_coord = !header ? coord : coord - offsets.select1(header) - 1;
                if (local_coord > std::numeric_limits<uint64_t>::max() / num_headers(col)) {
                    throw std::runtime_error(fmt::format("Local coordinate {} is too large for the "
                            "given {} headers in column {} to encode both local_coord and header_id"
                            " in a single 64-bit integer", local_coord, num_headers(col), col));
                }
                coord = local_coord * num_headers(col) + header;
            }
        }
    }
}

bool CoordToHeader::load(const std::string &filename_base) {
    std::unique_ptr<std::ifstream> in
        = utils::open_ifstream(utils::make_suffix(filename_base, kExtension));
    if (!in)
        return false;

    try {
        uint64_t num_columns = load_number(*in);
        headers_.resize(num_columns);
        coord_offsets_.resize(num_columns);

        for (uint64_t i = 0; i < num_columns; ++i) {
            load_string_vector(*in, &headers_[i]);
            coord_offsets_[i].load(*in);
        }
        return true;

    } catch (...) {
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
