#include "annotation/coord_to_accession.hpp"

#include <sstream>

#include <tsl/hopscotch_map.h>

#include "common/serialization.hpp"
#include "common/utils/file_utils.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/base/sequence_graph.hpp"

namespace mtg {
namespace annot {

using Tuple = CoordToAccession::Tuple;

CoordToAccession::CoordToAccession(std::vector<std::vector<std::string>> &&headers,
                                   std::vector<std::vector<uint64_t>> &&num_kmers)
      : seq_id_labels_(std::move(headers)), seq_delims_(num_kmers.size()) {
    assert(seq_id_labels_.size() == num_kmers.size());
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t i = 0; i < num_kmers.size(); ++i) {
        if (num_kmers[i].empty())
            continue;
        auto &offsets = num_kmers[i];
        assert(std::find(offsets.begin(), offsets.end(), 0) == offsets.end());
        std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
        seq_delims_[i] = bit_vector_sd([&](const auto &callback) {
            for (uint64_t cur_coord : offsets) {
                callback(cur_coord - 1);
            }
        }, offsets.back(), offsets.size());
    }
}

std::vector<std::tuple<std::string, size_t, std::vector<Tuple>>>
CoordToAccession::rows_tuples_to_label_tuples(const std::vector<RowTuples> &rows_tuples, size_t min_count) const {
    // RowTuples = Vector<std::pair<Column, Tuple>>
    tsl::hopscotch_map<std::pair<Column, size_t>, std::vector<Tuple>> conv_coords;
    for (size_t i = 0; i < rows_tuples.size(); ++i) {
        const auto &row_tuples = rows_tuples[i];
        for (const auto &[col, tuple] : row_tuples) {
            const auto &delims = seq_delims_[col];
            for (uint64_t coord : tuple) {
                size_t seq_id = coord ? delims.rank1(coord - 1) : 0;
                uint64_t conv_coord = seq_id > 0 ? coord - delims.select1(seq_id) - 1 : coord;
                auto it = conv_coords.try_emplace(std::make_pair(col, seq_id), rows_tuples.size()).first;
                it.value()[i].emplace_back(conv_coord);
            }
        }
    }

    std::vector<std::tuple<std::string, size_t, std::vector<Tuple>>> result;
    result.reserve(conv_coords.size());
    for (auto &[key, tuples] : conv_coords) {
        auto [col, seq_id] = key;
        size_t count = std::count_if(tuples.begin(), tuples.end(),
                                     [](const auto &a) { return !a.empty(); });
        if (count >= min_count)
            result.emplace_back(seq_id_labels_[col][seq_id], count, std::move(tuples));
    }

    return result;
}

bool CoordToAccession::load(const std::string &filename_base) {
    try {
        std::unique_ptr<std::ifstream> in
            = utils::open_ifstream(utils::make_suffix(filename_base, kExtension));
        if (!in->good())
            return false;

        uint64_t num_columns = load_number(*in);
        seq_id_labels_.resize(num_columns);
        seq_delims_.resize(num_columns);

        for (uint64_t i = 0; i < num_columns; ++i) {
            load_string_vector(*in, &seq_id_labels_[i]);
            seq_delims_[i].load(*in);
        }
        return true;

    } catch (...) {
        return false;
    }
}

void CoordToAccession::serialize(const std::string &filename_base) const {
    auto fname = utils::make_suffix(filename_base, kExtension);
    std::ofstream out = utils::open_new_ofstream(fname);
    if (!out)
        throw std::ios_base::failure("Couldn't open file " + fname + " for writing");
    serialize_number(out, seq_id_labels_.size());
    for (size_t i = 0; i < seq_id_labels_.size(); ++i) {
        serialize_string_vector(out, seq_id_labels_[i]);
        seq_delims_[i].serialize(out);
    }
}

} // namespace annot
} // namespace mtg

