#include "coord_to_accession.hpp"

#include <sstream>

#include <tsl/hopscotch_map.h>

#include "common/serialization.hpp"
#include "common/utils/file_utils.hpp"
#include "common/threads/threading.hpp"

namespace mtg::graph {

using Tuple = CoordToAccession::Tuple;

CoordToAccession::CoordToAccession(const std::vector<std::vector<std::pair<std::string, uint64_t>>> &accessions,
                                   const std::vector<std::string> &col_names)
      : seq_id_labels_(col_names.size()), seq_delims_(col_names.size()) {
    assert(accessions.size() == col_names.size());
    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t i = 0; i < col_names.size(); ++i) {
        std::vector<uint64_t> delims;
        uint64_t cur_coord = 0;
        for (const auto &[seq_id, num_kmers] : accessions[i]) {
            assert(num_kmers);
            cur_coord += num_kmers;
            delims.emplace_back(cur_coord);
            seq_id_labels_[i].emplace_back(seq_id);
        }
        if (delims.size()) {
            seq_delims_[i] = bit_vector_sd([&](const auto &callback) {
                for (uint64_t cur_coord : delims) {
                    callback(cur_coord - 1);
                }
            }, delims.back(), delims.size());
        }
    }
}

CoordToAccession::CoordToAccession(const std::vector<std::string> &fnames) {
    for (size_t i = 0; i < fnames.size(); ++i) {
        CoordToAccession rt;
        if (!rt.load(fnames[i]))
            throw std::runtime_error("Cannot load CoordToAccession mapping from file "
                                     + utils::make_suffix(fnames[i], kExtension));
        seq_id_labels_.insert(seq_id_labels_.end(),
                              std::make_move_iterator(rt.seq_id_labels_.begin()),
                              std::make_move_iterator(rt.seq_id_labels_.end()));
        seq_delims_.insert(seq_delims_.end(),
                           std::make_move_iterator(rt.seq_delims_.begin()),
                           std::make_move_iterator(rt.seq_delims_.end()));
    }
    // remove the files
    for (size_t i = 0; i < fnames.size(); ++i) {
        std::filesystem::remove(utils::make_suffix(fnames[i], kExtension));
    }
}

std::vector<std::tuple<std::string, size_t, std::vector<Tuple>>>
CoordToAccession::rows_tuples_to_label_tuples(const std::vector<RowTuples> &rows_tuples) const {
    // RowTuples = Vector<std::pair<Column, Tuple>>
    tsl::hopscotch_map<std::pair<Column, size_t>, std::vector<Tuple>> conv_coords;
    for (size_t i = 0; i < rows_tuples.size(); ++i) {
        const auto &row_tuples = rows_tuples[i];
        for (const auto &[col, tuple] : row_tuples) {
            const auto &delims = seq_delims_[col];
            for (uint64_t coord : tuple) {
                uint64_t seq_id = coord ? delims.rank1(coord - 1) : 0;
                uint64_t conv_coord = seq_id > 0 ? coord - delims.select1(seq_id) - 1 : coord;
                auto &tuple = conv_coords[std::make_pair(col, seq_id)];
                tuple.resize(rows_tuples.size());
                tuple[i].emplace_back(conv_coord);
            }
        }
    }

    std::vector<std::tuple<std::string, size_t, std::vector<Tuple>>> result;
    result.reserve(conv_coords.size());
    for (auto &[key, tuples] : conv_coords) {
        auto [col, seq_id] = key;
        size_t count = std::count_if(tuples.begin(), tuples.end(),
                                     [](const auto &a) { return !a.empty(); });
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
    std::ofstream out = utils::open_new_ofstream(utils::make_suffix(filename_base, kExtension));
    serialize_number(out, seq_id_labels_.size());
    for (size_t i = 0; i < seq_id_labels_.size(); ++i) {
        serialize_string_vector(out, seq_id_labels_[i]);
        seq_delims_[i].serialize(out);
    }
}

bool CoordToAccession::is_compatible(const SequenceGraph &graph, bool verbose) const {
    assert(seq_id_labels_.size() == seq_delims_.size());
    const auto *dbg = dynamic_cast<const DeBruijnGraph*>(&graph);
    if (!dbg) {
        if (verbose)
            std::cerr << "Incompatible: graph is not a DeBruijnGraph" << std::endl;
        return false;
    }
    return true;
}

} // namespace mtg::graph
