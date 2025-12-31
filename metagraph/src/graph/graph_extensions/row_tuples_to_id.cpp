#include "row_tuples_to_id.hpp"

#include <sstream>

#include <tsl/hopscotch_map.h>

#include "common/serialization.hpp"
#include "common/utils/file_utils.hpp"
#include "common/threads/threading.hpp"

namespace mtg::graph {

using Tuple = RowTuplesToId::Tuple;

RowTuplesToId::RowTuplesToId(const std::vector<std::string> &fai_infiles, size_t k)
      : seq_id_labels_(fai_infiles.size()), seq_delims_(fai_infiles.size()) {
    size_t num_labels = fai_infiles.size();

    #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    for (size_t i = 0; i < num_labels; ++i) {
        std::vector<uint64_t> delims;
        std::ifstream fin(fai_infiles[i]);
        std::string line;
        std::string seq_id;
        uint64_t len;
        uint64_t cur_coord = 0;
        while (std::getline(fin, line)) {
            std::istringstream sin(line);
            sin >> seq_id >> len;
            if (len < k)
                continue;
            cur_coord += len - k + 1;
            seq_id_labels_[i].emplace_back(std::move(seq_id));
            delims.emplace_back(cur_coord);
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

std::vector<std::tuple<std::string, size_t, std::vector<Tuple>>> RowTuplesToId
::rows_tuples_to_label_tuples(const std::vector<RowTuples> &rows_tuples) const {
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

bool RowTuplesToId::load(const std::string &filename_base) {
    const auto rowtuples_filename
            = utils::make_suffix(filename_base, kRowTuplesExtension);
    try {
        std::unique_ptr<std::ifstream> in = utils::open_ifstream(rowtuples_filename);
        if (!in->good())
            return false;

        uint64_t num_labels = load_number(*in);

        seq_id_labels_.resize(num_labels);
        seq_delims_.resize(num_labels);

        for (uint64_t i = 0; i < num_labels; ++i) {
            load_string_vector(*in, &seq_id_labels_[i]);
            seq_delims_[i].load(*in);
        }

        return true;

    } catch (...) {
        std::cerr << "ERROR: Cannot load sequence delimiters from file "
                  << rowtuples_filename << std::endl;
        return false;
    }
}

void RowTuplesToId::serialize(const std::string &filename_base) const {
    const auto fname = utils::make_suffix(filename_base, kRowTuplesExtension);

    std::ofstream out = utils::open_new_ofstream(fname);
    serialize_number(out, seq_id_labels_.size());
    for (size_t i = 0; i < seq_id_labels_.size(); ++i) {
        serialize_string_vector(out, seq_id_labels_[i]);
        seq_delims_[i].serialize(out);
    }
}

bool RowTuplesToId::is_compatible(const SequenceGraph &graph, bool verbose) const {
    const auto *dbg = dynamic_cast<const DeBruijnGraph*>(&graph);
    if (!dbg) {
        if (verbose)
            std::cerr << "Incompatible: graph is not a DeBruijnGraph" << std::endl;
        return false;
    }
    if (seq_id_labels_.size() != seq_delims_.size()) {
        if (verbose)
            std::cerr << "Incompatible: mapping corrupted -- different number of label sets and delimiter sets"
                      << std::endl;
        return false;
    }
    return true;
}

} // namespace mtg::graph
