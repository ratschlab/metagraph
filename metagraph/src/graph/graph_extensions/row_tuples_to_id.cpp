#include "row_tuples_to_id.hpp"

#include <sstream>

#include <tsl/hopscotch_map.h>

#include "common/serialization.hpp"
#include "common/utils/file_utils.hpp"
#include "common/threads/threading.hpp"

namespace mtg::graph {


RowTuplesToId::RowTuplesToId(const std::vector<std::string> &fai_infiles, size_t k)
      : seq_id_labels_(fai_infiles.size()), seq_delims_(fai_infiles.size()), k_(k) {
    size_t num_labels = fai_infiles.size();

    #pragma omp parallel for schedule(dynamic)
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
            if (len < k_)
                continue;
            cur_coord += len - k_ + 1;
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

std::vector<std::tuple<std::string, size_t, std::vector<SmallVector<uint64_t>>>> RowTuplesToId
::rows_tuples_to_label_tuples(const std::vector<RowTuples> &rows_tuples) const {
    using KmerCoords = std::vector<SmallVector<uint64_t>>;
    tsl::hopscotch_map<std::uint64_t, std::tuple<std::string, size_t, KmerCoords>> conv_coords(seq_id_labels_.size());
    for (size_t i = 0; i < rows_tuples.size(); ++i) {
        const auto &row_tuples = rows_tuples[i];
        for (const auto &[col, tuple] : row_tuples) {
            assert(col < seq_id_labels_.size());
            assert(col < seq_delims_.size());
            const auto &labels = seq_id_labels_[col];
            const auto &delims = seq_delims_[col];
            for (uint64_t coord : tuple) {
                uint64_t seq_id = coord ? delims.rank1(coord - 1) : 0;
                uint64_t conv_coord = seq_id > 0 ? coord - delims.select1(seq_id) - 1 : coord;
                auto it = conv_coords.find(seq_id);
                if (it == conv_coords.end()) {
                    it = conv_coords.try_emplace(seq_id, std::make_tuple(labels[seq_id], size_t(0), KmerCoords())).first;
                    std::get<2>(it.value()).resize(rows_tuples.size());
                }
                auto &[label, count, cur_tuples] = it.value();
                cur_tuples[i].emplace_back(conv_coord);
            }
        }
    }

    std::vector<std::tuple<std::string, size_t, std::vector<SmallVector<uint64_t>>>> result;
    result.reserve(conv_coords.size());
    for (auto it = conv_coords.begin(); it != conv_coords.end(); ++it) {
        auto &[label, count, cur_tuples] = it.value();
        count = std::count_if(cur_tuples.begin(), cur_tuples.end(),
                              [](const auto &a) { return !a.empty(); });
        result.emplace_back(std::move(it.value()));
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

        k_ = load_number(*in);

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
    serialize_number(out, k_);
}

bool RowTuplesToId::is_compatible(const SequenceGraph &graph, bool verbose) const {
    const auto *dbg = dynamic_cast<const DeBruijnGraph*>(&graph);
    if (!dbg) {
        if (verbose)
            std::cerr << "Incompatible: graph is not a DeBruijnGraph" << std::endl;
        return false;
    }
    if (dbg->get_k() != k_) {
        if (verbose)
            std::cerr << "Incompatible: graph k=" << dbg->get_k()
                      << " != extension k=" << k_ << std::endl;
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
