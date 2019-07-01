#include "weighted_graph.hpp"
#include "dbg_succinct.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_bitmap.hpp"

#include "utils.hpp"


template <class DBG, typename Weights>
bool WeightedDBG<DBG, Weights>::load(const std::string &filename) {

    if constexpr (std::is_same<DBG, DBGSuccinct>::value) {
        if (!DBG::load_without_mask(filename))
            return false;
    } else if (!DBG::load(filename)) {
        return false;
    }

    const auto weights_filename = utils::remove_suffix(filename, DBG::file_extension())
                                        + DBG::file_extension()
                                        + kWeightsExtension;
    try {
        std::ifstream instream(weights_filename, std::ios::binary);
        this->weights_.load(instream);
        return DBG::num_nodes() + 1 == this->weights_.size();
    } catch (...) {
        std::cerr << "ERROR: Cannot load graph weights from file "
                  << weights_filename << std::endl;
        return false;
    }
}


template <class DBG, typename Weights>
void WeightedDBG<DBG, Weights>::serialize(const std::string &filename) const {

    DBG::serialize(filename);

    std::ofstream outstream(utils::remove_suffix(filename, DBG::file_extension())
                                + DBG::file_extension()
                                + kWeightsExtension,
                            std::ios::binary);

    this->weights_.serialize(outstream);
}


template class WeightedDBG<DBGSuccinct>;
template class WeightedDBG<DBGHashOrdered>;
template class WeightedDBG<DBGHashString>;
template class WeightedDBG<DBGBitmap>;
