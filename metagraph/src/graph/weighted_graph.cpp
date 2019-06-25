#include "weighted_graph.hpp"
#include "dbg_succinct.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_bitmap.hpp"

#include "utils.hpp"


template <class T, typename Weights>
bool WeightedMixin<T, Weights>::load(const std::string &filename) {

    if constexpr (std::is_same<T, DBGSuccinct>::value) {
        if (!T::load_without_mask(filename))
            return false;
    } else if (!T::load(filename)) {
        return false;
    }

    try {
        std::ifstream instream(utils::remove_suffix(filename, T::file_extension())
                                    + T::file_extension()
                                    + kWeightsExtension,
                               std::ios::binary);
        this->weights_.load(instream);
        return T::num_nodes() + 1 == this->weights_.size();
    } catch (...) {
        std::cerr << "ERROR: Cannot load graph weights from file "
                  << filename + kWeightsExtension << std::endl;
        return false;
    }
}


template <class T, typename Weights>
void WeightedMixin<T, Weights>::serialize(const std::string &filename) const {

    T::serialize(filename);

    std::ofstream outstream(utils::remove_suffix(filename, T::file_extension())
                                + T::file_extension()
                                + kWeightsExtension,
                            std::ios::binary);

    this->weights_.serialize(outstream);
}


template class WeightedMixin<DBGSuccinct>;
template class WeightedMixin<DBGHashOrdered>;
template class WeightedMixin<DBGHashString>;
template class WeightedMixin<DBGBitmap>;
