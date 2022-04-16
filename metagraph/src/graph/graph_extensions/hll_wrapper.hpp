#ifndef __HLL_WRAPPER_HPP__
#define __HLL_WRAPPER_HPP__

#include "graph/representation/base/sequence_graph.hpp"
#include "annotation/binary_matrix/hll/hll_matrix.hpp"


namespace mtg {
namespace graph {

template <class HLL = annot::binmat::HLLMatrix<>>
class HLLWrapper : public SequenceGraph::GraphExtension {
  public:
    HLLWrapper() {}

    template <typename... Args>
    HLLWrapper(Args&&... args) : hll_(std::forward<Args>(args)...) {}

    bool load(const std::string &filename_base) {
        std::string fname = utils::make_suffix(filename_base, kHLLExtension);
        std::ifstream fin(fname, std::ios::binary);
        return hll_.load(fin);
    }

    void serialize(const std::string &filename_base) const {
        std::string fname = utils::make_suffix(filename_base, kHLLExtension);
        std::ofstream fout(fname, std::ios::binary);
        hll_.serialize(fout);
    }

    bool is_compatible(const SequenceGraph &, bool = true) const { return true; }

    const HLL& data() const { return hll_; }

  private:
    HLL hll_;
    static constexpr auto kHLLExtension = ".hll";
};

} // namespace graph
} // namespace mtg

#endif // __MER_DISTANCES_HPP__
