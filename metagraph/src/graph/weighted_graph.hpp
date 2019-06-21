#ifndef __WEIGHTED_GRAPH_HPP__
#define __WEIGHTED_GRAPH_HPP__

#include "sequence_graph.hpp"
#include "dbg_succinct.hpp"

template <typename Weights = sdsl::int_vector<>>
class Weighted {
  public:
    using weight = typename Weights::value_type;
    using node_index = typename DeBruijnGraph::node_index;

    virtual uint64_t num_weights() const = 0;
    virtual void set_weights(const Weights &weights) = 0;
    virtual weight get_weight(node_index i) const = 0;
};


template <class T, typename Weights = sdsl::int_vector<>>
class WeightedMixin : public T, public Weighted<Weights> {

    //TODO: better if WeightedMixin<DeBruijnGraph> were possible
    static_assert(std::is_base_of<DeBruijnGraph, T>::value);
    static_assert(!std::is_same<DeBruijnGraph, T>::value);

    using T::T;

  public:

    using W = Weighted<Weights>;
    using typename W::weight;
    using typename W::node_index;

    virtual uint64_t num_weights() const { return T::num_nodes(); };
    virtual void set_weights(const Weights &weights) { weights_ = std::move(weights); };
    virtual weight get_weight(node_index i) const { return weights_[i]; };

    virtual bool load(const std::string &filename_base);
    virtual void serialize(const std::string &filename_base) const;

  private:
    Weights weights_;

    static constexpr auto kWeightsExtension = ".weights";
};

typedef WeightedMixin<DBGSuccinct> WeightedDBGSuccinct;

#endif // __WEIGHTED_GRAPH_HPP__
