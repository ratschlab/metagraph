#ifndef __WEIGHTED_GRAPH_HPP__
#define __WEIGHTED_GRAPH_HPP__

#include "sequence_graph.hpp"
#include "dbg_succinct.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_bitmap.hpp"


template <typename Weights = sdsl::int_vector<>>
class IWeighted {
  public:
    using weight = typename Weights::value_type;
    using node_index = typename DeBruijnGraph::node_index;

    virtual uint64_t num_weights() const = 0;
    virtual void set_weights(const Weights &weights) = 0;
    virtual weight get_weight(node_index i) const = 0;
};


template <class T, typename Weights = sdsl::int_vector<>>
class Weighted : public T, public IWeighted<Weights> {
    static_assert(std::is_base_of<DeBruijnGraph, T>::value);

  public:
    using W = IWeighted<Weights>;
    using typename W::weight;
    using typename W::node_index;

    using T::T;

    virtual uint64_t num_weights() const { return this->num_nodes(); };
    virtual void set_weights(const Weights &weights) { weights_ = std::move(weights); };
    virtual weight get_weight(node_index i) const { return weights_[i]; };

  protected:
    Weights weights_;
};


template <class T, typename Weights = sdsl::int_vector<>>
class WeightedMixin : public Weighted<T, Weights> {
  public:

    using Weighted<T, Weights>::Weighted;

    using W = Weighted<T, Weights>;
    using typename W::weight;
    using typename W::node_index;

    virtual bool load(const std::string &filename_base);
    virtual void serialize(const std::string &filename_base) const;

  private:
    static constexpr auto kWeightsExtension = ".weights";
};

typedef WeightedMixin<DeBruijnGraph> WeightedDBG;
typedef WeightedMixin<DBGSuccinct> WeightedDBGSuccinct;
typedef WeightedMixin<DBGHashOrdered> WeightedDBGHashOrdered;
typedef WeightedMixin<DBGHashString> WeightedDBGHashString;
typedef WeightedMixin<DBGBitmap> WeightedDBGBitmap;

#endif // __WEIGHTED_GRAPH_HPP__
