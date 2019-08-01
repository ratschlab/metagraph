#ifndef __WEIGHTED_GRAPH_HPP__
#define __WEIGHTED_GRAPH_HPP__

#include "sequence_graph.hpp"
#include "dbg_succinct.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_hash_string.hpp"
#include "dbg_bitmap.hpp"


template <typename index, typename Weights = sdsl::int_vector<>>
class IWeighted {
  public:
    using weight = typename Weights::value_type;

    virtual void set_weights(Weights&& weights) = 0;
    virtual weight get_weight(index i) const = 0;
};


template <class DBG, typename Weights = sdsl::int_vector<>>
class IWeightedDBG : public DBG,
                     public IWeighted<typename DBG::node_index, Weights> {
    static_assert(std::is_base_of<DeBruijnGraph, DBG>::value);

  public:
    using typename DBG::node_index;
    using IW = IWeighted<node_index, Weights>;
    using typename IW::weight;

    using DBG::DBG;

    virtual void set_weights(Weights&& weights) { weights_ = std::move(weights); };
    virtual weight get_weight(node_index i) const { return weights_[i]; };

    virtual bool load(const std::string &filename_base) = 0;
    virtual void serialize(const std::string &filename_base) const = 0;

  protected:
    Weights weights_;
};


template <class DBG, typename Weights = sdsl::int_vector<>>
class WeightedDBG : public IWeightedDBG<DBG, Weights> {
  public:
    using IW = IWeightedDBG<DBG, Weights>;
    using typename IW::weight;
    using typename IW::node_index;

    using IW::IWeightedDBG;

    virtual bool load(const std::string &filename_base);
    virtual void serialize(const std::string &filename_base) const;

  private:
    static constexpr auto kWeightsExtension = ".weights";
};

typedef WeightedDBG<DBGSuccinct> WeightedDBGSuccinct;
typedef WeightedDBG<DBGHashOrdered> WeightedDBGHashOrdered;
typedef WeightedDBG<DBGHashString> WeightedDBGHashString;
typedef WeightedDBG<DBGBitmap> WeightedDBGBitmap;

#endif // __WEIGHTED_GRAPH_HPP__
