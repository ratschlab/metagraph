#ifndef __MASKED_DBG_ALIGNER_HPP__
#define __MASKED_DBG_ALIGNER_HPP__

#include "dbg_aligner.hpp"

#include <tsl/hopscotch_map.h>

#include "graph/annotated_dbg.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/annotated_graph_algorithm.hpp"


namespace sdsl {
namespace util {

template <class t_int_vec>
typename t_int_vec::size_type next_zero(const t_int_vec& v, uint64_t idx)
{
    uint64_t pos = idx>>6;
    uint64_t node = ~v.data()[pos];
    node >>= (idx&0x3F);
    if (node) {
        return idx+bits::lo(node);
    } else {
        ++pos;
        while ((pos<<6) < v.bit_size()) {
            node = ~v.data()[pos];
            if (node) {
                return (pos<<6)|bits::lo(node);
            }
            ++pos;
        }
        return v.bit_size();
    }

}

} // namespace util
} // namespace sdsl


template <class BaseSeeder = ExactSeeder<>>
class LabeledSeeder : public BaseSeeder {
  public:
    typedef typename BaseSeeder::Seed Seed;
    typedef typename Seed::node_index node_index;

    LabeledSeeder(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : BaseSeeder(anno_graph.get_graph(), config), anno_graph_(&anno_graph) {}

    LabeledSeeder(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : BaseSeeder(graph, config), anno_graph_(nullptr) {
        throw std::runtime_error("Not implemented");
    }

    void call_seeds(std::function<void(Seed&&)> callback) const {
        assert(anno_graph_);
        const auto &config = BaseSeeder::get_config();

        BaseSeeder::call_seeds([&](Seed&& seed) {
            const auto label_signatures = anno_graph_->get_top_label_signatures(
                seed.get_sequence(),
                anno_graph_->get_annotation().num_labels()
            );

            for (const auto &[label, signature] : label_signatures) {
                size_t i = sdsl::util::next_bit(signature, 0);
                while (i < signature.size()) {
                    size_t next = sdsl::util::next_zero(signature, i + 1);

                    const char *begin = seed.get_query().data() + i;
                    size_t size = seed.get_query().size() - seed.size() + (next - i);

                    Seed cur_seed(std::string_view(begin, size),
                                  std::vector<node_index>(seed.begin() + i,
                                                          seed.begin() + next),
                                  config.match_score(std::string_view(begin, size)),
                                  seed.get_clipping() + i,
                                  seed.get_orientation(),
                                  i <= seed.get_offset() ? seed.get_offset() - i : 0);

                    cur_seed.set_label(label);
                    callback(std::move(cur_seed));

                    i = sdsl::util::next_bit(signature, next + 1);
                }
            }
        });
    }

  private:
    const AnnotatedDBG *anno_graph_ = nullptr;
};

template <class BaseExtender = DefaultColumnExtender<>>
class MaskedColumnExtender : public BaseExtender {
  public:
    typedef BaseExtender base_extender;
    typedef typename BaseExtender::DBGAlignment DBGAlignment;
    typedef typename BaseExtender::node_index node_index;
    typedef typename BaseExtender::score_t score_t;

    MaskedColumnExtender(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : BaseExtender(anno_graph.get_graph(), config), anno_graph_(&anno_graph) {}

    MaskedColumnExtender(const DeBruijnGraph &graph, const DBGAlignerConfig &config)
          : BaseExtender(graph, config), anno_graph_(nullptr) {
        throw std::runtime_error("Not implemented");
    }

    virtual
    void operator()(const DBGAlignment &path,
                    std::string_view query,
                    std::function<void(DBGAlignment&&, node_index)> callback,
                    bool orientation,
                    score_t min_path_score = std::numeric_limits<score_t>::min()) override {
        assert(anno_graph_);
        assert(path.get_label().size());

        auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(
            anno_graph_->get_graph_ptr()
        );

        assert(graph);

        std::function<bool(const node_index&)> check_node = [&](const node_index &node) {
            return anno_graph_->has_label(node, path.get_label());
        };

        if (graph->is_canonical_mode()) {
            check_node = [&](const node_index &node) {
                node_index canon_node;
                graph->map_to_nodes(graph->get_node_sequence(node), [&](auto i) {
                    canon_node = i;
                });

                return anno_graph_->has_label(canon_node, path.get_label());
            };
        }

        masked_graph_ = std::make_unique<MaskedDeBruijnGraph>(graph, std::move(check_node));

        assert(std::all_of(path.begin(), path.end(),
                           [&](auto node) { return masked_graph_->in_subgraph(node); }));

        BaseExtender::reset();
        BaseExtender::set_graph(*masked_graph_);
        BaseExtender::operator()(path, query,
                                 [&](DBGAlignment&& alignment, node_index node) {
            alignment.set_label(path.get_label());
            callback(std::move(alignment), node);
        }, orientation, min_path_score);
    }

  private:
    const AnnotatedDBG *anno_graph_ = nullptr;
    std::unique_ptr<MaskedDeBruijnGraph> masked_graph_;
};


template <class Seeder = ExactSeeder<>,
          class Extender = DefaultColumnExtender<>,
          class AlignmentCompare = std::less<Alignment<>>,
          template <class BaseSeeder> class LabeledSeeder = LabeledSeeder,
          template <class BaseExtender> class MaskedExtender = MaskedColumnExtender>
class MaskedDBGAligner : public DBGAligner<LabeledSeeder<Seeder>,
                                           MaskedExtender<Extender>,
                                           AlignmentCompare> {
  public:
    typedef DBGAligner<LabeledSeeder<Seeder>,
                       MaskedExtender<Extender>,
                       AlignmentCompare> Aligner;
    typedef typename Aligner::node_index node_index;
    typedef typename Aligner::DBGAlignment DBGAlignment;

    MaskedDBGAligner(const AnnotatedDBG &anno_graph, const DBGAlignerConfig &config)
          : Aligner(anno_graph.get_graph(), config), anno_graph_(anno_graph) {}

  private:
    virtual LabeledSeeder<Seeder> build_seeder() const override {
        return LabeledSeeder<Seeder>(anno_graph_, Aligner::get_config());
    }

    virtual MaskedExtender<Extender> build_extender() const override {
        return MaskedExtender<Extender>(anno_graph_, Aligner::get_config());
    }

    virtual void align_aggregate(std::string_view query,
                                 const std::function<void(DBGAlignment&&)> &callback,
                                 bool orientation,
                                 score_t min_path_score,
                                 const LabeledSeeder<Seeder> &seeder) const override {
        tsl::hopscotch_map<std::string, typename Aligner::AlignmentQueue> path_queues;
        size_t num_alternative_paths = Aligner::get_config().num_alternative_paths;

        Aligner::align(query,
                       [&](DBGAlignment&& alignment) {
                           path_queues.emplace(alignment.get_label(), num_alternative_paths);
                           path_queues[alignment.get_label()].emplace(std::move(alignment));
                       },
                       orientation,
                       min_path_score,
                       seeder);

        for (auto it = path_queues.begin(); it != path_queues.end(); ++it) {
            while (it->second.size()) {
                auto path = it.value().pop_top();
                assert(path.is_valid(Aligner::get_graph(), &Aligner::get_config()));
                assert(std::all_of(path.begin(), path.end(), [&](auto node) {
                    if (!anno_graph_.get_graph().is_canonical_mode())
                        return anno_graph_.has_label(node, path.get_label());

                    const auto &graph = anno_graph_.get_graph();
                    node_index canon_node;
                    graph.map_to_nodes(graph.get_node_sequence(node), [&](auto i) {
                        canon_node = i;
                    });

                    return anno_graph_.has_label(canon_node, path.get_label());
                }));
                callback(std::move(path));
            }
        }
    }

    const AnnotatedDBG &anno_graph_;
};

#endif // __MASKED_DBG_ALIGNER_HPP__
