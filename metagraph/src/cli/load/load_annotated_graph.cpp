#include "load_annotated_graph.hpp"

#include "graph/annotated_graph_algorithm.hpp"
#include "common/logger.hpp"
#include "common/utils/string_utils.hpp"
#include "common/threads/threading.hpp"
#include "load_graph.hpp"
#include "load_annotation.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;

using mtg::common::logger;


std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(std::shared_ptr<DeBruijnGraph> graph,
                                                       const Config &config) {
    auto annotation_temp = config.infbase_annotators.size()
            ? initialize_annotation(config.infbase_annotators.at(0), config, 0)
            : initialize_annotation(config.anno_type, config, graph->max_index());

    if (config.infbase_annotators.size()
            && !annotation_temp->load(config.infbase_annotators.at(0))) {
        logger->error("Cannot load annotations for graph {}, file corrupted",
                      config.infbase);
        exit(1);
    }

    // load graph
    auto anno_graph = std::make_unique<AnnotatedDBG>(std::move(graph),
                                                     std::move(annotation_temp));

    if (!anno_graph->check_compatibility()) {
        logger->error("Graph and annotation are not compatible");
        exit(1);
    }

    return anno_graph;
}

std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(const Config &config) {
    return initialize_annotated_dbg(load_critical_dbg(config.infbase), config);
}


std::unique_ptr<MaskedDeBruijnGraph>
mask_graph_from_labels(const AnnotatedDBG &anno_graph,
                       std::vector<std::string> &label_mask_in,
                       std::vector<std::string> &label_mask_out,
                       Config *config,
                       size_t num_threads) {
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());

    if (!graph.get())
        throw std::runtime_error("Masking only supported for DeBruijnGraph");

    // Remove non-present labels
    label_mask_in.erase(
        std::remove_if(label_mask_in.begin(), label_mask_in.end(),
                       [&](const auto &label) {
                           bool exists = anno_graph.label_exists(label);
                           if (!exists)
                               logger->trace("Removing mask-in label {}", label);

                           return !exists;
                       }),
        label_mask_in.end()
    );

    label_mask_out.erase(
        std::remove_if(label_mask_out.begin(), label_mask_out.end(),
                       [&](const auto &label) {
                           bool exists = anno_graph.label_exists(label);
                           if (!exists)
                               logger->trace("Removing mask-out label {}", label);

                           return !exists;
                       }),
        label_mask_out.end()
    );

    logger->trace("Masked in: {}", fmt::join(label_mask_in, " "));
    logger->trace("Masked out: {}", fmt::join(label_mask_out, " "));

    if (!config->filter_by_kmer) {
        return std::make_unique<MaskedDeBruijnGraph>(
            make_masked_graph_by_unitig_labels(
                anno_graph,
                label_mask_in,
                label_mask_out,
                num_threads,
                config->label_mask_in_fraction,
                config->label_mask_out_fraction,
                config->label_other_fraction,
                config->canonical
            )
        );
    }

    size_t num_mask_in = label_mask_in.size();
    size_t num_mask_out = label_mask_out.size();

    return std::make_unique<MaskedDeBruijnGraph>(
        graph,
        mask_nodes_by_node_label(
            anno_graph,
            label_mask_in,
            label_mask_out,
            [config,num_mask_in,num_mask_out,&anno_graph](auto index,
                                                          auto get_num_in_labels,
                                                          auto get_num_out_labels) {
                assert(index != DeBruijnGraph::npos);

                size_t num_in_labels = get_num_in_labels();

                if (num_in_labels < config->label_mask_in_fraction * num_mask_in)
                    return false;

                size_t num_out_labels = get_num_out_labels();

                if (num_out_labels < config->label_mask_out_fraction * num_mask_out)
                    return false;

                size_t num_total_labels = anno_graph.get_labels(index).size();

                return num_total_labels - num_in_labels - num_out_labels
                            <= config->label_other_fraction * num_total_labels;
            },
            num_threads,
            config->identity == Config::ASSEMBLE ? 1.0 : 0.05
        )
    );
}

void
call_masked_graphs(const AnnotatedDBG &anno_graph, Config *config,
                   const std::function<void(const MaskedDeBruijnGraph&,
                                            const std::string& /* header */)> &callback,
                   size_t num_parallel_graphs_masked,
                   size_t num_threads_per_graph) {
    if (config->label_mask_file.empty()) {
        callback(*mask_graph_from_labels(anno_graph,
                                         config->label_mask_in,
                                         config->label_mask_out,
                                         config,
                                         num_threads_per_graph),
                 config->header);
        return;
    }

    std::ifstream fin(config->label_mask_file);
    if (!fin.good()) {
        throw std::iostream::failure("Failed to read label mask file");
        exit(1);
    }

    ThreadPool thread_pool(num_parallel_graphs_masked);

    std::string line;
    while (std::getline(fin, line)) {
        thread_pool.enqueue([&](std::string line) {
            auto line_split = utils::split_string(line, "\t");
            if (line_split.size() <= 1)
                return;

            if (line_split.size() > 3)
                throw std::iostream::failure("Each line in mask file must have 2-3 fields.");

            if (config->enumerate_out_sequences)
                line_split[0] += ".";

            auto foreground_labels = utils::split_string(line_split[1], ",");
            auto background_labels = utils::split_string(
                line_split.size() == 3 ? line_split[2] : "",
                ","
            );

            callback(*mask_graph_from_labels(anno_graph,
                                             foreground_labels,
                                             background_labels,
                                             config,
                                             num_threads_per_graph),
                     line_split[0]);
        }, line);
    }
}

std::unique_ptr<MaskedDeBruijnGraph>
mask_graph(const AnnotatedDBG &anno_graph, Config *config) {
    assert(config->label_mask_file.empty());
    return mask_graph_from_labels(anno_graph,
                                  config->label_mask_in,
                                  config->label_mask_out,
                                  config,
                                  get_num_threads());
}

} // namespace cli
} // namespace mtg
