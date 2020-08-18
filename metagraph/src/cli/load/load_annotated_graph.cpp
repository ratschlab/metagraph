#include "load_annotated_graph.hpp"

#include "graph/annotated_graph_algorithm.hpp"
#include "common/algorithms.hpp"
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
                       std::vector<std::string> &label_mask_in_post,
                       std::vector<std::string> &label_mask_out_post,
                       const DifferentialAssemblyConfig &diff_config,
                       size_t num_threads,
                       const sdsl::int_vector<> *init_counts = nullptr) {
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());

    if (!graph.get())
        throw std::runtime_error("Masking only supported for DeBruijnGraph");

    std::vector<std::vector<std::string>*> label_sets {
        &label_mask_in, &label_mask_out,
        &label_mask_in_post, &label_mask_out_post
    };

    // Remove non-present labels
    for (auto *label_set : label_sets) {
        assert(label_set);
        label_set->erase(
            std::remove_if(label_set->begin(), label_set->end(),
                           [&](const std::string &label) {
                               bool exists = anno_graph.label_exists(label);
                               if (!exists)
                                   logger->trace("Removing label {}", label);

                               return !exists;
                           }),
            label_set->end()
        );

        std::sort(label_set->begin(), label_set->end());
    }

    for (auto *label_set : label_sets) {
        for (auto *other_label_set : label_sets) {
            if (label_set == other_label_set)
                continue;

            if (utils::count_intersection(label_set->begin(), label_set->end(),
                                          other_label_set->begin(), other_label_set->end()))
                logger->warn("Overlapping label sets");
        }
    }

    logger->trace("Masked in: {}", fmt::join(label_mask_in, " "));
    logger->trace("Masked in (post-processing): {}", fmt::join(label_mask_in_post, " "));
    logger->trace("Masked out: {}", fmt::join(label_mask_out, " "));
    logger->trace("Masked out (post-processing): {}", fmt::join(label_mask_out_post, " "));

    return std::make_unique<MaskedDeBruijnGraph>(mask_nodes_by_label(
        anno_graph,
        label_mask_in, label_mask_out,
        label_mask_in_post, label_mask_out_post,
        diff_config, num_threads, init_counts
    ));
}

DifferentialAssemblyConfig parse_diff_config(const std::string &config_str,
                                             bool canonical) {
    DifferentialAssemblyConfig diff_config;
    diff_config.add_complement = canonical;

    auto vals = utils::split_string(config_str, ",");
    auto it = vals.begin();
    if (it != vals.end()) {
        diff_config.label_mask_in_kmer_fraction = std::stof(*it);
        ++it;
    }
    if (it != vals.end()) {
        diff_config.label_mask_in_unitig_fraction = std::stof(*it);
        ++it;
    }
    if (it != vals.end()) {
        diff_config.label_mask_out_kmer_fraction = std::stof(*it);
        ++it;
    }
    if (it != vals.end()) {
        diff_config.label_mask_out_unitig_fraction = std::stof(*it);
        ++it;
    }
    if (it != vals.end()) {
        diff_config.label_mask_other_unitig_fraction = std::stof(*it);
        ++it;
    }

    assert(it == vals.end());

    logger->trace("Per-kmer mask in fraction: {}", diff_config.label_mask_in_kmer_fraction);
    logger->trace("Per-unitig mask in fraction: {}", diff_config.label_mask_in_unitig_fraction);
    logger->trace("Per-kmer mask out fraction: {}", diff_config.label_mask_out_kmer_fraction);
    logger->trace("Per-unitig mask out fraction: {}", diff_config.label_mask_out_unitig_fraction);
    logger->trace("Per-unitig other label fraction: {}", diff_config.label_mask_other_unitig_fraction);
    logger->trace("Include reverse complements: {}", diff_config.add_complement);

    return diff_config;
}

void call_masked_graphs(const AnnotatedDBG &anno_graph, Config *config,
                        const CallMaskedGraphHeader &callback,
                        size_t num_parallel_graphs_masked,
                        size_t num_threads_per_graph) {
    assert(!config->label_mask_file.empty());

    std::ifstream fin(config->label_mask_file);
    if (!fin.good()) {
        throw std::iostream::failure("Failed to read label mask file");
        exit(1);
    }

    ThreadPool thread_pool(num_parallel_graphs_masked);
    std::vector<std::string> shared_foreground_labels;
    std::vector<std::string> shared_background_labels;

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        std::unique_ptr<sdsl::int_vector<>> shared_counts;
        if (line[0] == '@') {
            logger->trace("Counting shared k-mers");

            // shared in and out labels
            auto line_split = utils::split_string(line, "\t", false);
            if (line_split.size() <= 1 || line_split.size() > 3)
                throw std::iostream::failure("Each line in mask file must have 2-3 fields.");

            // sync all assembly jobs before clearing current shared_counts
            thread_pool.join();

            shared_foreground_labels = utils::split_string(line_split[1], ",");
            shared_background_labels = utils::split_string(
                line_split.size() == 3 ? line_split[2] : "",
                ","
            );

            continue;
        }

        thread_pool.enqueue([&](std::string line) {
            auto line_split = utils::split_string(line, "\t", false);
            if (line_split.size() <= 2 || line_split.size() > 4)
                throw std::iostream::failure("Each line in mask file must have 3-4 fields.");

            auto diff_config = parse_diff_config(line_split[1], config->canonical);

            if (config->enumerate_out_sequences)
                line_split[0] += ".";

            auto foreground_labels = utils::split_string(line_split[2], ",");
            auto background_labels = utils::split_string(
                line_split.size() == 4 ? line_split[3] : "",
                ","
            );

            callback(*mask_graph_from_labels(anno_graph,
                                             foreground_labels, background_labels,
                                             shared_foreground_labels,
                                             shared_background_labels,
                                             diff_config, num_threads_per_graph),
                     line_split[0]);
        }, line);
    }

    thread_pool.join();
}

} // namespace cli
} // namespace mtg
