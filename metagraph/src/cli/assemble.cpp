#include "assemble.hpp"

#include <fmt/format.h>

#include "common/algorithms.hpp"
#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/masked_graph.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using mtg::graph::AnnotatedDBG;
using mtg::graph::DeBruijnGraph;
using mtg::graph::MaskedDeBruijnGraph;
using mtg::graph::DifferentialAssemblyConfig;


void clean_label_set(const AnnotatedDBG &anno_graph,
                     std::vector<std::string> &label_set) {
    label_set.erase(std::remove_if(label_set.begin(), label_set.end(),
        [&](const std::string &label) {
            bool exists = anno_graph.label_exists(label);
            if (!exists)
                logger->trace("Removing label {}", label);

            return !exists;
        }
    ), label_set.end());

    std::sort(label_set.begin(), label_set.end());
    auto end = std::unique(label_set.begin(), label_set.end());
    for (auto it = end; it != label_set.end(); ++it) {
        logger->trace("Removing duplicate label {}", *it);
    }
    label_set.erase(end, label_set.end());
}


std::unique_ptr<MaskedDeBruijnGraph>
mask_graph_from_labels(const AnnotatedDBG &anno_graph,
                       const std::vector<std::string> &label_mask_in,
                       const std::vector<std::string> &label_mask_out,
                       const std::vector<std::string> &label_mask_in_post,
                       const std::vector<std::string> &label_mask_out_post,
                       const DifferentialAssemblyConfig &diff_config,
                       size_t num_threads,
                       const sdsl::int_vector<> *init_counts = nullptr) {
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());

    if (!graph.get())
        throw std::runtime_error("Masking only supported for DeBruijnGraph");

    std::vector<const std::vector<std::string>*> label_sets {
        &label_mask_in, &label_mask_out,
        &label_mask_in_post, &label_mask_out_post
    };

    for (const auto *label_set : label_sets) {
        for (const auto *other_label_set : label_sets) {
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

typedef std::function<void(const graph::MaskedDeBruijnGraph&,
                           const std::string& /* header */)> CallMaskedGraphHeader;

void call_masked_graphs(const AnnotatedDBG &anno_graph, Config *config,
                        const CallMaskedGraphHeader &callback,
                        size_t num_parallel_graphs_masked = 1,
                        size_t num_threads_per_graph = 1) {
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

            clean_label_set(anno_graph, shared_foreground_labels);
            clean_label_set(anno_graph, shared_background_labels);

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

            clean_label_set(anno_graph, foreground_labels);
            clean_label_set(anno_graph, background_labels);

            callback(*mask_graph_from_labels(anno_graph,
                                             foreground_labels, background_labels,
                                             shared_foreground_labels,
                                             shared_background_labels,
                                             diff_config, num_threads_per_graph),
                     line_split[0]);
        }, std::move(line));
    }

    thread_pool.join();
}


int assemble(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(files.size() == 1);
    assert(config->outfbase.size());

    Timer timer;
    logger->trace("Graph loading...");

    auto graph = load_critical_dbg(files.at(0));

    logger->trace("Graph loaded in {} sec", timer.elapsed());

    if (config->infbase_annotators.size()) {
        assert(config->label_mask_file.size());
        auto anno_graph = initialize_annotated_dbg(graph, *config);

        logger->trace("Generating masked graphs...");

        std::filesystem::remove(
            utils::remove_suffix(config->outfbase, ".gz", ".fasta") + ".fasta.gz"
        );

        std::mutex file_open_mutex;
        std::mutex write_mutex;

        size_t num_threads_per_traversal = get_num_threads() / config->parallel_assemblies;

        call_masked_graphs(*anno_graph, config,
            [&](const graph::MaskedDeBruijnGraph &graph, const std::string &header) {
                std::lock_guard<std::mutex> file_lock(file_open_mutex);
                seq_io::FastaWriter writer(config->outfbase, header,
                                           config->enumerate_out_sequences,
                                           get_num_threads() > 1, /* async write */
                                           true /* append */);

                if (config->unitigs || config->min_tip_size > 1) {
                    graph.call_unitigs([&](const auto &unitig, auto&&) {
                                           std::lock_guard<std::mutex> lock(write_mutex);
                                           writer.write(unitig);
                                       },
                                       num_threads_per_traversal,
                                       config->min_tip_size,
                                       config->kmers_in_single_form);
                } else {
                    graph.call_sequences([&](const auto &seq, auto&&) {
                                             std::lock_guard<std::mutex> lock(write_mutex);
                                             writer.write(seq);
                                         },
                                         num_threads_per_traversal,
                                         config->kmers_in_single_form);
                }
            },
            config->parallel_assemblies,
            num_threads_per_traversal
        );

        return 0;
    }

    logger->trace("Extracting sequences from graph...");

    timer.reset();

    if (config->to_gfa) {
        if (!config->unitigs) {
            logger->error("Flag '--unitigs' must be set for GFA output");
            exit(1);
        }

        logger->trace("Writing graph to GFA...");

        std::ofstream gfa_file(utils::remove_suffix(config->outfbase, ".gfa") + ".gfa");
        std::mutex str_mutex;

        gfa_file << "H\tVN:Z:1.0" << std::endl;
        size_t k = graph->get_k();
        size_t overlap = k - 1;
        graph->call_unitigs(
            [&](const auto &unitig, const auto &path) {
                std::ostringstream ostr;
                if (config->output_compacted) {
                    ostr << fmt::format("S\t{}\t{}\n", path.back(), unitig);
                    graph->adjacent_incoming_nodes(path.front(), [&](uint64_t node) {
                        ostr << fmt::format("L\t{}\t+\t{}\t+\t{}M\n",
                                            node, path.back(), overlap);
                    });
                } else {
                    for (size_t i = 0; i < path.size(); ++i) {
                        ostr << fmt::format("S\t{}\t{}\n", path[i],
                                            std::string_view(unitig.data() + i, k));
                        if (i)
                            ostr << fmt::format("L\t{}\t+\t{}\t+\t{}M\n",
                                                path[i - 1], path[i], overlap);
                    }
                    graph->adjacent_incoming_nodes(path.front(), [&](uint64_t node) {
                        ostr << fmt::format("L\t{}\t+\t{}\t+\t{}M\n",
                                            node, path.front(), overlap);
                    });
                }
                std::lock_guard<std::mutex> lock(str_mutex);
                gfa_file << ostr.str();
            },
            get_num_threads(),
            config->min_tip_size
        );

        return 0;
    }

    seq_io::FastaWriter writer(config->outfbase, config->header,
                               config->enumerate_out_sequences,
                               get_num_threads() > 1);
    std::mutex write_mutex;

    if (config->unitigs || config->min_tip_size > 1) {
        graph->call_unitigs([&](const auto &unitig, auto&&) {
                                std::lock_guard<std::mutex> lock(write_mutex);
                                writer.write(unitig);
                            },
                            get_num_threads(),
                            config->min_tip_size,
                            config->kmers_in_single_form);
    } else {
        graph->call_sequences([&](const auto &contig, auto&&) {
                                  std::lock_guard<std::mutex> lock(write_mutex);
                                  writer.write(contig);
                              },
                              get_num_threads(),
                              config->kmers_in_single_form);
    }

    logger->trace("Sequences extracted in {} sec", timer.elapsed());

    return 0;
}

} // namespace cli
} // namespace mtg
