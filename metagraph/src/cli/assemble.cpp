#include "assemble.hpp"

#include <json/json.h>
#include <spdlog/fmt/fmt.h>
#include <tsl/hopscotch_set.h>

#include "common/algorithms.hpp"
#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/threads/threading.hpp"
#include "seq_io/sequence_io.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotated_graph.hpp"
#include "graph/annotated_graph_algorithm.hpp"
#include "graph/representation/masked_graph.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using mtg::graph::DeBruijnGraph;
using mtg::graph::MaskedDeBruijnGraph;
using mtg::graph::AnnotatedDBG;
using mtg::graph::DifferentialAssemblyConfig;


void check_labels(const tsl::hopscotch_set<std::string> &label_set,
                  const AnnotatedDBG &anno_graph) {
    bool detected_missing_labels = false;
    for (const std::string &label : label_set) {
        if (!anno_graph.label_exists(label)) {
            detected_missing_labels = true;
            logger->error("Label {} is not found in annotation", label);
        }
    }

    if (detected_missing_labels)
        exit(1);
}

DifferentialAssemblyConfig diff_assembly_config(const Json::Value &experiment,
                                                const DeBruijnGraph &graph) {
    DifferentialAssemblyConfig diff_config;
    diff_config.add_complement = graph.get_mode() == DeBruijnGraph::CANONICAL;
    diff_config.label_mask_in_kmer_fraction = experiment.get("in_min_fraction", 1.0).asDouble();
    diff_config.label_mask_in_unitig_fraction = experiment.get("unitig_in_min_fraction", 0.0).asDouble();
    diff_config.label_mask_out_kmer_fraction = experiment.get("out_max_fraction", 0.0).asDouble();
    diff_config.label_mask_out_unitig_fraction = experiment.get("unitig_out_max_fraction", 1.0).asDouble();
    diff_config.label_mask_other_unitig_fraction = experiment.get("unitig_other_max_fraction", 1.0).asDouble();

    logger->trace("Per-kmer mask in fraction:\t\t{}", diff_config.label_mask_in_kmer_fraction);
    logger->trace("Per-unitig mask in fraction:\t\t{}", diff_config.label_mask_in_unitig_fraction);
    logger->trace("Per-kmer mask out fraction:\t\t{}", diff_config.label_mask_out_kmer_fraction);
    logger->trace("Per-unitig mask out fraction:\t\t{}", diff_config.label_mask_out_unitig_fraction);
    logger->trace("Per-unitig other label fraction:\t{}", diff_config.label_mask_other_unitig_fraction);
    logger->trace("Include reverse complements:\t\t{}", diff_config.add_complement);

    return diff_config;
}

typedef std::function<void(const graph::MaskedDeBruijnGraph&,
                           const std::string& /* header */)> CallMaskedGraphHeader;

void call_masked_graphs(const AnnotatedDBG &anno_graph,
                        Config *config,
                        const CallMaskedGraphHeader &callback) {
    if (!std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr()).get())
        throw std::runtime_error("Masking only supported for DeBruijnGraph");

    assert(!config->assembly_config_file.empty());

    std::ifstream fin(config->assembly_config_file);
    if (!fin.good())
        throw std::iostream::failure("Failed to read assembly config JSON from " + config->assembly_config_file);

    size_t num_threads = std::max(1u, get_num_threads());

    Json::Value diff_json;
    fin >> diff_json;

    tsl::hopscotch_set<std::string> foreground_labels;
    tsl::hopscotch_set<std::string> background_labels;
    tsl::hopscotch_set<std::string> shared_foreground_labels;
    tsl::hopscotch_set<std::string> shared_background_labels;

    for (const Json::Value &group : diff_json["groups"]) {
        if (group["shared_labels"]) {
            shared_foreground_labels.clear();
            shared_background_labels.clear();

            for (const Json::Value &in_label : group["shared_labels"]["in"]) {
                shared_foreground_labels.emplace(in_label.asString());
            }

            for (const Json::Value &out_label : group["shared_labels"]["out"]) {
                shared_background_labels.emplace(out_label.asString());
            }

            check_labels(shared_foreground_labels, anno_graph);
            check_labels(shared_background_labels, anno_graph);
        }

        if (!group["experiments"])
            throw std::runtime_error("Missing experiments in group");

        for (const Json::Value &experiment : group["experiments"]) {
            DifferentialAssemblyConfig diff_config = diff_assembly_config(
                experiment, anno_graph.get_graph()
            );

            foreground_labels.clear();
            background_labels.clear();

            for (const Json::Value &in_label : experiment["in"]) {
                foreground_labels.emplace(in_label.asString());
            }

            for (const Json::Value &out_label : experiment["out"]) {
                background_labels.emplace(out_label.asString());
            }

            check_labels(foreground_labels, anno_graph);
            check_labels(background_labels, anno_graph);

            std::string exp_name = experiment["name"].asString()
                                    + (config->enumerate_out_sequences ? "." : "");

            logger->trace("Running assembly: {}", exp_name);

            callback(*mask_nodes_by_label(anno_graph,
                                          foreground_labels,
                                          background_labels,
                                          shared_foreground_labels,
                                          shared_background_labels,
                                          diff_config, num_threads), exp_name);
        }
    }
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
        config->infbase = files.at(0);

        assert(config->assembly_config_file.size());
        auto anno_graph = initialize_annotated_dbg(graph, *config);

        logger->trace("Generating masked graphs...");

        std::filesystem::remove(
            utils::remove_suffix(config->outfbase, ".gz", ".fasta") + ".fasta.gz"
        );

        std::mutex write_mutex;

        size_t num_threads = std::max(1u, get_num_threads());

        call_masked_graphs(*anno_graph, config,
            [&](const graph::MaskedDeBruijnGraph &graph, const std::string &header) {
                seq_io::FastaWriter writer(config->outfbase, header,
                                           config->enumerate_out_sequences,
                                           num_threads > 1, /* async write */
                                           "a" /* append mode */);

                if (config->unitigs || config->min_tip_size > 1) {
                    graph.call_unitigs([&](const std::string &unitig, auto&&) {
                                           std::lock_guard<std::mutex> lock(write_mutex);
                                           writer.write(unitig);
                                       },
                                       num_threads, config->min_tip_size,
                                       config->kmers_in_single_form);
                } else {
                    graph.call_sequences([&](const std::string &seq, auto&&) {
                                             std::lock_guard<std::mutex> lock(write_mutex);
                                             writer.write(seq);
                                         },
                                         num_threads, config->kmers_in_single_form);
                }
            }
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

        std::ofstream gfa_file(utils::make_suffix(config->outfbase, ".gfa"));
        std::mutex str_mutex;

        gfa_file << "H\tVN:Z:1.0" << std::endl;
        size_t k = graph->get_k();
        size_t overlap = k - 1;
        graph->call_unitigs(
            [&](const auto &unitig, const auto &path) {
                std::ostringstream ostr;
                if (config->output_compacted) {
                    ostr << fmt::format("S\t{}\t{}\n", path.back(), unitig);
                    graph->call_incoming_kmers(path.front(), [&](uint64_t node, char c) {
                        if (c != graph::boss::BOSS::kSentinel) {
                            ostr << fmt::format("L\t{}\t+\t{}\t+\t{}M\n",
                                                node, path.back(), overlap);
                        }
                    });
                } else {
                    for (size_t i = 0; i < path.size(); ++i) {
                        ostr << fmt::format("S\t{}\t{}\n", path[i],
                                            std::string_view(unitig.data() + i, k));
                        if (i)
                            ostr << fmt::format("L\t{}\t+\t{}\t+\t{}M\n",
                                                path[i - 1], path[i], overlap);
                    }
                    graph->call_incoming_kmers(path.front(), [&](uint64_t node, char c) {
                        if (c != graph::boss::BOSS::kSentinel) {
                            ostr << fmt::format("L\t{}\t+\t{}\t+\t{}M\n",
                                                node, path.front(), overlap);
                        }
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
