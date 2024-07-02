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
#include "load/load_annotation.hpp"
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

DifferentialAssemblyConfig diff_assembly_config(const Json::Value &experiment) {
    DifferentialAssemblyConfig diff_config;
    diff_config.count_kmers = experiment.get("count_kmers", false).asBool();
    diff_config.family_wise_error_rate = experiment.get("family_wise_error_rate", 0.05).asDouble();
    diff_config.test_by_unitig = experiment.get("test_by_unitig", false).asBool();
    diff_config.test_type = experiment.get("test_type", "nbinom_exact").asString();
    diff_config.clean = experiment.get("clean", false).asBool();
    diff_config.min_count = experiment.get("min_count", 0).asUInt64();
    diff_config.min_recurrence = experiment.get("min_recurrence", 1).asUInt64();
    diff_config.min_in_recurrence = experiment.get("min_in_recurrence", 0).asUInt64();
    diff_config.min_out_recurrence = experiment.get("min_out_recurrence", 0).asUInt64();
    diff_config.max_in_recurrence = experiment.get("max_in_recurrence", std::numeric_limits<uint64_t>::max()).asUInt64();
    diff_config.max_out_recurrence = experiment.get("max_out_recurrence", std::numeric_limits<uint64_t>::max()).asUInt64();
    diff_config.assemble_shared = experiment.get("assemble_shared", false).asBool();
    return diff_config;
}

typedef std::function<void(const graph::DeBruijnGraph&,
                           const std::string& /* header */)> CallMaskedGraphHeader;

void call_masked_graphs(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                        Config *config,
                        const CallMaskedGraphHeader &callback) {
    assert(!config->assembly_config_file.empty());

    std::ifstream fin(config->assembly_config_file);
    if (!fin.good())
        throw std::iostream::failure("Failed to read assembly config JSON from " + config->assembly_config_file);

    for (const auto &name : config->fnames) {
        if (parse_annotation_type(name) != Config::ColumnCompressed) {
            throw std::runtime_error("Multiple annotation files must be ColumnCompressed");
        }
    }

    size_t num_threads = std::max(1u, get_num_threads());

    Json::Value diff_json;
    fin >> diff_json;

    for (const Json::Value &group : diff_json["groups"]) {
        if (!group["experiments"])
            throw std::runtime_error("Missing experiments in group");

        for (const Json::Value &experiment : group["experiments"]) {
            tsl::hopscotch_set<std::string> foreground_labels;
            tsl::hopscotch_set<std::string> background_labels;

            DifferentialAssemblyConfig diff_config = diff_assembly_config(experiment);
            diff_config.outfbase = config->outfbase;

            std::string exp_name = experiment["name"].asString()
                                    + (config->enumerate_out_sequences ? "." : "");

            for (const Json::Value &in_label : experiment["in"]) {
                foreground_labels.emplace(in_label.asString());
            }

            for (const Json::Value &out_label : experiment["out"]) {
                background_labels.emplace(out_label.asString());
            }

            logger->trace("Running assembly: {}", exp_name);

            auto filenames = (config->infbase_annotators.size() > 1) ? config->infbase_annotators : config->fnames;

            std::shared_ptr<DeBruijnGraph> ingraph;
            std::shared_ptr<DeBruijnGraph> outgraph;

            if (diff_config.test_by_unitig) {
                auto [cur_ingraph, cur_outgraph, pvals, tmp_file] =
                    graph::mask_nodes_by_label_dual<std::vector<uint64_t>>(graph_ptr,
                                        filenames,
                                        foreground_labels,
                                        background_labels,
                                        diff_config, num_threads,
                                        config->tmp_dir,
                                        config->parallel_nodes);
                std::swap(ingraph, cur_ingraph);
                std::swap(outgraph, cur_outgraph);

                std::ofstream fout_all(config->outfbase + ".all.pvals", ios::out | ios::app);
                for (uint64_t pval : pvals) {
                    static_assert(sizeof(double) == sizeof(pval));
                    double pval_d;
                    memcpy(&pval_d, &pval, sizeof(pval));
                    fout_all << pval_d << "\n";
                }
            } else {
                auto [cur_ingraph, cur_outgraph, pvals, tmp_file] =
                    graph::mask_nodes_by_label_dual<sdsl::int_vector_buffer<64>>(graph_ptr,
                                        filenames,
                                        foreground_labels,
                                        background_labels,
                                        diff_config, num_threads,
                                        config->tmp_dir,
                                        config->parallel_nodes);
                std::swap(ingraph, cur_ingraph);
                std::swap(outgraph, cur_outgraph);

                std::ofstream fout_all(config->outfbase + ".all.pvals", ios::out | ios::app);
                for (uint64_t pval : pvals) {
                    static_assert(sizeof(double) == sizeof(pval));
                    double pval_d;
                    memcpy(&pval_d, &pval, sizeof(pval));
                    fout_all << pval_d << "\n";
                }
                pvals.close();
            }

            sdsl::bit_vector mask;
            callback(*ingraph, exp_name + "in.");
            if (diff_config.assemble_shared) {
                if (auto masked_graph = std::dynamic_pointer_cast<MaskedDeBruijnGraph>(ingraph)) {
                    std::unique_ptr<bitmap_vector> inmask(static_cast<bitmap_vector*>(masked_graph->release_mask()));
                    mask = const_cast<sdsl::bit_vector&&>(inmask->data());
                    inmask.reset();
                }
            }
            ingraph.reset();

            callback(*outgraph, exp_name + "out.");
            if (diff_config.assemble_shared) {
                if (auto masked_graph = std::dynamic_pointer_cast<MaskedDeBruijnGraph>(outgraph)) {
                    std::unique_ptr<bitmap_vector> outmask(static_cast<bitmap_vector*>(masked_graph->release_mask()));
                    outmask->add_to(&mask);
                    outmask.reset();
                    mask.flip();
                    auto masked_graph_null = std::make_shared<MaskedDeBruijnGraph>(
                        graph_ptr, std::make_unique<bitmap_vector>(std::move(mask)), true,
                        graph_ptr->get_mode() == DeBruijnGraph::PRIMARY ? DeBruijnGraph::PRIMARY : DeBruijnGraph::BASIC
                    );
                    callback(*masked_graph_null, exp_name + "shared.");
                }
            }
            outgraph.reset();
        }
    }
}

int assemble(Config *config) {
    assert(config);

    if (config->infbase_annotators.size()) {
        assert(config->assembly_config_file.size());
        if (!std::filesystem::exists(config->assembly_config_file)) {
            logger->error("Differential assembly config does not exist\n{}",
                          config->assembly_config_file);
            exit(1);
        }
    }

    const auto &files = config->fnames;

    assert(config->outfbase.size());

    Timer timer;
    logger->trace("Graph loading...");

    auto graph = load_critical_dbg(files.at(0));

    logger->trace("Graph loaded in {} sec", timer.elapsed());

    if (config->infbase_annotators.size() > 1 || files.size() > 1) {
        logger->trace("Generating masked graphs...");
        std::mutex write_mutex;
        size_t num_threads = std::max(1u, get_num_threads());

        std::filesystem::remove(
            utils::remove_suffix(config->outfbase, ".gz", ".fasta") + ".fasta.gz"
        );

        std::filesystem::remove(config->outfbase + ".pvals");
        std::filesystem::remove(config->outfbase + ".all.pvals");

        auto graph_callback = [&](const graph::DeBruijnGraph &graph, const std::string &header) {
            seq_io::FastaWriter writer(config->outfbase,
                                       header,
                                       config->enumerate_out_sequences,
                                       num_threads > 1, /* async write */
                                       "a" /* append mode */);

            if (config->unitigs || config->min_tip_size > 1) {
                graph.call_unitigs(
                    [&](const std::string &unitig, auto&& path) {
                        std::lock_guard<std::mutex> lock(write_mutex);
                        writer.write(unitig);
                    },
                    num_threads, config->min_tip_size,
                    config->kmers_in_single_form
                );
            } else {
                graph.call_sequences(
                    [&](const std::string &seq, auto&& path) {
                        std::lock_guard<std::mutex> lock(write_mutex);
                        writer.write(seq);
                    },
                    num_threads, config->kmers_in_single_form
                );
            }
        };

        config->fnames.erase(config->fnames.begin(), config->fnames.begin() + 1);
        call_masked_graphs(graph, config, graph_callback);

        return 0;
    }

    assert(files.size() == 1);

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
