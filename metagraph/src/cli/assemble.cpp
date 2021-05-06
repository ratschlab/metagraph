#include "assemble.hpp"

#include <spdlog/fmt/fmt.h>

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


namespace mtg {
namespace cli {

using mtg::common::logger;
using mtg::graph::DeBruijnGraph;
using mtg::graph::MaskedDeBruijnGraph;
using mtg::graph::AnnotatedDBG;
using mtg::graph::DifferentialAssemblyConfig;


void check_and_sort_labels(const AnnotatedDBG &anno_graph,
                     std::vector<std::string> &label_set) {
    bool detected_missing_labels = false;
    for (const std::string &label : label_set) {
        if (!anno_graph.label_exists(label)) {
            detected_missing_labels = true;
            logger->trace("Label {} is not found in annotation", label);
        }
    }

    if (detected_missing_labels)
        exit(1);

    std::sort(label_set.begin(), label_set.end());
}


std::unique_ptr<MaskedDeBruijnGraph>
mask_graph_from_labels(const AnnotatedDBG &anno_graph,
                       const std::vector<std::string> &label_mask_in,
                       const std::vector<std::string> &label_mask_out,
                       const std::vector<std::string> &label_mask_in_post,
                       const std::vector<std::string> &label_mask_out_post,
                       const DifferentialAssemblyConfig &diff_config,
                       size_t num_threads) {
    auto graph = std::dynamic_pointer_cast<const DeBruijnGraph>(anno_graph.get_graph_ptr());

    if (!graph.get())
        throw std::runtime_error("Masking only supported for DeBruijnGraph");

    std::vector<const std::vector<std::string>*> label_sets {
        &label_mask_in, &label_mask_out,
        &label_mask_in_post, &label_mask_out_post
    };

    bool has_overlap = false;
    for (const auto *label_set : label_sets) {
        for (const auto *other_label_set : label_sets) {
            if (label_set == other_label_set)
                continue;

            if (utils::count_intersection(label_set->begin(), label_set->end(),
                                          other_label_set->begin(), other_label_set->end())) {
                has_overlap = true;
                break;
            }
        }

        if (has_overlap)
            break;
    }

    if (has_overlap)
        logger->warn("Overlapping label sets");

    logger->trace("Masked in: {}", fmt::join(label_mask_in, " "));
    logger->trace("Masked in (post-processing): {}", fmt::join(label_mask_in_post, " "));
    logger->trace("Masked out: {}", fmt::join(label_mask_out, " "));
    logger->trace("Masked out (post-processing): {}", fmt::join(label_mask_out_post, " "));

    return std::make_unique<MaskedDeBruijnGraph>(mask_nodes_by_label(
        anno_graph,
        label_mask_in, label_mask_out,
        label_mask_in_post, label_mask_out_post,
        diff_config, num_threads
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

void call_masked_graphs(const AnnotatedDBG &anno_graph,
                        Config *config,
                        const std::vector<std::string> &lines,
                        const CallMaskedGraphHeader &callback,
                        size_t num_threads) {
    std::vector<std::string> shared_foreground_labels;
    std::vector<std::string> shared_background_labels;

    for (const std::string &line : lines) {
        assert(line.size() && line[0] != '#');

        auto line_split = utils::split_string(line, "\t", false);

        if (line[0] == '@') {
            logger->trace("Counting shared k-mers");

            // shared in and out labels
            if (line_split.size() <= 1 || line_split.size() > 3)
                throw std::iostream::failure("Each line in mask file must have 2-3 fields.");

            shared_foreground_labels = utils::split_string(line_split[1], ",");
            shared_background_labels = utils::split_string(
                line_split.size() == 3 ? line_split[2] : "",
                ","
            );

            check_and_sort_labels(anno_graph, shared_foreground_labels);
            check_and_sort_labels(anno_graph, shared_background_labels);

            continue;
        }

        if (line_split.size() <= 2 || line_split.size() > 4)
            throw std::iostream::failure("Each line in mask file must have 3-4 fields.");

        DifferentialAssemblyConfig diff_config = parse_diff_config(
            line_split[1],
            anno_graph.get_graph().get_mode() == DeBruijnGraph::CANONICAL
        );

        if (config->enumerate_out_sequences)
            line_split[0] += ".";

        auto foreground_labels = utils::split_string(line_split[2], ",");
        auto background_labels = utils::split_string(
            line_split.size() == 4 ? line_split[3] : "",
            ","
        );

        check_and_sort_labels(anno_graph, foreground_labels);
        check_and_sort_labels(anno_graph, background_labels);

        callback(*mask_graph_from_labels(anno_graph,
                                         foreground_labels, background_labels,
                                         shared_foreground_labels,
                                         shared_background_labels,
                                         diff_config, num_threads),
                 line_split[0]);
    }
}

std::pair<std::vector<std::string>, unsigned int>
parse_diff_file(const std::string &fname) {
    std::ifstream fin(fname);
    if (!fin.good()) {
        throw std::iostream::failure("Failed to read label mask file");
        exit(1);
    }

    std::string line;
    std::vector<std::string> lines;
    unsigned int num_experiments = 0;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        num_experiments += (line[0] != '@');
        lines.emplace_back(std::move(line));
    }

    return { std::move(lines), num_experiments };
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

        assert(config->label_mask_file.size());
        auto anno_graph = initialize_annotated_dbg(graph, *config);

        logger->trace("Generating masked graphs...");

        std::filesystem::remove(
            utils::remove_suffix(config->outfbase, ".gz", ".fasta") + ".fasta.gz"
        );

        std::mutex write_mutex;

        assert(!config->label_mask_file.empty());

        auto [lines, num_experiments] = parse_diff_file(config->label_mask_file);
        if (!num_experiments)
            return 0;

        size_t num_threads = std::max(1u, get_num_threads());

        call_masked_graphs(*anno_graph, config, lines,
            [&](const graph::MaskedDeBruijnGraph &graph, const std::string &header) {
                seq_io::FastaWriter writer(config->outfbase, header,
                                           config->enumerate_out_sequences,
                                           num_threads > 1, /* async write */
                                           "a" /* append mode */);

                if (config->unitigs || config->min_tip_size > 1) {
                    graph.call_unitigs([&](const auto &unitig, auto&&) {
                                           std::lock_guard<std::mutex> lock(write_mutex);
                                           writer.write(unitig);
                                       },
                                       num_threads,
                                       config->min_tip_size,
                                       config->kmers_in_single_form);
                } else {
                    graph.call_sequences([&](const auto &seq, auto&&) {
                                             std::lock_guard<std::mutex> lock(write_mutex);
                                             writer.write(seq);
                                         },
                                         num_threads,
                                         config->kmers_in_single_form);
                }
            },
            num_threads
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
