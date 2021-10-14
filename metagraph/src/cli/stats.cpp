#include "stats.hpp"

#include <random>

#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/serialization.hpp"
#include "common/threads/threading.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss.hpp"
#include "graph/graph_extensions/node_weights.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "config/config.hpp"
#include "load/load_graph.hpp"
#include "load/load_annotation.hpp"


namespace mtg {
namespace cli {

using mtg::common::logger;
using mtg::common::get_verbose;

typedef annot::MultiLabelEncoded<std::string> Annotator;

void print_degree_stats(const graph::DeBruijnGraph &graph,
                        size_t sample_length,
                        size_t num_samples = 1000) {
    std::cout << "====================== FORK STATS ======================" << std::endl;
    std::random_device rd;
    std::mt19937_64 gen(42);
    std::uniform_int_distribution<uint64_t> distrib(1, graph.max_index());
    double mean_average_outdegree = 0;
    #pragma omp parallel for num_threads(get_num_threads())
    for (size_t i = 0; i < num_samples; ++i) {
        graph::DeBruijnGraph::node_index start;

        #pragma omp critical
        start = distrib(gen);

        std::string seq = graph.get_node_sequence(start);

        size_t dist = 0;
        size_t num_nodes = 0;
        while (dist < sample_length) {
            ++dist;
            std::vector<std::pair<graph::DeBruijnGraph::node_index, char>> nodes;
            graph.call_outgoing_kmers(start, [&](graph::DeBruijnGraph::node_index next, char c) {
                if (c != graph::boss::BOSS::kSentinel)
                    nodes.emplace_back(next, c);
            });
            num_nodes += nodes.size();
            if (nodes.empty()) {
                break;
            } else {
                std::uniform_int_distribution<size_t> out_distrib(0, nodes.size() - 1);
                size_t idx;

                #pragma omp critical
                idx = out_distrib(gen);

                start = nodes[idx].first;
                seq += nodes[idx].second;
            }
        }

        double average_outdegree = static_cast<double>(num_nodes) / dist;

        #pragma omp critical
        {
            std::cout << average_outdegree << "\t" << seq << std::endl;
            mean_average_outdegree += average_outdegree;
        }
    }

    std::cout << "sequence length: " << sample_length << "\n";
    std::cout << "number of samples: " << num_samples << "\n";
    std::cout << "mean avg. outdegree: " << mean_average_outdegree / num_samples << "\n";
}


void print_boss_stats(const graph::boss::BOSS &boss_graph,
                      bool count_dummy,
                      size_t num_threads) {
    std::cout << "====================== BOSS STATS ======================" << std::endl;
    std::cout << "k: " << boss_graph.get_k() + 1 << std::endl;
    std::cout << "nodes (k-1): " << boss_graph.num_nodes() << std::endl;
    std::cout << "edges ( k ): " << boss_graph.num_edges() << std::endl;
    std::cout << "state: " << Config::state_to_string(boss_graph.get_state()) << std::endl;

    assert(boss_graph.rank_W(boss_graph.num_edges(), boss_graph.alph_size) == 0);
    std::cout << "W stats: {'" << boss_graph.decode(0) << "': "
              << boss_graph.rank_W(boss_graph.num_edges(), 0);
    for (int i = 1; i < boss_graph.alph_size; ++i) {
        std::cout << ", '" << boss_graph.decode(i) << "': "
                  << boss_graph.rank_W(boss_graph.num_edges(), i)
                        + boss_graph.rank_W(boss_graph.num_edges(), i + boss_graph.alph_size);
    }
    std::cout << "}" << std::endl;

    assert(boss_graph.get_F(0) == 0);
    std::cout << "F stats: {'";
    for (int i = 1; i < boss_graph.alph_size; ++i) {
        std::cout << boss_graph.decode(i - 1) << "': "
                  << boss_graph.get_F(i) - boss_graph.get_F(i - 1)
                  << ", '";
    }
    std::cout << boss_graph.decode(boss_graph.alph_size - 1) << "': "
              << boss_graph.num_edges() - boss_graph.get_F(boss_graph.alph_size - 1)
              << "}" << std::endl;

    if (count_dummy) {
        uint64_t num_source_dummy_edges
            = boss_graph.mark_source_dummy_edges(NULL, num_threads);
        uint64_t num_sink_dummy_edges = boss_graph.mark_sink_dummy_edges(NULL);

        std::cout << "dummy source edges: " << num_source_dummy_edges << std::endl;
        std::cout << "dummy sink edges: " << num_sink_dummy_edges << std::endl;
        std::cout << "real edges: "
                  << boss_graph.num_edges() - num_source_dummy_edges - num_sink_dummy_edges
                  << std::endl;
    }
    std::cout << "indexed suffix length: " << boss_graph.get_indexed_suffix_length() << std::endl;
    std::cout << "========================================================" << std::endl;
}

void print_stats(const graph::DeBruijnGraph &graph, bool print_counts_hist) {
    std::cout << "====================== GRAPH STATS =====================" << std::endl;
    std::cout << "k: " << graph.get_k() << std::endl;
    std::cout << "nodes (k): " << graph.num_nodes() << std::endl;
    std::cout << "mode: " << Config::graphmode_to_string(graph.get_mode()) << std::endl;

    if (auto weights = graph.get_extension<graph::NodeWeights>()) {
        double sum_weights = 0;
        uint64_t num_non_zero_weights = 0;
        std::vector<uint64_t> hist;
        auto add_to_hist = [&](uint64_t c) {
            assert(c && "All k-mers in graph must have non-zero weights");
            while (c >= hist.size()) {
                hist.push_back(0);
            }
            hist[c]++;
        };
        if (const auto *dbg_succ = dynamic_cast<const graph::DBGSuccinct*>(&graph)) {
            // In DBGSuccinct some of the nodes may be masked out
            // TODO: Fix this by using non-contiguous indexing in graph
            //       so that mask of dummy edges does not change indexes.
            for (uint64_t i = 1; i <= dbg_succ->get_boss().num_edges(); ++i) {
                if (uint64_t weight = (*weights)[i]) {
                    sum_weights += weight;
                    num_non_zero_weights++;
                    if (print_counts_hist)
                        add_to_hist(weight);
                }
            }
        } else {
            if (!weights->is_compatible(graph)) {
                logger->error("Node weights are not compatible with graph");
                exit(1);
            }
            graph.call_nodes([&](auto i) {
                if (uint64_t weight = (*weights)[i]) {
                    sum_weights += weight;
                    num_non_zero_weights++;
                    if (print_counts_hist)
                        add_to_hist(weight);
                }
            });
        }
        std::cout << "nnz weights: " << num_non_zero_weights << std::endl;
        std::cout << "avg weight: " << static_cast<double>(sum_weights) / num_non_zero_weights << std::endl;
        if (print_counts_hist) {
            std::cout << "weights histogram:\n";
            if (hist.size() > 1u && hist[1])
                std::cout << fmt::format("{}:{}", 1, hist[1]);
            for (size_t i = 2; i < hist.size(); i++) {
                if (hist[i])
                    std::cout << fmt::format(",{}:{}", i, hist[i]);
            }
            std::cout << std::endl;
        }

        if (get_verbose()) {
            if (const auto *dbg_succ = dynamic_cast<const graph::DBGSuccinct*>(&graph)) {
                // In DBGSuccinct some of the nodes may be masked out
                // TODO: Fix this by using non-contiguous indexing in graph
                //       so that mask of dummy edges does not change indexes.
                for (uint64_t i = 1; i <= dbg_succ->get_boss().num_edges(); ++i) {
                    if (uint64_t weight = (*weights)[i])
                        std::cout << weight << " ";
                }
            } else {
                graph.call_nodes([&](auto i) { std::cout << (*weights)[i] << " "; });
            }
            std::cout << std::endl;
        }
    }

    std::cout << "========================================================" << std::endl;
}

template <class KmerHasher>
void print_bloom_filter_stats(const kmer::KmerBloomFilter<KmerHasher> *kmer_bloom) {
    if (!kmer_bloom)
        return;

    std::cout << "====================== BLOOM STATS =====================" << std::endl;
    std::cout << "Size (bits):\t" << kmer_bloom->size() << std::endl
              << "Num hashes:\t" << kmer_bloom->num_hash_functions() << std::endl;
    std::cout << "========================================================" << std::endl;
}

template <class Matrix>
void print_anchor_stats(const Matrix& m) {
    uint64_t num_anchors = m.anchor().num_set_bits();
    if (num_anchors != 0) {
        std::cout << "num anchors: " << m.anchor().num_set_bits() << std::endl;
    } else {
        std::cout << "Please specify the anchor file via '-i anchors_file' to get "
                     "anchor stats" << std::endl;
    }
}

void print_brwt_stats(const annot::binmat::BRWT& brwt) {
    std::cout << "=================== Multi-BRWT STATS ===================" << std::endl;
    std::cout << "num nodes: " << brwt.num_nodes() << std::endl;
    std::cout << "avg arity: " << brwt.avg_arity() << std::endl;
    std::cout << "shrinkage: " << brwt.shrinking_rate() << std::endl;
    if (get_verbose()) {
        std::cout << "==================== Multi-BRWT TREE ===================" << std::endl;
        brwt.print_tree_structure(std::cout);
    }
}

void print_stats(const Annotator &annotation) {
    std::cout << "=================== ANNOTATION STATS ===================" << std::endl;
    std::cout << "labels:  " << annotation.num_labels() << std::endl;
    std::cout << "objects: " << annotation.num_objects() << std::endl;
    std::cout << "density: " << static_cast<double>(annotation.num_relations())
                                    / annotation.num_objects()
                                    / annotation.num_labels() << std::endl;
    std::cout << "representation: "
              << utils::split_string(annotation.file_extension(), ".").at(0) << std::endl;

    using namespace annot::binmat;

    const BinaryMatrix *mat = &annotation.get_matrix();

#define CHECK_IF_DIFFED_AND_PRINT_STATS(RD_TYPE, NAME) \
    if (const auto *rd = dynamic_cast<const RD_TYPE *>(mat)) { \
        std::cout << "=================== DIFF ANNOTATION ====================" << std::endl; \
        print_anchor_stats(*rd); \
        std::cout << "underlying matrix: " NAME << std::endl; \
        mat = &rd->diffs(); \
    }

    CHECK_IF_DIFFED_AND_PRINT_STATS(RowDiff<ColumnMajor>, "ColumnMajor");
    CHECK_IF_DIFFED_AND_PRINT_STATS(RowDiff<RowSparse>, "RowSparse");
    CHECK_IF_DIFFED_AND_PRINT_STATS(RowDiff<BRWT>, "Multi-BRWT");

    CHECK_IF_DIFFED_AND_PRINT_STATS(typename annot::IntRowDiffBRWTAnnotator::binary_matrix_type, "Multi-BRWT");

    CHECK_IF_DIFFED_AND_PRINT_STATS(typename annot::RowDiffCoordAnnotator::binary_matrix_type, "ColumnMajor");
    CHECK_IF_DIFFED_AND_PRINT_STATS(typename annot::RowDiffBRWTCoordAnnotator::binary_matrix_type, "Multi-BRWT");

    if (const auto *mat_coord = dynamic_cast<const annot::matrix::MultiIntMatrix *>(mat)) {
        std::cout << "================== COORDINATES STATS ===================" << std::endl;
        std::cout << "coordinates: " << mat_coord->num_attributes() << std::endl;
        mat = &mat_coord->get_binary_matrix();
    } else if (const auto *int_mat = dynamic_cast<const annot::matrix::IntMatrix *>(mat)) {
        mat = &int_mat->get_binary_matrix();
    }

    if (const auto *rbmat = dynamic_cast<const RainbowMatrix *>(mat)) {
        std::cout << "================= RAINBOW MATRIX STATS =================" << std::endl;
        std::cout << "distinct rows: " << rbmat->num_distinct_rows() << std::endl;
        if (const auto *rb_brwt = dynamic_cast<const Rainbow<BRWT> *>(mat))
            mat = &rb_brwt->get_reduced_matrix();
    }

    if (const auto *brwt = dynamic_cast<const BRWT *>(mat))
        print_brwt_stats(*brwt);

    std::cout << "========================================================" << std::endl;
}

int print_stats(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    for (const auto &file : files) {
        std::shared_ptr<graph::DeBruijnGraph> graph;

        graph = load_critical_dbg(file);
        graph->load_extension<graph::NodeWeights>(file);

        logger->info("Statistics for graph '{}'", file);

        if (config->print_degree_stats) {
            print_degree_stats(*graph, config->distance);
            return 0;
        }

        print_stats(*graph, config->print_counts_hist);

        if (auto dbg_succ = dynamic_cast<graph::DBGSuccinct*>(graph.get())) {
            const auto &boss_graph = dbg_succ->get_boss();

            print_boss_stats(boss_graph, config->count_dummy, get_num_threads());

            if (config->print_graph_internal_repr) {
                logger->info("Printing internal representation");
                boss_graph.print_internal_representation(std::cout);
            }
            print_bloom_filter_stats(dbg_succ->get_bloom_filter());
        }

        if (config->print_graph)
            std::cout << *graph;
    }

    for (const std::string &file : config->infbase_annotators) {
        std::unique_ptr<annot::MultiLabelEncoded<std::string>> annotation
                = initialize_annotation(file, *config);

        if (config->print_column_names) {
            annot::LabelEncoder<std::string> label_encoder;

            logger->info("Scanning annotation '{}'", file);

            try {
                std::ifstream instream(file, std::ios::binary);

                // TODO: make this more reliable
                if (dynamic_cast<const annot::ColumnCompressed<> *>(annotation.get())) {
                    // Column compressed dumps the number of rows first
                    // skipping it...
                    load_number(instream);
                }

                if (!label_encoder.load(instream))
                    throw std::ios_base::failure("");

            } catch (...) {
                logger->error("Cannot read label encoder from file '{}'", file);
                exit(1);
            }

            std::cout << "Number of columns: " << label_encoder.size() << std::endl;
            for (size_t c = 0; c < label_encoder.size(); ++c) {
                std::cout << label_encoder.decode(c) << '\n';
            }

            continue;
        }

        if (!annotation->load(file)) {
            logger->error("Cannot load annotations from file '{}'", file);
            exit(1);
        }

        using RowDiffCol = annot::binmat::RowDiff<annot::binmat::ColumnMajor>;
        if (auto *rd = dynamic_cast<const RowDiffCol *>(&annotation->get_matrix())) {
            std::string anchors_file = utils::make_suffix(config->infbase,
                                                          annot::binmat::kRowDiffAnchorExt);
            if (!config->infbase.empty() && std::filesystem::exists(anchors_file)) {
                const_cast<RowDiffCol *>(rd)->load_anchor(anchors_file);
            }
        }

        logger->info("Statistics for annotation '{}'", file);
        print_stats(*annotation);
    }

    return 0;
}

int compare(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    assert(files.size());

    logger->info("Loading graph                '{}'", files.at(0));
    auto graph = load_critical_dbg(files.at(0));

    for (size_t f = 1; f < files.size(); ++f) {
        logger->info("Loading graph for comparison '{}'", files[f]);
        auto second = load_critical_dbg(files[f]);
        if (*graph == *second) {
            logger->info("Graphs are identical");
        } else {
            logger->info("Graphs are not identical");
        }
    }

    return 0;
}

} // namespace cli
} // namespace mtg
