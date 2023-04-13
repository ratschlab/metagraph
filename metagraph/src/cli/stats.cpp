#include "stats.hpp"

#include <ips4o/ips4o.hpp>

#include "common/algorithms.hpp"
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

void print_annotation_stats(const std::string &fname, const Config &config) {
    if (config.print_column_names) {
        annot::LabelEncoder<std::string> label_encoder;

        logger->info("Scanning annotation '{}'", fname);

        try {
            std::ifstream instream(fname, std::ios::binary);

            // TODO: make this more reliable
            if (parse_annotation_type(fname) == Config::ColumnCompressed) {
                // Column compressed dumps the number of rows first
                // skipping it...
                load_number(instream);
            }

            if (!label_encoder.load(instream))
                throw std::ios_base::failure("");

        } catch (...) {
            logger->error("Cannot read label encoder from file '{}'", fname);
            exit(1);
        }

        std::cout << "Number of columns: " << label_encoder.size() << std::endl;
        for (size_t c = 0; c < label_encoder.size(); ++c) {
            std::cout << label_encoder.decode(c) << '\n';
        }

        return;
    }

    logger->info("Statistics for annotation '{}'", fname);

    std::unique_ptr<annot::MultiLabelEncoded<std::string>> anno_p
            = initialize_annotation(fname, config);
    auto &annotation = *anno_p;

    if (!annotation.load(fname)) {
        logger->error("Cannot load annotations from file '{}'", fname);
        exit(1);
    }

    using RowDiffCol = annot::binmat::RowDiff<annot::binmat::ColumnMajor>;
    if (auto *rd = dynamic_cast<const RowDiffCol *>(&annotation.get_matrix())) {
        std::string anchors_fname = utils::make_suffix(config.infbase,
                                                       annot::binmat::kRowDiffAnchorExt);
        if (!config.infbase.empty() && std::filesystem::exists(anchors_fname)) {
            const_cast<RowDiffCol *>(rd)->load_anchor(anchors_fname);
        }
    }

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
    CHECK_IF_DIFFED_AND_PRINT_STATS(RowDiff<RowFlat<>>, "RowFlat");
    CHECK_IF_DIFFED_AND_PRINT_STATS(RowDiff<RowSparse>, "RowSparse");
    CHECK_IF_DIFFED_AND_PRINT_STATS(RowDiff<BRWT>, "Multi-BRWT");

    CHECK_IF_DIFFED_AND_PRINT_STATS(typename annot::IntRowDiffBRWTAnnotator::binary_matrix_type, "Multi-BRWT");

    CHECK_IF_DIFFED_AND_PRINT_STATS(typename annot::RowDiffCoordAnnotator::binary_matrix_type, "ColumnMajor");
    CHECK_IF_DIFFED_AND_PRINT_STATS(typename annot::RowDiffBRWTCoordAnnotator::binary_matrix_type, "Multi-BRWT");

    const auto &coords_fname = utils::remove_suffix(fname, annot::ColumnCompressed<>::kExtension)
                                                        + annot::ColumnCompressed<>::kCoordExtension;

    if (const auto *mat_coord = dynamic_cast<const annot::matrix::MultiIntMatrix *>(mat)) {
        std::cout << "================== COORDINATES STATS ===================" << std::endl;
        std::cout << "coordinates: " << mat_coord->num_attributes() << std::endl;
        mat = &mat_coord->get_binary_matrix();
    } else if (parse_annotation_type(fname) == Config::ColumnCompressed
                                && std::filesystem::exists(coords_fname)) {
        std::unique_ptr<std::ifstream> in = utils::open_ifstream(coords_fname);
        if (!*in) {
            logger->error("Could not open file with coordinates {}", coords_fname);
            exit(1);
        }

        uint64_t num_coords = 0;
        try {
            for (size_t j = 0; j < annotation.num_labels(); ++j) {
                {
                    bit_vector_smart delims;
                    delims.load(*in);
                    num_coords += delims.size() - delims.num_set_bits();
                }
                if (j + 1 < annotation.num_labels()) {
                    sdsl::int_vector<> coords;
                    coords.load(*in);
                }
            }
            std::cout << "================== COORDINATES STATS ===================" << std::endl;
            std::cout << "coordinates: " << num_coords << std::endl;
        } catch (...) {
            logger->error("Couldn't read coordinates from {}", coords_fname);
            exit(1);
        }
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

    if (config.print_counts_hist || config.count_quantiles.size()) {
        std::cout << "===================== COUNTS STATS =====================" << std::endl;
        if (!dynamic_cast<const annot::ColumnCompressed<> *>(&annotation)) {
            logger->error("Printing statistics for counts is currently only supported for"
                          " ColumnCompressed ({}) annotations with counts",
                          Config::annotype_to_string(Config::ColumnCompressed));
            exit(1);
        }
        // the annotation is not needed anymore -> free memory
        anno_p.reset();

        std::cout << fmt::format("Column-index\tLabel\tNum-counts");

        for (double q : config.count_quantiles) {
            if (q < 0. or q > 1.) {
                logger->error("Count quantiles must be in interval [0, 1], got {}", q);
                exit(1);
            }
            std::cout << fmt::format("\tQuantile({})", q);
        }
        if (config.print_counts_hist)
            std::cout << fmt::format("\tHistogram(count:multiplicity[,...])");

        std::cout << std::endl;

        annot::ColumnCompressed<>::load_column_values({ fname },
            [&](size_t j, const std::string &label, sdsl::int_vector<>&& counts) {

                // compute count histogram
                std::vector<std::pair<uint64_t, uint64_t>> count_hist_v;
                {
                    tsl::hopscotch_map<uint64_t, uint64_t> counts_hist;
                    for (auto c : counts) {
                        counts_hist[c]++;
                    }
                    count_hist_v.assign(counts_hist.begin(), counts_hist.end());
                }

                ips4o::parallel::sort(count_hist_v.begin(), count_hist_v.end(),
                                      utils::LessFirst(), get_num_threads());

                std::cout << fmt::format("{}\t{}\t{}", j, label, counts.size());

                for (double q : config.count_quantiles) {
                    if (count_hist_v.size()) {
                        std::cout << fmt::format("\t{}", utils::get_quantile(count_hist_v, q));
                    } else {
                        std::cout << "\tnan";
                    }
                }

                if (config.print_counts_hist) {
                    std::cout << "\t";
                    if (count_hist_v.size())
                        std::cout << fmt::format("{}:{}", count_hist_v[0].first, count_hist_v[0].second);
                    for (size_t i = 2; i < count_hist_v.size(); i++) {
                        std::cout << fmt::format(",{}:{}", count_hist_v[i].first, count_hist_v[i].second);
                    }
                }

                std::cout << std::endl;
            }
        );
    }

    std::cout << "========================================================" << std::endl;
}

int print_stats(Config *config) {
    assert(config);

    const auto &files = config->fnames;

    for (const auto &file : files) {
        if (utils::ends_with(file, ".annodbg")) {
            print_annotation_stats(file, *config);
            continue;
        }

        std::shared_ptr<graph::DeBruijnGraph> graph;

        graph = load_critical_dbg(file);
        graph->load_extension<graph::NodeWeights>(file);

        logger->info("Statistics for graph '{}'", file);

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
        print_annotation_stats(file, *config);
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
