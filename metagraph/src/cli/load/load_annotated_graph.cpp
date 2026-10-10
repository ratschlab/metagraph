#include "load_annotated_graph.hpp"

#include <cassert>
#include <filesystem>
#include <mutex>
#include <random>

#include "annotation/binary_matrix/multi_brwt/brwt.hpp"
#include "annotation/binary_matrix/column_sparse/column_major.hpp"
#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/binary_matrix/row_sparse/row_sparse.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/coord_to_header.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/annotated_dbg.hpp"
#include "common/logger.hpp"
#include "common/threads/threading.hpp"
#include "common/unix_tools.hpp"
#include "common/utils/file_utils.hpp"
#include "common/utils/string_utils.hpp"
#include "cli/config/config.hpp"
#include "cli/pattern.hpp"
#include "load_graph.hpp"
#include "load_annotation.hpp"


namespace mtg {
namespace cli {

using namespace mtg::graph;
using mtg::common::logger;

namespace {

// The graphs whose mask build_mask_at_load built. Weak references: a graph freed and another
// allocated at its address must not pass for it (the unit tests load many graphs).
std::mutex built_masks_mutex;
std::vector<std::weak_ptr<const DBGSuccinct>> built_masks;

// The graphs whose mask check_mask_at_load found invalid (weak references, as above)
std::vector<std::weak_ptr<const DBGSuccinct>> invalid_masks;

// Whether |registry| holds |graph| (or the PRIMARY graph its CanonicalDBG wraps); prunes the
// expired entries. The caller holds built_masks_mutex.
bool registered(std::vector<std::weak_ptr<const DBGSuccinct>> *registry,
                const DeBruijnGraph &graph) {
    const DeBruijnGraph *base = &graph;
    if (const auto *canonical = dynamic_cast<const CanonicalDBG*>(&graph))
        base = &canonical->get_graph();
    bool found = false;
    for (auto it = registry->begin(); it != registry->end(); ) {
        if (auto held = it->lock()) {
            found |= held.get() == base;
            ++it;
        } else {
            it = registry->erase(it);
        }
    }
    return found;
}

} // namespace

void build_mask_at_load(const std::shared_ptr<DeBruijnGraph> &graph, bool stdout_reserved) {
    const auto progress = stdout_reserved ? spdlog::level::trace : spdlog::level::info;
    auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(graph);
    if (!dbg_succ) {
        logger->warn("--pattern-build-mask: the graph is not a succinct graph: it has no "
                     "dummy-edge mask to build (and the pattern search does not serve it)");
        return;
    }
    if (dbg_succ->get_mask()) {
        logger->log(progress, "--pattern-build-mask: the graph has its dummy-edge mask (its "
                    ".edgemask file): nothing to build");
        return;
    }
    // (the graph was loaded without one: none beside it, or one that could not be opened,
    // which DBGSuccinct::load names in a warning)
    logger->log(progress, "--pattern-build-mask: building the dummy-edge mask in memory "
                          "(the graph was loaded without a .edgemask)...");
    try {
        const DummyMaskCounts counts = mask_dummy_edges(dbg_succ.get(), get_num_threads());
        logger->log(progress, "--pattern-build-mask: dummy-edge mask built in {:.3f} s with {} "
                    "threads: {} edges, {} source dummies (the main dummy edge included), {} "
                    "sink dummies, {} k-mers; it lives in memory only (transform --mask-dummy "
                    "writes it to the graph's .edgemask once)",
                    counts.seconds, get_num_threads(), counts.edges, counts.source_dummy,
                    counts.sink_dummy, counts.kmers);
    } catch (const std::exception &e) {
        logger->error("--pattern-build-mask: the dummy-edge mask could not be built: {}",
                      e.what());
        exit(1);
    }
    std::lock_guard<std::mutex> lock(built_masks_mutex);
    built_masks.emplace_back(dbg_succ);
}

bool mask_built_at_load(const DeBruijnGraph &graph) {
    std::lock_guard<std::mutex> lock(built_masks_mutex);
    return registered(&built_masks, graph);
}

MaskSentinelSample sample_mask_sentinels(const DBGSuccinct &graph, uint64_t samples) {
    MaskSentinelSample result;
    const bit_vector *mask = graph.get_mask();
    if (!mask)
        return result;

    const boss::BOSS &boss = graph.get_boss();
    const auto &W = boss.get_W();
    const auto plain = static_cast<boss::BOSS::TAlphabet>(boss::BOSS::kSentinelCode);
    const auto marked = static_cast<boss::BOSS::TAlphabet>(boss::BOSS::kSentinelCode
                                                            + boss.alph_size);
    // the occurrences of $ and of its marked form in W[1..num_edges] (position 0 holds a
    // placeholder, no edge): two ranks each
    const uint64_t last = W.size() - 1;
    const uint64_t first_plain = W.rank(plain, 0);
    const uint64_t first_marked = W.rank(marked, 0);
    const uint64_t num_plain = W.rank(plain, last) - first_plain;
    const uint64_t num_marked = W.rank(marked, last) - first_marked;
    result.sentinel_edges = num_plain + num_marked;

    // the edge of the r-th occurrence (0-based) of either form, plain ones first
    auto edge_of = [&](uint64_t r) {
        return r < num_plain ? W.select(plain, first_plain + r + 1)
                             : W.select(marked, first_marked + (r - num_plain) + 1);
    };
    auto look_up = [&](uint64_t edge) {
        assert(edge >= 1 && edge < W.size() && W[edge] % boss.alph_size == 0);
        ++result.looked_up;
        result.valid += (*mask)[edge];
    };
    if (result.sentinel_edges <= samples) {
        // all of them: the full check
        for (uint64_t r = 0; r < result.sentinel_edges; ++r) {
            look_up(edge_of(r));
        }
        return result;
    }
    // the main dummy edge, then the draws: the generator the standard specifies, not a
    // distribution of the library (whose draws differ between implementations)
    if (W[1] % boss.alph_size == 0)
        look_up(1);
    std::mt19937_64 rng(boss.num_edges());
    for (uint64_t i = 0; i < samples; ++i) {
        look_up(edge_of(rng() % result.sentinel_edges));
    }
    return result;
}

uint64_t check_mask_at_load(const std::shared_ptr<DeBruijnGraph> &graph, bool stdout_reserved,
                            bool progress) {
    auto dbg_succ = std::dynamic_pointer_cast<DBGSuccinct>(graph);
    if (!dbg_succ || !dbg_succ->get_mask() || mask_built_at_load(*graph))
        return 0;

    Timer timer;
    const MaskSentinelSample sample = sample_mask_sentinels(*dbg_succ);
    logger->log(stdout_reserved || !progress ? spdlog::level::trace : spdlog::level::info,
                "Dummy-edge mask checked for the pattern search in {:.3f} s: {} of {} edges "
                "with W = $ looked up, {} of them marked valid", timer.elapsed(),
                sample.looked_up, sample.sentinel_edges, sample.valid);
    if (!sample.valid)
        return 0;

    logger->warn("The dummy-edge mask (.edgemask) marks dummy edges with W = $ valid ({} of the "
                 "{} looked up; a mask written by `metagraph extend` on a masked graph, or a "
                 "stale one): the pattern search answers mask_invalid. Remedy: rebuild the mask "
                 "with `metagraph transform --mask-dummy --force <graph>.dbg` and restart",
                 sample.valid, sample.looked_up);
    std::lock_guard<std::mutex> lock(built_masks_mutex);
    invalid_masks.emplace_back(dbg_succ);
    return sample.valid;
}

bool mask_invalid_at_load(const DeBruijnGraph &graph) {
    std::lock_guard<std::mutex> lock(built_masks_mutex);
    return registered(&invalid_masks, graph);
}

PatternPreparation pattern_preparation(const Config &config) {
    // where the pattern search is served on a graph loaded from |config|: the CLI and a
    // single-graph server (a graph list's graphs are prepared by the server, which loads them
    // itself). The CLI's stdout carries its answers, and the logger writes info lines there
    PatternPreparation prep;
    const bool cli = config.identity == Config::PATTERN;
    const bool serves_pattern
            = cli || (config.identity == Config::SERVER_QUERY && config.fnames.empty());
    prep.check_mask = serves_pattern;
    prep.build_mask = config.pattern_build_mask;
    prep.notices = serves_pattern;
    prep.sample_fraction = serves_pattern;
    prep.stdout_reserved = cli;
    return prep;
}

void prepare_pattern_graph(const std::shared_ptr<DeBruijnGraph> &graph, const std::string &path,
                           const PatternPreparation &prep) {
    const bool cli = prep.stdout_reserved;
    // a mask read from its file is checked once per load, before the graph is shared
    // (mask_invalid); one built at load is not
    if (prep.check_mask)
        check_mask_at_load(graph, cli, prep.progress);
    if (prep.build_mask) {
        build_mask_at_load(graph, cli);
    } else if (prep.notices) {
        const auto *dbg_succ = dynamic_cast<const DBGSuccinct*>(graph.get());
        const std::string mask_path = utils::remove_suffix(path, DBGSuccinct::kExtension)
                + DBGSuccinct::kDummyMaskExtension;
        if (dbg_succ && !dbg_succ->get_mask() && std::filesystem::exists(mask_path)) {
            // a mask that is there but could not be opened (permissions): not "no mask",
            // whose remedy transform would refuse
            logger->log(cli ? spdlog::level::warn : spdlog::level::info,
                        "The dummy-edge mask {} exists but could not be opened "
                        "(permissions?): the graph was loaded without it, and the pattern "
                        "search counts upper bounds with estimates (counting: "
                        "upper_bound). For exact counts: make it readable and restart, or "
                        "`metagraph transform --mask-dummy --force {}`, or "
                        "--pattern-build-mask (builds it in memory at every start)",
                        mask_path, path);
        } else if (dbg_succ && !dbg_succ->get_mask()) {
            // the operator learns at start-up, not from the first answer (in the server's
            // log; on the CLI's stderr, as a warning)
            logger->log(cli ? spdlog::level::warn : spdlog::level::info,
                        "The graph has no dummy-edge mask (.edgemask): the pattern search "
                        "counts upper bounds with estimates (counting: upper_bound). For "
                        "exact counts: `metagraph transform --mask-dummy {}` once (writes "
                        "the .edgemask beside the graph), or --pattern-build-mask (builds "
                        "it in memory at every start)",
                        path);
        }
    }
    if (prep.sample_fraction) {
        // without a mask (none, or not built): the dummy fraction of the estimates
        // (counts are bounds, DESIGN-pattern-search.md §4.4), sampled once here rather
        // than in the first request
        sample_dummy_fraction_at_load(graph, cli || !prep.progress);
    }
}

std::shared_future<std::shared_ptr<DeBruijnGraph>>
async_load_critical_dbg(const std::string &path, const PatternPreparation &prep) {
    return std::async(std::launch::async, [path, prep]() -> std::shared_ptr<DeBruijnGraph> {
        auto graph = load_critical_dbg(path);
        // here, in the loading thread: while the annotation loads, and before the graph is
        // shared with anyone, so that nothing ever sees it without the mask
        prepare_pattern_graph(graph, path, prep);
        return graph;
    }).share();
}

std::shared_future<std::shared_ptr<DeBruijnGraph>> async_load_critical_dbg(const Config &config) {
    return async_load_critical_dbg(config.infbase, pattern_preparation(config));
}

namespace {

std::unique_ptr<annot::CoordToHeader>
load_coord_to_header(const annot::MultiLabelAnnotation<std::string> &annotation,
                     const Config &config) {
    std::unique_ptr<annot::CoordToHeader> coord_to_header;

    if (dynamic_cast<const annot::matrix::MultiIntMatrix *>(&annotation.get_matrix())
            && !config.no_coord_mapping && config.identity != Config::ANNOTATE) {
        // Load CoordToHeader mapping if exists
        auto cth_fname = utils::remove_suffix(config.infbase_annotators.at(0),
                                              annotation.file_extension())
                                        + annot::CoordToHeader::kExtension;
        if (std::filesystem::exists(cth_fname)) {
            coord_to_header = std::make_unique<annot::CoordToHeader>();
            if (!coord_to_header->load(cth_fname)) {
                exit(1);
            }
            logger->trace("CoordToHeader mapping loaded successfully from {}. "
                          "All queries will be performed against individual sequences.",
                          cth_fname);
        } else {
            const auto anno_basename = utils::remove_suffix(config.infbase_annotators.at(0),
                                                            annotation.file_extension());
            logger->warn("No CoordToHeader mapping found at '{}'. Coords output will use "
                         "file-level positions (e.g., '<file_37.fa>:0-1086-1090') instead of "
                         "per-sequence positions (e.g., '<seq_9>:0-1-5'), and the per-target "
                         "sequence length needed to compute the fraction of the target covered "
                         "('kmers_in_target' for query, 'nt_length' for align) will be omitted. "
                         "To enable per-sequence reporting, run once against the final annotation "
                         "with all input FASTAs:\n"
                         "    metagraph annotate --anno-filename --index-header-coords "
                         "-i {} -o {} <input_fastas>\n"
                         "Pass '--no-coord-mapping' to suppress this warning.",
                         cth_fname, config.infbase, anno_basename);
        }
    }

    return coord_to_header;
}

// Build the AnnotatedDBG from a graph future. Annotation and CTH loads run
// in parallel with the graph load; the graph is awaited only when needed
// (for row_diff set_graph and for the final PRIMARY wrap).
// |shared_coord_to_header|: the mapping of this annotation that another index holds, used
// instead of loading the .seqs (initialize_annotated_dbg); null loads it.
std::unique_ptr<AnnotatedDBG>
build_annotated_dbg(std::shared_future<std::shared_ptr<DeBruijnGraph>> graph_future,
                    const Config &config,
                    size_t max_chunks_open,
                    std::shared_ptr<const annot::CoordToHeader> shared_coord_to_header = nullptr) {
    // Construct the annotation. The bulk load is graph-independent and runs
    // in parallel with the graph load.
    std::unique_ptr<annot::MultiLabelAnnotation<std::string>> annotation;
    std::shared_ptr<const annot::CoordToHeader> coord_to_header;
    if (!config.infbase_annotators.size()) {
        // No annotators configured: build an empty annotation sized to the graph.
        auto graph = graph_future.get();
        annotation = initialize_annotation(config.anno_type, config,
                                           graph->max_index(), max_chunks_open);
    } else {
        annotation = initialize_annotation(config.infbase_annotators.at(0),
                                           config, 0, max_chunks_open);

        bool loaded = false;
        if (auto *cc = dynamic_cast<annot::ColumnCompressed<>*>(annotation.get())) {
            loaded = cc->merge_load(config.infbase_annotators);
        } else {
            if (config.infbase_annotators.size() > 1) {
                logger->warn("Cannot merge annotations of this type. Only the first"
                             " file {} will be loaded.", config.infbase_annotators.at(0));
            }
            loaded = annotation->load(config.infbase_annotators.at(0));
        }
        if (!loaded)
            exit(1);

        // The mapping another index of this annotation holds is shared (a column per label,
        // as the .seqs of this annotation has: a mapping of another annotation is not, and the
        // .seqs is loaded instead); else the CTH load (graph-independent) overlaps with the
        // graph load.
        if (shared_coord_to_header
                && shared_coord_to_header->num_columns() != annotation->num_labels()) {
            logger->warn("The record mapping offered for {} has {} columns, the annotation {} "
                         "labels: not the mapping of this annotation; its .seqs is loaded",
                         config.infbase_annotators.at(0), shared_coord_to_header->num_columns(),
                         annotation->num_labels());
            shared_coord_to_header = nullptr;
        }
        if (shared_coord_to_header) {
            coord_to_header = std::move(shared_coord_to_header);
        } else {
            coord_to_header = load_coord_to_header(*annotation, config);
        }

        using namespace annot::matrix;
        BinaryMatrix &matrix = const_cast<BinaryMatrix &>(annotation->get_matrix());
        if (IRowDiff *row_diff = dynamic_cast<IRowDiff*>(&matrix)) {
            // row_diff side files (no graph needed).
            if (auto *row_diff_column = dynamic_cast<RowDiff<ColumnMajor> *>(&matrix)) {
                row_diff_column->load_anchor(config.infbase + kRowDiffAnchorExt);
                row_diff_column->load_fork_succ(config.infbase + kRowDiffForkSuccExt);
            }
            // Attach the graph (the only annotation step that needs it).
            auto graph = graph_future.get();
            const DeBruijnGraph *base_graph = graph.get();
            if (auto *canonical = dynamic_cast<const CanonicalDBG *>(graph.get()))
                base_graph = &canonical->get_graph();
            row_diff->set_graph(base_graph);
        }
    }

    // Final assembly: await graph (may already be resolved) and wrap PRIMARY.
    auto graph = graph_future.get();
    if (graph->get_mode() == DeBruijnGraph::PRIMARY) {
        graph = std::make_shared<CanonicalDBG>(graph);
        logger->trace("Primary graph wrapped into canonical");
    }
    auto anno_graph = std::make_unique<AnnotatedDBG>(std::move(graph), std::move(annotation),
                                                     false, std::move(coord_to_header));
    if (!anno_graph->check_compatibility()) {
        logger->error("Graph and annotation are not compatible");
        exit(1);
    }
    return anno_graph;
}

} // namespace

std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(std::shared_ptr<DeBruijnGraph> graph,
                                                       const Config &config,
                                                       size_t max_chunks_open) {
    std::promise<std::shared_ptr<DeBruijnGraph>> ready;
    ready.set_value(std::move(graph));
    return build_annotated_dbg(ready.get_future().share(), config, max_chunks_open);
}

std::unique_ptr<AnnotatedDBG> initialize_annotated_dbg(const Config &config) {
    return build_annotated_dbg(async_load_critical_dbg(config), config, kDefaultMaxChunksOpen);
}

std::unique_ptr<AnnotatedDBG>
initialize_annotated_dbg(const Config &config, const PatternPreparation &prep,
                         std::shared_ptr<const annot::CoordToHeader> coord_to_header) {
    return build_annotated_dbg(async_load_critical_dbg(config.infbase, prep), config,
                               kDefaultMaxChunksOpen, std::move(coord_to_header));
}

std::pair<std::shared_future<std::shared_ptr<DeBruijnGraph>>,
          std::future<std::unique_ptr<AnnotatedDBG>>>
load_graph_with_async_annotation(const Config &config) {
    auto graph_future = async_load_critical_dbg(config);
    std::future<std::unique_ptr<AnnotatedDBG>> anno_dbg_future;
    if (config.infbase_annotators.size()) {
        anno_dbg_future = std::async(std::launch::async, [graph_future, config] {
            return build_annotated_dbg(graph_future, config, kDefaultMaxChunksOpen);
        });
    } else {
        // No annotation requested: deliver an immediately-ready null future.
        std::promise<std::unique_ptr<AnnotatedDBG>> ready;
        ready.set_value(nullptr);
        anno_dbg_future = ready.get_future();
    }
    return std::make_pair(std::move(graph_future), std::move(anno_dbg_future));
}


} // namespace cli
} // namespace mtg
