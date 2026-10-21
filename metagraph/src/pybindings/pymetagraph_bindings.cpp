// In-process Python bindings for querying an annotated metagraph index.
//
// Exposes AnnotatedDBG::get_top_labels -- the primitive behind `metagraph query`
// and the server's /search endpoint -- for batches of sequences, without
// running a server. Results are plain (label, k-mer count) pairs, so the
// bindings carry no assumptions about how callers turn them into decisions.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <optional>
#include <future>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "cli/load/load_graph.hpp"
#include "common/threads/threading.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/alignment/alignment.hpp"
#include "graph/alignment_redone/aln_query.hpp"
#include "graph/alignment_redone/aln_seeder.hpp"

namespace py = pybind11;

namespace {

// Loading an annotation goes through the CLI's Config (it carries the
// annotation type and loading options), so build a minimal `server_query`
// command line for it. "server_query" is the CLI identity that accepts -i/-a
// without positional input files. Do not pass -p here: Config would mutate the
// process-global thread count. Each Python object owns its query/align pool.
mtg::cli::Config make_config(const std::string &graph,
                             const std::string &annotation) {
    std::vector<std::string> args = {
        "pymetagraph", "server_query",
        "-i", graph,
        "-a", annotation,
    };
    std::vector<char*> argv;
    argv.reserve(args.size());
    for (auto &arg : args) {
        argv.push_back(arg.data());
    }
    return mtg::cli::Config(static_cast<int>(argv.size()), argv.data());
}

// ThreadPool deliberately rethrows task failures on the worker thread. Capture
// them inside each task so Python receives an exception after every task drains.
template <typename Callback>
void run_batch(ThreadPool &pool, size_t size, const Callback &callback) {
    std::vector<std::exception_ptr> failures(size);
    std::vector<std::shared_future<void>> futures;
    futures.reserve(size);
    try {
        for (size_t i = 0; i < size; ++i) {
            futures.push_back(pool.enqueue([&, i]() {
                try { callback(i); }
                catch (...) { failures[i] = std::current_exception(); }
            }));
        }
    } catch (...) {
        // Even an enqueue failure must not leave workers using stack storage.
        for (const auto &future : futures) future.wait();
        throw;
    }
    for (const auto &future : futures) future.get();
    for (const auto &failure : failures) {
        if (failure) std::rethrow_exception(failure);
    }
}

} // namespace

class AnnotatedGraph {
  public:
    using Hits = std::vector<std::pair<std::string, size_t>>;

    AnnotatedGraph(const std::string &graph,
                   const std::string &annotation,
                   unsigned int threads) {
        if (!threads)
            throw std::invalid_argument("threads must be at least 1");
        graph_ = mtg::cli::load_critical_dbg(graph);
        if (!annotation.empty()) {
            auto config = make_config(graph, annotation);
            anno_dbg_ = mtg::cli::initialize_annotated_dbg(graph_, config);
        }
        if (graph_->get_mode() == mtg::graph::DeBruijnGraph::PRIMARY)
            graph_ = std::make_shared<mtg::graph::CanonicalDBG>(graph_);
        thread_pool_ = std::make_unique<ThreadPool>(threads);
    }

    size_t k() const { return graph_->get_k(); }

    std::vector<std::string> labels() const {
        return anno_dbg_ ? anno_dbg_->get_annotator().get_label_encoder().get_labels()
                         : std::vector<std::string>();
    }

    // Top labels of each sequence with their numbers of matching k-mers.
    // A sequence without hits (or shorter than k) gets an empty list.
    std::vector<Hits> query(const std::vector<std::string> &seqs,
                            size_t num_top_labels,
                            double discovery_fraction,
                            double presence_fraction) const {
        if (!anno_dbg_)
            throw std::invalid_argument("query requires an annotation");
        if (!std::isfinite(discovery_fraction) || !std::isfinite(presence_fraction)
                || discovery_fraction < 0 || discovery_fraction > 1
                || presence_fraction < 0 || presence_fraction > 1) {
            throw std::invalid_argument(
                "discovery_fraction and presence_fraction must be in [0, 1]");
        }
        std::vector<Hits> results(seqs.size());
        run_batch(*thread_pool_, seqs.size(), [&](size_t i) {
            if (seqs[i].size() >= k())
                results[i] = anno_dbg_->get_top_labels(seqs[i], num_top_labels,
                                                       discovery_fraction,
                                                       presence_fraction);
        });
        return results;
    }

    std::vector<std::vector<std::string>> align(
            const std::vector<std::string> &seqs, size_t seed_length,
            size_t max_alternative_alignments, double min_exact_match,
            int32_t xdrop, double rel_score_cutoff, bool only_forwards,
            bool connect_anchors, bool extend_chains, size_t max_num_seeds_per_locus,
            double max_num_nodes_per_seq_char, double max_ram_per_alignment) const {
        namespace ar = mtg::graph::align_redone;
        if (!max_alternative_alignments || xdrop < 0
                || !std::isfinite(min_exact_match) || min_exact_match < 0 || min_exact_match > 1
                || !std::isfinite(rel_score_cutoff) || rel_score_cutoff < 0 || rel_score_cutoff > 1
                || !std::isfinite(max_num_nodes_per_seq_char) || max_num_nodes_per_seq_char < 0
                || !std::isfinite(max_ram_per_alignment) || max_ram_per_alignment < 0)
            throw std::invalid_argument("invalid alignment limits or score thresholds");

        ar::DBGAlignerConfig config{};
        config.min_seed_length = std::min(seed_length ? seed_length : size_t(19), k());
        config.max_seed_length = std::numeric_limits<size_t>::max();
        config.num_alternative_paths = max_alternative_alignments;
        config.min_exact_match = min_exact_match;
        config.max_num_seeds_per_locus = max_num_seeds_per_locus;
        config.max_nodes_per_seq_char = max_num_nodes_per_seq_char;
        config.max_ram_per_alignment = max_ram_per_alignment;
        config.xdrop = xdrop;
        config.rel_score_cutoff = rel_score_cutoff;
        config.forward_and_reverse_complement = !only_forwards;
        config.gap_opening_penalty = -6;
        config.gap_extension_penalty = -2;
        config.left_end_bonus = config.right_end_bonus = 5;
        config.alignment_match_score = 2;
        config.alignment_mm_transition_score = config.alignment_mm_transversion_score = 3;
        config.set_scoring_matrix();

        std::vector<std::vector<std::string>> results(seqs.size());
        run_batch(*thread_pool_, seqs.size(), [&](size_t i) {
            ar::Query query(*graph_, seqs[i]);
            std::optional<ar::AnnotationBuffer> annotations;
            std::unique_ptr<ar::ExactSeeder> seeder;
            std::unique_ptr<ar::Extender> extender;
            if (anno_dbg_) {
                annotations.emplace(*graph_, anno_dbg_->get_annotator());
                seeder = std::make_unique<ar::LabeledSeeder>(*annotations, query, config);
                extender = std::make_unique<ar::LabeledExtender>(*annotations, query, config);
            } else {
                seeder = std::make_unique<ar::ExactSeeder>(query, config);
                extender = std::make_unique<ar::Extender>(query, config);
            }
            auto callback = [&](ar::Alignment&& path) {
                if (!connect_anchors
                        && (path.get_spelling().find('$') != std::string_view::npos
                            || std::find(path.get_path().begin(), path.get_path().end(),
                                         mtg::graph::DeBruijnGraph::npos) != path.get_path().end()))
                    throw std::invalid_argument(
                        "Unconnected chains cannot be serialized; use connect_anchors=True");
                mtg::graph::align::Alignment displayed(
                    path.get_seed(), std::vector(path.get_path()), std::string(path.get_spelling()),
                    ar::score_match(path, config),
                    mtg::graph::align::Cigar(path.get_cigar().to_string()), 0, path.get_orientation());
                const auto *encoder = anno_dbg_ ? &anno_dbg_->get_annotator().get_label_encoder() : nullptr;
                if (annotations) {
                    const auto &cols = annotations->get_cached_column_set(path.get_label_classes()[0]);
                    displayed.label_columns.assign(cols.begin(), cols.end());
                    displayed.label_coordinates = path.get_coords();
                }
                auto result = displayed.to_json(k(), !results[i].empty(), {}, {}, encoder,
                    anno_dbg_ ? anno_dbg_->get_coord_to_header() : nullptr, false);
                result["sequence"] = seqs[i];
                result["strand"] = path.get_orientation() ? "-" : "+";
                result["query_start"] = Json::UInt64(path.get_orientation()
                    ? path.get_end_clipping() : path.get_clipping());
                result["query_end"] = Json::UInt64(seqs[i].size() - (path.get_orientation()
                    ? path.get_clipping() : path.get_end_clipping()));
                Json::StreamWriterBuilder writer;
                writer["indentation"] = "";
                results[i].push_back(Json::writeString(writer, result));
            };
            if (extend_chains)
                ar::align_query(query, *seeder, *extender, callback, connect_anchors);
            else
                ar::align_query(query, *seeder, callback, connect_anchors);
        });
        return results;
    }

  private:
    std::shared_ptr<mtg::graph::DeBruijnGraph> graph_;
    std::unique_ptr<mtg::graph::AnnotatedDBG> anno_dbg_;
    std::unique_ptr<ThreadPool> thread_pool_;
};

PYBIND11_MODULE(_pymetagraph_core, m) {
    m.doc() = "In-process queries against an annotated metagraph index";

    py::class_<AnnotatedGraph>(m, "AnnotatedGraph")
        .def(py::init<const std::string&, const std::string&, unsigned int>(),
             py::arg("graph"), py::arg("annotation"), py::arg("threads") = 1,
             py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("k", &AnnotatedGraph::k)
        .def_property_readonly("labels", &AnnotatedGraph::labels)
        .def("query", &AnnotatedGraph::query,
             py::arg("seqs"), py::arg("num_top_labels") = 1,
             py::arg("discovery_fraction") = 0.0, py::arg("presence_fraction") = 0.0,
             py::call_guard<py::gil_scoped_release>())
        .def("align", &AnnotatedGraph::align,
             py::arg("seqs"), py::arg("seed_length") = 0,
             py::arg("max_alternative_alignments") = 1, py::arg("min_exact_match") = 0.7,
             py::arg("xdrop") = 27,
             py::arg("rel_score_cutoff") = 0.95,
             py::arg("only_forwards") = false, py::arg("connect_anchors") = true,
             py::arg("extend_chains") = true,
             py::arg("max_num_seeds_per_locus") = 1000,
             py::arg("max_num_nodes_per_seq_char") = 5.0,
             py::arg("max_ram_per_alignment") = 200.0,
             py::call_guard<py::gil_scoped_release>());
}
