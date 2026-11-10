// In-process Python bindings for querying an annotated metagraph index.
//
// Exposes, for batches of sequences and without running a server:
//  - AnnotatedDBG::get_top_labels, the primitive behind `metagraph query` and
//    the server's /search endpoint,
//  - sequence-to-graph alignment as in `metagraph align` (labeled with the
//    annotation, or unlabeled as in the server's /align endpoint),
//  - the seed-only label mapping of `metagraph align --align-map-only`.
// Results are plain values, so the bindings carry no assumptions about how
// callers turn them into decisions.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <future>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <shared_mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "cli/align.hpp"
#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "cli/load/load_graph.hpp"
#include "common/threads/threading.hpp"
#include "graph/alignment/aligner_labeled.hpp"
#include "graph/alignment/dbg_aligner.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/graph_extensions/node_first_cache.hpp"
#include "graph/representation/canonical_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"

namespace py = pybind11;

namespace {

// Loading an annotation goes through the CLI's Config (it carries the
// annotation type and loading options), so build a minimal `server_query`
// command line for it. "server_query" is the CLI identity that accepts -i/-a
// without positional input files.
mtg::cli::Config make_config(const std::string &graph,
                             const std::string &annotation,
                             unsigned int threads) {
    std::vector<std::string> args = {
        "pymetagraph", "server_query",
        "-i", graph,
        "-a", annotation,
        "-p", std::to_string(std::max(1u, threads)),
    };
    std::vector<char*> argv;
    argv.reserve(args.size());
    for (auto &arg : args) {
        argv.push_back(arg.data());
    }
    return mtg::cli::Config(static_cast<int>(argv.size()), argv.data());
}

// Config sets MetaGraph's process-global (non-atomic) thread count, which the
// loaders and some annotation queries read. Constructions write it exclusively,
// queries read it under a shared lock.
std::shared_mutex num_threads_mutex;

// Run task(i) for every i in indices on the pool and wait for all of them.
// ThreadPool terminates the process if a task throws on a worker, so catch
// inside each task and rethrow the first exception once all have finished.
void run_in_pool(ThreadPool &pool, const std::vector<size_t> &indices,
                 const std::function<void(size_t)> &task) {
    std::mutex error_mutex;
    std::exception_ptr error;
    std::vector<std::shared_future<void>> futures;
    futures.reserve(indices.size());
    for (size_t i : indices) {
        futures.push_back(pool.enqueue([&, i]() {
            try {
                task(i);
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!error)
                    error = std::current_exception();
            }
        }));
    }
    for (auto &future : futures) {
        future.wait();
    }
    if (error)
        std::rethrow_exception(error);
}

// Split [0, n) into at most num_chunks contiguous ranges of similar size.
std::vector<std::pair<size_t, size_t>> split_range(size_t n, size_t num_chunks) {
    std::vector<std::pair<size_t, size_t>> ranges;
    num_chunks = std::max<size_t>(1, std::min(n, num_chunks));
    for (size_t c = 0, begin = 0; c < num_chunks; ++c) {
        size_t end = begin + (n - begin) / (num_chunks - c);
        ranges.emplace_back(begin, end);
        begin = end;
    }
    return ranges;
}

} // namespace

// One alignment of a query to the graph. Query coordinates refer to the query
// in the aligned orientation (its reverse complement if orientation is true).
struct GraphAlignment {
    std::vector<std::string> labels;
    int64_t score;
    int64_t max_score;
    std::string cigar;
    bool orientation;
    size_t query_begin;
    size_t query_end;
};

// The labels whose seeds cover the most query characters on the better strand.
// Query coordinates refer to the forward query.
struct LabelMapping {
    std::vector<std::string> labels;
    size_t num_matches;
    bool orientation;
    size_t query_begin;
    size_t query_end;
};

class AnnotatedGraph {
  public:
    using Hits = std::vector<std::pair<std::string, size_t>>;

    AnnotatedGraph(const std::string &graph,
                   const std::string &annotation,
                   unsigned int threads) {
        std::unique_lock<std::shared_mutex> lock(num_threads_mutex);
        auto config = make_config(graph, annotation, threads);
        // a file that cannot be loaded, or a record mapping (.seqs) that is not the
        // annotation's, raises (RuntimeError) instead of ending the Python process
        graph_ = mtg::cli::load_critical_dbg(config.infbase, mtg::cli::LoadFailure::THROW);
        anno_dbg_ = mtg::cli::initialize_annotated_dbg(graph_, config,
                                                       mtg::cli::kDefaultMaxChunksOpen,
                                                       mtg::cli::LoadFailure::THROW);
        num_threads_ = std::max(1u, threads);
        thread_pool_ = std::make_unique<ThreadPool>(num_threads_);
        config_ = std::make_unique<mtg::cli::Config>(config);

        // as in `metagraph align`: drop the dummy k-mer mask, so that suffix
        // seeds can match dummy nodes (full k-mers never do, so the label query
        // is unaffected)
        if (auto dbg_succ = std::dynamic_pointer_cast<mtg::graph::DBGSuccinct>(graph_))
            dbg_succ->reset_mask();

        // as in `metagraph align` (with a cache shared by all threads): align
        // primary graphs through their canonical wrapper, and cache backward
        // traversals of succinct graphs
        aln_graph_ = graph_;
        if (aln_graph_->get_mode() == mtg::graph::DeBruijnGraph::PRIMARY)
            aln_graph_ = std::make_shared<mtg::graph::CanonicalDBG>(aln_graph_);
        if (auto dbg_succ = std::dynamic_pointer_cast<mtg::graph::DBGSuccinct>(graph_))
            aln_graph_->add_extension(std::make_shared<mtg::graph::NodeFirstCache>(*dbg_succ));
    }

    size_t k() const { return graph_->get_k(); }

    std::vector<std::string> labels() const {
        return anno_dbg_->get_annotator().get_label_encoder().get_labels();
    }

    // Top labels of each sequence with their numbers of matching k-mers, by
    // decreasing count. A sequence without hits (or shorter than k) gets an
    // empty list. No num_top_labels means no limit, as in `metagraph query`.
    std::vector<Hits> query(const std::vector<std::string> &seqs,
                            std::optional<size_t> num_top_labels,
                            double discovery_fraction,
                            double presence_fraction) const {
        if (!std::isfinite(discovery_fraction) || !std::isfinite(presence_fraction)
                || discovery_fraction < 0 || discovery_fraction > 1
                || presence_fraction < 0 || presence_fraction > 1) {
            throw std::invalid_argument(
                "discovery_fraction and presence_fraction must be in [0, 1]");
        }
        const size_t max_labels = num_top_labels.value_or(std::numeric_limits<size_t>::max());
        std::shared_lock<std::shared_mutex> lock(num_threads_mutex);
        std::vector<Hits> results(seqs.size());
        std::vector<size_t> indices;
        for (size_t i = 0; i < seqs.size(); ++i) {
            if (seqs[i].size() >= k())
                indices.push_back(i);
        }
        run_in_pool(*thread_pool_, indices, [&](size_t i) {
            results[i] = anno_dbg_->get_top_labels(seqs[i], max_labels,
                                                   discovery_fraction,
                                                   presence_fraction);
            // get_top_labels sorts only when it truncates to num_top_labels
            std::ranges::stable_sort(results[i], std::greater<>(), &Hits::value_type::second);
        });

        return results;
    }

    // Alignments of each sequence to the graph, best first, as in `metagraph
    // align`. With labeled=true, alignments are restricted to and report the
    // labels of the annotation (`metagraph align -a`); with labeled=false they
    // have no labels, as with the server's /align endpoint. Arguments left
    // unset keep the defaults of `metagraph align`.
    std::vector<std::vector<GraphAlignment>>
    align(const std::vector<std::string> &seqs,
          bool labeled,
          std::optional<double> min_exact_match,
          std::optional<size_t> max_alternative_alignments,
          std::optional<double> max_num_nodes_per_seq_char,
          std::optional<size_t> min_seed_length) const {
        auto aligner_config = make_aligner_config(min_exact_match, max_alternative_alignments,
                                                  max_num_nodes_per_seq_char, min_seed_length);
        const auto &encoder = anno_dbg_->get_annotator().get_label_encoder();
        std::shared_lock<std::shared_mutex> lock(num_threads_mutex);
        std::vector<std::vector<GraphAlignment>> results(seqs.size());
        for_each_batch(seqs, [&](const std::vector<mtg::graph::align::IDBGAligner::Query> &batch) {
            std::unique_ptr<mtg::graph::align::IDBGAligner> aligner;
            if (labeled) {
                aligner = std::make_unique<mtg::graph::align::LabeledAligner<>>(
                        *aln_graph_, aligner_config, anno_dbg_->get_annotator());
            } else {
                aligner = std::make_unique<mtg::graph::align::DBGAligner<>>(
                        *aln_graph_, aligner_config);
            }
            aligner->align_batch(batch,
                [&](const std::string &header, mtg::graph::align::AlignmentResults&& paths) {
                    auto &out = results[std::stoull(header)];
                    for (const auto &path : paths) {
                        GraphAlignment a;
                        for (auto column : path.label_columns) {
                            a.labels.push_back(encoder.decode(column));
                        }
                        a.score = path.get_score();
                        a.max_score = aligner_config.match_score(paths.get_query(path.get_orientation()))
                                + aligner_config.left_end_bonus + aligner_config.right_end_bonus;
                        a.cigar = path.get_cigar().to_string();
                        a.orientation = path.get_orientation();
                        a.query_begin = path.get_clipping();
                        a.query_end = path.get_clipping() + path.get_query_view().size();
                        out.push_back(std::move(a));
                    }
                }
            );
        });
        return results;
    }

    // Seed-only mapping of each sequence to labels, without extending the seeds
    // into alignments (`metagraph align --align-map-only`). Seeds are found and
    // filtered exactly as for labeled alignment. A sequence without a mapping
    // gets None.
    std::vector<std::optional<LabelMapping>>
    map(const std::vector<std::string> &seqs,
        std::optional<double> min_exact_match,
        std::optional<size_t> min_seed_length) const {
        auto aligner_config = make_aligner_config(min_exact_match, std::nullopt, std::nullopt,
                                                  min_seed_length);
        const auto &encoder = anno_dbg_->get_annotator().get_label_encoder();
        std::shared_lock<std::shared_mutex> lock(num_threads_mutex);
        std::vector<std::optional<LabelMapping>> results(seqs.size());
        for_each_batch(seqs, [&](const std::vector<mtg::graph::align::IDBGAligner::Query> &batch) {
            mtg::graph::align::LabeledAligner<> aligner(*aln_graph_, aligner_config,
                                                        anno_dbg_->get_annotator());
            aligner.map_batch(batch,
                [&](const std::string &header,
                    std::optional<mtg::graph::align::LabeledAligner<>::LabelMapping>&& mapping) {
                    if (!mapping)
                        return;
                    LabelMapping m;
                    for (auto column : mapping->labels) {
                        m.labels.push_back(encoder.decode(column));
                    }
                    m.num_matches = mapping->num_matches;
                    m.orientation = mapping->orientation;
                    m.query_begin = mapping->query_begin;
                    m.query_end = mapping->query_end;
                    results[std::stoull(header)] = std::move(m);
                }
            );
        });
        return results;
    }

  private:
    std::shared_ptr<mtg::graph::DeBruijnGraph> graph_;
    std::shared_ptr<mtg::graph::DeBruijnGraph> aln_graph_;
    std::unique_ptr<mtg::graph::AnnotatedDBG> anno_dbg_;
    std::unique_ptr<ThreadPool> thread_pool_;
    size_t num_threads_;
    std::unique_ptr<mtg::cli::Config> config_;

    // The aligner configuration of `metagraph align`, with the given overrides
    // validated as the CLI validates its options
    mtg::graph::align::DBGAlignerConfig
    make_aligner_config(std::optional<double> min_exact_match,
                        std::optional<size_t> max_alternative_alignments,
                        std::optional<double> max_num_nodes_per_seq_char,
                        std::optional<size_t> min_seed_length) const {
        mtg::cli::Config config(*config_);
        if (min_exact_match) {
            if (!std::isfinite(*min_exact_match) || *min_exact_match < 0 || *min_exact_match > 1)
                throw std::invalid_argument("min_exact_match must be in [0, 1]");
            config.alignment_min_exact_match = *min_exact_match;
        }
        if (max_alternative_alignments) {
            if (*max_alternative_alignments == 0)
                throw std::invalid_argument("max_alternative_alignments must be positive");
            config.alignment_num_alternative_paths = *max_alternative_alignments;
        }
        if (max_num_nodes_per_seq_char) {
            if (!(*max_num_nodes_per_seq_char > 0))
                throw std::invalid_argument("max_num_nodes_per_seq_char must be positive");
            config.alignment_max_nodes_per_seq_char = *max_num_nodes_per_seq_char;
        }
        if (min_seed_length) {
            if (*min_seed_length == 0)
                throw std::invalid_argument("min_seed_length must be positive");
            config.alignment_min_seed_length = *min_seed_length;
        }
        return mtg::cli::initialize_aligner_config(config, *aln_graph_);
    }

    // Split the sequences into one batch per thread and run task(batch) on the
    // pool. Queries are named by their index in seqs, which the callbacks use to
    // place results.
    void for_each_batch(const std::vector<std::string> &seqs,
                        const std::function<void(const std::vector<mtg::graph::align::IDBGAligner::Query>&)> &task) const {
        // empty sequences have no alignments and are left out
        std::vector<size_t> nonempty;
        for (size_t i = 0; i < seqs.size(); ++i) {
            if (seqs[i].size())
                nonempty.push_back(i);
        }
        auto ranges = split_range(nonempty.size(), num_threads_);
        std::vector<size_t> indices(ranges.size());
        std::iota(indices.begin(), indices.end(), 0);
        run_in_pool(*thread_pool_, indices, [&](size_t c) {
            auto [begin, end] = ranges[c];
            std::vector<mtg::graph::align::IDBGAligner::Query> batch;
            batch.reserve(end - begin);
            for (size_t i = begin; i < end; ++i) {
                batch.emplace_back(std::to_string(nonempty[i]), seqs[nonempty[i]]);
            }
            task(batch);
        });
    }
};

PYBIND11_MODULE(_pymetagraph_core, m) {
    m.doc() = "In-process queries and alignment against an annotated metagraph index";

    py::class_<GraphAlignment>(m, "GraphAlignment")
        .def_readonly("labels", &GraphAlignment::labels)
        .def_readonly("score", &GraphAlignment::score)
        .def_readonly("max_score", &GraphAlignment::max_score)
        .def_readonly("cigar", &GraphAlignment::cigar)
        .def_readonly("orientation", &GraphAlignment::orientation)
        .def_readonly("query_begin", &GraphAlignment::query_begin)
        .def_readonly("query_end", &GraphAlignment::query_end)
        .def("__repr__", [](const GraphAlignment &a) {
            return "GraphAlignment(labels=" + py::repr(py::cast(a.labels)).cast<std::string>()
                + ", score=" + std::to_string(a.score) + ", cigar='" + a.cigar
                + "', orientation=" + (a.orientation ? "True" : "False")
                + ", query_begin=" + std::to_string(a.query_begin)
                + ", query_end=" + std::to_string(a.query_end) + ")";
        });

    py::class_<LabelMapping>(m, "LabelMapping")
        .def_readonly("labels", &LabelMapping::labels)
        .def_readonly("num_matches", &LabelMapping::num_matches)
        .def_readonly("orientation", &LabelMapping::orientation)
        .def_readonly("query_begin", &LabelMapping::query_begin)
        .def_readonly("query_end", &LabelMapping::query_end)
        .def("__repr__", [](const LabelMapping &a) {
            return "LabelMapping(labels=" + py::repr(py::cast(a.labels)).cast<std::string>()
                + ", num_matches=" + std::to_string(a.num_matches)
                + ", orientation=" + (a.orientation ? "True" : "False")
                + ", query_begin=" + std::to_string(a.query_begin)
                + ", query_end=" + std::to_string(a.query_end) + ")";
        });

    py::class_<AnnotatedGraph>(m, "AnnotatedGraph")
        .def(py::init<const std::string&, const std::string&, unsigned int>(),
             py::arg("graph"), py::arg("annotation"), py::arg("threads") = 1,
             py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("k", &AnnotatedGraph::k)
        .def_property_readonly("labels", &AnnotatedGraph::labels)
        .def("query", &AnnotatedGraph::query,
             // same defaults as `metagraph query`
             py::arg("seqs"), py::arg("num_top_labels") = py::none(),
             py::arg("discovery_fraction") = 0.7, py::arg("presence_fraction") = 0.0,
             py::call_guard<py::gil_scoped_release>())
        .def("align", &AnnotatedGraph::align,
             // unset arguments keep the defaults of `metagraph align`
             py::arg("seqs"), py::arg("labeled") = true,
             py::arg("min_exact_match") = py::none(),
             py::arg("max_alternative_alignments") = py::none(),
             py::arg("max_num_nodes_per_seq_char") = py::none(),
             py::arg("min_seed_length") = py::none(),
             py::call_guard<py::gil_scoped_release>())
        .def("map", &AnnotatedGraph::map,
             py::arg("seqs"), py::arg("min_exact_match") = py::none(),
             py::arg("min_seed_length") = py::none(),
             py::call_guard<py::gil_scoped_release>());

    // Test hook: a task that throws on a worker thread must surface as a Python
    // exception rather than terminate the process
    m.def("_raise_in_worker", [](unsigned int threads) {
        ThreadPool pool(std::max(1u, threads));
        run_in_pool(pool, { 0, 1, 2, 3 }, [](size_t i) {
            if (i == 2)
                throw std::runtime_error("error in worker");
        });
    }, py::arg("threads") = 2, py::call_guard<py::gil_scoped_release>());
}
