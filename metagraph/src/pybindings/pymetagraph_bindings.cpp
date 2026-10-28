// In-process Python bindings for querying an annotated metagraph index.
//
// Exposes AnnotatedDBG::get_top_labels -- the primitive behind `metagraph query`
// and the server's /search endpoint -- for batches of sequences, without
// running a server. Results are plain (label, k-mer count) pairs, so the
// bindings carry no assumptions about how callers turn them into decisions.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <exception>
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

} // namespace

class AnnotatedGraph {
  public:
    using Hits = std::vector<std::pair<std::string, size_t>>;

    AnnotatedGraph(const std::string &graph,
                   const std::string &annotation,
                   unsigned int threads) {
        auto config = make_config(graph, annotation, threads);
        graph_ = mtg::cli::load_critical_dbg(config.infbase);
        anno_dbg_ = mtg::cli::initialize_annotated_dbg(graph_, config);
        thread_pool_ = std::make_unique<ThreadPool>(std::max(1u, get_num_threads()));
    }

    size_t k() const { return graph_->get_k(); }

    std::vector<std::string> labels() const {
        return anno_dbg_->get_annotator().get_label_encoder().get_labels();
    }

    // Top labels of each sequence with their numbers of matching k-mers.
    // A sequence without hits (or shorter than k) gets an empty list.
    std::vector<Hits> query(const std::vector<std::string> &seqs,
                            size_t num_top_labels,
                            double discovery_fraction,
                            double presence_fraction) const {
        if (discovery_fraction < 0 || discovery_fraction > 1
                || presence_fraction < 0 || presence_fraction > 1) {
            throw std::invalid_argument(
                "discovery_fraction and presence_fraction must be in [0, 1]");
        }
        std::vector<Hits> results(seqs.size());
        std::vector<std::shared_future<void>> futures;
        futures.reserve(seqs.size());
        for (size_t i = 0; i < seqs.size(); ++i) {
            if (seqs[i].size() < k())
                continue;
            futures.push_back(thread_pool_->enqueue([&, i]() {
                results[i] = anno_dbg_->get_top_labels(seqs[i], num_top_labels,
                                                       discovery_fraction,
                                                       presence_fraction);
            }));
        }
        // Wait for every task before rethrowing, as they write into results
        std::exception_ptr error;
        for (auto &future : futures) {
            try {
                future.get();
            } catch (...) {
                if (!error)
                    error = std::current_exception();
            }
        }
        if (error)
            std::rethrow_exception(error);

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
             py::call_guard<py::gil_scoped_release>());
}
