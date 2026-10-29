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
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <shared_mutex>
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

} // namespace

class AnnotatedGraph {
  public:
    using Hits = std::vector<std::pair<std::string, size_t>>;

    AnnotatedGraph(const std::string &graph,
                   const std::string &annotation,
                   unsigned int threads) {
        std::unique_lock<std::shared_mutex> lock(num_threads_mutex);
        auto config = make_config(graph, annotation, threads);
        graph_ = mtg::cli::load_critical_dbg(config.infbase);
        anno_dbg_ = mtg::cli::initialize_annotated_dbg(graph_, config);
        thread_pool_ = std::make_unique<ThreadPool>(std::max(1u, threads));
    }

    size_t k() const { return graph_->get_k(); }

    std::vector<std::string> labels() const {
        return anno_dbg_->get_annotator().get_label_encoder().get_labels();
    }

    // Top labels of each sequence with their numbers of matching k-mers, by
    // decreasing count. A sequence without hits (or shorter than k) gets an
    // empty list.
    std::vector<Hits> query(const std::vector<std::string> &seqs,
                            size_t num_top_labels,
                            double discovery_fraction,
                            double presence_fraction) const {
        if (!std::isfinite(discovery_fraction) || !std::isfinite(presence_fraction)
                || discovery_fraction < 0 || discovery_fraction > 1
                || presence_fraction < 0 || presence_fraction > 1) {
            throw std::invalid_argument(
                "discovery_fraction and presence_fraction must be in [0, 1]");
        }
        std::shared_lock<std::shared_mutex> lock(num_threads_mutex);
        std::vector<Hits> results(seqs.size());
        std::vector<size_t> indices;
        for (size_t i = 0; i < seqs.size(); ++i) {
            if (seqs[i].size() >= k())
                indices.push_back(i);
        }
        run_in_pool(*thread_pool_, indices, [&](size_t i) {
            results[i] = anno_dbg_->get_top_labels(seqs[i], num_top_labels,
                                                   discovery_fraction,
                                                   presence_fraction);
            // get_top_labels sorts only when it truncates to num_top_labels
            std::ranges::stable_sort(results[i], std::greater<>(), &Hits::value_type::second);
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
