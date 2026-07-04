#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <future>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "cli/load/load_graph.hpp"
#include "common/threads/threading.hpp"
#include "graph/annotated_dbg.hpp"

namespace py = pybind11;

namespace {

// Converts Python kwargs into a fake argv so metagraph's existing CLI argument
// parser (mtg::cli::Config) can be reused unchanged instead of reimplementing
// config parsing/validation here. argv[0] is a dummy program name and argv[1]
// is fixed to "server_query" -- the one identity whose validation rules allow
// loading a graph/annotation from -i/-a with no positional FASTA files, since
// it's designed for exactly this "load once, serve many queries" use case
// (plain "query" always demands input files, even reading from stdin otherwise).
std::vector<std::string> kwargs_to_argv(const py::kwargs &kwargs) {
    static const std::unordered_map<std::string, std::string> kFlagAliases = {
        { "input", "-i" },
        { "annotator", "-a" },
        { "threads", "-p" },
        { "num_top_labels", "--num-top-labels" },
        { "discovery_fraction", "--min-kmers-fraction-label" },
        { "presence_fraction", "--min-kmers-fraction-graph" },
    };

    std::vector<std::string> args = { "pymetagraph", "server_query" };

    for (auto item : kwargs) {
        auto key = item.first.cast<std::string>();
        auto it = kFlagAliases.find(key);
        if (it == kFlagAliases.end()) {
            throw std::invalid_argument("pymetagraph: unknown keyword argument '" + key + "'");
        }
        args.push_back(it->second);
        args.push_back(py::str(item.second).cast<std::string>());
    }

    return args;
}

} // namespace

// An alignment-like result. Only the fields Readfish's decision logic reads
// (`ctg`, `r_st`, `r_en`, `strand`) are populated; r_st/r_en/strand are dummy
// placeholders since pymetagraph reports label presence, not real coordinates.
struct Alignment {
    std::string ctg;
    int64_t r_st = 0;
    int64_t r_en = 0;
    int strand = 1;

    Alignment() = default;
    Alignment(std::string ctg, int64_t r_st, int64_t r_en, int strand)
        : ctg(std::move(ctg)), r_st(r_st), r_en(r_en), strand(strand) {}
};

// Loads a prebuilt metagraph graph + single annotation file and answers
// label-presence queries for it. Read-only after construction and safe to
// query concurrently from multiple threads (AnnotatedDBG::get_top_labels has
// no mutable per-call state), unlike metagraph's aligner classes.
class Index {
  public:
    explicit Index(const py::kwargs &kwargs) {
        auto arg_strings = kwargs_to_argv(kwargs);
        std::vector<char*> argv;
        argv.reserve(arg_strings.size());
        for (auto &s : arg_strings)
            argv.push_back(s.data());

        mtg::cli::Config config(static_cast<int>(argv.size()), argv.data());

        if (config.infbase.empty())
            throw std::invalid_argument("pymetagraph: missing required argument 'input'");
        if (config.infbase_annotators.size() != 1) {
            throw std::invalid_argument(
                "pymetagraph: exactly one 'annotator' file is required");
        }

        graph_ = mtg::cli::load_critical_dbg(config.infbase);
        anno_dbg_ = mtg::cli::initialize_annotated_dbg(graph_, config);

        num_top_labels_ = config.num_top_labels;
        discovery_fraction_ = config.discovery_fraction;
        presence_fraction_ = config.presence_fraction;

        thread_pool_ = std::make_unique<ThreadPool>(std::max(1u, get_num_threads()));
    }

    Alignment query(const std::string &seq) const {
        if (seq.empty())
            return Alignment("*", 0, 0, 1);

        auto top_labels = anno_dbg_->get_top_labels(seq, num_top_labels_,
                                                     discovery_fraction_,
                                                     presence_fraction_);
        if (top_labels.empty())
            return Alignment("*", 0, 0, 1);

        return Alignment(top_labels.front().first, 0, static_cast<int64_t>(seq.size()), 1);
    }

    std::vector<Alignment> query_batch(const std::vector<std::string> &seqs) const {
        std::vector<Alignment> results(seqs.size());
        std::vector<std::shared_future<void>> futures;
        futures.reserve(seqs.size());

        for (size_t i = 0; i < seqs.size(); ++i) {
            futures.push_back(thread_pool_->enqueue([&, i]() {
                results[i] = query(seqs[i]);
            }));
        }
        for (auto &f : futures)
            f.get();

        return results;
    }

  private:
    std::shared_ptr<mtg::graph::DeBruijnGraph> graph_;
    std::unique_ptr<mtg::graph::AnnotatedDBG> anno_dbg_;
    std::unique_ptr<ThreadPool> thread_pool_;
    size_t num_top_labels_ = 1;
    double discovery_fraction_ = 0.7;
    double presence_fraction_ = 0.0;
};

PYBIND11_MODULE(_pymetagraph_core, m) {
    py::class_<Alignment>(m, "Alignment")
        .def(py::init<>())
        .def(py::init<std::string, int64_t, int64_t, int>(),
             py::arg("ctg"), py::arg("r_st"), py::arg("r_en"), py::arg("strand"))
        .def_readonly("ctg", &Alignment::ctg)
        .def_readonly("r_st", &Alignment::r_st)
        .def_readonly("r_en", &Alignment::r_en)
        .def_readonly("strand", &Alignment::strand);

    py::class_<Index>(m, "Index")
        .def(py::init<const py::kwargs&>())
        .def("query", &Index::query, py::arg("seq"),
             py::call_guard<py::gil_scoped_release>())
        .def("query_batch", &Index::query_batch, py::arg("seqs"),
             py::call_guard<py::gil_scoped_release>());
}
