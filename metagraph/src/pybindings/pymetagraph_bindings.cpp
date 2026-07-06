#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstdint>
#include <future>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "cli/align.hpp"
#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "cli/load/load_graph.hpp"
#include "common/threads/threading.hpp"
#include "graph/alignment_redone/aln_match.hpp"
#include "graph/alignment_redone/aln_query.hpp"
#include "graph/alignment_redone/aln_seeder.hpp"
#include "graph/alignment_redone/annotation_buffer.hpp"
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
//
// "method" is consumed here rather than forwarded to Config: it selects which
// pymetagraph binding code path runs (label-presence vs. real alignment), not
// a CLI flag. "connect_anchors"/"extend_chains" are consumed specially too,
// since their CLI flags are bare negations (`--align-no-...`) with no value
// syntax, unlike the rest of kFlagAliases.
std::vector<std::string> kwargs_to_argv(const py::kwargs &kwargs, std::string *method) {
    static const std::unordered_map<std::string, std::string> kFlagAliases = {
        { "input", "-i" },
        { "annotator", "-a" },
        { "threads", "-p" },
        { "num_top_labels", "--num-top-labels" },
        { "discovery_fraction", "--min-kmers-fraction-label" },
        { "presence_fraction", "--min-kmers-fraction-graph" },
        { "seed_length", "--align-min-seed-length" },
        { "max_alternative_alignments", "--align-alternative-alignments" },
        { "max_num_nodes_per_seq_char", "--align-max-nodes-per-seq-char" },
        { "min_exact_match", "--align-min-exact-match" },
        { "max_num_seeds_per_locus", "--align-max-num-seeds-per-locus" },
    };
    static const std::unordered_map<std::string, std::string> kNegationFlagAliases = {
        { "connect_anchors", "--align-no-connect-anchors" },
        { "extend_chains", "--align-no-extend-chains" },
    };

    std::vector<std::string> args = { "pymetagraph", "server_query" };
    *method = "query";

    for (auto item : kwargs) {
        auto key = item.first.cast<std::string>();

        if (key == "method") {
            *method = py::str(item.second).cast<std::string>();
            continue;
        }

        auto neg_it = kNegationFlagAliases.find(key);
        if (neg_it != kNegationFlagAliases.end()) {
            if (!item.second.cast<bool>())
                args.push_back(neg_it->second);
            continue;
        }

        auto it = kFlagAliases.find(key);
        if (it == kFlagAliases.end()) {
            throw std::invalid_argument("pymetagraph: unknown keyword argument '" + key + "'");
        }
        args.push_back(it->second);
        args.push_back(py::str(item.second).cast<std::string>());
    }

    if (*method != "query" && *method != "align") {
        throw std::invalid_argument(
            "pymetagraph: unknown method '" + *method + "' (expected 'query' or 'align')");
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
// queries against it, using one of two unrelated underlying methods picked
// at construction time via the "method" kwarg:
//  - "query" (default): AnnotatedDBG::get_top_labels, the same label-presence
//    primitive backing the CLI's `query` command. Read-only, no mutable
//    per-call state, safe to call concurrently from multiple threads.
//  - "align": the alignment_redone seed-and-extend pipeline, the same one
//    backing the CLI's `align` command and the server's `/align` endpoint
//    (see cli/align.cpp::align_to_graph, cli/server.cpp::process_align_request).
//    Reports real r_st/r_en/strand instead of dummy placeholders. Each call
//    builds its own AnnotationBuffer, since that type caches mutable
//    per-query state and is not meant to be shared across concurrent queries
//    (mirrored from how align_to_graph gives each worker its own buffer).
class Index {
  public:
    explicit Index(const py::kwargs &kwargs) {
        std::string method_str;
        auto arg_strings = kwargs_to_argv(kwargs, &method_str);
        method_ = (method_str == "align") ? Method::kAlign : Method::kQuery;

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

        if (method_ == Method::kAlign) {
            connect_anchors_ = config.alignment_connect_anchors;
            extend_chains_ = config.alignment_extend_chains;
            aligner_config_ = mtg::cli::initialize_aligner_config(config, anno_dbg_->get_graph());
        }

        thread_pool_ = std::make_unique<ThreadPool>(std::max(1u, get_num_threads()));
    }

    Alignment query(const std::string &seq) const {
        if (seq.empty())
            return Alignment("*", 0, 0, 1);

        return method_ == Method::kAlign ? query_by_alignment(seq) : query_by_label(seq);
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
    enum class Method { kQuery, kAlign };

    Alignment query_by_label(const std::string &seq) const {
        auto top_labels = anno_dbg_->get_top_labels(seq, num_top_labels_,
                                                     discovery_fraction_,
                                                     presence_fraction_);
        if (top_labels.empty())
            return Alignment("*", 0, 0, 1);

        return Alignment(top_labels.front().first, 0, static_cast<int64_t>(seq.size()), 1);
    }

    Alignment query_by_alignment(const std::string &seq) const {
        namespace align_redone = mtg::graph::align_redone;

        const auto &graph = anno_dbg_->get_graph();
        align_redone::Query aln_query(graph, seq);
        align_redone::AnnotationBuffer anno_buffer(graph, anno_dbg_->get_annotator());
        align_redone::LabeledSeeder seeder(anno_buffer, aln_query, *aligner_config_);

        std::vector<align_redone::Alignment> paths;
        auto aln_callback = [&](align_redone::Alignment &&aln) {
            paths.emplace_back(std::move(aln));
        };

        if (extend_chains_) {
            align_redone::LabeledExtender extender(anno_buffer, aln_query, *aligner_config_);
            align_redone::align_query(aln_query, seeder, extender, aln_callback, connect_anchors_);
        } else {
            align_redone::align_query(aln_query, seeder, aln_callback, connect_anchors_);
        }

        if (paths.empty())
            return Alignment("*", 0, 0, 1);

        const auto &best = *std::max_element(paths.begin(), paths.end(),
            [&](const align_redone::Alignment &a, const align_redone::Alignment &b) {
                return align_redone::score_match(a, *aligner_config_)
                     < align_redone::score_match(b, *aligner_config_);
            });

        std::string label = "*";
        if (!best.get_label_classes().empty()) {
            auto label_class = best.get_label_classes().front();
            if (label_class != align_redone::Anchor::nannot) {
                const auto &columns = anno_buffer.get_cached_column_set(label_class);
                if (!columns.empty())
                    label = anno_buffer.get_annotator().get_label_encoder().decode(columns.front());
            }
        }

        // Query coordinates are always reported on the original (forward) read,
        // mirroring the reverse-complement clip positions back when needed.
        auto qlen = static_cast<int64_t>(seq.size());
        int64_t r_st, r_en;
        if (!best.get_orientation()) {
            r_st = static_cast<int64_t>(best.get_clipping());
            r_en = qlen - static_cast<int64_t>(best.get_end_clipping());
        } else {
            r_st = static_cast<int64_t>(best.get_end_clipping());
            r_en = qlen - static_cast<int64_t>(best.get_clipping());
        }
        int strand = best.get_orientation() ? -1 : 1;

        return Alignment(label, r_st, r_en, strand);
    }

    Method method_ = Method::kQuery;
    std::shared_ptr<mtg::graph::DeBruijnGraph> graph_;
    std::unique_ptr<mtg::graph::AnnotatedDBG> anno_dbg_;
    std::unique_ptr<ThreadPool> thread_pool_;
    size_t num_top_labels_ = 1;
    double discovery_fraction_ = 0.7;
    double presence_fraction_ = 0.0;

    // Only populated when method_ == Method::kAlign.
    std::optional<mtg::graph::align_redone::DBGAlignerConfig> aligner_config_;
    bool connect_anchors_ = true;
    bool extend_chains_ = true;
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
