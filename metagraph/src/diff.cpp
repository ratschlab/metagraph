#include <cmath>
#include <random>

#include <gflags/gflags.h>
#include <progress_bar.hpp>
#include <sdsl/enc_vector.hpp>

#include "annotation//binary_matrix/multi_brwt/brwt.hpp"
#include "annotation//binary_matrix/row_flat/flat_matrix.hpp"
#include "cli/load/load_annotation.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "graph/annotated_dbg.hpp"
#include "seq_io/sequence_io.hpp"

DEFINE_string(format, "diff", "Format to compress in: diff or zip");
DEFINE_string(o, "/tmp/graph_750" + FLAGS_format, "Output file");
DEFINE_string(graph, "/tmp/diffs/graph_subset_750_primary.dbg", "Input graph file");
DEFINE_string(annotation, "/tmp/graph_subset_750_rowflat.flat.annodbg", "Input annotation file");

// Note: if testing with many chunks use 'ulimit -n <max_files>' to increase the maxium
// number of files the system allows you to open. On mac the default is only 256!

using namespace mtg;

std::shared_ptr<DeBruijnGraph> build_graph(const std::string &filename) {
    std::vector<std::string> sequences;
    seq_io::read_fasta_file_critical(
            filename,
            [&](seq_io::kseq_t *stream) { sequences.emplace_back(stream->seq.s); }, true);

    size_t k = 12;

    BOSSConstructor constructor(k - 1);
    constructor.add_sequences(std::move(sequences));
    auto graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor));
    dynamic_cast<DBGSuccinct *>(graph.get())->mask_dummy_kmers(1, false);

    return graph;
}

std::unique_ptr<annotate::MultiLabelEncoded<std::string>>
load_annotations(const std::vector<std::string> &anno_files) {
    std::unique_ptr<annotate::MultiLabelEncoded<std::string>> annotation
            = cli::initialize_annotation(anno_files[0]);

    if (!annotation->merge_load(anno_files)) {
        std::cerr << "Can't load annotation from " << std::getenv("ANNO") << std::endl;
        exit(1);
    }

    return annotation;
}

std::unique_ptr<annotate::MultiLabelEncoded<std::string>>
load_annotation(const std::string &anno_file) {
    std::unique_ptr<annotate::MultiLabelEncoded<std::string>> annotation
            = cli::initialize_annotation(anno_file);

    if (!annotation->load(anno_file)) {
        std::cerr << "Can't load annotation from " << std::getenv("ANNO") << std::endl;
        exit(1);
    }
    std::cout << "Loaded annotation " << anno_file << std::endl;
    return annotation;
}

std::shared_ptr<DBGSuccinct> load_graph(const std::string &name) {
    std::shared_ptr<DBGSuccinct> result = std::make_shared<DBGSuccinct>(20);
    if (!result->load(name)) {
        std::cerr << "Cannot load graph " << name << std::endl;
        std::exit(1);
    }
    std::cout << "Loaded graph with k-mer length " << result->get_k() << ", type "
              << result->get_state() << " number of nodes: " << result->num_nodes()
              << std::endl;
    return result;
}

void store_diff(const AnnotatedDBG& annotated_dbg) {
    std::vector<uint64_t> diff;
    std::vector<bool> boundary;
    std::vector<bool> full;
    uint64_t max_id = 0;
    const DeBruijnGraph &graph = annotated_dbg.get_graph();
    ProgressBar progress_bar(graph.num_nodes(), "Building Diff-anno");
    graph.call_sequences([&](const std::string &, const std::vector<uint64_t> &path) {
      std::vector<uint64_t> anno_ids;

      for (const uint64_t node_id : path) {
          anno_ids.push_back(annotated_dbg.graph_to_anno_index(node_id));
      }
      std::vector<Vector<uint64_t>> rows
              = annotated_dbg.get_annotation().get_matrix().get_rows(anno_ids);

      for (uint32_t i = 0; i < rows.size() - 1; ++i) {
          uint64_t idx1 = 0;
          uint64_t idx2 = 0;
          while (idx1 < rows[i].size() && idx2 < rows[i + 1].size()) {
              if (rows[i][idx1] == rows[i + 1][idx2]) {
                  idx1++;
                  idx2++;
              } else {
                  if (rows[i][idx1] < rows[i + 1][idx2]) {
                      diff.push_back(rows[i][idx1]);
                      boundary.push_back(0);
                      idx1++;
                  } else {
                      diff.push_back(rows[i + 1][idx2]);
                      boundary.push_back(0);
                      idx2++;
                  }
                  if (max_id < diff.back()) {
                      max_id = diff.back();
                  }
              }
          }
          while (idx1 < rows[i].size()) {
              diff.push_back(rows[i][idx1]);
              boundary.push_back(0);
              idx1++;
          }
          while (idx2 < rows[i + 1].size()) {
              diff.push_back(rows[i + 1][idx2]);
              boundary.push_back(0);
              idx2++;
          }
          boundary.push_back(1);
          full.push_back(0);
      }

      for (const auto v : rows[rows.size() - 1]) {
          diff.push_back(v);
          boundary.push_back(0);
      }
      full.push_back(1);
      boundary.push_back(1);

      progress_bar += path.size();
    });

    std::ofstream f(FLAGS_o);
    if (max_id > 0) {
        uint8_t width = std::ceil(std::log2(max_id));
        std::cout << "Max id is " << max_id << " storing in " << (int)width << " bits."
                  << std::endl;
        sdsl::int_vector<64> sdiff(diff.size(), 0, width);
        for (uint64_t i = 0; i < diff.size(); ++i) {
            sdiff[i] = diff[i];
        }
        sdsl::enc_vector<> cdiff(sdiff);
        cdiff.serialize(f);
    }
    sdsl::bit_vector sboundary(boundary.size());
    for (uint64_t i = 0; i < diff.size(); ++i) {
        sboundary[i] = boundary[i];
    }
    sdsl::sd_vector<> cboundary(sboundary);
    cboundary.serialize(f);

    sdsl::bit_vector sfull(full.size());
    for (uint64_t i = 0; i < diff.size(); ++i) {
        sfull[i] = full[i];
    }
    sdsl::sd_vector<> cfull(sfull);
    cfull.serialize(f);

    f.close();
}

void store_zip(const AnnotatedDBG& annotated_dbg) {
    std::ofstream out(FLAGS_o);
    const DeBruijnGraph &graph = annotated_dbg.get_graph();
    ProgressBar progress_bar(graph.num_nodes(), "Building Diff-anno");
    graph.call_sequences([&](const std::string &, const std::vector<uint64_t> &path) {
      std::vector<uint64_t> anno_ids;

      for (const uint64_t node_id : path) {
          anno_ids.push_back(annotated_dbg.graph_to_anno_index(node_id));
      }
      std::vector<Vector<uint64_t>> rows
              = annotated_dbg.get_annotation().get_matrix().get_rows(anno_ids);

      for (uint32_t i = 0; i < rows.size(); ++i) {
          sdsl::bit_vector row(annotated_dbg.get_annotation().get_matrix().num_columns());
          for (const auto v : rows[i]) {
              row[v] = 1;
          }
          row.serialize(out);
      }

      progress_bar += path.size();
    });

    out.close();
}

int main(int argc, char* argv[]) {
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    spdlog::set_level(spdlog::level::trace);

    std::unique_ptr<annotate::MultiLabelEncoded<std::string>> annotation
            = load_annotation(FLAGS_annotation);
    std::shared_ptr<DeBruijnGraph> graph = load_graph(FLAGS_graph);
    AnnotatedDBG annotated_dbg(graph, std::move(annotation), false);
    if (FLAGS_format == "diff") {
        store_diff(annotated_dbg);
    } else if (FLAGS_format == "zip") {
        store_zip(annotated_dbg);
    } else {
        std::cerr << "Invalid -format, only 'zip' and 'diff' supported.";
        std::exit(1);
    }
}
