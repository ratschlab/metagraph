#include "row_diff_builder.hpp"

#include <progress_bar.hpp>

#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/sorted_sets/sorted_set_disk.hpp"
#include "common/threads/threading.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "graph/annotated_dbg.hpp"

namespace mtg {
namespace annot {

using namespace mtg::annot::binmat;
using mtg::common::logger;

constexpr uint64_t chunk_size = 1 << 20;

void traverse_anno_chunked(
        const std::string &log_header,
        uint64_t num_rows,
        const std::string &pred_succ_fprefix,
        const std::vector<std::unique_ptr<annot::ColumnCompressed<>>> &col_annotations,
        const std::function<void()> &before_chunk,
        const CallOnes &call_ones,
        const std::function<void(node_index start, uint64_t size)> &after_chunk) {
    if (col_annotations.empty())
        return;

    const uint32_t num_threads = get_num_threads();

    sdsl::int_vector_buffer succ(pred_succ_fprefix + ".succ", std::ios::in, 1024 * 1024);
    sdsl::int_vector_buffer pred(pred_succ_fprefix + ".pred", std::ios::in, 1024 * 1024);
    sdsl::int_vector_buffer<1> pred_boundary(pred_succ_fprefix + ".pred_boundary",
                                             std::ios::in, 1024 * 1024);

    assert(succ.size() == num_rows);
    assert(static_cast<uint64_t>(std::count(pred_boundary.begin(), pred_boundary.end(), 0))
                   == pred.size());
    std::vector<uint64_t> succ_chunk;
    std::vector<uint64_t> pred_chunk_idx;

    auto pred_boundary_it = pred_boundary.begin();
    auto pred_it = pred.begin();
    ProgressBar progress_bar(num_rows, log_header, std::cerr, !common::get_verbose());
    for (node_index chunk = 0; chunk < num_rows; chunk += chunk_size) {
        succ_chunk.resize(std::min(chunk_size, num_rows - chunk));
        pred_chunk_idx.resize(std::min(chunk_size, num_rows - chunk) + 1);
        pred_chunk_idx[0] = 0;

        std::vector<uint64_t> pred_chunk;
        pred_chunk.reserve(succ_chunk.size() * 1.1);

        before_chunk();

        for (node_index idx = chunk, i = 0; i < succ_chunk.size(); ++idx, ++i) {
            succ_chunk[i] = succ[idx] == 0
                            ? std::numeric_limits<uint64_t>::max()
                            : graph::AnnotatedSequenceGraph::graph_to_anno_index(succ[idx]);
            pred_chunk_idx[i + 1] = pred_chunk_idx[i];
            while (*pred_boundary_it == 0) {
                ++pred_chunk_idx[i + 1];
                pred_chunk.push_back(
                        graph::AnnotatedSequenceGraph::graph_to_anno_index(*pred_it));
                ++pred_it;
                ++pred_boundary_it;
            }
            ++pred_boundary_it;
        }

        assert(pred_chunk.size() == pred_chunk_idx.back());

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (uint64_t l_idx = 0; l_idx < col_annotations.size(); ++l_idx) {
            if (col_annotations[l_idx]->num_labels()
                && col_annotations[l_idx]->num_objects() != num_rows) {
                logger->error(
                        "Graph and annotation are incompatible. Graph has {} nodes, "
                        "annotation has {} entries",
                        num_rows, col_annotations[l_idx]->num_objects());
                std::exit(1);
            }
            for (uint64_t l_idx2 = 0; l_idx2 < col_annotations[l_idx]->num_labels(); ++l_idx2) {
                const std::unique_ptr<bit_vector> &source_col
                        = col_annotations[l_idx]->get_matrix().data()[l_idx2];
                source_col->call_ones_in_range(
                        chunk, chunk + succ_chunk.size(), [&](node_index row_idx) {
                          call_ones(*source_col, row_idx, row_idx - chunk, l_idx,
                                    l_idx2, succ_chunk, pred_chunk, pred_chunk_idx);
                        });
            }
        }
        after_chunk(chunk, succ_chunk.size());

        progress_bar += succ_chunk.size();
    }
    assert(pred_boundary_it == pred_boundary.end());
}


void convert_batch_to_row_diff(const std::string &graph_fname,
                               const std::vector<std::string> &source_files,
                               const std::filesystem::path &dest_dir) {
    if (source_files.empty())
        return;

    uint32_t num_threads = get_num_threads();

    std::vector<std::unique_ptr<annot::ColumnCompressed<>>> sources;
    for (const auto &fname : source_files) {
        auto anno = std::make_unique<annot::ColumnCompressed<>>() ;
        anno->merge_load({fname});
        sources.push_back(std::move(anno));
    }
    logger->trace("Done loading {} annotations", sources.size());

    sdsl::rrr_vector rterminal;

    // if we just generated anchor nodes, attempt a greedy anchor optimization
    if (!std::filesystem::exists(graph_fname + "terminal")) {
        logger->trace("Performing anchor optimization");
        sdsl::bit_vector terminal;
        std::ifstream f(graph_fname + ".terminal.unopt", std::ios::binary);
        terminal.load(f);
        f.close();

        // total number of set bits in the original and sparsified rows
        std::vector<std::atomic<uint32_t>> orig_ones(chunk_size);
        std::vector<std::atomic<uint32_t>> sparse_ones(chunk_size);

        traverse_anno_chunked(
                "Anchor opt", terminal.size(), graph_fname, sources,
                [&]() {
                  std::fill(orig_ones.begin(), orig_ones.end(), 0);
                  std::fill(sparse_ones.begin(), sparse_ones.end(), 0);
                },
                [&](const bit_vector &source_col, node_index row_idx, node_index chunk_idx,
                    uint64_t, uint64_t, const std::vector<uint64_t> &succ_chunk,
                    const std::vector<uint64_t> &, const std::vector<uint64_t> &) {
                  if (terminal[row_idx])
                      return;

                  // if successor is not set, add it to the diff
                  if (!source_col[succ_chunk[chunk_idx]])
                      sparse_ones[row_idx]++;

                  orig_ones[row_idx]++;
                },
                [&](node_index chunk_start, uint64_t chunk_size) {
                  for (uint64_t i = 0; i < chunk_size; ++i) {
                      if (sparse_ones[i] >= 2 && sparse_ones[i] > orig_ones[i] / 2) {
                          terminal[i + chunk_start] = 1;
                      }
                  }
                });

        // save the optimized terminal bit vector, and delete the unoptimized one
        std::ofstream fterm(graph_fname + ".terminal", std::ios::binary);
        rterminal = sdsl::rrr_vector(terminal);
        rterminal.serialize(fterm);
        fterm.close();
        std::filesystem::remove(graph_fname + ".terminal.unopt");
    } else {
        std::ifstream f(graph_fname + ".terminal", std::ios::binary);
        rterminal.load(f);
        f.close();
    }

    // accumulate the indices for the set bits in each column into a #SortedSetDisk
    using SSD = common::SortedSetDisk<uint64_t>;
    std::vector<std::vector<std::unique_ptr<SSD>>> targets(sources.size());
    std::vector<std::vector<uint64_t>> targets_size(sources.size());
    const std::filesystem::path tmp_path = utils::create_temp_dir(
            std::filesystem::path(dest_dir).remove_filename(), "col");
    std::filesystem::remove_all(tmp_path);
    logger->trace("Using temporary directory {}", tmp_path);

    // stores the set rows for each of the sources, per chunk
    std::vector<std::vector<std::vector<uint64_t>>> set_rows(sources.size());

    for (uint64_t i = 0; i < sources.size(); ++i) {
        if (sources[i]->num_labels() == 0)
            continue;
        uint64_t num_elements
                = std::max((uint64_t)2,  //CWQ needs a buffer size of at least 2
                           std::min((uint64_t)1'000'000, sources[i]->num_relations()));
        targets_size[i].assign(sources[i]->num_labels(), 0U);
        set_rows[i].resize(sources[i]->num_labels());
        for (uint64_t j = 0; j < sources[i]->num_labels(); ++j) {
            const std::filesystem::path tmp_dir
                    = tmp_path/fmt::format("{}/col_{}_{}", i / 100, i, j);
            std::filesystem::create_directories(tmp_dir);
            auto sorted_set = std::make_unique<SSD>(num_threads, num_elements, tmp_dir,
                                                    std::numeric_limits<uint64_t>::max());
            targets[i].push_back(std::move(sorted_set));
        }
    }

    traverse_anno_chunked(
            "Compute diffs", rterminal.size(), graph_fname, sources,
            [&]() {
              for (uint32_t source_idx = 0; source_idx < sources.size(); ++source_idx) {
                  for (uint32_t col_idx = 0; col_idx < set_rows[source_idx].size();
                       ++col_idx) {
                      set_rows[source_idx][col_idx].resize(0);
                  }
              }
            },
            [&](const bit_vector &source_col, node_index row_idx, node_index chunk_idx,
                uint64_t source_idx, uint64_t col_idx, const std::vector<uint64_t> &succ_chunk,
                const std::vector<uint64_t> &pred_chunk,
                const std::vector<uint64_t> &pred_chunk_idx) {

              // check successor node and add current node if it's either terminal
              // or if its successor is 0
              if (rterminal[row_idx] || !source_col[succ_chunk[chunk_idx]])
                  set_rows[source_idx][col_idx].push_back(row_idx);

              // check non-terminal predecessor nodes and add them if they are zero
              for (uint64_t p_idx = pred_chunk_idx[chunk_idx];
                   p_idx < pred_chunk_idx[chunk_idx + 1]; ++p_idx) {
                  if (!source_col[pred_chunk[p_idx]] && !rterminal[pred_chunk[p_idx]])
                      set_rows[source_idx][col_idx].push_back(pred_chunk[p_idx]);
              }
            },
            [&](node_index /* start */, uint64_t /* chunk_size */) {
              for (size_t source_idx = 0; source_idx < sources.size(); ++source_idx) {
                  for (size_t col_idx = 0; col_idx < set_rows[source_idx].size(); ++col_idx) {
                      targets[source_idx][col_idx]->insert(
                              set_rows[source_idx][col_idx].begin(),
                              set_rows[source_idx][col_idx].end());
                      targets_size[source_idx][col_idx]
                              += set_rows[source_idx][col_idx].size();
                  }
              }
            });

    std::vector<LabelEncoder<std::string>> label_encoders;
    std::for_each(sources.begin(), sources.end(), [&](auto& source) {
        label_encoders.push_back(source->get_label_encoder());
    });

    // free memory occupied by sources
    std::vector<std::unique_ptr<ColumnCompressed<>>>().swap(sources);

    logger->trace("Generating row_diff columns...");
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint32_t l_idx = 0; l_idx < targets.size(); ++l_idx) {
        std::vector<std::unique_ptr<bit_vector>> columns(targets[l_idx].size());
        for (uint64_t l_idx2 = 0; l_idx2 < targets[l_idx].size(); ++l_idx2) {
            auto call_ones = [&](const std::function<void(uint64_t)>& call) {
                auto &queue = targets[l_idx][l_idx2]->data(true);
                for (auto &it = queue.begin(); it != queue.end(); ++it) {
                    call(*it);
                }
            };
            columns[l_idx2] = std::make_unique<bit_vector_sd>(call_ones, rterminal.size(),
                                                              targets_size[l_idx][l_idx2]);
        }
        ColumnMajor matrix(std::move(columns));
        auto diff_annotation = std::make_unique<RowDiff<ColumnMajor>>(
                nullptr, std::move(matrix), graph_fname + ".terminal");
        RowDiffAnnotator annotator(std::move(diff_annotation), label_encoders[l_idx]);
        auto fname = std::filesystem::path(source_files[l_idx])
                .filename()
                .replace_extension()
                .replace_extension(RowDiffAnnotator::kExtension);
        auto fpath = dest_dir/fname;
        annotator.serialize(fpath);
        logger->trace("Serialized {}", fpath);
    }
    logger->trace("Removing temp directory: {}", tmp_path);
    std::filesystem::remove_all(tmp_path);
}

} // namespace annot
} // namespace mtg
