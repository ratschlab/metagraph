#include "row_diff_builder.hpp"

#include <ips4o.hpp>
#include <omp.h>
#include <progress_bar.hpp>

#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/threads/threading.hpp"
#include "common/file_merger.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "graph/annotated_dbg.hpp"

constexpr uint64_t BLOCK_SIZE = 1 << 20;
constexpr uint64_t ROW_REDUCTION_WIDTH = 16;

namespace mtg {
namespace annot {

using namespace mtg::annot::binmat;
using mtg::common::logger;

using anchor_bv_type = RowDiff<ColumnMajor>::anchor_bv_type;


void build_successor(const std::string &graph_fname,
                     const std::string &outfbase,
                     uint32_t max_length,
                     uint32_t num_threads) {
    if (std::filesystem::exists(outfbase + ".succ")
        && std::filesystem::exists(outfbase + ".pred")
        && std::filesystem::exists(outfbase + ".pred_boundary")
        && std::filesystem::exists(outfbase + kRowDiffAnchorExt + ".unopt")) {
        logger->trace("Using existing pred/succ/anchors.unopt files in {}.*", outfbase);
        return;

    }
    logger->trace("Building and writing successor, predecessor and anchor files to {}.*",
                  outfbase);

    graph::DBGSuccinct graph(2);
    logger->trace("Loading graph...");
    if (!graph.load(graph_fname)) {
        logger->error("Cannot load graph from {}", graph_fname);
        std::exit(1);
    }

    using graph::boss::BOSS;
    const BOSS &boss = graph.get_boss();
    sdsl::bit_vector terminal;
    sdsl::bit_vector dummy;
    boss.row_diff_traverse(num_threads, max_length, &terminal, &dummy);

    uint64_t num_rows = graph.num_nodes();

    // terminal uses BOSS edges as indices, so we need to map it to annotation indices
    sdsl::bit_vector term(num_rows, 0);
    for (BOSS::edge_index i = 1; i < terminal.size(); ++i) {
        if (terminal[i]) {
            uint64_t graph_idx = graph.boss_to_kmer_index(i);
            uint64_t anno_index
                = graph::AnnotatedSequenceGraph::graph_to_anno_index(graph_idx);
            assert(anno_index < num_rows);
            term[anno_index] = 1;
        }
    }
    logger->trace("Number of anchors before anchor optimization: {}",
                  sdsl::util::cnt_one_bits(term));

    std::ofstream fterm(outfbase + kRowDiffAnchorExt + ".unopt", ios::binary);
    anchor_bv_type(term).serialize(fterm);
    term = sdsl::bit_vector();
    fterm.close();
    logger->trace("Anchor nodes written to {}.unopt", outfbase + kRowDiffAnchorExt);

    // create the succ file, indexed using annotation indices
    uint32_t width = sdsl::bits::hi(graph.num_nodes()) + 1;
    sdsl::int_vector_buffer succ(outfbase + ".succ", std::ios::out, 1024 * 1024, width);
    sdsl::int_vector_buffer<1> pred_boundary(outfbase + ".pred_boundary", std::ios::out, 1024 * 1024);
    sdsl::int_vector_buffer pred(outfbase + ".pred", std::ios::out, 1024 * 1024, width);

    ProgressBar progress_bar(graph.num_nodes(), "Compute successors", std::cerr,
                             !common::get_verbose());

    // traverse BOSS table in parallel processing blocks of size |BLOCK_SIZE|
    for (uint64_t start = 1; start <= graph.num_nodes(); start += BLOCK_SIZE) {
        std::vector<std::vector<uint64_t>> pred_buf(num_threads);
        std::vector<std::vector<uint64_t>> succ_buf(num_threads);
        std::vector<std::vector<bool>> pred_boundary_buf(num_threads);

        // use static scheduling to make threads process ordered contiguous blocks
        #pragma omp parallel for num_threads(num_threads) schedule(static)
        for (uint64_t i = start; i < std::min(start + BLOCK_SIZE, graph.num_nodes() + 1); ++i) {
            const size_t r = omp_get_thread_num();
            BOSS::edge_index boss_idx = graph.kmer_to_boss_index(i);
            if (dummy[boss_idx] || terminal[boss_idx]) {
                succ_buf[r].push_back(num_rows);
            } else {
                const BOSS::TAlphabet d = boss.get_W(boss_idx) % boss.alph_size;
                assert(d && "must not be dummy");
                uint64_t next = graph.boss_to_kmer_index(boss.fwd(boss_idx, d));
                assert(next);
                succ_buf[r].push_back(
                    graph::AnnotatedSequenceGraph::graph_to_anno_index(next)
                );
            }
            // ignore predecessors if boss_idx is not the last outgoing
            // edge (bc. we only traverse the last outgoing at a bifurcation)
            if (!dummy[boss_idx] && boss.get_last(boss_idx)) {
                BOSS::TAlphabet d = boss.get_node_last_value(boss_idx);
                BOSS::edge_index back_idx = boss.bwd(boss_idx);
                boss.call_incoming_to_target(back_idx, d,
                    [&](BOSS::edge_index pred) {
                        // terminal and dummy predecessors are ignored
                        if (!terminal[pred] && !dummy[pred]) {
                            uint64_t node_index = graph.boss_to_kmer_index(pred);
                            pred_buf[r].push_back(
                                graph::AnnotatedSequenceGraph::graph_to_anno_index(node_index)
                            );
                            pred_boundary_buf[r].push_back(0);
                        }
                    }
                );
            }
            pred_boundary_buf[r].push_back(1);
            ++progress_bar;
        }
        for (uint32_t i = 0; i < num_threads; ++i) {
            for (uint64_t v : succ_buf[i]) {
                succ.push_back(v);
            }
            for (uint64_t v : pred_buf[i]) {
                pred.push_back(v);
            }
            for (bool v : pred_boundary_buf[i]) {
                pred_boundary.push_back(v);
            }
        }
    }
    succ.close();
    pred.close();
    pred_boundary.close();

    logger->trace("Pred/succ nodes written to {}.pred/succ", outfbase);
}


/**
 * Callback invoked by #traverse_anno_chunked for each set bit in the annotation matrix.
 * @param source_col the column for which the callback was invoked, in bit_vector format
 * @param row_idx the row in which the bit is set
 * @param row_idx_chunk relative index of the row in the current chunk
 * @param source_idx index of the source file for the current column
 * @param coll_idx index of the column in the current source file (typically, the source
 *        files contain a single column each, but that's not a requirement)
 * @param succ successor of #row_idx in the row-diff path
 * @param pred_begin begin of the predecessor values
 * @param pred_end end of the predecessor values
 */
using CallOnes = std::function<void(const bit_vector &source_col,
                                    uint64_t row_idx,
                                    uint64_t row_idx_chunk,
                                    size_t source_idx,
                                    size_t col_idx,
                                    uint64_t succ,
                                    const uint64_t *pred_begin,
                                    const uint64_t *pred_end)>;

/**
 * Traverses a group of column compressed annotations (loaded in memory) in chunks of
 * BLOCK_SIZE rows at a time and invokes #call_ones for each set bit.
 * @param log_header label to be displayed in the progress bar
 * @param num_rows number of rows in the annotation
 * @param pred_succ_fprefix prefix for the pred/succ files containg the predecessors and
 * the successor for each node
 * @param col_annotations the annotations to transform
 * @param before_chunk callback to invoke before a chunk is traversed
 * @param call_ones callback to invoke on a set bit
 * @param after_chunk callback to invoke after a chunk is traversed
 */
void traverse_anno_chunked(
        const std::string &log_header,
        uint64_t num_rows,
        const std::string &pred_succ_fprefix,
        const std::vector<annot::ColumnCompressed<>> &col_annotations,
        const std::function<void(uint64_t chunk_size)> &before_chunk,
        const CallOnes &call_ones,
        const std::function<void()> &after_chunk) {
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
    std::vector<uint64_t> pred_chunk;
    std::vector<uint64_t> pred_chunk_idx;

    auto pred_boundary_it = pred_boundary.begin();
    auto pred_it = pred.begin();

    ProgressBar progress_bar(num_rows, log_header, std::cerr, !common::get_verbose());

    for (uint64_t chunk = 0; chunk < num_rows; chunk += BLOCK_SIZE) {
        uint64_t block_size = std::min(BLOCK_SIZE, num_rows - chunk);

        before_chunk(block_size);

        succ_chunk.resize(block_size);
        pred_chunk.resize(0);
        pred_chunk_idx.resize(block_size + 1);
        pred_chunk_idx[0] = 0;

        for (uint64_t i = 0; i < block_size; ++i) {
            succ_chunk[i] = succ[chunk + i];
            pred_chunk_idx[i + 1] = pred_chunk_idx[i];
            while (*pred_boundary_it == 0) {
                ++pred_chunk_idx[i + 1];
                pred_chunk.push_back(*pred_it);
                ++pred_it;
                ++pred_boundary_it;
            }
            ++pred_boundary_it;
        }

        assert(pred_chunk.size() == pred_chunk_idx.back());

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t l_idx = 0; l_idx < col_annotations.size(); ++l_idx) {
            for (size_t j = 0; j < col_annotations[l_idx].num_labels(); ++j) {
                const std::unique_ptr<bit_vector> &source_col
                        = col_annotations[l_idx].get_matrix().data()[j];
                source_col->call_ones_in_range(chunk, chunk + block_size,
                    [&](uint64_t i) {
                        call_ones(*source_col, i, i - chunk, l_idx, j,
                                  succ_chunk[i - chunk],
                                  pred_chunk.data() + pred_chunk_idx[i - chunk],
                                  pred_chunk.data() + pred_chunk_idx[i - chunk + 1]);
                    }
                );
            }
        }
        after_chunk();

        progress_bar += succ_chunk.size();
    }
    assert(pred_boundary_it == pred_boundary.end());
}


void convert_batch_to_row_diff(const std::string &pred_succ_fprefix,
                               const std::string &anchors_fname,
                               const std::vector<std::string> &source_files,
                               const std::filesystem::path &dest_dir,
                               const std::string &row_reduction_fname) {
    if (source_files.empty())
        return;

    uint32_t num_threads = get_num_threads();

    anchor_bv_type anchor;
    {
        std::ifstream f(anchors_fname, std::ios::binary);
        anchor.load(f);
    }

    std::vector<annot::ColumnCompressed<>> sources(source_files.size());

    #pragma omp parallel for num_threads(num_threads)
    for (size_t i = 0; i < source_files.size(); ++i) {
        sources[i].load(source_files[i]);

        if (sources[i].num_labels() && sources[i].num_objects() != anchor.size()) {
            logger->error("Anchor vector {} and annotation {} are incompatible."
                          " Vector size: {}, number of rows: {}",
                          anchors_fname, source_files[i], anchor.size(), sources[i].num_objects());
            std::exit(1);
        }
    }
    logger->trace("Done loading {} annotations", sources.size());

    const std::filesystem::path tmp_path = utils::create_temp_dir(
            std::filesystem::path(dest_dir).remove_filename(), "col");
    logger->trace("Using temporary directory {}", tmp_path);

    // stores the row indices that were set because of differences to incoming/outgoing
    // edges, for each of the sources, per chunk. set_rows_fwd is already sorted
    std::vector<std::vector<std::vector<uint64_t>>> set_rows_bwd(sources.size());
    std::vector<std::vector<std::vector<uint64_t>>> set_rows_fwd(sources.size());
    std::vector<std::vector<uint64_t>> col_size(sources.size());
    std::vector<std::vector<uint64_t>> chunks(sources.size());

    auto tmp_dir = [&](uint32_t s, uint32_t c) {
      return tmp_path/fmt::format("{}/col_{}_{}", s / 100, s, c);
    };
    auto tmp_file = [&](uint32_t source_idx, uint32_t col_idx, uint32_t chunk) {
      return tmp_dir(source_idx, col_idx)/fmt::format("chunk_{}", chunk);
    };

    constexpr uint64_t BUF_SIZE = 1'000'000;

    for (size_t i = 0; i < sources.size(); ++i) {
        if (sources[i].num_labels() == 0)
            continue;

        col_size[i].assign(sources[i].num_labels(), 0U);
        chunks[i].assign(sources[i].num_labels(), 0);
        set_rows_fwd[i].resize(sources[i].num_labels());
        set_rows_bwd[i].resize(sources[i].num_labels());

        for (size_t c = 0; c < sources[i].num_labels(); ++c) {
            set_rows_bwd[i][c].reserve(std::min(BUF_SIZE, sources[i].num_relations()));
            set_rows_fwd[i][c].reserve(std::min(BUF_SIZE, sources[i].num_relations()));
            std::filesystem::create_directories(tmp_dir(i, c));
        }
    }

    // total number of set bits in the original rows
    std::vector<uint32_t> row_nbits_batch;
    const std::string temp_row_reduction_fname = row_reduction_fname + ".temp";
    sdsl::int_vector_buffer source_row_nbits(temp_row_reduction_fname, std::ios::out,
                                             1024 * 1024, ROW_REDUCTION_WIDTH);
    traverse_anno_chunked(
            "Compute diffs", anchor.size(), pred_succ_fprefix, sources,
            [&](uint64_t chunk_size) {
                row_nbits_batch.assign(chunk_size, 0);
            },
            [&](const bit_vector &source_col, uint64_t row_idx, uint64_t chunk_idx,
                    size_t source_idx, size_t j,
                    uint64_t succ,
                    const uint64_t *pred_begin, const uint64_t *pred_end) {

                __atomic_add_fetch(&row_nbits_batch[chunk_idx], 1, __ATOMIC_RELAXED);

                // check successor node and add current node if it's either an anchor
                // or if its successor is 0
                if (anchor[row_idx] || !source_col[succ]) {
                    auto &v = set_rows_fwd[source_idx][j];
                    v.push_back(row_idx);
                    if (v.size() == v.capacity()) {
                        std::ofstream f(tmp_dir(source_idx, j)/"chunk_sorted",
                                        std::ios::binary | std::ios::app);
                        f.write(reinterpret_cast<char *>(v.data()), v.size() * sizeof(uint64_t));
                        col_size[source_idx][j] += v.size();
                        v.resize(0);
                    }
                }

                // check non-anchor predecessor nodes and add them if they are zero
                for (const uint64_t *pred_p = pred_begin; pred_p < pred_end; ++pred_p) {
                    if (!source_col[*pred_p] && !anchor[*pred_p]) {
                        auto &v = set_rows_bwd[source_idx][j];
                        v.push_back(*pred_p);
                        if (v.size() == v.capacity()) {
                            ips4o::parallel::sort(v.begin(), v.end(), std::less<uint64_t>(), num_threads);
                            std::ofstream f(tmp_file(source_idx, j, chunks[source_idx][j]),
                                            std::ios::binary);
                            f.write(reinterpret_cast<char *>(v.data()), v.size() * sizeof(uint64_t));
                            col_size[source_idx][j] += v.size();
                            chunks[source_idx][j]++;
                            v.resize(0);
                        }
                    }
                }
            },
            [&]() {
                __atomic_thread_fence(__ATOMIC_ACQUIRE);
                for (uint32_t nbits : row_nbits_batch) {
                    source_row_nbits.push_back(nbits);
                }
            });

    #pragma omp parallel for num_threads(num_threads)
    for (size_t s = 0; s < sources.size(); ++s) {
        for (size_t j = 0; j < set_rows_fwd[s].size(); ++j) {
            auto &fwd = set_rows_fwd[s][j];
            std::ofstream f(tmp_dir(s, j)/"chunk_sorted",
                            std::ios::binary | std::ios::app);
            f.write(reinterpret_cast<char *>(fwd.data()),
                    fwd.size() * sizeof(uint64_t));
            f.close();

            auto &bwd = set_rows_bwd[s][j];
            if (!bwd.empty()) {
                ips4o::parallel::sort(bwd.begin(), bwd.end(), std::less<uint64_t>(), num_threads);
                f.open(tmp_file(s, j, chunks[s][j]), std::ios::binary);
                f.write(reinterpret_cast<char *>(bwd.data()),
                        bwd.size() * sizeof(uint64_t));
                chunks[s][j]++;
            }
            col_size[s][j] += fwd.size() + bwd.size();
        }
    }

    std::vector<LabelEncoder<std::string>> label_encoders;
    for (const auto &source : sources) {
        label_encoders.push_back(source.get_label_encoder());
    }

    // free memory occupied by sources
    sources.clear();

    if (std::filesystem::exists(row_reduction_fname)) {
        logger->trace("Merging row reduction vectors {} and {}",
                      row_reduction_fname, temp_row_reduction_fname);

        sdsl::int_vector_buffer result(row_reduction_fname, std::ios::in | std::ios::out, 1024 * 1024);
        if (result.size() != source_row_nbits.size()) {
            logger->error("Incompatible sizes of '{}': {} and '{}': {}",
                          row_reduction_fname, result.size(),
                          temp_row_reduction_fname, source_row_nbits.size());
            exit(1);
        }
        for (uint64_t i = 0; i < result.size(); ++i) {
            result[i] += source_row_nbits[i];
        }
        source_row_nbits.close();
        std::filesystem::remove(temp_row_reduction_fname);

    } else {
        source_row_nbits.close();
        std::filesystem::rename(temp_row_reduction_fname, row_reduction_fname);
    }

    std::vector<std::unique_ptr<RowDiffAnnotator>> row_diff(label_encoders.size());

    logger->trace("Generating row_diff columns...");
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint32_t l_idx = 0; l_idx < label_encoders.size(); ++l_idx) {
        std::vector<std::unique_ptr<bit_vector>> columns(label_encoders[l_idx].size());
        for (size_t j = 0; j < label_encoders[l_idx].size(); ++j) {
            auto call_ones = [&](const std::function<void(uint64_t)>& call) {
                std::vector<std::string> filenames;
                for (uint32_t chunk = 0; chunk < chunks[l_idx][j]; ++chunk) {
                    filenames.push_back(tmp_file(l_idx, j, chunk));
                }
                filenames.push_back(tmp_dir(l_idx, j)/"chunk_sorted");
                common::merge_files<uint64_t>(filenames, [&](const uint64_t &v) { call(v); });
            };
            columns[j] = std::make_unique<bit_vector_sd>(call_ones, anchor.size(),
                                                         col_size[l_idx][j]);
        }

        ColumnMajor matrix(std::move(columns));
        auto diff_annotation = std::make_unique<RowDiff<ColumnMajor>>(
                nullptr, std::move(matrix));
        row_diff[l_idx] = std::make_unique<RowDiffAnnotator>(std::move(diff_annotation),
                                                             label_encoders[l_idx]);
        auto fname = std::filesystem::path(source_files[l_idx])
                .filename()
                .replace_extension()
                .replace_extension(RowDiffAnnotator::kExtension);
        auto fpath = dest_dir/fname;
        row_diff[l_idx]->serialize(fpath);
        logger->trace("Serialized {}", fpath);
    }
    logger->trace("Removing temp directory: {}", tmp_path);
    std::filesystem::remove_all(tmp_path);

    sdsl::int_vector_buffer row_reduction(row_reduction_fname, std::ios::in | std::ios::out, 1024 * 1024);
    uint64_t num_larger_rows = 0;
    ProgressBar progress_bar(row_reduction.size(), "Update row reduction",
                             std::cerr, !common::get_verbose());
    for (uint64_t chunk = 0; chunk < row_reduction.size(); chunk += BLOCK_SIZE) {
        row_nbits_batch.assign(std::min(BLOCK_SIZE, row_reduction.size() - chunk), 0);

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t l_idx = 0; l_idx < row_diff.size(); ++l_idx) {
            const auto &columns = row_diff[l_idx]->get_matrix().diffs().data();
            for (size_t j = 0; j < columns.size(); ++j) {
                const std::unique_ptr<bit_vector> &source_col = columns[j];
                source_col->call_ones_in_range(chunk, chunk + row_nbits_batch.size(),
                    [&](uint64_t i) {
                        __atomic_add_fetch(&row_nbits_batch[i - chunk], 1, __ATOMIC_RELAXED);
                    }
                );
            }
        }

        __atomic_thread_fence(__ATOMIC_ACQUIRE);

        for (uint64_t i = 0; i < row_nbits_batch.size(); ++i) {
            row_reduction[chunk + i] -= row_nbits_batch[i];
            // check if the row reduction is negative
            if (row_reduction[chunk + i] >> (row_reduction.width() - 1))
                num_larger_rows++;
        }

        progress_bar += row_nbits_batch.size();
    }

    logger->trace("Rows with negative row reduction: {} in vector {}",
                  num_larger_rows, row_reduction_fname);
}

void optimize_anchors_in_row_diff(const std::string &graph_fname,
                                  const std::filesystem::path &dest_dir,
                                  const std::string &row_reduction_extension) {
    if (std::filesystem::exists(graph_fname + kRowDiffAnchorExt)) {
        logger->info("Found optimized anchors {}", graph_fname + kRowDiffAnchorExt);
        return;
    }

    logger->trace("Optimizing anchors");

    std::vector<std::string> filenames;
    for (const auto &p : std::filesystem::directory_iterator(dest_dir)) {
        auto path = p.path();
        if (utils::ends_with(path, row_reduction_extension)) {
            logger->info("Found row reduction vector {}", path);
            filenames.push_back(path);
        }
    }

    if (!filenames.size()) {
        logger->error("Didn't find any row reduction vectors in {} to merge", dest_dir);
        exit(1);
    }

    if (filenames.size() > 1u)
        logger->trace("Merging row reduction vectors");

    while (filenames.size() > 1u) {
        std::vector<std::string> filenames_new((filenames.size() + 1) / 2);

        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (size_t t = 1; t < filenames.size(); t += 2) {
            // compute sum of t-1 and t.
            sdsl::int_vector_buffer first(filenames[t - 1], std::ios::in, 1024 * 1024);
            sdsl::int_vector_buffer second(filenames[t], std::ios::in, 1024 * 1024);
            if (first.size() != second.size() || first.width() != second.width()) {
                logger->error("Sizes of row reduction vectors are incompatible, {}: {}, {}: {}",
                              filenames[t - 1], first.size(), filenames[t], second.size());
                exit(1);
            }
            filenames_new[t / 2] = fmt::format("{}.merged", filenames[t]);
            sdsl::int_vector_buffer sum(filenames_new[t / 2], std::ios::out, 1024 * 1024, first.width());

            for (uint64_t i = 0; i < first.size(); ++i) {
                sum.push_back(first[i] + second[i]);
            }
        }

        if (filenames.size() % 2)
            filenames_new.back() = filenames.back();

        logger->trace("Merged {} row reduction vectors into {}",
                      filenames.size(), filenames_new.size());

        filenames.swap(filenames_new);
    }

    assert(filenames.size() == 1u && "All must be merged into one vector");

    anchor_bv_type old_anchors;
    auto original_anchors_fname = graph_fname + kRowDiffAnchorExt + ".unopt";
    {
        std::ifstream f(original_anchors_fname, ios::binary);
        old_anchors.load(f);
    }
    uint64_t num_anchors_old = old_anchors.num_set_bits();

    sdsl::int_vector_buffer row_reduction(filenames[0], std::ios::in, 1024 * 1024);

    if (old_anchors.size() != row_reduction.size()) {
        logger->error(
            "Original anchors {}: {} and row reduction vector {}: {} are incompatible",
            original_anchors_fname, old_anchors.size(), filenames[0], row_reduction.size());
        exit(1);
    }
    sdsl::bit_vector anchors = old_anchors.convert_to<sdsl::bit_vector>();
    for (uint64_t i = 0; i < row_reduction.size(); ++i) {
        // check if the reduction is negative
        if (row_reduction[i] >> (row_reduction.width() - 1))
            anchors[i] = true;
    }
    anchor_bv_type ranchors(std::move(anchors));

    logger->trace("Number of anchors increased after optimization from {} to {}",
                  num_anchors_old, ranchors.num_set_bits());

    std::ofstream f(graph_fname + kRowDiffAnchorExt, ios::binary);
    ranchors.serialize(f);
    logger->trace("Serialized optimized anchors to {}", graph_fname + kRowDiffAnchorExt);
}

} // namespace annot
} // namespace mtg
