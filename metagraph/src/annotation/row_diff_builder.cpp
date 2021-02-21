#include "row_diff_builder.hpp"

#include <omp.h>
#include <progress_bar.hpp>

#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/threads/threading.hpp"
#include "common/elias_fano_file_merger.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/bit_vector_sd.hpp"
#include "graph/annotated_dbg.hpp"

const uint64_t BLOCK_SIZE = 1 << 25;
const uint64_t ROW_REDUCTION_WIDTH = 32;
const uint32_t MAX_NUM_FILES_OPEN = 2000;


namespace mtg {
namespace annot {

using namespace mtg::annot::binmat;
using mtg::common::logger;
namespace fs = std::filesystem;

using anchor_bv_type = RowDiff<ColumnMajor>::anchor_bv_type;
template <typename T>
using Encoder = common::EliasFanoEncoder<T>;


void build_successor(const std::string &graph_fname,
                     const std::string &outfbase,
                     uint32_t max_length,
                     uint32_t num_threads) {
    if (fs::exists(outfbase + ".succ")
        && fs::exists(outfbase + ".pred")
        && fs::exists(outfbase + ".pred_boundary")
        && fs::exists(outfbase + kRowDiffAnchorExt + ".unopt")) {
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
    sdsl::int_vector_buffer succ(outfbase + ".succ", std::ios::out, BLOCK_SIZE, width);
    sdsl::int_vector_buffer<1> pred_boundary(outfbase + ".pred_boundary", std::ios::out, BLOCK_SIZE);
    sdsl::int_vector_buffer pred(outfbase + ".pred", std::ios::out, BLOCK_SIZE, width);

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

void read_next_block(sdsl::int_vector_buffer<>::iterator *succ_it_p,
                     sdsl::int_vector_buffer<1>::iterator *pred_boundary_it_p,
                     sdsl::int_vector_buffer<>::iterator *pred_it_p,
                     uint64_t block_size,
                     std::array<std::vector<uint64_t>, 3> *out) {
    auto &succ_it = *succ_it_p;
    auto &pred_boundary_it = *pred_boundary_it_p;
    auto &pred_it = *pred_it_p;

    std::vector<uint64_t> &succ_chunk = out->at(0);
    std::vector<uint64_t> &pred_chunk_idx = out->at(1);
    std::vector<uint64_t> &pred_chunk = out->at(2);

    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
            succ_chunk.resize(block_size);
            for (uint64_t i = 0; i < block_size; ++i, ++succ_it) {
                succ_chunk[i] = *succ_it;
            }
        }

        #pragma omp section
        {
            // read predecessor offsets
            pred_chunk_idx.resize(block_size + 1);
            pred_chunk_idx[0] = 0;
            for (uint64_t i = 1; i <= block_size; ++i) {
                // find where the last predecessor for the node ends
                pred_chunk_idx[i] = pred_chunk_idx[i - 1];
                while (*pred_boundary_it == 0) {
                    ++pred_chunk_idx[i];
                    ++pred_boundary_it;
                }
                ++pred_boundary_it;
            }

            // read all predecessors for the block
            pred_chunk.resize(pred_chunk_idx.back());
            for (uint64_t i = 0; i < pred_chunk.size(); ++i, ++pred_it) {
                pred_chunk[i] = *pred_it;
            }
        }
    }
}

/**
 * Traverses a group of column compressed annotations (loaded in memory) in chunks of
 * BLOCK_SIZE rows at a time and invokes #call_ones for each set bit.
 * @param num_rows number of rows in the annotation
 * @param pred_succ_fprefix prefix for the pred/succ files containg the predecessors and
 * the successor for each node
 * @param col_annotations the annotations to transform
 * @param before_chunk callback to invoke before a chunk is traversed
 * @param call_ones callback to invoke on a set bit
 * @param after_chunk callback to invoke after a chunk is traversed
 */
void traverse_anno_chunked(
        uint64_t num_rows,
        const std::string &pred_succ_fprefix,
        const std::vector<annot::ColumnCompressed<>> &col_annotations,
        const std::function<void(uint64_t chunk_size)> &before_chunk,
        const CallOnes &call_ones,
        const std::function<void(uint64_t chunk_begin)> &after_chunk) {
    if (col_annotations.empty())
        return;

    const uint32_t num_threads = get_num_threads();

    sdsl::int_vector_buffer succ(pred_succ_fprefix + ".succ", std::ios::in, BLOCK_SIZE);
    sdsl::int_vector_buffer pred(pred_succ_fprefix + ".pred", std::ios::in, BLOCK_SIZE);
    sdsl::int_vector_buffer<1> pred_boundary(pred_succ_fprefix + ".pred_boundary",
                                             std::ios::in, BLOCK_SIZE);

    assert(succ.size() == num_rows);
    assert(static_cast<uint64_t>(std::count(pred_boundary.begin(), pred_boundary.end(), 0))
                   == pred.size());

    auto succ_it = succ.begin();
    auto pred_boundary_it = pred_boundary.begin();
    auto pred_it = pred.begin();

    ThreadPool async_reader(1, 1);
    // start reading the first block
    uint64_t next_block_size = std::min(BLOCK_SIZE, num_rows);
    std::array<std::vector<uint64_t>, 3> context;
    std::array<std::vector<uint64_t>, 3> context_other;
    async_reader.enqueue(read_next_block,
                         &succ_it, &pred_boundary_it, &pred_it, next_block_size,
                         &context_other);

    ProgressBar progress_bar(num_rows, "Compute diffs", std::cerr, !common::get_verbose());

    for (uint64_t chunk = 0; chunk < num_rows; chunk += BLOCK_SIZE) {
        uint64_t block_size = next_block_size;
        next_block_size = std::min(BLOCK_SIZE, num_rows - (chunk + block_size));

        before_chunk(block_size);

        // finish reading this block
        async_reader.join();
        context.swap(context_other);
        std::vector<uint64_t> &succ_chunk = context[0];
        std::vector<uint64_t> &pred_chunk_idx = context[1];
        std::vector<uint64_t> &pred_chunk = context[2];

        // start reading next block
        async_reader.enqueue(read_next_block,
                             &succ_it, &pred_boundary_it, &pred_it, next_block_size,
                             &context_other);

        assert(succ_chunk.size() == block_size);
        assert(pred_chunk.size() == pred_chunk_idx.back());
        // process the current block
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t l_idx = 0; l_idx < col_annotations.size(); ++l_idx) {
            for (size_t j = 0; j < col_annotations[l_idx].num_labels(); ++j) {
                const bit_vector &source_col
                        = *col_annotations[l_idx].get_matrix().data()[j];
                source_col.call_ones_in_range(chunk, chunk + block_size,
                    [&](uint64_t i) {
                        call_ones(source_col, i, i - chunk, l_idx, j,
                                  succ_chunk[i - chunk],
                                  pred_chunk.data() + pred_chunk_idx[i - chunk],
                                  pred_chunk.data() + pred_chunk_idx[i - chunk + 1]);
                    }
                );
            }
        }

        after_chunk(chunk);

        progress_bar += succ_chunk.size();
    }

    if (succ_it != succ.end() || pred_it != pred.end()
                              || pred_boundary_it != pred_boundary.end()) {
        logger->error("Buffers were not read to the end, they might be corrupted");
        exit(1);
    }
}


void convert_batch_to_row_diff(const std::string &pred_succ_fprefix,
                               const std::string &anchors_fname,
                               const std::vector<std::string> &source_files,
                               const fs::path &col_out_dir,
                               const fs::path &swap_dir,
                               const std::string &row_reduction_fname,
                               uint64_t buf_size_bytes,
                               bool compute_row_reduction) {
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
        if (!sources[i].load(source_files[i])) {
            logger->error("Can't load source annotations from {}", source_files[i]);
            std::exit(1);
        }

        if (sources[i].num_labels() && sources[i].num_objects() != anchor.size()) {
            logger->error("Anchor vector {} and annotation {} are incompatible."
                          " Vector size: {}, number of rows: {}",
                          anchors_fname, source_files[i], anchor.size(), sources[i].num_objects());
            std::exit(1);
        }
    }
    logger->trace("Done loading {} annotations", sources.size());

    const fs::path tmp_path = utils::create_temp_dir(swap_dir, "col");

    // stores the row indices that were set because of differences to incoming/outgoing
    // edges, for each of the sources, per chunk. set_rows_fwd is already sorted
    std::vector<std::vector<std::vector<uint64_t>>> set_rows_bwd(sources.size());
    std::vector<std::vector<std::vector<uint64_t>>> set_rows_fwd(sources.size());
    std::vector<std::vector<uint64_t>> row_diff_bits(sources.size());
    std::vector<std::vector<uint64_t>> num_chunks(sources.size());

    auto tmp_file = [&](size_t s, size_t j, size_t chunk) {
        return tmp_path/fmt::format("{}/col_{}_{}/chunk_{}", s / 100, s, j, chunk);
    };
    auto dump_chunk_to_disk = [&](const std::vector<uint64_t> &v,
                                  size_t s, size_t j, size_t chunk) {
        assert(std::is_sorted(v.begin(), v.end()) && "all bits in chunks must be sorted");
        Encoder<uint64_t>::append(v, tmp_file(s, j, chunk));
        row_diff_bits[s][j] += v.size();
    };

    // In the first stage, only one buffer is created per column (`bwd`).
    // In the last stage, two buffers (`fwd` and `bwd`) are created per column.
    const uint64_t buf_size = compute_row_reduction
                                ? buf_size_bytes / sizeof(uint64_t)
                                : buf_size_bytes / sizeof(uint64_t) / 2;

    #pragma omp parallel for num_threads(num_threads)
    for (size_t s = 0; s < sources.size(); ++s) {
        set_rows_fwd[s].resize(sources[s].num_labels());
        set_rows_bwd[s].resize(sources[s].num_labels());
        row_diff_bits[s].assign(sources[s].num_labels(), 0);
        // The first chunk will contain forward bits, all sorted.
        // The other ones (added later) will contain chunks of sorted pred bits.
        num_chunks[s].assign(sources[s].num_labels(), 1);

        for (size_t j = 0; j < sources[s].num_labels(); ++j) {
            fs::create_directories(tmp_file(s, j, 0).parent_path());
            uint64_t original_nbits = sources[s].get_matrix().data()[j]->num_set_bits();
            set_rows_bwd[s][j].reserve(std::min(buf_size, original_nbits));

            if (!compute_row_reduction) {
                // make sure the first chunk exists even if empty
                dump_chunk_to_disk({}, s, j, 0);
                set_rows_fwd[s][j].reserve(std::min(buf_size, original_nbits));
            }
        }
    }

    ThreadPool async_writer(1, 1);
    sdsl::int_vector_buffer row_reduction;
    const bool new_reduction_vector = !fs::exists(row_reduction_fname);
    if (compute_row_reduction) {
        if (new_reduction_vector) {
            // create an empty vector
            sdsl::int_vector_buffer(row_reduction_fname,
                                    std::ios::out, BLOCK_SIZE, ROW_REDUCTION_WIDTH);
        }

        row_reduction = sdsl::int_vector_buffer(row_reduction_fname,
                                                std::ios::in | std::ios::out, BLOCK_SIZE);

        if (!new_reduction_vector && row_reduction.size() != anchor.size()) {
            logger->error("Incompatible sizes of '{}': {} and '{}': {}",
                          row_reduction_fname, row_reduction.size(),
                          anchors_fname, anchor.size());
            exit(1);
        }
    }

    // total number of set bits in the original rows
    std::vector<uint32_t> row_nbits_block;
    // buffer for writing previous block while populating next row_nbits_block
    std::vector<uint32_t> row_nbits_block_other;

    traverse_anno_chunked(
            anchor.size(), pred_succ_fprefix, sources,
            [&](uint64_t chunk_size) {
                if (compute_row_reduction)
                    row_nbits_block.assign(chunk_size, 0);
            },
            [&](const bit_vector &source_col, uint64_t row_idx, uint64_t chunk_idx,
                    size_t source_idx, size_t j,
                    uint64_t succ,
                    const uint64_t *pred_begin, const uint64_t *pred_end) {

                // add current bit if this node is an anchor
                // or if the successor has zero bit
                if (anchor[row_idx] || !source_col[succ]) {
                    // no reduction, we must keep the bit
                    // Push to the buffer only if it's the final stage of the
                    // row-diff transform, where we construct the full row-diff
                    // columns.
                    if (!compute_row_reduction) {
                        auto &v = set_rows_fwd[source_idx][j];
                        v.push_back(row_idx);

                        if (v.size() == v.capacity()) {
                            // dump chunk to disk
                            dump_chunk_to_disk(v, source_idx, j, 0);
                            v.resize(0);
                        }
                    }
                } else if (compute_row_reduction) {
                    // reduction (no bit in row-diff)
                    __atomic_add_fetch(&row_nbits_block[chunk_idx], 1, __ATOMIC_RELAXED);
                }

                // check non-anchor predecessor nodes and add them if they are zero
                for (const uint64_t *pred_p = pred_begin; pred_p < pred_end; ++pred_p) {
                    if (!source_col[*pred_p] && !anchor[*pred_p]) {
                        auto &v = set_rows_bwd[source_idx][j];
                        v.push_back(*pred_p);

                        if (v.size() == v.capacity()) {
                            std::sort(v.begin(), v.end());
                            dump_chunk_to_disk(v, source_idx, j, num_chunks[source_idx][j]++);
                            v.resize(0);
                        }
                    }
                }
            },
            [&](uint64_t block_begin) {
                if (!compute_row_reduction)
                    return;

                __atomic_thread_fence(__ATOMIC_ACQUIRE);

                async_writer.join();
                row_nbits_block_other.swap(row_nbits_block);

                async_writer.enqueue([&,block_begin]() {
                    for (size_t i = 0; i < row_nbits_block_other.size(); ++i) {
                        if (new_reduction_vector) {
                            row_reduction.push_back(row_nbits_block_other[i]);
                        } else {
                            row_reduction[block_begin + i] += row_nbits_block_other[i];
                        }
                    }
                });
            });

    #pragma omp parallel for num_threads(num_threads)
    for (size_t s = 0; s < sources.size(); ++s) {
        for (size_t j = 0; j < sources[s].num_labels(); ++j) {
            auto &fwd = set_rows_fwd[s][j];
            assert(fwd.empty() || !compute_row_reduction);
            if (fwd.size())
                dump_chunk_to_disk(fwd, s, j, 0);

            auto &bwd = set_rows_bwd[s][j];
            if (bwd.size()) {
                std::sort(bwd.begin(), bwd.end());
                dump_chunk_to_disk(bwd, s, j, num_chunks[s][j]++);
            }
        }
    }

    async_writer.join();

    set_rows_fwd.clear(); // free up memory
    set_rows_bwd.clear();

    std::vector<LabelEncoder<std::string>> label_encoders;
    for (const auto &source : sources) {
        label_encoders.push_back(source.get_label_encoder());
    }

    // free memory occupied by original columns
    sources.clear();

    std::vector<std::unique_ptr<RowDiffColumnAnnotator>> row_diff(label_encoders.size());

    const uint32_t files_open_per_thread
            = MAX_NUM_FILES_OPEN / std::max((uint32_t)1, num_threads);
    if (files_open_per_thread < 3) {
        logger->error("Can't merge with less than 3 files per thread open. "
                      "Max num files open: {}. Current number of threads: {}.",
                      MAX_NUM_FILES_OPEN, num_threads);
        exit(1);
    }

    logger->trace("Generating row_diff columns...");
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint32_t l_idx = 0; l_idx < label_encoders.size(); ++l_idx) {
        std::vector<std::unique_ptr<bit_vector>> columns(label_encoders[l_idx].size());

        for (size_t j = 0; j < label_encoders[l_idx].size(); ++j) {
            auto call_ones = [&](const std::function<void(uint64_t)> &call) {
                std::vector<std::string> filenames;
                // skip chunk with fwd bits which have already been counted if stage 1
                for (uint32_t chunk = compute_row_reduction ? 1 : 0;
                                    chunk < num_chunks[l_idx][j]; ++chunk) {
                    filenames.push_back(tmp_file(l_idx, j, chunk));
                }

                // if there are too many chunks, merge them into larger ones
                // TODO: move this pre-merging to common::merge_files
                //       and implement it for SortedSetDisk too.
                while (filenames.size() > files_open_per_thread) {
                    // chunk 0 stores fwd bits and hence not merged
                    std::vector<std::string> new_chunks = { filenames.at(0) };
                    size_t i = 1;
                    while (i < filenames.size()) {
                        std::vector<std::string> to_merge;
                        while (i < filenames.size()
                                && to_merge.size() + 1 < files_open_per_thread) {
                            to_merge.push_back(filenames[i++]);
                        }

                        if (to_merge.size() < 2) {
                            // nothing to merge
                            new_chunks.push_back(to_merge.at(0));
                            continue;
                        }

                        assert(to_merge.size() < files_open_per_thread);

                        new_chunks.push_back(to_merge.at(0) + "_");
                        std::vector<uint64_t> buf;
                        buf.reserve(1024 * 1024);

                        common::merge_files<uint64_t>(to_merge, [&](uint64_t i) {
                            buf.push_back(i);
                            if (buf.size() == buf.capacity()) {
                                Encoder<uint64_t>::append(buf, new_chunks.back());
                                buf.resize(0);
                            }
                        });
                        if (buf.size()) {
                            Encoder<uint64_t>::append(buf, new_chunks.back());
                        }

                        // remove merged chunks
                        for (const std::string &chunk_file : to_merge) {
                            fs::remove(chunk_file);
                        }
                    }
                    filenames.swap(new_chunks);
                }

                common::merge_files<uint64_t>(filenames, call);
            };
            columns[j] = std::make_unique<bit_vector_sd>(call_ones, anchor.size(),
                                                         row_diff_bits[l_idx][j]);
        }

        row_diff[l_idx] = std::make_unique<RowDiffColumnAnnotator>(
                std::make_unique<RowDiff<ColumnMajor>>(nullptr, ColumnMajor(std::move(columns))),
                std::move(label_encoders[l_idx]));

        if (!compute_row_reduction) {
            auto fpath = col_out_dir/fs::path(source_files[l_idx])
                                    .filename()
                                    .replace_extension()
                                    .replace_extension(RowDiffColumnAnnotator::kExtension);

            row_diff[l_idx]->serialize(fpath);

            logger->trace("Serialized {}", fpath);
        }
    }

    utils::remove_temp_dir(tmp_path);

    if (!compute_row_reduction)
        return;

    uint64_t num_larger_rows = 0;
    ProgressBar progress_bar(row_reduction.size(), "Update row reduction",
                             std::cerr, !common::get_verbose());
    for (uint64_t chunk = 0; chunk < row_reduction.size(); chunk += BLOCK_SIZE) {
        row_nbits_block.assign(std::min(BLOCK_SIZE, row_reduction.size() - chunk), 0);

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t l_idx = 0; l_idx < row_diff.size(); ++l_idx) {
            for (const auto &col_ptr : row_diff[l_idx]->get_matrix().diffs().data()) {
                col_ptr->call_ones_in_range(chunk, chunk + row_nbits_block.size(),
                    [&](uint64_t i) {
                        __atomic_add_fetch(&row_nbits_block[i - chunk], 1, __ATOMIC_RELAXED);
                    }
                );
            }
        }

        __atomic_thread_fence(__ATOMIC_ACQUIRE);

        async_writer.join();
        row_nbits_block_other.swap(row_nbits_block);

        async_writer.enqueue([&,chunk]() {
            for (uint64_t i = 0; i < row_nbits_block_other.size(); ++i) {
                row_reduction[chunk + i] -= row_nbits_block_other[i];
                // check if the row reduction is negative
                if (row_reduction[chunk + i] >> (row_reduction.width() - 1))
                    num_larger_rows++;
            }
            progress_bar += row_nbits_block_other.size();
        });
    }

    async_writer.join();

    logger->trace("Rows with negative row reduction: {} in vector {}",
                  num_larger_rows, row_reduction_fname);
}

void optimize_anchors_in_row_diff(const std::string &graph_fname,
                                  const fs::path &dir,
                                  const std::string &row_reduction_extension) {
    if (fs::exists(graph_fname + kRowDiffAnchorExt)) {
        logger->info("Found optimized anchors {}", graph_fname + kRowDiffAnchorExt);
        return;
    }

    logger->trace("Optimizing anchors");

    std::vector<sdsl::int_vector_buffer<>> row_reduction;
    for (const auto &p : fs::directory_iterator(dir)) {
        auto path = p.path();
        if (utils::ends_with(path, row_reduction_extension)) {
            logger->info("Found row reduction vector {}", path);
            row_reduction.emplace_back(path, std::ios::in, BLOCK_SIZE);

            if (row_reduction.back().size() != row_reduction.front().size()) {
                logger->error("Row reduction vectors have different sizes");
                exit(1);
            }
            if (row_reduction.back().width() != row_reduction.front().width()) {
                logger->error("Row reduction vectors have different width: {} != {}",
                              row_reduction.back().width(),
                              row_reduction.front().width());
                exit(1);
            }
        }
    }

    if (!row_reduction.size()) {
        logger->error("Didn't find any row reduction vectors in {} to merge", dir);
        exit(1);
    }

    sdsl::bit_vector anchors;
    uint64_t num_anchors_old;
    {
        anchor_bv_type old_anchors;
        auto original_anchors_fname = graph_fname + kRowDiffAnchorExt + ".unopt";

        std::ifstream f(original_anchors_fname, ios::binary);
        old_anchors.load(f);

        if (old_anchors.size() != row_reduction.front().size()) {
            logger->error(
                "Original anchors {}: {} and row reduction vectors: {} are incompatible",
                original_anchors_fname, old_anchors.size(), row_reduction.front().size());
            exit(1);
        }
        num_anchors_old = old_anchors.num_set_bits();
        anchors = old_anchors.convert_to<sdsl::bit_vector>();
    }

    std::vector<uint64_t> buf;
    static_assert(sizeof(uint64_t) * 8 > ROW_REDUCTION_WIDTH);
    const uint64_t MINUS_BIT = 1llu << (row_reduction.front().width() - 1);

    ProgressBar progress_bar(anchors.size(), "Assign extra anchors", std::cerr,
                             !common::get_verbose());
    for (uint64_t i = 0; i < anchors.size(); i += BLOCK_SIZE) {
        // adjust if the last block
        buf.assign(std::min(BLOCK_SIZE, anchors.size() - i), 0);

        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (size_t j = 0; j < row_reduction.size(); ++j) {
            for (uint64_t t = 0; t < buf.size(); ++t) {
                __atomic_add_fetch(&buf[t], row_reduction[j][i + t], __ATOMIC_RELAXED);
            }
        }

        __atomic_thread_fence(__ATOMIC_ACQUIRE);

        for (uint64_t t = 0; t < buf.size(); ++t) {
            // check if the reduction is negative
            if (buf[t] & MINUS_BIT)
                anchors[i + t] = true;
        }

        progress_bar += buf.size();
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
