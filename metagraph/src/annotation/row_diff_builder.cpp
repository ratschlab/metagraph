#include "row_diff_builder.hpp"

#include <omp.h>
#include <progress_bar.hpp>

#include "annotation/binary_matrix/row_diff/row_diff.hpp"
#include "annotation/int_matrix/row_diff/int_row_diff.hpp"
#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "common/threads/threading.hpp"
#include "common/elias_fano/elias_fano_merger.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/bit_vector_sdsl.hpp"
#include "graph/annotated_dbg.hpp"

const uint64_t BLOCK_SIZE = 1 << 25;
const uint64_t BUFFER_SIZE = 1024 * 1024; // 1 MiB
const uint64_t ROW_REDUCTION_WIDTH = 32;
const uint32_t MAX_NUM_FILES_OPEN = 2000;


namespace mtg {
namespace annot {

using namespace mtg::annot::binmat;
using mtg::common::logger;
using mtg::graph::boss::BOSS;
namespace fs = std::filesystem;

using anchor_bv_type = RowDiff<ColumnMajor>::anchor_bv_type;
using rd_succ_bv_type = RowDiff<ColumnMajor>::fork_succ_bv_type;
template <typename T>
using Encoder = mtg::elias_fano::EliasFanoEncoderBuffered<T>;


std::vector<annot::ColumnCompressed<>>
load_columns(const std::vector<std::string> &source_files, uint64_t *num_rows) {
    *num_rows = 0;

    std::vector<annot::ColumnCompressed<>> sources(source_files.size());

    #pragma omp parallel for num_threads(get_num_threads())
    for (size_t i = 0; i < source_files.size(); ++i) {
        if (!sources[i].load(source_files[i])) {
            logger->error("Can't load source annotations from {}", source_files[i]);
            std::exit(1);
        }

        if (!sources[i].num_labels())
            continue;

        #pragma omp critical
        {
            if (!*num_rows) {
                *num_rows = sources[i].num_objects();
            } else if (*num_rows != sources[i].num_objects()) {
                logger->error("Annotations have different number of rows");
                std::exit(1);
            }
        }
    }
    logger->trace("Done loading {} annotations", sources.size());

    return sources;
}

template <class Callback>
void load_coordinates(const std::vector<std::string> &source_files,
                      const std::vector<annot::ColumnCompressed<>> &sources,
                      Callback callback,
                      size_t num_threads) {
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t i = 0; i < source_files.size(); ++i) {
        if (!sources[i].num_labels())
            continue;

        const auto &coords_fname = utils::remove_suffix(source_files[i],
                                                        ColumnCompressed<>::kExtension)
                                        + ColumnCompressed<>::kCoordExtension;
        std::ifstream in(coords_fname, std::ios::binary);
        if (!in) {
            logger->error("Could not open file with coordinates {}", coords_fname);
            exit(1);
        }

        sdsl::int_vector<> coords;
        bit_vector_smart delims;
        for (size_t j = 0; j < sources[i].num_labels(); ++j) {
            try {
                delims.load(in);
                coords.load(in);
            } catch (...) {
                logger->error("Couldn't read coordinates from {}", coords_fname);
                exit(1);
            }
            callback(i, std::move(coords), std::move(delims));
        }
    }
}

void count_labels_per_row(const std::vector<std::string> &source_files,
                          const std::string &row_count_fname,
                          bool with_coordinates) {
    if (source_files.empty())
        return;

    uint64_t num_rows;
    std::vector<annot::ColumnCompressed<>> sources = load_columns(source_files, &num_rows);
    if (!num_rows) {
        logger->warn("Input annotations have no rows");
        return;
    }

    std::vector<std::vector<bit_vector_smart>> delims(source_files.size());
    if (with_coordinates) {
        load_coordinates(source_files, sources,
            [&](size_t s, sdsl::int_vector<>&&, bit_vector_smart&& d) {
                delims[s].push_back(std::move(d));
            },
            get_num_threads()
        );
        logger->trace("Done loading coordinates");
    }

    // initializing the row count vector
    const bool new_vector = !fs::exists(row_count_fname);
    logger->trace("Row count vector: {}", row_count_fname);
    if (new_vector) {
        // create an empty vector
        sdsl::int_vector_buffer<>(row_count_fname,
                                  std::ios::out, BUFFER_SIZE, ROW_REDUCTION_WIDTH);
        logger->trace("Initialized new row count vector {}", row_count_fname);
    } else {
        logger->trace("Row count vector {} already exists and will be updated",
                      row_count_fname);
    }
    sdsl::int_vector_buffer<> row_count(row_count_fname,
                                        std::ios::in | std::ios::out, BUFFER_SIZE);

    // total number of set bits in the original rows
    std::vector<uint32_t> row_count_block;
    // buffer for writing previous block while populating next row_count_block
    std::vector<uint32_t> row_count_block_other;

    ThreadPool async_writer(1, 1);

    ProgressBar progress_bar(num_rows, "Count row labels", std::cerr, !common::get_verbose());
    uint64_t next_block_size = std::min(BLOCK_SIZE, num_rows);

    for (uint64_t block_begin = 0; block_begin < num_rows; block_begin += BLOCK_SIZE) {
        uint64_t block_size = next_block_size;
        next_block_size = std::min(BLOCK_SIZE, num_rows - (block_begin + block_size));

        row_count_block.assign(block_size, 0);

        // process the current block
        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (size_t l_idx = 0; l_idx < sources.size(); ++l_idx) {
            for (size_t j = 0; j < sources[l_idx].num_labels(); ++j) {
                const bit_vector &source_col
                        = *sources[l_idx].get_matrix().data()[j];
                source_col.call_ones_in_range(block_begin, block_begin + block_size,
                    [&](uint64_t i) {
                        if (with_coordinates) {
                            uint64_t rk = source_col.rank1(i);
                            uint32_t num_coords = delims[l_idx][j].select1(rk + 1)
                                                    - delims[l_idx][j].select1(rk);
                            __atomic_add_fetch(&row_count_block[i - block_begin], num_coords, __ATOMIC_RELAXED);
                        } else {
                            __atomic_add_fetch(&row_count_block[i - block_begin], 1, __ATOMIC_RELAXED);
                        }
                    }
                );
            }
        }

        __atomic_thread_fence(__ATOMIC_ACQUIRE);

        async_writer.join();
        row_count_block_other.swap(row_count_block);

        async_writer.enqueue([&,block_begin]() {
            for (size_t i = 0; i < row_count_block_other.size(); ++i) {
                uint32_t c = row_count_block_other[i];
                if (new_vector) {
                    row_count.push_back(c);
                } else if (c) {
                    row_count[block_begin + i] += c;
                }
            }
        });

        progress_bar += block_size;
    }

    async_writer.join();
}


// convert index
inline uint64_t to_row(graph::DeBruijnGraph::node_index i) {
    return graph::AnnotatedSequenceGraph::graph_to_anno_index(i);
}
inline graph::DeBruijnGraph::node_index to_node(uint64_t i) {
    return graph::AnnotatedSequenceGraph::anno_to_graph_index(i);
}

template <class Callback>
void sum_and_call_counts(const fs::path &dir,
                         const std::string &file_extension,
                         const std::string &counts_name,
                         const Callback &callback) {
    std::vector<sdsl::int_vector_buffer<>> vectors;
    for (const auto &p : fs::directory_iterator(dir)) {
        auto path = p.path();
        if (utils::ends_with(path, file_extension)) {
            logger->trace("Found count vector {}", path);
            vectors.emplace_back(path, std::ios::in, BUFFER_SIZE);

            if (vectors.back().size() != vectors.front().size()) {
                logger->error("Count vectors have different sizes");
                exit(1);
            }
            if (vectors.back().width() != 32) {
                logger->error("Count vectors must have width 32, but {} has {}",
                              path, vectors.back().width());
                exit(1);
            }
        }
    }

    if (!vectors.size()) {
        logger->error("Didn't find any count vectors to merge in {}", dir);
        exit(1);
    }

    std::vector<int32_t> buf;

    ProgressBar progress_bar(vectors.front().size(), "Sum " + counts_name,
                             std::cerr, !common::get_verbose());

    for (uint64_t i = 0; i < vectors.front().size(); i += BLOCK_SIZE) {
        // adjust if the last block
        buf.assign(std::min(BLOCK_SIZE, vectors.front().size() - i), 0);

        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
        for (size_t j = 0; j < vectors.size(); ++j) {
            for (uint64_t t = 0; t < buf.size(); ++t) {
                __atomic_add_fetch(&buf[t], (int32_t)vectors[j][i + t], __ATOMIC_RELAXED);
            }
        }

        __atomic_thread_fence(__ATOMIC_ACQUIRE);

        std::for_each(buf.begin(), buf.end(), callback);

        progress_bar += buf.size();
    }
}

rd_succ_bv_type route_at_forks(const graph::DBGSuccinct &graph,
                               const std::string &rd_succ_filename,
                               const std::string &count_vectors_dir,
                               const std::string &row_count_extension) {
    logger->trace("Assigning row-diff successors at forks...");

    rd_succ_bv_type rd_succ;

    bool optimize_forks = false;
    for (const auto &p : fs::directory_iterator(count_vectors_dir)) {
        if (utils::ends_with(p.path(), row_count_extension))
            optimize_forks = true;
    }

    if (optimize_forks) {
        logger->trace("RowDiff successors will be set to the adjacent nodes with"
                      " the largest number of labels");

        const bit_vector &last = graph.get_boss().get_last();
        graph::DeBruijnGraph::node_index graph_idx = to_node(0);

        std::vector<uint32_t> outgoing_counts;

        sdsl::bit_vector rd_succ_bv(last.size(), false);

        sum_and_call_counts(count_vectors_dir, row_count_extension, "row counts",
            [&](int32_t count) {
                // TODO: skip single outgoing
                outgoing_counts.push_back(count);
                if (last[graph.kmer_to_boss_index(graph_idx)]) {
                    // pick the node with the largest count
                    size_t max_pos = std::max_element(outgoing_counts.rbegin(),
                                                      outgoing_counts.rend())
                                     - outgoing_counts.rbegin();
                    rd_succ_bv[graph.kmer_to_boss_index(graph_idx - max_pos)] = true;
                    outgoing_counts.resize(0);
                }
                graph_idx++;
            }
        );

        if (graph_idx != graph.num_nodes() + 1) {
            logger->error("Size the count vectors is incompatible with the"
                          " graph: {} != {}", graph_idx - 1, graph.num_nodes());
            exit(1);
        }

        rd_succ = rd_succ_bv_type(std::move(rd_succ_bv));

    } else {
        logger->warn("No count vectors could be found in {}. The last outgoing"
                     " edges will be selected for assigning RowDiff successors",
                     count_vectors_dir);
    }

    std::ofstream f(rd_succ_filename, ios::binary);
    rd_succ.serialize(f);
    logger->trace("RowDiff successors are assigned for forks and written to {}",
                  rd_succ_filename);
    return rd_succ;
}

void build_pred_succ(const std::string &graph_fname,
                     const std::string &outfbase,
                     const std::string &count_vectors_dir,
                     const std::string &row_count_extension,
                     uint32_t num_threads) {
    if (fs::exists(outfbase + ".succ")
            && fs::exists(outfbase + ".succ_boundary")
            && fs::exists(outfbase + ".pred")
            && fs::exists(outfbase + ".pred_boundary")
            && fs::exists(outfbase + kRowDiffForkSuccExt)) {
        logger->trace("Using existing pred/succ files in {}.*", outfbase);
        return;
    }
    logger->trace("Building and writing successor and predecessor files to {}.*",
                  outfbase);

    graph::DBGSuccinct graph(2);
    logger->trace("Loading graph...");
    if (!graph.load(graph_fname)) {
        logger->error("Cannot load graph from {}", graph_fname);
        std::exit(1);
    }

    // assign row-diff successors at forks
    rd_succ_bv_type rd_succ = route_at_forks(graph, outfbase + kRowDiffForkSuccExt,
                                             count_vectors_dir, row_count_extension);

    const BOSS &boss = graph.get_boss();

    sdsl::bit_vector dummy = boss.mark_all_dummy_edges(num_threads);

    // create the succ/pred files, indexed using annotation indices
    uint32_t width = sdsl::bits::hi(graph.num_nodes()) + 1;
    sdsl::int_vector_buffer<> succ(outfbase + ".succ", std::ios::out, BUFFER_SIZE, width);
    sdsl::int_vector_buffer<1> succ_boundary(outfbase + ".succ_boundary", std::ios::out, BUFFER_SIZE);
    sdsl::int_vector_buffer<> pred(outfbase + ".pred", std::ios::out, BUFFER_SIZE, width);
    sdsl::int_vector_buffer<1> pred_boundary(outfbase + ".pred_boundary", std::ios::out, BUFFER_SIZE);

    ProgressBar progress_bar(graph.num_nodes(), "Compute succ/pred", std::cerr,
                             !common::get_verbose());

    const uint64_t BS = 1'000'000;
    // traverse BOSS table in parallel processing blocks of size |BS|
    // use static scheduling to make threads process ordered contiguous blocks
    #pragma omp parallel for ordered num_threads(num_threads) schedule(dynamic)
    for (uint64_t start = 1; start <= graph.num_nodes(); start += BS) {
        std::vector<uint64_t> succ_buf;
        std::vector<bool> succ_boundary_buf;
        std::vector<uint64_t> pred_buf;
        std::vector<bool> pred_boundary_buf;

        for (uint64_t i = start; i < std::min(start + BS, graph.num_nodes() + 1); ++i) {
            BOSS::edge_index boss_idx = graph.kmer_to_boss_index(i);
            if (!dummy[boss_idx]) {
                const BOSS::TAlphabet d = boss.get_W(boss_idx) % boss.alph_size;
                assert(d && "must not be dummy");
                BOSS::edge_index next = boss.fwd(boss_idx, d);
                assert(next);
                if (!dummy[next]) {
                    while (rd_succ.size() && !rd_succ[next]) {
                        next--;
                        assert(!boss.get_last(next));
                    }
                    succ_buf.push_back(to_row(graph.boss_to_kmer_index(next)));
                    succ_boundary_buf.push_back(0);
                }
                // compute predecessors only for row-diff successors
                if (rd_succ.size() ? rd_succ[boss_idx] : boss.get_last(boss_idx)) {
                    BOSS::TAlphabet d = boss.get_node_last_value(boss_idx);
                    BOSS::edge_index back_idx = boss.bwd(boss_idx);
                    boss.call_incoming_to_target(back_idx, d,
                        [&](BOSS::edge_index pred) {
                            // dummy predecessors are ignored
                            if (!dummy[pred]) {
                                uint64_t node_index = graph.boss_to_kmer_index(pred);
                                pred_buf.push_back(to_row(node_index));
                                pred_boundary_buf.push_back(0);
                            }
                        }
                    );
                }
            }
            succ_boundary_buf.push_back(1);
            pred_boundary_buf.push_back(1);
            ++progress_bar;
        }

        #pragma omp ordered
        {
            // append to the files on disk
            for (uint64_t v : succ_buf) { succ.push_back(v); }
            for (bool v : succ_boundary_buf) { succ_boundary.push_back(v); }
            for (uint64_t v : pred_buf) { pred.push_back(v); }
            for (bool v : pred_boundary_buf) { pred_boundary.push_back(v); }
        }
    }

    logger->trace("Pred/succ nodes written to {}.pred/succ", outfbase);
}

void assign_anchors(const std::string &graph_fname,
                    const std::string &outfbase,
                    const std::filesystem::path &count_vectors_dir,
                    uint32_t max_length,
                    const std::string &row_reduction_extension,
                    uint32_t num_threads) {
    std::string anchor_filename = outfbase + kRowDiffAnchorExt;
    if (fs::exists(anchor_filename)) {
        logger->trace("Using existing anchors {}", anchor_filename);
        return;
    }

    graph::DBGSuccinct graph(2);
    logger->trace("Loading graph...");
    if (!graph.load(graph_fname)) {
        logger->error("Cannot load graph from {}", graph_fname);
        std::exit(1);
    }
    const BOSS &boss = graph.get_boss();
    const uint64_t num_rows = graph.num_nodes();

    bool optimize_anchors = false;
    for (const auto &p : fs::directory_iterator(count_vectors_dir)) {
        if (utils::ends_with(p.path(), row_reduction_extension))
            optimize_anchors = true;
    }

    sdsl::bit_vector anchors_bv(boss.get_last().size(), false);

    if (optimize_anchors) {
        logger->trace("Making every row with negative reduction an anchor...");

        uint64_t i = 0;
        sum_and_call_counts(count_vectors_dir, row_reduction_extension, "row reduction",
            [&](int32_t count) {
                // check if the reduction is negative
                if (count < 0)
                    anchors_bv[graph.kmer_to_boss_index(to_node(i))] = true;
                i++;
            }
        );

        if (i != num_rows) {
            logger->error("Reduction vectors are incompatible with the graph size:"
                          " {} != {}", i, num_rows);
            exit(1);
        }

        logger->trace("Number of initial anchors (rows with negative reduction): {}",
                      sdsl::util::cnt_one_bits(anchors_bv));
    } else {
        logger->warn("Didn't find any precomputed row reduction vectors in {}."
                     " Anchors will be assigned from scratch.", count_vectors_dir);
    }

    // assign extra anchors and restrict the length of row-diff paths
    logger->trace("Assigning required anchors...");
    {
        rd_succ_bv_type rd_succ;
        const std::string &rd_succ_fname = outfbase + kRowDiffForkSuccExt;
        std::ifstream f(rd_succ_fname, ios::binary);
        if (!rd_succ.load(f)) {
            logger->error("Couldn't load row-diff successor bitmap from {}", rd_succ_fname);
            exit(1);
        }

        if (rd_succ.size()) {
            logger->trace("Assigning anchors for RowDiff successors {}...", rd_succ_fname);
            boss.row_diff_traverse(num_threads, max_length, rd_succ, &anchors_bv);
        } else {
            logger->warn("Assigning anchors without chosen RowDiff successors."
                         " The last outgoing edges will be used for routing.");
            boss.row_diff_traverse(num_threads, max_length, boss.get_last(), &anchors_bv);
        }
    }

    // anchors_bv uses BOSS edges as indices, so we need to map it to annotation indices
    {
        sdsl::bit_vector anchors(num_rows, false);
        for (BOSS::edge_index i = 1; i < anchors_bv.size(); ++i) {
            if (anchors_bv[i]) {
                uint64_t graph_idx = graph.boss_to_kmer_index(i);
                assert(to_row(graph_idx) < num_rows);
                anchors[to_row(graph_idx)] = 1;
            }
        }
        anchors_bv = std::move(anchors);
    }

    anchor_bv_type anchors(std::move(anchors_bv));
    logger->trace("Final number of anchors in row-diff: {}", anchors.num_set_bits());

    std::ofstream f(anchor_filename, ios::binary);
    anchors.serialize(f);
    logger->trace("Serialized anchors to {}", anchor_filename);
}


/**
 * Callback invoked by #traverse_anno_chunked for each set bit in the annotation matrix.
 * @param source_col the column for which the callback was invoked, in bit_vector format
 * @param row_idx the row in which the bit is set
 * @param row_idx_chunk relative index of the row in the current chunk
 * @param source_idx index of the source file for the current column
 * @param coll_idx index of the column in the current source file (typically, the source
 *        files contain a single column each, but that's not a requirement)
 * @param succ_begin begin of the successor values
 * @param succ_end end of the successor values
 * @param pred_begin begin of the predecessor values
 * @param pred_end end of the predecessor values
 */
using CallOnes = std::function<void(const bit_vector &source_col,
                                    uint64_t row_idx,
                                    uint64_t row_idx_chunk,
                                    size_t source_idx,
                                    size_t col_idx,
                                    const uint64_t *succ,
                                    const uint64_t *pred_begin, const uint64_t *pred_end)>;

void read_next_block(sdsl::int_vector_buffer<>::iterator &it,
                     sdsl::int_vector_buffer<1>::iterator &boundary_it,
                     uint64_t block_size,
                     std::vector<uint64_t> &chunk,
                     std::vector<uint64_t> &chunk_idx) {
    // read offsets
    chunk_idx.resize(block_size + 1);
    chunk_idx[0] = 0;
    for (uint64_t i = 1; i <= block_size; ++i) {
        // find where the last element for the node ends
        chunk_idx[i] = chunk_idx[i - 1];
        while (*boundary_it == 0) {
            ++chunk_idx[i];
            ++boundary_it;
        }
        ++boundary_it;
    }
    // read all elements from the block
    chunk.resize(chunk_idx.back());
    for (uint64_t i = 0; i < chunk.size(); ++i, ++it) {
        chunk[i] = *it;
    }
}

void read_next_blocks(sdsl::int_vector_buffer<>::iterator *succ_it_p,
                      sdsl::int_vector_buffer<1>::iterator *succ_boundary_it_p,
                      sdsl::int_vector_buffer<>::iterator *pred_it_p,
                      sdsl::int_vector_buffer<1>::iterator *pred_boundary_it_p,
                      uint64_t block_size,
                      std::array<std::vector<uint64_t>, 4> *out) {
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        read_next_block(*succ_it_p, *succ_boundary_it_p, block_size,
                        out->at(0), out->at(1));
        #pragma omp section
        read_next_block(*pred_it_p, *pred_boundary_it_p, block_size,
                        out->at(2), out->at(3));
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

    sdsl::int_vector_buffer<> succ(pred_succ_fprefix + ".succ", std::ios::in, BUFFER_SIZE);
    sdsl::int_vector_buffer<1> succ_boundary(pred_succ_fprefix + ".succ_boundary",
                                             std::ios::in, BUFFER_SIZE);
    sdsl::int_vector_buffer<> pred(pred_succ_fprefix + ".pred", std::ios::in, BUFFER_SIZE);
    sdsl::int_vector_buffer<1> pred_boundary(pred_succ_fprefix + ".pred_boundary",
                                             std::ios::in, BUFFER_SIZE);
    auto succ_it = succ.begin();
    auto succ_boundary_it = succ_boundary.begin();
    auto pred_it = pred.begin();
    auto pred_boundary_it = pred_boundary.begin();

    ThreadPool async_reader(1, 1);
    // start reading the first block
    uint64_t next_block_size = std::min(BLOCK_SIZE, num_rows);
    std::array<std::vector<uint64_t>, 4> context;
    std::array<std::vector<uint64_t>, 4> context_other;
    async_reader.enqueue(read_next_blocks,
                         &succ_it, &succ_boundary_it,
                         &pred_it, &pred_boundary_it,
                         next_block_size, &context_other);

    ProgressBar progress_bar(num_rows, "Compute diffs", std::cerr, !common::get_verbose());

    for (uint64_t chunk = 0; chunk < num_rows; chunk += BLOCK_SIZE) {
        uint64_t block_size = next_block_size;
        next_block_size = std::min(BLOCK_SIZE, num_rows - (chunk + block_size));

        before_chunk(block_size);

        // finish reading this block
        async_reader.join();
        context.swap(context_other);
        std::vector<uint64_t> &succ_chunk = context[0];
        std::vector<uint64_t> &succ_chunk_idx = context[1];
        std::vector<uint64_t> &pred_chunk = context[2];
        std::vector<uint64_t> &pred_chunk_idx = context[3];

        // start reading next block
        async_reader.enqueue(read_next_blocks,
                             &succ_it, &succ_boundary_it,
                             &pred_it, &pred_boundary_it,
                             next_block_size, &context_other);

        assert(succ_chunk.size() == succ_chunk_idx.back());
        assert(pred_chunk.size() == pred_chunk_idx.back());
        // process the current block
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t l_idx = 0; l_idx < col_annotations.size(); ++l_idx) {
            for (size_t j = 0; j < col_annotations[l_idx].num_labels(); ++j) {
                const bit_vector &source_col
                        = *col_annotations[l_idx].get_matrix().data()[j];
                source_col.call_ones_in_range(chunk, chunk + block_size,
                    [&](uint64_t i) {
                        assert(succ_chunk_idx[i - chunk + 1] >= succ_chunk_idx[i - chunk]);
                        assert(succ_chunk_idx[i - chunk + 1] <= succ_chunk_idx[i - chunk] + 1);
                        const uint64_t *succ = succ_chunk_idx[i - chunk + 1]
                                                > succ_chunk_idx[i - chunk]
                                                ? succ_chunk.data() + succ_chunk_idx[i - chunk]
                                                : NULL;
                        call_ones(source_col, i, i - chunk, l_idx, j, succ,
                                  pred_chunk.data() + pred_chunk_idx[i - chunk],
                                  pred_chunk.data() + pred_chunk_idx[i - chunk + 1]);
                    }
                );
            }
        }

        after_chunk(chunk);

        progress_bar += block_size;
    }

    if (succ_it != succ.end()
            || succ_boundary_it != succ_boundary.end()
            || pred_it != pred.end()
            || pred_boundary_it != pred_boundary.end()) {
        logger->error("Buffers were not read to the end, they might be corrupted");
        exit(1);
    }
}


template <typename T = uint64_t>
void convert_batch_to_row_diff(const std::string &pred_succ_fprefix,
                               const std::vector<std::string> &source_files,
                               const fs::path &col_out_dir,
                               const fs::path &swap_dir,
                               const std::string &row_reduction_fname,
                               uint64_t buf_size_bytes,
                               bool compute_row_reduction);

void convert_batch_to_row_diff_coord(const std::string &pred_succ_fprefix,
                                     const std::vector<std::string> &source_files,
                                     const fs::path &col_out_dir,
                                     const fs::path &swap_dir,
                                     const std::string &row_reduction_fname,
                                     uint64_t buf_size_bytes,
                                     bool compute_row_reduction,
                                     size_t num_coords_per_seq);

void convert_batch_to_row_diff(const std::string &pred_succ_fprefix,
                               const std::vector<std::string> &source_files,
                               const fs::path &col_out_dir,
                               const fs::path &swap_dir,
                               const std::string &row_reduction_fname,
                               uint64_t buf_size_bytes,
                               bool compute_row_reduction,
                               bool with_values,
                               bool with_coordinates,
                               size_t num_coords_per_seq) {
    if (with_values) {
        convert_batch_to_row_diff<std::pair<uint64_t, uint64_t>>(
                pred_succ_fprefix, source_files, col_out_dir, swap_dir,
                row_reduction_fname, buf_size_bytes, compute_row_reduction);
    } else if (with_coordinates) {
        convert_batch_to_row_diff_coord(
                pred_succ_fprefix, source_files, col_out_dir, swap_dir,
                row_reduction_fname, buf_size_bytes, compute_row_reduction,
                num_coords_per_seq);
    } else {
        convert_batch_to_row_diff<uint64_t>(
                pred_succ_fprefix, source_files, col_out_dir, swap_dir,
                row_reduction_fname, buf_size_bytes, compute_row_reduction);
    }
}

// 'T' is either a row index, or a pair (row index, value)
template <typename T>
void convert_batch_to_row_diff(const std::string &pred_succ_fprefix,
                               const std::vector<std::string> &source_files,
                               const fs::path &col_out_dir,
                               const fs::path &swap_dir,
                               const std::string &row_reduction_fname,
                               uint64_t buf_size_bytes,
                               bool compute_row_reduction) {
    if (source_files.empty())
        return;

    const uint32_t num_threads = get_num_threads();

    uint64_t num_rows;
    std::vector<annot::ColumnCompressed<>> sources = load_columns(source_files, &num_rows);
    if (!num_rows) {
        logger->warn("Input annotations have no rows");
        return;
    }

    constexpr bool with_values = utils::is_pair_v<T>;
    std::vector<std::vector<sdsl::int_vector<>>> values;
    if (with_values) {
        values.resize(source_files.size());
        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 0; i < source_files.size(); ++i) {
            if (!sources[i].num_labels())
                continue;

            const auto &values_fname = utils::remove_suffix(source_files[i],
                                                            ColumnCompressed<>::kExtension)
                                            + ColumnCompressed<>::kCountExtension;
            std::ifstream instream(values_fname, std::ios::binary);
            if (!instream) {
                logger->error("Could not open file with column values {}", values_fname);
                exit(1);
            }

            values[i].resize(sources[i].num_labels());
            for (size_t j = 0; j < values[i].size(); ++j) {
                try {
                    values[i][j].load(instream);
                } catch (...) {
                    logger->error("Couldn't read column values from {}", values_fname);
                    exit(1);
                }
            }
        }
        logger->trace("Done loading column values");
    }

    anchor_bv_type anchor;
    if (!compute_row_reduction) {
        const std::string anchors_fname = pred_succ_fprefix + kRowDiffAnchorExt;
        std::ifstream f(anchors_fname, std::ios::binary);
        if (!anchor.load(f)) {
            logger->error("Can't load anchors from {}", anchors_fname);
            exit(1);
        }
        if (anchor.size() != num_rows) {
            logger->error("Anchor vector {} is incompatible with annotations."
                          " Vector size: {}, number of rows: {}",
                          anchors_fname, anchor.size(), num_rows);
            exit(1);
        }
    }

    const fs::path tmp_path = utils::create_temp_dir(swap_dir, "col");

    // stores the row indices that were set because of differences to incoming/outgoing
    // edges, for each of the sources, per chunk. set_rows_fwd is already sorted
    std::vector<std::vector<std::vector<T>>> set_rows_bwd(sources.size());
    std::vector<std::vector<std::vector<T>>> set_rows_fwd(sources.size());
    std::vector<std::vector<uint64_t>> row_diff_bits(sources.size());
    std::vector<std::vector<uint64_t>> num_relations_anchored(sources.size());
    std::vector<std::vector<uint64_t>> num_chunks(sources.size());

    auto tmp_file = [&](size_t s, size_t j, size_t chunk) {
        return tmp_path/fmt::format("{}/col_{}_{}/chunk_{}", s / 100, s, j, chunk);
    };
    auto dump_chunk_to_disk = [&](const std::vector<T> &v,
                                  size_t s, size_t j, size_t chunk) {
        assert(std::is_sorted(v.begin(), v.end()) && "all bits in chunks must be sorted");
        Encoder<T>::append_block(v, tmp_file(s, j, chunk));
        row_diff_bits[s][j] += v.size();
    };

    // In the first stage, only one buffer is created per column (`bwd`).
    // In the last stage, two buffers (`fwd` and `bwd`) are created per column.
    const uint64_t buf_size = compute_row_reduction
                                ? buf_size_bytes / sizeof(T)
                                : buf_size_bytes / sizeof(T) / 2;

    #pragma omp parallel for num_threads(num_threads)
    for (size_t s = 0; s < sources.size(); ++s) {
        set_rows_fwd[s].resize(sources[s].num_labels());
        set_rows_bwd[s].resize(sources[s].num_labels());
        row_diff_bits[s].assign(sources[s].num_labels(), 0);
        num_relations_anchored[s].assign(sources[s].num_labels(), 0);
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
    sdsl::int_vector_buffer<> row_reduction;
    const bool new_reduction_vector = !fs::exists(row_reduction_fname);
    if (compute_row_reduction) {
        logger->trace("Row reduction vector: {}", row_reduction_fname);
        if (new_reduction_vector) {
            // create an empty vector
            sdsl::int_vector_buffer<>(row_reduction_fname,
                                      std::ios::out, BUFFER_SIZE, ROW_REDUCTION_WIDTH);
            logger->trace("Initialized new row reduction vector");
        } else {
            logger->trace("Row reduction vector already exists and will be updated");
        }

        row_reduction = sdsl::int_vector_buffer<>(row_reduction_fname,
                                                  std::ios::in | std::ios::out, BUFFER_SIZE);

        if (!new_reduction_vector && row_reduction.size() != num_rows) {
            logger->error("Count vector {} is incompatible with annotations."
                          " Vector size: {}, number of rows: {}",
                          row_reduction_fname, row_reduction.size(), num_rows);
            exit(1);
        }
    }

    // total number of set bits in the original rows
    std::vector<uint32_t> row_nbits_block;
    // buffer for writing previous block while populating next row_nbits_block
    std::vector<uint32_t> row_nbits_block_other;

    // get bit at position |i| or its value
    auto get_value = [&](const bit_vector &col,
                         size_t s, size_t j, uint64_t i) -> uint64_t {
        if (!with_values) {
            return col[i];
        } else {
            if (uint64_t rk = col.conditional_rank1(i)) {
                return values[s][j][rk - 1];
            } else {
                return 0;
            }
        }
    };

    traverse_anno_chunked(
            num_rows, pred_succ_fprefix, sources,
            [&](uint64_t chunk_size) {
                if (compute_row_reduction)
                    row_nbits_block.assign(chunk_size, 0);
            },
            [&](const bit_vector &source_col, uint64_t row_idx, uint64_t chunk_idx,
                    size_t source_idx, size_t j,
                    const uint64_t *succ,
                    const uint64_t *pred_begin, const uint64_t *pred_end) {

                // get bits for these positions (or values, hence uint64_t)
                uint64_t curr_value;
                if constexpr(with_values) {
                    curr_value = values[source_idx][j][source_col.rank1(row_idx) - 1];
                } else {
                    curr_value = 1;
                }

                if (compute_row_reduction) {
                    if (succ && curr_value == get_value(source_col, source_idx, j, *succ)) {
                        // reduction (zero diff)
                        __atomic_add_fetch(&row_nbits_block[chunk_idx], 1, __ATOMIC_RELAXED);
                    }
                } else {
                    bool is_anchor = anchor[row_idx];
                    // add current bit if this node is an anchor
                    // or if the successor has zero diff
                    uint64_t succ_value = is_anchor ? 0 : get_value(source_col, source_idx, j, *succ);
                    if (succ_value != curr_value) {
                        // no reduction, we must keep the bit
                        auto &v = set_rows_fwd[source_idx][j];
                        if constexpr(with_values) {
                            v.emplace_back(row_idx, curr_value - succ_value);
                        } else {
                            v.push_back(row_idx);
                        }
                        if (is_anchor)
                            num_relations_anchored[source_idx][j]++;

                        if (v.size() == v.capacity()) {
                            // dump chunk to disk
                            dump_chunk_to_disk(v, source_idx, j, 0);
                            v.resize(0);
                        }
                    }
                }

                // check non-anchor predecessor nodes and add them if they are zero
                for (const uint64_t *pred_p = pred_begin; pred_p < pred_end; ++pred_p) {
                    if (curr_value && !source_col[*pred_p] && (compute_row_reduction || !anchor[*pred_p])) {
                        auto &v = set_rows_bwd[source_idx][j];
                        if constexpr(with_values) {
                            v.emplace_back(*pred_p, -curr_value);
                        } else {
                            v.push_back(*pred_p);
                        }

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

            logger->trace("Number of relations for column {} reduced from {}"
                          " to {}, of them stored in anchors: {}",
                          sources[s].get_label_encoder().decode(j),
                          sources[s].get_matrix().data()[j]->num_set_bits(),
                          row_diff_bits[s][j], num_relations_anchored[s][j]);
        }
    }

    async_writer.join();

    set_rows_fwd.clear(); // free up memory
    set_rows_bwd.clear();
    anchor = anchor_bv_type();

    std::vector<LabelEncoder<std::string>> label_encoders;
    for (const auto &source : sources) {
        label_encoders.push_back(source.get_label_encoder());
    }

    // free memory occupied by original columns
    sources.clear();

    std::vector<std::vector<std::unique_ptr<bit_vector>>> diff_columns(label_encoders.size());

    const uint32_t chunks_open_per_thread
            = MAX_NUM_FILES_OPEN / std::max((uint32_t)1, num_threads) / (2 + with_values);
    if (chunks_open_per_thread < 3) {
        logger->error("Can't merge with less than 3 chunks per thread open. "
                      "Max num files open: {}. Current number of threads: {}.",
                      MAX_NUM_FILES_OPEN, num_threads);
        exit(1);
    }

    logger->trace("Generating row_diff columns...");
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint32_t l_idx = 0; l_idx < label_encoders.size(); ++l_idx) {
        std::vector<std::unique_ptr<bit_vector>> columns(label_encoders[l_idx].size());

        for (size_t j = 0; j < label_encoders[l_idx].size(); ++j) {
            if constexpr(with_values) {
                // diff values may be negative, hence we need wider integers
                values[l_idx][j] = sdsl::int_vector<>(row_diff_bits[l_idx][j], 0,
                                                      std::min(values[l_idx][j].width() + 1, 64));
            }

            auto call_ones = [&](const std::function<void(uint64_t)> &call) {
                std::vector<std::string> filenames;
                // skip chunk with fwd bits which have already been counted if stage 1
                for (uint32_t chunk = compute_row_reduction ? 1 : 0;
                                    chunk < num_chunks[l_idx][j]; ++chunk) {
                    filenames.push_back(tmp_file(l_idx, j, chunk));
                }

                const bool remove_chunks = true;
                uint64_t r = 0;
                elias_fano::merge_files<T>(filenames, [&](T v) {
                    call(utils::get_first(v));
                    if constexpr(with_values) {
                        assert(v.second && "zero diffs must have been skipped");
                        values[l_idx][j][r++] = matrix::encode_diff(v.second);
                    }
                }, remove_chunks, chunks_open_per_thread);
            };
            columns[j] = std::make_unique<bit_vector_smart>(call_ones, num_rows,
                                                            row_diff_bits[l_idx][j]);
        }

        if (compute_row_reduction) {
            diff_columns[l_idx] = std::move(columns);

        } else {
            if (!columns.size())
                continue;

            auto fpath = col_out_dir/fs::path(source_files[l_idx]).filename();
            if constexpr(with_values) {
                fpath.replace_extension().replace_extension(ColumnCompressed<>::kExtension);
                ColumnCompressed<>(std::move(columns), label_encoders[l_idx]).serialize(fpath);
                logger->trace("Serialized {}", fpath);

                fpath.replace_extension().replace_extension(ColumnCompressed<>::kCountExtension);
                std::ofstream outstream(fpath, std::ios::binary);
                for (size_t j = 0; j < label_encoders[l_idx].size(); ++j) {
                    values[l_idx][j].serialize(outstream);
                }
                logger->trace("Serialized {}", fpath);
                values[l_idx].clear();

            } else {
                fpath.replace_extension().replace_extension(RowDiffColumnAnnotator::kExtension);
                RowDiffColumnAnnotator(
                        std::make_unique<RowDiff<ColumnMajor>>(nullptr, ColumnMajor(std::move(columns))),
                        std::move(label_encoders[l_idx])).serialize(fpath);
                logger->trace("Serialized {}", fpath);
            }
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
        for (size_t l_idx = 0; l_idx < diff_columns.size(); ++l_idx) {
            for (const auto &col_ptr : diff_columns[l_idx]) {
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

void convert_batch_to_row_diff_coord(const std::string &pred_succ_fprefix,
                                     const std::vector<std::string> &source_files,
                                     const fs::path &col_out_dir,
                                     const fs::path &swap_dir,
                                     const std::string &row_reduction_fname,
                                     uint64_t buf_size_bytes,
                                     bool compute_row_reduction,
                                     size_t num_coords_per_seq) {
    if (source_files.empty())
        return;

    const uint32_t num_threads = get_num_threads();
    using T = std::pair<uint64_t, uint64_t>;

    uint64_t num_rows;
    std::vector<annot::ColumnCompressed<>> sources = load_columns(source_files, &num_rows);
    if (!num_rows) {
        logger->warn("Input annotations have no rows");
        return;
    }

    std::vector<std::vector<sdsl::int_vector<>>> coords(source_files.size());
    std::vector<std::vector<bit_vector_smart>> delims(source_files.size());
    load_coordinates(source_files, sources,
        [&](size_t s, sdsl::int_vector<>&& c, bit_vector_smart&& d) {
            coords[s].push_back(std::move(c));
            delims[s].push_back(std::move(d));
        },
        num_threads
    );
    logger->trace("Done loading coordinates");

    // get bit at position |i| or its value
    auto get_value = [&](const bit_vector &col,
                         size_t s, size_t j, uint64_t i) {
        std::vector<uint64_t> result;
        if (uint64_t rk = col.conditional_rank1(i)) {
            uint64_t t = delims[s][j].select1(rk);
            while (!delims[s][j][++t]) {
                result.push_back(coords[s][j][t - rk]);
            }
        }
        std::sort(result.begin(), result.end());
        return result;
    };

    auto get_diff = [&](std::vector<uint64_t> curr,
                        std::vector<uint64_t>&& next) {
        if (num_coords_per_seq > 1) {
            // skip starting coordinates -- they can always be reconstructed
            // if all sequences have the same length
            size_t i = 0;
            for (size_t j = 0; j < next.size(); ++j) {
                if (next[j] % num_coords_per_seq)
                    next[i++] = next[j];
            }
            next.resize(i);
        }
        assert(std::is_sorted(curr.begin(), curr.end()));
        assert(std::is_sorted(next.begin(), next.end()));
        for (uint64_t &c : curr) {
            c++;
        }
        std::vector<uint64_t> diff;
        diff.reserve(curr.size() + next.size());
        std::set_symmetric_difference(curr.begin(), curr.end(),
                                      next.begin(), next.end(),
                                      std::back_inserter(diff));
        return diff;
    };

    if (compute_row_reduction) {
        ThreadPool async_writer(1, 1);
        sdsl::int_vector_buffer<> row_reduction;
        const bool new_reduction_vector = !fs::exists(row_reduction_fname);
        logger->trace("Row reduction vector: {}", row_reduction_fname);
        if (new_reduction_vector) {
            // create an empty vector
            sdsl::int_vector_buffer<>(row_reduction_fname,
                                      std::ios::out, BUFFER_SIZE, ROW_REDUCTION_WIDTH);
            logger->trace("Initialized new row reduction vector");
        } else {
            logger->trace("Row reduction vector already exists and will be updated");
        }

        row_reduction = sdsl::int_vector_buffer<>(row_reduction_fname,
                                                  std::ios::in | std::ios::out, BUFFER_SIZE);

        if (!new_reduction_vector && row_reduction.size() != num_rows) {
            logger->error("Count vector {} is incompatible with annotations."
                          " Vector size: {}, number of rows: {}",
                          row_reduction_fname, row_reduction.size(), num_rows);
            exit(1);
        }

        // total number of set bits in the original rows
        std::vector<uint32_t> row_nbits_block;
        // buffer for writing previous block while populating next row_nbits_block
        std::vector<uint32_t> row_nbits_block_other;

        traverse_anno_chunked(
                num_rows, pred_succ_fprefix, sources,
                [&](uint64_t chunk_size) {
                    row_nbits_block.assign(chunk_size, 0);
                },
                [&](const bit_vector &source_col, uint64_t row_idx, uint64_t chunk_idx,
                        size_t s, size_t j,
                        const uint64_t *succ, const uint64_t *, const uint64_t *) {
                    if (!succ)
                        return;

                    // get annotated coordinates for this k-mer
                    const auto curr_value = get_value(source_col, s, j, row_idx);
                    const auto diff = get_diff(curr_value, get_value(source_col, s, j, *succ));
                    // reduction (zero diff)
                    __atomic_add_fetch(&row_nbits_block[chunk_idx],
                                       curr_value.size() - diff.size(),
                                       __ATOMIC_RELAXED);
                },
                [&](uint64_t block_begin) {
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
                }
        );

        async_writer.join();
        return;
    }

    anchor_bv_type anchor;
    const std::string anchors_fname = pred_succ_fprefix + kRowDiffAnchorExt;
    std::ifstream f(anchors_fname, std::ios::binary);
    if (!anchor.load(f)) {
        logger->error("Can't load anchors from {}", anchors_fname);
        exit(1);
    }
    if (anchor.size() != num_rows) {
        logger->error("Anchor vector {} is incompatible with annotations."
                      " Vector size: {}, number of rows: {}",
                      anchors_fname, anchor.size(), num_rows);
        exit(1);
    }

    const fs::path tmp_path = utils::create_temp_dir(swap_dir, "col");

    // stores the row indices that were set because of differences to incoming/outgoing
    // edges, for each of the sources, per chunk. set_rows_fwd is already sorted
    std::vector<std::vector<std::vector<T>>> set_rows_bwd(sources.size());
    std::vector<std::vector<std::vector<T>>> set_rows_fwd(sources.size());
    std::vector<std::vector<uint64_t>> row_diff_bits(sources.size());
    std::vector<std::vector<uint64_t>> row_diff_coords(sources.size());
    std::vector<std::vector<uint64_t>> num_coords_anchored(sources.size());
    std::vector<std::vector<uint64_t>> num_chunks(sources.size());

    auto tmp_file = [&](size_t s, size_t j, size_t chunk) {
        return tmp_path/fmt::format("{}/col_{}_{}/chunk_{}", s / 100, s, j, chunk);
    };
    auto dump_chunk_to_disk = [&](const std::vector<T> &v,
                                  size_t s, size_t j, size_t chunk) {
        assert(std::is_sorted(v.begin(), v.end()) && "all bits in chunks must be sorted");
        Encoder<T>::append_block(v, tmp_file(s, j, chunk));
    };

    // two buffers (`fwd` and `bwd`) are created per column
    const uint64_t buf_size = buf_size_bytes / sizeof(T) / 2;

    #pragma omp parallel for num_threads(num_threads)
    for (size_t s = 0; s < sources.size(); ++s) {
        set_rows_fwd[s].resize(sources[s].num_labels());
        set_rows_bwd[s].resize(sources[s].num_labels());
        row_diff_bits[s].assign(sources[s].num_labels(), 0);
        row_diff_coords[s].assign(sources[s].num_labels(), 0);
        num_coords_anchored[s].assign(sources[s].num_labels(), 0);
        // The first chunk will contain forward bits, all sorted.
        // The other ones (added later) will contain chunks of sorted pred bits.
        num_chunks[s].assign(sources[s].num_labels(), 1);

        for (size_t j = 0; j < sources[s].num_labels(); ++j) {
            fs::create_directories(tmp_file(s, j, 0).parent_path());
            uint64_t original_nbits = sources[s].get_matrix().data()[j]->num_set_bits();
            set_rows_bwd[s][j].reserve(std::min(buf_size, original_nbits));

            // make sure the first chunk exists even if empty
            dump_chunk_to_disk({}, s, j, 0);
            set_rows_fwd[s][j].reserve(std::min(buf_size, original_nbits));
        }
    }

    // We use this dummy index for an optimization where we don't store diff
    // for non-anchor k-mers with no coordinates.
    uint64_t DUMMY_COORD = -1;

    traverse_anno_chunked(
            num_rows, pred_succ_fprefix, sources,
            [&](uint64_t) {},
            [&](const bit_vector &source_col, uint64_t row_idx, uint64_t,
                    size_t s, size_t j,
                    const uint64_t *succ,
                    const uint64_t *pred_begin, const uint64_t *pred_end) {
                // get annotated coordinates for this k-mer
                const auto curr_value = get_value(source_col, s, j, row_idx);

                bool is_anchor = anchor[row_idx];

                if (is_anchor)
                    num_coords_anchored[s][j] += curr_value.size();

                const auto diff = is_anchor
                    ? curr_value
                    : get_diff(curr_value, get_value(source_col, s, j, *succ));

                if (diff.size()) {
                    // must write the coordinates/diff
                    auto &v = set_rows_fwd[s][j];

                    for (uint64_t coord : diff) {
                        assert((!v.size() || v.back() != std::make_pair(row_idx, coord))
                               && "coordinates must be unique and can't repeat");
                        v.emplace_back(row_idx, coord);
                        row_diff_coords[s][j]++;

                        if (v.size() == v.capacity()) {
                            // dump chunk to disk
                            dump_chunk_to_disk(v, s, j, 0);
                            v.resize(0);
                        }
                    }
                    row_diff_bits[s][j]++;
                }

                // check non-anchor predecessor nodes and add them if they are zero
                for (const uint64_t *pred_p = pred_begin; pred_p < pred_end; ++pred_p) {
                    if (curr_value.size() && !source_col[*pred_p] && !anchor[*pred_p]) {
                        auto &v = set_rows_bwd[s][j];
                        // indicate that there are no coordinates for the predecessor
                        v.emplace_back(*pred_p, DUMMY_COORD);
                        row_diff_bits[s][j]++;

                        if (v.size() == v.capacity()) {
                            std::sort(v.begin(), v.end());
                            dump_chunk_to_disk(v, s, j, num_chunks[s][j]++);
                            v.resize(0);
                        }
                    }
                }
            },
            [&](uint64_t) {}
    );

    #pragma omp parallel for num_threads(num_threads)
    for (size_t s = 0; s < sources.size(); ++s) {
        for (size_t j = 0; j < sources[s].num_labels(); ++j) {
            auto &fwd = set_rows_fwd[s][j];
            if (fwd.size())
                dump_chunk_to_disk(fwd, s, j, 0);

            auto &bwd = set_rows_bwd[s][j];
            if (bwd.size()) {
                std::sort(bwd.begin(), bwd.end());
                dump_chunk_to_disk(bwd, s, j, num_chunks[s][j]++);
            }
            logger->trace("Number of coordinates for column {} reduced from {}"
                          " to {}, number of coordinates stored in anchors: {}",
                          sources[s].get_label_encoder().decode(j), coords[s][j].size(),
                          row_diff_coords[s][j], num_coords_anchored[s][j]);
        }
    }

    set_rows_fwd.clear(); // free up memory
    set_rows_bwd.clear();
    anchor = anchor_bv_type();

    std::vector<LabelEncoder<std::string>> label_encoders;
    for (const auto &source : sources) {
        label_encoders.push_back(source.get_label_encoder());
    }

    // free memory occupied by original columns
    sources.clear();
    coords.clear();
    delims.clear();

    std::vector<std::vector<std::unique_ptr<bit_vector>>> diff_columns(label_encoders.size());

    const uint32_t chunks_open_per_thread
            = MAX_NUM_FILES_OPEN / std::max((uint32_t)1, num_threads) / 3;
    if (chunks_open_per_thread < 3) {
        logger->error("Can't merge with less than 3 chunks per thread open. "
                      "Max num files open: {}. Current number of threads: {}.",
                      MAX_NUM_FILES_OPEN, num_threads);
        exit(1);
    }

    logger->trace("Generating row_diff columns...");
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (uint32_t l_idx = 0; l_idx < label_encoders.size(); ++l_idx) {
        if (!label_encoders[l_idx].size())
            continue;

        std::vector<std::unique_ptr<bit_vector>> columns(label_encoders[l_idx].size());

        auto fpath_coord = col_out_dir/fs::path(source_files[l_idx]).filename();
        fpath_coord.replace_extension().replace_extension(ColumnCompressed<>::kCoordExtension);
        std::ofstream out_coord(fpath_coord, std::ios::binary);

        for (size_t j = 0; j < label_encoders[l_idx].size(); ++j) {
            // diff values may be negative, hence we need wider integers
            sdsl::int_vector<> diff_coords(row_diff_coords[l_idx][j]);
            auto coords_it = diff_coords.begin();
            sdsl::bit_vector diff_delims(diff_coords.size() + row_diff_bits[l_idx][j] + 1, 0);
            auto delims_it = diff_delims.begin();

            auto call_ones = [&](const std::function<void(uint64_t)> &call) {
                std::vector<std::string> filenames;
                // skip chunk with fwd bits which have already been counted if stage 1
                for (uint32_t chunk = 0; chunk < num_chunks[l_idx][j]; ++chunk) {
                    filenames.push_back(tmp_file(l_idx, j, chunk));
                }

                uint64_t last = -1;
                elias_fano::merge_files<T>(filenames, [&](T pair) {
                    const auto &[i, coord] = pair;
                    if (i != last) {
                        call(i);
                        last = i;
                        *delims_it++ = 1;
                    }
                    if (coord != DUMMY_COORD) {
                        *coords_it++ = coord;
                        ++delims_it;
                    }
                }, true, chunks_open_per_thread);
                *delims_it++ = 1;
                assert(coords_it == diff_coords.end());
                assert(delims_it == diff_delims.end());
            };
            columns[j] = std::make_unique<bit_vector_smart>(call_ones, num_rows,
                                                            row_diff_bits[l_idx][j]);
            bit_vector_smart(std::move(diff_delims)).serialize(out_coord);
            sdsl::util::bit_compress(diff_coords);
            diff_coords.serialize(out_coord);
        }
        logger->trace("Serialized {}", fpath_coord);
        auto fpath = col_out_dir/fs::path(source_files[l_idx]).filename();
        fpath.replace_extension().replace_extension(ColumnCompressed<>::kExtension);
        ColumnCompressed<>(std::move(columns), label_encoders[l_idx]).serialize(fpath);
        logger->trace("Serialized {}", fpath);
    }

    utils::remove_temp_dir(tmp_path);
}

} // namespace annot
} // namespace mtg
