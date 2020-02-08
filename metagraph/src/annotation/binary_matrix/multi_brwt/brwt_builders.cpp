#include "brwt_builders.hpp"

#include <omp.h>
#include <progress_bar.hpp>

#include "common/algorithms.hpp"
#include "common/logger.hpp"
#include "common/vectors/int_vector_algorithm.hpp"

const uint64_t kBlockSize = 1'000'000 / 64 * 64;

// Each block is a multiple of 64 bits for thread safety
static_assert((kBlockSize % 64) == 0);


BRWTBottomUpBuilder::Partitioner
BRWTBottomUpBuilder::get_basic_partitioner(size_t arity) {
    assert(arity > 1u);

    return [arity](const VectorPtrs &vectors) {
        if (!vectors.size())
            return Partition(0);

        Partition partition((vectors.size() - 1) / arity + 1);
        for (size_t i = 0; i < vectors.size(); ++i) {
            partition[i / arity].push_back(i);
        }
        return partition;
    };
}

// check if the bit vector is sparser than the threshold for its type
bool is_sparser(const bit_vector &bv,
                double threshold_stat,
                double threshold_sd,
                double threshold_rrr,
                double threshold_other) {
    if (bv.size() < 64u)
        return false;

    double density = static_cast<double>(bv.num_set_bits()) / bv.size();

    if (auto *adaptive = dynamic_cast<const bit_vector_adaptive *>(&bv)) {
        switch (adaptive->representation_tag()) {
            case bit_vector_adaptive::STAT_VECTOR:
                return density < threshold_stat;
            case bit_vector_adaptive::SD_VECTOR:
                return density < threshold_sd;
            case bit_vector_adaptive::RRR_VECTOR:
                return density < threshold_rrr;
        }
    }

    if (dynamic_cast<const bit_vector_stat *>(&bv))
        return density < threshold_stat;

    if (dynamic_cast<const bit_vector_sd *>(&bv))
        return density < threshold_sd;

    return density < threshold_other;
}

// TODO: move to int_vector_algorithm.hpp
void compute_or(const std::vector<const bit_vector *> &columns,
                sdsl::bit_vector *result,
                ThreadPool &thread_pool) {
    uint64_t size = columns.at(0)->size();

    assert(result);
    assert(result->size() == size);

    const uint64_t block_size = std::max(kBlockSize, size / 100 / 64 * 64);

    // Each block is a multiple of 64 bits for thread safety
    assert(!(block_size & 0x3F));

    std::vector<std::future<void>> results;

    for (uint64_t i = 0; i < size; i += block_size) {
        results.push_back(thread_pool.enqueue([&](uint64_t begin, uint64_t end) {
            std::fill(result->data() + (begin >> 6),
                      result->data() + ((end + 63) >> 6),
                      0);

            for (const auto &col_ptr : columns) {
                assert(col_ptr);
                assert(col_ptr->size() == size);

                if (is_sparser(*col_ptr, 0, 0.2, 0.01, 0.01)) {
                    col_ptr->call_ones_in_range(begin, end,
                        [result](uint64_t k) { (*result)[k] = true; }
                    );
                } else {
                    assert((begin & 0x3F) == 0);

                    uint64_t i;
                    for (i = begin; i + 64 <= end; i += 64) {
                        result->data()[i / 64] |= col_ptr->get_int(i, 64);
                    }
                    if (i < end)
                        result->data()[i / 64] |= col_ptr->get_int(i, end - i);
                }
            }

        }, i, std::min(i + block_size, size)));
    }

    std::for_each(results.begin(), results.end(), [](auto &res) { res.wait(); });
}

void compute_subindex(const bit_vector &column,
                      const sdsl::bit_vector &reference,
                      uint64_t begin, uint64_t end,
                      uint64_t reference_rank_begin,
                      sdsl::bit_vector *subindex);

sdsl::bit_vector generate_subindex(const bit_vector &column,
                                   const sdsl::bit_vector &reference,
                                   uint64_t reference_num_set_bits,
                                   ThreadPool &thread_pool) {
    assert(column.size() == reference.size());

    sdsl::select_support_scan_offset<> reference_select(&reference);

    // no shrinkage if vectors are the same
    if (column.num_set_bits() == reference_num_set_bits)
        return sdsl::bit_vector(reference_num_set_bits, true);

    sdsl::bit_vector subindex(reference_num_set_bits, false);

    const uint64_t block_size
      = std::max(reference_num_set_bits / 100,
                 static_cast<uint64_t>(1'000)) & ~0x3Full;
    // Each block is a multiple of 64 bits for thread safety
    assert(!(block_size & 0x3F));

    std::vector<std::future<void>> futures;

    uint64_t block_end = 0;

    for (uint64_t offset = 0; offset < reference_num_set_bits; offset += block_size) {
        // ........... [1     1       1    1  1] ............. //
        //              ^                     ^                //
        //          (rank = 1 mod 64)     (rank = 0 mod 64)    //
        //              |                     |                //
        //              |_______ BLOCK _______|                //

        uint64_t block_begin = block_end;

        if (offset + block_size >= reference_num_set_bits) {
            block_end = reference.size();
        } else {
            block_end = reference_select.select_offset(block_size + 1, block_begin);
            assert(count_ones(reference, block_begin, block_end) == block_size);
        }

        futures.push_back(thread_pool.enqueue(
            [&](uint64_t begin, uint64_t end, uint64_t offset_rank) {
                compute_subindex(column, reference, begin, end, offset_rank, &subindex);
            },
            block_begin, block_end, offset
        ));
    }

    std::for_each(futures.begin(), futures.end(), [](auto &f) { f.wait(); });

    return subindex;
}

sdsl::bit_vector generate_subindex(const bit_vector &column,
                                   const bit_vector_stat &reference) {
    assert(column.size() == reference.size());

    uint64_t shrinked_size = reference.num_set_bits();

    // no shrinkage if vectors are the same
    if (column.num_set_bits() == shrinked_size)
        return sdsl::bit_vector(shrinked_size, true);

    sdsl::bit_vector subindex(shrinked_size, false);

    compute_subindex(column, reference.data(), 0, reference.size(), 0, &subindex);

    return subindex;
}

void compute_subindex(const bit_vector &column,
                      const sdsl::bit_vector &reference,
                      uint64_t begin, uint64_t end,
                      uint64_t offset,
                      sdsl::bit_vector *subindex) {
    assert(column.size() == reference.size());
    assert(begin <= end);
    assert(end <= reference.size());

    if (begin == end)
        return;

    uint64_t popcount = column.rank1(end - 1)
                        - (begin ? column.rank1(begin - 1) : 0);

    // check if all zeros
    if (!popcount)
        return;

    uint64_t i = begin;
    uint64_t rank = offset;

    // check if all ones
    if (popcount == end - begin) {
        for ( ; i + 64 <= end; i += 64, rank += 64) {
            subindex->set_int(rank, 0xFFFF, 64);
        }
        if (begin < end) {
            subindex->set_int(rank, sdsl::bits::lo_set[end - begin], end - begin);
        }
        return;
    }


#if __BMI2__
    if (is_sparser(column, 0, 0.02, 0.01, 0.01)) {
#endif
        // sparse or without __BMI2__
        column.call_ones_in_range(begin, end, [&](uint64_t j) {
            assert(j >= i);
            assert(reference[j]);

            rank += count_ones(reference, i, j);

            assert(rank < subindex->size());

            (*subindex)[rank++] = true;
            i = j + 1;
        });
#if __BMI2__
    } else {
        // dense
        for ( ; i + 64 <= end; i += 64) {
            uint64_t mask = reference.get_int(i, 64);
            int popcount = sdsl::bits::cnt(mask);
            if (uint64_t a = column.get_int(i, 64))
                subindex->set_int(rank, _pext_u64(a, mask), popcount);
            rank += popcount;
        }
        if (i < end) {
            uint64_t mask = reference.get_int(i, end - i);
            int popcount = sdsl::bits::cnt(mask);
            if (uint64_t a = column.get_int(i, end - i))
                subindex->set_int(rank, _pext_u64(a, mask), popcount);
        }
    }
#endif
}

BRWT BRWTBottomUpBuilder::concatenate(std::vector<BRWT>&& submatrices,
                                      sdsl::bit_vector *buffer,
                                      ThreadPool &thread_pool) {
    assert(submatrices.size());

    if (submatrices.size() == 1)
        return std::move(submatrices[0]);

    BRWT parent;

    VectorPtrs index_columns(submatrices.size());
    for (size_t i = 0; i < submatrices.size(); ++i) {
        index_columns[i] = submatrices[i].nonzero_rows_.get();
    }

    // build an aggregated parent index column
    compute_or(index_columns, buffer, thread_pool);

    // set column assignments
    uint64_t num_columns = 0;
    Partition partition;
    for (const BRWT &submatrix : submatrices) {
        partition.push_back(utils::arange(num_columns, submatrix.num_columns()));
        num_columns += submatrix.num_columns();
    }
    parent.assignments_ = RangePartition(std::move(partition));

    // set child nodes
    parent.child_nodes_.resize(submatrices.size());

    std::vector<std::future<void>> results;

    // compress the parent index vector
    results.push_back(thread_pool.enqueue([&]() {
        parent.nonzero_rows_.reset(new bit_vector_smart(*buffer));
    }));

    uint64_t subindex_size = sdsl::util::cnt_one_bits(*buffer);

    for (size_t i = 0; i < submatrices.size(); ++i) {
        // generate an index column for the child node
        sdsl::bit_vector subindex
            = generate_subindex(*submatrices[i].nonzero_rows_, *buffer,
                                subindex_size, thread_pool);

        // the full index vector is not needed anymore, the subindex
        // will be used instead
        submatrices[i].nonzero_rows_.reset();

        results.push_back(thread_pool.enqueue(
            [&,i,subindex(std::move(subindex))]() {
                // compress the subindex vector and set it to the child node
                // all in a single thread
                submatrices[i].nonzero_rows_
                    = std::make_unique<bit_vector_small>(std::move(subindex));
                parent.child_nodes_[i].reset(new BRWT(std::move(submatrices[i])));
            }
        ));
    }

    std::for_each(results.begin(), results.end(), [](auto &res) { res.wait(); });

    return parent;
}

template <typename T>
std::vector<T> subset(std::vector<T> *vector,
                      const std::vector<uint64_t> indexes) {
    assert(vector);

    std::vector<T> result;
    result.reserve(indexes.size());
    for (auto j : indexes) {
        result.push_back(std::move((*vector)[j]));
    }
    return result;
}

BRWT BRWTBottomUpBuilder::build(std::vector<std::unique_ptr<bit_vector>>&& columns,
                                Partitioner partitioner,
                                size_t num_nodes_parallel,
                                size_t num_threads) {
    if (!columns.size())
        return BRWT();

    std::vector<BRWT> nodes(columns.size());

    // construct leaves
    for (size_t i = 0; i < columns.size(); ++i) {
        nodes[i].assignments_ = RangePartition(Partition(1, { 0 }));
        nodes[i].nonzero_rows_ = std::move(columns[i]);
    }

    return merge(std::move(nodes), partitioner, num_nodes_parallel, num_threads);
}

BRWT BRWTBottomUpBuilder::merge(std::vector<BRWT>&& nodes,
                                Partitioner partitioner,
                                size_t num_nodes_parallel,
                                size_t num_threads) {
    if (!nodes.size())
        return BRWT();

    num_nodes_parallel = std::min(num_nodes_parallel, nodes.size());

    num_threads = std::max(num_nodes_parallel, num_threads);

    ThreadPool thread_pool(num_threads, 100'000 * num_threads);

    // initialize buffers for merging columns
    // these may be huge, so we keep only a limited number of them
    std::vector<sdsl::bit_vector> buffers;
    for (size_t i = 0; i < std::max(num_nodes_parallel, size_t(1)); ++i) {
        buffers.emplace_back(nodes.at(0).nonzero_rows_->size());
    }

    uint64_t num_columns_total = 0;
    Partition current_partition;
    for (const BRWT &node : nodes) {
        current_partition.push_back(utils::arange(num_columns_total, node.num_columns()));
        num_columns_total += node.num_columns();
    }

    for (size_t level = 1; nodes.size() > 1; ++level) {
        mg::common::logger->trace("BRWT construction: level {}", level);

        VectorPtrs columns(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) {
            columns[i] = nodes[i].nonzero_rows_.get();
        }
        auto groups = partitioner(columns);

        assert(groups.size() > 0);
        assert(groups.size() < nodes.size());

        std::vector<BRWT> parent_nodes(groups.size());

        Partition partition = std::move(current_partition);
        current_partition.assign(groups.size(), {});

        ProgressBar progress_bar(groups.size(), "Building BRWT level",
                                 std::cerr, !utils::get_verbose());

        #pragma omp parallel for num_threads(num_nodes_parallel) schedule(dynamic)
        for (size_t g = 0; g < groups.size(); ++g) {
            // compute columns assignments for the parent
            for (size_t i : groups[g]) {
                std::copy(partition[i].begin(), partition[i].end(),
                          std::back_inserter(current_partition[g]));
            }
            // merge submatrices
            parent_nodes[g] = concatenate(subset(&nodes, groups[g]),
                                          &buffers.at(omp_get_thread_num()),
                                          thread_pool);
            ++progress_bar;
        }

        nodes = std::move(parent_nodes);
    }

    // All submatrices must be merged into one
    assert(nodes.size() == 1u);
    assert(current_partition.size() == 1u);

    // get the root node in Multi-BRWT
    BRWT root = std::move(nodes.at(0));
    // compress the index vector
    root.nonzero_rows_ = std::make_unique<bit_vector_small>(
        root.nonzero_rows_->convert_to<bit_vector_small>()
    );
    // update the column arrangement to be consistent with the initial
    // order 1,2,...,m
    const auto &column_arrangement = current_partition.at(0);
    std::vector<size_t> submatrix_sizes;
    for (size_t g = 0; g < root.assignments_.num_groups(); ++g) {
        submatrix_sizes.push_back(root.assignments_.group_size(g));
    }
    root.assignments_ = RangePartition(column_arrangement, submatrix_sizes);

    return root;
}


void BRWTOptimizer::relax(BRWT *brwt_matrix, uint64_t max_arity, size_t num_threads) {
    assert(brwt_matrix);

    std::deque<BRWT*> parents;

    // traverse the tree and collect the nodes
    // order: leaves first, root last
    brwt_matrix->BFT([&](const BRWT &node) {
        if (node.child_nodes_.size())
            parents.push_front(const_cast<BRWT*>(&node));
    });

    ProgressBar progress_bar(parents.size(), "Relax Multi-BRWT",
                             std::cerr, !utils::get_verbose());

    for (BRWT *parent : parents) {
        assert(parent->child_nodes_.size());

        for (int g = parent->child_nodes_.size() - 1; g >= 0; --g) {
            const auto *node = dynamic_cast<const BRWT*>(parent->child_nodes_[g].get());

            if (node && should_prune(*node)
                     && parent->child_nodes_.size() - 1
                            + node->child_nodes_.size() <= max_arity) {
                // remove this node and reassign all its children directly to its parent
                reassign(g, parent, num_threads);
            }
        }

        ++progress_bar;
    }
}

bool BRWTOptimizer::should_prune(const BRWT &node) {
    // we don't remove leaves in BRWT
    if (!node.child_nodes_.size())
        return false;

    // we relax only those child nodes, that have only BRWT nodes
    // because nodes of other types can't be reassigned
    for (const auto &child_node : node.child_nodes_) {
        if (!dynamic_cast<const BRWT *>(child_node.get()))
            return false;
    }

    // we don't remove nodes if it doesn't reduce the size
    return pruning_delta(node) <= 0;
}

void BRWTOptimizer::reassign(size_t node_rank, BRWT *parent, size_t num_threads) {
    assert(parent);

    BRWT &node = dynamic_cast<BRWT&>(*parent->child_nodes_.at(node_rank));

    std::vector<uint64_t> column_arrangement;
    std::vector<size_t> group_sizes;
    for (size_t g = 0; g < parent->assignments_.num_groups(); ++g) {
        if (g == node_rank)
            continue;  // skip the columns assigned to |node|

        size_t group_size = parent->child_nodes_[g]->num_columns();
        group_sizes.push_back(group_size);
        for (size_t r = 0; r < group_size; ++r) {
            column_arrangement.push_back(parent->assignments_.get(g, r));
        }
    }
    // reassign the node's columns
    for (size_t g = 0; g < node.assignments_.num_groups(); ++g) {
        size_t group_size = node.child_nodes_[g]->num_columns();
        group_sizes.push_back(group_size);
        for (size_t r = 0; r < group_size; ++r) {
            column_arrangement.push_back(
                parent->assignments_.get(node_rank, node.assignments_.get(g, r))
            );
        }
    }

    auto children_reassigned = std::move(node.child_nodes_);
    const auto index_column = node.nonzero_rows_->convert_to<sdsl::bit_vector>();

    // remove the |node| from the list of parent's children
    // after erasing, |node| will be deleted and thus must not be accessed
    parent->child_nodes_.erase(parent->child_nodes_.begin() + node_rank);
    // insert new children to |parent|
    parent->child_nodes_.resize(parent->child_nodes_.size()
                                    + children_reassigned.size());
    parent->assignments_ = RangePartition(column_arrangement, group_sizes);

    // the column assignments are updated
    // nonzero_rows_ stays unchanged in |parent|
    // now we need to insert the child nodes and update them accordingly
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t t = 0; t < children_reassigned.size(); ++t) {

        std::unique_ptr<BRWT> grand_child {
            dynamic_cast<BRWT*>(children_reassigned[t].release())
        };
        assert(grand_child);

        sdsl::bit_vector subindex(index_column.size(), false);

        uint64_t child_i = 0;
        uint64_t w = 0;
        call_ones(index_column, [&](auto i) {
            if (child_i % 64 == 0) {
                w = grand_child->nonzero_rows_->get_int(child_i,
                    std::min(static_cast<uint64_t>(64),
                             grand_child->nonzero_rows_->size() - child_i)
                );
            }

            if (w & 1)
                subindex[i] = true;

            w >>= 1;
            child_i++;
        });

        grand_child->nonzero_rows_
            = std::make_unique<bit_vector_small>(std::move(subindex));

        parent->child_nodes_[parent->child_nodes_.size()
                        - children_reassigned.size() + t] = std::move(grand_child);
    }
}

double BRWTOptimizer::pruning_delta(const BRWT &node) {
    assert(node.child_nodes_.size());

    double delta = 0;

    for (const auto &submatrix : node.child_nodes_) {
        auto *brwt_child = dynamic_cast<const BRWT*>(submatrix.get());
        assert(brwt_child);

        assert(brwt_child->num_rows() <= node.num_rows());
        assert(brwt_child->nonzero_rows_->size() <= node.num_rows());

        // updated vector
        delta += predict_size<bit_vector_small>(
                        node.num_rows(),
                        brwt_child->nonzero_rows_->num_set_bits());

        // old index vector
        delta -= predict_size<bit_vector_small>(
                        brwt_child->nonzero_rows_->size(),
                        brwt_child->nonzero_rows_->num_set_bits());
    }

    // removed index vector
    delta -= predict_size<bit_vector_small>(
                    node.nonzero_rows_->size(),
                    node.nonzero_rows_->num_set_bits());

    return delta;
}
