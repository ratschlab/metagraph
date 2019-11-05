#include "BRWT_builders.hpp"

#include <omp.h>
#include <progress_bar.hpp>

#include "utils.hpp"

const uint64_t kBlockSize = 1'000'000 / 64 * 64;

// Each block is a multiple of 64 bits for thread safety
static_assert((kBlockSize % 64) == 0);

BRWT
BRWTBuilder::initialize(NodeBRWT&& node, bit_vector&& nonzero_rows) {
    BRWT brwt;

    brwt.assignments_ = utils::RangePartition(node.column_arrangement,
                                              node.group_sizes);
    brwt.nonzero_rows_ = nonzero_rows.convert_to<bit_vector_small>();

    brwt.child_nodes_ = std::move(node.child_nodes);

    return brwt;
}

BRWTBottomUpBuilder::Partitioner
BRWTBottomUpBuilder::get_basic_partitioner(size_t arity) {
    assert(arity > 1u);

    return [arity](const VectorsPtr &vectors) {
        if (!vectors.size())
            return Partition(0);

        Partition partition((vectors.size() - 1) / arity + 1);
        for (size_t i = 0; i < vectors.size(); ++i) {
            partition[i / arity].push_back(i);
        }
        return partition;
    };
}

void compute_or(const std::vector<std::unique_ptr<bit_vector>> &columns,
                sdsl::bit_vector *buffer,
                ThreadPool &thread_pool) {
    uint64_t size = columns.at(0)->size();

    assert(buffer);
    assert(buffer->size() == size);

    auto &merged_result = *buffer;

    const uint64_t block_size = std::max(kBlockSize, size / 100 / 64 * 64);

    // Each block is a multiple of 64 bits for thread safety
    assert(!(block_size & 0x3F));

    std::vector<std::future<void>> results;

    for (uint64_t i = 0; i < size; i += block_size) {
        results.push_back(thread_pool.enqueue([&](uint64_t begin, uint64_t end) {
            std::fill(merged_result.data() + (begin >> 6),
                      merged_result.data() + ((end + 63) >> 6),
                      0);

            for (const auto &col_ptr : columns) {
                assert(col_ptr.get());
                assert(col_ptr->size() == size);
                col_ptr->call_ones_in_range(
                    begin, end,
                    [&merged_result](uint64_t k) { merged_result[k] = true; }
                );
            }

        }, i, std::min(i + block_size, size)));
    }

    std::for_each(results.begin(), results.end(), [](auto &res) { res.wait(); });
}

sdsl::bit_vector generate_subindex(const bit_vector &column,
                                   const sdsl::bit_vector &reference) {
    assert(column.size() == reference.size());

    uint64_t shrinked_size = sdsl::util::cnt_one_bits(reference);

    // no shrinkage if vectors are the same
    if (column.num_set_bits() == shrinked_size)
        return sdsl::bit_vector(shrinked_size, true);

    sdsl::bit_vector subindex(shrinked_size, false);

    uint64_t rank = 0;
    uint64_t offset = 0;

    column.call_ones([&](auto i) {
        assert(i >= offset);
        assert(reference[i]);

        rank += count_ones(reference, offset, i);

        assert(rank < subindex.size());

        subindex[rank++] = true;
        offset = i + 1;
    });

    return subindex;
}

std::pair<BRWTBottomUpBuilder::NodeBRWT, std::unique_ptr<bit_vector>>
BRWTBottomUpBuilder::merge(std::vector<NodeBRWT>&& nodes,
                           std::vector<std::unique_ptr<bit_vector>>&& index,
                           sdsl::bit_vector *buffer,
                           ThreadPool &thread_pool) {
    assert(nodes.size());
    assert(nodes.size() == index.size());

    if (nodes.size() == 1)
        return { std::move(nodes[0]), std::move(index[0]) };

    NodeBRWT parent;

    // build the parent aggregated column
    compute_or(index, buffer, thread_pool);

    for (size_t i = 0; i < nodes.size(); ++i) {
        // define column assignments for parent
        parent.column_arrangement.insert(parent.column_arrangement.end(),
                                         nodes[i].column_arrangement.begin(),
                                         nodes[i].column_arrangement.end());
        parent.group_sizes.push_back(nodes[i].column_arrangement.size());
    }

    // initialize child nodes
    parent.child_nodes.resize(nodes.size());

    std::vector<std::future<void>> results;

    for (size_t i = 0; i < nodes.size(); ++i) {
        results.push_back(thread_pool.enqueue([&](size_t i) {
            std::iota(nodes[i].column_arrangement.begin(),
                      nodes[i].column_arrangement.end(), 0);

            // shrink index column
            bit_vector_small shrinked_index(
                generate_subindex(std::move(*index[i]), *buffer)
            );
            index[i].reset();

            parent.child_nodes.at(i).reset(
                new BRWT(initialize(std::move(nodes[i]),
                                    std::move(shrinked_index)))
            );

        }, i));
    }

    std::unique_ptr<bit_vector_smart> parent_index_compressed;

    results.push_back(thread_pool.enqueue([&]() {
        parent_index_compressed.reset(new bit_vector_smart(*buffer));
    }));

    std::for_each(results.begin(), results.end(), [](auto &res) { res.wait(); });

    return { std::move(parent), std::move(parent_index_compressed) };
}

template <typename T>
std::vector<T> subset(std::vector<T> *vector,
                      const std::vector<BRWT::Column> indexes) {
    assert(vector);

    std::vector<T> result;
    result.reserve(indexes.size());
    for (auto j : indexes) {
        result.push_back(std::move((*vector)[j]));
    }
    return result;
}

BRWT BRWTBottomUpBuilder::build(VectorsPtr&& columns,
                                Partitioner partitioner,
                                size_t num_nodes_parallel,
                                size_t num_threads) {
    if (!columns.size())
        return BRWT();

    std::vector<NodeBRWT> nodes(columns.size());

    // construct leaves
    for (size_t i = 0; i < columns.size(); ++i) {
        nodes[i].column_arrangement = { i };
        nodes[i].group_sizes = { 1 };
    }

    return merge(std::move(columns),
                 std::move(nodes),
                 partitioner,
                 num_nodes_parallel,
                 num_threads);
}

BRWT BRWTBottomUpBuilder::merge(VectorsPtr&& columns,
                                std::vector<NodeBRWT>&& nodes,
                                Partitioner partitioner,
                                size_t num_nodes_parallel,
                                size_t num_threads) {
    num_threads = std::max(num_nodes_parallel, num_threads);

    ThreadPool thread_pool(num_threads, 100'000 * num_threads);

    // initialize buffer vectors for merging columns
    std::stack<sdsl::bit_vector> buffers;
    for (size_t i = 0; i < std::max(num_nodes_parallel, size_t(1)); ++i) {
        buffers.emplace(columns.at(0)->size());
    }

    size_t level = 0;

    while (nodes.size() > 1) {
        if (utils::get_verbose())
            std::cerr << "BRWT construction: level " << ++level << std::endl;

        auto groups = partitioner(columns);

        assert(groups.size() > 0);
        assert(groups.size() < columns.size());

        std::vector<NodeBRWT> parent_nodes(groups.size());
        std::vector<std::unique_ptr<bit_vector>> parent_columns(groups.size());

        ProgressBar progress_bar(groups.size(), "Building BRWT level",
                                 std::cerr, !utils::get_verbose());

        #pragma omp parallel for num_threads(num_nodes_parallel) schedule(dynamic)
        for (size_t g = 0; g < groups.size(); ++g) {
            const auto &group = groups[g];
            assert(group.size());

            // get available buffer vector
            sdsl::bit_vector buffer;
            #pragma omp critical
            {
                assert(buffers.size());
                buffer = std::move(buffers.top());
                buffers.pop();
            }

            auto parent = merge(subset(&nodes, group),
                                subset(&columns, group),
                                &buffer,
                                thread_pool);

            // insert buffer back to stack
            #pragma omp critical
            buffers.push(std::move(buffer));

            parent_nodes[g] = std::move(parent.first);
            parent_columns[g] = std::move(parent.second);

            ++progress_bar;
        }

        nodes = std::move(parent_nodes);
        columns = std::move(parent_columns);
    }

    return initialize(std::move(nodes.at(0)),
                      std::move(*columns.at(0)));
}

BRWT BRWTBottomUpBuilder::merge(std::vector<BRWT>&& brwts,
                                Partitioner partitioner,
                                size_t num_nodes_parallel,
                                size_t num_threads) {
    if (!brwts.size())
        return BRWT();

    std::vector<NodeBRWT> nodes(brwts.size());
    VectorsPtr columns;

    // initialize nodes
    uint64_t num_cols = 0;

    for (size_t i = 0; i < brwts.size(); ++i) {
        auto &brwt = brwts[i];
        auto &node = nodes[i];

        assert(brwt.child_nodes_.size() == brwt.assignments_.num_groups());

        for (size_t j = 0; j < brwt.child_nodes_.size(); ++j) {
            node.group_sizes.push_back(brwt.child_nodes_[j]->num_columns());

            for (size_t r = 0; r < brwt.child_nodes_[j]->num_columns(); ++r) {
                node.column_arrangement.push_back(num_cols + brwt.assignments_.get(j, r));
            }
        }

        num_cols += brwt.num_columns();

        node.child_nodes = std::move(brwt.child_nodes_);

        columns.emplace_back(new bit_vector_small(std::move(brwt.nonzero_rows_)));
    }

    return merge(std::move(columns),
                 std::move(nodes),
                 partitioner,
                 num_nodes_parallel,
                 num_threads);
}

void BRWTOptimizer::relax(BRWT *brwt_matrix,
                          uint64_t max_arity,
                          size_t num_threads) {
    assert(brwt_matrix);

    std::deque<BRWT*> parents;

    // traverse the tree and collect the nodes
    // order: leaves first, root last
    brwt_matrix->BFT([&](const BRWT &node) {
        if (node.child_nodes_.size())
            parents.push_front(const_cast<BRWT*>(&node));
    });

    ProgressBar progress_bar(parents.size(), "Optimizing BRWT",
                             std::cerr, !utils::get_verbose());

    while (!parents.empty()) {
        auto &parent = *parents.front();
        parents.pop_front();

        assert(parent.child_nodes_.size());

        NodeBRWT updated_parent;

        for (size_t g = 0; g < parent.child_nodes_.size(); ++g) {
            // we leave the column arrangement unchanged
            auto num_columns = parent.child_nodes_[g]->num_columns();
            for (size_t r = 0; r < num_columns; ++r) {
                updated_parent.column_arrangement.push_back(
                    parent.assignments_.get(g, r)
                );
            }
            // column arrangement in the parrent is set
            // now we need to initialize group_sizes and child_nodes
            add_submatrix(
                std::move(parent.child_nodes_[g]),
                &updated_parent,
                max_arity - std::min(max_arity,
                        static_cast<uint64_t>(updated_parent.group_sizes.size())
                            + parent.child_nodes_.size() - g - 1),
                num_threads
            );
        }

        parent = BRWTBuilder::initialize(std::move(updated_parent),
                                         std::move(parent.nonzero_rows_));

        ++progress_bar;
    }
}

void BRWTOptimizer::add_submatrix(std::unique_ptr<BinaryMatrix>&& submatrix,
                                  NodeBRWT *parent,
                                  uint64_t max_delta_arity,
                                  size_t num_threads) {
    assert(parent);

    const auto *brwt_cptr = dynamic_cast<const BRWT*>(submatrix.get());

    // we don't relax other submatrix representations
    bool prune = brwt_cptr != NULL;

    // we don't remove leaves in BRWT
    if (prune && !brwt_cptr->child_nodes_.size())
        prune = false;

    // check if the arity will not be exceeded
    if (prune && brwt_cptr->child_nodes_.size() > max_delta_arity)
        prune = false;

    // we relax only those child nodes, that have only BRWT nodes
    // because nodes of other types can't be reassigned
    if (prune) {
        for (const auto &child_node : brwt_cptr->child_nodes_) {
            if (!dynamic_cast<const BRWT *>(child_node.get())) {
                prune = false;
                break;
            }
        }
    }

    // we don't remove nodes if it doesn't reduce the space requirements
    if (!prune || pruning_delta(*brwt_cptr) > 0) {
        // leave this submatrix unchanged
        parent->group_sizes.push_back(submatrix->num_columns());
        parent->child_nodes.push_back(std::move(submatrix));
    } else {
        // at this point we are going to remove this node and
        // reassign all its children to its parent directly
        reassign(std::unique_ptr<BRWT>(dynamic_cast<BRWT*>(submatrix.release())),
                 parent,
                 num_threads);
    }
}

void BRWTOptimizer::reassign(std::unique_ptr<BRWT>&& node,
                             NodeBRWT *parent,
                             size_t num_threads) {
    assert(node);
    assert(parent);

    // TODO: assignemnts?

    const auto node_index = node->nonzero_rows_.convert_to<bit_vector_stat>();

    auto offset = parent->group_sizes.size();

    parent->group_sizes.resize(parent->group_sizes.size() + node->child_nodes_.size());
    parent->child_nodes.resize(parent->child_nodes.size() + node->child_nodes_.size());

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t i = 0; i < node->child_nodes_.size(); ++i) {

        auto submatrix = std::move(node->child_nodes_[i]);

        assert(dynamic_cast<BRWT*>(submatrix.get()));

        std::unique_ptr<BRWT> grand_child {
            dynamic_cast<BRWT*>(submatrix.release())
        };
        sdsl::bit_vector subindex(node_index.size(), false);

        uint64_t child_i = 0;
        node_index.call_ones([&](auto i) {
            if (grand_child->nonzero_rows_[child_i])
                subindex[i] = true;

            child_i++;
        });

        grand_child->nonzero_rows_ = bit_vector_small(subindex);

        parent->group_sizes[offset + i] = grand_child->num_columns();
        parent->child_nodes[offset + i] = std::move(grand_child);
    }
}

double BRWTOptimizer::pruning_delta(const BRWT &node) {
    assert(node.child_nodes_.size());

    double delta = 0;

    for (const auto &submatrix : node.child_nodes_) {
        auto *brwt_child = dynamic_cast<const BRWT*>(submatrix.get());
        assert(brwt_child);

        assert(brwt_child->num_rows() <= node.num_rows());
        assert(brwt_child->nonzero_rows_.size() <= node.num_rows());

        // updated vector
        delta += predict_size<decltype(brwt_child->nonzero_rows_)>(
                        node.num_rows(),
                        brwt_child->nonzero_rows_.num_set_bits());

        // old index vector
        delta -= predict_size<decltype(brwt_child->nonzero_rows_)>(
                        brwt_child->nonzero_rows_.size(),
                        brwt_child->nonzero_rows_.num_set_bits());
    }

    // removed index vector
    delta -= predict_size<decltype(node.nonzero_rows_)>(
                    node.nonzero_rows_.size(),
                    node.nonzero_rows_.num_set_bits());

    return delta;
}
