#include "BRWT_builders.hpp"

#include <omp.h>
#include <progress_bar.hpp>

#include "common/algorithms.hpp"
#include "common/vectors/int_vector_algorithm.hpp"

const uint64_t kBlockSize = 1'000'000 / 64 * 64;

// Each block is a multiple of 64 bits for thread safety
static_assert((kBlockSize % 64) == 0);


BRWT
BRWTBuilder::initialize(NodeBRWT&& node, bit_vector&& nonzero_rows) {
    BRWT brwt;

    brwt.assignments_ = RangePartition(node.column_arrangement,
                                       node.group_sizes);
    brwt.nonzero_rows_ = std::make_unique<bit_vector_small>(
        nonzero_rows.convert_to<bit_vector_small>()
    );

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

        columns.push_back(std::move(brwt.nonzero_rows_));
    }

    return merge(std::move(columns),
                 std::move(nodes),
                 partitioner,
                 num_nodes_parallel,
                 num_threads);
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

            if (node && better_split(*node)
                     && parent->child_nodes_.size() - 1
                            + node->child_nodes_.size() <= max_arity) {
                // remove this node and reassign all its children directly to its parent
                reassign(g, parent, num_threads);
            }
        }

        ++progress_bar;
    }
}

bool BRWTOptimizer::better_split(const BRWT &node) {
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
        call_ones(index_column, [&](auto i) {
            // TODO: use get_int() for faster access
            if ((*grand_child->nonzero_rows_)[child_i])
                subindex[i] = true;

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
