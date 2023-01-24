#include "brwt_builders.hpp"

#include <omp.h>
#include <progress_bar.hpp>

#include "common/algorithms.hpp"
#include "common/logger.hpp"
#include "common/utils/file_utils.hpp"
#include "common/vectors/vector_algorithm.hpp"


namespace mtg {
namespace annot {
namespace binmat {

using mtg::common::logger;


BRWTBottomUpBuilder::Partitioner
BRWTBottomUpBuilder::get_basic_partitioner(size_t arity) {
    assert(arity > 0u);

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

BRWT BRWTBottomUpBuilder::concatenate(std::vector<BRWT>&& submatrices,
                                      sdsl::bit_vector *buffer,
                                      ThreadPool &thread_pool) {
    assert(submatrices.size());

    if (submatrices.size() == 1)
        return std::move(submatrices[0]);

    BRWT parent;

    VectorPtrs index_columns(submatrices.size());
    uint64_t sum_num_set_bits = 0;
    for (size_t i = 0; i < submatrices.size(); ++i) {
        index_columns[i] = submatrices[i].nonzero_rows_.get();
        sum_num_set_bits += submatrices[i].nonzero_rows_->num_set_bits();
    }

    // build an aggregated parent index column
    if (sum_num_set_bits * 64 * 3 < buffer->size()) {
        // work with bits
        return concatenate_sparse(std::move(submatrices), buffer, thread_pool);
    }
    // work with uncompressed bitmap stored in buffer

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

    std::vector<std::shared_future<void>> results;

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
            [&](size_t i, sdsl::bit_vector &subindex) {
                // compress the subindex vector and set it to the child node
                // all in a single thread
                submatrices[i].nonzero_rows_
                    = std::make_unique<bit_vector_smallrank>(std::move(subindex));
                parent.child_nodes_[i].reset(new BRWT(std::move(submatrices[i])));
            }, i, std::move(subindex)
        ));
    }

    std::for_each(results.begin(), results.end(), [](auto &res) { res.wait(); });

    return parent;
}

BRWT BRWTBottomUpBuilder::concatenate_sparse(std::vector<BRWT>&& submatrices,
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

    // work with bits

    // build an aggregated parent index column
    parent.nonzero_rows_ = compute_or(index_columns, buffer->data(), thread_pool);

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

    std::vector<std::shared_future<void>> results;

    for (size_t i = 0; i < submatrices.size(); ++i) {
        // generate an index column for the child node
        sdsl::bit_vector subindex
            = generate_subindex(*submatrices[i].nonzero_rows_,
                                *parent.nonzero_rows_, thread_pool);

        // the full index vector is not needed anymore, the subindex
        // will be used instead
        submatrices[i].nonzero_rows_.reset();

        // compress the subindex vector and set it to the child node
        // all in a single thread
        submatrices[i].nonzero_rows_
            = std::make_unique<bit_vector_smallrank>(std::move(subindex));
        parent.child_nodes_[i].reset(new BRWT(std::move(submatrices[i])));
    }

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

// linkage[i] is a vector of ids of clusters merged into cluster 'i'
// linkage[c] = {} for each c < num_columns
BRWT BRWTBottomUpBuilder::build(
        const std::function<void(const CallColumn &)> &get_columns,
        const std::vector<std::vector<uint64_t>> &linkage,
        const std::filesystem::path &tmp_path,
        size_t num_nodes_parallel,
        size_t num_threads) {

    if (!linkage.size()) {
        logger->warn("Passed no linkage rules. Assembling Multi-BRWT without internal nodes...");

        std::vector<std::unique_ptr<bit_vector>> columns;

        std::mutex mu;
        uint64_t num_rows = 0;
        get_columns([&](uint64_t i, std::unique_ptr<bit_vector>&& column) {
            std::unique_lock<std::mutex> lock(mu);

            uint64_t size = column->size();
            if (!num_rows)
                num_rows = size;

            if (size != num_rows) {
                logger->error("Can't merge columns of different size");
                exit(1);
            }

            while (i >= columns.size()) {
                columns.emplace_back();
            }
            assert(!columns[i]);
            columns[i] = std::move(column);
        });

        return build(std::move(columns), get_basic_partitioner(columns.size()),
                     num_nodes_parallel, num_threads);
    }

    std::function<void(BRWT&& node, uint64_t id)> dump_node;
    std::function<BRWT(uint64_t id)> get_node;

    if (!tmp_path.empty()) {
        // keep all temp nodes in |tmp_dir|
        std::filesystem::path tmp_dir = utils::create_temp_dir(tmp_path, "brwt");
        dump_node = [tmp_dir](BRWT&& node, uint64_t id) {
            std::ofstream out(tmp_dir/std::to_string(id), std::ios::binary);
            node.serialize(out);
            node = BRWT();
        };
        get_node = [tmp_dir](uint64_t id) {
            BRWT node;
            auto filename = tmp_dir/std::to_string(id);
            std::unique_ptr<std::ifstream> in = utils::open_ifstream(filename, utils::with_mmap());
            if (!node.load(*in)) {
                logger->error("Can't load temp BRWT node {}", filename);
                exit(1);
            }
            std::filesystem::remove(filename);
            return node;
        };
    } else {
        // keep all temp nodes in memory
        auto temp_nodes = std::make_shared<std::vector<BRWT>>(linkage.size());
        dump_node = [temp_nodes](BRWT&& node, uint64_t id) {
            temp_nodes->at(id) = std::move(node);
        };
        get_node = [temp_nodes](uint64_t id) {
            return std::move(temp_nodes->at(id));
        };
    }

    std::vector<bool> done(linkage.size(), false);
    std::mutex done_mu;

    ProgressBar progress_bar(linkage.size(), "Building BRWT",
                             std::cerr, !common::get_verbose());

    uint64_t num_leaves = 0;
    uint64_t num_rows = 0;

    // prepare leaves
    get_columns([&](uint64_t i, std::unique_ptr<bit_vector>&& column) {
        uint64_t size = column->size();

        BRWT node;
        node.assignments_ = RangePartition(Partition(1, { 0 }));
        node.nonzero_rows_ = std::make_unique<bit_vector_smart>(
                                column->convert_to<bit_vector_smart>());

        dump_node(std::move(node), i);

        std::unique_lock<std::mutex> lock(done_mu);
        done[i] = true;

        if (!num_rows)
            num_rows = size;

        if (size != num_rows) {
            logger->error("Can't merge columns of different size");
            exit(1);
        }

        ++progress_bar;
        num_leaves++;
    });

    for (size_t i = 0; i < num_leaves; ++i) {
        if (linkage[i].size()) {
            logger->error("Parsed more columns than expected");
            exit(1);
        }
    }

    std::condition_variable done_cond;

    num_nodes_parallel = std::min(num_nodes_parallel, done.size());
    // initialize buffers for merging columns
    // these may be huge, so we keep only a few of them
    std::vector<sdsl::bit_vector> buffers;
    for (size_t i = 0; i < std::max(num_nodes_parallel, size_t(1)); ++i) {
        buffers.emplace_back(num_rows);
    }

    num_threads = std::max(num_nodes_parallel, num_threads);

    ThreadPool thread_pool(num_threads, 100'000 * num_threads);

    std::vector<std::vector<uint64_t>> stored_columns(linkage.size());

    #pragma omp parallel for num_threads(num_nodes_parallel) schedule(dynamic)
    for (size_t i = num_leaves; i < linkage.size(); ++i) {
        if (!linkage[i].size()) {
            logger->error("Invalid linkage: no rule for node {}", i);
            exit(1);
        }

        std::vector<BRWT> children(linkage[i].size());
        for (size_t r = 0; r < children.size(); ++r) {
            size_t j = linkage[i][r];
            {
                std::unique_lock<std::mutex> lock(done_mu);
                done_cond.wait(lock, [&]() { return done[j]; });
            }
            children[r] = get_node(j);
        }

        // merge submatrices
        BRWT node = concatenate(std::move(children),
                                &buffers.at(omp_get_thread_num()),
                                thread_pool);

        // the child nodes will be serialized separately in order not
        // to load them together with the parent next time
        std::vector<std::unique_ptr<BRWT>> child_nodes;
        for (size_t r = 0; r < node.child_nodes_.size(); ++r) {
            child_nodes.push_back(std::make_unique<BRWT>());
        }
        // replace the child nodes with dummy empty matrices
        std::swap(child_nodes, node.child_nodes_);
        // serialize the node without its child nodes
        dump_node(std::move(node), i);

        // compute column assignments for the parent
        for (size_t j : linkage[i]) {
            if (stored_columns[j].empty()) {
                // the child j is a leaf
                stored_columns[i].push_back(j);
            } else {
                // the child j is an internal node
                for (size_t c : stored_columns[j]) {
                    stored_columns[i].push_back(c);
                }
            }
        }

        std::unique_lock<std::mutex> lock(done_mu);
        if (done[i]) {
            logger->error("Invalid linkage: multiple rules for node {}", i);
            exit(1);
        }

        done[i] = true;
        done_cond.notify_all();

        for (size_t r = 0; r < children.size(); ++r) {
            dump_node(std::move(*child_nodes[r]), linkage[i][r]);
        }

        ++progress_bar;
    }

    buffers.clear();

    // All submatrices must be merged into one
    if (stored_columns.back().size() != num_leaves) {
        logger->error("Invalid linkage: all must merge into a single root");
        exit(1);
    }

    logger->trace("All {} index bitmaps have been constructed", linkage.size());
    logger->trace("Assembling Multi-BRWT...");

    std::vector<std::unique_ptr<BRWT>> nodes(linkage.size());

    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (size_t i = 0; i < linkage.size(); ++i) {
        nodes[i] = std::make_unique<BRWT>(get_node(i));
    }

    // make the root node in Multi-BRWT
    BRWT &root = *nodes.back();

    // compress the index vector
    root.nonzero_rows_ = std::make_unique<bit_vector_smallrank>(
        root.nonzero_rows_->convert_to<bit_vector_smallrank>()
    );
    // update the column arrangement to be consistent with the initial
    // order 1,2,...,m
    const auto &column_arrangement = stored_columns.back();
    std::vector<size_t> submatrix_sizes;
    for (size_t g = 0; g < root.assignments_.num_groups(); ++g) {
        submatrix_sizes.push_back(root.assignments_.group_size(g));
    }
    root.assignments_ = RangePartition(column_arrangement, submatrix_sizes);

    // put all child nodes in place
    for (size_t i = num_leaves; i < linkage.size(); ++i) {
        nodes[i]->child_nodes_.clear();
        for (size_t j : linkage[i]) {
            nodes[i]->child_nodes_.push_back(std::move(nodes[j]));
        }
    }

    return std::move(root);
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
        logger->trace("BRWT construction: level {}", level);

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
                                 std::cerr, !common::get_verbose());

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
    root.nonzero_rows_ = std::make_unique<bit_vector_smallrank>(
        root.nonzero_rows_->convert_to<bit_vector_smallrank>()
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
                             std::cerr, !common::get_verbose());

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
            = std::make_unique<bit_vector_smallrank>(std::move(subindex));

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
        delta += bit_vector_smallrank::predict_size(
                        node.num_rows(),
                        brwt_child->nonzero_rows_->num_set_bits());

        // old index vector
        delta -= bit_vector_smallrank::predict_size(
                        brwt_child->nonzero_rows_->size(),
                        brwt_child->nonzero_rows_->num_set_bits());
    }

    // removed index vector
    delta -= bit_vector_smallrank::predict_size(
                    node.nonzero_rows_->size(),
                    node.nonzero_rows_->num_set_bits());

    return delta;
}

} // namespace binmat
} // namespace annot
} // namespace mtg
