#ifndef __BRWT_DISK_HPP__
#define __BRWT_DISK_HPP__

#include <vector>
#include <unordered_map>
#include <memory>

#include "common/logger.hpp"

#include "common/vectors/bit_vector_adaptive.hpp"
#include "common/range_partition.hpp"
#include "common/ifstream_with_name_and_offset.hpp"

#include "annotation/binary_matrix/base/binary_matrix.hpp"

using mtg::common::logger;

namespace mtg {
namespace annot {
namespace binmat {

using NodeDepth = size_t;

class NonZeroRows {    
    // mkokot altertavielly unique_ptr could be skipped and instead somehow cleaning of underlying data could be performed
    // additional field to inform if loaded may be needed in this class
    // without unique_ptr there would be one indirection less which may or may not matter for performance
    std::unique_ptr<bit_vector> data; //maybe loaded or not
    std::streampos file_pos;
    size_t in_file_size;
    NodeDepth depth; // 0 - means root
public:
    NonZeroRows(std::streampos file_pos, size_t in_file_size, NodeDepth depth):
        file_pos(file_pos), in_file_size(in_file_size), depth(depth) {        
    }

    bit_vector* get_bit_vector() { return data.get(); }

    NodeDepth get_depth() const { return depth; }

    size_t required_size() const { return in_file_size; }

    bool is_in_memory() const { return data.get(); }

    void load_from_disk(const std::string& fname) {
        if (data) 
            return; 
        
        std::ifstream in(fname);
        if(!in)
        {
            std::cerr << "Error: cannot open file " << fname << "\n"; // mkokot TODO: move to logger
            exit(1); // mkokot TODO: return false? throw exception?
        }
        in.seekg(file_pos);

        data = std::make_unique<bit_vector_smallrank>();
        data->load(in);

        read_bytes += in_file_size;

        return;
    }

    void release_memory() {
        data.reset();        
    }

    // mkokot TODO: debug members, consider removing them    
    uint64_t disk_load_was_needed{};
    uint64_t disk_load_was_not_needed{};        
    double load_disk_time{};
    uint64_t read_bytes{};
    uint64_t slice_rows_cals{};
};

// mkokot
// the idea is to release memory of nodes close to the bottom of a three if currently used space is too much
// it may be better to extend it to give some heuristic value for a node, e.g. basing on how often it is accessed
class BRWT_DiskManager {
    size_t allowed_memory; //in bytes
    size_t used_memory = 0; //in_bytes
    std::string fname;

    // nodes that are currenty allocated, greater to first access nodes closer to the bottom of a tree
    
    // mkokot, just checking if it is better to remove oldest added nodes
    std::map<NodeDepth, std::vector<NonZeroRows*>, std::greater<>> nodes_in_memory;
    //std::map<NodeDepth, std::queue<NonZeroRows*>, std::greater<>> nodes_in_memory;
    
    std::map<NodeDepth, std::vector<NonZeroRows*>> tree_structure; // mkokot TODO: consider removing, this is only for debug purposes
    size_t tot_calls_get_rows{}, tot_calls_get_row{}, tot_calls_slice_rows{}, tot_calls_get{}; // mkokot TODO: consider removing, this is only for debug purposes
    //

    bool allocate_node_if_enough_space(NonZeroRows* node) {
        auto memory_left = allowed_memory - used_memory;
        if (node->required_size() <= memory_left) {        
            auto start = std::chrono::high_resolution_clock::now(); // mkokot dbg

            node->load_from_disk(fname);
            
            auto end =  std::chrono::high_resolution_clock::now(); // mkokot dbg
            node->load_disk_time += std::chrono::duration<double>(end - start).count(); // mkokot dbg
            node->disk_load_was_needed++; // mkokot dbg
            
            used_memory += node->required_size();
            nodes_in_memory[node->get_depth()].push_back(node);
            //nodes_in_memory[node->get_depth()].push(node);
            return true;
        }
        return false;
    }

public:
    BRWT_DiskManager(size_t allowed_memory, const std::string& fname):
        allowed_memory(allowed_memory),
        fname(fname) {
    }

    ~BRWT_DiskManager() {
        //printStats(); // mkokot, I guess calling this in destructor may cause the issues with correctness,
        //maybe it should be called at the first BRWT_Disk dtor call, like some flag put here (stats_was_printed)
    }
    void printStats() {
        static bool was_already_printed = false;
        if (was_already_printed)
            return;
        was_already_printed = true;

        // mkokot dbg, only for debuging purposes
        std::cerr << "BRWT_DiskManager stats\n"; 

        double tot_time{};
        uint64_t tot_load_was_needed{}, tot_load_was_not_needed{};
        size_t tot_req_size{};
        size_t tot_slice_rows_cals{};
        uint64_t tot_read_bytes{};
        uint64_t tot_nodes{}; 

        for (auto& [depth, nodes] : tree_structure)
        {
            std::cerr << "depth: " << depth << "\n";
            double sum_time{};
            uint64_t sum_load_was_needed{}, sum_load_was_not_needed{};
            size_t sum_req_size{};
            size_t sum_slice_rows_cals{};
            uint64_t sum_read_bytes{};
            for (auto& node : nodes)
            {
                // mkokot, this is very detailed (for each node), so probably not very informative
                std::cerr 
                    << node->load_disk_time << "\t" 
                    << node->disk_load_was_needed 
                    << "\t" << node->disk_load_was_not_needed 
                    << "\t" << node->required_size() 
                    << "\t" << node->read_bytes << "\n";

                sum_time += node->load_disk_time;
                sum_load_was_needed += node->disk_load_was_needed;
                sum_load_was_not_needed += node->disk_load_was_not_needed;
                sum_req_size += node->required_size();
                sum_slice_rows_cals += node->slice_rows_cals;
                sum_read_bytes += node->read_bytes;
            }

            std::cerr << "tot nodes: " << nodes.size() << "\n";
            std::cerr << "sum time: " << sum_time << "\n";
            std::cerr << "sum_load_was_needed: " << sum_load_was_needed << "\n";
            std::cerr << "sum_load_was_not_needed: " << sum_load_was_not_needed << "\n";
            std::cerr << "sum_req_size: " << sum_req_size << "\n";
            std::cerr << "sum_slice_rows_cals: " << sum_slice_rows_cals << "\n";
            std::cerr << "sum_read_bytes: " << sum_read_bytes << "\n";

            tot_nodes += nodes.size();
            tot_time += sum_time;
            tot_load_was_needed += sum_load_was_needed;
            tot_load_was_not_needed += sum_load_was_not_needed;
            tot_req_size += sum_req_size;
            tot_slice_rows_cals += sum_slice_rows_cals;
            tot_read_bytes += sum_read_bytes;
        }

        std::cerr << "Summary:\n";
        std::cerr << "tot nodes: " << tot_nodes << "\n";
        std::cerr << "tot time: " << tot_time << "\n";
        std::cerr << "tot_load_was_needed: " << tot_load_was_needed << "\n";
        std::cerr << "tot_load_was_not_needed: " << tot_load_was_not_needed << "\n";
        std::cerr << "tot_req_size: " << tot_req_size << "\n";
        std::cerr << "tot_slice_rows_cals: " << tot_slice_rows_cals << "\n";
        std::cerr << "tot_read_bytes: " << tot_read_bytes << "\n";
        std::cerr << "tot_calls_get_rows: " << tot_calls_get_rows << "\n";
        std::cerr << "tot_calls_get_row: " << tot_calls_get_row << "\n";
        std::cerr << "tot_calls_get: " << tot_calls_get << "\n";        
    }

    // mkokot TODO: for debuging purposes consider removing
    void register_node(uint32_t depth, NonZeroRows* node) {
        tree_structure[depth].push_back(node);
    }

    void assure_node_loaded(NonZeroRows* node, bool count_as_new_access) {
        if (node->is_in_memory()) {
            if (count_as_new_access)
                node->disk_load_was_not_needed++;
            return;
        }

        assert(used_memory <= allowed_memory);
        
        while (!allocate_node_if_enough_space(node)) {
            //we need to free some memory
            if (nodes_in_memory.empty()) { // but there it is not possible, lets force extend allowed memory
                assert(used_memory == 0);
                allowed_memory = node->required_size();
                if(!allocate_node_if_enough_space(node)) // must succeed
                {
                    std::cerr << "Error: this should not happen\n"; // mkokot TODO: better message, maybe throw exception and/or use logger
                    exit(1);
                }
                return;                
            }

            //get one of the allocated nodes close to the bottom
            auto it_to_node_to_free = nodes_in_memory.begin();
            
            auto node_to_free = it_to_node_to_free->second.back();            
            node_to_free->release_memory();
            used_memory -= node_to_free->required_size();

            //remove from nodes_in_memory
            it_to_node_to_free->second.pop_back();            
            if (it_to_node_to_free->second.empty()) 
                nodes_in_memory.erase(it_to_node_to_free);

        }
    }

    // mkokot TODO: consider removing, this is only for debug purposes
    void notify_get_rows_called() { tot_calls_get_rows++; }

    // mkokot TODO: consider removing, this is only for debug purposes
    void notify_get_row_called() { tot_calls_get_row++; }

    // mkokot TODO: consider removing, this is only for debug purposes
    void notify_get_slice_rows_called() { tot_calls_slice_rows++; }

    // mkokot TODO: consider removing, this is only for debug purposes
    void notify_get_called() { tot_calls_get++; }
};


// The Multi-BRWT compressed binary matrix representation
class BRWT_Disk : public BinaryMatrix {
    // friend class BRWTBuilder; // mkokot
    // friend class BRWTBottomUpBuilder;
    // friend class BRWTOptimizer;

    typedef uint32_t Child;

    std::shared_ptr<BRWT_DiskManager> disk_manager; // only root should allocate, using shared_ptr here is not elegant, because there is only one owner -> the root
  public:
    BRWT_Disk() = default;

    ~BRWT_Disk()
    {
        disk_manager->printStats();
    }
    
    uint64_t num_columns() const override { return assignments_.size(); }
    //uint64_t num_rows() const override { return nonzero_rows_->size(); }

    uint64_t num_rows() const override { return get_bit_vector()->size(); }

    bool get(Row row, Column column) const override;
    SetBitPositions get_row(Row row) const override;
    std::vector<SetBitPositions> get_rows(const std::vector<Row> &rows) const override;
    std::vector<Row> get_column(Column column) const override;
    // get all selected rows appended with -1 and concatenated
    std::vector<Column> slice_rows(const std::vector<Row> &rows) const override;
    // query row and get ranks of each set bit in its column
    Vector<std::pair<Column, uint64_t>> get_column_ranks(Row row) const;
    std::vector<Vector<std::pair<Column, uint64_t>>>
    get_column_ranks(const std::vector<Row> &rows) const;

    bool load(std::istream &in) override;
    void serialize(std::ostream &out) const override;

    // number of ones in the matrix
    uint64_t num_relations() const override;

    // internal stats functions
    double avg_arity() const;
    uint64_t num_nodes() const;
    double shrinking_rate() const;
    uint64_t total_column_size() const;
    uint64_t total_num_set_bits() const;

    void print_tree_structure(std::ostream &os) const;

    void set_brwt_max_anno_mem(size_t _brwt_max_anno_mem) { brwt_max_anno_mem = _brwt_max_anno_mem; }
  private:

    BRWT_Disk(std::shared_ptr<BRWT_DiskManager>& disk_manager):
        disk_manager(disk_manager) 
    {

    }

    size_t brwt_max_anno_mem = 0;

    mutable bool count_as_new_access = true; // mkokot, TODO: only for collecting some stats, may be removed

    bit_vector* get_bit_vector()
    {
        auto node = nonzero_rows_.get();
        disk_manager->assure_node_loaded(node, count_as_new_access);
        count_as_new_access = false;
        return node->get_bit_vector();
    }
    
    
    const bit_vector* get_bit_vector() const
    {
        const auto node = nonzero_rows_.get();
        disk_manager->assure_node_loaded(node, count_as_new_access);
        count_as_new_access = false;
        return node->get_bit_vector();
    }

    bool load_impl(std::istream &in, NodeDepth depth);

    // breadth-first traversal
    void BFT(std::function<void(const BRWT_Disk &node)> callback) const;
    // helper function for querying rows in batches
    template <typename T>
    std::vector<T> slice_rows(const std::vector<Row> &rows) const;

    // assigns columns to the child nodes
    RangePartition assignments_;
    //std::unique_ptr<bit_vector> nonzero_rows_;
    std::unique_ptr<NonZeroRows> nonzero_rows_;
    // generally, these child matrices can be abstract BinaryMatrix instances
    std::vector<std::unique_ptr<BRWT_Disk>> child_nodes_;
};

} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __BRWT_DISK_HPP__
