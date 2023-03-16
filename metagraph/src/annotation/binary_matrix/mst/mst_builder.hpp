#ifndef __MST_BUILDER__
#define __MST_BUILDER__

#include <vector>
#include <fstream>
#include <filesystem>

#include "annotation/binary_matrix/base/binary_matrix.hpp"
#include "common/logger.hpp"
#include "common/unix_tools.hpp"
#include "common/serialization.hpp"

// mkokot_TODO: refactor code in this file and move implementation to cpp
// also make it more general: i.e. not assuming that input annotation is row_disk

namespace mtg {
namespace annot {
namespace binmat {
namespace mst {

using mtg::common::logger;

enum class MSTStage {
    BUILD_SIGNATURES,
    COMPUTE_BINS,
    BUILD_KRUSKAL_INPUT,
    RUN_KRUSKAL,
    COMPUTE_SPARSIFIED_ROWS
};

using hamming_dist_t = uint64_t;

inline hamming_dist_t hamming_dist(const BinaryMatrix::SetBitPositions &v1,
                            const BinaryMatrix::SetBitPositions &v2) {
    hamming_dist_t res {};
    size_t i1 = 0;
    size_t i2 = 0;
    while (i1 < v1.size() && i2 < v2.size()) {
        if (v1[i1] == v2[i2]) {
            ++i1, ++i2;
        } else {
            ++res;
            if (v1[i1] < v2[i2])
                ++i1;
            else
                ++i2;
        }
    }

    res += v1.size() - i1;
    res += v2.size() - i2;
    return res;
}

//mkokot_TODO: remove this
inline std::string fix_fname(const std::string &fname) {
    if (std::filesystem::exists(fname)) {
        size_t name_iter = 0;
        while (true) {
            std::string new_name = fname + std::to_string(name_iter++);
            if (!std::filesystem::exists(new_name)) {
                logger->warn("File {} already exists, new name will be used: {}", fname,
                             new_name);
                return new_name;
            }
        }
    }
    return fname;
}

class SetBitsPerRow {
public:
    SetBitsPerRow(size_t size) : set_bits_per_row_(size) {}
    SetBitsPerRow(const std::string &fname) { Load(fname); }

    void Load(const std::string &fname) {
        std::ifstream in(fname, std::ios::binary);
        if (!in) {
            logger->error("Cannot open file {}", fname);
            throw std::ostream::failure("Bad stream");
        }

        // serialize set_bits_per_row_
        size_t num_rows = load_number(in);
        logger->info("num_rows: {}", num_rows);

        set_bits_per_row_.resize(num_rows);
        for (size_t row_id = 0; row_id < num_rows; ++row_id) {
            set_bits_per_row_[row_id] = load_number(in);
        }
        logger->info("set bits per row loaded");
    }

    void Serialize(std::string fname) {
        fname = fix_fname(fname);

        std::ofstream out(fname, std::ios::binary);

        serialize_number(out, set_bits_per_row_.size());
        for (auto set_bits : set_bits_per_row_)
            serialize_number(out, set_bits);
    }

    size_t &operator[](int64_t i) { return set_bits_per_row_[i]; }

    const std::vector<size_t> &Get() const { return set_bits_per_row_; }
private:
    std::vector<size_t> set_bits_per_row_;
};

class MinHash {
public:
    MinHash(size_t num_cols, size_t num_rows, size_t num_min_hashes = 200)
        : num_min_hashes_(num_min_hashes),
          orderings_(num_min_hashes, std::vector<uint32_t>(num_cols)),
          row_signatures_(num_rows,
                         std::vector<uint32_t>(num_min_hashes_,
                                               std::numeric_limits<uint32_t>::max())),
          set_bits_per_row_(num_rows) {
        std::default_random_engine eng;

        auto progress_bar
                = std::make_unique<ProgressBar>(num_min_hashes_, "Create orderings",
                                                std::cerr, !common::get_verbose());
        for (auto &x : orderings_) {
            std::iota(x.begin(), x.end(), 0ul);
            std::shuffle(x.begin(), x.end(), eng);
            ++(*progress_bar);
        }
    }

    std::vector<uint32_t> GetSignature(const BinaryMatrix::SetBitPositions &row) const {
        std::vector<uint32_t> signature(num_min_hashes_, std::numeric_limits<uint32_t>::max());
        for (auto set_bit_pos : row) {
            for (size_t k = 0; k < num_min_hashes_; ++k) {
                if (orderings_[k][set_bit_pos] < signature[k])
                    signature[k] = orderings_[k][set_bit_pos];
            }
        }
        return signature;
    }

    MinHash(const std::string &fname_sginatures, const std::string &fname_set_bits_per_row)
        : set_bits_per_row_(fname_set_bits_per_row) {
        loadSignatures(fname_sginatures);
    }
    SetBitsPerRow &GetSetBitsPerRow() { return set_bits_per_row_; };

    const std::vector<std::vector<uint32_t>> &GetRowSignatures() const {
        return row_signatures_;
    }
    void AddRows(const std::vector<BinaryMatrix::SetBitPositions> &rows) {
        size_t set_bits {};
        #pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic) reduction(+:set_bits)
        for (size_t row_id = 0; row_id < rows.size(); ++row_id) {
            set_bits_per_row_[row_id + offset_] = rows[row_id].size();
            set_bits += rows[row_id].size();
            for (auto set_bit_pos : rows[row_id]) {
                for (size_t k = 0; k < num_min_hashes_; ++k) {
                    if (orderings_[k][set_bit_pos] < row_signatures_[offset_ + row_id][k])
                        row_signatures_[offset_ + row_id][k] = orderings_[k][set_bit_pos];
                }
            }
            //++(*progress_bar);
        }
        tot_set_bits_ += set_bits;
        offset_ += rows.size();
    }

    size_t GetTotSetBits() const { return tot_set_bits_; }

    void SerializeSetBitsPerRow(std::string fname) { set_bits_per_row_.Serialize(fname); }
    void SerializeSignatures(std::string fname) {
        fname = fix_fname(fname);

        std::ofstream out(fname, std::ios::binary);

        // num_rows
        serialize_number(out, row_signatures_.size());
        // num_cols
        serialize_number(out, orderings_[0].size());
        // num_min_hashes_
        serialize_number(out, orderings_.size());
        // serialize signatures
        for (auto &signature : row_signatures_) {
            for (auto &val : signature)
                serialize_number32(out, val);
        }
    }
private:
    void loadSignatures(const std::string &fname) {
        std::ifstream in(fname, std::ios::binary);
        if (!in) {
            throw std::ostream::failure("Bad stream");
        }

        size_t num_rows = load_number(in);
        /*size_t num_cols =*/load_number(in);
        num_min_hashes_ = load_number(in);
        row_signatures_.resize(num_rows);
        for (size_t row_id = 0; row_id < num_rows; ++row_id) {
            row_signatures_[row_id].resize(num_min_hashes_);
            for (size_t sig_id = 0; sig_id < num_min_hashes_; ++sig_id)
                row_signatures_[row_id][sig_id] = load_number32(in);
        }
    }

    size_t num_min_hashes_;

    //mkokot_TODO: I will just generate random permutations, but it may be done better with regular hash functions
    std::vector<std::vector<uint32_t>> orderings_;
    std::vector<std::vector<uint32_t>> row_signatures_;
    SetBitsPerRow set_bits_per_row_;
    size_t tot_set_bits_ {};
    size_t offset_ = 0;


};

class LSH {
public:
    LSH(const std::vector<std::vector<uint32_t>> &rows_signatures,
        size_t num_bands = 40,
        size_t single_hash_table_size = 1ull << 26)
        : rows_signatures_(rows_signatures),
          num_bands_(num_bands),
          single_hash_table_size_(single_hash_table_size) {
        assert(rows_signatures[0].size() % num_bands == 0);
    }

    std::vector<std::vector<size_t>> AssignToBins() {
        logger->info("num_bands: {}", num_bands_);
        size_t num_hashes_per_band = rows_signatures_[0].size() / num_bands_;

        auto num_rows = rows_signatures_.size();


        std::vector<std::vector<std::list<size_t>>> hash_tables(num_bands_);
        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t band_id = 0; band_id < num_bands_; ++band_id)
            hash_tables[band_id] = std::vector<std::list<size_t>>(single_hash_table_size_);

        auto progress_bar = std::make_unique<ProgressBar>(num_rows, "Compute LSH", std::cerr,
                                                          !common::get_verbose());

        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t row_id = 0; row_id < num_rows; ++row_id) {
            for (size_t band_id = 0; band_id < num_bands_; ++band_id) {
                auto band_start = band_id * num_hashes_per_band;
                auto band_end = band_start + num_hashes_per_band;
                std::vector<uint32_t> band(rows_signatures_[row_id].begin() + band_start,
                                           rows_signatures_[row_id].begin() + band_end);

                uint64_t hash_val = hash_vec(band, band_id) % single_hash_table_size_;
                #pragma omp critical
                hash_tables[band_id][hash_val].push_back(row_id);
            }
            ++(*progress_bar);
        }


        std::vector<std::vector<size_t>> res;

        size_t max_of_all_bands_bucket_size {};
        for (size_t band_id = 0; band_id < num_bands_; ++band_id) {
            logger->info("band {}", band_id);

            size_t tot_size {};
            size_t max_bucket_size {};
            size_t num_non_empty_buckets {};


            for (size_t i = 0; i < single_hash_table_size_; ++i) {
                auto &bucket = hash_tables[band_id][i];

                if (!bucket.empty())
                    res.emplace_back(bucket.begin(), bucket.end());

                tot_size += bucket.size();
                if (!bucket.empty())
                    num_non_empty_buckets++;
                if (bucket.size() > max_bucket_size)
                    max_bucket_size = bucket.size();
            }
            logger->info("tot size: {}", tot_size);
            logger->info("num_non_empty_buckets: {}", num_non_empty_buckets);
            logger->info("max_bucket_size: {}", max_bucket_size);

            if (max_bucket_size > max_of_all_bands_bucket_size)
                max_of_all_bands_bucket_size = max_bucket_size;
        }

        size_t max_res_bucket {};
        for (auto &x : res) {
            if (x.size() > max_res_bucket)
                max_res_bucket = x.size();
        }

        logger->info("max_of_all_bands_bucket_size: {}", max_of_all_bands_bucket_size);
        logger->info("max_res_bucket: {}", max_res_bucket);
        return res;
    }
private:
    std::size_t hash_vec(const std::vector<uint32_t> &vec, size_t modifier) {
        std::size_t seed = vec.size() + modifier;
        for (auto x : vec) {
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    };

    const std::vector<std::vector<uint32_t>> &rows_signatures_;
    size_t num_bands_;
    size_t single_hash_table_size_;
};

struct edge {
    uint64_t idx1, idx2;
    hamming_dist_t dist;

    edge(uint64_t idx1,
         uint64_t idx2,
         hamming_dist_t dist = std::numeric_limits<hamming_dist_t>::max())
        : idx1(idx1), idx2(idx2), dist(dist) {}
};

class CandidateCollector {
public:
    CandidateCollector(const std::vector<size_t> &num_set_bits_per_row, size_t max_neighbours)
        : num_set_bits_per_row(num_set_bits_per_row),
          candidates(this->num_set_bits_per_row.size()),
          max_neighbours(max_neighbours)
    {}

    void AddCandidate(size_t u,
                      size_t v,
                      const BinaryMatrix::SetBitPositions &full_u,
                      const BinaryMatrix::SetBitPositions &full_v) {
        if (v < u)
            std::swap(u, v);

        if (alreadyContained(u, v))
            return;

        auto new_dist = hamming_dist(full_u, full_v);
        #pragma omp critical
        if (haveSpace(u)) {
            std::lock_guard<std::mutex> lck(candidates_mtx);
            candidates[u].emplace_back(v, new_dist);
        } else {
            auto worst_id = getWorst(u);

            if (candidates[u][worst_id].dist > get_best_possible_hamming(u, v)) {
                if (new_dist < candidates[u][worst_id].dist) {
                    candidates[u][worst_id].dist = new_dist;
                    candidates[u][worst_id].idx = v;
                }
            }
        }
    }

    std::vector<edge> GetEdges() {
        std::vector<edge> res;
        for (size_t u = 0; u < candidates.size(); ++u)
            for (auto &[v, dist] : candidates[u])
                res.emplace_back(u, v, dist);
        return res;
    }

private:
    const std::vector<size_t> &num_set_bits_per_row;
    struct elem_t {
        size_t idx {};
        hamming_dist_t dist {};
        elem_t(size_t idx, hamming_dist_t dist) : idx(idx), dist(dist) {}
    };
    std::vector<std::vector<elem_t>> candidates;
    std::mutex candidates_mtx;
    size_t max_neighbours;

    bool alreadyContained(size_t u, size_t v) {
        assert(u < v);
        assert(u < candidates.size());
        std::lock_guard<std::mutex> lck(candidates_mtx);
        auto &tmp = candidates[u];
        return std::find_if(tmp.begin(), tmp.end(),
                            [v](const auto &x) { return x.idx == v; })
                != tmp.end();
    }

    size_t getWorst(size_t u) {
        auto &tmp = candidates[u];
        assert(tmp.size() == max_neighbours);
        hamming_dist_t worst_hamming = tmp[0].dist;
        size_t worst_idx = 0;
        for (size_t i = 1; i < tmp.size(); ++i)
            if (tmp[i].dist > worst_hamming) {
                worst_hamming = tmp[i].dist;
                worst_idx = i;
            }
        return worst_idx;
    }
    hamming_dist_t get_best_possible_hamming(size_t u, size_t v) {
        auto sbu = num_set_bits_per_row[u];
        auto sbv = num_set_bits_per_row[v];
        return sbu < sbv ? sbv - sbu : sbu - sbv;
    }

    bool haveSpace(size_t u) { return candidates[u].size() < max_neighbours; }
};

class BatchedKruskal {
    class DisjointSet {
        struct elem {
            uint64_t parent;
        };

        std::vector<elem> s;

      public:
        DisjointSet(size_t size) : s(size) {
            for (size_t i = 0; i < s.size(); ++i)
                s[i].parent = i;
        }
        uint64_t Find(uint64_t x) {
            if (x != s[x].parent) {
                s[x].parent = Find(s[x].parent);
                return s[x].parent;
            }
            return x;
        }

        void Union(uint64_t x, uint64_t y) {
            x = Find(x);
            y = Find(y);
            if (x == y)
                return;
            s[y].parent = x;
        }

        std::vector<size_t> GetConnectedComponents() const {
            std::vector<size_t> res;
            for (size_t i = 0; i < s.size(); ++i)
                if (s[i].parent == i)
                    res.push_back(i);
            return res;
        }
    };

    std::vector<size_t> num_set_bits_per_row;
    std::vector<edge> selected_edges;
    DisjointSet disjoint_set;

    // mkokot_TODO: remove times?
    double prepare_edges_time = {};
    double compute_costs_time = {};
    double sort_edges_time = {};
    double selecting_edges_time = {};

  public:
    BatchedKruskal(std::vector<size_t> &&num_set_bits_per_row)
        : num_set_bits_per_row(std::move(num_set_bits_per_row)),
          disjoint_set(this->num_set_bits_per_row.size()) {}
    void AddEdges(std::vector<edge> &edges) {
        ips4o::sort(edges.begin(), edges.end(),
                    [](const auto &e1, const auto &e2) { return e1.dist < e2.dist; });
        // this is not stable and may result in having slightly different results per each run
        // I mean MST cost will be the same, but the tree roots may be different
        // std::sort(edges.begin(), edges.end(), [](const auto &e1, const auto &e2) { return e1.dist < e2.dist; });

        uint64_t edges_pos = 0;

        while (edges_pos < edges.size()) {
            auto edge = edges[edges_pos++];

            // this edge is not helpful
            if (edge.dist > num_set_bits_per_row[edge.idx1]
                && edge.dist > num_set_bits_per_row[edge.idx2])
                continue;

            if (disjoint_set.Find(edge.idx1) != disjoint_set.Find(edge.idx2)) {
                disjoint_set.Union(edge.idx1, edge.idx2);
                selected_edges.push_back(edge);
            }
        }
    }
    using LoadRow = std::function<BinaryMatrix::SetBitPositions(size_t /*row_id*/)>;
    void AddRows(std::vector<size_t> &batch, LoadRow load_row) {
        size_t num_rows = num_set_bits_per_row.size();
        size_t max_edges_in_MST = num_rows - 1;

        // we have full mst
        if (selected_edges.size() >= max_edges_in_MST) {
            logger->warn("Full MST is already created!");
            return;
        }
        // filter out ids that are already connected
        std::vector<edge> kruskal_input;

        Timer timer;

        // mkokot_TODO: move to param?
        size_t max_considered = 10;
        std::sort(batch.begin(), batch.end(), [this](size_t u, size_t v) {
            return num_set_bits_per_row[u] < num_set_bits_per_row[v];
        });
        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 0; i < batch.size(); ++i) {
            size_t added = 0;
            for (size_t j = i + 1; j < batch.size() && added < max_considered; ++j) {
                auto u = batch[i];
                auto v = batch[j];
                assert(u < num_rows);
                assert(v < num_rows);
                ++added;
                if (disjoint_set.Find(u) != disjoint_set.Find(v)) {
                    #pragma omp critical
                    kruskal_input.emplace_back(u, v);
                }
            }
        }

        prepare_edges_time += timer.elapsed();

        if (kruskal_input.empty())
            return;

        // compute costs
        timer.reset();
        #pragma omp parallel for num_threads(get_num_threads())
        for (size_t i = 0; i < kruskal_input.size(); ++i) {
            kruskal_input[i].dist = hamming_dist(load_row(kruskal_input[i].idx1),
                                                 load_row(kruskal_input[i].idx2));
        }
        compute_costs_time += timer.elapsed();

        timer.reset();
        ips4o::sort(kruskal_input.begin(), kruskal_input.end(),
                    [](const auto &e1, const auto &e2) { return e1.dist < e2.dist; });
        sort_edges_time += timer.elapsed();
        uint64_t edges_pos = 0;

        timer.reset();
        while (selected_edges.size() < max_edges_in_MST && edges_pos < kruskal_input.size()) {
            auto edge = kruskal_input[edges_pos++];

            // this edge is not helpful
            if (edge.dist > num_set_bits_per_row[edge.idx1]
                && edge.dist > num_set_bits_per_row[edge.idx2])
                continue;

            if (disjoint_set.Find(edge.idx1) != disjoint_set.Find(edge.idx2)) {
                disjoint_set.Union(edge.idx1, edge.idx2);
                selected_edges.push_back(edge);
            }
        }
        selecting_edges_time += timer.elapsed();
    }

    void PrintTimes() {
        logger->info("prepare_edges_time: {}", prepare_edges_time);
        logger->info("compute_costs_time: {}", compute_costs_time);
        logger->info("sort_edges_time: {}", sort_edges_time);
        logger->info("selecting_edges_time: {}", selecting_edges_time);
    }
    size_t GetTotCost() const {
        size_t res {};
        for (const auto &e : selected_edges)
            res += e.dist;
        return res;
    }

    size_t GetTotSetBitsInput() const {
        size_t res {};
        for (auto x : num_set_bits_per_row)
            res += x;
        return res;
    }

    size_t GetNumEdgesInResult() const { return selected_edges.size(); }

    std::vector<size_t> GetConnectedComponents() const {
        return disjoint_set.GetConnectedComponents();
    }

    const std::vector<edge> &GetResult() const { return selected_edges; }

    const std::vector<size_t> &GetSetBitsPerRow() const { return num_set_bits_per_row; }
};

struct parent {
    uint64_t idx;
    hamming_dist_t dist;

  public:
    parent() = default;
    parent(uint64_t idx, hamming_dist_t dist) : idx(idx), dist(dist) {}
};

size_t FindNodeWithLowestVal(const std::vector<std::list<parent>> &Graph,
                             uint64_t start_node,
                             const std::vector<size_t> &node_values,
                             std::vector<uint8_t> &visited) {
    uint64_t cur_min_id = start_node;
    std::queue<size_t> Q;

    Q.push(start_node);
    // std::vector<uint8_t> visited(node_values.size());
    while (!Q.empty()) {
        auto id = Q.front();
        Q.pop();
        if (visited[id]) {
            logger->error("Double visiting node!");
        }
        visited[id] = true;
        if (node_values[id] < node_values[cur_min_id])
            cur_min_id = id;

        for (const auto &e : Graph[id]) {
            if (!visited[e.idx])
                Q.push(e.idx);
        }
    }
    return cur_min_id;
}

// mkokot_TODO: clean this up

// node values are number of set bits per appropriate row
// the idea is to split trees that if a node/row self-cost (number of ones) is better than the MST transformation use this node as a root
std::vector<parent> ComputeTrees(const std::vector<edge> &MST,
                                 const std::vector<size_t> &connected_components,
                                 const std::vector<size_t> &set_bits_per_row,
                                 size_t &num_sparsified_set_bits) {
    num_sparsified_set_bits = 0;
    std::vector<parent> parents(set_bits_per_row.size());
    if (set_bits_per_row.size() == 1) {
        assert(MST.empty());
        parents[0].idx = 0;
        parents[0].dist = set_bits_per_row[0];
        num_sparsified_set_bits = set_bits_per_row[0];
        return parents;
    }

    std::vector<std::list<parent>> G(set_bits_per_row.size());
    for (const auto &edge : MST) {
        G[edge.idx1].emplace_back(edge.idx2, edge.dist);

        G[edge.idx2].emplace_back(edge.idx1, edge.dist);
    }

    logger->info("connected_components.size(): {}", connected_components.size());

    auto progress_bar
            = std::make_unique<ProgressBar>(connected_components.size(), "Traverse trees",
                                            std::cerr, !common::get_verbose());

    std::vector<uint8_t> done(set_bits_per_row.size(), false);
    std::vector<uint8_t> done_for_find_root(set_bits_per_row.size(), false);
    std::vector<uint64_t> roots;
    double time_find_lowest {};
    double time_rest {};
    // turns out this is very fast and does not need to be parallelized, but it is to be
    // considered if uncomment this for (size_t node_id : connected_components) {
    //#pragma omp parallel for num_threads(get_num_threads()) schedule(dynamic)
    //reduction(+:time_find_lowest) reduction(+:time_rest)
    for (size_t i = 0; i < connected_components.size(); ++i) {
        size_t node_id = connected_components[i];

        auto start = std::chrono::high_resolution_clock::now();
        auto root_node
                = FindNodeWithLowestVal(G, node_id, set_bits_per_row, done_for_find_root);
        time_find_lowest += std::chrono::duration<double>(
                                    std::chrono::high_resolution_clock::now() - start)
                                    .count();

        start = std::chrono::high_resolution_clock::now();
        std::queue<uint64_t> Q;
        parents[root_node].idx = root_node;
        parents[root_node].dist
                = set_bits_per_row[root_node]; // self dist is cost of plain representation
        #pragma omp critical
        { roots.push_back(root_node); }
        Q.push(root_node);

        // mkokot_TODO: remove?
        if (done[root_node]) {
            logger->error("Something wrong, probably with connected components");
        }
        done[root_node] = true;

        while (!Q.empty()) {
            auto v = Q.front();
            Q.pop();
            for (const auto &x : G[v]) {
                auto child_id = x.idx;
                auto dist = x.dist;
                if (!done[child_id]) {
                    done[child_id] = true;
                    parents[child_id].idx = v;
                    parents[child_id].dist = dist;
                    Q.push(child_id);
                }
            }
        }
        time_rest += std::chrono::duration<double>(
                             std::chrono::high_resolution_clock::now() - start)
                             .count();
        ++(*progress_bar);
    }

    logger->info("time_find_lowest: {}", time_find_lowest);
    logger->info("time_rest: {}", time_rest);

    size_t cost_before_splitting_tree {};
    for (size_t node_id = 0; node_id < set_bits_per_row.size(); ++node_id) {
        cost_before_splitting_tree += parents[node_id].dist;
    }

    logger->info("cost_before_splitting_tree: {}", cost_before_splitting_tree);

    logger->info("roots.size() before opt: {}", roots.size());
    for (size_t node_id = 0; node_id < set_bits_per_row.size(); ++node_id) {
        if (parents[node_id].idx == node_id) // already is a root
            continue;

        if (set_bits_per_row[node_id] < parents[node_id].dist) { // convert node into new root
            parents[node_id].dist = set_bits_per_row[node_id];
            parents[node_id].idx = node_id;
            roots.push_back(node_id);
        }
    }

    logger->info("roots.size() after opt: {}", roots.size());

    // now all should be done
    if (!std::all_of(done.begin(), done.end(), [](bool v) { return v; })) {
        logger->error("Not all done, but should be");
    }

    // now all should be done
    if (!std::all_of(done_for_find_root.begin(), done_for_find_root.end(),
                     [](bool v) { return v; })) {
        logger->error("Not all done_for_find_root, but should be");
    }

    size_t cost_after_splitting_tree {};
    for (size_t node_id = 0; node_id < set_bits_per_row.size(); ++node_id) {
        cost_after_splitting_tree += parents[node_id].dist;
    }

    uint64_t tot_ones_in_rows {};
    for (const auto &node : set_bits_per_row) {
        tot_ones_in_rows += node;
    }

    logger->info("tot_ones_in_rows: {}", tot_ones_in_rows);

    logger->info("cost_after_splitting_tree: {}", cost_after_splitting_tree);

    num_sparsified_set_bits = cost_after_splitting_tree;


    // mkokot_TODO: lets try to cut these trees so that the depth <= max_depth
    {
        std::vector<std::list<size_t>> G(set_bits_per_row.size());

        // mkokot_TODO: make this a parameter
        const size_t max_depth = 100;
        logger->info("max_depth: {}", max_depth);
        std::queue<size_t> roots;
        for (size_t node_id = 0; node_id < set_bits_per_row.size(); ++node_id) {
            if (parents[node_id].idx == node_id) {
                roots.emplace(node_id);
                continue;
            }
            G[parents[node_id].idx].push_back(node_id);
        }

        size_t roots_after_max_depth = roots.size();

        while (!roots.empty()) {
            size_t root = roots.front();
            roots.pop();
            std::queue<std::pair<size_t, size_t>> Q; // node_id, dist
            Q.emplace(root, 0);
            while (!Q.empty()) {
                auto [node, dist] = Q.front();
                Q.pop();
                assert(dist <= max_depth);
                if (dist == max_depth) { // make all childs new roots
                    for (auto child : G[node]) {
                        roots.emplace(child);
                        ++roots_after_max_depth;
                        parents[child].idx = child;
                        parents[child].dist = set_bits_per_row[child];
                    }
                } else { // dist is OK, traverse childs
                    for (auto child : G[node])
                        Q.emplace(child, dist + 1);
                }
            }
        }

        size_t cost_after_max_depth {};
        for (size_t node_id = 0; node_id < set_bits_per_row.size(); ++node_id) {
            cost_after_max_depth += parents[node_id].dist;
        }

        logger->info("roots_after_max_depth: {}", roots_after_max_depth);
        logger->info("cost_after_max_depth: {}", cost_after_max_depth);

        // split trees again
        size_t roots_after_split_trees_again {};
        for (size_t node_id = 0; node_id < set_bits_per_row.size(); ++node_id) {
            if (parents[node_id].idx == node_id) // already is a root
                continue;

            if (set_bits_per_row[node_id]
                < parents[node_id].dist) { // convert node into new root
                parents[node_id].dist = set_bits_per_row[node_id];
                parents[node_id].idx = node_id;
                ++roots_after_split_trees_again;
            }
        }

        logger->info("number of new roots_after_split_trees_again: {}",
                     roots_after_split_trees_again);

        size_t cost_afters_split_trees_again {};
        for (size_t node_id = 0; node_id < set_bits_per_row.size(); ++node_id) {
            cost_afters_split_trees_again += parents[node_id].dist;
        }
        logger->info("cost_afters_split_trees_again: {}", cost_afters_split_trees_again);

        num_sparsified_set_bits = cost_afters_split_trees_again;
    }

    // mkokot_TODO: turs out it is worth to reverse again, but I don't understand why
    // so it must be investigated, maybe should be done in loop
    // the code below reverses all the remaining trees so that the minimal row becames a root    
    /*
            G.clear();
            G.resize(set_bits_per_row.size());
            for(size_t node_id = 0; node_id < parents.size() ; ++node_id) {
                auto u = node_id;
                auto v = parents[u].idx;
                if (u == v) { //skip self-conections, they indicate roots
                    continue;
                }
                auto dist = parents[u].dist;

                G[u].emplace_back(v, dist);
                G[v].emplace_back(u, dist);
            }

            std::fill(done.begin(), done.end(), false);

            size_t num_changed_roots{};
            size_t num_keeped_roots{};
            for(size_t old_root : roots) {
                auto new_root = FindNodeWithLowestVal(G, old_root, set_bits_per_row);

                if(new_root != old_root)
                    ++num_changed_roots;
                else
                    ++num_keeped_roots;

                std::queue<uint64_t> Q;
                parents[new_root].idx = new_root;
                parents[new_root].dist = set_bits_per_row[new_root]; // self dist is cost
       of plain representation roots.push_back(new_root); Q.push(new_root);

                if (done[old_root]) {
                    logger->error("Something wrong, probably with computing roots");
                }

                done[new_root] = true;

                while (!Q.empty()) {
                    auto v = Q.front();
                    Q.pop();
                    for (const auto &x : G[v]) {
                        auto child_id = x.idx;
                        auto dist = x.dist;
                        if (!done[child_id]) {
                            done[child_id] = true;
                            parents[child_id].idx = v;
                            parents[child_id].dist = dist;
                            Q.push(child_id);
                        }
                    }
                }
            }

            logger->info("num_changed_roots: {}", num_changed_roots);
            logger->info("num_keeped_roots: {}", num_keeped_roots);

            //now all should be done
            if(!std::all_of(done.begin(), done.end(), [](bool v) { return v;})) {
                logger->error("Not all done, but should be");
            }

            size_t cost_after_reverting_tree {};
            for (size_t node_id = 0; node_id < set_bits_per_row.size(); ++node_id) {
                cost_after_reverting_tree += parents[node_id].dist;
            }

            logger->info("cost_after_reverting_tree: {}", cost_after_reverting_tree);
    */

    return parents;
}

struct xored_row {
    uint64_t parent_idx;
    BinaryMatrix::SetBitPositions set_bits;
};

BinaryMatrix::SetBitPositions bit_xor(const BinaryMatrix::SetBitPositions &v1,
                                      const BinaryMatrix::SetBitPositions &v2) {
    BinaryMatrix::SetBitPositions res;

    size_t i1 = 0;
    size_t i2 = 0;
    while (i1 < v1.size() && i2 < v2.size()) {
        if (v1[i1] == v2[i2]) {
            ++i1, ++i2;
        } else {
            if (v1[i1] < v2[i2])
                res.push_back(v1[i1++]);
            else
                res.push_back(v2[i2++]);
        }
    }
    while (i1 < v1.size())
        res.push_back(v1[i1++]);
    while (i2 < v2.size())
        res.push_back(v2[i2++]);

    return res;
}

class SparsifyRows {
    uint64_t tot_set_bits {};
    sdsl::int_vector<64> parents;

  public:
    const sdsl::int_vector<64> &GetParents() const { return parents; }
    SparsifyRows(size_t num_rows) : parents(num_rows) {}

    std::vector<xored_row> GetSparsifiedRows(const std::vector<parent> &parents_batch,
                                             RowDisk &row_disk,
                                             size_t offset) {
        // as a result I should get a batch of rows, where I have new set_bits, parent_row_id

        std::vector<xored_row> results(parents_batch.size());

        for (uint32_t i = 0; i < parents_batch.size(); ++i) {
            auto parent_idx = parents_batch[i].idx;
            auto node_id = i + offset;

            #pragma omp critical
            parents[node_id] = parent_idx;

            results[i].parent_idx = parent_idx;

            if (node_id == parent_idx) // root
                results[i].set_bits = row_disk.get_row(node_id);
            else
                results[i].set_bits
                        = bit_xor(row_disk.get_row(parent_idx), row_disk.get_row(node_id));

            if (results[i].set_bits.size() != parents_batch[i].dist) {
                logger->info("Error");
               
                logger->info("node_id: {}, {}", node_id, i);
                logger->info("parents_batch[i].dist: {}", parents_batch[i].dist);
                logger->info("parents_batch[i].idx: {}", parents_batch[i].idx);
                logger->info("results[i].set_bits.size(): {}",
                             results[node_id].set_bits.size());
                //   logger->info("recomputed hamming: {}",
                //                mkokot_kmeans::hamming_dist(node, parent));
            } else {
                // logger->info("OK");
            }
            #pragma omp critical
            tot_set_bits += results[i].set_bits.size();
        }
        //            logger->info("GetSparsifiedRows, tot_set_bits: {}", tot_set_bits);
        return results;
    }

    size_t GetTotSetBits() const { return tot_set_bits; }
};

} // namespace mst
} // namespace binmat
} // namespace annot
} // namespace mtg

#endif // __MST_BUILDER__