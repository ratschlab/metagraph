#include "annotated_graph_algorithm.hpp"

#include "common/utils/string_utils.hpp"
#include "common/logger.hpp"
#include "common/vectors/vector_algorithm.hpp"
#include "common/vectors/bitmap.hpp"
#include "graph/representation/masked_graph.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "differential_tests.hpp"

namespace mtg {
namespace graph {

using mtg::common::logger;

typedef AnnotatedDBG::node_index node_index;
typedef AnnotatedDBG::Annotator Annotator;
typedef AnnotatedDBG::Annotator::Label Label;
using Column = annot::binmat::BinaryMatrix::Column;

typedef std::function<size_t()> LabelCountCallback;

constexpr std::memory_order MO_RELAXED = std::memory_order_relaxed;


/**
 * Return an int_vector<>, bit_vector of lengths anno_graph.get_graph().max_index() * 2
 * and anno_graph.get_graph().max_index(), respectively, and a bit_vector of length
 * equal to the total number of labels.
 * For an index i, the int_vector at indices 2*i and 2*i + 1 represent the
 * number of labels in labels_in and labels_out which the k-mer of index i is
 * annotated with, respectively. The width of the int_vector<> is computed to be
 * wide enough to contain counts up to num_labels.
 * The returned bit_vector is a k-mer mask indicating those k-mers annotated
 * with at least one in-label. If add_out_labels_to_mask is true, then it also indicates
 * which those k-mers with at least one out-label.
 * The second bit vector indicates which of the labels correspond to an "other"
 * label (i.e., not in- or out-).
 */


using ValueCallback = std::function<void(uint64_t /* row index */, uint64_t /* row value */)>;
using ValueGenerator = std::function<void(const ValueCallback&)>;
using ColumnCallback = std::function<void(uint64_t /* offset */, const Label& /* annotation column */, const ValueGenerator&)>;
using ColumnGenerator = std::function<void(const ColumnCallback&)>;

template <typename T>
std::string outstring(std::vector<T> out_array);
template <typename T>
std::string outstring(std::vector<T> out_array){
    std::string out = "";
    for (auto &item: out_array) { std::string x = std::to_string(item); out = out + x + ", "; }
    return out;
}

void kmer_distribution_table(const ColumnGenerator &generate_columns,
                             const tsl::hopscotch_set<Label> &labels_in,
                             const tsl::hopscotch_set<Label> &labels_out,
                             size_t max_index,
                             MaskedDeBruijnGraph &masked_graph, // std::shared_ptr<const DeBruijnGraph> graph_ptr,
                             const DifferentialAssemblyConfig &config,
                             size_t num_threads);  // For kmer distribution.

std::tuple<sdsl::int_vector_buffer<>, sdsl::bit_vector, sdsl::bit_vector, std::tuple<size_t, size_t>>
construct_diff_label_count_vector(const ColumnGenerator &generate_columns,
                                  const std::vector<std::string> &files,
                                  const tsl::hopscotch_set<Label> &labels_in,
                                  const tsl::hopscotch_set<Label> &labels_out,
                                  const DifferentialAssemblyConfig &config,
                                  size_t max_index,
                                  size_t num_labels,
                                  size_t num_threads,
                                  bool add_out_labels_to_mask);

// Regions of a graph mask which should be kept (i.e., masked in)
typedef std::vector<std::pair<size_t, size_t>> Intervals;

// Returns a vector of kept regions given a unitig and its corresponding path
typedef std::function<Intervals(const std::string&, const std::vector<node_index>&)> GetKeptIntervals;

// Assemble unitigs from the masked graph, then use get_kept_intervals to generate
// regions which should be masked in. Update the graph mask accordingly.
void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const GetKeptIntervals &get_kept_intervals,
                                   size_t num_threads);

// Given an initial mask and counts, generate a masked graph. If add_complement
// is true, then add the reverse complements of all nodes to the graph as well.
std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                          sdsl::int_vector_buffer<> &counts,
                          sdsl::bit_vector&& mask,
                          bool add_complement,
                          size_t num_threads);

// Helper function for calls in interface
std::shared_ptr<MaskedDeBruijnGraph>
mask_nodes_by_label(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                    const AnnotatedDBG *anno_graph,
                    sdsl::int_vector_buffer<>&& counts,
                    sdsl::bit_vector&& mask,
                    sdsl::bit_vector&& other_mask,
                    size_t num_labels,
                    const tsl::hopscotch_set<Label> &labels_in,
                    const tsl::hopscotch_set<Label> &labels_out,
                    const tsl::hopscotch_set<Label> &labels_in_round2,
                    const tsl::hopscotch_set<Label> &labels_out_round2,
                    std::tuple<size_t, size_t> total_kmers,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads);

std::shared_ptr<MaskedDeBruijnGraph> // This version of mask_nodes_by_label is for when the columns are already loaded.
mask_nodes_by_label(const AnnotatedDBG &anno_graph,
                    const tsl::hopscotch_set<Label> &labels_in,
                    const tsl::hopscotch_set<Label> &labels_out,
                    const tsl::hopscotch_set<Label> &labels_in_round2,
                    const tsl::hopscotch_set<Label> &labels_out_round2,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads,
                    size_t num_parallel_files) {
    num_parallel_files = std::min(num_threads, num_parallel_files);

    auto graph_ptr = std::static_pointer_cast<const DeBruijnGraph>(
        anno_graph.get_graph_ptr()
    );

    bool check_other = config.label_mask_other_unitig_fraction != 1.0;
    bool unitig_mode = check_other || labels_in_round2.size() || labels_out_round2.size()
            || config.label_mask_in_unitig_fraction != 0.0
            || config.label_mask_out_unitig_fraction != 1.0
            || config.label_mask_other_unitig_fraction != 1.0;

    bool add_complement = graph_ptr->get_mode() == DeBruijnGraph::CANONICAL
        && (config.add_complement || unitig_mode);

    logger->trace("Generating initial mask");

    auto count_vector = construct_diff_label_count_vector(
        [&](const ColumnCallback &column_callback) {
            if (config.count_kmers) {
                throw std::runtime_error("AnnotatedDBG with counts not supported yet");
            } else {
                size_t offset = 0;
                for (const auto &label : labels_in) {
                    column_callback(offset, label, [&](const ValueCallback &value_callback) {
                        anno_graph.call_annotated_nodes(label, [&](node_index i) {
                            i = AnnotatedDBG::graph_to_anno_index(i);
                            value_callback(i, 1);
                        });
                    });
                    ++offset;
                }
                for (const auto &label : labels_out) {
                    column_callback(offset, label, [&](const ValueCallback &value_callback) {
                        anno_graph.call_annotated_nodes(label, [&](node_index i) {
                            i = AnnotatedDBG::graph_to_anno_index(i);
                            value_callback(i, 1);
                        });
                    });
                    ++offset;
                }
            };
        },
        { }, labels_in, labels_out, config, graph_ptr->max_index(),
        std::max(labels_in.size(), labels_out.size()), num_parallel_files,
        add_complement
    );

    auto &[counts, init_mask, other_labels, total_kmers] = count_vector; //ERROR  only 3 names provided for structured binding, std::vector<long unsigned int, std::allocator<long unsigned int> > >’ decomposes into 4 elements
    sdsl::bit_vector other_mask(init_mask.size() * check_other, false);
    return mask_nodes_by_label(graph_ptr, &anno_graph,
                               std::move(counts), std::move(init_mask),
                               std::move(other_mask),
                               anno_graph.get_annotator().num_labels(),
                               labels_in, labels_out,
                               labels_in_round2, labels_out_round2, total_kmers,
                               config, num_threads);
}



std::shared_ptr<MaskedDeBruijnGraph>
mask_nodes_by_label(std::shared_ptr<const DeBruijnGraph> graph_ptr, // Myrthe: this version i don't need to change at first
                    const std::vector<std::string> &files,
                    const tsl::hopscotch_set<Label> &labels_in,
                    const tsl::hopscotch_set<Label> &labels_out,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads,
                    size_t num_parallel_files) {
    bool check_other = config.label_mask_other_unitig_fraction != 1.0;
    bool unitig_mode = check_other
            || config.label_mask_in_unitig_fraction != 0.0
            || config.label_mask_out_unitig_fraction != 1.0
            || config.label_mask_other_unitig_fraction != 1.0;

    bool add_complement = graph_ptr->get_mode() == DeBruijnGraph::CANONICAL
        && (config.add_complement || unitig_mode);

    logger->trace("Generating initial mask");



    auto count_vector = construct_diff_label_count_vector(
            [&](const ColumnCallback &column_callback) {
            if (config.count_kmers) {
                assert(files.size() > 0);
                annot::ColumnCompressed<>::load_columns_and_values(files,  // TODO instead of looping per column, m_buffersize. Give an extra argument to the load_cols_and_vals.
                    [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> && column, auto&& column_values) {
                        column_callback(offset, label, [&](const ValueCallback &value_callback) {
                            assert(std::accumulate(column_values.begin(), column_values.end(), 0) > 0);
                            assert(column->num_set_bits() > 0);
                            assert(column->num_set_bits() == column_values.size());
                            int normalization_constant = 1e9;  //Normalization. This constant should be the same for every column. Calculate it from  the max width column. with the max value over all columns being the width of that column devidied by its column_values_sum, //(int) 2**16 / width * column_values_sum;
                            auto total_kmers = std::accumulate(column_values.begin(), column_values.end(), 0);
                            auto normalization_factor = normalization_constant / total_kmers;
                            int min_col_width = std::pow(2, 16); //column_values = std::min((int) std::ceil((*column_values)*normalization_factor), min_column_width); // limit to a certain number of bits.
                            call_ones(*column, [&](uint64_t i) {
                                // OR TODO instead of looping per column, m_buffersize. Give an extra argument to the call_ones for the range (is already implemented to be possible). But then we have to keep the threads in sync.
                                auto column_value = column_values[column->rank1(i)-1];
                                column_value = std::min((int) std::ceil(column_value*normalization_factor), min_col_width); // normalize value
                                value_callback(i, column_value);
                            });
                        });
                    }
                );
            } else {
                annot::ColumnCompressed<>::merge_load(files,
                    [&](uint64_t offset, const Label &label, std::unique_ptr<bit_vector> && column) {
                        column_callback(offset, label, [&](const ValueCallback &value_callback) {
                            assert(column->num_set_bits() > 0);
                            call_ones(*column, [&](uint64_t i) {
                                value_callback(i, 1);
                            });
                        });
                    }
                );
            }
        },
        files, labels_in, labels_out, config, graph_ptr->max_index(),
        std::max(labels_in.size(), labels_out.size()), num_parallel_files,
        add_complement
    );


    auto &[counts, init_mask, other_labels, total_kmers] = count_vector;

    if (config.evaluate_assembly){     // Myrthe temporary: evaluate what k-mers overlap. Run with one thread

        logger->trace("Evaluate confusion table");
        assert(counts.size() > 1);
        assert( std::accumulate(counts.begin(), counts.end(), 0) > 0);
        int true_positive = 0; int true_negative = 0; int false_positive = 0; int false_negative =0;
        for (size_t i = 0; i < counts.size() ; i+=2){
            auto &[counts, init_mask, other_labels, total_kmers] = count_vector;
            uint64_t true_count = counts[i];
            uint64_t found_count = counts[i + 1];
            if (true_count == found_count){
                true_positive += true_count; // if 0 nothing will be added.
            }else if (true_count > found_count){
                false_negative += (true_count - found_count);
                true_positive += found_count; // add to true_positive
            }else{false_positive += (found_count - true_count);
                true_positive += true_count;} // TODO :Why are there no false positives?
        };
        // true_negatives: total number of possible k-mers? Total number of k-mers in the original foreground samples from which the result had been calculated ?
        auto precision = true_positive / ( true_positive - true_negative); // accuracy/precision = Ratio of true positives to total predicted positives.
        auto recall = true_positive / (true_positive + false_negative);  //recall / sensitivity:  ratio of true positives to total (actual) positives in the data.
        //auto specificity = true_negative / (true_negative + false_positive); // Specificity is the Ratio of true negatives to total negatives in the data.
        std::cout<< "precision" << precision << std::endl << std::flush;
        std::cout<< "recall" << recall << std::endl << std::flush;
        std::cout << std::accumulate(counts.begin(), counts.end(), 0) << std::flush;
        assert(true_positive + false_negative + false_positive == std::accumulate(counts.begin(), counts.end(), 0));
    } // in mask probably only has sequences


    sdsl::bit_vector other_mask(init_mask.size() * check_other, false);
    if (check_other && sdsl::util::cnt_one_bits(other_labels)) {
        bool parallel = num_parallel_files > 1;
        size_t j = 0;
        std::atomic_thread_fence(std::memory_order_release);
        for (size_t i = 0; i < files.size(); ++i) {
            size_t num_labels_per_file = annot::ColumnCompressed<>::read_label_encoder(files[i]).size();
            if (count_ones(other_labels, j, j + num_labels_per_file)) {
                annot::ColumnCompressed<>::merge_load({ files[i] },
                    [&](size_t offset, const auto &, auto&& column) {
                        auto &[counts, init_mask, other_labels, total_kmers] = count_vector;
                        if (other_labels[j + offset]) {
                            call_ones(init_mask, [&](size_t i) { // check how many set bits are also in other labels.
                                if ((*column)[AnnotatedDBG::graph_to_anno_index(i)])
                                    set_bit(other_mask.data(), i, parallel, MO_RELAXED);
                            });
                        }
                    },
                    num_parallel_files
                );
            }

            j += num_labels_per_file;
        }
        std::atomic_thread_fence(std::memory_order_acquire);
    }


    return mask_nodes_by_label(graph_ptr, nullptr,
                               std::move(counts), std::move(init_mask),
                               std::move(other_mask),
                               files.size(),
                               labels_in, labels_out, {}, {},
                               total_kmers,
                               config, num_threads);
}

std::shared_ptr<MaskedDeBruijnGraph>
mask_nodes_by_label(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                    const AnnotatedDBG *anno_graph,
                    sdsl::int_vector_buffer<>&& counts,
                    sdsl::bit_vector&& init_mask,
                    sdsl::bit_vector&& other_mask,
                    size_t num_labels,
                    const tsl::hopscotch_set<Label> &labels_in,
                    const tsl::hopscotch_set<Label> &labels_out,
                    const tsl::hopscotch_set<Label> &labels_in_round2,
                    const tsl::hopscotch_set<Label> &labels_out_round2,
                    std::tuple<size_t, size_t> total_kmers,
                    const DifferentialAssemblyConfig &config,
                    size_t num_threads) {
    // in and out counts are stored interleaved in the counts vector
    assert(counts.size() == init_mask.size() * 2);

    bool parallel = num_threads > 1;
    bool check_other = config.label_mask_other_unitig_fraction != 1.0;
    bool unitig_mode = check_other || labels_in_round2.size() || labels_out_round2.size()
            || config.label_mask_in_unitig_fraction != 0.0
            || config.label_mask_out_unitig_fraction != 1.0
            || config.label_mask_other_unitig_fraction != 1.0;

    size_t num_in_labels = labels_in.size() + labels_in_round2.size();
    size_t num_out_labels = labels_out.size() + labels_out_round2.size();

    bool add_complement = graph_ptr->get_mode() == DeBruijnGraph::CANONICAL
            && (config.add_complement || unitig_mode);

    auto masked_graph = make_initial_masked_graph(graph_ptr, counts, std::move(init_mask),
                                                  add_complement, num_threads);

    auto mask_or = [&](sdsl::bit_vector &a, const sdsl::bit_vector &b,
                       const std::vector<node_index> &id_map) {
        call_ones(b, [&](size_t i) {
            if (id_map[i])
                set_bit(a.data(), id_map[i], parallel, MO_RELAXED);
        });
    };

    // check all other labels and round 2 labels
    if (anno_graph && (check_other || labels_in_round2.size() || labels_out_round2.size())) {
        sdsl::bit_vector union_mask
                = static_cast<const bitmap_vector &>(masked_graph->get_mask()).data();
        std::mutex vector_backup_mutex;
        std::atomic_thread_fence(std::memory_order_release);

        auto count_merge
                = [&](sdsl::bit_vector
                              &a, // Myrthe ANSWER: what does count_merge do? Is a hack, I do not have to touch this now. Count_add would be better name. If you have a dense column, and include it in foreground or background. Wate time, Instead first create mask.  You only access dense columns, for instance a reference.
                      const sdsl::bit_vector &b, const std::vector<node_index> &id_map,
                      size_t offset = 0) {
                      call_ones(b, [&](size_t i) {
                          if (id_map[i]) {
                              set_bit(a.data(), id_map[i], parallel, MO_RELAXED);
                              { // here want no multiple threads working at once
                                  std::lock_guard<std::mutex> lock(vector_backup_mutex);
                                  counts[id_map[i] * 2 + offset] += 1;//atomic_fetch_and_add(counts, id_map[i] * 2 + offset, 1, vector_backup_mutex, MO_RELAXED);
                              }
                          }
                      });
                  };

        logger->trace("Checking shared and other labels");
        masked_graph->call_sequences(
                [&](const std::string &contig, const std::vector<node_index> &path) {
                    for (const auto &[label, sig] :
                         anno_graph->get_top_label_signatures(contig, num_labels)) {
                        bool found_in = labels_in.count(label);
                        bool found_out = labels_out.count(label);
                        bool found_in_round2 = labels_in_round2.count(label);
                        bool found_out_round2 = labels_out_round2.count(label);
                        if (!found_in && !found_out && !found_in_round2
                            && !found_out_round2 && check_other) {
                            mask_or(other_mask, sig, path);
                        }

                        if (found_in_round2)
                            count_merge(union_mask, sig, path);

                        if (found_out_round2)
                            count_merge(union_mask, sig, path, 1);
                    }
                },
                num_threads);

        std::atomic_thread_fence(std::memory_order_acquire);

        masked_graph->set_mask(new bitmap_vector(std::move(union_mask)));
    }


    // Filter unitigs from masked graph based on filtration criteria
    logger->trace("Filtering out background");
    if (config.count_kmers) {
        auto &[in_total_kmers, out_total_kmers] = total_kmers;

        if (config.test_by_unitig == false) { // Statistical testing part when k-mer counts are included.
            auto total_hypotheses = counts.size()/2; // the total number of hypotheses tested.
            auto statistical_model = DifferentialTest(config.family_wise_error_rate, total_hypotheses, config.pvalue,
                                                      in_total_kmers, out_total_kmers); // similar to the width of the counts vector, the size of the log-factorial table should be the maximum joint coverage over the in_labels resp. out_labels. Convert the width from bits to the number of stored values.
            const auto &in_mask = static_cast<const bitmap_vector &>(masked_graph->get_mask()).data();
            sdsl::bit_vector mask = in_mask;
            size_t total_nodes = masked_graph->num_nodes();
            size_t kept_nodes = 0;
            call_ones(in_mask, [&](node_index node) {
                //{ // here want no multiple threads working at once
                    //std::lock_guard<std::mutex> lock(vector_backup_mutex);
                    uint64_t in_sum = counts[node * 2];
                    uint64_t out_sum = counts[node * 2 + 1];
                //}
                if (std::get<0>(statistical_model.likelihood_ratio_test(in_sum, out_sum)))
                    ++kept_nodes;
                else
                    mask[node] = false;
            });

            masked_graph->set_mask(new bitmap_vector(std::move(mask)));

            logger->trace("Kept {} out of {} nodes", kept_nodes, total_nodes);

        } else {
            masked_graph->likelihood_ratios.resize(counts.size()/2+1);
            std::fill(masked_graph->likelihood_ratios.begin(), masked_graph->likelihood_ratios.end(), 0);
            std::atomic<uint64_t> total_unitigs(0);
            masked_graph->call_unitigs([&](const std::string &, const std::vector<node_index> &) {
                total_unitigs.fetch_add(1, MO_RELAXED);
            });
            auto statistical_model = DifferentialTest(config.family_wise_error_rate, total_unitigs, // total number of hypotheses equals the number of unitigs
                                               config.pvalue, in_total_kmers, out_total_kmers);
            { // here want no multiple threads working at once TODO luck thread.
             //   std::lock_guard<std::mutex> lock(vector_backup_mutex);
            update_masked_graph_by_unitig(*masked_graph, // intitial mask
                [&](const std::string &, const std::vector<node_index> &path) -> Intervals { // return a set of intervals to keep in the graph
                    int in_kmer_count_unitig = 0;
                    int out_kmer_count_unitig = 0;
                    for (size_t i = 0; i < path.size(); ++i) {
                            in_kmer_count_unitig += counts[path[i] * 2];
                            out_kmer_count_unitig += counts[path[i] * 2 + 1];
                    }
                    auto [significance, likelihood_ratio] = statistical_model.likelihood_ratio_test(in_kmer_count_unitig, out_kmer_count_unitig);

                    if (significance){// if significant add the likelihood ratio or the effect (in_kmer_count_unitig - out_kmer_count_unitig) to the likelihood_ratios vector. This is a temporary solution
                        masked_graph->likelihood_ratios.at(path[0]) = in_kmer_count_unitig - out_kmer_count_unitig; //likelihood_ratio;
                        masked_graph->likelihood_ratios.at(path[path.size()-1]) = in_kmer_count_unitig - out_kmer_count_unitig; //likelihood_ratio;
                        return { std::make_pair(0, path.size()) };
                    }
                    else
                        return {};
                },
                num_threads);
            }
            }
            return masked_graph;
    }


    size_t min_label_in_count = std::ceil(config.label_mask_in_kmer_fraction
                                    * num_in_labels);
    size_t max_label_out_count = std::floor(config.label_mask_out_kmer_fraction
                                    * num_out_labels);

    if (config.label_mask_in_unitig_fraction == 0.0
            && config.label_mask_out_unitig_fraction == 1.0
            && config.label_mask_other_unitig_fraction == 1.0) {
        logger->trace("Filtering by node");
        size_t total_nodes = masked_graph->num_nodes();
        const auto &in_mask = static_cast<const bitmap_vector &>(masked_graph->get_mask()).data();
        sdsl::bit_vector mask = in_mask;

        // TODO: make this part multithreaded
        size_t kept_nodes = 0;
        //{ // here want no multiple threads working at once TODO lock thread
         //   std::lock_guard<std::mutex> lock(vector_backup_mutex);
        call_ones(in_mask, [&](node_index node) { // TODO alternatively, input a mutex argument to the call_back.
                uint64_t in_count = counts[node * 2];
                uint64_t out_count = counts[node * 2 + 1];
            if (in_count >= min_label_in_count && out_count <= max_label_out_count) {
                ++kept_nodes;
            } else {
                mask[node] = false;
            }
        //}
        });



        masked_graph->set_mask(new bitmap_vector(std::move(mask)));

        logger->trace("Kept {} out of {} nodes", kept_nodes, total_nodes);

        return masked_graph;
    }

    logger->trace("Filtering by unitig");

    update_masked_graph_by_unitig(*masked_graph,
        [&](const std::string &, const std::vector<node_index> &path) -> Intervals {
            // return a set of intervals to keep in the graph
            size_t in_kmer_count = 0;

            size_t begin = path.size();
            size_t end = 0;
            for (size_t i = 0; i < path.size(); ++i) {
                if (counts[path[i] * 2] >= min_label_in_count) {
                    if (begin == path.size())
                        begin = i;

                    end = std::max(end, i + 1);

                    ++in_kmer_count;
                }
            }

            if (begin >= end)
                return {};

            size_t size = end - begin;
            size_t label_in_cutoff = std::ceil(config.label_mask_in_unitig_fraction * size);
            if (in_kmer_count < label_in_cutoff)
                return {};

            size_t out_kmer_count = 0;
            size_t other_kmer_count = 0;
            size_t label_out_cutoff = std::floor(config.label_mask_out_unitig_fraction * size);
            size_t other_cutoff = std::floor(config.label_mask_other_unitig_fraction * size);


            for (size_t i = begin; i < end; ++i) {
                if (counts[path[i] * 2 + 1] > max_label_out_count
                        && ++out_kmer_count > label_out_cutoff) {
                    return {};
                }

                if (check_other && other_mask[path[i]]
                        && ++other_kmer_count > other_cutoff) {
                    return {};
                }
            }

            return { std::make_pair(begin, end) };
        },
        num_threads
    );

    return masked_graph;
}


/**
 * Helpers
 */

std::shared_ptr<MaskedDeBruijnGraph>
make_initial_masked_graph(std::shared_ptr<const DeBruijnGraph> graph_ptr,
                          sdsl::int_vector_buffer<> &counts,
                          // TODO: Why make it a template class? int_vector_buffer and make it a template class.
                          sdsl::bit_vector&& mask,
                          bool add_complement,
                          size_t num_threads) {
    // counts is a double-length vector storing the in-label and out-label
    // counts interleaved
    assert(counts.size() == mask.size() * 2);

    auto masked_graph = std::make_shared<MaskedDeBruijnGraph>(
        graph_ptr,
        add_complement
            ? std::make_unique<bitmap_vector>(mask)
            : std::make_unique<bitmap_vector>(std::move(mask)),
        true
    );

    logger->trace("Constructed masked graph with {} nodes", masked_graph->num_nodes());

    if (add_complement) {
        logger->trace("Adding reverse complements");
        std::mutex vector_backup_mutex;
        std::atomic_thread_fence(std::memory_order_release);
        masked_graph->call_sequences([&](const std::string &seq, const std::vector<node_index> &path) {
            std::string rc_seq = seq;
            std::vector<DeBruijnGraph::node_index> rc_path = path;
            reverse_complement_seq_path(*graph_ptr, rc_seq, rc_path);

            auto it = rc_path.rbegin();
            for (size_t i = 0; i < path.size(); ++i, ++it) {
                if (*it) {
                    { // here want no multiple threads working at once
                        std::lock_guard<std::mutex> lock(vector_backup_mutex);
                        uint64_t in_count = counts[path[i] * 2]; //atomic_fetch(counts, path[i] * 2, vector_backup_mutex, MO_RELAXED);
                        uint64_t out_count = counts[path[i] * 2 + 1]; //atomic_fetch(counts, path[i] * 2 + 1, vector_backup_mutex, MO_RELAXED);
                        counts[*it * 2] += in_count; //atomic_fetch_and_add(counts, *it * 2, in_count, vector_backup_mutex, MO_RELAXED);
                        counts[*it * 2 + 1] += out_count; //atomic_fetch_and_add(counts, *it * 2 + 1, out_count, vector_backup_mutex, MO_RELAXED);
                        set_bit(mask.data(), *it, in_count, MO_RELAXED);
                    }
                }
            }
        }, num_threads, true);

        std::atomic_thread_fence(std::memory_order_acquire);

        masked_graph->set_mask(new bitmap_vector(std::move(mask)));

        logger->trace("Updated masked graph has {} nodes", masked_graph->num_nodes());
    }

    return masked_graph;
}


std::tuple<sdsl::int_vector_buffer<>, sdsl::bit_vector, sdsl::bit_vector, std::tuple<size_t, size_t>>
construct_diff_label_count_vector(const ColumnGenerator &generate_columns,
                                  const std::vector<std::string> &files,
                                  const tsl::hopscotch_set<Label> &labels_in,
                                  const tsl::hopscotch_set<Label> &labels_out,
                                  const DifferentialAssemblyConfig &config,
                                  size_t max_index,
                                  size_t num_labels,
                                  size_t num_threads,
                                  bool add_out_labels_to_mask) {
    bool parallel = (num_threads > 1);
    logger->trace("Allocating mask vector");

    size_t width;

    if (config.count_kmers){ // calculate the sum of column widths to derive an upper bound for the required counts vector width.
        assert(labels_in.size() > 0);
        size_t sum_widths_in = 0; size_t sum_widths_out = 0;
        for (size_t i = 0; i < files.size(); ++i) {
            const auto &values_fname = utils::remove_suffix(files[i],
                 annot::ColumnCompressed<>::kExtension) + annot::ColumnCompressed<>::kCountExtension;
            sdsl::int_vector_size_type col_size; uint8_t col_width;
            std::ifstream in_stream(values_fname);
            sdsl::int_vector<>::read_header(col_size, col_width, in_stream);
            auto label_encoder = annot::ColumnCompressed<>::read_label_encoder(files[i]);
            for (size_t c = 0; c < label_encoder.size(); ++c) {
                auto label = label_encoder.decode(c);
                if (labels_in.count(label)) sum_widths_in += std::pow(2, col_width);
                if (labels_out.count(label)) sum_widths_out += std::pow(2, col_width);
            }
        }
        width = std::min((size_t) 32, (size_t) std::ceil(std::log(std::max(sum_widths_in, sum_widths_out)))); // TODO, check if this is really correct. this should be 8.5 -> 9 rather then 7.
    } else{
        width = (sdsl::bits::hi(num_labels) + 1) * (1 + add_out_labels_to_mask);
    }
    assert(width);

    sdsl::bit_vector indicator(max_index + 1, false);
    size_t length = indicator.size() * 2;

    logger->trace("Allocating count vector"); // the in and out counts are stored interleaved

    std::string filename = std::tmpnam(nullptr);
    sdsl::int_vector_buffer<> counts(filename, std::ios::out, 1024*1024, width);
    for (size_t i = 0; i < length; ++i) counts.push_back(0); // set the size of the count vector using push backs

    logger->trace("done");

    sdsl::bit_vector other_labels(num_labels, false);
    std::mutex vector_backup_mutex; // a mutex variable declared before the multi-threaded part starts

    logger->trace("Populating count vector");
    std::atomic_thread_fence(std::memory_order_release);


    size_t in_total_kmers = 0; size_t out_total_kmers = 0;

    generate_columns([&](uint64_t offset, const Label &label, const ValueGenerator &value_generator) {
        uint8_t col_indicator = static_cast<bool>(labels_in.count(label));
        if (labels_out.count(label))
            col_indicator |= 2;

        if (!col_indicator)
            set_bit(other_labels.data(), offset, parallel, MO_RELAXED);

        ValueCallback add_in = [&](Column i, uint64_t value) {
            i = AnnotatedDBG::anno_to_graph_index(i); // indices in the annotation matrix have an offset of 1 compared to those in the graph. Since the mask should be compatible with the graph, we have to convert the index.
            assert(i != DeBruijnGraph::npos);
            set_bit(indicator.data(), i, parallel, MO_RELAXED); // set the corresponding bit to true in the mask, with the 'indicator' representing the mask.
            { // here want no multiple threads working at once
                std::lock_guard<std::mutex> lock(vector_backup_mutex);
                counts[i * 2] += value; //atomic_fetch_and_add(counts, i * 2, value, vector_backup_mutex, MO_RELAXED);
            }
            in_total_kmers += value;
        };

        ValueCallback add_out = [&](Column i, uint64_t value) {
            i = AnnotatedDBG::anno_to_graph_index(i);
            assert(i != DeBruijnGraph::npos);
            if (add_out_labels_to_mask)
                set_bit(indicator.data(), i, parallel, MO_RELAXED);
            { // here want no multiple threads working at once
                std::lock_guard<std::mutex> lock(vector_backup_mutex);
                counts[i * 2 + 1] += value; // atomic_fetch_and_add(counts, i * 2 + 1, value, vector_backup_mutex, MO_RELAXED);
            }
            out_total_kmers += value;
        };

        ValueCallback add_both = [&](Column i, uint64_t value) { // theoretically a label should not be contained in both the in- and out-labels, but this is just in case.
            i = AnnotatedDBG::anno_to_graph_index(i);
            assert(i != DeBruijnGraph::npos);
            set_bit(indicator.data(), i, parallel, MO_RELAXED);
            { // here want no multiple threads working at once
                std::lock_guard<std::mutex> lock(vector_backup_mutex);
                counts[i * 2] += value; // atomic_fetch_and_add(counts, i * 2, value, vector_backup_mutex, MO_RELAXED);
                counts[i * 2 + 1] += value; // atomic_fetch_and_add(counts, i * 2 + 1, value, vector_backup_mutex, MO_RELAXED);
            }
            out_total_kmers += value;
            in_total_kmers += value;
        };

        switch (col_indicator) {
            case 1: { value_generator(add_in); } break;
            case 2: { value_generator(add_out); } break;
            case 3: { value_generator(add_both); } break;
            default: {}
        }
    });

    std::atomic_thread_fence(std::memory_order_acquire);
    logger->trace("done");
    assert( std::accumulate(counts.begin(), counts.end(), 0) > 0); // assert that the sum of the count vector is greater than 0

    return std::make_tuple(std::move(counts), std::move(indicator), std::move(other_labels),
                           std::make_tuple(std::move(in_total_kmers), std::move(out_total_kmers)));
}



void kmer_distribution_table(const ColumnGenerator &generate_columns,
                                  const tsl::hopscotch_set<Label> &labels_in,
                                  const tsl::hopscotch_set<Label> &labels_out,
                                  size_t max_index,
                                  MaskedDeBruijnGraph &masked_graph, // std::shared_ptr<const DeBruijnGraph> graph_ptr,
                             const DifferentialAssemblyConfig &config,
                             size_t num_threads) {

    logger->trace("Allocating count vector");
    std::vector<std::vector<int>> count_matrix(max_index*2 + 2); // outer array 'max_index + 1'*2, with arrays of initially zero lenght or of max size in labels/out_labels.

    logger->trace("Populating count matrix");
    generate_columns([&](uint64_t offset, const Label &label, const ValueGenerator &value_generator) {
        std::ignore = offset;// to overcome the error that 'offset' is set but not used.
        uint8_t col_indicator = static_cast<bool>(labels_in.count(label));
        if (labels_out.count(label))
            col_indicator |= 2;

        ValueCallback add_in = [&](Column i, uint64_t value) {
            count_matrix[i * 2].push_back(value);
        };

        ValueCallback add_out = [&](Column i, uint64_t value) {
            count_matrix[i * 2 + 1].push_back(value);
        };

        ValueCallback add_both = [&](Column i, uint64_t value) { // theoretically a label should not be contained in both the in- and out-labels, but this is just in case.
            count_matrix[i * 2].push_back(value);
            count_matrix[i * 2 + 1].push_back(value);

        };

        switch (col_indicator) {
            case 1: { value_generator(add_in); } break;
            case 2: { value_generator(add_out); } break;
            case 3: { value_generator(add_both); } break;
            default: {}
        }

        // Calculate the product of the poisson probabilities and compare with negative binom
    });

    if (config.test_by_unitig){         // unitig mode.
        std::vector<std::vector<int>> unitig_count_matrix;
        update_masked_graph_by_unitig(masked_graph, // intitial mask
            [&](const std::string &, const std::vector<node_index> &path)
                    -> Intervals {
                // return a set of intervals to keep in the graph
                std::vector<int> in_kmer_unitig(labels_in.size(), 0);
                std::vector<int> out_kmer_unitig(labels_out.size(), 0);
                for (size_t i = 0; i < path.size(); ++i) { // first strategy is summing counts per unitig. // The update_mask.. function takes a path. Count vector per path. QUESTION: where is the count vector created, such that unitig k-mers can be accessed sequentially?
                    assert(count_matrix.size() >= path[i] * 2);
                    if (count_matrix[path[i] * 2].size() != in_kmer_unitig.size() or count_matrix[path[i] * 2 + 1].size() != out_kmer_unitig.size()){
                        std::cout << "sizes" << std::endl;
                        std::cout << std::to_string(count_matrix[path[i] * 2].size()) << std::endl << std::flush;
                        std::cout << std::to_string(in_kmer_unitig.size()) << std::endl << std::flush;
                        std::cout << std::to_string(count_matrix[path[i] * 2 + 1].size()) << std::endl << std::flush;
                        std::cout << std::to_string(out_kmer_unitig.size()) << std::endl << std::flush;
                        std::cout << std::to_string(path[i] * 2) << std::endl << std::flush; // TODO: do I have to convert the index? It was not done in the original code
                    }else{
                        assert(count_matrix[path[i] * 2].size() == in_kmer_unitig.size());
                        assert(count_matrix[path[i] * 2 + 1].size() == out_kmer_unitig.size());
                        std::transform(in_kmer_unitig.begin(), in_kmer_unitig.end(), count_matrix[path[i] * 2].begin(), in_kmer_unitig.begin(), std::plus<int>()); // sum the values of count_matrix[x] to in_kmer_count.
                        std::transform(out_kmer_unitig.begin(), out_kmer_unitig.end(), count_matrix[path[i] * 2 + 1].begin(), out_kmer_unitig.begin(), std::plus<int>());
                    }
                }
                unitig_count_matrix.push_back(in_kmer_unitig);
                unitig_count_matrix.push_back(out_kmer_unitig); // add to matrix
                return { };
                },
            num_threads);
        count_matrix = unitig_count_matrix;
    }

    int matrix_sum = 0;
    for (size_t i = 0; i < count_matrix.size(); i+=2){
        matrix_sum += std::accumulate(count_matrix[i].begin(), count_matrix[i].end(), 0);
    }
    assert(matrix_sum > 0);

    logger->trace("Writing k-mer count matrix to table");

    std::ofstream out_stream("kmer_dist_table.py");
    out_stream << "count_matrix_in = [" ;
    for (size_t i = 0; i < std::min(count_matrix.size(), (size_t) 50000000); i+=2){
        if (count_matrix[i+1].size() > 3 or count_matrix[i].size() > 3){ // labels_out.size()/3
                std::vector<int> zeros_vector(labels_in.size()-count_matrix[i].size(), 0);
                count_matrix[i].insert(count_matrix[i].end(), zeros_vector.begin(), zeros_vector.end());
                out_stream << "[" + outstring(count_matrix[i]) << "]," <<std::endl;
            }
    }
    out_stream << "]" << std::endl;
    out_stream << "count_matrix_out = [" ;
    for (size_t i = 1; i < std::min(count_matrix.size(), (size_t) 50000000); i+=2){
        if (count_matrix[i].size() > 3 or count_matrix[i-1].size() > 3){ // labels_out.size()/3
            std::vector<int> zeros_vector(labels_out.size()-count_matrix[i].size(), 0);
            count_matrix[i].insert(count_matrix[i].end(), zeros_vector.begin(), zeros_vector.end());
            out_stream << "[" + outstring(count_matrix[i]) << "]," <<std::endl;        }
    }
    out_stream << "]"<< std::endl;
    out_stream.close();
}

void update_masked_graph_by_unitig(MaskedDeBruijnGraph &masked_graph,
                                   const GetKeptIntervals &get_kept_intervals,
                                   size_t num_threads) { // traverses unitigs and takes a callback that does the filtering.
    std::atomic<uint64_t> kept_unitigs(0);
    std::atomic<uint64_t> total_unitigs(0);
    std::atomic<uint64_t> num_kept_nodes(0);
    bool parallel = num_threads > 1;

    sdsl::bit_vector mask = static_cast<const bitmap_vector&>(masked_graph.get_mask()).data();

    std::atomic_thread_fence(std::memory_order_release);

    masked_graph.call_unitigs([&](const std::string &unitig,
                                  const std::vector<node_index> &path) {
        total_unitigs.fetch_add(1, MO_RELAXED);

        size_t last = 0;
        for (const auto &pair : get_kept_intervals(unitig, path)) {
            const auto &[begin, end] = pair;
            kept_unitigs.fetch_add(1, MO_RELAXED);
            num_kept_nodes.fetch_add(end - begin, MO_RELAXED);
            for ( ; last < begin; ++last) {
                unset_bit(mask.data(), path[last], parallel, MO_RELAXED);
            }
            last = end;
        }

        for ( ; last < path.size(); ++last) {
            unset_bit(mask.data(), path[last], parallel, MO_RELAXED);
        }

    }, num_threads);
    std::atomic_thread_fence(std::memory_order_acquire);

    masked_graph.set_mask(new bitmap_vector(std::move(mask)));

    logger->trace("Kept {} out of {} unitigs with average length {} and a total of {} nodes",
                  kept_unitigs, total_unitigs,
                  static_cast<double>(num_kept_nodes + kept_unitigs * (masked_graph.get_k() - 1))
                      / kept_unitigs,
                  num_kept_nodes);
}

//%TODO A post-filtering step is done to filter out short contigs, percentage/cutoff , that is set to by default.  Short unitigs are removed afterwards
//        % TODO For each contig the confidence score is output


} // namespace graph
} // namespace mtg


