#include "dbg_aligner.hpp"

#include <string>
#include <random>

#define NUMERICAL_DIGITS "0123456789"

DBGAligner::DBGAligner(DeBruijnGraph *dbg, Annotator *annotation,
                       size_t num_top_paths, bool verbose,
                       float sw_threshold, float re_seeding_threshold,
                       float insertion_penalty, float deletion_penalty,
                       size_t num_threads) :
                            AnnotatedDBG(dbg, annotation, num_threads),
                            num_top_paths_(num_top_paths),
                            verbose_(verbose),
                            sw_threshold_(sw_threshold),
                            re_seeding_threshold_(re_seeding_threshold),
                            insertion_penalty_(insertion_penalty),
                            deletion_penalty_(deletion_penalty) {
   match_score_ = 2;
    // Substitution score for each pair of nucleotides.
    // Transition and transversion mutations have different score values.
    sub_score_ = {
        {'a', {{'a', match_score_}, {'t', -2}, {'c', -2}, {'g', -1}}},
        {'g', {{'g', match_score_}, {'t', -2}, {'c', -2}, {'a', -1}}},
        {'c', {{'c', match_score_}, {'a', -2}, {'g', -2}, {'t', -1}}},
        {'t', {{'t', match_score_}, {'a', -2}, {'g', -2}, {'c', -1}}},
        {'$', {{'$', -2}}}};

    // score for the following chars respectively:
    // A  C  G  T   N (or other ambiguous code)
    int8_t cssw_score_matrix_[5*5] = {
    sub_score_['a']['a'], sub_score_['a']['c'], sub_score_['a']['g'], sub_score_['a']['t'], sub_score_['$']['$'],
    sub_score_['c']['a'], sub_score_['c']['c'], sub_score_['c']['g'], sub_score_['c']['t'], sub_score_['$']['$'],
    sub_score_['g']['a'], sub_score_['g']['c'], sub_score_['g']['g'], sub_score_['g']['t'], sub_score_['$']['$'],
    sub_score_['t']['a'], sub_score_['t']['c'], sub_score_['t']['g'], sub_score_['t']['t'], sub_score_['$']['$'],
    sub_score_['$']['$'], sub_score_['$']['$'], sub_score_['$']['$'], sub_score_['$']['$'], sub_score_['$']['$']};

    int8_t cssw_translation_matrix_[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  //   A     C            G
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  //             T
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  //   a     c            g
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  //             t
    4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4};

    cssw_aligner_.ReBuild(cssw_score_matrix_, 5, cssw_translation_matrix_, 128);
    cssw_aligner_.SetGapPenalty(std::abs(deletion_penalty_), std::abs(deletion_penalty_));
}

namespace { // Helper function.
template <typename T>
float std_dev(const std::vector<T> &list, T mean) {
    float std_dev = 0;
    for (auto it = list.begin(); it != list.end(); ++it)
        std_dev += std::pow((*it) - mean, 2.0);
    return std::pow(std_dev / list.size(), 0.5);
}
}  // namespace

std::vector<DBGAligner::AlignedPath> DBGAligner::align(const std::string &sequence) const {
    if (sequence.size() < graph_->get_k())
        return std::vector<AlignedPath>();

    std::vector<AlignedPath> partial_paths;
    uint32_t print_frequency = 50;
    std::map<DPAlignmentKey, DPAlignmentValue> dp_alignment;
    BoundedPriorityQueue<AlignedPath> queue(num_top_paths_);
    AlignedPath initial_path(std::begin(sequence), std::begin(sequence));
    exact_map(initial_path, sequence);
    queue.push(std::move(initial_path));
    while (!queue.empty()) {
        auto path = std::move(queue.top());
        queue.pop();

        if (path.get_query_it() + graph_->get_k() > std::end(sequence)) {
            partial_paths.push_back(path);
            break;
        }

        if (path.get_total_score() < 0) {
            trim(path);
            partial_paths.push_back(path);
            queue = BoundedPriorityQueue<AlignedPath>(num_top_paths_);
            dp_alignment.clear();
            path = AlignedPath(path.get_query_it() + 1, path.get_query_it() + 1);
        }

        inexact_map(path, queue, dp_alignment, sequence);

        if (verbose_ && path.size() % print_frequency == 0) {
            auto shadow_queue(queue);
            std::cout << "Loss mean and std of paths (top path size: "
                      << path.size() << "): ";
            std::vector<float> score_values;
            score_values.reserve(queue.size());
            while (!shadow_queue.empty()) {
                score_values.push_back(shadow_queue.top().get_total_score());
                shadow_queue.pop();
            }
            float mean = std::accumulate(score_values.begin(), score_values.end(), 0);
            mean /= score_values.size();
            std::cout << mean << ", " << std_dev<float>(score_values, mean) << std::endl;
        }
    }

    if (partial_paths.empty())
        return std::vector<AlignedPath>();

    std::vector<AlignedPath> complete_path;
    for (auto partial_path = std::begin(partial_paths);
         partial_path != std::end(partial_paths); ++ partial_path) {
        // TODO: Find shortest path between the end node of a partial path
        // and the first node of the next partial path.
        complete_path.push_back(*partial_path);
    }
    // TODO: Stitch partial paths together and return the corresponding path.
    // Insertions is the parts in query that aren't matched and for deletion,
    // we have to find approximation or exact distance using labeling schemes.
    return complete_path;
}

void DBGAligner::exact_map(AlignedPath &path, const std::string &sequence) const {
    node_index seed = path.size() > 0 ? path.back() : graph_->npos;
    graph_->extend_from_seed(path.get_query_it(), std::end(sequence),
            [&](node_index node) {
            if (node == graph_->npos)
                return;
            if (path.size() > 0) {
                path.extend(node, annotator_->get(node),
                            *(path.get_query_it() + graph_->get_k() - 1), match_score_);
            } else {
                path.seed(node, annotator_->get(node), graph_->get_node_sequence(node),
                          match_score_ * graph_->get_k());
            }},
            [&]() { return (path.size() > 0 &&
                    (graph_->outdegree(path.back()) > 1)); },
            seed);
}

void DBGAligner::inexact_map(AlignedPath& path,
                             BoundedPriorityQueue<AlignedPath> &queue,
                             std::map<DPAlignmentKey, DPAlignmentValue> &dp_alignment,
                             const std::string &sequence) const {
    if (path.get_query_it() + graph_->get_k() > std::end(sequence))
        return;
    if (path.size() == 0) {
        // To do inexact seeding, a path with incremented query_it is pushed to the queue.
        path.set_query_begin_it(path.get_query_it() + 1);
        path.set_query_it(path.get_query_it() + 1);
        exact_map(path, sequence);
        queue.push(std::move(path));
        return;
    }
    graph_->call_outgoing_kmers(path.back(), [&](node_index node, char extension) {
        AlignedPath alternative_path(path);
        alternative_path.set_query_it(path.get_query_it());

        if (queue.size() > 0 && queue.back() < alternative_path) {
            alternative_path.extend(node, annotator_->get(node), extension,
                single_char_score(*(alternative_path.get_query_it() +
                    graph_->get_k() - 1), extension));
        } else {
            alternative_path.extend(node, annotator_->get(node), extension);
            alternative_path.update_total_score(
                ssw_score(alternative_path));
        }
        exact_map(alternative_path, sequence);

        DPAlignmentKey alternative_key{.node = alternative_path.back(),
                           .query_begin_it = alternative_path.get_query_begin_it(),
                           .query_it = alternative_path.get_query_it()};
        DPAlignmentValue alternative_value = {.score = alternative_path.get_total_score()};
        auto dp_alignment_it = dp_alignment.find(alternative_key);
        if (dp_alignment_it == dp_alignment.end()) {
            dp_alignment[alternative_key] = alternative_value;
            queue.push(std::move(alternative_path));
        }
        else if (dp_alignment_it->second.score < alternative_path.get_total_score()) {
            dp_alignment_it->second = alternative_value;
            queue.push(std::move(alternative_path));
        }
        });
}

float DBGAligner::single_char_score(char char_in_query, char char_in_graph) const {
    try {
        return sub_score_.at(std::tolower(char_in_query)).at(std::tolower(char_in_graph));
    }
    catch (const std::out_of_range&) {
        return sub_score_.at('$').at('$');
    }
}

namespace {
// helper function to print an alignment based on SSW library.
static void PrintAlignment(const StripedSmithWaterman::Alignment &alignment){
    std::cout << "===== SSW result =====" << std::endl;
    std::cout << "Best Smith-Waterman score:\t" << alignment.sw_score << std::endl
       << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << std::endl
       << "Reference start:\t" << alignment.ref_begin << std::endl
       << "Reference end:\t" << alignment.ref_end << std::endl
       << "Query start:\t" << alignment.query_begin << std::endl
       << "Query end:\t" << alignment.query_end << std::endl
       << "Next-best reference end:\t" << alignment.ref_end_next_best << std::endl
       << "Number of mismatches:\t" << alignment.mismatches << std::endl
       << "Cigar: " << alignment.cigar_string << std::endl;
    std::cout << "======================" << std::endl;
}
}  // namespace

bool DBGAligner::cssw_align(const AlignedPath &path, StripedSmithWaterman::Alignment& alignment) const {
    auto ref = path.get_sequence();
    std::string query(path.get_query_begin_it(), path.get_query_it() + graph_->get_k() - 1);

    int32_t maskLen = strlen(query.c_str()) / 2;
    maskLen = maskLen < 15 ? 15 : maskLen;

    StripedSmithWaterman::Filter filter;
    return cssw_aligner_.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);
}

float DBGAligner::ssw_score(const AlignedPath &path) const {
    StripedSmithWaterman::Alignment alignment;
    if (verbose_) {
        PrintAlignment(alignment);
    }
    if (!cssw_align(path, alignment)) {
        std::cout << "Failure in SSW calculation. Computing SSW inefficiently.\n";
        return whole_path_score(path);
    }
    std::string query(path.get_query_begin_it(), path.get_query_it() + graph_->get_k() - 1);

    // Calculate mismatch score for either X or S in cigar
    // (mismatched or not counted due to clipping).
    // TODO: Count mismatches more accurately based on different values in sub_score_.
    // TODO: Count insertion, deletions based on cigar.
    float score = sub_score_.at('$').at('$') * alignment.mismatches
                  + match_score_ *
                    (alignment.query_end - alignment.query_begin + 1.0)
                  + sub_score_.at('$').at('$') *
                    (query.size() - 1.0 - alignment.query_end);

    return score;
}

float DBGAligner::whole_path_score(const AlignedPath &path) const {
    auto path_sequence = path.get_sequence();
    std::string query_sequence(path.get_query_begin_it(), path.get_query_it() + graph_->get_k() - 1);
    float alignment_score[query_sequence.size() + 1][path_sequence.size() + 1] = {};
    for (size_t i = 1; i <= query_sequence.size(); ++i) {
        for (size_t j = 1; j <= path_sequence.size(); ++j) {
            alignment_score[i][j] = std::max({alignment_score[i-1][j-1] +
                    single_char_score(query_sequence[i-1], path_sequence[j-1]),
                alignment_score[i][j-1] - insertion_penalty_, alignment_score[i-1][j] - deletion_penalty_});
        }
    }
    return alignment_score[query_sequence.size()][path_sequence.size()];
}

void DBGAligner::trim(AlignedPath &path) const {
    StripedSmithWaterman::Alignment alignment;
    if (verbose_) {
        PrintAlignment(alignment);
    }
    if (!cssw_align(path, alignment)) {
        std::cout << "Failure in SSW calculation. Computing SSW inefficiently.\n";
        return;
    }
    // Returns in case of no clipping in CSSW cigar.
    if (alignment.cigar_string.back() != 'S')
        return;
    auto num_starting_pos = alignment.cigar_string.find_last_not_of(NUMERICAL_DIGITS,
                                alignment.cigar_string.size() - 2) + 1;
    auto trim_length = std::atoi(alignment.cigar_string.substr(num_starting_pos,
                                 alignment.cigar_string.size() - num_starting_pos - 1).c_str());
    path.trim_with_score(trim_length, alignment.sw_score);
    path.set_cigar(alignment.cigar_string);
}
