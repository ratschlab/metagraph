#include "benchmark_graph_helpers.hpp"

#include <vector>
#include <string>

#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"
#include "cli/query.hpp"
#include "common/utils/string_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "seq_io/sequence_io.hpp"


namespace mg {

std::shared_ptr<DeBruijnGraph> bm::build_graph(const std::string &filename) {
    std::vector<std::string> sequences;
    read_fasta_file_critical(filename,
                             [&](kseq_t *stream) {
                                 sequences.emplace_back(stream->seq.s);
                             },
                             true);

    size_t k = 12;

    BOSSConstructor constructor(k - 1);
    constructor.add_sequences(std::move(sequences));
    auto graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor));
    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);

    return graph;
}

template <class Annotation>
std::unique_ptr<AnnotatedDBG> bm::build_anno_graph(const std::string &filename) {
    std::vector<std::string> sequences;
    std::vector<std::string> labels;
    read_fasta_file_critical(filename,
                             [&](kseq_t *stream) {
                                 for (const auto &label : utils::split_string(stream->name.s, "|")) {
                                     sequences.emplace_back(stream->seq.s);
                                     labels.push_back(label);
                                 }
                             },
                             true);

    size_t k = 12;

    BOSSConstructor constructor(k - 1);
    constructor.add_sequences(std::vector<std::string>(sequences));
    auto graph = std::make_shared<DBGSuccinct>(new BOSS(&constructor));
    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);

    uint64_t max_index = graph->max_index();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<Annotation>(max_index)
    );

    for (size_t i = 0; i < sequences.size(); ++i) {
        anno_graph->annotate_sequence(std::string(sequences[i]), { labels[i] });
    }

    return anno_graph;
}

template std::unique_ptr<AnnotatedDBG> bm::build_anno_graph<annotate::ColumnCompressed<>>(const std::string &filename);
template std::unique_ptr<AnnotatedDBG> bm::build_anno_graph<annotate::RowCompressed<>>(const std::string &filename);


std::unique_ptr<AnnotatedDBG> bm::build_query_graph(const AnnotatedDBG &anno_graph,
                                                    const std::string &query_filename) {
    return construct_query_graph(
        anno_graph,
        [&](std::function<void(const std::string&)> call_sequence) {
            read_fasta_file_critical(
                query_filename,
                [&](kseq_t *stream) { call_sequence(stream->seq.s); }
            );
        },
        0.0,
        1
    );
}

} //namespace mg
