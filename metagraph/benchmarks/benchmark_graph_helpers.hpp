#include <vector>
#include <string>

#include "annotation/annotation_converters.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "common/utils/string_utils.hpp"
#include "graph/annotated_dbg.hpp"
#include "graph/representation/succinct/dbg_succinct.hpp"
#include "graph/representation/succinct/boss_construct.hpp"
#include "seq_io/sequence_io.hpp"

std::unique_ptr<AnnotatedDBG> build_anno_graph(const std::string &filename) {
    std::vector<std::string> sequences;
    std::vector<std::string> labels;
    read_fasta_file_critical(filename,
                             [&](kseq_t *stream) {
                                 std::string name(stream->name.s);
                                 for (const auto &label : utils::split_string(name, "|")) {
                                     sequences.push_back(stream->seq.s);
                                     labels.push_back(label);
                                 }
                             },
                             true);

    size_t k = 12;

    BOSSConstructor constructor(k - 1);
    constructor.add_sequences(sequences);
    std::shared_ptr<DeBruijnGraph> graph { new DBGSuccinct(new BOSS(&constructor)) };
    dynamic_cast<DBGSuccinct*>(graph.get())->mask_dummy_kmers(1, false);

    uint64_t max_index = graph->max_index();

    auto anno_graph = std::make_unique<AnnotatedDBG>(
        graph,
        std::make_unique<annotate::ColumnCompressed<>>(max_index)
    );

    for (size_t i = 0; i < sequences.size(); ++i) {
        anno_graph->annotate_sequence(std::string(sequences[i]), { labels[i] });
    }

    return anno_graph;
}
