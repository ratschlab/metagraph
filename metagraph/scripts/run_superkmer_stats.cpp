#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/representation/base/annotation.hpp"

#include "cli/load/load_annotation.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"
#include "annotation/representation/row_compressed/annotate_row_compressed.hpp"

#include "string_utils.hpp"

#include "annotation/representation/annotation_matrix/static_annotators_def.hpp"
#include "annotation/binary_matrix/base/binary_matrix.hpp"

int main (int argc, char *argv[]){
    
    using namespace mtg::graph;
    using namespace mtg::annot;
    using namespace mtg::annot::matrix;
    using namespace mtg::cli;

    if(argc <3){
    	std::cerr<<"missing input files!\n" ;//<< "graph path, annotation path\n";
	return EXIT_FAILURE;
    }
    std::string graph_path = argv[1];
    std::string anno_path = argv[2];

    std::shared_ptr<mtg::graph::DBGSSHash> graph_ptr = std::make_shared<mtg::graph::DBGSSHash>(31);
    graph_ptr->load(graph_path);
    
    std::string sk_mask_path = utils::remove_suffix(graph_path, graph_ptr->kExtension) + "_sk_mask";

    // Warning: Unused variable...
    //Config::AnnotationType config = parse_annotation_type(anno_path);
    //assert(config == Config::AnnotationType::RowFlat);

    std::unique_ptr<RowFlatAnnotator> anno_ptr = std::make_unique<RowFlatAnnotator>();
    anno_ptr->load(anno_path);
    std::unique_ptr<AnnotatedDBG> anno_graph = std::make_unique<AnnotatedDBG>(graph_ptr, std::move(anno_ptr), false);
    
    //graph_ptr->superkmer_stats(anno_graph);
    graph_ptr->superkmer_bv(anno_graph, sk_mask_path);

    return 0;
}