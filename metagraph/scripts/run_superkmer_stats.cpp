#include "graph/representation/hash/dbg_sshash.hpp"
#include "graph/annotated_dbg.hpp"
#include "annotation/representation/base/annotation.hpp"

#include "cli/load/load_annotation.hpp"
#include "annotation/representation/column_compressed/annotate_column_compressed.hpp"



int main (){
    
    using namespace mtg::graph;
    using namespace mtg::annot;
    using namespace mtg::annot::matrix;
    using namespace mtg::cli;

    std::string graph_path = "/home/marianna/Documents/Masterthesis/new_clone/metagraph/metagraph/build/graph_subset_750.sshashdbg";
    std::string anno_path = "/home/marianna/Documents/Masterthesis/new_clone/metagraph/metagraph/build/graph_subset_750_annotation.column.annodbg";

    std::shared_ptr<mtg::graph::DBGSSHash> graph_ptr = std::make_shared<mtg::graph::DBGSSHash>(31);
    graph_ptr->load(graph_path);
    
    Config::AnnotationType config = parse_annotation_type(anno_path);
    /*
    Types of annotations:
    ColumnCompressed = 1,
        RowCompressed,
        BRWT,
        BinRelWT,
        RowDiff,
        RowDiffBRWT,
        RowDiffRowFlat,
        RowDiffRowSparse,
        RowDiffDisk,
        RowFlat,
        RowSparse,
        RBFish,
        RbBRWT,
        IntBRWT,
        IntRowDiffBRWT,
        IntRowDiffDisk,
        ColumnCoord,
        BRWTCoord,
        RowDiffCoord,
        RowDiffBRWTCoord,
        RowDiffDiskCoord,
        */
    assert(config == Config::AnnotationType::ColumnCompressed);
    std::unique_ptr<ColumnCompressed<>> anno_ptr =  std::make_unique<ColumnCompressed<>> (); 
    anno_ptr->load(anno_path);

    std::unique_ptr<AnnotatedDBG> anno_graph = std::make_unique<AnnotatedDBG>(graph_ptr, std::move(anno_ptr), false);
    graph_ptr->superkmer_statistics(anno_graph);

    return 0;
}
