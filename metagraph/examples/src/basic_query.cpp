/**
 * This example shows how to:
 * 1. Load a de Bruijn graph from disk
 * 2. Load annotations for the graph
 * 3. Query sequences against the annotated graph
 * 4. Process and display results
 *
 * Prerequisites:
 * Before running this example, you need a graph and annotation.
 * See the README.md for instructions on building test data.
 *
 * Usage:
 *   ./basic_query <graph.dbg> <annotation.annodbg> <query.fa>
 *
 * Example:
 *   ./basic_query ../data/graphs/test_graph.dbg \
 *                 ../data/graphs/test_graph.column.annodbg \
 *                 ../data/test_query.fa
 */

#include <iostream>
#include <memory>
#include <string>
#include <map>

// IMPORTANT: When using ExternalProject_Add, metagraph.hpp should be the FIRST
// MetaGraph header included. It contains compile-time configuration definitions
// that affect how other MetaGraph headers are interpreted.
#include "metagraph.hpp"
#include "seq_io/sequence_io.hpp"
#include "cli/load/load_annotation.hpp"

using namespace mtg;
using namespace mtg::cli;
using namespace mtg::graph;

void print_usage(const char* program_name) {
    std::cerr << "Usage: " << program_name 
              << " <graph.dbg> <annotation.annodbg> <query.fa>\n"
              << "\n"
              << "Example:\n"
              << "  " << program_name 
              << " graph.dbg graph.column.annodbg query.fa\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        print_usage(argv[0]);
        return 1;
    }

    std::string graph_path = argv[1];
    std::string annotation_path = argv[2];
    std::string query_path = argv[3];

    std::cout << "Loading graph from: " << graph_path << std::endl;
    
    // Load the de Bruijn graph
    auto graph = load_critical_dbg(graph_path);
    if (!graph) {
        std::cerr << "Error: Failed to load graph from " << graph_path << std::endl;
        return 1;
    }
    
    std::cout << "Graph loaded successfully" << std::endl;
    std::cout << "  k: " << graph->get_k() << std::endl;
    std::cout << "  Number of nodes: " << graph->num_nodes() << std::endl;
    
    std::cout << "\nLoading annotation from: " << annotation_path << std::endl;
    
    // Initialize annotation using parse_annotation_type and initialize_annotation
    auto anno_type = parse_annotation_type(annotation_path);
    auto annotation = initialize_annotation(anno_type);
    
    if (!annotation->load(annotation_path)) {
        std::cerr << "Error: Failed to load annotations from " << annotation_path << std::endl;
        return 1;
    }
    
    // Create an AnnotatedDBG
    auto anno_graph = std::make_unique<AnnotatedDBG>(graph, std::move(annotation));
    
    std::cout << "Annotation loaded successfully" << std::endl;
    std::cout << "  Number of labels: " << anno_graph->get_annotator().num_labels() << std::endl;
    
    // Display some label names (first few)
    const auto& label_encoder = anno_graph->get_annotator().get_label_encoder();
    size_t num_labels_to_show = std::min(size_t(5), label_encoder.size());
    if (num_labels_to_show > 0) {
        std::cout << "  Sample labels:" << std::endl;
        for (size_t i = 0; i < num_labels_to_show; ++i) {
            std::cout << "    - " << label_encoder.decode(i) << std::endl;
        }
    }
    
    // Now perform a simple query using the lower-level API
    std::cout << "\nQuerying sequences from: " << query_path << std::endl;
    
    size_t num_sequences = 0;
    double min_kmers_fraction = 0.7;  // Require 70% of k-mers to match
    
    seq_io::read_fasta_file_critical(query_path, 
        [&](seq_io::kseq_t* read_stream) {
            num_sequences++;
            
            std::string header = read_stream->name.s;
            std::string sequence = read_stream->seq.s;
            
            std::cout << "\nSequence: " << header << std::endl;
            std::cout << "  Length: " << sequence.length() << " bp" << std::endl;
            
            // Map the sequence to the graph
            std::vector<DeBruijnGraph::node_index> mapping;
            anno_graph->get_graph().map_to_nodes_sequentially(sequence,
                [&](DeBruijnGraph::node_index node) { mapping.push_back(node); });
            
            // Count how many k-mers matched
            size_t total_kmers = mapping.size();
            size_t matched_kmers = 0;
            for (const auto& node : mapping) {
                if (node != DeBruijnGraph::npos) {
                    matched_kmers++;
                }
            }
            
            double match_fraction = total_kmers > 0 
                ? static_cast<double>(matched_kmers) / total_kmers 
                : 0.0;
            
            std::cout << "  K-mers matched: " << matched_kmers << "/" << total_kmers 
                      << " (" << (match_fraction * 100) << "%)" << std::endl;
            
            if (match_fraction >= min_kmers_fraction && matched_kmers > 0) {
                // Get annotations for matched nodes
                std::map<std::string, size_t> label_counts;
                
                for (const auto& node : mapping) {
                    if (node != DeBruijnGraph::npos) {
                        // Convert graph node index to annotation index
                        auto anno_index = AnnotatedDBG::graph_to_anno_index(node);
                        auto labels = anno_graph->get_annotator().get_labels(anno_index);
                        for (const auto& label : labels) {
                            label_counts[label]++;
                        }
                    }
                }
                
                if (!label_counts.empty()) {
                    std::cout << "  Matching labels:" << std::endl;
                    for (const auto& [label, count] : label_counts) {
                        std::cout << "    - " << label << " (in " << count << " k-mers)" << std::endl;
                    }
                } else {
                    std::cout << "  No label annotations found" << std::endl;
                }
            } else {
                std::cout << "  Match fraction below threshold (" 
                          << (min_kmers_fraction * 100) << "%)" << std::endl;
            }
        }
    );
    
    std::cout << "\nExample completed successfully!" << std::endl;
    return 0;
}
