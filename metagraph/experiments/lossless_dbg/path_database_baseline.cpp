//
//  path_database_baseline.cpp
//  PathDatabase
//
//  Created by Jan Studený on 21/03/2019.
//

#include "path_database_baseline.hpp"

PathDatabaseBaseline::PathDatabaseBaseline(std::shared_ptr<const DeBruijnGraph> graph) : PathDatabase(graph), graph(*graph_) {

}

PathDatabaseBaseline::PathDatabaseBaseline(const vector<string> &raw_reads, size_t k_kmer) : PathDatabase(raw_reads,k_kmer), graph(*graph_) {}
