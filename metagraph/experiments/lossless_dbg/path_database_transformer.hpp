//
// Created by Jan Studený on 2019-05-20.
//

#ifndef METAGRAPH_PATH_DATABASE_TRANSFORMER_HPP
#define METAGRAPH_PATH_DATABASE_TRANSFORMER_HPP

#include "path_database_dynamic.hpp"
#include "path_database_wavelet.hpp"

class PathDatabaseTransformer {
public:
    PathDatabaseDynamicCore<>& pd;

    void freeze(PathDatabaseWavelet<>& frozen, PathDatabaseDynamicCore<>& ) {

    }
    IncomingTable<>

};

PathDatabaseWavelet<> freeze(PathDatabaseDynamicCore<>& pd) {

}

#endif //METAGRAPH_PATH_DATABASE_TRANSFORMER_HPP
