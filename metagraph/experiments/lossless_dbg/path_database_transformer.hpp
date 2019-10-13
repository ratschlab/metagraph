//
// Created by Jan Studený on 2019-05-20.
//

#ifndef __PATH_DATABASE_TRANSFORMER_HPP__
#define __PATH_DATABASE_TRANSFORMER_HPP__

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

#endif // __PATH_DATABASE_TRANSFORMER_HPP__
