#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>


uint64_t kFromFile(std::string infbase) {
    uint64_t k = 0;
    std::ifstream instream((infbase + ".F.dbg").c_str()); 
    std::string line;
    size_t mode = 0;
    while (std::getline(instream, line)) {
        if (strcmp(line.c_str(), ">F") == 0 || strcmp(line.c_str(), ">p") == 0) {
            mode = 1;
        } else if (strcmp(line.c_str(), ">k") == 0) {
            mode = 2;
        } else {
            if (mode == 2) {
                k = strtoul(line.c_str(), NULL, 10);
            } else if (mode == 1) {
                continue;
            } else {
                fprintf(stderr, "ERROR: input file corrupted\n");
                exit(1);
            }
        }
    }
    instream.close();
    return k;
};



#endif
