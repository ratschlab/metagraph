#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <tclap/CmdLine.h>
#include <gtest/gtest.h>

using TCLAP::ValueArg;
using TCLAP::MultiArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::UnlabeledMultiArg;
using TCLAP::ValuesConstraint;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#pragma clang diagnostic ignored "-Wcomma"



#define main main_sampler
#include "playground/sampler/main.cpp"
#undef main

#define main main_compressor
#include "playground/compressor/main.cpp"
#undef main

#define main main_decompressor
#include "playground/decompressor/main.cpp"
#undef main

#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"


int main(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Toolbox for experiments", ' ', "");

    std::vector<std::string> regimes {
            "sample",
            "compress",
            "decompress"
    };

    ValuesConstraint<std::string> regime_constraint(regimes);
    UnlabeledValueArg<std::string> regime_arg("regime", "Regime", true, "", &regime_constraint, cmd);
    cmd.parse(std::min(argc, 2), argv);

    char ** new_argv = (char**) malloc((argc-1)*sizeof(char*));
    new_argv[0] = argv[0];
    for(int i=2;i<argc;i++) {
        new_argv[i-1] = argv[i];
    }

    std::string regime = regime_arg.getValue();
    regime_arg.reset();

    if (regime == "sample") {
        main_sampler(argc-1,new_argv);
    }
    else if (regime == "compress") {
        main_compressor(argc-1,new_argv);
    }
    else if (regime == "decompress") {
        main_decompressor(argc-1,new_argv);
    }
}
