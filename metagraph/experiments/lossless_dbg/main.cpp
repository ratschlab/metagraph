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

#include "playground/sampler/sampler.hpp"
#include "playground/compressor/compressor.hpp"
#include "playground/decompressor/decompressor.hpp"
#include "playground/statistics/statistics.hpp"
#include "playground/tester/tester.hpp"
#include "playground/tidbits/dump_routing_table.hpp"

int main(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Toolbox for experiments", ' ', "");

    std::vector<std::string> regimes {
            "sample",
            "compress",
            "decompress",
            "statistics",
            "identity",
            "dumprt",
    };

    ValuesConstraint<std::string> regime_constraint(regimes);
    UnlabeledValueArg<std::string> regime_arg("regime", "Regime", true, "", &regime_constraint, cmd);
    cmd.parse(std::min(argc, 2), argv);

    std::vector<char*> new_argv(argc - 1);
    new_argv[0] = argv[0];
    for (int i = 2; i < argc; i++) {
        new_argv[i - 1] = argv[i];
    }

    std::string regime = regime_arg.getValue();
    regime_arg.reset();

    if (regime == "sample") {
        main_sampler(new_argv.size(), new_argv.data());
    }
    else if (regime == "compress") {
        main_compressor(new_argv.size(), new_argv.data());
    }
    else if (regime == "decompress") {
        main_decompressor(new_argv.size(), new_argv.data());
    }
    else if (regime == "statistics") {
        main_statistics(new_argv.size(), new_argv.data());
    }
    else if (regime == "identity") {
        main_tester(new_argv.size(), new_argv.data());
    }
    else if (regime == "dumprt") {
        main_dumprt(new_argv.size(),new_argv.data());
    }
    else if (regime == "test") {
       // main_compressor()
    }
}
