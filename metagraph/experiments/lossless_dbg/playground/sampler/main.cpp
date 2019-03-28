//
// Created by Jan Studený on 2019-03-11.
//
#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <tclap/CmdLine.h>
#include <random>

using TCLAP::ValueArg;
using TCLAP::MultiArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::UnlabeledMultiArg;
using TCLAP::ValuesConstraint;


#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#pragma clang diagnostic ignored "-Wcomma"

#include "samplers.hpp"
#include "utilities.hpp"

#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"



using namespace std;
using namespace std::string_literals;



int main(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Compress reads",' ', "0.1");
    TCLAP::ValueArg<std::string> referenceArg("r",
                                         "reference",
                                         "Path to the reference from which the sampler will sample (FASTA file)",
                                         true,
                                         "",
                                         "string");
    TCLAP::ValueArg<double>          coverageArg("c",
                                              "coverage",
                                              "Coverage of reads",
                                              false,
                                              1.0,
                                              "double");
    TCLAP::ValueArg<int>          read_lengthArg("l",
                                                 "read-length",
                                                 "Length of individual reads",
                                                 false,
                                                 100,
                                                 "int");
    TCLAP::ValueArg<int>          subsampleArg("z",
                                                "subsample-size",
                                                "Subsample the reference to 'subsample-size' base pairs.",
                                                false,
                                                -1,
                                                "int");
    TCLAP::ValueArg<int>          seedArg("s",
                                                 "seed",
                                                 "Seed for a random generator",
                                                 false,
                                                 12345,
                                                 "int");
    TCLAP::ValueArg<std::string> outputArg("o",
                                          "output",
                                          "Filename where to store sampled reads.",
                                          false,
                                          "output.fasta",
                                          "string");
    cmd.add(referenceArg);
    cmd.add(coverageArg);
    cmd.add(outputArg);
    cmd.add(read_lengthArg);
    cmd.add(subsampleArg);
    cmd.parse(argc, argv);
    auto reference_filename = referenceArg.getValue();
    auto coverage = coverageArg.getValue();
    auto output_filename = outputArg.getValue();
    auto read_length = read_lengthArg.getValue();
    auto seed = seedArg.getValue();
    auto subsample_size = subsampleArg.getValue();
    auto generator = std::mt19937(seed);
    auto reference = read_reads_from_fasta(reference_filename)[0];
    Sampler sampler = subsampleArg.isSet() ? Sampler(reference,generator) :
                                             SubSampler(reference,subsample_size,generator);
    auto reads = sampler.sample_coverage(read_length,coverage);
    write_reads_to_fasta(reads,output_filename);
    return 0;
}
