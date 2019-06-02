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


int main_sampler(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Compress reads",' ', "0.1");
    TCLAP::ValueArg<std::string> referenceArg("r",
                                         "reference",
                                         "Path to the reference from which the sampler will sample (FASTA file)",
                                         true,
                                         "",
                                         "string",cmd);
    TCLAP::ValueArg<double>          coverageArg("c",
                                              "coverage",
                                              "Coverage of reads",
                                              false,
                                              1.0,
                                              "double",cmd);
    TCLAP::ValueArg<double>          error_probabilityArg("p",
                                                 "error-probability",
                                                 "Probability of a substitution error",
                                                 false,
                                                 0.0,
                                                 "double",cmd);
    TCLAP::ValueArg<int>          read_lengthArg("l",
                                                 "read-length",
                                                 "Length of individual reads",
                                                 false,
                                                 100,
                                                 "int64_t",cmd);
    TCLAP::ValueArg<int>          subsampleArg("z",
                                                "subsample-size",
                                                "Subsample the reference to 'subsample-size' base pairs.",
                                                false,
                                                -1,
                                                "int64_t",cmd);
    TCLAP::ValueArg<int>          seedArg("s",
                                                 "seed",
                                                 "Seed for a random generator",
                                                 false,
                                                 12345,
                                                 "int64_t",cmd);
    TCLAP::ValueArg<std::string> outputArg("o",
                                          "output",
                                          "Filename where to store sampled reads.",
                                          false,
                                          "output.fasta",
                                          "string",cmd);
    cmd.parse(argc, argv);
    auto reference_filename = referenceArg.getValue();
    auto coverage = coverageArg.getValue();
    auto output_filename = outputArg.getValue();
    auto read_length = read_lengthArg.getValue();
    auto error_probability = error_probabilityArg.getValue();
    auto seed = seedArg.getValue();
    auto subsample_size = subsampleArg.getValue();
    auto generator = std::mt19937(seed);
    auto reference = read_reads_from_fasta(reference_filename)[0];
    NoisySampler sampler = subsampleArg.isSet() ?  SubSampler(reference,subsample_size,generator,error_probability) :
                                              NoisySampler(reference,generator,error_probability);
    auto reads = sampler.sample_coverage(read_length,coverage);
    write_reads_to_fasta(reads,output_filename);
    return 0;
}
