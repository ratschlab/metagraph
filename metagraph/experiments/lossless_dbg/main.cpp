#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <ProgressBar.hpp>
#include <tclap/CmdLine.h>
#include <gtest/gtest.h>

using TCLAP::ValueArg;
using TCLAP::MultiArg;
using TCLAP::UnlabeledValueArg;
using TCLAP::UnlabeledMultiArg;
using TCLAP::ValuesConstraint;

using json = nlohmann::json;
#define _DNA_GRAPH 1

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#pragma clang diagnostic ignored "-Wcomma"

#include "dbg_succinct.hpp"
#include "sequence_graph.hpp"
#include "dbg_succinct_construct.hpp"
#include "dbg_hash.hpp"


#include "compressed_reads.hpp"
#include "samplers.hpp"
#include "utils.hpp"
#include "tests.hpp"

#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"



using namespace std;
using namespace std::string_literals;
namespace fs = std::filesystem;

// debugging functions
struct d_t {
    template<typename T> d_t & operator,(const T & x) {
        std::cerr << ' ' <<  x;
        return *this;
    }
};

#define D(args ...) { d_t, "|", __LINE__, "|", #args, ":", args, "\n"; }

using node_index = SequenceGraph::node_index;

string local_file(string filename) {
    return fs::path(__FILE__).parent_path() / filename;
}

// TODO: remove these constants
string HUMAN_REFERENCE_FILENAME = local_file("genomic_data/GCF_000001405.38_GRCh38.p12_genomic.fna");
string HUMAN_CHROMOSOME_10_SAMPLE = local_file("genomic_data/human_chromosome_10_sample.fasta");
string HUMAN_CHROMOSOME_10_FILENAME = local_file("genomic_data/human_chromosome_10.fasta");
string HUMAN_CHROMOSOME_10_STRIPPED_N_FILENAME = local_file("genomic_data/human_chromosome_10_n_trimmed.fasta");
string JSON_OUTPUT_FILE = local_file("statistics.json");
const int CHROMOSOME_NUMBER = 10;
const int READ_LENGTH = 100;
const double READ_COVERAGE = 0.00001;
const int test_seed = 3424;
const int DEFAULT_K_KMER = 21;
// todo change to proper thing
#define x first
#define y second
#define all(x) begin(x),end(x)

string get_human_chromosome(int chromosome_number,bool five_letter_alphabet=true) {
    int current_chromosome=1;
    string result;
    read_fasta_file_critical(
                             HUMAN_REFERENCE_FILENAME,
                             [&](kseq_t* chromosome) {
                                 if (chromosome->comment.s == ("Homo sapiens chromosome "s + to_string(chromosome_number) + ", GRCh38.p12 Primary Assembly"s)) {
                                     result = chromosome->seq.s;
                                 }
                                 current_chromosome++;
                             });
    if (five_letter_alphabet) {
        transform(all(result),result.begin(),::toupper);
    }
    return result;
}

void to_be_determined() {
    auto chromosome = get_human_chromosome(CHROMOSOME_NUMBER);
    auto sampler = Sampler(chromosome,test_seed);
    auto reads = sampler.sample_coverage(READ_LENGTH, READ_COVERAGE);
    transform_to_fasta(HUMAN_CHROMOSOME_10_SAMPLE,reads);
    auto compressed_reads = CompressedReads(reads);
}
void code_to_violate_assertion() {
    auto reads = read_reads_from_fasta(HUMAN_CHROMOSOME_10_SAMPLE);
    auto compressed_reads = CompressedReads(reads);
}
void save_human_chromosome() {
    auto chromosome = get_human_chromosome(CHROMOSOME_NUMBER);
    transform_to_fasta(HUMAN_CHROMOSOME_10_FILENAME,{chromosome});
}

void get_statistics() {
    auto chromosome_cleaned
        = read_reads_from_fasta(HUMAN_CHROMOSOME_10_STRIPPED_N_FILENAME)[0];

    DBGHash graph(DEFAULT_K_KMER);

    graph.add_sequence(chromosome_cleaned);

    int kmers_count = graph.num_nodes();

    auto pb = ProgressBar(kmers_count, 70, '=', ' ', 10000);

    map<int,int> kmer_outgoing_edges_statistics;
    
    vector<node_index> outgoing_edges;

    #pragma omp parallel for reduction(map_reduction:kmer_outgoing_edges_statistics)
    for (size_t i = 1; i <= graph.num_nodes(); ++i) {
        // const auto &kmer = graph.get_node_sequence(i);
        graph.adjacent_outgoing_nodes(i, &outgoing_edges);

        kmer_outgoing_edges_statistics[outgoing_edges.size()]++;

        outgoing_edges.clear();

        ++pb;
        pb.display();
    }

    json statistics(kmer_outgoing_edges_statistics);

    ofstream myfile;

    myfile.open(JSON_OUTPUT_FILE);

    myfile << statistics.dump(4) << endl;
    cout << statistics.dump(4) << endl;
}


int main(int argc, char *argv[]) {
    TCLAP::CmdLine cmd("Compress reads",' ', "0.1");
    TCLAP::ValueArg<std::string> nameArg("r",
                                         "reference",
                                         "path to human reference",
                                         false,
                                         HUMAN_CHROMOSOME_10_STRIPPED_N_FILENAME,
                                         "string");
    cmd.add(nameArg);
    cmd.parse( argc, argv );
    HUMAN_CHROMOSOME_10_STRIPPED_N_FILENAME = nameArg.getValue();
    //get_statistics();
    //save_human_chromosome();
    //playground_dbg();
    //to_be_determined();
    //    code_to_violate_assertion();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
