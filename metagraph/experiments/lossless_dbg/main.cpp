#include <utility>
#include <iostream>
#include <map>
#include <filesystem>
#include <nlohmann/json.hpp>
#include <ProgressBar.hpp>

using json = nlohmann::json;
#define _DNA_GRAPH 1

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#pragma clang diagnostic ignored "-Wcomma"

#include "dbg_succinct.hpp"
#include "sequence_graph.hpp"
#include "sequence_io.hpp"
#include "dbg_succinct_construct.hpp"
#include "dbg_hash.hpp"

#pragma clang diagnostic pop
#include <gtest/gtest.h>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"


using namespace std;
using namespace std::string_literals;
namespace fs = std::filesystem;

// debugging functions
struct d_t {
    template<typename T> d_ & operator,(const T & x) {
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
const string HUMAN_REFERENCE_FILENAME = local_file("genomic_data/GCF_000001405.38_GRCh38.p12_genomic.fna");
const string HUMAN_CHROMOSOME_10_SAMPLE = local_file("genomic_data/human_chromosome_10_sample.fasta");
const string HUMAN_CHROMOSOME_10_FILENAME = local_file("genomic_data/human_chromosome_10.fasta");
const string HUMAN_CHROMOSOME_10_STRIPPED_N_FILENAME = local_file("genomic_data/human_chromosome_10_n_trimmed.fasta");
const string JSON_OUTPUT_FILE = local_file("statistics.json");
const int CHROMOSOME_NUMBER = 10;
const int READ_LENGTH = 100;
const double READ_COVERAGE = 0.00001;
const int test_seed = 3424;
const int DEFAULT_K_KMER = 21;
// todo change to proper thing
#define x first
#define y second
#define all(x) begin(x),end(x)

//class ProgressBar {
//public:
//    int skips_before_report;
//    int current_time = 0;
//    ProgressBar(int skips_before_report) : skips_before_report(skips_before_report) {};
//    void invoke(function<void(void)> fn) {
//        current_time++;
//        if ((current_time%skips_before_report) == 0) {
//            fn();
//        }
//    }
//};


// openmp reductions

void reduce_maps(\
                 std::map<int, int>& output, \
                 std::map<int, int>& input)
{
    for (auto& X : input) {
        output[X.first] += X.second;
    }
}

#pragma omp declare reduction(map_reduction : \
std::map<int, int> : \
reduce_maps(omp_out, omp_in)) \
initializer(omp_priv(omp_orig))


node_index k_mer_to_node(DBGSuccinct& de_Bruijn_Graph,
                         const string &k_mer) {
    node_index result;
    de_Bruijn_Graph.map_to_nodes(k_mer, [&](node_index node) { result = node; },
        [](){return true;}
    );
    return result;
}

class CompressedReads {
    using bifurcation_choices_t = vector<char>;
    using kmer_t = string;
    using compressed_read_t = pair<kmer_t, bifurcation_choices_t>;

public:
    // TODO: read about lvalues, rvalues, etc
    // http://thbecker.net/articles/rvalue_references/section_01.html
    // const value
    // CompressedReads(const vector<string> &raw_reads)
    // non-const pointer to modify
    // CompressedReads(vector<string> *raw_reads)
    CompressedReads(const vector<string>& raw_reads, int k_kmer = defult_k_k_mer)
          : k_k_mer(k_kmer),
            read_length(raw_reads[0].length()),
            graph(dbg_succ_graph_constructor(raw_reads,k_kmer)) {
        for(auto & read : raw_reads) {
            compressed_reads.insert(align_read(read));
        }
    }

    int compressed_size_without_reference() {
        // returns size in bits
        int size = 0;
        for(auto &read : reads) {
            // *2 for two bit encoding
            size += read.x.size() * 2;
            size += read.y.size() * 2;
        }
        // to be really fair
        // size += sizeof(read_length);
        return size;
    }

    double bits_per_symbol(bool include_reference = false) {
        assert(!include_reference);
        return static_cast<double>(compressed_size_without_reference())
                    / (compressed_reads.size()*read_length);
    }
    vector<string> get_reads() {
        vector<string> reads;
        for(auto& compressed_read : compressed_reads) {
            reads.push_back(decompress_read(compressed_read));
        }
        return reads;
    }
    
private:
    compressed_read_t align_read(const string &read) {
        auto k_mer = read.substr(0,k_k_mer);
        auto bifurcation_choices = bifurcation_choices_t();
        auto node = graph.kmer_to_node(read);
        // for all other characters
        for (auto &character : read.substr(k_k_mer)) {
            vector<node_index> outnodes;
            graph.adjacent_incoming_nodes(node, &outnodes);
            if (outnodes.size() > 1) {
                bifurcation_choices.push_back(character);
            }
            assert(!outnodes.empty());
            
        }
        return {k_mer,bifurcation_choices};
    }
    static DBG_succ* dbg_succ_graph_constructor(vector<string> raw_reads, int k_k_mer) {
        auto graph_constructor = DBGSuccConstructor(k_k_mer);
        cerr << "Starting building the graph" << endl;
        for(auto &read : raw_reads) {
            assert(read.size() >= k_k_mer);
            graph_constructor.add_sequence(read);
        }
        return new DBG_succ(&graph_constructor);
        //graph = DBGSuccinct(stupid_old_representation);
    }
    
    string decompress_read(compressed_read_t compressed_read) {
        string read;
        auto& [starting_kmer,bifurcation_choices] = compressed_read;
        auto current_bifurcation_choice = bifurcation_choices.begin();
        auto node = k_mer_to_node(graph,starting_kmer);
        char next_char;
        while(read.size() < read_length) {
            int outgoing_degree = 0;
            graph.call_outgoing_kmers(node,[&](node_index next_node, char character) {
                outgoing_degree++;
                if (outgoing_degree>1) {
                    if (character == *current_bifurcation_choice) {
                        //initial guess was wrong
                        node = next_node;
                        next_char = character;
                    }
                }
                else {
                    //initial guess
                    node = next_node;
                    next_char = character;
                }
            });
            if (outgoing_degree>1) {
                current_bifurcation_choice++;//prepare next choice
            }
            read += next_char;
        }
        return read;
    }
    
    
    multiset<compressed_read_t> compressed_reads;
    int read_length;
    int k_k_mer;
    static const int defult_k_k_mer = 21;
    DBGSuccinct graph;
};

class SamplerConvenient {
public:
    virtual string sample(int length) = 0;
    virtual int reference_size() = 0;
    virtual vector<string> sample_coverage(int length, double coverage) {
        int count = ceil(reference_size()*coverage/length);
        return sample(length, count);
    }
    virtual vector<string> sample(int length, int count) {
        auto res = vector<string>();
        for(int i=0;i<count;i++) {
            res.push_back(sample(length));
        }
        return res;
    }
};

class Sampler : public SamplerConvenient {
public:
    Sampler(string reference, unsigned int seed) : reference(std::move(reference)) {
        generator = std::mt19937(seed); //Standard mersenne_twister_engine seeded with rd()
    };
    string sample(int length) override {
        std::uniform_int_distribution<> dis(0, reference.length()-1-length);
        return reference.substr(dis(generator),length);
    }
    int reference_size() override {
        return reference.size();
    }
private:
    string reference;
    std::mt19937 generator;
};

class DeterministicSampler : public SamplerConvenient {
public:
    DeterministicSampler(vector<string> samples, int reference_size) : _reference_size(reference_size), samples(std::move(samples)) {};
    string sample(int length) override {
        string sample = samples[current_sample];
        assert(length==sample.length());
        current_sample = (current_sample + 1) % samples.size();
        return sample;
    }
    int reference_size() override {
        return _reference_size;
    }
    vector<string> samples;
    int _reference_size;
    int current_sample = 0;
};

void transform_to_fasta(const string &filename,vector<string> reads) {
    ofstream myfile;
    myfile.open (filename);
    for(auto& read : reads) {
        myfile << ">" << endl;
        myfile << read << endl;
    }
    myfile.close();
}

vector<string> read_reads_from_fasta(const string &filename) {
    vector<string> result;
    read_fasta_file_critical(
                             filename,
                             [&](kseq_t* read) {
                                 result.push_back(read->seq.s);
                             });
    return result;
}

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

TEST(SamplerTest,SampleNoRandom) {
    auto sampler = Sampler("AAAAAAAAA",test_seed);
    ASSERT_EQ(sampler.sample(2),"AA");
}

TEST(SamplerTest,SampleNormal) {
    auto sampler = Sampler("ADFAGADFDS",test_seed);
    ASSERT_EQ(sampler.sample(4),"ADFD");
}
TEST(SamplerTest,SampleCoverage) {
    auto sequence = "ADFAGADFDS"s;
    auto sampler = Sampler(sequence,test_seed);
    auto reads = sampler.sample_coverage(sequence.length()/2, 1);
    ASSERT_EQ(reads.size(), 2);
}

TEST(CompressingReads,GetChromosomeWorks) {
    auto chromosome = get_human_chromosome(CHROMOSOME_NUMBER);
    EXPECT_EQ(chromosome.length(), 133'797'422);
    EXPECT_EQ(chromosome.substr(0,10),"NNNNNNNNNN");
}

TEST(CompressedReads,IdentityTest1) {
    set<string> reads = {"ATGCGATCGATATGCGAGA",
                         "ATGCGATCGAGACTACGAG",
                         "GTACGATAGACATGACGAG",
                         "ACTGACGAGACACAGATGC"};
    auto compressed_reads = CompressedReads(vector<string>(all(reads)));
    auto decompressed_reads = compressed_reads.get_reads();
    set<string> decompressed_read_set = set<string>(all(decompressed_reads));
    ASSERT_EQ(reads, decompressed_read_set);
}

static void playground_dbg() {
    //    DBGSuccinct dbgSuccinct(k_k_mer);
    //    dbgSuccinct.add_sequence("ATAGAGAGAGAGAGAGAG");
    //
    //    auto node = k_mer_to_node(dbgSuccinct,"ATAG");
    //    auto next_node = dbgSuccinct.traverse(node,'A');
    //    cout << dbgSuccinct.get_path_sequence({next_node}) << endl;
    //    cout << "Hello metagraph!" << endl;
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

    // openmp doesn't work with maps
    //#pragma omp parallel for reduction(map_reduction:kmer_outgoing_edges_statistics)
    //for(auto it = graph.indices_.begin(); it != graph.indices_.end(); it++)
    vector<node_index> outgoing_edges;

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
    get_statistics();
    //save_human_chromosome();
    //playground_dbg();
    //to_be_determined();
    //    code_to_violate_assertion();
    //    ::testing::InitGoogleTest(&argc, argv);
    //    return RUN_ALL_TESTS();
}
