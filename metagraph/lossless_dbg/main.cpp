#include <utility>

#include <iostream>
#define _DNA_GRAPH 1
#ifdef SELF_STANDING
#include <metagraph/dbg_succinct.hpp>
#include <metagraph/sequence_graph.hpp>
#include <metagraph/sequence_io.hpp>
#include <gtest/gtest.h>
#else
#include "dbg_succinct.hpp"
#include "sequence_graph.hpp"
#include "sequence_io.hpp"
#include <gtest/gtest.h>
#endif

using node_index = SequenceGraph::node_index;

using namespace std;
using namespace std::string_literals;

const int k_k_mer = 4;
const string HUMAN_REFERENCE_FILENAME = "/Users/janstudeny/Library/Mobile Documents/com~apple~CloudDocs/Active projects/Master Thesis/genomic-data/GCF_000001405.38_GRCh38.p12_genomic.fna";
const int CHROMOSOME_NUMBER = 10;
const int READ_LENGTH = 100;
const int READ_COVERAGE = 10;
const int test_seed = 3424;
// todo change to proper thing
#define x first
#define y second
#define all(x) begin(x),end(x)



node_index k_mer_to_node(DBGSuccinct& de_Bruijn_Graph, const string &k_mer) {
    node_index result;
    de_Bruijn_Graph.map_to_nodes(k_mer,[&](node_index node) {result = node;},[](){return true;});
    return result;
}

class CompressedReads {
    using bifurcation_choices_t = vector<char>;
    using k_mer_t = string;
    using compressed_read_t = pair<k_mer_t,bifurcation_choices_t>;

public:
    CompressedReads(vector<string> raw_reads) :  read_length(raw_reads[0].length()), graph(k_k_mer) {
        build_graph(raw_reads);
        for(auto & read : raw_reads) {
            reads.insert(align_read(read));
        }
    }
    int compressed_size_without_reference() {
        // returns size in bits
        int size = 0;
        for(auto & read : reads) {
            // *2 for two bit encoding
            size += read.x.size()*2;
            size += read.y.size()*2;
        }
        // to be really fair
        // size += sizeof(read_length);
        return size;
    }
    double bits_per_symbol(bool include_reference=false) {
        assert(!include_reference);
        return static_cast<double>(compressed_size_without_reference())/(reads.size()*read_length);
    }
private:
    compressed_read_t align_read(const string &read) {
        auto k_mer = read.substr(0,k_k_mer);
        auto bifurcation_choices = bifurcation_choices_t();
        auto node = k_mer_to_node(graph,k_mer);
        // for all other characters
        for(auto& character : read.substr(k_k_mer)) {
            vector<node_index> outnodes;
            graph.adjacent_incoming_nodes(node, &outnodes);
            if (outnodes.size() > 1) {
                bifurcation_choices.push_back(character);
            }
            assert(outnodes.size() != 0);
        }
        return {k_mer,bifurcation_choices};
    }
    void build_graph(vector<string> raw_reads) {
        for(auto & read : raw_reads) {
            assert(read.size() >= k_k_mer);
            graph.add_sequence(read);
        }
    }

    multiset<compressed_read_t> reads;
    int read_length;
    size_t k_k_mer = 21;
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

TEST(CompressingReads,GetChromosomeWorks) {
    auto chromosome = get_human_chromosome(CHROMOSOME_NUMBER);
    EXPECT_EQ(chromosome.length(), 133'797'422);
    EXPECT_EQ(chromosome.substr(0,10),"NNNNNNNNNN");
}

static void playground_dbg() {
    DBGSuccinct dbgSuccinct(k_k_mer);
    dbgSuccinct.add_sequence("ATAGAGAGAGAGAGAGAG");

    auto node = k_mer_to_node(dbgSuccinct,"ATAG");
    auto next_node = dbgSuccinct.traverse(node,'A');
    cout << dbgSuccinct.get_node_sequence(next_node) << endl;
    cout << "Hello metagraph!" << endl;
}
void to_be_determined() {
    auto chromosome = get_human_chromosome(CHROMOSOME_NUMBER);
    auto sampler = Sampler(chromosome,test_seed);
    auto reads = sampler.sample_coverage(READ_LENGTH, READ_COVERAGE);
    auto compressed_reads = CompressedReads(reads);
}


int main(int argc,char**argv) {
    //playground_dbg();
    to_be_determined();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}