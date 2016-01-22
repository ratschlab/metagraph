#include <iostream>
#include <assert.h>

#include <vector>
#include <string>
#include <config.hpp>

#include <dbg_succinct_libmaus.hpp>
#include <MockDatabase.cpp>
#include <RocksdbImpl.cpp>

/**
 * Tests that the mock database code is working properly.
 */
void test_MockDatabase() {
    dbg_database::IDatabase *db = new dbg_database::MockDatabase();

    std::string kmer = "TCGA";
    std::string annotation = "annotationA";

    db->annotate_kmer(kmer, annotation);
    // TODO fix.
    // assert(annotation == db->get_annotation(kmer));

    auto datastore = ((dbg_database::MockDatabase *) db)->getDatastore();
    assert(1 == datastore.size());
}

/*
 * Runs DBG graph construction on the mock database and makes sure
 * that the database is behaving properly.
 */
void test_DatabaseIntegration() {
    CFG config;
    config.k = 3;
    DBG_succ* graph = new DBG_succ(config.k, config);
    dbg_database::IDatabase *mockdb = new dbg_database::MockDatabase();

    seqan::String<seqan::Dna5> seq = "TCGA";
    seqan::CharString annotation = "annotationA";
    graph->add_annotation_for_seq(mockdb, seq, annotation);

    auto datastore = ((dbg_database::MockDatabase *) mockdb)->getDatastore();

    assert(2 == datastore.size());
    assert(datastore["TCG"] == "annotationA");
    assert(datastore["CGA"] == "annotationA");

    delete graph;
    delete mockdb;
}

/**
 * Runs the DBG graph construction on the RocksDB implementation
 */
void test_RocksDbIntegration() {
    CFG config;
    config.k = 3;
    DBG_succ* graph = new DBG_succ(config.k, config);
    dbg_database::IDatabase *therocksdb = new dbg_database::RocksdbImpl("/tmp/dbg-test-annotation-db");

    seqan::String<seqan::Dna5> seq = "TCGA";
    seqan::CharString annotation = "annotationA";
    graph->add_annotation_for_seq(therocksdb, seq, annotation);

    therocksdb->get_annotation("TCG");

    graph->add_annotation_for_seq(therocksdb, "AACTGATGGGTAATA", annotation);
    graph->add_annotation_for_seq(therocksdb, "AACTGATGGG", "annotationB");

    auto retrieved_annotations = therocksdb->get_annotation("GAT");
    assert("annotationA" == retrieved_annotations[0]);
    assert("annotationB" == retrieved_annotations[1]);

    delete graph;
    delete therocksdb;
}

int main() {
    test_MockDatabase();
    test_RocksDbIntegration();
}
