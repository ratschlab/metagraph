#include <iostream>
#include <assert.h>

#include <vector>
#include <string>
#include <config.hpp>

#include <dbg_succinct_libmaus.hpp>
#include <MockDatabase.cpp>

/**
 * Tests that the mock database code is working properly.
 */
void test_MockDatabase() {
    dbg_database::IDatabase *db = new dbg_database::MockDatabase();

    std::string kmer = "TCGA";
    std::string annotation = "annotationA";

    db->annotate_kmer(kmer, annotation);
    assert(annotation == db->get_annotation(kmer));

    auto datastore = ((dbg_database::MockDatabase *) db)->getDatastore();
    assert(1 == datastore.size());
}

/*
 * Runs DBG graph construction on the given database and makes sure
 * that the database is behaving properly.
 */
void test_DatabaseIntegration() {
    CFG config;
    config.k = 3;
    DBG_succ* graph = new DBG_succ(config.k, config);
    dbg_database::IDatabase *db = new dbg_database::MockDatabase();

    seqan::String<seqan::Dna5> seq = "TCGA";
    seqan::CharString annotation = "annotationA";
    graph->add_annotation_for_seq(db, seq, annotation);

    auto datastore = ((dbg_database::MockDatabase *) db)->getDatastore();

    assert(2 == datastore.size());
    assert(datastore["TCG"] == "annotationA");
    assert(datastore["CGA"] == "annotationA");
}

int main() {
    test_MockDatabase();
    test_DatabaseIntegration();
}
