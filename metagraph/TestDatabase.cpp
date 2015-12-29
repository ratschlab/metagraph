#include "MockDatabase.hpp"
#include <iostream>
#include <assert.h>

using namespace models;
using namespace std;
using namespace dbg_database;

void test_MockDatabase() {
    MockDatabase db;

    // test putting
    Kmer kmer;
    kmer.name = "TCGA";
    models::Annotation annotation;
    vector<models::Tag> tags;
    models::Tag tag1;
    models::Tag tag2;
    models::Tag tag3;
    tag1.name = "tag1";
    tag2.name = "tag2";
    tag3.name = "tag3";
    tags.push_back(tag1);
    tags.push_back(tag2);
    tags.push_back(tag3);
    annotation.tags = tags;
    assert(1 == db.annotate_kmer(kmer, annotation));

    // test getting
    Annotation gotten_annotation = db.get_annotation(kmer);
    auto sometags = gotten_annotation.tags;
    assert(3 == sometags.size());
    assert("tag1" == sometags[0].name);

    // test putting more than one thing.
    Kmer kmer2;
    kmer2.name = "AAAA";
    assert(2 == db.annotate_kmer(kmer2, annotation));

    // test adding new annotations
    Annotation new_annotation;
    vector<Tag> new_tags;
    new_tags.push_back((Tag) {"tag1"});
    new_tags.push_back((Tag) {"new_tag1"});
    new_tags.push_back((Tag) {"new_tag2"});
    new_tags.push_back((Tag) {"new_tag3"});
    new_annotation.tags = new_tags;
    db.annotate_kmer(kmer, new_annotation);
    assert(6 == db.get_annotation(kmer).tags.size());
}

int main() {
    test_MockDatabase();
}
