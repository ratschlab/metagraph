#include <random>

#include "gtest/gtest.h"

#include "annotate.hpp"
#include "hashers.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";


std::vector<sdsl::uint256_t> generate_kmers(size_t num) {
    std::vector<sdsl::uint256_t> kmers(num);
    int mod = pow(num, 0.25);
    for (size_t i = 0; i < kmers.size(); ++i) {
        *(reinterpret_cast<uint64_t*>(&kmers[i]))     = rand() % mod;
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 1) = rand() % mod;
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 2) = rand() % mod;
        *(reinterpret_cast<uint64_t*>(&kmers[i]) + 3) = rand() % mod;
    }
    return kmers;
}

TEST(Annotate, RandomTestNoFalseNegative) {
    //create annotation
    annotate::BloomFilter bloom(7, 1000);
    annotate::ExactFilter exact;
    //generate a bunch of kmers
    auto kmers = generate_kmers(1000);
    size_t total = 0, fp = 0;
    for (size_t i = 0; i < kmers.size(); ++i) {
        if (i < kmers.size() / 2) {
            bloom.insert(&kmers[i], &kmers[i] + 1);
            exact.insert(&kmers[i], &kmers[i] + 1);
            ASSERT_TRUE(bloom.find(&kmers[i], &kmers[i] + 1));
            ASSERT_TRUE(exact.find(&kmers[i], &kmers[i] + 1));
        } else {
            if (exact.find(&kmers[i], &kmers[i] + 1)) {
                ASSERT_TRUE(bloom.find(&kmers[i], &kmers[i] + 1));
            }
            if (!exact.find(&kmers[i], &kmers[i] + 1)) {
                total++;
                if (bloom.find(&kmers[i], &kmers[i] + 1)) {
                    fp++;
                }
            }
        }
    }
    EXPECT_TRUE(total > fp) << "Total: " << total << " FP: " << fp << std::endl;
}

TEST(Annotate, RandomHashAnnotator) {
    annotate::HashAnnotation<annotate::BloomFilter> bloomhash(7);
    annotate::HashAnnotation<annotate::ExactFilter> exacthash;
    size_t num_bits = 5;
    std::vector<size_t> bounds(num_bits);
    std::iota(bounds.begin(), bounds.end(), 0);
    for (size_t i = 0; i < num_bits; ++i) {
        bloomhash.append_bit(1000);
        exacthash.append_bit();
    }
    ASSERT_EQ(bloomhash.size(), num_bits);
    ASSERT_EQ(exacthash.size(), num_bits);
    auto kmers = generate_kmers(1000);
    for (size_t i = 0; i < kmers.size(); ++i) {
        size_t pick_bits = 0;
        if (i < kmers.size()) {
            //insert into random bit positions
            pick_bits = rand() % (1u << num_bits);
            ASSERT_TRUE(pick_bits < (1u << num_bits));
            for (size_t j = 0; j < num_bits; ++j) {
                if (((1lu << j) | pick_bits) != pick_bits)
                    continue;

                //insert
                auto testbloom = bloomhash.insert(&kmers[i], &kmers[i] + 1, j);
                auto testexact = exacthash.insert(&kmers[i], &kmers[i] + 1, j);
                //check if it's there
                testbloom = bloomhash.find(&kmers[i], &kmers[i] + 1, j);
                testexact = exacthash.find(&kmers[i], &kmers[i] + 1, j);

                //test OR
                auto testbloom_merged = testbloom;
                annotate::merge_or(testbloom_merged, testexact);
                ASSERT_TRUE(annotate::equal(testbloom, testbloom_merged));

                //test bit
                ASSERT_TRUE(annotate::test_bit(testbloom, j));
                ASSERT_TRUE(annotate::test_bit(testexact, j));

                //test AND
                auto testbloom_and = testbloom;
                annotate::merge_and(testbloom_merged, testexact);
                ASSERT_TRUE(annotate::equal(testexact, testbloom_and));
            }
        }
        auto testbloom = bloomhash.find(&kmers[i], &kmers[i] + 1);
        auto testexact = exacthash.find(&kmers[i], &kmers[i] + 1);
        ASSERT_TRUE(annotate::equal(testbloom, annotate::merge_or(testbloom, testexact)));
    }
}

TEST(Annotate, HashIterator) {
    std::string test_string;
    for (size_t i = 0; i < 8; ++i) {
        test_string += std::string("$NATGC");
    }
    ASSERT_EQ(48llu, test_string.length());

    size_t num_hash_functions = 5;
    size_t kmer_size = 20;

    annotate::HashIterator hash_it(test_string, num_hash_functions, kmer_size);
    auto pos = hash_it.pos();
    ASSERT_EQ(num_hash_functions, hash_it.size());
    ASSERT_EQ(0llu, pos);

    auto hashes = annotate::hash_murmur(
            test_string,
            num_hash_functions,
            kmer_size
    );

    ASSERT_EQ(test_string.length() - kmer_size + 1, hashes.size());

    uint64_t bigint[2];
    for (size_t i = 0; i + kmer_size <= test_string.length(); ++i) {
        auto hash = annotate::hash_murmur(
                std::string(&test_string[i], kmer_size),
                num_hash_functions,
                kmer_size
        );
        ASSERT_EQ(1llu, hash.size());
        for (uint32_t j = 0; j < num_hash_functions; ++j) {
            ASSERT_NE('\0', *(&test_string[i] + kmer_size - 1));
            annotate::Murmur3Hasher(&test_string[i], kmer_size, j, &bigint[0]);
            ASSERT_EQ(bigint[0], (*hash_it)[j]);
            ASSERT_EQ(bigint[0], hashes[i][j]);
            ASSERT_EQ(bigint[0], hash[0][j]);
        }
        ++hash_it;
    }
    EXPECT_EQ(hashes.size(), hash_it.pos());
}

TEST(Annotate, ntHash) {

    //TODO: if N in string, ntHashIterator fails
    std::string test_string("NATGCA");

    size_t num_hash = 5;

    ntHashIterator hash_nt(test_string, num_hash, test_string.length());
    ASSERT_NE(hash_nt, hash_nt.end());

    auto hashes = annotate::hash(hash_nt, num_hash);
    ASSERT_EQ(1llu, hashes.size());

}

TEST(Annotate, HashIteratorInsert) {
    std::string test_string;
    for (size_t i = 0; i < 8; ++i) {
        test_string += std::string("$NATGC");
    }
    ASSERT_EQ(48llu, test_string.length());

    size_t num_hash_functions = 5;
    size_t kmer_size = 20;

    annotate::HashIterator hash_it(test_string, num_hash_functions, kmer_size);
    ntHashIterator hash_nt_it(test_string, num_hash_functions, kmer_size);
    ASSERT_EQ(num_hash_functions, hash_it.size());

    auto hashes = annotate::hash_murmur(
            test_string,
            num_hash_functions,
            kmer_size
    );

    ASSERT_EQ(test_string.length() - kmer_size + 1, hashes.size());

    annotate::HashAnnotation<annotate::BloomFilter> bloomhash(num_hash_functions);
    annotate::HashAnnotation<annotate::ExactFilter> exacthash;
    annotate::HashAnnotation<annotate::BloomFilter> bloomhash_it(num_hash_functions);
    annotate::HashAnnotation<annotate::ExactFilter> exacthash_it;
    annotate::HashAnnotation<annotate::ExactFilter> exacthash_it2;
    annotate::HashAnnotation<annotate::ExactFilter> exacthash_nt_it;

    bloomhash.append_bit(1000);
    exacthash.append_bit(1000);
    bloomhash_it.append_bit(1000);
    exacthash_it.append_bit(1000);
    exacthash_it2.append_bit(1000);
    exacthash_nt_it.append_bit(1000);

    while (hash_it != hash_it.end()) {
        exacthash_it2.insert(annotate::MultiHash(*hash_it, num_hash_functions), 0);
        ++hash_it;
    }

    while (hash_nt_it != hash_nt_it.end()) {
        exacthash_nt_it.insert(annotate::MultiHash(*hash_nt_it, num_hash_functions), 0);
        ++hash_nt_it;
    }

    for (size_t i = 0; i + kmer_size <= test_string.length(); ++i) {
        bloomhash.insert(&test_string[i], &test_string[i] + kmer_size, 0);
        exacthash.insert(&test_string[i], &test_string[i] + kmer_size, 0);
        bloomhash_it.insert(hashes[i], 0);
        exacthash_it.insert(hashes[i], 0);
    }
    EXPECT_TRUE(bloomhash == bloomhash_it);
    EXPECT_TRUE(exacthash == exacthash_it);
    EXPECT_TRUE(exacthash == exacthash_it2);

    for (size_t i = 0; i + kmer_size <= test_string.length(); ++i) {
        auto hash = annotate::hash_murmur(std::string(&test_string[i], kmer_size), num_hash_functions, kmer_size)[0];
        ASSERT_TRUE(exacthash.find(hash)[0]);
        ASSERT_TRUE(exacthash_it.find(hash)[0]);
        ASSERT_TRUE(exacthash_it2.find(hash)[0]);
        ASSERT_TRUE(bloomhash.find(hash)[0]);
        ASSERT_TRUE(bloomhash_it.find(hash)[0]);
        //ASSERT_TRUE(bloomhash_it2.find(hash)[0]);
    }
}


TEST(ColorCompressed, EmptyConstructor) {
    annotate::ColorCompressed annotation(5);
    EXPECT_EQ(0u, annotation.get(0).size());
    EXPECT_EQ(0u, annotation.get(1).size());
    EXPECT_EQ(0u, annotation.get(2).size());
    EXPECT_EQ(0u, annotation.get(3).size());
    EXPECT_EQ(0u, annotation.get(4).size());
}

TEST(ColorCompressed, add_label) {
    annotate::ColorCompressed annotation(5);
    annotation.add_label(0, "0");
    annotation.add_label(1, "0");
    annotation.add_label(2, "1");
    annotation.add_label(1, "2");

    EXPECT_EQ(std::set<std::string>({"0"}), annotation.get(0));
    EXPECT_EQ(std::set<std::string>({"0", "2"}), annotation.get(1));
    EXPECT_EQ(std::set<std::string>({"1"}), annotation.get(2));
    EXPECT_EQ(std::set<std::string>({}), annotation.get(3));
    EXPECT_EQ(std::set<std::string>({}), annotation.get(4));
}

TEST(ColorCompressed, set_label) {
    annotate::ColorCompressed annotation(5);
    annotation.set_label(0, { "Label0", "Label2", "Label8" });
    annotation.set_label(2, { "Label1", "Label2" });
    annotation.set_label(4, { "Label8" });

    EXPECT_EQ(std::set<std::string>({ "Label0", "Label2", "Label8" }), annotation.get(0));
    EXPECT_EQ(std::set<std::string>({}), annotation.get(1));
    EXPECT_EQ(std::set<std::string>({ "Label1", "Label2" }), annotation.get(2));
    EXPECT_EQ(std::set<std::string>({}), annotation.get(3));
    EXPECT_EQ(std::set<std::string>({ "Label8" }), annotation.get(4));
}

TEST(ColorCompressed, Serialization) {
    {
        annotate::ColorCompressed annotation(5);
        annotation.set_label(0, { "Label0", "Label2", "Label8" });
        annotation.set_label(2, { "Label1", "Label2" });
        annotation.set_label(4, { "Label8" });

        annotation.serialize(test_dump_basename + "_color_compressed");
    }
    {
        annotate::ColorCompressed annotation(5);
        ASSERT_FALSE(annotation.load(test_dump_basename + "_bad_file"));
        ASSERT_TRUE(annotation.load(test_dump_basename + "_color_compressed"));

        EXPECT_EQ(std::set<std::string>({ "Label0", "Label2", "Label8" }), annotation.get(0));
        EXPECT_EQ(std::set<std::string>({}), annotation.get(1));
        EXPECT_EQ(std::set<std::string>({ "Label1", "Label2" }), annotation.get(2));
        EXPECT_EQ(std::set<std::string>({}), annotation.get(3));
        EXPECT_EQ(std::set<std::string>({ "Label8" }), annotation.get(4));
    }
}

TEST(ColorCompressed, has_label) {
    std::unique_ptr<annotate::AnnotationCategory<std::set<std::string>>> annotation(
        new annotate::ColorCompressed(5)
    );
    annotation->set_label(0, {"Label0", "Label2", "Label8"});
    annotation->set_label(2, {"Label1", "Label2"});
    annotation->set_label(4, {"Label8"});

    EXPECT_FALSE(annotation->has_label(0, { "Label0", "Label1",
                                            "Label2", "Label4",
                                            "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_label(0, { "Label0",
                                            "Label2", "Label4",
                                            "Label8" }));

    EXPECT_TRUE(annotation->has_label(0, { "Label0", "Label2", "Label8" }));
    EXPECT_TRUE(annotation->has_label(0, { "Label0", "Label8" }));
    EXPECT_TRUE(annotation->has_label(0, { "Label2" }));
    EXPECT_TRUE(annotation->has_label(0, {}));

    EXPECT_FALSE(annotation->has_label(1, { "Label0", "Label1",
                                            "Label2", "Label4",
                                            "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_label(1, { "Label0",
                                            "Label2", "Label4",
                                            "Label8" }));

    EXPECT_FALSE(annotation->has_label(1, { "Label0", "Label2", "Label8" }));
    EXPECT_FALSE(annotation->has_label(1, { "Label0", "Label8" }));
    EXPECT_FALSE(annotation->has_label(1, { "Label2" }));
    EXPECT_TRUE(annotation->has_label(1, {}));

    EXPECT_FALSE(annotation->has_label(2, { "Label0", "Label1",
                                            "Label2", "Label4",
                                            "Label5", "Label8" }));

    EXPECT_FALSE(annotation->has_label(2, { "Label0",
                                            "Label2", "Label4",
                                            "Label8" }));

    EXPECT_FALSE(annotation->has_label(2, { "Label1", "Label2", "Label8" }));
    EXPECT_FALSE(annotation->has_label(2, { "Label1", "Label8" }));
    EXPECT_TRUE(annotation->has_label(2, { "Label1", "Label2" }));
    EXPECT_TRUE(annotation->has_label(2, { "Label2" }));
    EXPECT_TRUE(annotation->has_label(2, {}));
}

TEST(ColorCompressed, set_label_cache) {
    size_t graph_half_size = 500;
    annotate::ColorCompressed annotation(graph_half_size * 2);
    for (size_t i = 0; i < graph_half_size; ++i) {
        annotation.add_label(i, "Label1");
    }
    for (size_t i = graph_half_size; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, "Label2");
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get(i).size());
    }
}

TEST(ColorCompressed, set_label_random) {
    size_t graph_half_size = 500;
    annotate::ColorCompressed annotation(graph_half_size * 2);

    std::vector<std::string> colors { "Label1", "Label2" };

    for (size_t i = 0; i < 2 * graph_half_size; ++i) {
        annotation.add_label(i, colors[i % 2]);
    }
    for (size_t i = 0; i < 2 * graph_half_size; i+= 100) {
        ASSERT_EQ(1u, annotation.get(i).size());
    }
}
