#include <stdio.h>

#include "gtest/gtest.h"

#include "utils.hpp"


TEST(get_filetype, VCF) {
    EXPECT_EQ("VCF", utils::get_filetype("file.VCF"));
    EXPECT_EQ("VCF", utils::get_filetype("file.vcf"));
    EXPECT_EQ("VCF", utils::get_filetype("file.VCF.gz"));
    EXPECT_EQ("VCF", utils::get_filetype("file.vcf.gz"));
}

TEST(get_filetype, FASTA) {
    EXPECT_EQ("FASTA", utils::get_filetype("file.FASTA"));
    EXPECT_EQ("FASTA", utils::get_filetype("file.fasta"));
    EXPECT_EQ("FASTA", utils::get_filetype("file.FASTA.gz"));
    EXPECT_EQ("FASTA", utils::get_filetype("file.fasta.gz"));
}

TEST(get_filetype, FASTQ) {
    EXPECT_EQ("FASTQ", utils::get_filetype("file.fq"));
    EXPECT_EQ("FASTQ", utils::get_filetype("file.FQ"));
    EXPECT_EQ("FASTQ", utils::get_filetype("file.fastq"));
    EXPECT_EQ("FASTQ", utils::get_filetype("file.fq.gz"));
    EXPECT_EQ("FASTQ", utils::get_filetype("file.FQ.gz"));
}

TEST(get_filetype, IncorrectFiletype) {
    EXPECT_EQ("", utils::get_filetype("fq"));
    EXPECT_EQ("", utils::get_filetype("fasta"));
    EXPECT_EQ("", utils::get_filetype("bz2"));
    EXPECT_EQ("", utils::get_filetype("file.gz"));
    EXPECT_EQ("", utils::get_filetype("fasta.gz"));
}

bool colexicographically_greater(const std::string &s1, const std::string &s2) {
    return utils::colexicographically_greater(s1, s2);
}


TEST(colexicographically_greater, basics) {
    EXPECT_FALSE(colexicographically_greater("AAAA", "AAAA"));

    EXPECT_TRUE(colexicographically_greater("BAAA", "AAAA"));
    EXPECT_TRUE(colexicographically_greater("AAAB", "AAAA"));

    EXPECT_FALSE(colexicographically_greater("AAAA", "BAAA"));
    EXPECT_FALSE(colexicographically_greater("AAAA", "AAAB"));
}

TEST(colexicographically_greater, different_lengths) {
    EXPECT_TRUE(colexicographically_greater("AAAAA", "AAAA"));
    EXPECT_FALSE(colexicographically_greater("AAAAA", "BBBB"));

    EXPECT_TRUE(colexicographically_greater("AAAAB", "AAAA"));
    EXPECT_TRUE(colexicographically_greater("ABAAA", "AAAA"));
    EXPECT_TRUE(colexicographically_greater("BAAAA", "AAAA"));
}


TEST(seq_equal, basics) {
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("ABAAACD"), 0));
    EXPECT_FALSE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 0));
    EXPECT_FALSE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 1));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 2));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 3));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 4));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 5));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 6));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 7));
    EXPECT_TRUE(utils::seq_equal(std::string("ABAAACD"), std::string("BAAAACD"), 100));
}

TEST(ThreadPool, SingleThreadEmpty) {
    utils::ThreadPool pool(1);
    pool.join();
}

TEST(ThreadPool, SingleThread) {
    utils::ThreadPool pool(1);
    std::vector<size_t> result;
    std::mutex mu;

    pool.enqueue([&](size_t i) {
        mu.lock();
        result.push_back(i);
        mu.unlock();
    }, 1);
    pool.enqueue([&](size_t i) {
        mu.lock();
        result.push_back(i);
        mu.unlock();
    }, 1);
    pool.enqueue([&](size_t i) {
        mu.lock();
        result.push_back(i);
        mu.unlock();
    }, 1);
    pool.enqueue([&](size_t i) {
        mu.lock();
        result.push_back(i);
        mu.unlock();
    }, 1);
    pool.enqueue([&](size_t i) {
        mu.lock();
        result.push_back(i);
        mu.unlock();
    }, 1);

    pool.join();

    ASSERT_EQ(5u, result.size());
    for (int value : result) {
        ASSERT_EQ(1, value);
    }
}

TEST(ThreadPool, MultiThreadTwo) {
    utils::ThreadPool pool(2);
    size_t num_tasks = 10000;

    std::vector<size_t> result;
    std::mutex mu;
    for (size_t t = 0; t < num_tasks; ++t) {
        pool.enqueue([&](size_t i) {
            size_t count = 3;
            for (size_t k = 0; k < num_tasks + i; ++k) {
                count += count * 3 - count / 2 + 1;
            }
            mu.lock();
            result.push_back(count);
            mu.unlock();
        }, t);
    }

    pool.join();

    ASSERT_EQ(num_tasks, result.size());
}

TEST(ThreadPool, MultiThreadFour) {
    utils::ThreadPool pool(4);
    size_t num_tasks = 10000;

    std::vector<size_t> result;
    std::mutex mu;
    for (size_t t = 0; t < num_tasks; ++t) {
        pool.enqueue([&](size_t i) {
            size_t count = 3;
            for (size_t k = 0; k < num_tasks + i; ++k) {
                count += count * 3 - count / 2 + 1;
            }
            mu.lock();
            result.push_back(count);
            mu.unlock();
        }, t);
    }

    pool.join();

    ASSERT_EQ(num_tasks, result.size());
}

TEST(ThreadPool, MultiThread) {
    for (size_t i = 2; i < 20; ++i) {
        utils::ThreadPool pool(i);
        std::vector<size_t> result;
        std::mutex mu;
        for (size_t t = 0; t < 1000; ++t) {
            pool.enqueue([&](size_t i) {
                mu.lock();
                result.push_back(i);
                mu.unlock();
            }, 1);
        }

        pool.join();

        ASSERT_EQ(1000u, result.size());
        for (int value : result) {
            ASSERT_EQ(1, value);
        }
    }
}

TEST(ThreadPool, MultiThreadFuture) {
    for (size_t i = 2; i < 20; ++i) {
        utils::ThreadPool pool(i);

        std::vector<std::future<size_t>> result;
        for (size_t t = 0; t < 1000; ++t) {
            result.emplace_back(pool.enqueue([&](size_t i) { return i; }, 1));
        }

        ASSERT_EQ(1000u, result.size());
        for (auto &value : result) {
            ASSERT_EQ(1u, value.get());
        }
    }
}

TEST(ThreadPool, MultiThreadException) {
    for (size_t i = 2; i < 20; ++i) {
        try {
            utils::ThreadPool pool(i);

            std::vector<std::future<size_t>> result;
            for (size_t t = 0; t < 1000; ++t) {
                result.emplace_back(pool.enqueue([&](size_t i) {
                    throw std::runtime_error("error");
                    return i;
                }, 1));
            }

            ASSERT_EQ(1000u, result.size());
            for (auto &value : result) {
                ASSERT_EQ(1u, value.get());
            }
        } catch (...) {
            continue;
        }
    }
}
