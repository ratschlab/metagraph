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

TEST(ThreadPool, EmptyConstructor) {
    utils::ThreadPool pool(1);
    utils::ThreadPool pool2(2);
    utils::ThreadPool pool4(4);
    utils::ThreadPool pool20(20);
}

TEST(ThreadPool, EmptyTasks) {
    for (size_t i = 1; i < 20; ++i) {
        utils::ThreadPool pool(i);
        for (size_t t = 0; t < 1000; ++t) {
            pool.enqueue([]() {});
        }
    }
}

TEST(ThreadPool, EmptyJoin) {
    utils::ThreadPool pool(1);
    utils::ThreadPool pool2(2);
    utils::ThreadPool pool4(4);
    utils::ThreadPool pool20(20);
    pool.join();
    pool2.join();
    pool4.join();
    pool20.join();
    pool.join();
    pool2.join();
    pool4.join();
    pool20.join();
    pool.join();
    pool2.join();
    pool4.join();
    pool20.join();
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

TEST(ThreadPool, DummyEmptyJoin) {
    utils::ThreadPool pool(0);
    pool.join();
    pool.join();
    pool.join();
}

TEST(ThreadPool, DummyPoolRun) {
    utils::ThreadPool pool(0);
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

TEST(ThreadPool, DummyFuture) {
    utils::ThreadPool pool(0);

    std::vector<std::future<size_t>> result;
    for (size_t t = 0; t < 1000; ++t) {
        result.emplace_back(pool.enqueue([&](size_t i) { return i; }, 1));
    }

    ASSERT_EQ(1000u, result.size());
    for (auto &value : result) {
        ASSERT_EQ(1u, value.get());
    }
}


TEST(AsyncActivity, RunUniqueOnly) {
    utils::AsyncActivity async;

    ASSERT_EQ(1, async.run_unique([](int i) { return i; }, 1));
    ASSERT_EQ('a', async.run_unique([]() { return 'a'; }));
    ASSERT_EQ("abc", async.run_unique([](char a, char b, char c) {
                                         return std::string(1, a) + b + c;
                                      }, 'a', 'b', 'c'));
}

TEST(AsyncActivity, RunAsyncOnly) {
    utils::AsyncActivity async;

    ASSERT_EQ(1, async.run_async([](int i) { return i; }, 1));
    ASSERT_EQ('a', async.run_async([]() { return 'a'; }));
    ASSERT_EQ("abc", async.run_async([](char a, char b, char c) {
                                        return std::string(1, a) + b + c;
                                     }, 'a', 'b', 'c'));
}

TEST(AsyncActivity, RunBoth) {
    utils::AsyncActivity async;

    ASSERT_EQ(1, async.run_async([](int i) { return i; }, 1));
    ASSERT_EQ('a', async.run_async([]() { return 'a'; }));
    ASSERT_EQ("abc", async.run_async([](char a, char b, char c) {
                                        return std::string(1, a) + b + c;
                                     }, 'a', 'b', 'c'));

    ASSERT_EQ(1, async.run_unique([](int i) { return i; }, 1));
    ASSERT_EQ('a', async.run_unique([]() { return 'a'; }));
    ASSERT_EQ("abc", async.run_unique([](char a, char b, char c) {
                                         return std::string(1, a) + b + c;
                                      }, 'a', 'b', 'c'));
}

TEST(AsyncActivity, RunBothParallel) {
    utils::AsyncActivity async;

    std::string result;

    std::thread t1([&](std::string *result) {
                    async.run_async(
                        [](std::string *result_) {
                            std::this_thread::sleep_for(std::chrono::milliseconds(200));
                            result_->push_back('a');
                            return true;
                        }, result);
                    }, &result);
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    async.run_unique([](std::string *result_) { result_->push_back('b'); }, &result);
    ASSERT_EQ(2u, result.size()) << result;
    EXPECT_EQ('b', result.back()) << result;

    std::thread t2([&](std::string *result) {
                    async.run_async(
                        [](std::string *result_) {
                            std::this_thread::sleep_for(std::chrono::milliseconds(200));
                            result_->push_back('a');
                            return true;
                        }, result);
                    }, &result);
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    async.run_unique([](std::string *result_) { result_->push_back('c'); }, &result);
    ASSERT_EQ(4u, result.size()) << result;
    EXPECT_EQ('c', result.back()) << result;

    t1.join();
    t2.join();
}

TEST(Misc, SplitString) {
    EXPECT_EQ(std::vector<std::string>(),
              utils::split_string("", ""));

    EXPECT_EQ(std::vector<std::string>({ "123" }),
              utils::split_string("123", ""));

    EXPECT_EQ(std::vector<std::string>({ "23" }),
              utils::split_string("123", "1"));

    EXPECT_EQ(std::vector<std::string>({ "1", "3" }),
              utils::split_string("123", "2"));

    EXPECT_EQ(std::vector<std::string>({ "1", "34", "67" }),
              utils::split_string("12342672", "2"));

    EXPECT_EQ(std::vector<std::string>({ "1", "34", "72" }),
              utils::split_string("126342672", "26"));

    EXPECT_EQ(std::vector<std::string>({ "126342672" }),
              utils::split_string("126342672", ""));
}

TEST(Misc, JoinStrings) {
    EXPECT_EQ("", utils::join_strings({ "" }, ""));
    EXPECT_EQ("23", utils::join_strings({ "23" }, ":"));
    EXPECT_EQ("1:3", utils::join_strings(
                std::vector<std::string>({ "1", "3" }), ":"));
    EXPECT_EQ("1:34:67", utils::join_strings(
                std::vector<std::string>({ "1", "34", "67" }), ":"));
    EXPECT_EQ("1::34::72", utils::join_strings(
                std::vector<std::string>({ "1", "34", "72" }), "::"));
    EXPECT_EQ("13472", utils::join_strings(
                std::vector<std::string>({ "1", "34", "72" }), ""));

    EXPECT_EQ(std::vector<std::string>({ "123" }),
              utils::split_string("123", ""));

    EXPECT_EQ(std::vector<std::string>({ "23" }),
              utils::split_string("123", "1"));
}
