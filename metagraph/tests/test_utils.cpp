#include "gtest/gtest.h"
#include "test_helpers.hpp"

#include "common/vector.hpp"
#include "common/utils/string_utils.hpp"
#include "common/utils/file_utils.hpp"
#include "common/algorithms.hpp"
#include "common/vectors/bitmap_mergers.hpp"
#include "common/threads/threading.hpp"
#include "common/vectors/bit_vector_adaptive.hpp"
#include "seq_io/formats.hpp"


namespace {

using namespace mtg;
using namespace mtg::seq_io;

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";

const std::vector<sdsl::bit_vector> vectors {
    { 0, 1, 1, 0, 0, 1, 0 },
    { 0, 1, 1, 1, 1, 1, 0 }
};

const std::vector<sdsl::bit_vector> matrix {
    { 0, 0 },
    { 1, 1 },
    { 1, 1 },
    { 0, 1 },
    { 0, 1 },
    { 1, 1 },
    { 0, 0 }
};

const std::vector<std::vector<uint64_t>> indices {
    { 1, 0 }, { 1, 1 },
    { 2, 0 }, { 2, 1 },
    { 3, 1 },
    { 4, 1 },
    { 5, 0 }, { 5, 1 }
};


std::vector<bit_vector_small> generate_rows() {
    std::vector<bit_vector_small> rows;
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(rows),
        [](const auto &a) {
            return bit_vector_small(a);
        });

    return rows;
}

void check_rows(utils::RowsFromColumnsTransformer&& rct) {
    ASSERT_EQ(2u, rct.columns());
    ASSERT_EQ(7u, rct.rows());
    ASSERT_EQ(8u, rct.values_left());
    // ASSERT_EQ(std::vector<uint64_t>({ 3, 5 }), rct.num_set_bits());

    uint64_t i = 0;
    rct.call_rows<Vector<uint64_t>>([&](const auto &row_indices) {
        sdsl::bit_vector bv(rct.columns());
        for (auto j : row_indices) {
            bv[j] = 1;
        }
        EXPECT_EQ(matrix[i++], bv) << i;
    });

    ASSERT_EQ(7u, i);
}

void check_indices(utils::RowsFromColumnsTransformer&& rct) {
    ASSERT_EQ(2u, rct.columns());
    ASSERT_EQ(7u, rct.rows());
    ASSERT_EQ(8u, rct.values_left());
    // ASSERT_EQ(std::vector<uint64_t>({ 3, 5 }), rct.num_set_bits());

    uint64_t counter = 0;
    for (uint64_t i = 0; i < 8; ++i) {
        rct.call_next([&](auto row, auto column) {
            counter++;
            EXPECT_EQ(std::vector<uint64_t>({ row, column }),
                      indices[i]) << i;
        });
    }

    EXPECT_EQ(8u, counter);
    EXPECT_EQ(0u, rct.values_left());
}

TEST(Utils, RowsFromColumnsTransformerCallRowsColumns) {
    auto rows = generate_rows();
    utils::RowsFromColumnsTransformer rct(rows);
    check_rows(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallIndicesColumns) {
    auto rows = generate_rows();
    utils::RowsFromColumnsTransformer rct(rows);
    check_indices(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallRowsConcat) {
    bit_vector_small vectors_small{ 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0 };
    auto rct = utils::RowsFromColumnsTransformer(vectors_small, 7);
    check_rows(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallIndicesConcat) {
    bit_vector_small vectors_small{ 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0 };
    auto rct = utils::RowsFromColumnsTransformer(vectors_small, 7);
    check_indices(std::move(rct));
}

TEST(Utils, TempFileInitialize) {
    utils::TempFile tmp;
}

TEST(Utils, TempFileOpenWrite) {
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ofstream().good());
    }
}

TEST(Utils, TempFileOpenWriteRead) {
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ofstream().good());
        ASSERT_TRUE(tmp.ifstream().good());
    }
}

TEST(Utils, TempFileCheckStateFlow) {
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ifstream().good());
        EXPECT_DEATH(tmp.ofstream().good(), "Can't write after reading");
    }
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ofstream().good());
        ASSERT_TRUE(tmp.ofstream().good());
        ASSERT_TRUE(tmp.ifstream().good());
        ASSERT_TRUE(tmp.ifstream().good());
        EXPECT_DEATH(tmp.ofstream().good(), "Can't write after reading");
    }
}

TEST(Utils, TempFileReadWritePairs) {
    {
        utils::TempFile tmp[2];

        ASSERT_TRUE(tmp[0].ofstream().good());
        ASSERT_TRUE(tmp[0].ifstream().good());

        ASSERT_TRUE(tmp[1].ofstream().good());
        ASSERT_TRUE(tmp[1].ifstream().good());

        EXPECT_TRUE(tmp[0].ifstream().peek() == std::ifstream::traits_type::eof());
        EXPECT_TRUE(tmp[1].ifstream().peek() == std::ifstream::traits_type::eof());
    }

    utils::TempFile tmp[2];

    tmp[0].ofstream() << "test0 string0" << std::endl;
    tmp[1].ofstream() << "test1 string1" << std::endl;
    tmp[0].ofstream() << "test0 string0" << std::endl;
    tmp[1].ofstream() << "test1 string1" << std::endl;

    std::ifstream& in0 = tmp[0].ifstream();
    std::ifstream& in1 = tmp[1].ifstream();

    std::string temp;
    in0 >> temp; EXPECT_EQ("test0", temp);
    in1 >> temp; EXPECT_EQ("test1", temp);
    in0 >> temp; EXPECT_EQ("string0", temp);
    in1 >> temp; EXPECT_EQ("string1", temp);
    in0 >> temp; EXPECT_EQ("test0", temp);
    in1 >> temp; EXPECT_EQ("test1", temp);
    in0 >> temp; EXPECT_EQ("string0", temp);
    in1 >> temp; EXPECT_EQ("string1", temp);
    EXPECT_FALSE(in0 >> temp);
    EXPECT_FALSE(in1 >> temp);
    EXPECT_TRUE(in0.peek() == std::ifstream::traits_type::eof());
    EXPECT_TRUE(in1.peek() == std::ifstream::traits_type::eof());
}

template <typename T>
std::vector<T> uniqueize(const std::vector<T> &input) {
    auto vector = input;
    vector.erase(std::unique(vector.begin(), vector.end()), vector.end());
    return vector;
}

using utils::sample_indexes;

TEST(Utils, sample_indexes) {
    std::mt19937 gen;
    gen.seed(14);

    ASSERT_EQ(0u, sample_indexes(10, 0, gen).size());
    ASSERT_EQ(5u, sample_indexes(10, 5, gen).size());
    ASSERT_EQ(6u, sample_indexes(10, 6, gen).size());
    ASSERT_EQ(9u, sample_indexes(10, 9, gen).size());
    ASSERT_EQ(10u, sample_indexes(10, 10, gen).size());
    ASSERT_EQ(10u, sample_indexes(10, 11, gen).size());
    ASSERT_EQ(10u, sample_indexes(10, 1000, gen).size());

    ASSERT_EQ(0u, sample_indexes(0, 0, gen).size());
    ASSERT_EQ(0u, sample_indexes(0, 1000, gen).size());

    auto generated = sample_indexes(10, 0, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 5, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 6, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 9, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 10, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 11, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 1000, gen);
    ASSERT_EQ(generated, uniqueize(generated));

    generated = sample_indexes(100'000, 1'000, gen);
    EXPECT_TRUE(*std::min_element(generated.begin(), generated.end()) < 100'000u / 5);
    EXPECT_TRUE(*std::max_element(generated.begin(), generated.end()) > 4 * 100'000u / 5)
        << *std::max_element(generated.begin(), generated.end());

    ASSERT_EQ(2'000'000u, sample_indexes(100'000'000, 2'000'000, gen).size());
    ASSERT_EQ(12'000'000u, sample_indexes(100'000'000, 12'000'000, gen).size());
}


TEST(get_filetype, VCF) {
    EXPECT_EQ("VCF", file_format("file.VCF"));
    EXPECT_EQ("VCF", file_format("file.vcf"));
    EXPECT_EQ("VCF", file_format("file.VCF.gz"));
    EXPECT_EQ("VCF", file_format("file.vcf.gz"));
}

TEST(get_filetype, FASTA) {
    EXPECT_EQ("FASTA", file_format("file.FASTA"));
    EXPECT_EQ("FASTA", file_format("file.fasta"));
    EXPECT_EQ("FASTA", file_format("file.FASTA.gz"));
    EXPECT_EQ("FASTA", file_format("file.fasta.gz"));
}

TEST(get_filetype, FASTQ) {
    EXPECT_EQ("FASTQ", file_format("file.fq"));
    EXPECT_EQ("FASTQ", file_format("file.FQ"));
    EXPECT_EQ("FASTQ", file_format("file.fastq"));
    EXPECT_EQ("FASTQ", file_format("file.fq.gz"));
    EXPECT_EQ("FASTQ", file_format("file.FQ.gz"));
}

TEST(get_filetype, IncorrectFiletype) {
    EXPECT_EQ("", file_format("fq"));
    EXPECT_EQ("", file_format("fasta"));
    EXPECT_EQ("", file_format("bz2"));
    EXPECT_EQ("", file_format("file.gz"));
    EXPECT_EQ("", file_format("fasta.gz"));
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
    ThreadPool pool(1);
    ThreadPool pool2(2);
    ThreadPool pool4(4);
    ThreadPool pool20(20);
}

TEST(ThreadPool, EmptyTasks) {
    for (size_t i = 1; i < 20; ++i) {
        ThreadPool pool(i);
        for (size_t t = 0; t < 1000; ++t) {
            pool.enqueue([]() {});
        }
    }
}

TEST(ThreadPool, EmptyJoin) {
    ThreadPool pool(1);
    ThreadPool pool2(2);
    ThreadPool pool4(4);
    ThreadPool pool20(20);
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
    ThreadPool pool(1);
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
    ThreadPool pool(2);
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
    ThreadPool pool(4);
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
        ThreadPool pool(i);
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
        ThreadPool pool(i);

        std::vector<std::shared_future<size_t>> result;
        for (size_t t = 0; t < 1000; ++t) {
            result.emplace_back(pool.enqueue([&](size_t i) { return i; }, 1));
        }

        ASSERT_EQ(1000u, result.size());
        for (auto &value : result) {
            ASSERT_EQ(1u, value.get());
        }
    }
}

void throw_from_worker() {
    for (size_t i = 2; i < 20; ++i) {
        try {
            ThreadPool pool(i);

            std::vector<std::shared_future<size_t>> result;
            for (size_t t = 0; t < 1000; ++t) {
                result.emplace_back(pool.enqueue([&](size_t i) {
                    // This exception will be thrown in a worker thread, which
                    // should kill the main thread and the whole process as well
                    throw std::runtime_error("test exception");
                    return i;
                }, 1));
            }

        } catch (std::exception &e) {
            ASSERT_TRUE(false); // all exceptions must be thrown in workers
        }
    }
    std::this_thread::sleep_for(std::chrono::seconds(1));
}

TEST(ThreadPool, MultiThreadException) {
    // check that a worker throws even if we don't
    // wait for the result in the returned 'future'
    ASSERT_DEATH(throw_from_worker(), "");
}

TEST(ThreadPool, DummyEmptyJoin) {
    ThreadPool pool(0);
    pool.join();
    pool.join();
    pool.join();
}

TEST(ThreadPool, DummyPoolRun) {
    ThreadPool pool(0);
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
    ThreadPool pool(0);

    std::vector<std::shared_future<size_t>> result;
    for (size_t t = 0; t < 1000; ++t) {
        result.emplace_back(pool.enqueue([&](size_t i) { return i; }, 1));
    }

    ASSERT_EQ(1000u, result.size());
    for (auto &value : result) {
        ASSERT_EQ(1u, value.get());
    }
}


TEST(AsyncActivity, RunUniqueOnly) {
    AsyncActivity async;

    ASSERT_EQ(1, async.run_unique([](int i) { return i; }, 1));
    ASSERT_EQ('a', async.run_unique([]() { return 'a'; }));
    ASSERT_EQ("abc", async.run_unique([](char a, char b, char c) {
                                         return std::string(1, a) + b + c;
                                      }, 'a', 'b', 'c'));
}

TEST(AsyncActivity, RunAsyncOnly) {
    AsyncActivity async;

    ASSERT_EQ(1, async.run_async([](int i) { return i; }, 1));
    ASSERT_EQ('a', async.run_async([]() { return 'a'; }));
    ASSERT_EQ("abc", async.run_async([](char a, char b, char c) {
                                        return std::string(1, a) + b + c;
                                     }, 'a', 'b', 'c'));
}

TEST(AsyncActivity, RunBoth) {
    AsyncActivity async;

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
    AsyncActivity async;

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

    EXPECT_EQ(".1", utils::join_strings({ "", "1" }, "."));
    EXPECT_EQ("1.", utils::join_strings({ "1", "" }, "."));
    EXPECT_EQ(".1", utils::join_strings({ "", "1" }, ".", false));
    EXPECT_EQ("1.", utils::join_strings({ "1", "" }, ".", false));
    EXPECT_EQ("1", utils::join_strings({ "", "1" }, ".", true));
    EXPECT_EQ("1", utils::join_strings({ "1", "" }, ".", true));

    EXPECT_EQ("..", utils::join_strings({ "", "", "" }, "."));
    EXPECT_EQ("..3", utils::join_strings({ "", "", "3" }, "."));
    EXPECT_EQ("1..", utils::join_strings({ "1", "", "" }, "."));
    EXPECT_EQ(".2.", utils::join_strings({ "", "2", "" }, "."));
    EXPECT_EQ("1.2.", utils::join_strings({ "1", "2", "" }, "."));
    EXPECT_EQ("1..3", utils::join_strings({ "1", "", "3" }, "."));
    EXPECT_EQ(".2.3", utils::join_strings({ "", "2", "3" }, "."));
    EXPECT_EQ("1.2.3", utils::join_strings({ "1", "2", "3" }, "."));

    EXPECT_EQ("..", utils::join_strings({ "", "", "" }, ".", false));
    EXPECT_EQ("..3", utils::join_strings({ "", "", "3" }, ".", false));
    EXPECT_EQ("1..", utils::join_strings({ "1", "", "" }, ".", false));
    EXPECT_EQ(".2.", utils::join_strings({ "", "2", "" }, ".", false));
    EXPECT_EQ("1.2.", utils::join_strings({ "1", "2", "" }, ".", false));
    EXPECT_EQ("1..3", utils::join_strings({ "1", "", "3" }, ".", false));
    EXPECT_EQ(".2.3", utils::join_strings({ "", "2", "3" }, ".", false));
    EXPECT_EQ("1.2.3", utils::join_strings({ "1", "2", "3" }, ".", false));

    EXPECT_EQ("", utils::join_strings({ "", "", "" }, ".", true));
    EXPECT_EQ("3", utils::join_strings({ "", "", "3" }, ".", true));
    EXPECT_EQ("1", utils::join_strings({ "1", "", "" }, ".", true));
    EXPECT_EQ("2", utils::join_strings({ "", "2", "" }, ".", true));
    EXPECT_EQ("1.2", utils::join_strings({ "1", "2", "" }, ".", true));
    EXPECT_EQ("1.3", utils::join_strings({ "1", "", "3" }, ".", true));
    EXPECT_EQ("2.3", utils::join_strings({ "", "2", "3" }, ".", true));
    EXPECT_EQ("1.2.3", utils::join_strings({ "1", "2", "3" }, ".", true));
}

TEST(Misc, RemoveSuffix) {
    EXPECT_EQ("", utils::remove_suffix("", ""));
    EXPECT_EQ("", utils::remove_suffix("", "A"));
    EXPECT_EQ("", utils::remove_suffix("", "1AC"));
    EXPECT_EQ("", utils::remove_suffix("1", "1"));
    EXPECT_EQ("1.", utils::remove_suffix("1.1AC", "1AC"));
    EXPECT_EQ("1", utils::remove_suffix("1.1AC", ".1AC"));
    EXPECT_EQ("acb", utils::remove_suffix("acb.cde", ".cde"));
    EXPECT_EQ("acb.cde", utils::remove_suffix("acb.cde.efg", ".efg"));
    EXPECT_EQ("acb", utils::remove_suffix("acb.cde.efg", ".cde.efg"));
    EXPECT_EQ("acb.cde.efg", utils::remove_suffix("acb.cde.efg", ""));
    EXPECT_EQ("acb.cde.efg", utils::remove_suffix("acb.cde.efg", ".cde"));
    EXPECT_EQ("acb.cde.efg", utils::remove_suffix("acb.cde.efg", "Aacb.cde.efg"));
}

TEST(Misc, CountIntersectionEmpty) {
    {
        std::vector<int> first;
        std::vector<int> second;
        EXPECT_EQ(0u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 1, 2, 3, 4, 5 };
        std::vector<int> second;
        EXPECT_EQ(0u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first;
        std::vector<int> second { 1, 2, 3, 4, 5 };
        EXPECT_EQ(0u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 1, 2, 3, 4, 5 };
        std::vector<int> second { -5, -4, -3, -2, -1 };
        EXPECT_EQ(0u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { -1000 };
        std::vector<int> second { 1000 };
        EXPECT_EQ(0u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
}

TEST(Misc, CountIntersection) {
    {
        std::vector<int> first { 1, 2, 3, 4, 5 };
        std::vector<int> second { 1, 2, 3, 4, 5 };
        EXPECT_EQ(5u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 1, 2, 3, 4, 5 };
        std::vector<int> second { 0, 1, 2, 3, 4, 5, 6 };
        EXPECT_EQ(5u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 1, 2, 3, 4, 5, 6 };
        std::vector<int> second { 0, 1, 2, 3, 4, 5 };
        EXPECT_EQ(5u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 0, 1, 2, 3, 4, 5 };
        std::vector<int> second { -5, -4, -3, -2, -1, 0 };
        EXPECT_EQ(1u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
    {
        std::vector<int> first { 0, 1000000 };
        std::vector<int> second { -1000000, 0, 1000000 };
        EXPECT_EQ(2u, utils::count_intersection(first.begin(), first.end(),
                                                second.begin(), second.end()));
    }
}

TEST(Vector, ReserveInfinityCheckThrow) {
    Vector<int> vector;
    EXPECT_THROW(vector.reserve(1llu << 59), std::bad_alloc);
}

TEST(Vector, ResizeInfinityCheckThrow) {
    Vector<int> vector;
    EXPECT_THROW(vector.resize(1llu << 59), std::bad_alloc);
}

TEST(Deque, ResizeInfinityCheckThrow) {
    std::deque<int> array;
    EXPECT_THROW(array.resize(1llu << 60), std::bad_alloc);
}

TEST(Misc, get_quantile) {
    EXPECT_EQ(1, utils::get_quantile<int>({ {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1} }, 0.0));
    EXPECT_EQ(5, utils::get_quantile<int>({ {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1} }, 1.0));
    EXPECT_EQ(3, utils::get_quantile<int>({ {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1} }, 0.5));
    EXPECT_EQ(2, utils::get_quantile<int>({ {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1} }, 0.25));
    EXPECT_EQ(4, utils::get_quantile<int>({ {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 1} }, 0.75));

    EXPECT_EQ(1, utils::get_quantile<int>({ {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2} }, 0.0));
    EXPECT_EQ(5, utils::get_quantile<int>({ {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2} }, 1.0));
    EXPECT_EQ(3, utils::get_quantile<int>({ {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2} }, 0.5));
    EXPECT_EQ(2, utils::get_quantile<int>({ {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2} }, 0.25));
    EXPECT_EQ(4, utils::get_quantile<int>({ {1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2} }, 0.75));

    EXPECT_EQ(1, utils::get_quantile<int>({ {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 5} }, 0.0));
    EXPECT_EQ(5, utils::get_quantile<int>({ {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 5} }, 1.0));
    EXPECT_EQ(3, utils::get_quantile<int>({ {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 5} }, 0.5));
    EXPECT_EQ(2, utils::get_quantile<int>({ {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 5} }, 0.25));
    EXPECT_EQ(4, utils::get_quantile<int>({ {1, 2}, {2, 3}, {3, 5}, {4, 5}, {5, 5} }, 0.75));
}


template <typename T>
std::vector<T> insert_reference_impl(std::vector<T> vector,
                                     const std::vector<uint64_t> &new_pos,
                                     T value) {
    assert(std::is_sorted(new_pos.begin(), new_pos.end()));
    for (auto i : new_pos) {
        vector.insert(vector.begin() + i, value);
    }
    return vector;
}

bitmap_vector to_bitmap(size_t size, const std::vector<uint64_t> &set_bits_pos) {
    assert(std::is_sorted(set_bits_pos.begin(), set_bits_pos.end()));
    sdsl::bit_vector mask(size, false);
    for (auto i : set_bits_pos) {
        mask[i] = true;
    }
    return bitmap_vector(std::move(mask));
}

TEST(Misc, insert_zeros_to_empty) {
    for (const auto &new_pos : { std::vector<uint64_t>({}),
                                 std::vector<uint64_t>({ 0, }),
                                 std::vector<uint64_t>({ 0, 1, }),
                                 std::vector<uint64_t>({ 0, 1, 2, }),
                                 std::vector<uint64_t>({ 0, 1, 2, 3, }),
                                 std::vector<uint64_t>({ 0, 1, 2, 3, 4, }),
                                 std::vector<uint64_t>({ 0, 1, 2, 3, 4, 5, }) }) {
        std::vector<int> first;
        std::vector<int> second = first;
        auto expected = insert_reference_impl(first, new_pos, 0);
        utils::insert(&first, new_pos, 0);
        utils::insert(&second, to_bitmap(second.size() + new_pos.size(), new_pos), 0);

        EXPECT_EQ(expected, first);
        EXPECT_EQ(expected, second);
    }
}

TEST(Misc, insert_zeros_to_all_zeros) {
    for (const auto &new_pos : { std::vector<uint64_t>({}),
                                 std::vector<uint64_t>({ 0, }),
                                 std::vector<uint64_t>({ 0, 1, }),
                                 std::vector<uint64_t>({ 0, 1, 2, }),
                                 std::vector<uint64_t>({ 50, }),
                                 std::vector<uint64_t>({ 50, 51, }),
                                 std::vector<uint64_t>({ 50, 51, 52, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, 51, 52, 53, 54, 55, }) }) {
        std::vector<int> first(50, 0);
        std::vector<int> second = first;
        auto expected = insert_reference_impl(first, new_pos, 0);
        utils::insert(&first, new_pos, 0);
        utils::insert(&second, to_bitmap(second.size() + new_pos.size(), new_pos), 0);

        EXPECT_EQ(expected, first);
        EXPECT_EQ(expected, second);
    }
}

TEST(Misc, insert_zeros_to_all_ones) {
    for (const auto &new_pos : { std::vector<uint64_t>({}),
                                 std::vector<uint64_t>({ 0, }),
                                 std::vector<uint64_t>({ 0, 1, }),
                                 std::vector<uint64_t>({ 0, 1, 2, }),
                                 std::vector<uint64_t>({ 50, }),
                                 std::vector<uint64_t>({ 50, 51, }),
                                 std::vector<uint64_t>({ 50, 51, 52, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, }),
                                 std::vector<uint64_t>({ 0, 10, 20, 30, 40, 50, 51, 52, 53, 54, 55, }) }) {
        std::vector<int> first(50, 1);
        std::vector<int> second = first;
        auto expected = insert_reference_impl(first, new_pos, 0);
        utils::insert(&first, new_pos, 0);
        utils::insert(&second, to_bitmap(second.size() + new_pos.size(), new_pos), 0);

        EXPECT_EQ(expected, first);
        EXPECT_EQ(expected, second);
    }
}

TEST(Misc, drag_and_mark_segments) {
    EXPECT_EQ(std::vector<bool> {},
              utils::drag_and_mark_segments(std::vector<uint64_t> {}, 1, 0));
    EXPECT_EQ(std::vector<bool> {},
              utils::drag_and_mark_segments(std::vector<uint64_t> {}, 1, 1));
    EXPECT_EQ(std::vector<bool> {},
              utils::drag_and_mark_segments(std::vector<uint64_t> {}, 1, 100));

    EXPECT_EQ(std::vector<bool>({ 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 2, 0));
    EXPECT_EQ(std::vector<bool>({ 1 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 2, 1));
    EXPECT_EQ(std::vector<bool>({ 1 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 2, 100));

    EXPECT_EQ(std::vector<bool>({ 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 1, 0));
    EXPECT_EQ(std::vector<bool>({ 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 1, 1));
    EXPECT_EQ(std::vector<bool>({ 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 2 }), 1, 100));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 1, 1, 4, 6 }), 0, 0));
    EXPECT_EQ(std::vector<bool>({ 1, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 1, 1, 4, 6 }), 0, 1));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 1, 1, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 1, 1, 4, 6 }), 0, 4));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 1, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 1, 1, 4, 6 }), 0, 5));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 0, 1, 4, 6 }), 0, 0));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 0, 1, 4, 6 }), 0, 1));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 1, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 0, 1, 4, 6 }), 0, 2));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 1, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<uint64_t>({ 0, 1, 1, 4, 6 }), 0, 5));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 10, 0));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 1, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 10, 1));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 1, 1, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 10, 2));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 1, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 1, 10, 4, 6 }), 10, 5));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 6, 0));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 6, 1));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 0, 10, 4, 6 }), 6, 2));
    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 1, 10, 4, 6 }), 6, 50));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 6, 0, 10, 6, 4 }), 6, 0));
    EXPECT_EQ(std::vector<bool>({ 1, 0, 0, 1, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 6, 0, 10, 6, 4 }), 6, 1));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 0, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 6, 0, 10, 6, 4 }), 6, 2));
    EXPECT_EQ(std::vector<bool>({ 1, 1, 1, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 6, 1, 10, 6, 4 }), 6, 50));

    EXPECT_EQ(std::vector<bool>({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 6, 0, 10, 6, 4, 0, 0, 6, 0 }), 6, 0));
    EXPECT_EQ(std::vector<bool>({ 0, 1, 0, 0, 1, 0, 0, 0, 1, 0 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 6, 0, 10, 6, 4, 0, 0, 6, 0 }), 6, 1));
    EXPECT_EQ(std::vector<bool>({ 0, 1, 1, 0, 1, 1, 0, 0, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 6, 0, 10, 6, 4, 0, 0, 6, 0 }), 6, 2));
    EXPECT_EQ(std::vector<bool>({ 0, 1, 1, 1, 1, 1, 1, 0, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 6, 0, 10, 6, 4, 0, 0, 6, 0 }), 6, 3));
    EXPECT_EQ(std::vector<bool>({ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 }),
              utils::drag_and_mark_segments(std::vector<int>({ 0, 6, 0, 10, 6, 4, 0, 0, 6, 0 }), 6, 100));
}

} // namespace
