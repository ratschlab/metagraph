#include "gtest/gtest.h"

#define private public
#define protected public

#include "annotate_column_compressed.hpp"
#include "utils.hpp"
#include "bit_vector.hpp"

// Disable death tests
#ifndef _DEATH_TEST
#ifdef EXPECT_DEATH
#undef EXPECT_DEATH
#define EXPECT_DEATH(a, b) (void)0
#endif
#endif

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

utils::RowsFromColumnsTransformer generate_rct_file() {
    annotate::ColumnCompressed<> annotation(6);

    annotation.set_labels(1, { "Label0", "Label1" });
    annotation.set_labels(2, { "Label0", "Label1" });
    annotation.set_labels(3, { "Label1" });
    annotation.set_labels(4, { "Label1" });
    annotation.set_labels(5, { "Label0", "Label1" });

    annotation.dump_columns(test_dump_basename);

    utils::RowsFromColumnsTransformer rct(7, {
        test_dump_basename + ".0.raw.column.annodbg",
        test_dump_basename + ".1.raw.column.annodbg"
    });

    return rct;
}

std::vector<bit_vector_small> generate_rows() {
    std::vector<bit_vector_small> rows;
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(rows),
        [](const auto &a) {
            return bit_vector_small(a);
        });

    return rows;
}

std::vector<bit_vector_small const*>
generate_ptrs(const std::vector<bit_vector_small> &rows) {
    std::vector<bit_vector_small const*> rows_ptr;
    std::transform(rows.begin(), rows.end(), std::back_inserter(rows_ptr),
        [](auto &a) {
            return &a;
        });

    return rows_ptr;
}

void check_rows(utils::RowsFromColumnsTransformer&& rct) {
    ASSERT_EQ(2u, rct.columns());
    ASSERT_EQ(7u, rct.rows());
    ASSERT_EQ(8u, rct.values_left());
    // ASSERT_EQ(std::vector<uint64_t>({ 3, 5 }), rct.num_set_bits());

    uint64_t i = 0;
    utils::call_rows([&](auto&& row_indices) {
        sdsl::bit_vector bv(rct.columns());
        for (auto j : row_indices) {
            bv[j] = 1;
        }
        EXPECT_EQ(matrix[i++], bv) << i;
    }, std::move(rct));

    ASSERT_EQ(7u, i);
    EXPECT_EQ(0u, rct.values_left());
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

TEST(Utils, RowsFromColumnsTransformerCallRowsFile) {
    auto rct = generate_rct_file();
    check_rows(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallIndicesFile) {
    auto rct = generate_rct_file();
    check_indices(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallRowsColumns) {
    auto rows = generate_rows();
    utils::RowsFromColumnsTransformer rct(generate_ptrs(rows));
    check_rows(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallIndicesColumns) {
    auto rows = generate_rows();
    utils::RowsFromColumnsTransformer rct(generate_ptrs(rows));
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

TEST(Vector, ReserveInfinityCheckThrow) {
    Vector<int> vector;
    EXPECT_THROW(vector.reserve(2llu << 60), std::bad_alloc);
}

TEST(Vector, ResizeInfinityCheckThrow) {
    Vector<int> vector;
    EXPECT_THROW(vector.resize(2llu << 60), std::bad_alloc);
}
