#include "gtest/gtest.h"

#include <string>
#include <filesystem>

#include "seq_io/sequence_io.hpp"


namespace {

using namespace mtg;
using namespace mtg::seq_io;

const std::string test_data_dir = "../tests/data";
const std::string test_fasta = test_data_dir + "/transcripts_1000.fa";
const std::string dump_filename = test_data_dir + "/dump.fasta.gz";


TEST(FastaFile, iterator_read) {
    size_t num_records = 0;
    size_t total_size = 0;
    for (const auto &record : FastaParser(test_fasta)) {
        num_records++;
        total_size += record.seq.l;
    }
    EXPECT_EQ(1000u, num_records);
    EXPECT_EQ(1490627u, total_size);
}

TEST(FastaFile, full_iterator_read_empty) {
    {
        FastaWriter writer(dump_filename, "seq", true);
    }

    size_t num_records = 0;
    size_t total_size = 0;
    for (const auto &record : FastaParser(dump_filename)) {
        num_records++;
        total_size += record.seq.l;
    }
    EXPECT_EQ(0u, num_records);
    EXPECT_EQ(0u, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFile, full_iterator_read_1K) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 1'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    size_t num_records = 0;
    size_t total_size = 0;
    for (const auto &record : FastaParser(dump_filename)) {
        num_records++;
        total_size += record.seq.l;
    }
    EXPECT_EQ(1'000u, num_records);
    EXPECT_EQ(499'500u, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFile, full_iterator_read_100K) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 100'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    size_t num_records = 0;
    size_t total_size = 0;
    for (const auto &record : FastaParser(dump_filename)) {
        num_records++;
        total_size += record.seq.l;
    }
    EXPECT_EQ(100'000u, num_records);
    EXPECT_EQ(49'950'000u, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFile, full_iterator_read_100K_async) {
    {
        FastaWriter writer(dump_filename, "seq", true, true);
        for (size_t i = 0; i < 100'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    size_t num_records = 0;
    size_t total_size = 0;
    for (const auto &record : FastaParser(dump_filename)) {
        num_records++;
        total_size += record.seq.l;
    }
    EXPECT_EQ(100'000u, num_records);
    EXPECT_EQ(49'950'000u, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFile, full_iterator_read_100K_multithreaded) {
    {
        FastaWriter writer(dump_filename, "seq", true, true);
        #pragma omp parallel for num_threads(2)
        for (size_t i = 0; i < 100'000; ++i) {
            std::string seq(i % 1'000, 'A');

            #pragma omp critical
            writer.write(std::move(seq));
        }
    }

    size_t num_records = 0;
    size_t total_size = 0;
    for (const auto &record : FastaParser(dump_filename)) {
        num_records++;
        total_size += record.seq.l;
    }
    EXPECT_EQ(100'000u, num_records);
    EXPECT_EQ(49'950'000u, total_size);

    std::filesystem::remove(dump_filename);
}

// test that seek works fast so we can copy an iterator very quickly
TEST(FastaFile, full_iterator_read_100K_fast_copy) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 100'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    FastaParser fasta_parser(dump_filename);
    FastaParser::iterator copy;

    size_t total_size = 0;

    // on each iteration, copy the current iterator 50 times
    for (auto it = fasta_parser.begin(); it != fasta_parser.end(); ++it) {
        for (size_t i = 0; i < 20; ++i) {
            copy = it;
            total_size += copy->seq.l;
        }
    }

    EXPECT_EQ(49'950'000u * 20, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFile, twice_iterator_read_empty) {
    {
        FastaWriter writer(dump_filename, "seq", true);
    }

    FastaParser parser(dump_filename);
    auto begin = parser.begin();

    size_t num_records = 0;
    size_t total_size = 0;
    for (auto it = begin; it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(0u, num_records);
    EXPECT_EQ(0u, total_size);

    num_records = 0;
    total_size = 0;
    for (auto it = begin; it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(0u, num_records);
    EXPECT_EQ(0u, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFile, twice_full_iterator_read_1K) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 1'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    FastaParser parser(dump_filename);
    auto begin = parser.begin();

    size_t num_records = 0;
    size_t total_size = 0;
    for (auto it = begin; it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(1'000u, num_records);
    EXPECT_EQ(499'500u, total_size);

    num_records = 0;
    total_size = 0;
    for (auto it = begin; it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(1'000u, num_records);
    EXPECT_EQ(499'500u, total_size);

    std::filesystem::remove(dump_filename);
}


TEST(FastaFile, twice_iterator_read_empty_with_move) {
    {
        FastaWriter writer(dump_filename, "seq", true);
    }

    FastaParser parser(dump_filename);
    auto begin = parser.begin();

    size_t num_records = 0;
    size_t total_size = 0;
    for (auto it = begin; it != parser.end(); ++it ) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(0u, num_records);
    EXPECT_EQ(0u, total_size);

    num_records = 0;
    total_size = 0;
    for (auto it = std::move(begin); it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(0u, num_records);
    EXPECT_EQ(0u, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFile, twice_full_iterator_read_1K_with_move) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 1'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    FastaParser parser(dump_filename);
    auto begin = parser.begin();

    size_t num_records = 0;
    size_t total_size = 0;
    for (auto it = begin; it != parser.end(); ++it ) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(1'000u, num_records);
    EXPECT_EQ(499'500u, total_size);

    num_records = 0;
    total_size = 0;
    for (auto it = std::move(begin); it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(1'000u, num_records);
    EXPECT_EQ(499'500u, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFile, iterator_compare) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 1'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    FastaParser parser(dump_filename);

    for (auto it = parser.begin(); it != parser.end(); ++it) {
        auto copy = it;
        EXPECT_TRUE(copy == it);
        EXPECT_FALSE(copy != it);
        EXPECT_TRUE(++copy != it);
        EXPECT_FALSE(copy == it);
        EXPECT_TRUE(copy == ++it);
        EXPECT_FALSE(copy != it);
    }

    std::filesystem::remove(dump_filename);
}

TEST(FastaFile, iterator_read_1M_multiple_copies) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        writer.write(std::string(1'000'000, 'A'));
        writer.write(std::string(2'000'000, 'C'));
        writer.write(std::string(3'000'000, 'G'));
        writer.write(std::string(4'000'000, 'T'));
    }

    FastaParser parser(dump_filename);

    size_t num_records = 0;
    size_t total_size = 0;

    auto begin = parser.begin();
    auto first_copy = begin;
    auto second_copy = begin;

    for ( ; first_copy != parser.end(); ++first_copy) {
        num_records++;
        total_size += first_copy->seq.l;
    }

    EXPECT_EQ(4u, num_records);
    EXPECT_EQ(10'000'000u, total_size);

    for ( ; second_copy != parser.end(); ++second_copy) {
        num_records++;
        total_size += second_copy->seq.l;
    }

    EXPECT_EQ(8u, num_records);
    EXPECT_EQ(20'000'000u, total_size);

    ++begin;
    first_copy = begin;
    second_copy = begin;

    num_records = 0;
    total_size = 0;

    for ( ; first_copy != parser.end(); ++first_copy) {
        num_records++;
        total_size += first_copy->seq.l;
    }

    EXPECT_EQ(3u, num_records);
    EXPECT_EQ(9'000'000u, total_size);

    for ( ; second_copy != parser.end(); ++second_copy) {
        num_records++;
        total_size += second_copy->seq.l;
    }

    EXPECT_EQ(6u, num_records);
    EXPECT_EQ(18'000'000u, total_size);

    ++begin;
    first_copy = begin;
    second_copy = begin;

    num_records = 0;
    total_size = 0;

    for ( ; first_copy != parser.end(); ++first_copy) {
        num_records++;
        total_size += first_copy->seq.l;
    }

    EXPECT_EQ(2u, num_records);
    EXPECT_EQ(7'000'000u, total_size);

    for ( ; second_copy != parser.end(); ++second_copy) {
        num_records++;
        total_size += second_copy->seq.l;
    }

    EXPECT_EQ(4u, num_records);
    EXPECT_EQ(14'000'000u, total_size);

    ++begin;
    first_copy = begin;
    second_copy = begin;

    num_records = 0;
    total_size = 0;

    for ( ; first_copy != parser.end(); ++first_copy) {
        num_records++;
        total_size += first_copy->seq.l;
    }

    EXPECT_EQ(1u, num_records);
    EXPECT_EQ(4'000'000u, total_size);

    for ( ; second_copy != parser.end(); ++second_copy) {
        num_records++;
        total_size += second_copy->seq.l;
    }

    EXPECT_EQ(2u, num_records);
    EXPECT_EQ(8'000'000u, total_size);

    EXPECT_EQ(4'000'000u, begin->seq.l);

    std::filesystem::remove(dump_filename);
}


TEST(FastaFileWithCanonical, iterator_read) {
    size_t num_records = 0;
    size_t total_size = 0;
    for (const auto &record : FastaParser(test_fasta, true)) {
        num_records++;
        total_size += record.seq.l;
    }
    EXPECT_EQ(1000u * 2, num_records);
    EXPECT_EQ(1490627u * 2, total_size);
}

TEST(FastaFileWithCanonical, full_iterator_read_empty) {
    {
        FastaWriter writer(dump_filename, "seq", true);
    }

    size_t num_records = 0;
    size_t total_size = 0;
    for (const auto &record : FastaParser(dump_filename, true)) {
        num_records++;
        total_size += record.seq.l;
    }
    EXPECT_EQ(0u, num_records);
    EXPECT_EQ(0u, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFileWithCanonical, full_iterator_read_1K) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 1'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    size_t num_records = 0;
    size_t total_size = 0;
    for (const auto &record : FastaParser(dump_filename, true)) {
        num_records++;
        total_size += record.seq.l;
    }
    EXPECT_EQ(1'000u * 2, num_records);
    EXPECT_EQ(499'500u * 2, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFileWithCanonical, full_iterator_read_100K) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 100'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    size_t num_records = 0;
    size_t total_size = 0;
    for (const auto &record : FastaParser(dump_filename, true)) {
        num_records++;
        total_size += record.seq.l;
    }
    EXPECT_EQ(100'000u * 2, num_records);
    EXPECT_EQ(49'950'000u * 2, total_size);

    std::filesystem::remove(dump_filename);
}

// test that seek works fast so we can copy an iterator very quickly
TEST(FastaFileWithCanonical, full_iterator_read_100K_fast_copy) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 100'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    FastaParser fasta_parser(dump_filename, true);
    FastaParser::iterator copy;

    size_t total_size = 0;

    // on each iteration, copy the current iterator 50 times
    for (auto it = fasta_parser.begin(); it != fasta_parser.end(); ++it) {
        for (size_t i = 0; i < 20; ++i) {
            copy = it;
            total_size += copy->seq.l;
        }
    }

    EXPECT_EQ(49'950'000u * 20 * 2, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFileWithCanonical, twice_iterator_read_empty) {
    {
        FastaWriter writer(dump_filename, "seq", true);
    }

    FastaParser parser(dump_filename, true);
    auto begin = parser.begin();

    size_t num_records = 0;
    size_t total_size = 0;
    for (auto it = begin; it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(0u, num_records);
    EXPECT_EQ(0u, total_size);

    num_records = 0;
    total_size = 0;
    for (auto it = begin; it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(0u, num_records);
    EXPECT_EQ(0u, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFileWithCanonical, twice_full_iterator_read_1K) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 1'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    FastaParser parser(dump_filename, true);
    auto begin = parser.begin();

    size_t num_records = 0;
    size_t total_size = 0;
    for (auto it = begin; it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(1'000u * 2, num_records);
    EXPECT_EQ(499'500u * 2, total_size);

    num_records = 0;
    total_size = 0;
    for (auto it = begin; it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(1'000u * 2, num_records);
    EXPECT_EQ(499'500u * 2, total_size);

    std::filesystem::remove(dump_filename);
}


TEST(FastaFileWithCanonical, twice_iterator_read_empty_with_move) {
    {
        FastaWriter writer(dump_filename, "seq", true);
    }

    FastaParser parser(dump_filename, true);
    auto begin = parser.begin();

    size_t num_records = 0;
    size_t total_size = 0;
    for (auto it = begin; it != parser.end(); ++it ) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(0u, num_records);
    EXPECT_EQ(0u, total_size);

    num_records = 0;
    total_size = 0;
    for (auto it = std::move(begin); it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(0u, num_records);
    EXPECT_EQ(0u, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFileWithCanonical, twice_full_iterator_read_1K_with_move) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 1'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    FastaParser parser(dump_filename, true);
    auto begin = parser.begin();

    size_t num_records = 0;
    size_t total_size = 0;
    for (auto it = begin; it != parser.end(); ++it ) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(1'000u * 2, num_records);
    EXPECT_EQ(499'500u * 2, total_size);

    num_records = 0;
    total_size = 0;
    for (auto it = std::move(begin); it != parser.end(); ++it) {
        num_records++;
        total_size += it->seq.l;
    }
    EXPECT_EQ(1'000u * 2, num_records);
    EXPECT_EQ(499'500u * 2, total_size);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFileWithCanonical, iterator_compare) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        for (size_t i = 0; i < 1'000; ++i) {
            writer.write(std::string(i % 1'000, 'A'));
        }
    }

    FastaParser parser(dump_filename, true);

    for (auto it = parser.begin(); it != parser.end(); ++it) {
        auto copy = it;
        EXPECT_TRUE(copy == it);
        EXPECT_FALSE(copy != it);
        EXPECT_TRUE(++copy != it);
        EXPECT_FALSE(copy == it);
        EXPECT_TRUE(copy == ++it);
        EXPECT_FALSE(copy != it);
    }

    std::filesystem::remove(dump_filename);
}

TEST(FastaFileWithCanonical, iterator_read_1M_multiple_copies) {
    {
        FastaWriter writer(dump_filename, "seq", true);
        writer.write(std::string(1'000'000, 'A'));
        writer.write(std::string(2'000'000, 'C'));
        writer.write(std::string(3'000'000, 'G'));
        writer.write(std::string(4'000'000, 'T'));
    }

    FastaParser parser(dump_filename, true);

    size_t num_records = 0;
    size_t total_size = 0;

    auto begin = parser.begin();
    auto first_copy = begin;
    auto second_copy = begin;

    for ( ; first_copy != parser.end(); ++first_copy) {
        num_records++;
        total_size += first_copy->seq.l;
    }

    EXPECT_EQ(8u, num_records);
    EXPECT_EQ(20'000'000u, total_size);

    for ( ; second_copy != parser.end(); ++second_copy) {
        num_records++;
        total_size += second_copy->seq.l;
    }

    EXPECT_EQ(8u * 2, num_records);
    EXPECT_EQ(20'000'000u * 2, total_size);

    ++begin;
    first_copy = begin;
    second_copy = begin;

    num_records = 0;
    total_size = 0;

    for ( ; first_copy != parser.end(); ++first_copy) {
        num_records++;
        total_size += first_copy->seq.l;
    }

    EXPECT_EQ(7u, num_records);
    EXPECT_EQ(19'000'000u, total_size);

    for ( ; second_copy != parser.end(); ++second_copy) {
        num_records++;
        total_size += second_copy->seq.l;
    }

    EXPECT_EQ(7u * 2, num_records);
    EXPECT_EQ(19'000'000u * 2, total_size);

    ++begin;
    first_copy = begin;
    second_copy = begin;

    num_records = 0;
    total_size = 0;

    for ( ; first_copy != parser.end(); ++first_copy) {
        num_records++;
        total_size += first_copy->seq.l;
    }

    EXPECT_EQ(6u, num_records);
    EXPECT_EQ(18'000'000u, total_size);

    for ( ; second_copy != parser.end(); ++second_copy) {
        num_records++;
        total_size += second_copy->seq.l;
    }

    EXPECT_EQ(6u * 2, num_records);
    EXPECT_EQ(18'000'000u * 2, total_size);

    ++begin;
    first_copy = begin;
    second_copy = begin;

    num_records = 0;
    total_size = 0;

    for ( ; first_copy != parser.end(); ++first_copy) {
        num_records++;
        total_size += first_copy->seq.l;
    }

    EXPECT_EQ(5u, num_records);
    EXPECT_EQ(16'000'000u, total_size);

    for ( ; second_copy != parser.end(); ++second_copy) {
        num_records++;
        total_size += second_copy->seq.l;
    }

    EXPECT_EQ(5u * 2, num_records);
    EXPECT_EQ(16'000'000u * 2, total_size);

    EXPECT_EQ(2'000'000u, begin->seq.l);

    std::filesystem::remove(dump_filename);
}

TEST(FastaFromString, read_fasta_from_string) {
    std::string fasta_str = "";

    int nr_seqs = 20'000;
    for (int i = 0; i < nr_seqs; i++) {
        fasta_str.append(">hello\nAAA\n");
    }

    // precondition for test, want to test for larger fasta strings
    EXPECT_GT(fasta_str.size(), 65536);

    int seqs_cnt = 0;
    read_fasta_from_string(fasta_str,
                           [&](kseq_t*) {
                               seqs_cnt++;
                           });

    EXPECT_EQ(seqs_cnt, nr_seqs);
}

} // namespace
