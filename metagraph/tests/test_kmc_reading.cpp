#include "gtest/gtest.h"

#include <unordered_set>
#include <string>

#include "kmc_parser.hpp"

using kmc::read_kmers;

const std::string kTestDataDir = "../tests/data";
const std::string kTestDumpBasename = kTestDataDir + "/transcripts_1000_kmc_counters";
const std::string kTestDumpBothStrandsBasename = kTestDumpBasename + "_both_strands";


TEST(kmc_parser, OpenFile) {
    {
        bool exception_happened = false;
        try {
            ASSERT_FALSE(exception_happened);
            read_kmers(kTestDumpBasename + ".kmc_pre", [&](std::string&& string) {
                EXPECT_EQ(11u, string.size());
            });
            ASSERT_FALSE(exception_happened);
        } catch (...) {
            exception_happened = true;
        }
        ASSERT_FALSE(exception_happened);
    }
    {
        bool exception_happened = false;
        try {
            ASSERT_FALSE(exception_happened);
            read_kmers(kTestDumpBasename + ".kmc_pre", [&](std::string&& string) {
                EXPECT_EQ(11u, string.size());
            });
            ASSERT_FALSE(exception_happened);
        } catch (...) {
            exception_happened = true;
        }
        ASSERT_FALSE(exception_happened);
    }
}

TEST(kmc_parser, OpenBadFile) {
    bool exception_happened = false;
    try {
        ASSERT_FALSE(exception_happened);
        read_kmers(kTestDumpBothStrandsBasename + "_invalid", {});
        ASSERT_TRUE(exception_happened);
    } catch (...) {
        exception_happened = true;
    }
    ASSERT_TRUE(exception_happened);
}

TEST(kmc_parser, ReadKmers) {
    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        });
        EXPECT_EQ(469983u, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 0);
        EXPECT_EQ(469983u, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 1);
        EXPECT_EQ(469983u, kmers.size());
    }
}

TEST(kmc_parser, ReadKmersMinCountThreshold) {
    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 1);
        EXPECT_EQ(469983u, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 2);
        EXPECT_EQ(255127u, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 3);
        EXPECT_EQ(177441u, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 1000);
        EXPECT_EQ(0u, kmers.size());
    }
}

TEST(kmc_parser, ReadKmersBothStrands) {
    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBothStrandsBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        });
        EXPECT_EQ(401460u * 2, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBothStrandsBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 0);
        EXPECT_EQ(401460u * 2, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBothStrandsBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 1);
        EXPECT_EQ(401460u * 2, kmers.size());
    }
}

TEST(kmc_parser, ReadKmersBothStrandsMinCountThreshold) {
    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBothStrandsBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 1);
        EXPECT_EQ(401460u * 2, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBothStrandsBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 2);
        EXPECT_EQ(238473u * 2, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBothStrandsBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 3);
        EXPECT_EQ(173157u * 2, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestDumpBothStrandsBasename, [&](std::string&& string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(std::move(string));
        }, 1000);
        EXPECT_EQ(0u, kmers.size());
    }
}
