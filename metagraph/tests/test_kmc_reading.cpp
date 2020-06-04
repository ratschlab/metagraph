#include "gtest/gtest.h"

#include <unordered_set>
#include <string>

#include "seq_io/kmc_parser.hpp"


namespace {

using namespace mtg;

using mtg::seq_io::read_kmers;

const std::string kTestDataDir = "../tests/data";
// constructed with '-b' flag in KMC
const std::string kTestKMCDatabase = kTestDataDir + "/transcripts_1000_kmc_counters";
// constructed without '-b' flag in KMC
const std::string kTestKMCDatabaseCanonical = kTestKMCDatabase + "_both_strands";


TEST(kmc_parser, OpenFile) {
    for (const auto &file : { kTestKMCDatabase + ".kmc_pre",
                              kTestKMCDatabase + ".kmc_suf",
                              kTestKMCDatabaseCanonical + ".kmc_pre",
                              kTestKMCDatabaseCanonical + ".kmc_suf" }) {
        bool exception_happened = false;
        try {
            ASSERT_FALSE(exception_happened);
            read_kmers(file, [&](std::string_view string) {
                EXPECT_EQ(11u, string.size());
            }, false);
            ASSERT_FALSE(exception_happened);
            read_kmers(file, [&](std::string_view string) {
                EXPECT_EQ(11u, string.size());
            }, true);
            ASSERT_FALSE(exception_happened);
        } catch (...) {
            exception_happened = true;
        }
        ASSERT_FALSE(exception_happened);
    }
}

TEST(kmc_parser, OpenBadFile) {
    for (const auto &file : { kTestKMCDatabase + ".kmc_preeee",
                              kTestKMCDatabase + ".kmc_s",
                              kTestKMCDatabase + "_invalid",
                              kTestKMCDatabaseCanonical + ".kmc_preeee",
                              kTestKMCDatabaseCanonical + ".kmc_s",
                              kTestKMCDatabaseCanonical + "_invalid" }) {
        bool exception_happened = false;
        try {
            read_kmers(file, [&](std::string_view) {}, false);
            ASSERT_FALSE(exception_happened);
        } catch (...) {
            exception_happened = true;
        }
        ASSERT_TRUE(exception_happened);

        exception_happened = false;
        try {
            read_kmers(file, [&](std::string_view) {}, true);
            ASSERT_FALSE(exception_happened);
        } catch (...) {
            exception_happened = true;
        }
        ASSERT_TRUE(exception_happened);
    }
}

TEST(kmc_parser, ReadKmers) {
    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestKMCDatabase, [&](std::string_view string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(string);
        }, false);
        EXPECT_EQ(469983u, kmers.size());
    }

    for (size_t min_count : { 0, 1 }) {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestKMCDatabase, [&](std::string_view string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(string);
        }, false, min_count);
        EXPECT_EQ(469983u, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestKMCDatabase, [&](std::string_view string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(string);
        }, true);
        EXPECT_EQ(469983u, kmers.size());
    }

    for (size_t min_count : { 0, 1 }) {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestKMCDatabase, [&](std::string_view string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(string);
        }, true, min_count);
        EXPECT_EQ(469983u, kmers.size());
    }
}

TEST(kmc_parser, ReadKmersMinCountThreshold) {
    std::vector<std::pair<size_t, size_t>> values = {
        { 1, 469983u },
        { 2, 255127u },
        { 3, 177441u },
        { 1000, 0u },
    };

    for (bool call_both : { false, true }) {
        for (const auto &min_count__num_kmers : values) {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabase, [&](std::string_view string) {
                    EXPECT_EQ(11u, string.size());
                    kmers.emplace(string);
                },
                call_both,
                min_count__num_kmers.first
            );
            EXPECT_EQ(min_count__num_kmers.second, kmers.size());
        }
    }
}

TEST(kmc_parser, ReadKmersMaxCountThreshold) {
    for (bool call_both : { false, true }) {
        std::vector<std::pair<size_t, size_t>> values = {
            { 1, 469983u },
            { 2, 255127u },
            { 3, 177441u },
            { 1000, 0u },
        };

        for (const auto &min_count__num_kmers : values) {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabase, [&](std::string_view string) {
                    EXPECT_EQ(11u, string.size());
                    kmers.emplace(string);
                },
                call_both,
                min_count__num_kmers.first,
                min_count__num_kmers.first
            );
            EXPECT_EQ(0u, kmers.size());
        }

        for (const auto &min_count__num_kmers : values) {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabase, [&](std::string_view string) {
                    EXPECT_EQ(11u, string.size());
                    kmers.emplace(string);
                },
                call_both,
                1000,
                min_count__num_kmers.first
            );
            EXPECT_EQ(0u, kmers.size());
        }

        for (const auto &min_count__num_kmers : values) {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabase, [&](std::string_view string) {
                    EXPECT_EQ(11u, string.size());
                    kmers.emplace(string);
                },
                call_both,
                min_count__num_kmers.first,
                1000
            );
            EXPECT_EQ(min_count__num_kmers.second, kmers.size());
        }

        {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabase, [&](std::string_view string) {
                EXPECT_EQ(11u, string.size());
                kmers.emplace(string);
            }, call_both, 1, 2);
            EXPECT_EQ((469983u - 255127u), kmers.size());
        }

        {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabase, [&](std::string_view string) {
                EXPECT_EQ(11u, string.size());
                kmers.emplace(string);
            }, call_both, 2, 3);
            EXPECT_EQ((255127u - 177441u), kmers.size());
        }

        {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabase, [&](std::string_view string) {
                EXPECT_EQ(11u, string.size());
                kmers.emplace(string);
            }, call_both, 1, 3);
            EXPECT_EQ((469983u - 177441u), kmers.size());
        }
    }
}

TEST(kmc_parser, ReadKmersBothStrands) {
    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(string);
        }, false);
        EXPECT_EQ(401460u, kmers.size());
    }

    for (size_t min_count : { 0, 1 }) {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(string);
        }, false, min_count);
        EXPECT_EQ(401460u, kmers.size());
    }

    {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(string);
        }, true);
        EXPECT_EQ(401460u * 2, kmers.size());
    }

    for (size_t min_count : { 0, 1 }) {
        std::unordered_set<std::string> kmers;
        read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
            EXPECT_EQ(11u, string.size());
            kmers.emplace(string);
        }, true, min_count);
        EXPECT_EQ(401460u * 2, kmers.size());
    }
}

TEST(kmc_parser, ReadKmersBothStrandsMinCountThreshold) {
    std::vector<std::pair<size_t, size_t>> values = {
        { 1, 401460u },
        { 2, 238473u },
        { 3, 173157u },
        { 1000, 0u },
    };

    for (bool call_both : { false, true }) {
        for (const auto &min_count__num_kmers : values) {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
                    EXPECT_EQ(11u, string.size());
                    kmers.emplace(string);
                },
                call_both,
                min_count__num_kmers.first
            );
            EXPECT_EQ(min_count__num_kmers.second * (1 + call_both), kmers.size());
        }
    }
}

TEST(kmc_parser, ReadKmersBothStrandsMaxCountThreshold) {
    for (bool call_both : { false, true }) {
        std::vector<std::pair<size_t, size_t>> values = {
            { 1, 401460u },
            { 2, 238473u },
            { 3, 173157u },
            { 1000, 0u },
        };

        for (const auto &min_count__num_kmers : values) {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
                    EXPECT_EQ(11u, string.size());
                    kmers.emplace(string);
                },
                call_both,
                min_count__num_kmers.first,
                min_count__num_kmers.first
            );
            EXPECT_EQ(0u, kmers.size());
        }

        for (const auto &min_count__num_kmers : values) {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
                    EXPECT_EQ(11u, string.size());
                    kmers.emplace(string);
                },
                call_both,
                1000,
                min_count__num_kmers.first
            );
            EXPECT_EQ(0u, kmers.size());
        }

        for (const auto &min_count__num_kmers : values) {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
                    EXPECT_EQ(11u, string.size());
                    kmers.emplace(string);
                },
                call_both,
                min_count__num_kmers.first,
                1000
            );
            EXPECT_EQ(min_count__num_kmers.second * (call_both + 1), kmers.size());
        }

        {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
                EXPECT_EQ(11u, string.size());
                kmers.emplace(string);
            }, call_both, 1, 2);
            EXPECT_EQ((401460u - 238473u) * (call_both + 1), kmers.size());
        }

        {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
                EXPECT_EQ(11u, string.size());
                kmers.emplace(string);
            }, call_both, 2, 3);
            EXPECT_EQ((238473u - 173157u) * (call_both + 1), kmers.size());
        }

        {
            std::unordered_set<std::string> kmers;
            read_kmers(kTestKMCDatabaseCanonical, [&](std::string_view string) {
                EXPECT_EQ(11u, string.size());
                kmers.emplace(string);
            }, call_both, 1, 3);
            EXPECT_EQ((401460u - 173157u) * (call_both + 1), kmers.size());
        }
    }
}

} // namespace
