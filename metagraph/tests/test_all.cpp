#include "gtest/gtest.h"

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);

    // Earlier tests may leave library worker threads running in Release too.
    GTEST_FLAG_SET(death_test_style, "threadsafe");

    return RUN_ALL_TESTS();
}
