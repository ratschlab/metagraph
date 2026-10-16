#include "gtest/gtest.h"

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);

    // Earlier tests may leave background threads alive. Re-execute death tests
    // in a fresh process in Release builds as well as Debug builds.
    GTEST_FLAG_SET(death_test_style, "threadsafe");

    return RUN_ALL_TESTS();
}
