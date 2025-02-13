#include "gtest/gtest.h"

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);

#ifndef NDEBUG
    GTEST_FLAG_SET(death_test_style, "fast");
#endif

    return RUN_ALL_TESTS();
}
