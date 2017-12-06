#include "gtest/gtest.h"

#include "datatypes.hpp"


template <class T>
void test_bit_vector() {
    bit_vector *vector = new T();
    ASSERT_TRUE(vector);
    delete vector;

    vector = new T(0, 1);
    ASSERT_TRUE(vector);
    delete vector;

    vector = new T(10, 0);
    ASSERT_TRUE(vector);
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(0, (*vector)[i]);
        // EXPECT_EQ(0u, vector->select1(i));
        // EXPECT_EQ(0u, vector->rank1(i));
    }
    delete vector;

    vector = new T(10, 1);
    ASSERT_TRUE(vector);
    EXPECT_EQ(10u, vector->size());
    for (size_t i = 0; i < vector->size(); ++i) {
        EXPECT_EQ(1, (*vector)[i]);
        EXPECT_EQ(i, vector->select1(i));
        EXPECT_EQ(i + 1, vector->rank1(i));
    }
    delete vector;

    // vector = new T({ 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1 });
    // ASSERT_TRUE(vector);
    // for (size_t i = 0; i < vector->size(); ++i) {
    //     EXPECT_EQ(i + 1, vector->rank1(vector->select1(i)));
    //     // EXPECT_EQ(i, vector->select1(vector->rank1(i)));
    // }
    // delete vector;
}


TEST(bit_vector_dyn, methods) {
    test_bit_vector<bit_vector_dyn>();
}

TEST(bit_vector_stat, methods) {
    test_bit_vector<bit_vector_stat>();
}
