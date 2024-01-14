#ifndef BIT_VECTOR_TESTSUITE_H
#define BIT_VECTOR_TESTSUITE_H

#include "gmds/utils/BitVector.h"
#include "gtest/gtest.h"
using namespace gmds;

TEST(BitVectorTestSuite, DefaultConstructor)
{
	gmds::BitVector bv;
	EXPECT_EQ(bv.size(), 0);
	EXPECT_EQ(bv.top(), 0);
	EXPECT_EQ(bv.capacity(), gmds::GChunkSize);
	EXPECT_TRUE(bv.empty());
}

TEST(BitVectorTestSuite, CopyConstructor)
{
	gmds::BitVector bv;
	gmds::BitVector bv_copy(bv);
	EXPECT_EQ(bv_copy.size(), 0);
	EXPECT_EQ(bv_copy.top(), 0);
	EXPECT_EQ(bv_copy.capacity(), gmds::GChunkSize);
	EXPECT_TRUE(bv_copy.empty());
}

TEST(BitVectorTestSuite, Resize)
{
	gmds::BitVector bv;
	bv.resize(10);
	EXPECT_EQ(bv.capacity(), 10);
}

TEST(BitVectorTestSuite, FillAll)
{
	gmds::BitVector bv;
	bv.resize(5);
	bv.fillAll();
	EXPECT_EQ(bv.size(), 5);
	EXPECT_EQ(bv.top(), 5);
}

// Ajoutez d'autres tests en fonction de vos fonctions

#endif     // BIT_VECTOR_TEST_H
