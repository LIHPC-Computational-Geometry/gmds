/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <gmds/math/TransfiniteInterpolation.h>

/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
TEST(TransfiniteInterpolationTest, testUnitSquare) {
    std::vector<std::vector<math::Point> > pnts;
    auto I=5;
    auto J=5;
    pnts.resize(I);
    for(auto &line_p:pnts)
        line_p.resize(J);

    math::Point c00(0,0,0);
    math::Point cI0(1,0,0);
    math::Point c0J(0,1,0);
    math::Point cIJ(1,1,0);

    for(auto i=0;i<I;i++){
        auto w0 = (double)(I-1-i)/(double)(I-1);
        auto w1 = (double)(i)/(double)(I-1);
        pnts[i][0  ] = w0*c00+w1*cI0;
        pnts[i][J-1] = w0*c0J+w1*cIJ;
    }
    for(auto j=0;j<J;j++){
        auto w0 = (double)(J-1-j)/(double)(J-1);
        auto w1 = (double)(j)/(double)(J-1);
        pnts[0  ][j] = w0*c00+w1*c0J;
        pnts[I-1][j] = w0*cI0+w1*cIJ;
    }

    EXPECT_TRUE(math::TransfiniteInterpolation::compute(pnts));

    EXPECT_DOUBLE_EQ(pnts[1][1].X(),0.25);
    EXPECT_DOUBLE_EQ(pnts[1][1].Y(),0.25);
    EXPECT_DOUBLE_EQ(pnts[1][1].Z(),0.00);

    EXPECT_DOUBLE_EQ(pnts[2][2].X(),0.5);
    EXPECT_DOUBLE_EQ(pnts[2][2].Y(),0.5);
    EXPECT_DOUBLE_EQ(pnts[2][2].Z(),0.00);
}
/*----------------------------------------------------------------------------*/
TEST(TransfiniteInterpolationTest, testUnitSquare2) {
    auto N = 5;
    Array2D<math::Point> pnts(N, N);
    auto I=5;
    auto J=5;

    math::Point c00(0,0,0);
    math::Point cI0(1,0,0);
    math::Point c0J(0,1,0);
    math::Point cIJ(1,1,0);

    for(auto i=0;i<I;i++){
        auto w0 = (double)(I-1-i)/(double)(I-1);
        auto w1 = (double)(i)/(double)(I-1);
        pnts(i,0  ) = w0*c00+w1*cI0;
        pnts(i,J-1) = w0*c0J+w1*cIJ;
    }
    for(auto j=0;j<J;j++){
        auto w0 = (double)(J-1-j)/(double)(J-1);
        auto w1 = (double)(j)/(double)(J-1);
        pnts(0  ,j) = w0*c00+w1*c0J;
        pnts(I-1,j) = w0*cI0+w1*cIJ;
    }

    EXPECT_TRUE(math::TransfiniteInterpolation::computeQuad(pnts));

    EXPECT_DOUBLE_EQ(pnts(1,1).X(),0.25);
    EXPECT_DOUBLE_EQ(pnts(1,1).Y(),0.25);
    EXPECT_DOUBLE_EQ(pnts(1,1).Z(),0.00);

    EXPECT_DOUBLE_EQ(pnts(2,2).X(),0.5);
    EXPECT_DOUBLE_EQ(pnts(2,2).Y(),0.5);
    EXPECT_DOUBLE_EQ(pnts(2,2).Z(),0.00);
}

/*----------------------------------------------------------------------------*/
TEST(TransfiniteInterpolationTest, testUnitTriangle) {
    auto N = 5;
    TriArray<math::Point> P(N);

    for (auto i = 0; i < N; i++) {
        P(i, 0) = math::Point(i, 0, 0);
    }
    for (auto j = 1; j < N; j++) {
        P(0, j) = math::Point(0, j, 0);
    }

    P(1, 3) = math::Point(1, 3, 0);
    P(2, 2) = math::Point(2, 2, 0);
    P(3, 1) = math::Point(3, 1, 0);

    math::TransfiniteInterpolation::computeTri(P);

    EXPECT_EQ(P(1, 1), math::Point(1, 1, 0));
    EXPECT_EQ(P(1, 2), math::Point(1, 2, 0));
    EXPECT_EQ(P(2, 1), math::Point(2, 1, 0));

    EXPECT_EQ(P(2, 1, 1), math::Point(1, 1, 0));
    EXPECT_EQ(P(1, 2, 1), math::Point(1, 2, 0));
    EXPECT_EQ(P(1, 1, 2), math::Point(2, 1, 0));
}

/*----------------------------------------------------------------------------*/
TEST(TransfiniteInterpolationTest, testIcosehdronFaceTriangle) {
    auto N = 5;
    TriArray<math::Point> P(N);

    P(0,0) =math::Point(0.525731, -0.850651, 0);
    P(0,1) =math::Point(0.681718, -0.716567, -0.147621);
    P(0,2) =math::Point(0.809017, -0.5, -0.309017);
    P(0,3) =math::Point(0.864188, -0.238856, -0.442863);
    P(0,4) =math::Point(0.850651, 0, -0.525731);

    P(1,0) =math::Point(0.442863, -0.864188, -0.238856);
    P(2,0) =math::Point(0.309017, -0.809017, -0.5);
    P(3,0) =math::Point (0.147621, -0.681718, -0.716567);
    P(4,0) =math::Point (0, -0.525731, -0.850651);

    P(1,3)=math::Point(0.716567, -0.147621, -0.681718);
    P(2,2)=math::Point(0.5, -0.309017, -0.809017);
    P(3,1)=math::Point(0.238856, -0.442863, -0.864188);

    math::TransfiniteInterpolation::computeTri(P);
    EXPECT_TRUE(P(2, 1, 1).distance(math::Point(0.605127, -0.708494, -0.437874))<1e-4);
    EXPECT_TRUE(P(1, 2, 1).distance(math::Point(0.708494, -0.437874, -0.605127))<1e-4);
    EXPECT_TRUE(P(1, 1, 2).distance(math::Point(0.437874, -0.605127, -0.708494))<1e-4);

}
