#include "Vector.h"
#include <gtest/gtest.h>

TEST(vector, equality)
{
    std::vector<double> ent1(2, 0);
    ent1[0] = 0;
    ent1[1] = 10;
    std::vector<double> ent2(3, 200);
    ent2[1] = 314;
    ent2[2] = 503;
    std::vector<double> ent3(3, 20);
    ent3[1] = 30;
    ent3[2] = 50;

    Vector<double> v1(ent1);
    Vector<double> v2(ent2);
    Vector<double> v3(ent3);
    Vector<double> v4(ent3);

    EXPECT_NE(v1, v2);
    EXPECT_NE(v2, v3);
    EXPECT_EQ(v3, v4);
}

TEST(vector, multiplication_by_scalar)
{
    std::vector<double> ent1(2, 0);
    ent1[1] = 10;
    std::vector<double> res1(2, 0);
    res1[1] = 30;
    std::vector<double> ent2(3, 200);
    ent2[1] = 314;
    ent2[2] = 503;
    std::vector<double> res2(3, -200);
    res2[1] = -314;
    res2[2] = -503;

    Vector<double> v1(ent1);
    Vector<double> r1(res1);
    Vector<double> v2(ent2);
    Vector<double> r2(res2);

    EXPECT_EQ(v1 * 3, r1);
    EXPECT_EQ(v2 * (-1), r2);
}

TEST(vector, multiplication_by_scalar_in_place)
{
    std::vector<double> ent1(2, 0);
    ent1[1] = 10;
    std::vector<double> res1(2, 0);
    res1[1] = 30;
    std::vector<double> ent2(3, 200);
    ent2[1] = 314;
    ent2[2] = 503;
    std::vector<double> res2(3, -200);
    res2[1] = -314;
    res2[2] = -503;

    Vector<double> v1(ent1);
    Vector<double> r1(res1);
    Vector<double> v2(ent2);
    Vector<double> r2(res2);

    v1 *= 3;
    v2 *= -1;

    EXPECT_EQ(v1, r1);
    EXPECT_EQ(v2, r2);
}

TEST(vector, inner_multiplication)
{
    std::vector<double> ent1(3, 0);
    ent1[0] = -2;
    ent1[1] = 10;
    std::vector<double> ent2(3, 200);
    ent2[1] = 314;
    ent2[2] = 503;

    Vector<double> v1(ent1);
    Vector<double> v2(ent2);

    EXPECT_EQ(v1 * v2, 2740);
    EXPECT_EQ(v2 * v1, 2740);
    EXPECT_EQ(v1 * v1, 104);
}

TEST(vector, addition)
{
    std::vector<double> ent1(3, 0);
    ent1[0] = -2;
    ent1[1] = 10;
    std::vector<double> ent2(3, 200);
    ent2[1] = 314;
    ent2[2] = 503;
    std::vector<double> res1(3, 0);
    res1[0] = -6;
    res1[1] = 30;
    std::vector<double> res2(3, 198);
    res2[1] = 324;
    res2[2] = 503;

    Vector<double> v1(ent1);
    Vector<double> v2(ent2);
    Vector<double> r1(res1);
    Vector<double> r2(res2);

    EXPECT_EQ(v1 + v1 + v1, r1);
    EXPECT_EQ(v1 + v2, r2);
}

TEST(vector, addition_in_place)
{
    std::vector<double> ent1(3, 0);
    ent1[0] = -2;
    ent1[1] = 10;
    std::vector<double> ent2(3, 200);
    ent2[1] = 314;
    ent2[2] = 503;
    std::vector<double> res1(3, 0);
    res1[0] = -4;
    res1[1] = 20;
    std::vector<double> res2(3, 198);
    res2[1] = 324;
    res2[2] = 503;

    Vector<double> v1(ent1);
    Vector<double> v2(ent2);
    Vector<double> r1(res1);
    Vector<double> r2(res2);

    v2 += v1;
    v1 += v1;

    EXPECT_EQ(v1, r1);
    EXPECT_EQ(v2, r2);
}

TEST(vector, negation)
{
    std::vector<double> ent1(2, 0);
    ent1[1] = 10;
    std::vector<double> res1(2, 0);
    res1[1] = -10;
    std::vector<double> ent2(3, 200);
    ent2[1] = 314;
    ent2[2] = 503;
    std::vector<double> res2(3, -200);
    res2[1] = -314;
    res2[2] = -503;

    Vector<double> v1(ent1);
    Vector<double> r1(res1);
    Vector<double> v2(ent2);
    Vector<double> r2(res2);

    EXPECT_EQ(-v1, r1);
    EXPECT_EQ(-v2, r2);
}

TEST(vector, substraction)
{
    std::vector<double> ent1(3, 0);
    ent1[0] = -2;
    ent1[1] = 10;
    std::vector<double> ent2(3, 200);
    ent2[1] = 314;
    ent2[2] = 503;
    std::vector<double> res1(3, 0);
    std::vector<double> res2(3, 202);
    res2[1] = 304;
    res2[2] = 503;

    Vector<double> v1(ent1);
    Vector<double> r1(res1);
    Vector<double> v2(ent2);
    Vector<double> r2(res2);

    EXPECT_EQ(v1 - v1, r1);
    EXPECT_EQ(v2 - v1, r2);
}

TEST(vector, substraction_in_place)
{
    std::vector<double> ent1(3, 0);
    ent1[0] = -2;
    ent1[1] = 10;
    std::vector<double> ent2(3, 200);
    ent2[1] = 314;
    ent2[2] = 503;
    std::vector<double> res1(3, 0);
    std::vector<double> res2(3, 202);
    res2[1] = 304;
    res2[2] = 503;

    Vector<double> v1(ent1);
    Vector<double> r1(res1);
    Vector<double> v2(ent2);
    Vector<double> r2(res2);

    v2 -= v1;
    v1 -= v1;

    EXPECT_EQ(v1, r1);
    EXPECT_EQ(v2, r2);
}