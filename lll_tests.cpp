#include "lll.h"
#include <gtest/gtest.h>

TEST(lll, LLL)
{
    Vector<double> v1(3);
    Vector<double> v2(3);
    Vector<double> v3(3);

    v1(0) = 0;
    v1(1) = 10;
    v1(2) = 2;

    v2(0) = 200;
    v2(1) = 314;
    v2(2) = 503;

    v3(0) = 20;
    v3(1) = 30;
    v3(2) = 50;

    std::vector<Vector<double>> basis;
    basis.push_back(v1);
    basis.push_back(v2);
    basis.push_back(v3);

    Vector<double> res(3);
    res(0) = 0;
    res(1) = 0;
    res(2) = 1;

    std::vector<Vector<double>> lll_basis = lll(0.75, basis);

    EXPECT_EQ(lll_basis[0], res);

    v1(0) = 1;
    v1(1) = 1;
    v1(2) = 1;

    v2(0) = -1;
    v2(1) = 0;
    v2(2) = 2;

    v3(0) = 3;
    v3(1) = 5;
    v3(2) = 6;

    basis = std::vector<Vector<double>>();
    basis.push_back(v1);
    basis.push_back(v2);
    basis.push_back(v3);

    Vector<double> res1(3);
    res1(0) = 0;
    res1(1) = 1;
    res1(2) = 0;

    Vector<double> res2(3);
    res2(0) = 1;
    res2(1) = 0;
    res2(2) = 1;

    Vector<double> res3(3);
    res3(0) = -1;
    res3(1) = 0;
    res3(2) = 2;

    lll_basis = lll(0.75, basis);

    EXPECT_EQ(lll_basis[0], res1);
    EXPECT_EQ(lll_basis[1], res2);
    EXPECT_EQ(lll_basis[2], res3);
}
