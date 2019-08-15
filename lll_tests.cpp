// Tests for the lll algorithm
#include "lll.h"
#include <gtest/gtest.h>

TEST(lll, GramSchmidt)
{
    Vector<double> v1(3);
    v1(0) = 1;
    v1(1) = 1;
    v1(2) = 0;

    Vector<double> v2(3);
    v2(0) = 2;
    v2(1) = 2;
    v2(2) = 3;

    Vector<double> res2(3);
    res2(0) = 0;
    res2(1) = 0;
    res2(2) = 3;

    Vector<double> u1(4);
    u1(0) = 1;
    u1(1) = 0;
    u1(2) = 0;
    u1(3) = 0;

    Vector<double> u2(4);
    u2(0) = 2;
    u2(1) = 1;
    u2(2) = 0;
    u2(3) = 0;

    Vector<double> u3(4);
    u3(0) = 17;
    u3(1) = 12.3;
    u3(2) = 1;
    u3(3) = 0;

    Vector<double> u4(4);
    u4(0) = 2;
    u4(1) = 1;
    u4(2) = -12;
    u4(3) = 1;

    Vector<double> e1(4);
    e1(0) = 1;
    e1(1) = 0;
    e1(2) = 0;
    e1(3) = 0;

    Vector<double> e2(4);
    e2(0) = 0;
    e2(1) = 1;
    e2(2) = 0;
    e2(3) = 0;

    Vector<double> e3(4);
    e3(0) = 0;
    e3(1) = 0;
    e3(2) = 1;
    e3(3) = 0;

    Vector<double> e4(4);
    e4(0) = 0;
    e4(1) = 0;
    e4(2) = 0;
    e4(3) = 1;

    std::vector<Vector<double>> sub_basis;
    sub_basis.push_back(v1);
    sub_basis.push_back(v2);

    std::vector<Vector<double>> ubasis;
    ubasis.push_back(u1);
    ubasis.push_back(u2);
    ubasis.push_back(u3);
    ubasis.push_back(u4);

    std::vector<Vector<double>> standard_basis;
    standard_basis.push_back(e1);
    standard_basis.push_back(e2);
    standard_basis.push_back(e3);
    standard_basis.push_back(e4);

    std::vector<Vector<double>> gsd1 = gs(sub_basis);
    std::vector<Vector<double>> gsd2 = gs(ubasis);

    EXPECT_EQ(gsd1[0], v1);
    EXPECT_EQ(gsd1[1], res2);
    EXPECT_EQ(gsd2, standard_basis);
}

TEST(lll, LLL)
{
    Vector<mpq_class> v1(3);
    Vector<mpq_class> v2(3);
    Vector<mpq_class> v3(3);

    v1(0) = 0;
    v1(1) = 10;
    v1(2) = 2;

    v2(0) = 200;
    v2(1) = 314;
    v2(2) = 503;

    v3(0) = 20;
    v3(1) = 30;
    v3(2) = 50;

    std::vector<Vector<mpq_class>> basis;
    basis.push_back(v1);
    basis.push_back(v2);
    basis.push_back(v3);

    Vector<mpq_class> res(3);
    res(0) = 0;
    res(1) = 0;
    res(2) = 1;

    std::vector<Vector<mpq_class>> lll_basis = lll(0.75, basis);

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

    basis = std::vector<Vector<mpq_class>>();
    basis.push_back(v1);
    basis.push_back(v2);
    basis.push_back(v3);

    Vector<mpq_class> res1(3);
    res1(0) = 0;
    res1(1) = 1;
    res1(2) = 0;

    Vector<mpq_class> res2(3);
    res2(0) = 1;
    res2(1) = 0;
    res2(2) = 1;

    Vector<mpq_class> res3(3);
    res3(0) = -1;
    res3(1) = 0;
    res3(2) = 2;

    lll_basis = lll(0.75, basis);

    EXPECT_EQ(lll_basis[0], res1);
    EXPECT_EQ(lll_basis[1], res2);
    EXPECT_EQ(lll_basis[2], res3);

    Vector<mpq_class> u1(3);
    Vector<mpq_class> u2(3);
    Vector<mpq_class> u3(3);

    u1(0) = 69;
    u1(1) = 26;
    u1(2) = 15;

    u2(0) = 21;
    u2(1) = 96;
    u2(2) = 85;

    u3(0) = 333;
    u3(1) = -9476;
    u3(2) = -8694;

    std::vector<Vector<mpq_class>> new_basis;
    new_basis.push_back(u1);
    new_basis.push_back(u2);
    new_basis.push_back(u3);

    std::vector<Vector<mpq_class>> new_lll_basis = lll(0.75, new_basis);

    EXPECT_EQ(new_lll_basis[0], res);
}
