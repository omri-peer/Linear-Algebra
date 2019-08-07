#include "Manger.h"
#include <gtest/gtest.h>

TEST(manger, attacker)
{
    long n = 5897 * 6133;
    long e = 5;
    MangerAttacker ma(n, e);
    EXPECT_EQ(ma.find_plaintext(14508243), 123456);

    n = 13513 * 19423;
    e = 5;
    ma = MangerAttacker(n, e);
    EXPECT_EQ(ma.find_plaintext(57746787), 123456);
}
