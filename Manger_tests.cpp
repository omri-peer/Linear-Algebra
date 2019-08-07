#include "Manger.h"
#include <gtest/gtest.h>

TEST(manger, attacker)
{
    int n = 5897 * 6133;
    int e = 5;
    MangerAttacker ma(n, e);
    std::cout << ma.find_plaintext(14508243) << "\n";
}
