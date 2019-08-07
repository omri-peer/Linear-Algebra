#include "Manger.h"
#include <gtest/gtest.h>

TEST(manger, attacker)
{
    mpz_class n = 5897 * 6133;
    mpz_class e = 5;
    MangerAttacker ma1(n, e);
    EXPECT_EQ(ma1.find_plaintext(14508243), 123456);

    n = 13513 * 19423;
    e = 5;
    MangerAttacker ma2(n, e);
    EXPECT_EQ(ma2.find_plaintext(57746787), 123456);

    n = mpz_class("185652620748269215170218178196471833119") * mpz_class("295937400278117939395781398270220581033");
    e = 5;
    MangerAttacker ma3(n, e);
    EXPECT_EQ(ma3.find_plaintext(mpz_class("28678802168634497644363776")), mpz_class("123456"));
}
