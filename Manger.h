// Manger's attack
#pragma once

#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>

class MangerAttacker {
private:
    // the GMP library handles large integers and fractions with infinite precision
    mpz_class n; // public key
    mpz_class e; // public key
    mpz_class B; // B is the minimal number with the binary length of n

    // returns the quotient of the two given integers, rounded up to the nearest integer from above
    static mpz_class ceil_div(const mpz_class& a, const mpz_class& b)
    {
        mpz_class cand = a / b;         // the quotient rounded down (floor)
        if (cand * b == a) return cand; // if the division is actually accurate
        return cand + 1;                // round up instead
    }

    // given a ciphertext c and a coefficient f, returns whether or not (f*c^d mod n) is smaller than B.
    // An actual implementation of this function would probably involve some kind of side-channel attack,
    // and will not be given here. Instead, an implementation for the specific case: m = 123456 mod n is given.
    bool oracle(const mpz_class& f, const mpz_class& c)
    {
        return ((123456 * f) % n) < B;
    }

public:
    explicit MangerAttacker(mpz_class n, mpz_class e) : n(std::move(n)), e(std::move(e)), B((mpz_class)1 << this->n.get_str(2).length() / 8 * 8)
    {
    }

    // given a ciphertext c performs some oracle calls with appropriate analysis until finding the corresponding plaintext
    mpz_class find_plaintext(const mpz_class& c)
    {
        // Step 1
        // find f1 s.t. B/2 <= (f1/2)*m < B
        mpz_class f1 = 2;
        mpz_class j = 0;
        while (oracle(f1, c)) {
            f1 *= 2;
            j++;
        }

        // now B/2 <= (f1/2)*m < B
        //
        // Step 2
        // find f2 s.t. n <= f2*m < n + B
        mpz_class f2 = ((n + B) / B) * (f1 / 2);
        while (!oracle(f2, c)) {
            f2 += f1 / 2;
        }

        // now n <= f2*m < n + B
        // Step 3
        // we now know that n/f2 <= m <= (n+B)/f2, and can narrow it down (by roughly half at a time)
        mpz_class mmin = ceil_div(n, f2);
        mpz_class mmax = (n + B) / f2;
        mpz_class ftmp;
        mpz_class i;
        mpz_class f3;
        while (mmax - mmin > 0) { // until there is only one option
            ftmp = (2 * B) / (mmax - mmin);
            i = ftmp * mmin / n;
            f3 = ceil_div(i * n, mmin);

            if (!oracle(f3, c)) {
                mmin = ceil_div(i * n + B, f3);
            }
            else {
                mmax = (i * n + B) / f3;
            }
        }

        return mmin;
    }
};
