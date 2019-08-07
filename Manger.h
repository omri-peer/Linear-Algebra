#pragma once

#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include <iostream>

class MangerAttacker {
private:
    mpz_class n;
    mpz_class e;
    mpz_class B;

    static mpz_class ceil_div(const mpz_class& a, const mpz_class& b)
    {
        mpz_class cand = a / b;
        if (cand * b == a) return cand;
        return cand + 1;
    }

    bool oracle(const mpz_class& f, const mpz_class& c)
    {
        return ((123456 * f) % n) < B;
    }

public:
    explicit MangerAttacker(mpz_class n, mpz_class e) : n(std::move(n)), e(std::move(e)), B((mpz_class)1 << this->n.get_str(2).length() / 8 * 8)
    {
    }

    mpz_class find_plaintext(const mpz_class& c)
    {
        // Step 1
        mpz_class f1 = 2;
        mpz_class j = 0;
        while (oracle(f1, c)) {
            f1 *= 2;
            j++;
        }

        // Step 2
        mpz_class f2 = ((n + B) / B) * (f1 / 2);
        while (!oracle(f2, c)) {
            f2 += f1 / 2;
        }

        // Step 3
        mpz_class mmin = ceil_div(n, f2);
        mpz_class mmax = (n + B) / f2;
        mpz_class ftmp;
        mpz_class i;
        mpz_class f3;
        while (mmax - mmin > 0) {
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