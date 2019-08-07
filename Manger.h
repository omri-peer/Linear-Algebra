#pragma once

#include <cmath>
#include <iostream>

class MangerAttacker {
private:
    int n;
    int e;
    int k;
    int B;

public:
    explicit MangerAttacker(long n, long e) : n(n), e(e), k(ceil(log(n) / log(256))), B(pow(256, k - 1))
    {
    }
    bool oracle(long f, long c)
    {
        return (123456 * f) % n < B;
    }

    long find_plaintext(long c)
    {
        // Step 1
        long f1 = 2;
        long j = 0;
        while (oracle(f1, c)) {
            f1 *= 2;
            j++;
            if (j == 10) break;
        }
        // Step 2
        long f2 = floor((n + B) / (double)B) * (f1 / 2);
        while (!oracle(f2, c)) {
            f2 += f1 / 2;
        }
        // Step 3
        long mmin = ceil(n / (double)f2);
        long mmax = floor((n + B) / (double)f2);
        long ftmp;
        long i;
        long f3;
        while (mmax - mmin > 0) {
            ftmp = floor((2 * B) / (double)(mmax - mmin));
            i = floor(ftmp * mmin / (double)n);
            f3 = ceil(i * n / (double)mmin);

            if (!oracle(f3, c)) {
                mmin = ceil((i * n + B) / (double)f3);
            }
            else {
                mmax = floor((i * n + B) / (double)f3);
            }
        }

        return mmin;
    }
};