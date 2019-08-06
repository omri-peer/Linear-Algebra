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
        /*
        std::cout << f << ", " << c << "\n";
        unsigned int d = 463;
        unsigned int result = 1;
        unsigned int base = (int)(pow(f, e) * c) % n;
        std::cout << (pow(f, e) * c) << " " << (unsigned int)(pow(f, e) * c) << "\n";
        unsigned int exp = d;
        while (exp > 0) {
            if (exp % 2 == 1) {
                std::cout << "#@$$#@##@ " << result << " " << base << " " << result * base << " " << (result * base) % n << "\n";
                result = (result * base) % n;
            }
            exp = exp >> 1;
            base = (base * base) % n;
            std::cout << "result = " << result << ", base = " << base << "\n";
        }

        std::cout << "final result = " << result << "\n";*/
        return (123456 * f) % n < B;
    }

    long find_plaintext(int c)
    {
        std::cout << "n = " << n << ", e = " << e << ", k = " << k << ", B = " << B << "\n";
        // Step 1
        long f1 = 2;
        long j = 0;
        while (oracle(f1, c)) {
            f1 *= 2;
            j++;
            if (j == 10) break;
        }
        // Step 2
        std::cout << "f1 = " << f1 << "\n";
        long f2 = floor((n + B) / (double)B) * (f1 / 2);
        while (!oracle(f2, c)) {
            f2 += f1 / 2;
        }
        std::cout << "f2 = " << f2 << "\n";
        // Step 3
        long mmin = ceil(n / (double)f2);
        long mmax = floor((n + B) / (double)f2);
        long ftmp;
        long i;
        long f3;
        while (mmax - mmin > 0) {
            std::cout << "mmin = " << mmin << ", mmax = " << mmax << "\n";
            ftmp = floor((2 * B) / (double)(mmax - mmin));
            i = floor(ftmp * mmin / (double)n);
            std::cout << "i = " << i << "\n";
            f3 = ceil(i * n / (double)mmin);
            std::cout << "f3 = " << f3 << "\n";

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