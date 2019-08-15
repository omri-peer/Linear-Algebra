# Linear Algebra
Cryptography Attacks Workshop

This contains some basic linear algebra code, and implementations of the LLL algorithm and Manger's attack.

More specifically, two main classes implement the linear algebra functionality:
Matrix<T>, Vector<T> (implemented in Matrix.h and Vector.h).

Both classes support basic algebraic operations and some useful functionality,
and have some tests that make sure they indeed work (matrix_tests.cpp, vector_tests.cpp).

In addition, tests of the lll and Manger's algorithms are implemented in lll_tests.cpp and Manger_tests.cpp.

Note: the Matrix<T> class overloads the * operator with the naive implementation of matrices  mutiplication.
However, the class contains another, more efficient, implementation of multiplication, using an optimized version of Strassen's algorithm.

Build:
Navigate to the project's directoryת and run the commands:

cd cmake-build-debug
make

Test:
Navigate to the project's directory and run the commands:

cd cmake-build-debug
ctest -V
