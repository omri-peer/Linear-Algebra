#include "lll.h"

// Helping function for the lll algorithm
// recomputes both the orthogonal basis and the coefficient matrix needed in the implementation of the algorithm
// Is given a (usually linearly independent) basis on which it computes gram - schmidt, and a matrix in which to write the computed coefficients
std::vector<Vector<mpq_class>> update(unsigned int size, Matrix<mpq_class>& m, const std::vector<Vector<mpq_class>>& basis)
{
    std::vector<Vector<mpq_class>> ortho = gs(basis);

    for (unsigned int i = 0; i < size; i++) {
        for (unsigned int j = 0; j < size; j++) {
            m(i, j) = (basis[i] * ortho[j]) / (ortho[j] * ortho[j]);
        }
    }

    return ortho;
}

// The naive implementation of the LLL algorithm, as can be found in the internet
// Is given a (usually linearly independent) basis to reduce and a precision parameter delta, assumed to be in (0.25, 1)
std::vector<Vector<mpq_class>> lll(double delta, const std::vector<Vector<mpq_class>>& basis)
{
    Vector<mpq_class> temp;

    std::vector<Vector<mpq_class>> our_basis = basis;
    unsigned int size = basis.size();

    Matrix<mpq_class> m(size, size);
    std::vector<Vector<mpq_class>> ortho = update(size, m, our_basis);

    unsigned int k = 1;
    while (k <= size - 1) {
        for (int j = k - 1; j >= 0; --j) {
            if (abs(m(k, j)) > 0.5) {
                our_basis[k] -= (mpq_class)floor((mpf_class)(0.5 + m(k, j))) * our_basis[j];
                ortho = update(size, m, our_basis);
            }
        }
        if (ortho[k] * ortho[k] >= (delta - m(k, k - 1) * m(k, k - 1)) * (ortho[k - 1] * ortho[k - 1])) {
            k++;
        }
        else {
            temp = our_basis[k];
            our_basis[k] = our_basis[k - 1];
            our_basis[k - 1] = temp;

            ortho = update(size, m, our_basis);

            if (k - 1 > 1)
                k--;
            else
                k = 1;
        }
    }

    return our_basis;
}