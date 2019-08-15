// lll implementation
#include "lll.h"

// Helping function for the lll algorithm
// recomputes both the orthogonal basis and the coefficient matrix needed in the implementation of the algorithm
// Is given a (usually linearly independent) basis on which it computes gram - schmidt, and a matrix in which to write the computed coefficients
std::vector<Vector<mpq_class>> update(unsigned int size, Matrix<mpq_class>& m, const std::vector<Vector<mpq_class>>& basis)
{
    std::vector<Vector<mpq_class>> ortho_basis = gs(basis);

    for (unsigned int i = 0; i < size; ++i) {
        for (unsigned int j = 0; j < size; ++j) {
            m(i, j) = (basis[i] * ortho_basis[j]) / (ortho_basis[j] * ortho_basis[j]);
        }
    }

    return ortho_basis;
}

// The naive implementation of the LLL algorithm, as can be found in the internet
// Is given a (usually linearly independent) basis to reduce and a precision parameter delta, assumed to be in (0.25, 1)
std::vector<Vector<mpq_class>> lll(double delta, const std::vector<Vector<mpq_class>>& basis)
{
    Vector<mpq_class> temp;

    std::vector<Vector<mpq_class>> reduced_basis = basis;
    unsigned int size = basis.size();

    Matrix<mpq_class> m(size, size);
    std::vector<Vector<mpq_class>> ortho_basis = update(size, m, reduced_basis);

    unsigned int k = 1;
    while (k <= size - 1) {
        for (int j = k - 1; j >= 0; --j) {
            if (abs(m(k, j)) > 0.5) {
                reduced_basis[k] -= (mpq_class)floor((mpf_class)(0.5 + m(k, j))) * reduced_basis[j];
                ortho_basis = update(size, m, reduced_basis);
            }
        }
        if (ortho_basis[k] * ortho_basis[k] >= (delta - m(k, k - 1) * m(k, k - 1)) * (ortho_basis[k - 1] * ortho_basis[k - 1])) {
            ++k;
        }
        else {
            temp = reduced_basis[k];
            reduced_basis[k] = reduced_basis[k - 1];
            reduced_basis[k - 1] = temp;

            ortho_basis = update(size, m, reduced_basis);

            if (k - 1 > 1)
                --k;
            else
                k = 1;
        }
    }

    return reduced_basis;
}
