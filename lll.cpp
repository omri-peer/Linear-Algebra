#include "lll.h"

// Helping function to the lll algorithm
std::vector<Vector<mpq_class>> update(int size, Matrix<mpq_class>& m, const std::vector<Vector<mpq_class>>& basis)
{
    std::vector<Vector<mpq_class>> ortho = gs(basis);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            m(i, j) = (basis[i] * ortho[j]) / (ortho[j] * ortho[j]);
        }
    }

    return ortho;
}

std::vector<Vector<mpq_class>> lll(double delta, const std::vector<Vector<mpq_class>>& basis)
{
    Vector<mpq_class> temp;

    std::vector<Vector<mpq_class>> our_basis = basis;
    int size = basis.size();

    Matrix<mpq_class> m(size, size);
    std::vector<Vector<mpq_class>> ortho = update(size, m, our_basis);

    int k = 1;
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