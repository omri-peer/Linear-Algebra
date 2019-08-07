#pragma once

#include "Matrix.h"
#include "Vector.h"
#include <cmath>
#include <iostream>
#include <vector>

std::vector<Vector<double>> gs(const std::vector<Vector<double>>& basis)
{
    std::vector<Vector<double>> new_basis;
    for (int i = 0; i < basis.size(); i++) {
        new_basis.push_back(basis[i]);
        for (int j = 0; j < i; j++) {
            new_basis[i] -= Vector<double>(new_basis[j] * (basis[i] * new_basis[j] / (new_basis[j] * new_basis[j])));
        }
    }

    return new_basis;
}

// Helping function to the lll algorithm
std::vector<Vector<double>> update(int size, Matrix<double>& m, const std::vector<Vector<double>>& basis)
{
    std::vector<Vector<double>> ortho = gs(basis);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            m(i, j) = (basis[i] * ortho[j]) / (ortho[j] * ortho[j]);
        }
    }

    return ortho;
}

std::vector<Vector<double>> lll(double delta, const std::vector<Vector<double>>& basis)
{
    Vector<double> temp;

    std::vector<Vector<double>> our_basis = basis;
    int size = basis.size();

    Matrix<double> m(size, size);
    std::vector<Vector<double>> ortho = update(size, m, our_basis);

    int k = 1;
    while (k <= size - 1) {
        for (int j = k - 1; j >= 0; --j) {
            if (std::abs(m(k, j)) > 0.5) {
                our_basis[k] -= round(m(k, j)) * our_basis[j];
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