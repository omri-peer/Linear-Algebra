#include "Vector.h"

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
