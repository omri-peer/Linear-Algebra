#pragma once

#include "Matrix.h"
#include "Vector.h"
#include <cmath>
#include <iostream>
#include <vector>

// Helping function to the lll algorithm
std::vector<Vector<double>> update(int size, Matrix<double>& m, const std::vector<Vector<double>>& basis);

std::vector<Vector<double>> lll(double delta, const std::vector<Vector<double>>& basis);